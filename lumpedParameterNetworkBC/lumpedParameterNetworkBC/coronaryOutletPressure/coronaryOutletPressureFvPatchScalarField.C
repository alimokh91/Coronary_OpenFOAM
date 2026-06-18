#include "coronaryOutletPressureFvPatchScalarField.H"

#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "fvPatchFieldMacros.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "Time.H"
#include "IFstream.H"
#include "DynamicList.H"
#include "SortableList.H"
#include "IOmanip.H"

#include <cmath>

using namespace Foam;

//---------------------------------------------------------------------------//
//                              Constructors                                 //
//---------------------------------------------------------------------------//

coronaryOutletPressureFvPatchScalarField::
coronaryOutletPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{}


coronaryOutletPressureFvPatchScalarField::
coronaryOutletPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict)
{
    // Parameters
    dict.lookup("Ra")  >> Ra_;
    dict.lookup("Ram") >> Ram_;
    dict.lookup("Rv")  >> Rv_;
    dict.readIfPresent("Rvm", Rvm_);

    dict.lookup("Ca")  >> Ca_;
    dict.lookup("Cim") >> Cim_;

    dict.lookup("PimFile") >> PimFile_;
    dict.readIfPresent("PimPeriod",  PimPeriod_);
    dict.readIfPresent("PimScaling", PimScaling_);

    word schemeName("BDF1");
    dict.readIfPresent("diffScheme", schemeName);
    scheme_ = (schemeName == "BDF2" || schemeName == "bdf2") ? BDF2 : BDF1;

    dict.readIfPresent("useLaggedFlux", useLaggedFlux_);
    dict.readIfPresent("relaxFactor",   relaxFactor_);
    dict.readIfPresent("flipFluxSign",  flipFluxSign_);
    dict.readIfPresent("rho",           rho_);
    dict.readIfPresent("regularize",    regularize_);

    dict.readIfPresent("pimMode", pimMode_);
    dict.readIfPresent("transmuralOnCa", transmuralOnCa_);
    dict.readIfPresent("firstCycleRamp", firstCycleRamp_);

    // NEW: backflow penalty controls (opt-in)
    dict.readIfPresent("backflowStabilization", backflowStabilization_);
    dict.readIfPresent("penaltyFactor", penaltyFactor_);
    dict.readIfPresent("smoothingWidth", smoothingWidth_);
    dict.readIfPresent("U", UName_);

    haveWarmStart_ = dict.found("Po") || dict.found("Pio") || dict.found("Qo");
    if (haveWarmStart_)
    {
        dict.readIfPresent("Po",  Po0_);
        dict.readIfPresent("Pio", Pio0_);
        dict.readIfPresent("Qo",  Qo0_);
    }

    // --- self-calibrating sign defaults (OPTIONAL persistence) ---
    qSign_      = 1.0;
    signLocked_ = false;
    dict.readIfPresent("qSign",      qSign_);
    dict.readIfPresent("signLocked", signLocked_);

    // States
    Pn_      = haveWarmStart_ ? Po0_  : 0.0;
    Pnm1_    = Pn_;
    Pio_     = haveWarmStart_ ? Pio0_ : 0.0;
    Pio_m1_  = Pio_;
    Qn_      = haveWarmStart_ ? Qo0_  : 0.0;
    Qm1_     = Qn_;
    dtPrev_  = 0.0;
    firstIter_= true;
    pimHasPrev_ = false;

    readPimTable(PimFile_);

    // Start boundary at Po/ρ to avoid a big first-corrector jump
    this->operator==(Pn_/(rho_ + SMALL));
}


coronaryOutletPressureFvPatchScalarField::
coronaryOutletPressureFvPatchScalarField
(
    const coronaryOutletPressureFvPatchScalarField& rhs,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& m
)
:
    fixedValueFvPatchScalarField(rhs, p, iF, m)
{
    Ra_ = rhs.Ra_; Ram_ = rhs.Ram_; Rv_ = rhs.Rv_; Rvm_ = rhs.Rvm_;
    Ca_ = rhs.Ca_; Cim_ = rhs.Cim_;
    rho_ = rhs.rho_; regularize_ = rhs.regularize_;
    scheme_ = rhs.scheme_; useLaggedFlux_ = rhs.useLaggedFlux_;
    relaxFactor_ = rhs.relaxFactor_; flipFluxSign_ = rhs.flipFluxSign_;
    PimFile_ = rhs.PimFile_; PimPeriod_ = rhs.PimPeriod_;
    PimScaling_ = rhs.PimScaling_; pimMode_ = rhs.pimMode_;
    transmuralOnCa_ = rhs.transmuralOnCa_;
    firstCycleRamp_ = rhs.firstCycleRamp_;
    backflowStabilization_ = rhs.backflowStabilization_;
    penaltyFactor_ = rhs.penaltyFactor_;
    smoothingWidth_ = rhs.smoothingWidth_;
    UName_ = rhs.UName_;
    PimTable_ = rhs.PimTable_;
    Pn_ = rhs.Pn_; Pnm1_ = rhs.Pnm1_;
    Pio_ = rhs.Pio_; Pio_m1_ = rhs.Pio_m1_;
    Qn_ = rhs.Qn_; Qm1_ = rhs.Qm1_;
    dtPrev_ = rhs.dtPrev_; firstIter_ = rhs.firstIter_;
    haveWarmStart_ = rhs.haveWarmStart_;
    Po0_ = rhs.Po0_; Pio0_ = rhs.Pio0_; Qo0_ = rhs.Qo0_;
    pimPrev_ = rhs.pimPrev_; pimHasPrev_ = rhs.pimHasPrev_;

    // copy sign calibration state
    qSign_      = rhs.qSign_;
    signLocked_ = rhs.signLocked_;
}


coronaryOutletPressureFvPatchScalarField::
coronaryOutletPressureFvPatchScalarField
(
    const coronaryOutletPressureFvPatchScalarField& rhs,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(rhs, iF)
{
    Ra_ = rhs.Ra_; Ram_ = rhs.Ram_; Rv_ = rhs.Rv_; Rvm_ = rhs.Rvm_;
    Ca_ = rhs.Ca_; Cim_ = rhs.Cim_;
    rho_ = rhs.rho_; regularize_ = rhs.regularize_;
    scheme_ = rhs.scheme_; useLaggedFlux_ = rhs.useLaggedFlux_;
    relaxFactor_ = rhs.relaxFactor_; flipFluxSign_ = rhs.flipFluxSign_;
    PimFile_ = rhs.PimFile_; PimPeriod_ = rhs.PimPeriod_;
    PimScaling_ = rhs.PimScaling_; pimMode_ = rhs.pimMode_;
    transmuralOnCa_ = rhs.transmuralOnCa_;
    firstCycleRamp_ = rhs.firstCycleRamp_;
    backflowStabilization_ = rhs.backflowStabilization_;
    penaltyFactor_ = rhs.penaltyFactor_;
    smoothingWidth_ = rhs.smoothingWidth_;
    UName_ = rhs.UName_;
    PimTable_ = rhs.PimTable_;
    Pn_ = rhs.Pn_; Pnm1_ = rhs.Pnm1_;
    Pio_ = rhs.Pio_; Pio_m1_ = rhs.Pio_m1_;
    Qn_ = rhs.Qn_; Qm1_ = rhs.Qm1_;
    dtPrev_ = rhs.dtPrev_; firstIter_ = rhs.firstIter_;
    haveWarmStart_ = rhs.haveWarmStart_;
    Po0_ = rhs.Po0_; Pio0_ = rhs.Pio0_; Qo0_ = rhs.Qo0_;
    pimPrev_ = rhs.pimPrev_; pimHasPrev_ = rhs.pimHasPrev_;

    // copy sign calibration state
    qSign_      = rhs.qSign_;
    signLocked_ = rhs.signLocked_;
}

//---------------------------------------------------------------------------//
//                           Pim table reader                                //
//---------------------------------------------------------------------------//

void coronaryOutletPressureFvPatchScalarField::readPimTable(const fileName& path)
{
    PimTable_.clear();
    IFstream is(path);
    if (!is.good())
    {
        FatalIOErrorInFunction(is)
            << "Cannot open Pim file: " << path << nl
            << exit(FatalIOError);
    }

    DynamicList<std::pair<scalar, scalar>> tmp;
    token tok;
    bool insideOuterList = false;

    while (is.good())
    {
        is.read(tok);
        if (!is.good()) break;

        if (tok == token::BEGIN_LIST)
        {
            token next;
            is.read(next);
            if (next == token::BEGIN_LIST)
            {
                insideOuterList = true;
                scalar t(0), p(0);
                is >> t >> p;
                is.read(tok); // consume ')'
                tmp.append({t, p});
            }
            else if (next.isScalar() || next.isLabel())
            {
                scalar t = next.scalarToken();
                scalar p(0);
                is >> p;
                is.read(tok); // consume ')'
                tmp.append({t, p});
            }
            else
            {
                insideOuterList = true; // outer list wrapper
            }
        }
        else if (tok == token::END_LIST)
        {
            if (insideOuterList) break;
        }
        else if (tok.isScalar() || tok.isLabel())
        {
            scalar t = tok.scalarToken();
            scalar p(0); is >> p;
            tmp.append({t, p});
        }
    }

    if (tmp.size() < 2)
    {
        FatalErrorInFunction
            << "Pim file '" << path << "' must contain at least two (time pressure) pairs"
            << exit(FatalError);
    }

    std::sort
    (
        tmp.begin(), tmp.end(),
        [](const std::pair<scalar,scalar>& a, const std::pair<scalar,scalar>& b)
        { return a.first < b.first; }
    );

    PimTable_.assign(tmp.begin(), tmp.end());
}

//---------------------------------------------------------------------------//
//                        Pim interpolation                                  //
//---------------------------------------------------------------------------//

scalar coronaryOutletPressureFvPatchScalarField::interpolatePim(const scalar t) const
{
    if (PimTable_.empty()) return 0.0;

    const scalar t0 = PimTable_.front().first;
    const scalar t1 = PimTable_.back().first;
    scalar teff = t;

    if (PimPeriod_ > SMALL)
    {
        const scalar s = max(t - t0, 0.0);
        if (firstCycleRamp_)
        {
            // First cycle: play [t0, t0+T] once (startup ramp). Afterwards: loop the
            // SECOND cycle [t0+T, t0+2T] (settled waveform).
            if (s < PimPeriod_)
            {
                teff = t0 + s;
            }
            else
            {
                teff = t0 + PimPeriod_ + std::fmod(s - PimPeriod_, PimPeriod_);
            }
        }
        else
        {
            teff = t0 + std::fmod(s, PimPeriod_);
        }
    }
    else
    {
        if (teff <= t0) return PimScaling_ * PimTable_.front().second;
        if (teff >= t1) return PimScaling_ * PimTable_.back().second;
    }

    for (std::size_t i = 1; i < PimTable_.size(); ++i)
    {
        const scalar ta = PimTable_[i-1].first;
        const scalar tb = PimTable_[i].first;
        if (teff >= ta && teff <= tb)
        {
            const scalar ya = PimTable_[i-1].second;
            const scalar yb = PimTable_[i].second;
            const scalar a = (tb - ta > SMALL) ? (teff - ta)/(tb - ta) : 0.0;
            return PimScaling_ * ((1.0 - a)*ya + a*yb);
        }
    }

    return PimScaling_ * PimTable_.back().second;
}

//---------------------------------------------------------------------------//
//                           updateCoeffs                                    //
//---------------------------------------------------------------------------//

void coronaryOutletPressureFvPatchScalarField::updateCoeffs()
{
    if (updated()) return;

    const Time& runTime = this->db().time();
    const scalar dt = runTime.deltaTValue();
    const scalar t  = runTime.value();

    const surfaceScalarField& phi =
        this->db().lookupObject<surfaceScalarField>("phi");
    const label patchI = this->patch().index();
    const fvsPatchField<scalar>& phip = phi.boundaryField()[patchI];

    // OpenFOAM: sum(phi) > 0 => out of 3D domain (into LPN)
    scalar Qinst = gSum(phip);

    // If 'phi' is mass flux (kg/s), convert to volumetric using rho
    if (phi.dimensions()[0] != 0)        // mass exponent nonzero -> kg/s
    {
        Qinst /= max(rho_, SMALL);       // m^3/s
    }

    // Legacy per-patch flip (dictionary) if someone set it
    if (flipFluxSign_) Qinst = -Qinst;

    // --------- ONE-TIME PASSIVITY SELF-CHECK (no damping) ---------
    // If raising boundary pressure increases outflow at this patch (dQ*dP>0),
    // the 3D-0D interface is non-passive. Flip sign once so dQ/dP < 0.
    if (!signLocked_)
    {
        const scalar dP = Pn_ - Pnm1_;
        const scalar dQ = Qinst - Qm1_;
        if (mag(dP) > 5.0 && mag(dQ) > 1e-10)
        {
            if (dP*dQ > 0.0) qSign_ = -1.0;
            signLocked_ = true;
        }
    }
    Qinst *= qSign_;
    // --------------------------------------------------------------

    const scalar Qin = useLaggedFlux_ ? Qm1_ : Qinst;
    const scalar Pim = interpolatePim(t);

    scalar dPimdt = 0.0;
    if (pimHasPrev_) dPimdt = (Pim - pimPrev_) / max(dt, SMALL);

    const scalar invRam = Ram_ > SMALL ? 1.0/Ram_ : 0.0;
    const scalar invRv  = Rv_  > SMALL ? 1.0/Rv_  : 0.0;
    const scalar invRvm = Rvm_ > SMALL ? 1.0/Rvm_ : 0.0;

    const bool variableDt = (dtPrev_ > SMALL) && (mag(dt - dtPrev_) > VSMALL);
    const bool useBdf2Now = (scheme_ == BDF2) && !firstIter_ && !variableDt && (dtPrev_ > SMALL);

    scalar PoNew = Pn_;
    scalar PioNew = Pio_;

    if (!useBdf2Now)
    {
        // ---- BDF1 ----
        const scalar a11 = Ca_/dt + invRam + invRv + regularize_;
        const scalar a12 = -invRam;
        const scalar a21 =  invRam;
        const scalar a22 = Cim_/dt + invRam + invRvm + regularize_;

        scalar rhs1 = Ca_/dt * Pn_   + Qin;
        scalar rhs2 = Cim_/dt * Pio_ + invRvm * Pim;

        if (pimMode_ == "backPressure")
        {
            rhs1 += Pim * invRv;
        }
        else if (pimMode_ == "transmural")
        {
            rhs1 += Pim * invRv;      // venous back-pressure on Rv
            rhs2 += Cim_ * dPimdt;    // dPim/dt on Cim
            if (transmuralOnCa_) rhs1 += Ca_ * dPimdt;  // optional on Ca
        }

        //const scalar det = a11*a22 - a12*a21;
	const scalar A11 = max(a11, SMALL);
	const scalar A22 = max(a22, SMALL);
	const scalar det = A11*A22 - a12*a21;
	const scalar dPo_dQ = A22 / max(det, SMALL);
	if (dPo_dQ <= SMALL)
	{
    		FatalErrorInFunction
       		 << "Non-dissipative LPN coupling (dPo/dQ <= 0) at patch '"
        	 << patch().name() << "'. "
                 << "a11="<<a11<<", a22="<<a22<<", det="<<det << nl
                 << "Check signs/coefficients and optional entries."
                 << exit(FatalError);
	}
	if (mag(det) <= SMALL)
        {
            PoNew  = rhs1 / max(A11, VSMALL);
            PioNew = (rhs2 - invRam*PoNew) / max(A22, VSMALL);
        }
        else
        {
            PoNew  = ( rhs1*A22 - a12*rhs2)/det;
            PioNew = (-rhs1*a21 + A11*rhs2)/det;
        }
    }
    else
    {
        // ---- BDF2 ----  (used only when dt is constant)
        const scalar beta = 3.0/(2.0*dt);
        const scalar a11 = Ca_*beta + invRam + invRv + regularize_;
        const scalar a12 = -invRam;
        const scalar a21 =  invRam;
        const scalar a22 = Cim_*beta + invRam + invRvm + regularize_;

        scalar rhs1 = Ca_*beta*(4.0*Pn_  - Pnm1_)    + Qin;
        scalar rhs2 = Cim_*beta*(4.0*Pio_ - Pio_m1_) + invRvm * Pim;

        if (pimMode_ == "backPressure")
        {
            rhs1 += Pim * invRv;
        }
        else if (pimMode_ == "transmural")
        {
            rhs1 += Pim * invRv;

            // Proper BDF2 for dPim/dt using (n+1, n, n-1) with constant dt
            const scalar Pim_c   = Pim;
            const scalar Pim_n   = interpolatePim(max(t - dt, scalar(0)));
            const scalar Pim_nm1 = interpolatePim(max(t - 2.0*dt, scalar(0)));
            const scalar dPimdtBDF2 = (3.0*Pim_c - 4.0*Pim_n + Pim_nm1)/(2.0*dt);

            rhs2 += Cim_ * dPimdtBDF2;
            if (transmuralOnCa_) rhs1 += Ca_ * dPimdtBDF2;
        }

        const scalar det = a11*a22 - a12*a21;
        if (mag(det) <= SMALL)
        {
            PoNew  = rhs1 / max(a11, VSMALL);
            PioNew = (rhs2 - invRam*PoNew) / max(a22, VSMALL);
        }
        else
        {
            PoNew  = ( rhs1*a22 - a12*rhs2)/det;
            PioNew = (-rhs1*a21 + a11*rhs2)/det;
        }
    }

    // Apply: incompressible kinematic pressure p = P/ρ
    const scalar pk = PoNew / (rho_ + SMALL);
    if (!std::isfinite(pk))   // <— add std::
	{
    FatalErrorInFunction
        << "Non-finite patch pressure at '" << patch().name()
        << "' (pk=" << pk << ")"
        << exit(FatalError);
	}
    
    scalarField newVals(this->size(), pk);

    // --- Opt-in backflow penalty: on reversed faces add a kinematic dynamic
    //     pressure  mask*penaltyFactor*0.5*|U|^2  (no rho: field is kinematic).
    if (backflowStabilization_)
    {
        const fvPatchField<vector>& Up =
            this->db().lookupObject<volVectorField>(UName_).boundaryField()[patchI];

        scalarField mask(this->size(), Zero);
        if (smoothingWidth_ > SMALL)
        {
            const label nF = returnReduce(phip.size(), sumOp<label>());
            const scalar phiRef = gSum(mag(phip))/max(scalar(nF), scalar(1));
            const scalar sw = max(smoothingWidth_*phiRef, SMALL);
            forAll(mask, i) mask[i] = 0.5*(1.0 - Foam::tanh(phip[i]/sw));
        }
        else
        {
            forAll(mask, i) mask[i] = pos0(-phip[i] - SMALL);
        }

        forAll(newVals, i)
        {
            newVals[i] = pk + mask[i]*penaltyFactor_*0.5*magSqr(Up[i]);
        }
    }

    if (relaxFactor_ < 1.0 - SMALL)
    {
        newVals = relaxFactor_*newVals + (1.0 - relaxFactor_)*(*this);
    }
    this->operator==(newVals);

    // Bookkeeping
    Pnm1_   = Pn_;    Pn_    = PoNew;
    Pio_m1_ = Pio_;   Pio_   = PioNew;
    Qm1_    = Qn_;    Qn_    = Qinst;

    pimPrev_    = Pim;
    pimHasPrev_ = true;

    dtPrev_   = dt;
    firstIter_= false;

    fixedValueFvPatchScalarField::updateCoeffs();
}

//---------------------------------------------------------------------------//
//                                  write                                    //
//---------------------------------------------------------------------------//

void coronaryOutletPressureFvPatchScalarField::write(Ostream& os) const
{
    fixedValueFvPatchScalarField::write(os);
    os.writeKeyword("type") << this->type() << token::END_STATEMENT << nl;

    os.writeKeyword("Ra")  << Ra_  << token::END_STATEMENT << nl;
    os.writeKeyword("Ram") << Ram_ << token::END_STATEMENT << nl;
    os.writeKeyword("Rv")  << Rv_  << token::END_STATEMENT << nl;
    os.writeKeyword("Rvm") << Rvm_ << token::END_STATEMENT << nl;
    os.writeKeyword("Ca")  << Ca_  << token::END_STATEMENT << nl;
    os.writeKeyword("Cim") << Cim_ << token::END_STATEMENT << nl;

    os.writeKeyword("PimFile") << PimFile_ << token::END_STATEMENT << nl;
    os.writeKeyword("PimPeriod") << PimPeriod_ << token::END_STATEMENT << nl;
    os.writeKeyword("PimScaling") << PimScaling_ << token::END_STATEMENT << nl;

    os.writeKeyword("diffScheme") << (scheme_==BDF2 ? "BDF2" : "BDF1") << token::END_STATEMENT << nl;
    os.writeKeyword("useLaggedFlux") << useLaggedFlux_ << token::END_STATEMENT << nl;
    os.writeKeyword("relaxFactor") << relaxFactor_ << token::END_STATEMENT << nl;
    os.writeKeyword("flipFluxSign") << flipFluxSign_ << token::END_STATEMENT << nl;

    os.writeKeyword("rho") << rho_ << token::END_STATEMENT << nl;
    os.writeKeyword("regularize") << regularize_ << token::END_STATEMENT << nl;

    os.writeKeyword("pimMode") << pimMode_ << token::END_STATEMENT << nl;
    os.writeKeyword("transmuralOnCa") << transmuralOnCa_ << token::END_STATEMENT << nl;
    os.writeKeyword("firstCycleRamp") << firstCycleRamp_ << token::END_STATEMENT << nl;

    // NEW: backflow penalty controls
    os.writeKeyword("backflowStabilization") << backflowStabilization_ << token::END_STATEMENT << nl;
    os.writeKeyword("penaltyFactor")  << penaltyFactor_  << token::END_STATEMENT << nl;
    os.writeKeyword("smoothingWidth") << smoothingWidth_ << token::END_STATEMENT << nl;
    os.writeKeyword("U")              << UName_          << token::END_STATEMENT << nl;

    // write back current state (useful for warm restarts)
    os.writeKeyword("Po")  << Pn_  << token::END_STATEMENT << nl;
    os.writeKeyword("Pio") << Pio_ << token::END_STATEMENT << nl;
    os.writeKeyword("Qo")  << Qn_  << token::END_STATEMENT << nl;

    // persist sign calibration across restarts (optional)
    os.writeKeyword("qSign")      << qSign_      << token::END_STATEMENT << nl;
    os.writeKeyword("signLocked") << signLocked_ << token::END_STATEMENT << nl;

    writeEntry("value", os);
}

//---------------------------------------------------------------------------//
//                              Registration                                 //
//---------------------------------------------------------------------------//

namespace Foam
{
    makePatchTypeField(fvPatchScalarField, coronaryOutletPressureFvPatchScalarField);
}

