#include "windkesselOutletPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::windkesselOutletPressureFvPatchScalarField::
windkesselOutletPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    Rp_(1), Rd_(1), C_(1), L_(1), Pd_(0), rho_(1060),
    Pooo_(0), Poo_(0), Po_(0), Pn_(0),
    Qooo_(0), Qoo_(0), Qo_(0), Qn_(0),
    windkesselModel_(WK3),
    diffScheme_(secondOrder),
    timeIndex_(-1),
    lastDt_(-1),
    relaxFactor_(0.5),
    useLaggedFlux_(true),
    flipFluxSign_(false),
    backflowStabilization_(false),
    penaltyFactor_(1.0),
    smoothingWidth_(0.1),
    UName_("U")
{}


Foam::windkesselOutletPressureFvPatchScalarField::
windkesselOutletPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    Rp_(dict.lookupOrDefault<scalar>("Rp", 1)),
    Rd_(dict.lookupOrDefault<scalar>("Rd", 1)),
    C_(dict.lookupOrDefault<scalar>("C", 1)),
    L_(dict.lookupOrDefault<scalar>("L", 1)),
    Pd_(dict.lookupOrDefault<scalar>("Pd", 0)),
    rho_(dict.lookupOrDefault<scalar>("rho", 1060)),
    Pooo_(dict.lookupOrDefault<scalar>("Pooo", 0)),
    Poo_(dict.lookupOrDefault<scalar>("Poo", 0)),
    Po_(dict.lookupOrDefault<scalar>("Po", 0)),
    Pn_(dict.lookupOrDefault<scalar>("Pn", 0)),
    Qooo_(dict.lookupOrDefault<scalar>("Qooo", 0)),
    Qoo_(dict.lookupOrDefault<scalar>("Qoo", 0)),
    Qo_(dict.lookupOrDefault<scalar>("Qo", 0)),
    Qn_(dict.lookupOrDefault<scalar>("Qn", 0)),
    timeIndex_(-1),
    lastDt_(-1),
    relaxFactor_(dict.lookupOrDefault<scalar>("relaxFactor", 0.5)),
    useLaggedFlux_(dict.lookupOrDefault<bool>("useLaggedFlux", true)),
    flipFluxSign_(dict.lookupOrDefault<bool>("flipFluxSign", false)),
    backflowStabilization_(dict.lookupOrDefault<bool>("backflowStabilization", false)),
    penaltyFactor_(dict.lookupOrDefault<scalar>("penaltyFactor", 1.0)),
    smoothingWidth_(dict.lookupOrDefault<scalar>("smoothingWidth", 0.1)),
    UName_(dict.lookupOrDefault<word>("U", "U"))
{
    // Windkessel model
    word WKModelStr = dict.lookupOrDefault<word>("windkesselModel", "WK3");
    if      (WKModelStr == "Resistive")   windkesselModel_ = Resistive;
    else if (WKModelStr == "WK2")         windkesselModel_ = WK2;
    else if (WKModelStr == "WK3")         windkesselModel_ = WK3;
    else if (WKModelStr == "WK4Series")   windkesselModel_ = WK4Series;
    else if (WKModelStr == "WK4Parallel") windkesselModel_ = WK4Parallel;
    else
    {
        FatalErrorInFunction
            << "Unknown Windkessel Model: " << WKModelStr << nl
            << "Valid options: Resistive WK2 WK3 WK4Series WK4Parallel"
            << exit(FatalError);
    }

    // Differencing scheme (note: secondOrder assumes constant deltaT)
    word schemeStr = dict.lookupOrDefault<word>("diffScheme", "secondOrder");
    if      (schemeStr == "firstOrder")  diffScheme_ = firstOrder;
    else if (schemeStr == "secondOrder") diffScheme_ = secondOrder;
    else
    {
        FatalErrorInFunction
            << "Unknown differencing scheme: " << schemeStr << nl
            << "Valid options: firstOrder secondOrder"
            << exit(FatalError);
    }

    // Clamp relax factor
    relaxFactor_ = max(0.0, min(1.0, relaxFactor_));
}


Foam::windkesselOutletPressureFvPatchScalarField::
windkesselOutletPressureFvPatchScalarField
(
    const windkesselOutletPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    Rp_(ptf.Rp_), Rd_(ptf.Rd_), C_(ptf.C_), L_(ptf.L_), Pd_(ptf.Pd_), rho_(ptf.rho_),
    Pooo_(ptf.Pooo_), Poo_(ptf.Poo_), Po_(ptf.Po_), Pn_(ptf.Pn_),
    Qooo_(ptf.Qooo_), Qoo_(ptf.Qoo_), Qo_(ptf.Qo_), Qn_(ptf.Qn_),
    windkesselModel_(ptf.windkesselModel_),
    diffScheme_(ptf.diffScheme_),
    timeIndex_(ptf.timeIndex_),
    lastDt_(ptf.lastDt_),
    relaxFactor_(ptf.relaxFactor_),
    useLaggedFlux_(ptf.useLaggedFlux_),
    flipFluxSign_(ptf.flipFluxSign_),
    backflowStabilization_(ptf.backflowStabilization_),
    penaltyFactor_(ptf.penaltyFactor_),
    smoothingWidth_(ptf.smoothingWidth_),
    UName_(ptf.UName_)
{}


Foam::windkesselOutletPressureFvPatchScalarField::
windkesselOutletPressureFvPatchScalarField
(
    const windkesselOutletPressureFvPatchScalarField& wkpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(wkpsf, iF),
    Rp_(wkpsf.Rp_), Rd_(wkpsf.Rd_), C_(wkpsf.C_), L_(wkpsf.L_),
    Pd_(wkpsf.Pd_), rho_(wkpsf.rho_),
    Pooo_(wkpsf.Pooo_), Poo_(wkpsf.Poo_), Po_(wkpsf.Po_), Pn_(wkpsf.Pn_),
    Qooo_(wkpsf.Qooo_), Qoo_(wkpsf.Qoo_), Qo_(wkpsf.Qo_), Qn_(wkpsf.Qn_),
    windkesselModel_(wkpsf.windkesselModel_),
    diffScheme_(wkpsf.diffScheme_),
    timeIndex_(wkpsf.timeIndex_),
    lastDt_(wkpsf.lastDt_),
    relaxFactor_(wkpsf.relaxFactor_),
    useLaggedFlux_(wkpsf.useLaggedFlux_),
    flipFluxSign_(wkpsf.flipFluxSign_),
    backflowStabilization_(wkpsf.backflowStabilization_),
    penaltyFactor_(wkpsf.penaltyFactor_),
    smoothingWidth_(wkpsf.smoothingWidth_),
    UName_(wkpsf.UName_)
{}

// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

void Foam::windkesselOutletPressureFvPatchScalarField::updateCoeffs()
{
    if (updated()) return;

    const label  patchI = patch().index();
    const scalar dt     = db().time().deltaTValue();

    // Guard for variable dt when using BDF2
    const bool variableDt = (lastDt_ > SMALL && mag(dt - lastDt_) > SMALL);
    const diffSchemeType schemeToUse = (diffScheme_ == secondOrder && variableDt) ? firstOrder : diffScheme_;

    // --- Current outlet flux: always use corrected, conservative flux 'phi'
    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>("phi");

    scalar Qcur = gSum(phi.boundaryField()[patchI]);   // m^3/s
    if (phi.dimensions()[0] != 0)
    {
       Qcur /= max(rho_, SMALL);
    }
    // Optional sign flip (if your convention differs)
    if (flipFluxSign_) Qcur = -Qcur;

    // --- Patch-average kinematic pressure -> absolute (Pa) for ODE
    const fvPatchField<scalar>& pFld =
        db().lookupObject<volScalarField>("p").boundaryField()[patchI];

    const scalar area    = gSum(patch().magSf());
    const scalar PbarAbs = rho_ * (gSum(pFld * patch().magSf()) / (area + VSMALL));

    // --- Update histories once per (outer) time step
    if (db().time().timeIndex() != timeIndex_)
    {
        timeIndex_ = db().time().timeIndex();

        Pooo_ = Poo_;
        Poo_  = Po_;
        Po_   = Pn_;
        Pn_   = PbarAbs;

        Qooo_ = Qoo_;
        Qoo_  = Qo_;
        Qo_   = Qn_;
        Qn_   = Qcur;
    }

    // Choose which flux to use in the non-derivative terms
    const scalar Qin = useLaggedFlux_ ? Qn_ : Qcur;

    // --- Compute new absolute pressure (Pa)
    scalar PnAbs = PbarAbs; // default init

    if (db().time().timeIndex() > 1)
    {
        switch (windkesselModel_)
        {
            case Resistive:
            {
                PnAbs = Rd_*Qin + Pd_;
                break;
            }

            case WK2:
            {
                switch (diffScheme_)
                {
                    case firstOrder:
                    {
                        // Use lagged Qin in algebraic term; keep history for derivative
                        PnAbs = Rd_*Qin + Pd_ + Rd_*C_*(Po_/dt);
                        PnAbs /= (1.0 + Rd_*C_/dt) + SMALL;
                        break;
                    }
                    case secondOrder:
                    {
                        PnAbs = Rd_*Qin + Pd_
                              - Rd_*C_*((Poo_ - 4*Po_)/(2*dt));
                        PnAbs /= (1.0 + 3*Rd_*C_/(2*dt)) + SMALL;
                        break;
                    }
                }
                break;
            }

            case WK3:
            {
                switch (diffScheme_)
                {
                    case firstOrder:
                    {
                        PnAbs =
                              Rp_*Rd_*C_*((useLaggedFlux_ ? (Qn_ - Qo_) : (Qcur - Qo_))/dt)
                            + (Rp_ + Rd_)*Qin
                            + Pd_
                            + Rd_*C_*(Po_/dt);
                        PnAbs /= (1.0 + Rd_*C_/dt) + SMALL;
                        break;
                    }
                    case secondOrder:
                    {
                        PnAbs =
                              Rp_*Rd_*C_*((useLaggedFlux_ ? (3*Qn_ - 4*Qo_ + Qoo_) : (3*Qcur - 4*Qo_ + Qoo_))/(2*dt))
                            + (Rp_ + Rd_)*Qin
                            + Pd_
                            - Rd_*C_*((Poo_ - 4*Po_)/(2*dt));
                        PnAbs /= (1.0 + 3*Rd_*C_/(2*dt)) + SMALL;
                        break;
                    }
                }
                break;
            }

            case WK4Series:
            {
                switch (diffScheme_)
                {
                    case firstOrder:
                    {
                        PnAbs =
                              (Rp_ + Rd_)*Qin
                            + (L_ + Rp_*Rd_*C_)*((useLaggedFlux_ ? (Qn_ - Qo_) : (Qcur - Qo_))/dt)
                            + Rd_*C_*L_*((useLaggedFlux_ ? (Qn_ - 2*Qo_ + Qoo_) : (Qcur - 2*Qo_ + Qoo_))/(dt*dt))
                            + Pd_
                            + Rd_*C_*(Po_/dt);
                        PnAbs /= (1.0 + Rd_*C_/dt) + SMALL;
                        break;
                    }
                    case secondOrder:
                    {
                        PnAbs =
                              (Rp_ + Rd_)*Qin
                            + (L_ + Rp_*Rd_*C_)*((useLaggedFlux_ ? (3*Qn_ - 4*Qo_ + Qoo_) : (3*Qcur - 4*Qo_ + Qoo_))/(2*dt))
                            + Rd_*C_*L_*((useLaggedFlux_ ? (2*Qn_ - 5*Qo_ + 4*Qoo_ - Qooo_) : (2*Qcur - 5*Qo_ + 4*Qoo_ - Qooo_))/(dt*dt))
                            + Pd_
                            - Rd_*C_*((Poo_ - 4*Po_)/(2*dt));
                        PnAbs /= (1.0 + 3*Rd_*C_/(2*dt)) + SMALL;
                        break;
                    }
                }
                break;
            }

            case WK4Parallel:
            {
                switch (diffScheme_)
                {
                    case firstOrder:
                    {
                        PnAbs =
                              Rd_*Qin
                            + L_*(1 + Rd_/Rp_)*((useLaggedFlux_ ? (Qn_ - Qo_) : (Qcur - Qo_))/dt)
                            + Rd_*C_*L_*((useLaggedFlux_ ? (Qn_ - 2*Qo_ + Qoo_) : (Qcur - 2*Qo_ + Qoo_))/(dt*dt))
                            + Pd_
                            + ((L_ + Rp_*Rd_*C_)/Rp_)*(Po_/dt)
                            + (Rd_*C_*L_/Rp_)*((2*Po_ - Poo_)/(dt*dt));
                        PnAbs /= (1.0
                                 + (L_ + Rp_*Rd_*C_)/(Rp_*dt)
                                 + (Rd_*C_*L_)/(Rp_*dt*dt)) + SMALL;
                        break;
                    }
                    case secondOrder:
                    {
                        PnAbs =
                              Rd_*Qin
                            + L_*(1 + Rd_/Rp_)*((useLaggedFlux_ ? (3*Qn_ - 4*Qo_ + Qoo_) : (3*Qcur - 4*Qo_ + Qoo_))/(2*dt))
                            + Rd_*C_*L_*((useLaggedFlux_ ? (2*Qn_ - 5*Qo_ + 4*Qoo_ - Qooo_) : (2*Qcur - 5*Qo_ + 4*Qoo_ - Qooo_))/(dt*dt))
                            + Pd_
                            - ((L_ + Rp_*Rd_*C_)/Rp_)*((Poo_ - 4*Po_)/(2*dt))
                            - (Rd_*C_*L_/Rp_)*((-5*Po_ + 4*Poo_ - Pooo_)/(dt*dt));
                        PnAbs /= (1.0
                                 + 3*(L_ + Rp_*Rd_*C_)/(2*Rp_*dt)
                                 + 2*Rd_*C_*L_/(Rp_*dt*dt)) + SMALL;
                        break;
                    }
                }
                break;
            }

            default:
            {
                FatalErrorInFunction << "Unknown Windkessel Model!" << exit(FatalError);
            }
        }
    }

    // --- Under-relaxed write-back as kinematic pressure (m^2/s^2)
    const scalar pk = PnAbs / (rho_ + SMALL);

    // Base target: uniform Windkessel pressure on all faces
    scalarField target(patch().size(), pk);

    // --- Opt-in backflow penalty: on reversed faces add a kinematic dynamic
    //     pressure  mask*penaltyFactor*0.5*|U|^2  (no rho: field is kinematic).
    if (backflowStabilization_)
    {
        const fvPatchField<vector>& Up =
            db().lookupObject<volVectorField>(UName_).boundaryField()[patchI];

        const scalarField& phip = phi.boundaryField()[patchI];

        scalarField mask(patch().size(), Zero);
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

        forAll(target, i)
        {
            target[i] = pk + mask[i]*penaltyFactor_*0.5*magSqr(Up[i]);
        }
    }

    if (relaxFactor_ > 0 && relaxFactor_ < 1)
    {
        // mix with current boundary values (*this) for damping
        operator==( relaxFactor_*target + (1.0 - relaxFactor_)*(*this) );
    }
    else
    {
        operator==( target );
    }

    lastDt_ = dt;
    fixedValueFvPatchScalarField::updateCoeffs();
}

// * * * * * * * * * * * * * * * * Write * * * * * * * * * * * * * * * * * * //

void Foam::windkesselOutletPressureFvPatchScalarField::write(Ostream& os) const
{
    fixedValueFvPatchScalarField::write(os);
    os.writeKeyword("type") << this->type() << token::END_STATEMENT << nl;

    os.writeKeyword("Rd")  << Rd_  << token::END_STATEMENT << nl;
    os.writeKeyword("Pd")  << Pd_  << token::END_STATEMENT << nl;
    os.writeKeyword("rho") << rho_ << token::END_STATEMENT << nl;

    switch (windkesselModel_)
    {
        case Resistive:
            os.writeKeyword("windkesselModel") << "Resistive" << token::END_STATEMENT << nl;
            break;
        case WK2:
            os.writeKeyword("windkesselModel") << "WK2" << token::END_STATEMENT << nl;
            os.writeKeyword("C") << C_ << token::END_STATEMENT << nl;
            break;
        case WK3:
            os.writeKeyword("windkesselModel") << "WK3" << token::END_STATEMENT << nl;
            os.writeKeyword("Rp") << Rp_ << token::END_STATEMENT << nl;
            os.writeKeyword("C")  << C_  << token::END_STATEMENT << nl;
            break;
        case WK4Series:
            os.writeKeyword("windkesselModel") << "WK4Series" << token::END_STATEMENT << nl;
            os.writeKeyword("Rp") << Rp_ << token::END_STATEMENT << nl;
            os.writeKeyword("C")  << C_  << token::END_STATEMENT << nl;
            os.writeKeyword("L")  << L_  << token::END_STATEMENT << nl;
            break;
        case WK4Parallel:
            os.writeKeyword("windkesselModel") << "WK4Parallel" << token::END_STATEMENT << nl;
            os.writeKeyword("Rp") << Rp_ << token::END_STATEMENT << nl;
            os.writeKeyword("C")  << C_  << token::END_STATEMENT << nl;
            os.writeKeyword("L")  << L_  << token::END_STATEMENT << nl;
            break;
    }

    switch (diffScheme_)
    {
        case firstOrder:
            os.writeKeyword("diffScheme") << "firstOrder" << token::END_STATEMENT << nl;
            break;
        case secondOrder:
            os.writeKeyword("diffScheme") << "secondOrder" << token::END_STATEMENT << nl;
            break;
    }

    // New controls written for restart/reproducibility
    os.writeKeyword("relaxFactor")   << relaxFactor_   << token::END_STATEMENT << nl;
    os.writeKeyword("useLaggedFlux") << (useLaggedFlux_ ? "true" : "false") << token::END_STATEMENT << nl;
    os.writeKeyword("flipFluxSign")  << (flipFluxSign_  ? "true" : "false") << token::END_STATEMENT << nl;

    // NEW: backflow penalty controls
    os.writeKeyword("backflowStabilization") << (backflowStabilization_ ? "true" : "false") << token::END_STATEMENT << nl;
    os.writeKeyword("penaltyFactor")  << penaltyFactor_  << token::END_STATEMENT << nl;
    os.writeKeyword("smoothingWidth") << smoothingWidth_ << token::END_STATEMENT << nl;
    os.writeKeyword("U")              << UName_          << token::END_STATEMENT << nl;

    os.writeKeyword("Qooo") << Qooo_ << token::END_STATEMENT << nl;
    os.writeKeyword("Qoo")  << Qoo_  << token::END_STATEMENT << nl;
    os.writeKeyword("Qo")   << Qo_   << token::END_STATEMENT << nl;
    os.writeKeyword("Qn")   << Qn_   << token::END_STATEMENT << nl;

    os.writeKeyword("Pooo") << Pooo_ << token::END_STATEMENT << nl;
    os.writeKeyword("Poo")  << Poo_  << token::END_STATEMENT << nl;
    os.writeKeyword("Po")   << Po_   << token::END_STATEMENT << nl;
    os.writeKeyword("Pn")   << Pn_   << token::END_STATEMENT << nl;

    writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * Registration * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        windkesselOutletPressureFvPatchScalarField
    );
}

