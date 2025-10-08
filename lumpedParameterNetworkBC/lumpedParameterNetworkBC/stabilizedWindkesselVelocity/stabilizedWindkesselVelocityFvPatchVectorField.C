/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         |  stabilizedWindkesselVelocity
   \\    /   O peration     |
    \\  /    A nd           |  Smooth backflow damping for WK outlets
     \\/     M anipulation  |
-------------------------------------------------------------------------------
\*---------------------------------------------------------------------------*/

#include "stabilizedWindkesselVelocityFvPatchVectorField.H"
#include "fvCFD.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * *  Constructors  * * * * * * * * * * * * * * * * * *

stabilizedWindkesselVelocityFvPatchVectorField::
stabilizedWindkesselVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    zeroGradientFvPatchVectorField(p, iF),
    kappa_(5.0),
    Uref_(0.0),
    linearDamping_(false),
    alphaMax_(0.8),
    stateAware_(true),
    CoRef_(0.5),
    wkTau_(0.0),
    tauRef_(0.2)
{}


stabilizedWindkesselVelocityFvPatchVectorField::
stabilizedWindkesselVelocityFvPatchVectorField
(
    const stabilizedWindkesselVelocityFvPatchVectorField& rhs,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& m
)
:
    zeroGradientFvPatchVectorField(rhs, p, iF, m),
    kappa_(rhs.kappa_),
    Uref_(rhs.Uref_),
    linearDamping_(rhs.linearDamping_),
    alphaMax_(rhs.alphaMax_),
    stateAware_(rhs.stateAware_),
    CoRef_(rhs.CoRef_),
    wkTau_(rhs.wkTau_),
    tauRef_(rhs.tauRef_)
{}

// Clone with new internal field (no mapper)
stabilizedWindkesselVelocityFvPatchVectorField::
stabilizedWindkesselVelocityFvPatchVectorField
(
    const stabilizedWindkesselVelocityFvPatchVectorField& rhs,
    const DimensionedField<vector, volMesh>& iF
)
:
    zeroGradientFvPatchVectorField(rhs, iF),
    kappa_(rhs.kappa_),
    Uref_(rhs.Uref_),
    linearDamping_(rhs.linearDamping_),
    alphaMax_(rhs.alphaMax_),
    stateAware_(rhs.stateAware_),
    CoRef_(rhs.CoRef_),
    wkTau_(rhs.wkTau_),
    tauRef_(rhs.tauRef_)
{}


stabilizedWindkesselVelocityFvPatchVectorField::
stabilizedWindkesselVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    zeroGradientFvPatchVectorField(p, iF, dict),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 5.0)),
    Uref_(dict.lookupOrDefault<scalar>("Uref", 0.0)),
    linearDamping_(dict.lookupOrDefault<Switch>("linearDamping", false)),
    alphaMax_(dict.lookupOrDefault<scalar>("alphaMax", 0.8)),
    stateAware_(dict.lookupOrDefault<Switch>("stateAware", true)),
    CoRef_(dict.lookupOrDefault<scalar>("CoRef", 0.5)),
    wkTau_(dict.lookupOrDefault<scalar>("wkTau", 0.0)),
    tauRef_(dict.lookupOrDefault<scalar>("tauRef", 0.2))
{}


// * * * * * * * * * * * * * *  Read/Write  * * * * * * * * * * * * * * * * *

void stabilizedWindkesselVelocityFvPatchVectorField::read(const dictionary& dict)
{
    // NOTE: do NOT call any base read(dict) here on v2306.

    if (dict.found("kappa"))          { dict.lookup("kappa")          >> kappa_; }
    if (dict.found("Uref"))           { dict.lookup("Uref")           >> Uref_; }
    if (dict.found("linearDamping"))  { dict.lookup("linearDamping")  >> linearDamping_; }
    if (dict.found("alphaMax"))       { dict.lookup("alphaMax")       >> alphaMax_; }
    if (dict.found("stateAware"))     { dict.lookup("stateAware")     >> stateAware_; }
    if (dict.found("CoRef"))          { dict.lookup("CoRef")          >> CoRef_; }
    if (dict.found("wkTau"))          { dict.lookup("wkTau")          >> wkTau_; }
    if (dict.found("tauRef"))         { dict.lookup("tauRef")         >> tauRef_; }
}



void stabilizedWindkesselVelocityFvPatchVectorField::write(Ostream& os) const
{
    zeroGradientFvPatchVectorField::write(os);
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("Uref")  << Uref_  << token::END_STATEMENT << nl;
    os.writeKeyword("linearDamping") << linearDamping_ << token::END_STATEMENT << nl;
    os.writeKeyword("alphaMax") << alphaMax_ << token::END_STATEMENT << nl;
    os.writeKeyword("stateAware") << stateAware_ << token::END_STATEMENT << nl;
    os.writeKeyword("CoRef") << CoRef_ << token::END_STATEMENT << nl;
    os.writeKeyword("wkTau") << wkTau_ << token::END_STATEMENT << nl;
    os.writeKeyword("tauRef") << tauRef_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * *  Evaluate  * * * * * * * * * * * * * * * * * *

// Helper: smoothstep on [a,b]
static inline scalar smoothstep(const scalar a, const scalar b, const scalar x)
{
    if (x <= a) return 0.0;
    if (x >= b) return 1.0;
    const scalar t = (x - a)/(b - a + SMALL);
    return t*t*(3.0 - 2.0*t);
}

void stabilizedWindkesselVelocityFvPatchVectorField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    // zeroGradient behaviour first
    zeroGradientFvPatchVectorField::evaluate(commsType);

    const fvPatch& p = patch();

    // Face geometry
    const vectorField& Sf  = p.Sf();
    const scalarField  magSf = mag(Sf);
    vectorField nHat = Sf/(magSf + VSMALL); // outward normals

    // Current patch velocity from the owner cell (zeroGradient base)
    vectorField v = this->patchInternalField();

    // Estimate Uref if not provided: area-weighted mean *outflow* speed
    scalar UrefEff = Uref_;
    if (UrefEff <= SMALL)
    {
        scalar num = 0.0, den = 0.0;
        forAll(v, i)
        {
            const scalar vn = (v[i] & nHat[i]);
            if (vn > 0) { num += vn*magSf[i]; den += magSf[i]; }
        }
        UrefEff = (den > VSMALL ? max(num/den, 1e-3) : 0.05); // floor to avoid /0
    }

    // Optional: state-aware scaling
    scalar kappaScale = 1.0;
    if (stateAware_)
    {
        // Local Courant per face (approx): |vn| * dt * deltaCoeff
        const scalar dt = this->db().time().deltaTValue();
        const scalarField deltaCoeff = p.deltaCoeffs(); // ~ 1/delta
        scalar CoAvg = 0.0; scalar A = 0.0;

        forAll(v, i)
        {
            const scalar vn = mag(v[i] & nHat[i]);
            const scalar Co = vn * dt * deltaCoeff[i];
            CoAvg += Co * magSf[i];
            A += magSf[i];
        }
        CoAvg = (A > VSMALL ? CoAvg/A : 0.0);

        const scalar CoScale = min(1.0, (CoRef_ > SMALL ? CoAvg/CoRef_ : 0.0)); // [0,1]
        // Stronger damping as Co grows: (1 + CoScale) in [1,2]
        kappaScale *= (1.0 + CoScale);

        // Reduce damping for large WK time-scale (more compliant/longer reversal)
        if (wkTau_ > SMALL)
        {
            const scalar wkScale = 1.0/(1.0 + wkTau_/(tauRef_ + SMALL)); // (0,1]
            kappaScale *= max(0.1, wkScale); // keep some damping
        }
    }

    const scalar kappaEff = max(1e-6, kappa_ * kappaScale);

    // Apply smooth normal-only damping for backflow (vn < 0)
    forAll(v, i)
    {
        const vector vi = v[i];
        const vector nh = nHat[i];

        const scalar vn = (vi & nh);
        const vector vt = vi - vn*nh;

        scalar vnNew = vn;

        if (vn < 0.0) // only damp inflow at an outlet
        {
            if (linearDamping_)
            {
                // alpha rises smoothly from 0 to alphaMax_ as -vn goes 0 -> UrefEff
                const scalar a = alphaMax_ * smoothstep(0.0, UrefEff, -vn);
                vnNew = vn * (1.0 - a);
            }
            else
            {
                // Nonlinear (hyperbolic) damping: vn / (1 + kappaEff*|vn|/UrefEff)
                const scalar denom = 1.0 + kappaEff * (-vn) / (UrefEff + SMALL);
                vnNew = vn / denom;
            }
        }

        // Recompose velocity (preserve tangential)
        v[i] = vt + vnNew*nh;
    }

    // Assign the damped values back to the patch field
    this->operator==(v);
}

} // namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "addToRunTimeSelectionTable.H"
namespace Foam
{
    makePatchTypeField(fvPatchVectorField, stabilizedWindkesselVelocityFvPatchVectorField);
}

