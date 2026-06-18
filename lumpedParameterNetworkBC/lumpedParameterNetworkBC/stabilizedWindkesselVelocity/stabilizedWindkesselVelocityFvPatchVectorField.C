/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         |  stabilizedWindkesselVelocity
   \\    /   O peration     |
    \\  /    A nd           |  Directional backflow stabilization for WK outlets
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Ported to OpenFOAM v2506 from JieWangnk/OpenFOAM-WK (v12). License: GPLv3+
\*---------------------------------------------------------------------------*/

#include "stabilizedWindkesselVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::stabilizedWindkesselVelocityFvPatchVectorField::
stabilizedWindkesselVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(p, iF),
    phiName_("phi"),
    betaT_(0.2),
    betaN_(0.0),
    enableStabilization_(true),
    dampingFactor_(1.0),
    smoothingWidth_(0.1)
{
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = Zero;
}


Foam::stabilizedWindkesselVelocityFvPatchVectorField::
stabilizedWindkesselVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    directionMixedFvPatchVectorField(p, iF),
    phiName_(dict.getOrDefault<word>("phi", "phi")),
    betaT_(dict.getOrDefault<scalar>("betaT", 0.2)),
    betaN_(dict.getOrDefault<scalar>("betaN", 0.0)),
    enableStabilization_(dict.getOrDefault<Switch>("enableStabilization", true)),
    dampingFactor_(dict.getOrDefault<scalar>("dampingFactor", 1.0)),
    smoothingWidth_(dict.getOrDefault<scalar>("smoothingWidth", 0.1))
{
    fvPatchFieldBase::readDict(dict);

    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = Zero;

    if (!this->readValueEntry(dict, IOobjectOption::READ_IF_PRESENT))
    {
        // No "value" entry: initialise from the internal cell values
        fvPatchVectorField::operator=(this->patchInternalField());
    }
}


Foam::stabilizedWindkesselVelocityFvPatchVectorField::
stabilizedWindkesselVelocityFvPatchVectorField
(
    const stabilizedWindkesselVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    directionMixedFvPatchVectorField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    betaT_(ptf.betaT_),
    betaN_(ptf.betaN_),
    enableStabilization_(ptf.enableStabilization_),
    dampingFactor_(ptf.dampingFactor_),
    smoothingWidth_(ptf.smoothingWidth_)
{}


Foam::stabilizedWindkesselVelocityFvPatchVectorField::
stabilizedWindkesselVelocityFvPatchVectorField
(
    const stabilizedWindkesselVelocityFvPatchVectorField& ptf
)
:
    directionMixedFvPatchVectorField(ptf),
    phiName_(ptf.phiName_),
    betaT_(ptf.betaT_),
    betaN_(ptf.betaN_),
    enableStabilization_(ptf.enableStabilization_),
    dampingFactor_(ptf.dampingFactor_),
    smoothingWidth_(ptf.smoothingWidth_)
{}


Foam::stabilizedWindkesselVelocityFvPatchVectorField::
stabilizedWindkesselVelocityFvPatchVectorField
(
    const stabilizedWindkesselVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(ptf, iF),
    phiName_(ptf.phiName_),
    betaT_(ptf.betaT_),
    betaN_(ptf.betaN_),
    enableStabilization_(ptf.enableStabilization_),
    dampingFactor_(ptf.dampingFactor_),
    smoothingWidth_(ptf.smoothingWidth_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::stabilizedWindkesselVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    directionMixedFvPatchVectorField::autoMap(m);
}


void Foam::stabilizedWindkesselVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    directionMixedFvPatchVectorField::rmap(ptf, addr);
}


void Foam::stabilizedWindkesselVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (!enableStabilization_)
    {
        // Pure zero-gradient (no damping at all)
        refValue() = Zero;
        refGrad() = Zero;
        valueFraction() = Zero;

        directionMixedFvPatchVectorField::updateCoeffs();
        directionMixedFvPatchVectorField::evaluate();
        return;
    }

    // Patch flux (v2506 single-template-arg lookup)
    const scalarField& phip =
        patch().lookupPatchField<surfaceScalarField>(phiName_);

    // Backflow mask in [0,1]: 1 where flow re-enters (phi < 0)
    scalarField mask(patch().size(), Zero);
    if (smoothingWidth_ > SMALL)
    {
        const label nF = returnReduce(phip.size(), sumOp<label>());
        const scalar phiRef =
            gSum(mag(phip))/max(scalar(nF), scalar(1));
        const scalar sw = max(smoothingWidth_*phiRef, SMALL);

        forAll(mask, i)
        {
            mask[i] = 0.5*(1.0 - Foam::tanh(phip[i]/sw));
        }
    }
    else
    {
        forAll(mask, i)
        {
            mask[i] = pos0(-phip[i] - SMALL);
        }
    }

    // Directional projection tensors
    const vectorField n(patch().nf());
    const symmTensorField normalProj(sqr(n));            // n (x) n
    const symmTensorField tangProj(symmTensor::I - normalProj);

    const scalar effBetaT = betaT_*dampingFactor_;
    const scalar effBetaN = betaN_*dampingFactor_;

    refValue() = Zero;   // damped components pulled towards 0
    refGrad() = Zero;    // free components: zero-gradient
    valueFraction() = mask*(effBetaN*normalProj + effBetaT*tangProj);

    directionMixedFvPatchVectorField::updateCoeffs();
    directionMixedFvPatchVectorField::evaluate();
}


void Foam::stabilizedWindkesselVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchField<vector>::write(os);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    os.writeEntry("betaT", betaT_);
    os.writeEntry("betaN", betaN_);
    os.writeEntry("enableStabilization", enableStabilization_);
    os.writeEntry("dampingFactor", dampingFactor_);
    os.writeEntry("smoothingWidth", smoothingWidth_);
    fvPatchField<vector>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::stabilizedWindkesselVelocityFvPatchVectorField::operator=
(
    const fvPatchField<vector>& pvf
)
{
    tmp<vectorField> normalValue = transform(valueFraction(), refValue());
    tmp<vectorField> transformGradValue = transform(I - valueFraction(), pvf);
    fvPatchField<vector>::operator=(normalValue + transformGradValue);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        stabilizedWindkesselVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
