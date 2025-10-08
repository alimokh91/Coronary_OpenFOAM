/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         |  stabilizedWindkesselVelocity
   \\    /   O peration     |
    \\  /    A nd           |  Smooth backflow damping for WK outlets
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Author: Ali + ChatGPT
    License: GPLv3+
\*---------------------------------------------------------------------------*/

#ifndef stabilizedWindkesselVelocityFvPatchVectorField_H
#define stabilizedWindkesselVelocityFvPatchVectorField_H

#include "zeroGradientFvPatchFields.H"
#include "vectorField.H"
#include "Switch.H"

namespace Foam
{

class stabilizedWindkesselVelocityFvPatchVectorField
:
    public zeroGradientFvPatchVectorField
{
    // === user parameters (read from dictionary) ===
    scalar kappa_;          // nonlinearity strength (2..10 recommended)
    scalar Uref_;           // reference speed [m/s]; if <= 0, estimate per step
    Switch  linearDamping_; // if true: use smoothstep linear scaling
    scalar  alphaMax_;      // max linear damping (<=0.9), used only if linearDamping_=true

    // --- state-aware (optional) ---
    Switch  stateAware_;    // enable state-aware scaling
    scalar  CoRef_;         // reference Courant (~0.5)
    scalar  wkTau_;         // Windkessel time-scale [s] (e.g. Rd*C); 0 => disabled
    scalar  tauRef_;        // reference time-scale for scaling [s], e.g. 0.2

public:

    //- Runtime type information
    TypeName("stabilizedWindkesselVelocity");

    // Constructors
    stabilizedWindkesselVelocityFvPatchVectorField
    (
        const fvPatch&,
        const DimensionedField<vector, volMesh>&
    );

    stabilizedWindkesselVelocityFvPatchVectorField
    (
        const stabilizedWindkesselVelocityFvPatchVectorField&,
        const fvPatch&,
        const DimensionedField<vector, volMesh>&,
        const fvPatchFieldMapper&
    );

    // Clone with new internal field (no mapper needed)
    stabilizedWindkesselVelocityFvPatchVectorField
    (
        const stabilizedWindkesselVelocityFvPatchVectorField&,
        const DimensionedField<vector, volMesh>&
    );

    stabilizedWindkesselVelocityFvPatchVectorField
    (
        const fvPatch&,
        const DimensionedField<vector, volMesh>&,
        const dictionary&
    );

    //- Clone
    virtual tmp<fvPatchVectorField> clone() const
    {
        return tmp<fvPatchVectorField>
        (
            new stabilizedWindkesselVelocityFvPatchVectorField(*this)
        );
    }

    //- Clone with new internal field reference
    virtual tmp<fvPatchVectorField> clone
    (
        const DimensionedField<vector, volMesh>& iF
    ) const
    {
        return tmp<fvPatchVectorField>
        (
            new stabilizedWindkesselVelocityFvPatchVectorField(*this, iF)
        );
    }

    //- Mapping
    virtual void autoMap(const fvPatchFieldMapper& m)
    {
        zeroGradientFvPatchVectorField::autoMap(m);
    }
    virtual void rmap(const fvPatchVectorField& ptf, const labelList& addr)
    {
        zeroGradientFvPatchVectorField::rmap(ptf, addr);
    }

    //- Read dictionary
    virtual void read(const dictionary&);

    //- Evaluate (apply damping to the patch velocity)
    virtual void evaluate(const Pstream::commsTypes commsType=Pstream::commsTypes::blocking);

    //- Write
    virtual void write(Ostream&) const;
};

} // namespace Foam

#endif

