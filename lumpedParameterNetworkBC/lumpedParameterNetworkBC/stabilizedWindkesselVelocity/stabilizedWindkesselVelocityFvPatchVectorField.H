/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         |  stabilizedWindkesselVelocity
   \\    /   O peration     |
    \\  /    A nd           |  Directional backflow stabilization for WK outlets
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    Directional backflow-stabilized velocity boundary condition for Windkessel /
    lumped-parameter pressure outlets.

    Ported to OpenFOAM v2506 (ESI) from JieWangnk/OpenFOAM-WK (OpenFOAM v12).
    Derives from directionMixedFvPatchVectorField so the stabilization enters the
    implicit PIMPLE matrix. On flow reversal (phi < 0) the TANGENTIAL velocity
    component is damped towards zero (suppressing reentrant vortices) while the
    NORMAL component is left free (zero-gradient) so the Windkessel pressure can
    still drive the outflow. During outflow it reduces to a plain zero-gradient.

    valueFraction = mask(phi) * ( betaN*(n n) + betaT*(I - n n) )
    refValue      = 0
    refGrad       = 0

    where mask(phi) is a smooth tanh ramp (smoothingWidth > 0) or a hard
    Heaviside switch (smoothingWidth = 0) on the patch flux.

    License: GPLv3+
\*---------------------------------------------------------------------------*/

#ifndef stabilizedWindkesselVelocityFvPatchVectorField_H
#define stabilizedWindkesselVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "directionMixedFvPatchFields.H"
#include "Switch.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class stabilizedWindkesselVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class stabilizedWindkesselVelocityFvPatchVectorField
:
    public directionMixedFvPatchVectorField
{
    // Private data

        //- Flux field name
        word phiName_;

        //- Tangential backflow damping (0-1)
        scalar betaT_;

        //- Normal backflow damping (0-1); keep 0 to leave WK coupling free
        scalar betaN_;

        //- Master on/off switch
        Switch enableStabilization_;

        //- Overall multiplier on the damping
        scalar dampingFactor_;

        //- Smooth tanh ramp width (relative to mean |phi|); 0 = hard switch
        scalar smoothingWidth_;


public:

    //- Runtime type information
    TypeName("stabilizedWindkesselVelocity");


    // Constructors

        //- Construct from patch and internal field
        stabilizedWindkesselVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        stabilizedWindkesselVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        stabilizedWindkesselVelocityFvPatchVectorField
        (
            const stabilizedWindkesselVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        stabilizedWindkesselVelocityFvPatchVectorField
        (
            const stabilizedWindkesselVelocityFvPatchVectorField&
        );

        //- Construct as copy setting internal field reference
        stabilizedWindkesselVelocityFvPatchVectorField
        (
            const stabilizedWindkesselVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<vector>> clone() const
        {
            return fvPatchField<vector>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<vector>> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return fvPatchField<vector>::Clone(*this, iF);
        }


    // Member Functions

        //- True: this patch field is altered by assignment.
        virtual bool assignable() const { return true; }

        // Mapping functions

            virtual void autoMap(const fvPatchFieldMapper&);

            virtual void rmap
            (
                const fvPatchVectorField&,
                const labelList&
            );

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        virtual void operator=(const fvPatchField<vector>& pvf);
};

} // End namespace Foam

#endif

// ************************************************************************* //
