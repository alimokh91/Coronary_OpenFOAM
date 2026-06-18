#ifndef windkesselOutletPressureFvPatchScalarField_H
#define windkesselOutletPressureFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "volFieldsFwd.H"

namespace Foam
{

class windkesselOutletPressureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
public:

    // Windkessel types
    enum windkesselModelType { Resistive, WK2, WK3, WK4Series, WK4Parallel };
    enum diffSchemeType      { firstOrder, secondOrder };

private:

    // Parameters (absolute/SI for ODE)
    scalar Rp_, Rd_, C_, L_, Pd_, rho_;

    // Histories (absolute pressure in Pa; flux in m^3/s)
    scalar Pooo_, Poo_, Po_, Pn_;
    scalar Qooo_, Qoo_, Qo_, Qn_;

    // Model / scheme
    windkesselModelType windkesselModel_;
    diffSchemeType      diffScheme_;
    label               timeIndex_;

    // Track last time step to guard BDF2 against variable dt
    scalar             lastDt_;

    // --- NEW: stability/robustness controls (dictionary selectable)
    scalar relaxFactor_;     // 0..1 under-relaxation on patch value (kinematic)
    bool   useLaggedFlux_;   // true: use Qn_ instead of current sub-iter flux in ODE
    bool   flipFluxSign_;    // true: multiply outlet flux by -1 before using it

    // --- NEW: opt-in backflow pressure penalty (JieWangnk-style, kinematic)
    //     On reversed faces adds  mask*penaltyFactor*0.5*|U|^2  to the WK pressure
    bool   backflowStabilization_; // master on/off (default false -> unchanged)
    scalar penaltyFactor_;         // scale on the 0.5|U|^2 penalty (default 1)
    scalar smoothingWidth_;        // tanh ramp width rel. to mean |phi|; 0 = hard
    word   UName_;                 // velocity field name (default "U")

public:

    // Constructors
    TypeName("windkesselOutletPressure");

    windkesselOutletPressureFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    windkesselOutletPressureFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    windkesselOutletPressureFvPatchScalarField
    (
        const windkesselOutletPressureFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    windkesselOutletPressureFvPatchScalarField
    (
        const windkesselOutletPressureFvPatchScalarField&,
        const DimensionedField<scalar, volMesh>&
    );

    // Member functions
    virtual void updateCoeffs();
    virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif

