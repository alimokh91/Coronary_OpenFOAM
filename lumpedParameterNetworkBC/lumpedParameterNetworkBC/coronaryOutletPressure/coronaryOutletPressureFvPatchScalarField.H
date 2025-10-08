#ifndef coronaryOutletPressureFvPatchScalarField_H
#define coronaryOutletPressureFvPatchScalarField_H

#include <vector>
#include <utility>

#include "fixedValueFvPatchFields.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "fileName.H"
#include "word.H"
#include "fvPatchFieldMapper.H"
#include "dictionary.H"

namespace Foam
{

/**
 * \brief Coronary outlet pressure boundary condition with intramyocardial forcing.
 *
 * Two-state coronary LPN:
 *   - Resistances: Ra, Ram, Rv, (optional Rvm)
 *   - Compliances: Ca, Cim
 *   - Intramyocardial pressure Pim(t) from file (with period & scaling)
 *
 * Time integration: BDF1 / BDF2 (BDF2 used only when dt is constant).
 * Modes for Pim coupling:
 *   - "off"          : ignore Pim except via Rvm if provided
 *   - "backPressure" : Pim acts as back-pressure on Rv
 *   - "transmural"   : backPressure + transmural capacitors (Cim; Ca optional)
 *
 * The boundary writes **kinematic** pressure (m^2/s^2), converting absolute
 * pressure by division with rho.
 */
class coronaryOutletPressureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
public:

    // Time integration scheme
    enum diffScheme { BDF1, BDF2 };

private:

    // ---- Physical parameters ----
    scalar Ra_{0.0};
    scalar Ram_{0.0};
    scalar Rv_{0.0};
    scalar Rvm_{0.0};
    scalar Ca_{0.0};
    scalar Cim_{0.0};

    // ---- Fluid density & diagonal stabilization ----
    scalar rho_{1060.0};
    scalar regularize_{0.0};

    // ---- Numerics / controls ----
    diffScheme scheme_{BDF1};
    bool   useLaggedFlux_{true};
    scalar relaxFactor_{0.6};
    bool   flipFluxSign_{false};

    // ---- Intramyocardial pressure (Pim) controls ----
    fileName PimFile_{"constant/Pim.dat"};
    scalar   PimPeriod_{0.0};
    scalar   PimScaling_{1.0};
    word     pimMode_{"off"};        // "off" | "backPressure" | "transmural"
    bool     transmuralOnCa_{false}; // also add +Ca*dPim/dt if true

    // ---- Pim lookup table ----
    std::vector<std::pair<scalar, scalar>> PimTable_;  // (t, Pim[Pa])

    // ---- State variables (absolute Pa) ----
    scalar Pn_{0.0},   Pnm1_{0.0};   // arterial node
    scalar Pio_{0.0},  Pio_m1_{0.0}; // intramyocardial node
    scalar Qn_{0.0},   Qm1_{0.0};    // patch flux history (m^3/s)
    scalar dtPrev_{0.0};
    bool   firstIter_{true};

    // ---- Warm start ----
    bool   haveWarmStart_{false};
    scalar Po0_{0.0}, Pio0_{0.0}, Qo0_{0.0};

    // ---- For dPim/dt in BDF1 ----
    bool   pimHasPrev_{false};
    scalar pimPrev_{0.0};

    // ---- Flux sign self-calibration (one-time, per patch) ----
    scalar qSign_{1.0};
    bool   signLocked_{false};

    // ---- Internal helpers ----
    void   readPimTable(const fileName& path);
    scalar interpolatePim(const scalar t) const;

public:

    // Runtime type information
    TypeName("coronaryOutletPressure");

    // Constructors
    coronaryOutletPressureFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    coronaryOutletPressureFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    coronaryOutletPressureFvPatchScalarField
    (
        const coronaryOutletPressureFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    coronaryOutletPressureFvPatchScalarField
    (
        const coronaryOutletPressureFvPatchScalarField&,
        const DimensionedField<scalar, volMesh>&
    );

    // Core methods
    virtual void updateCoeffs();
    virtual void write(Ostream&) const;
};

} // namespace Foam

#endif

