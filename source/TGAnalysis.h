#pragma once

// Parent class
#include "BaseSolver.h"

// Maps
#include <maps/ThermodynamicsMap_CHEMKIN.h>
#include <maps/KineticsMap_CHEMKIN.h>
#include <maps/TransportPropertiesMap_CHEMKIN.h>
#include <maps/ThermodynamicsMap_Solid_CHEMKIN.h>
#include <maps/KineticsMap_Solid_CHEMKIN.h>

#include <math/external-ode-solvers/ODE_Parameters.h>

namespace BioSMOKE
{
class TGAnalysis : public virtual BaseSolver
{
  public: // clang-format off
    TGAnalysis( OpenSMOKE::ThermodynamicsMap_CHEMKIN &thermodynamicsMap, 
                OpenSMOKE::KineticsMap_CHEMKIN &kineticsMap,
                OpenSMOKE::TransportPropertiesMap_CHEMKIN &transportMap,
                OpenSMOKE::ThermodynamicsMap_Solid_CHEMKIN &thermodynamicsSolidMap,
                OpenSMOKE::KineticsMap_Solid_CHEMKIN &kineticsSolidMap, 
                OpenSMOKE::ODE_Parameters &ode_parameters,
                //OpenSMOKE::BioSMOKE_Options &biosmoke_options, 
                const double T0, 
                const double P0,
                const double rho0_solid,
                const std::vector<double> &omega0_gas, 
                const std::vector<double> &omega0_solid, 
                const double heating_rate); // clang-format on

    virtual void Solve(const double t0, const double tf);

    virtual int Equations(const double t, const std::vector<double> &y, std::vector<double> &dy);

    virtual int Print(const double t, const std::vector<double> &y);

  protected:
    double V_solid_;      // Current volume of the solid phase
    double T0_;           // Initial temperature
    double T_;            // Current temperature
    double P_;            // Current pressure
    double P0_;           // Initial pressure
    double heating_rate_; // Heating rate [K/s]

    double MW_solid_;                 // Current solid molecular weight
    double rho_solid_;                // Current solid density
    std::vector<double> omega_solid_; // Current solid composition
    std::vector<double> x_solid_;     // Current solid mole fractions

    double MW_gas_;                 // Current gas molecular weight
    double rho_gas_;                // Current gas density
    std::vector<double> omega_gas_; // Current gas composition
    std::vector<double> x_gas_;     // Current gas mole fractions
};
} // namespace BioSMOKE

#include "TGAnalysis.hpp"
