#pragma once

#include <maps/Maps_CHEMKIN>
#include <maps/ThermodynamicsMap_Solid_CHEMKIN.h>

namespace BioSMOKE
{

enum Analysis_Type
{
    THERMOGRAVIMETRIC_ANALYSIS,
    ONE_DIMENSIONAL_SPHERICAL_PARTICLE,
};

//! Base class for the BioSMOKE solver
/*!
    The purpose of this class is to provide a common interface for the BioSMOKE solver
*/

class BaseSolver
{

  public:
    /**
     * @brief Default constructor
     * @param thermodynamicsMap map containing the thermodynamic data for the gas phase
     * @param kineticsMap map containing the kinetic mechanism for the gas phase
     * @param transportMap map containing the transport properties for the gas phase
     * @param thermodynamicsSolidMap map containing the thermodynamic data for the solid phase
     * @param kineticsSolidMap map containing the kinetic mechanism for the solid phase
     * @param ode_parameters parameters governing the solution of the stiff ODE system
     * @param biosmoke_options options governing the output
     */

    BaseSolver(OpenSMOKE::ThermodynamicsMap_CHEMKIN &thermodynamicsMap, OpenSMOKE::KineticsMap_CHEMKIN &kineticsMap,
               OpenSMOKE::TransportPropertiesMap_CHEMKIN &transportMap,
               OpenSMOKE::ThermodynamicsMap_Solid_CHEMKIN &thermodynamicsSolidMap,
               OpenSMOKE::KineticsMap_Solid_CHEMKIN &kineticsSolidMap, OpenSMOKE::ODE_Parameters &ode_parameters,
               // OpenSMOKE::BioSMOKE_Options &biosmoke_options,
    );

    virtual ~BaseSolver() = 0;

    /**
     * @brief Solves the BioSMOKE solver
     * @param tf the final time of integration [s]
     */

    virtual void Solve(const double t0, const double tf) = 0;

    /**
     * @brief Ordinary differential equations corresponding to the symulation type
     * @param t current time [s]
     * @param y current solution
     * @param dy current unsteady terms
     */
    virtual int Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble &y,
                          OpenSMOKE::OpenSMOKEVectorDouble &dy) = 0;

    /**
     * @brief Writes the output (called at the end of each time step)
     * @param t current time [s]
     * @param y current solution
     */
    virtual int Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble &y) = 0;

    /**
     * @brief Returns the total number of equations
     * @return the total number of equations
     */
    unsigned int NumberOfEquations() const { return NE_; };

  protected:
    OpenSMOKE::ThermodynamicsMap_CHEMKIN &thermodynamicsMap_;
    OpenSMOKE::KineticsMap_CHEMKIN &kineticsMap_;
    OpenSMOKE::TransportPropertiesMap_CHEMKIN &transportMap_;
    OpenSMOKE::ThermodynamicsMap_Solid_CHEMKIN &thermodynamicsSolidMap_;
    OpenSMOKE::KineticsMap_Solid_CHEMKIN &kineticsSolidMap_;
    OpenSMOKE::ODE_Parameters &ode_parameters_;
    // OpenSMOKE::BioSMOKE_Options &biosmoke_options_;

  protected:
    Analysis_Type analysis_type_; // Type of analysis

    double T0_solid_;                  // Initial temperature of the solid phase [K]
    double P0_solid_;                  // Initial pressure of the solid phase [Pa]
    double rho0_solid_;                // Initial density of the solid phase [kg/m3]
    double MW0_solid_;                 // Initial molecular weight of the solid phase [kg/kmol]
    double V0_solid_;                  // Initial volume of the solid phase [m3]
    std::vector<double> omega0_solid_; // Initial composition of the solid phase [mass fractions]
    std::vector<double> x0_sold_;      // Initial composition of the solid phase [mole fractions]

    double T0_gas_;                  // Initial temperature of the gas phase [K]
    double P0_gas_;                  // Initial pressure of the gas phase [Pa]
    double rho0_gas_;                // Initial density of the gas phase [kg/m3]
    double MW0_gas_;                 // Initial molecular weight of the gas phase [kg/kmol]
    std::vector<double> omega0_gas_; // Initial composition of the gas phase [mass fractions]
    std::vector<double> x0_gas_;     // Initial composition of the gas phase [mole fractions]

    double final_time_; // Final time of the simulation [s]

    unsigned int NGS_; // Number of species in the gas phase [-]
    unsigned int NSS_; // Number of species in the solid phase [-]
    unsigned int NC_;  // Number of species in the system [-]
    unsigned int NE_;  // Number of equations [-]

    std::vector<double> y0_; // vector cibntaining the initial values for all the variables
    std::vector<double> yf_; // vector containing the final values for all the variables
};
} // namespace BioSMOKE

#include "BaseSolver.hpp"