namespace BioSMOKE
{

BaseSolver::BaseSolver(OpenSMOKE::ThermodynamicsMap_CHEMKIN &thermodynamicsMap,
                       OpenSMOKE::KineticsMap_CHEMKIN &kineticsMap,
                       OpenSMOKE::TransportPropertiesMap_CHEMKIN &transportMap,
                       OpenSMOKE::ThermodynamicsMap_Solid_CHEMKIN &thermodynamicsSolidMap,
                       OpenSMOKE::KineticsMap_Solid_CHEMKIN &kineticsSolidMap,
                       OpenSMOKE::ODE_Parameters &ode_parameters)
    : thermodynamicsMap_(thermodynamicsMap), kineticsMap_(kineticsMap), transportMap_(transportMap),
      thermodynamicsSolidMap_(thermodynamicsSolidMap), kineticsSolidMap_(kineticsSolidMap),
      ode_parameters_(ode_parameters)
{
}

BaseSolver::~BaseSolver() {}

} // namespace BioSMOKE