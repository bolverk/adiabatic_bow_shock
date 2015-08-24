#ifndef CUSTOM_FLUX_CALCULATOR_HPP
#define CUSTOM_FLUX_CALCULATOR_HPP 1

#include "source/newtonian/two_dimensional/simple_flux_calculator.hpp"

class CustomFluxCalculator: public FluxCalculator
{
public:

  CustomFluxCalculator(const RiemannSolver& rs);

  vector<Extensive> operator()
  (const Tessellation& tess,
   const vector<Vector2D>& point_velocities,
   const vector<ComputationalCell>& cells,
   const vector<Extensive>& /*extensive*/,
   const CacheData& cd,
   const EquationOfState& eos,
   const double /*time*/,
   const double /*dt*/) const;

private:
  const RiemannSolver& rs_;

  const Conserved calcHydroFlux
  (const Tessellation& tess,
   const vector<Vector2D>& point_velocities,
   const vector<ComputationalCell>& cells,
   const EquationOfState& eos,
   const size_t i) const;
};

#endif // CUSTOM_FLUX_CALCULATOR_HPP
