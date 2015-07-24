#ifndef TEMPERATURE_APPENDIX_HPP
#define TEMPERATURE_APPENDIX_HPP 1

#include "source/newtonian/two_dimensional/hdf5_diagnostics.hpp"

class TemperatureAppendix: public DiagnosticAppendix
{
public:

  TemperatureAppendix
  (double boltzmann_constant,
   double atomic_mass);

  string getName(void) const;

  vector<double> operator()(const hdsim& sim) const;

private:
  const double boltzmann_constant_;
  const double atomic_mass_;
};

#endif // TEMPERATURE_APPENDIX_HPP
