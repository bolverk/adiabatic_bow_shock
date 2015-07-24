#include "temperature_appendix.hpp"

TemperatureAppendix::TemperatureAppendix
(double boltzmann_constant,
 double atomic_mass):
  boltzmann_constant_(boltzmann_constant),
  atomic_mass_(atomic_mass) {}

string TemperatureAppendix::getName(void) const
{
  return "temperature";
}

vector<double> TemperatureAppendix::operator()(const hdsim& sim) const
{
  const vector<ComputationalCell>& cells = sim.getAllCells();
  vector<double> res(cells.size(),1);
  for(size_t i=0;i<cells.size();++i)
    res[i] = atomic_mass_*cells[i].pressure/boltzmann_constant_/cells[i].density;
  return res;
}
