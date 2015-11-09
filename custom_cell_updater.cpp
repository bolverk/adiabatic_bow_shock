#include "custom_cell_updater.hpp"

CustomCellUpdater::CustomCellUpdater(void) {}

vector<ComputationalCell> CustomCellUpdater::operator()
(const Tessellation& /*tess*/,
 const PhysicalGeometry& /*pg*/,
 const EquationOfState& eos,
 vector<Extensive>& extensives,
 const vector<ComputationalCell>& old,
 const CacheData& cd) const
{
  vector<ComputationalCell> res = old;
  for(size_t i=0;i<extensives.size();++i){
    if(old[i].stickers.find("obstacle")->second)
      continue;
    const double volume = cd.volumes[i];
    res[i].density = extensives[i].mass/volume;
    res[i].velocity = extensives[i].momentum / extensives[i].mass;
    const double total_energy = extensives[i].energy/extensives[i].mass;
    const double kinetic_energy =
      0.5*ScalarProd(res[i].velocity, res[i].velocity);
    const double energy = total_energy - kinetic_energy;
    res[i].pressure = eos.de2p(res[i].density, energy);
    for(boost::container::flat_map<std::string,double>::const_iterator it =
	  extensives[i].tracers.begin();
	it!=extensives[i].tracers.end();++it)
      res[i].tracers[it->first] = it->second/extensives[i].mass;
  }
  return res;
}
