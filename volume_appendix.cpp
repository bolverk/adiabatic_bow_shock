#include "volume_appendix.hpp"

VolumeAppendix::VolumeAppendix(void) {}

string VolumeAppendix::getName(void) const
{
  return "volumes";
}

vector<double> VolumeAppendix::operator()(const hdsim& sim) const
{
  const CacheData& cd = sim.getCacheData();
  vector<double> res(cd.volumes.size(),0);
  for(size_t i=0;i<cd.volumes.size();++i)
    res[i] = cd.volumes[i];
  return res;
}
