#include "wind.hpp"

Wind::Wind(const double specific_mass_loss,
	   const double speed,
	   const double t2k,
	   const double radius):
  specific_mass_loss_(specific_mass_loss),
  speed_(speed),
  t2k_(t2k),
  radius_(radius) {}

vector<Extensive> Wind::operator()
  (const Tessellation& tess,
   const PhysicalGeometry& /*pg*/,
   const CacheData& cd,
   const vector<ComputationalCell>& /*cells*/,
   const vector<Extensive>& /*fluxes*/,
   const vector<Vector2D>& /*point_velocities*/,
   const double /*time*/) const
{
  vector<Extensive> res(static_cast<size_t>(tess.GetPointNo()));
  for(size_t i=0;i<res.size();++i){
    res[i].mass = 0;
    res[i].momentum = Vector2D(0,0);
    res[i].energy = 0;
    const Vector2D r = tess.GetMeshPoint(static_cast<int>(i));
    if(abs(r)<radius_){
      res[i].mass = specific_mass_loss_*cd.volumes[i];
      res[i].momentum = res[i].mass*speed_*r/abs(r);
      res[i].energy = 
	0.5*ScalarProd(res[i].momentum,res[i].momentum)/res[i].mass*(1+t2k_);
    }
  }
  return res;
}
