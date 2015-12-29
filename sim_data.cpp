#include "sim_data.hpp"
#include "source/misc/vector_initialiser.hpp"

SimData::SimData(void):
  pg_(Vector2D(0,0), Vector2D(0,1)),
  outer_(Vector2D(0,-0.5),
	 Vector2D(2.0,4.5)),
  init_points_(cartesian_mesh(200,500,
			      outer_.getBoundary().first,
			      outer_.getBoundary().second)),
  tess_(init_points_, outer_),
  eos_(5./3.),
  point_motion_(),
  evc_(),
  rs_(),
  geom_source_(pg_.getAxis()),
  wind_source_
  (4*M_PI/1e4/(4.0*M_PI*pow(0.01,3)/3.),
   10,
   1e-3,
   0.01),
  source_
  (VectorInitialiser<SourceTerm*>
   (&geom_source_)
   (&wind_source_)
   ()),
  tsf_(0.3),
  fc_(rs_),
  eu_(),
  cu_(),
  sim_(tess_,
       outer_,
       pg_,
       calc_init_cond(tess_),
       eos_,
       point_motion_,
       evc_,
       source_,
       tsf_,
       fc_,
       eu_,
       cu_) {}

hdsim& SimData::getSim(void)
{
  return sim_;
}
