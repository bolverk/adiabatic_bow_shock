#include "calc_init_cond.hpp"
#include "source/tessellation/shape_2d.hpp"

vector<ComputationalCell> calc_init_cond(const Tessellation& tess)
{
  vector<ComputationalCell> res(static_cast<size_t>(tess.GetPointNo()));
  const Circle obstacle(Vector2D(0,0),0.01);
  for(size_t i=0;i<res.size();++i){
    const Vector2D r = tess.GetMeshPoint(static_cast<int>(i));
    res[i].density = 1;
    res[i].pressure = 1e-9;
    res[i].velocity = Vector2D(0,1);
    res[i].stickers["obstacle"] = false;
    if(obstacle(r)){
      res[i].velocity = Vector2D(0,0);
      res[i].stickers["obstacle"] = true;
    }
  }
  return res;
}
