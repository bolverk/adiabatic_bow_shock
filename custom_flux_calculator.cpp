#include "custom_flux_calculator.hpp"
#include "source/newtonian/common/hydrodynamics.hpp"

namespace {
    double calc_tracer_flux(const Edge& edge,
			  const Tessellation& tess,
			  const vector<ComputationalCell>& cells,
			  const string& name,
			  const Conserved& hf)
  {
    if(hf.Mass>0 &&
       edge.neighbors.first>0 &&
       edge.neighbors.first<tess.GetPointNo())
      return hf.Mass*
	cells[static_cast<size_t>(edge.neighbors.first)].tracers.find(name)->second;
    if(hf.Mass<0 &&
       edge.neighbors.second>0 &&
       edge.neighbors.second<tess.GetPointNo())
      return hf.Mass*
	cells[static_cast<size_t>(edge.neighbors.second)].tracers.find(name)->second;
    return 0;	      
  }
}

CustomFluxCalculator::CustomFluxCalculator
(const RiemannSolver& rs):
  rs_(rs) {}

vector<Extensive> CustomFluxCalculator::operator()
  (const Tessellation& tess,
   const vector<Vector2D>& point_velocities,
   const vector<ComputationalCell>& cells,
   const vector<Extensive>& /*extensive*/,
   const CacheData& /*cd*/,
   const EquationOfState& eos,
   const double /*time*/,
   const double /*dt*/) const
{
  vector<Extensive> res(tess.getAllEdges().size());
  for(size_t i=0;i<tess.getAllEdges().size();++i){
    const Conserved hydro_flux = 
      calcHydroFlux(tess, point_velocities,
		    cells, eos, i);
    res[i].mass = hydro_flux.Mass;
    res[i].momentum = hydro_flux.Momentum;
    res[i].energy = hydro_flux.Energy;
    for(boost::container::flat_map<std::string,double>::const_iterator it =
	  cells.front().tracers.begin();
	it!=cells.front().tracers.end();++it)
      res[i].tracers[it->first] = 
	calc_tracer_flux(tess.getAllEdges()[i],
			 tess,cells,it->first,hydro_flux);
  }
  return res;
}

const Conserved CustomFluxCalculator::calcHydroFlux
(const Tessellation& tess,
 const vector<Vector2D>& point_velocities,
 const vector<ComputationalCell>& cells,
 const EquationOfState& eos,
 const size_t i) const
{
  const Edge& edge = tess.GetEdge(static_cast<int>(i));
  const std::pair<bool,bool> flags 
    (edge.neighbors.first>=0 && edge.neighbors.first<tess.GetPointNo(),
     edge.neighbors.second>=0 && edge.neighbors.second<tess.GetPointNo());
  assert(flags.first || flags.second);
  if(!flags.first){
    const size_t right_index =
      static_cast<size_t>(edge.neighbors.second);
    const ComputationalCell& right_cell = cells[right_index];
    if(right_cell.stickers.find("obstacle")->second)
      return Conserved();
    const Vector2D p = Parallel(edge);
    const Primitive right = convert_to_primitive(right_cell,eos);
    //const Vector2D pos = tess.GetMeshPoint(edge.neighbors.second);
    const Primitive left = 
      CalcPrimitive(1,
		    1e-9,
		    Vector2D(0,1),
		    eos);
    const Vector2D n = remove_parallel_component
      (tess.GetMeshPoint(edge.neighbors.second) -
       edge.vertices.first, p);
    return rotate_solve_rotate_back
      (rs_, left, right, 0, n, p);
  }
  if(!flags.second){
    const size_t left_index = 
      static_cast<size_t>(edge.neighbors.first);
    const ComputationalCell& left_cell = cells[left_index];
    if(left_cell.stickers.find("obstacle")->second)
      return Conserved();
    const Primitive left = convert_to_primitive(left_cell, eos);
    const Vector2D p = Parallel(edge);
    //	const Vector2D pos = tess.GetMeshPoint(edge.neighbors.first);
    const Primitive right = 
      CalcPrimitive(1,
		    1e-9,
		    Vector2D(0,1),
		    eos);
    const Vector2D n = remove_parallel_component
      (edge.vertices.second - 
       tess.GetMeshPoint(edge.neighbors.first),p);
    return rotate_solve_rotate_back
      (rs_,left,right,0,n,p);
  }
  const size_t left_index =
    static_cast<size_t>(edge.neighbors.first);
  const size_t right_index =
    static_cast<size_t>(edge.neighbors.second);
  const ComputationalCell& left_cell = cells[left_index];
  const ComputationalCell& right_cell = cells[right_index];
  if(left_cell.stickers.find("obstacle")->second &&
     right_cell.stickers.find("obstacle")->second)
    return Conserved();
  const Vector2D p = Parallel(edge);
  const Vector2D n =
    tess.GetMeshPoint(edge.neighbors.second) -
    tess.GetMeshPoint(edge.neighbors.first);
  const double velocity = Projection
    (tess.CalcFaceVelocity
     (point_velocities[left_index],
      point_velocities[right_index],
      tess.GetCellCM(edge.neighbors.first),
      tess.GetCellCM(edge.neighbors.second),
      calc_centroid(edge)),n);
  if(left_cell.stickers.find("obstacle")->second){
    const Primitive right = convert_to_primitive(right_cell,eos);
    const Primitive left = reflect(right,p);
    return rotate_solve_rotate_back
      (rs_,left,right,velocity,n,p);
  }
  if(right_cell.stickers.find("obstacle")->second){
    const Primitive left = convert_to_primitive(left_cell,eos);
    const Primitive right = reflect(left,p);
    return rotate_solve_rotate_back
      (rs_,left,right,velocity,n,p);
  }
  const Primitive left =
    convert_to_primitive(left_cell,eos);
  const Primitive right = 
    convert_to_primitive(right_cell,eos);
  return rotate_solve_rotate_back
    (rs_,left,right,velocity,n,p);
}
