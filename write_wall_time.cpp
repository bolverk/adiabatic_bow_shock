#include "write_wall_time.hpp"
#include "source/misc/simple_io.hpp"

WriteWallTime::WriteWallTime(const string& fname):
  fname_(fname), start_(clock()) {}

void WriteWallTime::operator()(const hdsim& /*sim*/)
{
  write_number((clock()-start_)/CLOCKS_PER_SEC,fname_);
}
