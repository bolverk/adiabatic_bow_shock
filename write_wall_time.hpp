#ifndef WRITE_WALL_TIME_HPP
#define WRITE_WALL_TIME_HPP 1

#include <ctime>
#include "source/newtonian/test_2d/main_loop_2d.hpp"

class WriteWallTime: public DiagnosticFunction
{
public:

  WriteWallTime(const string& fname);

  void operator()(const hdsim& sim);

private:
const string fname_;
  const clock_t start_;
};

#endif // WRITE_WALL_TIME_HPP
