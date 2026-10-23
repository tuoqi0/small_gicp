#include <small_gicp/benchmark/benchmark_odom.hpp>

int main() {
  small_gicp::Stopwatch timer;
  timer.start();
  timer.stop();
  return timer.sec() < 0;
}
