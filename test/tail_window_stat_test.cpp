#include "tail_window_stat.hpp"

#include <gtest/gtest.h>
#include <vector>

using namespace gpw::sigproc;

TEST (Statistics, MeanMode) {
  std::vector<int>    x{1, 2, 2, 3, -5};
  std::vector<double> means{1, 1.5, 5. / 3., 7. / 3., 0.};
  std::vector<int>    modes{1, 1, 2, 2, 2};

  gpw::sigproc::tail_window_statistics<int> is{3};

  for (unsigned i = 0; i < 5; ++i) {
    is.update (x[i]);
    EXPECT_FLOAT_EQ (is.local_mean(), means[i]);
    EXPECT_EQ (is.local_mode(), modes[i]);
  }
}
