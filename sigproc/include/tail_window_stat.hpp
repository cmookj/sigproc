#include <unordered_map>
#include <vector>

namespace gpw::sigproc {

template <typename T> class tail_window_statistics {
public:
  tail_window_statistics<T> (const T ws)

      : count_values{0}
      , next_index{0}
      , sum_of_values{0}
      , buffer{std::vector<T> (ws)}
      , hist{std::unordered_map<T, int>()} {}

  void
  update (const T value) {
    next_index = (next_index + 1) % buffer.size();

    // Update sum and histogram for fast computation of algebraic mean and mode.
    sum_of_values += value;
    if (++count_values > static_cast<int> (buffer.size())) {
      sum_of_values -= buffer[next_index];
      if (hist[buffer[next_index]] > 0) {
        hist[buffer[next_index]]--;
      }
      count_values = static_cast<int> (buffer.size());
    }
    buffer[next_index] = value;
    hist[value]++;
  }

  double
  local_mean () const {
    return static_cast<double> (sum_of_values) /
           static_cast<double> (count_values);
  }

  int
  local_mode () const {
    int mode;
    int frequency{0};
    for (unsigned i = 0; i < count_values; ++i) {
      if (hist.find (buffer[i]) != hist.end() &&
          hist.at (buffer[i]) > frequency) {
        mode      = buffer[i];
        frequency = hist.at (buffer[i]);
      }
    }

    return mode;
  }

private:
  int count_values;
  int next_index;
  T   sum_of_values;

  std::vector<T>             buffer;
  std::unordered_map<T, int> hist;
};

}  // namespace gpw::sigproc
