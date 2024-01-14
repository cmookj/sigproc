// ===========================================================================
// File: "hampel.hpp"
//                        Created: 2023-05-03 08:00:13
//              Last modification: 2023-05-03 08:45:50
// Author: Changmook Chun
// e-mail: cmookj@duck.com
// ===========================================================================

#ifndef GPW_SP_HAMPEL_H
#define GPW_SP_HAMPEL_H

#include <algorithm>
#include <iterator>
#include <vector>

namespace gpw::sigproc {

template <
    class random_access_iter_t,
    class T = typename random_access_iter_t::value_type>
T median (random_access_iter_t first, random_access_iter_t last) {
    // Make a copy of the original data and perform some preparation.
    std::vector<T> s (first, last);
    const auto halfsize = s.size() / 2;
    const auto center_odd = s.begin() + halfsize;
    const auto center_even = s.begin() + halfsize - (s.size() % 2 == 0 ? 1 : 0);

    std::nth_element (s.begin(), center_odd, s.end());
    std::nth_element (s.begin(), center_even, s.end());
    return (*center_odd + *center_even) / 2;
}

// Given a sequence x[1], x[2], ... ,x[n] and a sliding window of length k,
// define point-to-point median using
//   m[i] = median(x[i-k], ... , x[i], ... , x[i+k])
//
// Near the sequence endpoints, the window is truncated.
// For i < k + 1,
//   m[i] = median(x[1], ... , x[i], ... , x[i+k]).
// For i > n - k,
//   m[i] = median(x[i-k], ... , x[i], ... , x[n]).
//
template <
    class random_access_iter_t,
    class T = typename random_access_iter_t::value_type>
std::vector<T> medians_windowed (
    random_access_iter_t first,
    random_access_iter_t last,
    const std::size_t window_length
) {
    // Length of the input sequence is too short.  Return an empty vector.
    if (last - first < 2 * window_length + 1)
        return {};

    // Make a copy of the original data and perform some preparation.
    std::vector<T> s (first, last);

    // To make the indices concise:
    const std::size_t k = window_length;
    const std::size_t n = s.size();

    // Calculate local medians
    std::vector<T> medians (n);
    for (std::size_t i = 0; i < k; ++i) {
        medians[i] = median (s.begin(), s.begin() + i + k + 1);
    }
    for (std::size_t i = k; i < n - k; ++i) {
        auto iter_i = s.begin() + i;
        medians[i] = median (iter_i - k, iter_i + k + 1);
    }
    for (std::size_t i = n - k; i < n; ++i) {
        medians[i] = median (s.begin() + i - k, s.end());
    }

    return medians;
}

// Given a sequence x[1], x[2], ... ,x[n] with local medians m[1], m[2], ... ,
// m[n] and a sliding window of length k, define point-to-point
// standard-deviation estimates using
//   s[i] = k * median(|x[i-k] - m[i]|, ... , |x[i+k] - m[i]|),
//   where k = 1.4826
//
// Near the sequence endpoints, the window is truncated.
// For i < k + 1,
//   s[i] = k * median(|x[1] - m[i]|, ... , |x[i+k] - m[i]|).
// For i > n - k,
//   s[i] = k * median(|x[i-k] - m[i]|, ... , |x[n] - m[i]|).
//
template <
    class random_access_iter_t,
    class T = typename random_access_iter_t::value_type>
std::vector<double> stdevs_windowed (
    random_access_iter_t first,
    random_access_iter_t last,
    const std::size_t window_length
) {
    // Length of the input sequence is too short.  Return an empty vector.
    if (last - first < 2 * window_length + 1)
        return {};

    // Local medians windowed
    auto medians = medians_windowed (first, last, window_length);

    // Make a copy of the original data and perform some preparation.
    std::vector<T> s (first, last);

    // To make the indices concise
    const std::size_t k = window_length;
    const std::size_t n = s.size();
    constexpr double kappa = 1.0 / 0.67449;

    // Calculate standard deviations.
    std::vector<double> stdevs (n);
    // s[i] = k * median(|x[1] - m[i]|, ... , |x[i+k] - m[i]|).
    for (std::size_t i = 0; i < k; ++i) {
        std::vector<double> absolute_diffs;

        const T m = medians[i];
        std::transform (
            s.begin(),
            s.begin() + i + k + 1,
            std::back_inserter (absolute_diffs),
            [m] (const T& v) { return abs (v - m); }
        );

        stdevs[i] =
            kappa * median (absolute_diffs.begin(), absolute_diffs.end());
    }
    // s[i] = k * median(|x[i-k] - m[i]|, ... , |x[i+k] - m[i]|)
    for (std::size_t i = k; i < n - k; ++i) {
        auto iter_i = s.begin() + i;
        T med = medians[i];

        std::vector<double> absolute_diffs;
        std::transform (
            iter_i - k,
            iter_i + k + 1,
            std::back_inserter (absolute_diffs),
            [med] (const T& v) { return abs (v - med); }
        );
        stdevs[i] =
            kappa * median (absolute_diffs.begin(), absolute_diffs.end());
    }
    // s[i] = k * median(|x[i-k] - m[i]|, ... , |x[n] - m[i]|).
    for (std::size_t i = n - k; i < n; ++i) {
        std::vector<double> absolute_diffs;
        const T m = medians[i];
        std::transform (
            s.begin() + i - k,
            s.end(),
            std::back_inserter (absolute_diffs),
            [m] (const T& v) { return abs (v - m); }
        );

        stdevs[i] =
            kappa * median (absolute_diffs.begin(), absolute_diffs.end());
    }

    return stdevs;
}

namespace hampel {

// Given a sequence x[1], x[2], ... ,x[n] and a sliding window of length k,
// define point-to-point median and standard-deviation estimates using
//   m[i] = median(x[i-k], ... , x[i], ... , x[i+k])
//   s[i] = k * median(|x[i-k] - m[i]|, ... , |x[i+k] - m[i]|),
//   where k = 1.4826
//
// The Hampel filter is defined as
//   x[i] = m[i],  if |x[i] - m[i]| > th * s[i], for a given threshold th;
//          x[i],  otherwise.
//
// Near the sequence endpoints, the window is truncated.
// For i < k + 1,
//   m[i] = median(x[1], ... , x[i], ... , x[i+k]),
//   s[i] = k * median(|x[1] - m[i]|, ... , |x[i+k] - m[i]|).
// For i > n - k,
//   m[i] = median(x[i-k], ... , x[i], ... , x[n]),
//   s[i] = k * median(|x[i-k] - m[i]|, ... , |x[n] - m[i]|).
//
template <
    class random_access_iter_t,
    class T = typename random_access_iter_t::value_type>
std::vector<T> filter (
    random_access_iter_t first,
    random_access_iter_t last,
    const std::size_t window_length,
    const double threshold = 3.
) {
    // Length of the input sequence is too short.  Return empty vector.
    if (last - first < 2 * window_length + 1)
        return {};

    auto medians = medians_windowed (first, last, window_length);
    auto stdevs = stdevs_windowed (first, last, window_length);

    // Filter ------------------------------------------------------------------
    // Make a copy of the original data and perform some preparation.
    std::vector<T> filtered (first, last);
    std::size_t idx = 0;
    for (auto& val : filtered) {
        if (abs (val - medians[idx]) > threshold * stdevs[idx]) {
            val = medians[idx];
        }

        idx++;
    }

    return filtered;
}

// Given a sequence x[1], x[2], ... ,x[n] and a sliding window of length k,
// define point-to-point median and standard-deviation estimates using
//   m[i] = median(x[i-k], ... , x[i], ... , x[i+k])
//   s[i] = k * median(|x[i-k] - m[i]|, ... , |x[i+k] - m[i]|),
//   where k = 1.4826
//
// x[i] is an outlier if |x[i] - m[i]| > th * s[i], for a given threshold th.
//
// Near the sequence endpoints, the window is truncated.
// For i < k + 1,
//   m[i] = median(x[1], ... , x[i], ... , x[i+k]),
//   s[i] = k * median(|x[1] - m[i]|, ... , |x[i+k] - m[i]|).
// For i > n - k,
//   m[i] = median(x[i-k], ... , x[i], ... , x[n]),
//   s[i] = k * median(|x[i-k] - m[i]|, ... , |x[n] - m[i]|).
//
template <
    class random_access_iter_t,
    class T = typename random_access_iter_t::value_type>
std::vector<bool> outlier_indices (
    random_access_iter_t first,
    random_access_iter_t last,
    const std::size_t window_length,
    const double threshold = 3.
) {
    // Length of the input sequence is too short.  Return empty vector.
    if (last - first < 2 * window_length + 1)
        return {};

    auto medians = medians_windowed (first, last, window_length);
    auto stdevs = stdevs_windowed (first, last, window_length);

    std::vector<bool> flags;
    flags.reserve (last - first);
    std::size_t idx = 0;
    std::transform (
        first,
        last,
        std::back_inserter (flags),
        [&medians, &stdevs, &idx, threshold] (const T& val) {
            bool is_outlier = false;
            if (abs (val - medians[idx]) > threshold * stdevs[idx]) {
                is_outlier = true;
            }

            idx++;
            return is_outlier;
        }
    );

    return flags;
}

} // namespace hampel
} // namespace gpw::sigproc

#endif
