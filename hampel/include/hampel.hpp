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
#include <vector>

namespace gpw::sp {

template <
    class random_access_iter_t,
    class T = typename random_access_iter_t::value_type>
bool hampel_identifier (
    random_access_iter_t first, random_access_iter_t last, double threshold = 3
) {
    // Make a copy of the original data and perform some preparation.
    std::vector<T> s (first, last);
    const auto halfsize = s.size() / 2;
    const auto center = s.begin() + halfsize;

    // Find median in the original sequence.
    std::nth_element (s.begin(), center, s.end());
    const T med = *center;

    // Find residuals with median value.
    s.resize (0);
    std::transform (first, last, std::back_inserter (s), [=] (T x) {
        return abs (x - med);
    });

    // Find a MAD and standard deviation.
    constexpr double k = 1.4826;
    std::nth_element (s.begin(), center, s.end()); // center = MAD
    const T sdev = *center * k;

    // Return whether the center value is an outlier or not.
    // True: the value is NOT an outlier.
    // False: the value is an outlier.
    return abs (*(first + halfsize) - med) <= abs (threshold * sdev);
}

template <
    class random_access_iter_t,
    class T = typename random_access_iter_t::value_type>
T hampel_filter (
    random_access_iter_t first, random_access_iter_t last, double threshold = 3
) {
    // Make a copy of the original data and perform some preparation.
    std::vector<T> s (first, last);
    const auto halfsize = s.size() / 2;
    const auto center = s.begin() + halfsize;

    // Find median in the original sequence.
    std::nth_element (s.begin(), center, s.end());
    const T med = *center;

    // Find residuals with median value.
    s.resize (0);
    std::transform (first, last, std::back_inserter (s), [=] (T x) {
        return abs (x - med);
    });

    // Find a MAD and standard deviation.
    constexpr double k = 1.4826;
    std::nth_element (s.begin(), center, s.end()); // center = MAD
    const T sdev = *center * k;

    // Return filtered value for central element in the window.
    return abs (*(first + halfsize) - med) <= abs (threshold * sdev)
               ? *(first + halfsize)
               : med;
}

} // namespace gpw::sp

#endif
