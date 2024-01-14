#include "hampel.hpp"
#include <gtest/gtest.h>
#include <rbm/rbm.hpp>

using namespace gpw::sigproc;
using namespace gpw::geometry;

TEST (Statistics, Median) {
    std::vector<int> x {2, 4, 6, 8, 10};
    EXPECT_EQ (median (x.begin(), x.end()), 6);

    EXPECT_EQ (median (x.begin(), x.begin() + 1), 2);
    EXPECT_EQ (median (x.begin(), x.begin() + 2), 3);
    EXPECT_EQ (median (x.begin(), x.begin() + 3), 4);
    EXPECT_EQ (median (x.begin(), x.begin() + 4), 5);
}

TEST (Statistics, Median_Windowed) {
    std::vector<int> ix {1, 3, 5, 7, 9, 10, 8, 6, 4, 2,
                         1, 3, 5, 7, 9, 10, 8, 6, 4, 2};

    std::vector<int> ixm {4, 5, 6, 7, 7, 7, 7, 6, 4, 4,
                          4, 4, 5, 7, 7, 7, 7, 7, 6, 5};

    EXPECT_EQ (medians_windowed (ix.begin(), ix.end(), 3), ixm);

    std::vector<double> x {
        0.062790519529313,
        0.125333233564304,
        0.187381314585725,
        0.248689887164855,
        0.309016994374947,
        0.368124552684678,
        0.425779291565073,
        0.481753674101715,
        0.535826794978997,
        0.587785252292473,
        0.637423989748690,
        0.684547105928689,
        0.728968627421412,
        0.770513242775789,
        0.809016994374947,
        0.844327925502015,
        0.876306680043864,
        0.904827052466020,
        0.929776485888251,
        0.951056516295154,
        0.968583161128631,
        0.982287250728689,
        0.992114701314478,
        0.998026728428272,
        1.0,
        0.998026728428272,
        0.992114701314478,
        0.982287250728689,
        0.968583161128631,
        0.951056516295154,
        0.929776485888251,
        0.904827052466019,
        0.876306680043863,
        0.844327925502015,
        0.809016994374947,
        0.770513242775789,
        0.728968627421411,
        0.684547105928689,
        0.637423989748690,
        0.587785252292473,
        0.535826794978997,
        0.481753674101715,
        0.425779291565073,
        0.368124552684678,
        0.309016994374948,
        0.248689887164855,
        0.187381314585725,
        0.125333233564305,
        0.062790519529314,
        1.224646799147353e-16,
        -0.062790519529313,
        -0.125333233564304,
        -0.187381314585725,
        -0.248689887164855,
        -0.309016994374947,
        -0.368124552684678,
        -0.425779291565072,
        -0.481753674101715,
        -0.535826794978996,
        -0.587785252292473,
        -0.637423989748690,
        -0.684547105928688,
        -0.728968627421412,
        -0.770513242775789,
        -0.809016994374947,
        -0.844327925502015,
        -0.876306680043864,
        -0.904827052466020,
        -0.929776485888251,
        -0.951056516295154,
        -0.968583161128631,
        -0.982287250728689,
        -0.992114701314478,
        -0.998026728428272,
        -1.0,
        -0.998026728428272,
        -0.992114701314478,
        -0.982287250728689,
        -0.968583161128631,
        -0.951056516295154,
        -0.929776485888252,
        -0.904827052466020,
        -0.876306680043863,
        -0.844327925502015,
        -0.809016994374948,
        -0.770513242775790,
        -0.728968627421412,
        -0.684547105928689,
        -0.637423989748690,
        -0.587785252292473,
        -0.535826794978997,
        -0.481753674101716,
        -0.425779291565073,
        -0.368124552684678,
        -0.309016994374948,
        -0.248689887164855,
        -0.187381314585725,
        -0.125333233564305,
        -0.062790519529313,
        -2.449293598294706e-16
    };

    std::vector<double> meds3 {
        0.156357274075014,  0.187381314585725,     0.218035600875290,
        0.248689887164855,  0.309016994374947,     0.368124552684678,
        0.425779291565073,  0.481753674101715,     0.535826794978997,
        0.587785252292473,  0.637423989748690,     0.684547105928689,
        0.728968627421412,  0.770513242775789,     0.809016994374947,
        0.844327925502015,  0.876306680043864,     0.904827052466020,
        0.929776485888251,  0.951056516295154,     0.968583161128631,
        0.982287250728689,  0.992114701314478,     0.992114701314478,
        0.992114701314478,  0.992114701314478,     0.992114701314478,
        0.982287250728689,  0.968583161128631,     0.951056516295154,
        0.929776485888251,  0.904827052466019,     0.876306680043863,
        0.844327925502015,  0.809016994374947,     0.770513242775789,
        0.728968627421411,  0.684547105928689,     0.637423989748690,
        0.587785252292473,  0.535826794978997,     0.481753674101715,
        0.425779291565073,  0.368124552684678,     0.309016994374948,
        0.248689887164855,  0.187381314585725,     0.125333233564305,
        0.062790519529314,  1.224646799147353e-16, -0.062790519529313,
        -0.125333233564304, -0.187381314585725,    -0.248689887164855,
        -0.309016994374947, -0.368124552684678,    -0.425779291565072,
        -0.481753674101715, -0.535826794978996,    -0.587785252292473,
        -0.637423989748690, -0.684547105928688,    -0.728968627421412,
        -0.770513242775789, -0.809016994374947,    -0.844327925502015,
        -0.876306680043864, -0.904827052466020,    -0.929776485888251,
        -0.951056516295154, -0.968583161128631,    -0.982287250728689,
        -0.992114701314478, -0.992114701314478,    -0.992114701314478,
        -0.992114701314478, -0.992114701314478,    -0.982287250728689,
        -0.968583161128631, -0.951056516295154,    -0.929776485888252,
        -0.904827052466020, -0.876306680043863,    -0.844327925502015,
        -0.809016994374948, -0.770513242775790,    -0.728968627421412,
        -0.684547105928689, -0.637423989748690,    -0.587785252292473,
        -0.535826794978997, -0.481753674101716,    -0.425779291565073,
        -0.368124552684678, -0.309016994374948,    -0.248689887164855,
        -0.187381314585725, -0.156357274075015,    -0.125333233564305,
        -0.094061876546809
    };

    std::vector<double> meds10 {
        0.368124552684678,  0.396951922124875,     0.425779291565073,
        0.453766482833394,  0.481753674101715,     0.508790234540356,
        0.535826794978997,  0.561806023635735,     0.587785252292473,
        0.612604621020581,  0.637423989748690,     0.684547105928689,
        0.728968627421412,  0.770513242775789,     0.809016994374947,
        0.844327925502015,  0.876306680043864,     0.904827052466020,
        0.929776485888251,  0.951056516295154,     0.951056516295154,
        0.951056516295154,  0.951056516295154,     0.951056516295154,
        0.951056516295154,  0.951056516295154,     0.951056516295154,
        0.951056516295154,  0.951056516295154,     0.951056516295154,
        0.929776485888251,  0.904827052466019,     0.876306680043863,
        0.844327925502015,  0.809016994374947,     0.770513242775789,
        0.728968627421411,  0.684547105928689,     0.637423989748690,
        0.587785252292473,  0.535826794978997,     0.481753674101715,
        0.425779291565073,  0.368124552684678,     0.309016994374948,
        0.248689887164855,  0.187381314585725,     0.125333233564305,
        0.062790519529314,  1.224646799147353e-16, -0.062790519529313,
        -0.125333233564304, -0.187381314585725,    -0.248689887164855,
        -0.309016994374947, -0.368124552684678,    -0.425779291565072,
        -0.481753674101715, -0.535826794978996,    -0.587785252292473,
        -0.637423989748690, -0.684547105928688,    -0.728968627421412,
        -0.770513242775789, -0.809016994374947,    -0.844327925502015,
        -0.876306680043864, -0.904827052466020,    -0.929776485888251,
        -0.951056516295154, -0.951056516295154,    -0.951056516295154,
        -0.951056516295154, -0.951056516295154,    -0.951056516295154,
        -0.951056516295154, -0.951056516295154,    -0.951056516295154,
        -0.951056516295154, -0.951056516295154,    -0.929776485888252,
        -0.904827052466020, -0.876306680043863,    -0.844327925502015,
        -0.809016994374948, -0.770513242775790,    -0.728968627421412,
        -0.684547105928689, -0.637423989748690,    -0.587785252292473,
        -0.561806023635735, -0.535826794978997,    -0.508790234540357,
        -0.481753674101716, -0.453766482833395,    -0.425779291565073,
        -0.396951922124875, -0.368124552684678,    -0.338570773529813,
        -0.309016994374948
    };

    auto meds3_calculated = medians_windowed (x.begin(), x.end(), 3);
    auto meds10_calculated = medians_windowed (x.begin(), x.end(), 10);

    for (size_t i = 0; i < meds3.size(); ++i)
        EXPECT_FLOAT_EQ (meds3[i], meds3_calculated[i]);

    for (size_t i = 0; i < meds10.size(); ++i)
        EXPECT_FLOAT_EQ (meds10[i], meds10_calculated[i]);
}

TEST (Statistics, Stdev_Windowed) {
    std::vector<double> x {
        0.062790519529313,
        0.125333233564304,
        0.187381314585725,
        0.248689887164855,
        0.309016994374947,
        0.368124552684678,
        0.425779291565073,
        0.481753674101715,
        0.535826794978997,
        0.587785252292473,
        0.637423989748690,
        0.684547105928689,
        0.728968627421412,
        0.770513242775789,
        0.809016994374947,
        0.844327925502015,
        0.876306680043864,
        0.904827052466020,
        0.929776485888251,
        0.951056516295154,
        0.968583161128631,
        0.982287250728689,
        0.992114701314478,
        0.998026728428272,
        1.0,
        0.998026728428272,
        0.992114701314478,
        0.982287250728689,
        0.968583161128631,
        0.951056516295154,
        0.929776485888251,
        0.904827052466019,
        0.876306680043863,
        0.844327925502015,
        0.809016994374947,
        0.770513242775789,
        0.728968627421411,
        0.684547105928689,
        0.637423989748690,
        0.587785252292473,
        0.535826794978997,
        0.481753674101715,
        0.425779291565073,
        0.368124552684678,
        0.309016994374948,
        0.248689887164855,
        0.187381314585725,
        0.125333233564305,
        0.062790519529314,
        1.224646799147353e-16,
        -0.062790519529313,
        -0.125333233564304,
        -0.187381314585725,
        -0.248689887164855,
        -0.309016994374947,
        -0.368124552684678,
        -0.425779291565072,
        -0.481753674101715,
        -0.535826794978996,
        -0.587785252292473,
        -0.637423989748690,
        -0.684547105928688,
        -0.728968627421412,
        -0.770513242775789,
        -0.809016994374947,
        -0.844327925502015,
        -0.876306680043864,
        -0.904827052466020,
        -0.929776485888251,
        -0.951056516295154,
        -0.968583161128631,
        -0.982287250728689,
        -0.992114701314478,
        -0.998026728428272,
        -1.0,
        -0.998026728428272,
        -0.992114701314478,
        -0.982287250728689,
        -0.968583161128631,
        -0.951056516295154,
        -0.929776485888252,
        -0.904827052466020,
        -0.876306680043863,
        -0.844327925502015,
        -0.809016994374948,
        -0.770513242775790,
        -0.728968627421412,
        -0.684547105928689,
        -0.637423989748690,
        -0.587785252292473,
        -0.535826794978997,
        -0.481753674101716,
        -0.425779291565073,
        -0.368124552684678,
        -0.309016994374948,
        -0.248689887164855,
        -0.187381314585725,
        -0.125333233564305,
        -0.062790519529313,
        -2.449293598294706e-16
    };

    std::vector<double> stdevs3 {
        0.091444424147802, 0.091992622576373, 0.136164975640656,
        0.177074100066164, 0.173112040851890, 0.168466787499742,
        0.163156672702487, 0.157202653057268, 0.150628226359493,
        0.143459338867767, 0.135724282905869, 0.127453585205868,
        0.118679886433063, 0.109437812368190, 0.099763837255280,
        0.089696139854485, 0.079274452767953, 0.068539905633398,
        0.057534862804213, 0.046302756156718, 0.034887913684390,
        0.023335384555521, 0.011690761324735, 0.008765184514776,
        0.008765184514776, 0.008765184514776, 0.011690761324734,
        0.023335384555521, 0.034887913684390, 0.046302756156718,
        0.057534862804213, 0.068539905633399, 0.079274452767953,
        0.089696139854485, 0.099763837255279, 0.109437812368189,
        0.118679886433063, 0.127453585205867, 0.135724282905869,
        0.143459338867768, 0.150628226359492, 0.157202653057268,
        0.163156672702487, 0.168466787499742, 0.173112040851890,
        0.177074100066164, 0.180337328704939, 0.182888848295604,
        0.184718589156012, 0.185819330134918, 0.184718589156012,
        0.182888848295603, 0.180337328704938, 0.177074100066165,
        0.173112040851890, 0.168466787499742, 0.163156672702487,
        0.157202653057267, 0.150628226359493, 0.143459338867768,
        0.135724282905869, 0.127453585205868, 0.118679886433063,
        0.109437812368190, 0.099763837255280, 0.089696139854485,
        0.079274452767953, 0.068539905633398, 0.057534862804213,
        0.046302756156719, 0.034887913684390, 0.023335384555521,
        0.011690761324735, 0.008765184514776, 0.008765184514776,
        0.008765184514776, 0.011690761324734, 0.023335384555521,
        0.034887913684390, 0.046302756156718, 0.057534862804213,
        0.068539905633398, 0.079274452767954, 0.089696139854486,
        0.099763837255279, 0.109437812368189, 0.118679886433062,
        0.127453585205868, 0.135724282905870, 0.143459338867768,
        0.150628226359492, 0.157202653057267, 0.163156672702487,
        0.168466787499743, 0.173112040851890, 0.177074100066163,
        0.180337328704939, 0.137807407437622, 0.092725966579640,
        0.092359294578006
    };

    std::vector<double> stdevs10 {
        0.248635716473921, 0.251371770311587, 0.262553143837599,
        0.288169021449789, 0.300661991925035, 0.311310611409869,
        0.336268713554377, 0.337972659886029, 0.327998671615335,
        0.353010088499417, 0.354168006594212, 0.326587537428542,
        0.297718176456296, 0.270912924073635, 0.254404169338932,
        0.219109001685347, 0.171697229256557, 0.138178046345638,
        0.101187961003379, 0.069638140712240, 0.068539905633398,
        0.068539905633398, 0.068539905633398, 0.068539905633398,
        0.068539905633398, 0.068539905633398, 0.068539905633398,
        0.068539905633398, 0.068539905633399, 0.069638140712240,
        0.101187961003379, 0.138178046345638, 0.171697229256557,
        0.219109001685348, 0.254404169338932, 0.270912924073635,
        0.297718176456296, 0.326587537428542, 0.354168006594211,
        0.380350736441824, 0.405032395698424, 0.428115577130903,
        0.445441992133202, 0.452688914941406, 0.458149281416230,
        0.461801541985801, 0.463631282846210, 0.463631282846210,
        0.461801541985801, 0.458149281416230, 0.461801541985801,
        0.463631282846210, 0.463631282846210, 0.461801541985801,
        0.458149281416230, 0.452688914941406, 0.445441992133201,
        0.428115577130903, 0.405032395698425, 0.380350736441825,
        0.354168006594211, 0.326587537428543, 0.297718176456296,
        0.270912924073636, 0.254404169338932, 0.219109001685347,
        0.171697229256557, 0.138178046345638, 0.101187961003379,
        0.069638140712240, 0.068539905633398, 0.068539905633398,
        0.068539905633398, 0.068539905633398, 0.068539905633398,
        0.068539905633398, 0.068539905633398, 0.068539905633398,
        0.068539905633398, 0.069638140712240, 0.101187961003379,
        0.138178046345638, 0.171697229256557, 0.219109001685348,
        0.254404169338932, 0.270912924073635, 0.297718176456295,
        0.326587537428542, 0.354168006594211, 0.380350736441824,
        0.370650554626400, 0.347946648156724, 0.356031162902723,
        0.345540887565907, 0.323101439745471, 0.313784899061979,
        0.296819883171202, 0.267970325785394, 0.258303021129658,
        0.256099784580199
    };

    auto stdevs3_calculated = stdevs_windowed (x.begin(), x.end(), 3);
    auto stdevs10_calculated = stdevs_windowed (x.begin(), x.end(), 10);

    for (size_t i = 0; i < stdevs3.size(); ++i)
        EXPECT_NEAR (stdevs3[i], stdevs3_calculated[i], 0.000001);

    for (size_t i = 0; i < stdevs10.size(); ++i)
        EXPECT_NEAR (stdevs10[i], stdevs10_calculated[i], 0.000001);
}

TEST (Hample, Filter) {
    std::vector<double> x {
        0.062790519529313,
        0.125333233564304,
        0.187381314585725,
        0.248689887164855,
        0.309016994374947,
        0.368124552684678,
        0.425779291565073,
        0.481753674101715,
        0.535826794978997,
        0.587785252292473,
        0.637423989748690,
        0.684547105928689,
        0.728968627421412,
        0.770513242775789,
        0.809016994374947,
        0.844327925502015,
        0.876306680043864,
        0.904827052466020,
        0.929776485888251,
        0.951056516295154,
        0.968583161128631,
        0.982287250728689,
        0.992114701314478,
        0.998026728428272,
        1.0,
        0.998026728428272,
        0.992114701314478,
        0.982287250728689,
        0.968583161128631,
        0.951056516295154,
        0.929776485888251,
        0.904827052466019,
        0.876306680043863,
        0.844327925502015,
        0.809016994374947,
        0.770513242775789,
        0.728968627421411,
        0.684547105928689,
        0.637423989748690,
        0.587785252292473,
        0.535826794978997,
        0.481753674101715,
        0.425779291565073,
        0.368124552684678,
        0.309016994374948,
        0.248689887164855,
        0.187381314585725,
        0.125333233564305,
        0.062790519529314,
        1.224646799147353e-16,
        -0.062790519529313,
        -0.125333233564304,
        -0.187381314585725,
        -0.248689887164855,
        -0.309016994374947,
        -0.368124552684678,
        -0.425779291565072,
        -0.481753674101715,
        -0.535826794978996,
        -0.587785252292473,
        -0.637423989748690,
        -0.684547105928688,
        -0.728968627421412,
        -0.770513242775789,
        -0.809016994374947,
        -0.844327925502015,
        -0.876306680043864,
        -0.904827052466020,
        -0.929776485888251,
        -0.951056516295154,
        -0.968583161128631,
        -0.982287250728689,
        -0.992114701314478,
        -0.998026728428272,
        -1.0,
        -0.998026728428272,
        -0.992114701314478,
        -0.982287250728689,
        -0.968583161128631,
        -0.951056516295154,
        -0.929776485888252,
        -0.904827052466020,
        -0.876306680043863,
        -0.844327925502015,
        -0.809016994374948,
        -0.770513242775790,
        -0.728968627421412,
        -0.684547105928689,
        -0.637423989748690,
        -0.587785252292473,
        -0.535826794978997,
        -0.481753674101716,
        -0.425779291565073,
        -0.368124552684678,
        -0.309016994374948,
        -0.248689887164855,
        -0.187381314585725,
        -0.125333233564305,
        -0.062790519529313,
        -2.449293598294706e-16
    };
    x[5] = 2.0;
    x[19] = -2.0;

    std::vector<double> y {
        0.062790519529313,
        0.125333233564304,
        0.187381314585725,
        0.248689887164855,
        0.309016994374947,
        0.425779291565073,
        0.425779291565073,
        0.481753674101715,
        0.535826794978997,
        0.587785252292473,
        0.637423989748690,
        0.684547105928689,
        0.728968627421412,
        0.770513242775789,
        0.809016994374947,
        0.844327925502015,
        0.876306680043864,
        0.904827052466020,
        0.929776485888251,
        0.929776485888251,
        0.968583161128631,
        0.982287250728689,
        0.992114701314478,
        0.998026728428272,
        1.0,
        0.998026728428272,
        0.992114701314478,
        0.982287250728689,
        0.968583161128631,
        0.951056516295154,
        0.929776485888251,
        0.904827052466019,
        0.876306680043863,
        0.844327925502015,
        0.809016994374947,
        0.770513242775789,
        0.728968627421411,
        0.684547105928689,
        0.637423989748690,
        0.587785252292473,
        0.535826794978997,
        0.481753674101715,
        0.425779291565073,
        0.368124552684678,
        0.309016994374948,
        0.248689887164855,
        0.187381314585725,
        0.125333233564305,
        0.062790519529314,
        1.224646799147353e-16,
        -0.062790519529313,
        -0.125333233564304,
        -0.187381314585725,
        -0.248689887164855,
        -0.309016994374947,
        -0.368124552684678,
        -0.425779291565072,
        -0.481753674101715,
        -0.535826794978996,
        -0.587785252292473,
        -0.637423989748690,
        -0.684547105928688,
        -0.728968627421412,
        -0.770513242775789,
        -0.809016994374947,
        -0.844327925502015,
        -0.876306680043864,
        -0.904827052466020,
        -0.929776485888251,
        -0.951056516295154,
        -0.968583161128631,
        -0.982287250728689,
        -0.992114701314478,
        -0.998026728428272,
        -1.0,
        -0.998026728428272,
        -0.992114701314478,
        -0.982287250728689,
        -0.968583161128631,
        -0.951056516295154,
        -0.929776485888252,
        -0.904827052466020,
        -0.876306680043863,
        -0.844327925502015,
        -0.809016994374948,
        -0.770513242775790,
        -0.728968627421412,
        -0.684547105928689,
        -0.637423989748690,
        -0.587785252292473,
        -0.535826794978997,
        -0.481753674101716,
        -0.425779291565073,
        -0.368124552684678,
        -0.309016994374948,
        -0.248689887164855,
        -0.187381314585725,
        -0.125333233564305,
        -0.062790519529313,
        -2.449293598294706e-16
    };

    auto filtered_x = hampel::filter (x.begin(), x.end(), 3, 3.0);

    for (std::size_t i = 0; i < x.size(); ++i) {
        EXPECT_FLOAT_EQ (filtered_x[i], y[i]);
    }

    auto flags = hampel::outlier_indices (x.begin(), x.end(), 3, 3.0);
    for (std::size_t i = 0; i < flags.size(); ++i) {
        if (i != 5 && i != 19)
            EXPECT_FALSE (flags[i]);
        else
            EXPECT_TRUE (flags[i]);
    }
}
