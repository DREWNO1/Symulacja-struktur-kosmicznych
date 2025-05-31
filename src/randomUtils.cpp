#include "RandomUtils.h" 
#include <algorithm> 




std::mt19937& get_random_engine() {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    return gen;
}

namespace RandomUtils {

double generateRandomNumberFromTo(double min_val, double max_val) {
    if (min_val > max_val) std::swap(min_val, max_val);
    std::uniform_real_distribution<double> distrib(min_val, max_val);
    return distrib(get_random_engine());
}

float generateRandomFloatFromTo(float min_val, float max_val) {
    if (min_val > max_val) std::swap(min_val, max_val);
    std::uniform_real_distribution<float> distrib(min_val, max_val);
    return distrib(get_random_engine());
}

} 