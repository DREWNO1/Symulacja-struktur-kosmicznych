#include "RandomUtils.h" // Zakładamy, że RandomUtils.h jest już zdefiniowany zgodnie z poprzednią sugestią
#include <algorithm> // Dla std::swap

// Definicja globalnego silnika losowego
// (Jeśli nie chcesz, aby był globalny, można go uczynić statycznym polem klasy RandomUtils
// lub przekazywać jako argument do funkcji)
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

} // namespace RandomUtils