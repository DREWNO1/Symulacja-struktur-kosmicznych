#ifndef RANDOMUTILS_H
#define RANDOMUTILS_H

#include <random> // Dla std::mt19937, std::uniform_real_distribution

// Globalny silnik losowy - deklaracja, definicja może być w .cpp lub jako inline
std::mt19937& get_random_engine();

namespace RandomUtils {

double generateRandomNumberFromTo(double min_val, double max_val);
float generateRandomFloatFromTo(float min_val, float max_val);

} // namespace RandomUtils

#endif // RANDOMUTILS_H