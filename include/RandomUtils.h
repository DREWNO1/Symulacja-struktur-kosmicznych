#ifndef RANDOMUTILS_H
#define RANDOMUTILS_H

#include <random> 

std::mt19937& get_random_engine();

namespace RandomUtils {

double generateRandomNumberFromTo(double min_val, double max_val);
float generateRandomFloatFromTo(float min_val, float max_val);

} 

#endif 