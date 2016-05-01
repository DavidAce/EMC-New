#ifndef RANDOMFUNCTIONS_H   // if x.h hasn't been included yet
#define RANDOMFUNCTIONS_H   //  #define this so the compiler knows it has been included
#include <random>

using namespace std;

typedef mt19937 RNGType;
//RNGType; //Forward declaration?
extern RNGType rng;
extern double uniform_double(RNGType *rn, const double lowerLimit, const double upperLimit);
extern int uniform_integer(RNGType *rn, const int lowerLimit, const int upperLimit);
extern double gaussian_truncated(RNGType *rn, const double lowerLimit, const double upperLimit, const double mean, const double std);

#endif