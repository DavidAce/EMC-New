#include <random>
#include "randomFunctions.hpp"

using namespace std;

double uniform_double(RNGType *rn, const double lowerLimit, const double upperLimit) {
	uniform_real_distribution<> dice(lowerLimit, upperLimit);
	return dice(*rn);
}
int uniform_integer(RNGType *rn, const int lowerLimit, const int upperLimit) {
	uniform_int_distribution<> dice(lowerLimit, upperLimit);
	return dice(*rn);
}

double gaussian_truncated(RNGType *rn, const double lowerLimit, const double upperLimit, const double mean, const double std) {
	normal_distribution<double> distribution(mean,std);
	double ul = fmax(lowerLimit, upperLimit);
	double ll = fmin(lowerLimit, upperLimit);
	double number;
	while (true) {
		number = distribution(*rn);
		if (number >= ll && number <= ul) {
			return number;
		}
	}
}
//class gaussian_truncated {
//	std::default_random_engine generator;
//	std::normal_distribution<double> distribution;
//	double min;
//	double max;
//public:
//	gaussian_truncated(double mean, double stddev, double min, double max) :
//		distribution(mean, stddev), min(min), max(max)
//	{}
//
//	double operator ()(double mean, double stddev, double min, double max) {
//		while (true) {
//			double number = this->distribution(generator);
//			if (number >= this->min && number <= this->max)
//				return number;
//		}
//	}
//};