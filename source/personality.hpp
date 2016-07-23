
#ifndef PERSONALITY_H   // if x.h hasn't been included yet...
#define PERSONALITY_H   //  #define this so the compiler knows it has been included
#include <bitset>
#include <iostream>
#include <random>
#include "constants.hpp"
#include "DNA.hpp"
#include <Eigen/Dense>
#include <Eigen/Core>

using namespace std;
using namespace constants;
using namespace Eigen;
class DNA;

class personality {
private:
public:

	personality(); //start generation count

	personality(bool); //TAKE CARE OF THIS!! NEEDS IN BOUNDARIES!
	double H;							//Fitness, or energy
	double t; 							//temperature
	double value;						//Actual fitting-value
	int born;							//Generation when DNA first emerged
	DNA genome;							//Contains binary and real representation of parameters
	int operator()() { //Return the bit at a.
		return 0;
	}
	friend ostream &operator<<(std::ostream &os, personality const &);

};
#endif
