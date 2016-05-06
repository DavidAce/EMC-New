#include <bitset>
#include <iostream>
#include <random>
#include "constants.hpp"
#include "randomFunctions.hpp"
#include "DNA.hpp"
#include "personality.hpp"
#include "../Eigen/Dense"
#include "../Eigen/Core"

using namespace std;
using namespace constants;
using namespace Eigen;


std::ostream &operator<<(std::ostream &os, personality const &guy) {
	for (int i = 0; i < nGenes; i++) {
		os << guy.H << endl;
	}
	return os;
}
personality::personality() : born(0) {
}

personality::personality(bool dna) : born(0), genome(dna) {

}

