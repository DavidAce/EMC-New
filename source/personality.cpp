#include <bitset>
#include <iostream>
#include <random>
#include "..\Eigen\Dense"
#include "..\Eigen\Core"
#include "constants.h"
#include "randomFunctions.h"
#include "DNA.h"
#include "personality.h"


using namespace std;
using namespace constants;
using namespace Eigen;


std::ostream &operator<<(std::ostream &os, personality const &guy) {
	for (int i = 0; i < nGenes; i++) {
		os << guy.H << endl;
	}
	return os;
}
personality::personality() : generation(0) {
}

personality::personality(bool dna) : generation(0), genome(dna) {

}

