#ifndef CONSTANTS_H   // if x.h hasn't been included yet
#define CONSTANTS_H   //  #define this so the compiler knows it has been included
#include "../Eigen/Dense"
#include "../Eigen/Core"
#include "mymath.h"
using namespace std;
using namespace Eigen;
namespace constants {
	//Evolutionary parameters
	const int M					= 8;				//Number of populations in a species
	const int N 				= 20;				//Number of individuals per population
	const int N_best			= 10;				//Number of individuals in "hall of fame". Best individuals of all time
	const int geneLength		= 12;				//Number of bits per gene (The number of possible values for a parameter is 2^geneLength-1)
	//const int nGenes 			= 12;				//Set this as the number of parameters in your model
	//const int genomeLength 		= nGenes*geneLength;//Number of bits for all genes
	extern int nGenes 			;				//Set this as the number of parameters in your model
	extern int genomeLength 	;					//Number of bits for all genes
	const int generations 		= (int)1e4;			//Number of generations to run the simulation
	const double Tmin 			= 0.0001;			//Minimum temperature of the ladder. Preferrably between close to 1
	const double Tmax 			= 50;				//Maximum temperature of the ladder. Preferrably around H_max
	const int seed 				= 8;				//Seed for the random number generator
	const double log_param 		= 1.3;
	const double log_const 		= 1/log(log_param);
	
	//Probabilities genetic operators
	const double qmig = 0.005;						//The migration probability vs evolution probability
	const double qm = 0.3;							//The mutation probability vs crossover (1-qm for crossover)
	const double qma = 0.005;						//elite mutation rate (qm*(1-qma) for regular mutation)
	const int r_num = 50;
	
	//Probabilities for smart copy crossover
	const double P0 = 0.1;							//If parents have the  SAME bit, reverse with probability p0
	const double P1 = 0.4;							//If parents have DIFFERENT bit, reverse good bit with probability p1
	const double P2 = 0.7;							//If parents have DIFFERENT bit, reverse bad  bit with probability p2
	
	//Probability matrix for smart copy crossover (Liang and Wong, 2000)
	const double p_matrix[] = { P0*P0 + (1 - P0)*(1 - P0), 2 * P0*(1 - P0), P1*(1 - P2) + P2*(1 - P1), P1*P2 + (1 - P1)*(1 - P2) };
}
using namespace constants;
class boundaries {
public:
	boundaries(){
		all_bounds.resize(nGenes, 2);				
		upper_bound.resize(nGenes, 1);
		lower_bound.resize(nGenes, 1);
	}
	MatrixXd all_bounds;							//Lower and upper boundary of parameters
	ArrayXd	upper_bound;							//upper boundary of parameters
	ArrayXd	lower_bound;							//lower boundary of parameters
};
extern boundaries bounds;

#endif
