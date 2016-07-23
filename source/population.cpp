#include <iostream>
#include <string>
#include <fstream>
#include <time.h>
#include <Eigen/Dense>
#include <Eigen/Core>
#include "constants.hpp"
#include "mymath.hpp"
#include "datafiles.hpp"
#include "minimization.hpp"
#include "DNA.hpp"
#include "evolution.hpp"
#include "population.hpp"
#include "species.hpp"



using namespace std;
using namespace constants;
using namespace Eigen;


std::ostream &operator<<(std::ostream &os, population const &pop) {
	for (int i = 0; i < N; i++) {
		os << pop.guys[i].H << endl;
	}
	return os;
}


void population::wakeUpGuys() {
	ArrayXd T(N);
	//Initialize some temperature ladder, here logarithmic.
	logspace( T,Tmax, Tmin, N);

	//linspace(T, Tmax, Tmin, N);
	//Initialize the genome sequence for all guys randomly
	for (int i = 0; i < N; i++) {
		guys[i].t = T(i); 		//Assign temperatures produced above
	}
}
void population::wakeUpBest() {
	int i;
	int j = 0;
	int copied_guys[N_best];
	double lowest_H;
	int lowest_i = 0; //The position i on the temperature ladder of the guy with lowest H
	while (j < N_best) {
		lowest_H = 1e10;
		//Find the best guy yet among guys
		for (i = 0; i < N; i++) {
			//Check if i is in skip-list
			if (isvalueinarray(i, copied_guys, N_best) == 1) { continue; }
			if (guys[i].H < lowest_H) {
				lowest_H = guys[i].H;
				lowest_i = i;
			}
		}
		//By now we should have a winner, copy him to bestguys and add him to skiplist
		copy(bestguys[N_best - j - 1], guys[lowest_i]);
		copied_guys[j] = lowest_i;
		j++;

	}
}


void population::getFitness4All(inData &in) {
	for (int i = 0; i < N; i++) {
		guys[i].H = fitnessTest(guys[i], in);
	}
}

void population::getFitness(personality& guy, inData &in) {
	guy.genome.update_parameters();
	guy.H = fitnessTest(guy, in);
}



void population::copy(personality &destination, personality &source) {
	//destination.born				= source.born;
	destination.H					= source.H;
	destination.t					= source.t;
	destination.genome.parameters	= source.genome.parameters;
	destination.value				= source.value;
	destination.genome.chromosomes  = source.genome.chromosomes;

	// std::copy(destination.genome.chromosomes, destination.genome.chromosomes + nGenes, source.genome.chromosomes);
}

void population::copy(DNA &destination, DNA &source) {
	destination.chromosomes  = source.chromosomes;

	// std::copy(destination.chromosomes, destination.chromosomes + nGenes, source.chromosomes);
}

void population::wakeUpPop(inData &in) {
	wakeUpGuys();
	getFitness4All(in);
	for (int i = 0; i < N; i++) {
		copy(newguys[i], guys[i]);
	}
	wakeUpBest();
}

