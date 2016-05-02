// PID.cpp : Defines the entry point for the console application.
//
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "../Eigen/Dense"
#include "../Eigen/Core"
#include <string.h>
#include <omp.h>
#include <random>
#include <iostream>
#include <string>
#include <fstream>
#include <bitset>

#include "constants.h"
#include "randomFunctions.h"
#include "mymath.h"
#include "DNA.h"
#include "personality.h"
#include "population.h"
#include "species.h"
#include "datafiles.h"
#include "evolution.h"
using namespace Eigen;
using namespace std;
using namespace constants;


RNGType rng;
boundaries bounds;

int main(int argc, const char **argv) {
	//Start up some files and folder for saving out
	Eigen::initParallel();
	omp_set_num_threads(8);
	species sp(argc,argv);
	rng.seed(6);
	//Start algorithm
	sp.count.evolution_tic = clock();				//Start timer
	while (sp.count.generation < generations &&  sp.champion_fitness() > 0.001) {
		sp.print_progress();
		if (uniform_double(&rng, 0, 1) < qmig) {
			migration(sp);
		}
		else {
			#pragma omp parallel for
			for (int i = 0; i < M; i++) {
				evolve(sp.pop[i], sp.in); 			//Evolve all the populations
			}
		}
	}
	
	sp.count.evolution_toc = clock();				//End timer
	sp.count.evolution_time += (double)(sp.count.evolution_toc - sp.count.evolution_tic) / CLOCKS_PER_SEC; //Collect computation time

	sp.out.print_to_file(sp);
	//Print final parameters
	cout << endl << "Best Parameters: " 
		 << sp.pop[sp.champion_number()].bestguys[N_best - 1].genome.parameters.transpose() << endl;
	//Print timing to console
	printf("\nTotal time:		%.3f seconds\n", sp.count.evolution_time);
	return 0;
}