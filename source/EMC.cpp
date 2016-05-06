// PID.cpp : Defines the entry point for the console application.
//
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../Eigen/Dense"
#include "../Eigen/Core"
#include <string.h>
#include <omp.h>
#include <random>
#include <iostream>
#include <string>
#include <fstream>
#include <bitset>
#include <chrono>
#include "constants.hpp"
#include "randomFunctions.hpp"
#include "mymath.hpp"
#include "DNA.hpp"
#include "personality.hpp"
#include "population.hpp"
#include "species.hpp"
#include "datafiles.hpp"
#include "evolution.hpp"
using namespace Eigen;
using namespace std;
using namespace constants;
using namespace std::chrono;

RNGType rng;
boundaries bounds;
int main(int argc, const char **argv) {
	//Start up some files and folder for saving out
	Eigen::initParallel();
	omp_set_num_threads(8);
	species sp(argc,argv);
	rng.seed(6);
	//Start algorithm
	sp.count.simulation_tic = high_resolution_clock::now();
	// sp.count.evolution_tic = clock();				//Start timer
	#pragma omp parallel
	while (sp.count.generation < generations &&  sp.champion_fitness() > lowest_H) {
		#pragma omp single nowait
		{
		sp.print_progress();
		if (uniform_double(&rng, 0, 1) < qmig) {
			migration(sp);
		}
		}
		#pragma omp for nowait
		for (int i = 0; i < M; i++) {
			evolve(sp.pop[i], sp.in); 			//Evolve all the populations
			//cout << sp.pop[0] <<endl;

		}
		

		
	}

	//sp.count.evolution_time += (double)(sp.count.evolution_toc - sp.count.evolution_tic) / CLOCKS_PER_SEC; //Collect computation time

	sp.out.print_to_file(sp);
	//Print final parameters
	cout << endl << "Best Parameters: " 
		 << sp.pop[sp.champion_number()].bestguys[N_best - 1].genome.parameters.transpose() << endl;
	//Print timing to console
	sp.count.simulation_toc = high_resolution_clock::now();
	printf("\nTotal time:		%.3f seconds\n", std::chrono::duration<double>(sp.count.simulation_toc - sp.count.simulation_tic).count());
	return 0;
}