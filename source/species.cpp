#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <time.h>
#include <Eigen/Dense>
#include <Eigen/Core>
#include "constants.hpp"
#include "mymath.hpp"
#include "minimization.hpp"
#include "DNA.hpp"
#include "datafiles.hpp"
#include "population.hpp"
#include "evolution.hpp"
#include "species.hpp"


using namespace std;
using namespace constants;
using namespace Eigen;

double species::champion_value() {
	int best_index = 0;
	double best_H = pop[0].bestguys[N_best - 1].H;
	for (int i = 1; i < M; i++) {
		if (pop[i].bestguys[N_best - 1].H < best_H) {
			best_index = i;
			best_H = pop[i].bestguys[N_best - 1].H;
		}
	}
	return pop[best_index].bestguys[N_best - 1].value;
}

double species::champion_fitness() {
	double best_H = pop[0].bestguys[N_best - 1].H;
	for (int i = 1; i < M; i++) {
		if (pop[i].bestguys[N_best - 1].H < best_H) {
			best_H = pop[i].bestguys[N_best - 1].H;
		}
	}
	return best_H;
}

int species::champion_number() {
	int best_index = 0;
	double best_H = pop[0].bestguys[N_best - 1].H;
	for (int i = 1; i < M; i++) {
		if (pop[i].bestguys[N_best - 1].H < best_H) {
			best_index = i;
			best_H = pop[i].bestguys[N_best - 1].H;
		}
	}
	return best_index;
}

void species::print_progress(){
	count.generation++;
	if (mod(count.generation, 1000) == 0) {
		cout << fixed << setprecision(9);
		cout << "\rGeneration... " << setw(7) << count.generation << " | Best Fitness H: ";
		for (int m = 0; m < M; m++) {
			cout << setw(10) << pop[m].bestguys[N_best - 1].H << " ";
			//cout << setw(9) << pop[m].guys[N - 1].H << " ";
		}
		cout << flush;
		/*cout << " | H best: " << setw(7) << pop.bestguys[N_best - 1].H;
		cout << " | H current: " << setw(7) << pop.guys[N - 1].H;
		cout << " | Value: " << setw(7) << pop.bestguys[N_best - 1].value;*/
	}
}

void species::wakeUpAll() {
	for (int i = 0; i < M; i++) {
		pop[i].wakeUpPop(in);
		pop[i].population_number = i;
	}
}
void species::zero_counters() {
	count.evolution_time = 0;
	count.simulation_time = 0;
	count.store_counter = 0;
	count.generation = 0;
}

void species::copy(personality &destination, personality &source) {
	destination.born = source.born;
	destination.H = source.H;
	destination.t = source.t;
	destination.genome.parameters = source.genome.parameters;
	destination.value = source.value;
	destination.genome.chromosomes = source.genome.chromosomes;
	// std::copy(destination.genome.chromosomes, destination.genome.chromosomes + nGenes, source.genome.chromosomes);
	//std::copy(destination.genome.chromosomes.begin(), destination.genome.chromosomes.end(), source.genome.chromosomes);

}


species::species(const int argc, const char **argv) : in(argc, argv)
{	
	zero_counters();
	wakeUpAll();
}
