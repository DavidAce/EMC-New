#ifndef SPECIES_H   // if x.h hasn't been included yet...
#define SPECIES_H   //  #define this so the compiler knows it has been included
#include <iostream>
#include <vector>
#include <chrono>
#include <Eigen/Dense>
#include <Eigen/Core>
#include "personality.hpp"
#include "population.hpp"
#include "datafiles.hpp"
using namespace std::chrono;

typedef std::chrono::high_resolution_clock Clock;
class counters {
public:
	int store_counter;
	int generation;
	double simulation_time;
	high_resolution_clock::time_point simulation_tic;
	high_resolution_clock::time_point simulation_toc;
	double evolution_time;
	high_resolution_clock::time_point evolution_tic;
	high_resolution_clock::time_point evolution_toc;
};



class species{
private:
	void zero_counters();
public:
	species(const int, const char **) ;
	inData in;				//Experimental data for minimization
	population pop[M];		//Array of separate populations to evolve independently
	outData out;			//Experimental data for minimization
	counters count;
	double champion_fitness();
	double champion_value();
	int champion_number();
	void print_progress();
	void copy(personality &, personality &);
	void wakeUpAll();
};


#endif
