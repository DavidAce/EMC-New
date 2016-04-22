#ifndef SPECIES_H   // if x.h hasn't been included yet...
#define SPECIES_H   //  #define this so the compiler knows it has been included
#include <iostream>
#include <vector>
#include "..\Eigen\Dense"
#include "..\Eigen\Core"
#include "personality.h"
#include "population.h"
#include "datafiles.h"

class counters {
public:
	int store_counter;
	int generation;
	double simulation_time;
	clock_t simulation_tic;
	clock_t simulation_toc;
	double evolution_time;
	clock_t evolution_tic;
	clock_t evolution_toc;
};



class species{
private:
	void zero_counters();
public:
	species(const int, const char **) ;
	inData in;				//Experimental data for minimization
	population pop[M];	//Array of separate populations to evolve independently
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
