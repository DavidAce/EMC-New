
#ifndef EVOLUTION_H
#define EVOLUTION_H  

#include "datafiles.hpp"
#include "personality.hpp"
//#include "population.h"

class population;	//Forward declaration
class personality;	//Forward declaration
class inData;		//Forward declaration
extern void roulette_select(personality [], int [], double *, double );
extern void bitselector_smartCopy(population &, int [], int []);
extern void mutation(population &, inData&);
extern void mutation_elite(population &, inData&);
extern void crossover(population &, inData&);
extern void crossover_elite(population &, inData&);
extern void crossover_smartCopy(population &, inData&);
extern void crossover_snooker(population &, inData&);
extern void exchange(population &);
extern void migration(species &);
extern void insertguy(population &, int , int );
extern void find_elite(population &);
extern void evolve(population &, inData &);

#endif