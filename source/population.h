
#ifndef POPULATION_H   // if x.h hasn't been included yet...
#define POPULATION_H   //  #define this so the compiler knows it has been included
#include <iostream>
#include "..\Eigen\Dense"
#include "..\Eigen\Core"
#include "personality.h"
using namespace std;
using namespace constants;
using namespace Eigen;
class inData;		//Forward declaration
class personality;	//Forward declaration
//class DNA;			//Forward declaration

class population{
private:
	void wakeUpGuys();
	void wakeUpBest();
	void wakeUpSnookerGuys();
	void getFitness4All(inData &);
public:
	population() :generation(0) {};
	void wakeUpPop (inData &);
	personality guys[N]; //Make an array of N guys
	personality newguys[N]; //Make a temporary array of N guinneapigs
	personality bestguys[N_best]; //Make an array of N/10 good performers
	personality snookerGuys[r_num];// (r_num, personality(true));//Make an array of r_num snooker guys
	int generation;				  //Number of generations for this population
	int population_number;		  //Which number this instance is in a species

	void copy(personality&, personality&);
	void copy(DNA&, DNA&);
	void getFitness(personality&, inData&);
	int operator()() { //Return the bit at a.
		return 0;
	}
	friend ostream &operator<<(std::ostream &os, population const &);


};

#endif
