
#ifndef POPULATION_H   // if x.h hasn't been included yet...
#define POPULATION_H   //  #define this so the compiler knows it has been included
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>
#include "personality.hpp"
using namespace std;
using namespace constants;
using namespace Eigen;
class inData;		//Forward declaration
class personality;	//Forward declaration
//class DNA;			//Forward declaration

class paramLine{
	private:
		ArrayXd o;
		ArrayXd d;
	public:
		paramLine(): o(nGenes), d(nGenes){}
		//Always use "Through" first, to set "o" and "d".
		void Through(ArrayXd &p1, ArrayXd &p2){ //Points in parameter space p1 and p2
			for(int i = 0;i < nGenes; i++){
				o(i) = p1(i);
				d(i) = p2(i) - p1(i);
			}
		}
		ArrayXd pointAt(double r){
			ArrayXd v(nGenes);
			for(int i = 0;i < nGenes; i++){
				v(i) = o(i)+d(i)*r;
			}
			return v;
		}
		double distance(ArrayXd &p1, ArrayXd &p2){
			return sqrt((p2-p1).square().sum());
		}
		double line_max(ArrayXd &Bu ) { //Get the largest r within upper boundary Bu and lower boundary Bl
		return (Bu - o).cwiseQuotient(d).minCoeff();
		//return fmin((Bu - o).cwiseQuotient(d).minCoeff(),
		//			(Bl - o).cwiseQuotient(d).minCoeff());
	}
	double line_min(ArrayXd &Bl) { //Get the largest r within upper boundary Bu and lower boundary Bl
		return -(Bl - o).cwiseQuotient(-d).minCoeff();
	}
};

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
	//ParametrizedLine<double, Dynamic> line; //Line for doing the snooker crossover
	paramLine line;
	void copy(personality&, personality&);
	void copy(DNA&, DNA&);
	void getFitness(personality&, inData&);
	int operator()() { //Return the bit at a.
		return 0;
	}
	friend ostream &operator<<(std::ostream &os, population const &);


};

#endif
