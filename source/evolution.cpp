#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <Eigen/Dense>
#include <Eigen/Core>
#include "../Eigen/Geometry"
#include "constants.hpp"
#include "randomFunctions.hpp"
#include "mymath.hpp"
#include "minimization.hpp"
#include "DNA.hpp"
#include "personality.hpp"
#include "population.hpp"
#include "species.hpp"
#include "datafiles.hpp"
#include "evolution.hpp"

using namespace std;
using namespace Eigen;
using namespace constants;

inline void roulette_select(personality guys[], ArrayXi &selected, long double &Z, double s) {
	//Selected 0 is the guy with good fitness (lowest), selected 1 is random

	//double total_H = 0;
	ArrayXd roulette(N);
	double lucky_number;
	//Make a roulette wheel
	Z = 0;
	for (int i = 0; i < N; i++) {
		Z += static_cast<long double>(exp(-guys[i].H / s));
		roulette(i) = static_cast <double> (Z); //Cumulative sum - Low fitness gives large area on roulette
	}
	lucky_number = uniform_double(&rng, 0 , static_cast<double>(Z));
	for (int i = 0; i < N; i++) {
		if (lucky_number <= roulette(i)) {
			selected(0) = i;
			break;
		}
	}
	selected(1) = selected(0);
	while (selected(1) == selected(0)) {
		selected(1) = uniform_integer(&rng, 0, N - 1);
	}
}
inline void bitselector_smartCopy(population &pop, ArrayXi &selected, ArrayXd &n_ab_XY, ArrayXd &n_ab_YX) {
	//n_ab_XY(0): n_11 = guys: SAME | newguys: SAME
	//n_ab_XY(1): n_12 = guys: SAME | newguys: DIFF
	//n_ab_XY(2): n_21 = guys: DIFF | newguys: SAME
	//n_ab_XY(3): n_22 = guys: DIFF | newguys: DIFF

	//n_ab_YX(0): n_11 = newguys: SAME | guys: SAME
	//n_ab_YX(1): n_12 = newguys: SAME | guys: DIFF
	//n_ab_YX(2): n_21 = newguys: DIFF | guys: SAME
	//n_ab_YX(3): n_22 = newguys: DIFF | guys: DIFF
	
	//#pragma omp parallel for 
	for (int i = 0; i < genomeLength; i++) {
		if (pop.guys[selected(0)].genome(i) == pop.guys[selected(1)].genome(i)) {
			//If parents loci are the same: Copy the values
			pop.newguys[selected(0)].genome.copy_loci(i, pop.guys[selected(0)].genome(i));
			pop.newguys[selected(1)].genome.copy_loci(i, pop.guys[selected(1)].genome(i));

			//Independently reverse with probability p0
			for (int j = 0; j < 2; j++) {
				if (uniform_double(&rng, 0, 1) < P0) {
					pop.newguys[selected[j]].genome.flip_loci(i);
				}
			}

			if (pop.newguys[selected(0)].genome(i) == pop.newguys[selected(1)].genome(i)) {
				//#pragma omp atomic
				n_ab_XY(0)++; //Offsprings have identical bit
				//#pragma omp atomic
				n_ab_YX(0)++; //Parents have identical bit
			}
			else {
				//#pragma omp atomic
				n_ab_XY(1)++; //Offsprings have different bit
				//#pragma omp atomic
				n_ab_YX(2)++; //Parents have different bit
			}

		}
		else {
			//If parents loci are different: Copy the values
			pop.newguys[selected(0)].genome.copy_loci(i, pop.guys[selected(0)].genome(i));
			pop.newguys[selected(1)].genome.copy_loci(i, pop.guys[selected(1)].genome(i));
			
			//Reverse with probabilities p1 and p2 respectively
			if (uniform_double(&rng, 0, 1) < P1) {
				pop.newguys[selected(0)].genome.flip_loci(i);
			}
			if (uniform_double(&rng, 0, 1) < P2) {
				pop.newguys[selected(1)].genome.flip_loci(i);
			}

			if (pop.newguys[selected(0)].genome(i) == pop.newguys[selected(1)].genome(i)) {
				//#pragma omp atomic
				n_ab_XY(2)++; //Offsprings have identical bit
				//#pragma omp atomic
				n_ab_YX(1)++; //Offsprings have identical bit
			}
			else {
				//#pragma omp atomic
				n_ab_XY(3)++; //Offsprings have different bit
				//#pragma omp atomic
				n_ab_YX(3)++; //Parents have different bit

			}

		}
	}

}
void mutation(population &pop, inData& in) {
	int mutantGenes;	//Number of points to mutate
	int mutant;			//Which guy to mutate
	double dH;
	//#pragma omp parallel for private(mutantGenes, mutant,dH)
	for (int i = 0; i < N; i++) {
		mutantGenes =  uniform_integer(&rng, 1, genomeLength - 1);
		ArrayXi loci(mutantGenes);
		mutant = i;// uniform_integer(&rng, 0, N - 1);
		rndChoice(loci.data(), mutantGenes, genomeLength);				//Choose locus to mutate
		pop.newguys[mutant].genome.flip_loci(loci);	//Flip bits
		pop.getFitness(pop.newguys[mutant],in);							//Get Fitness score
		
		//Perform Metropolis
		dH = pop.newguys[mutant].H - pop.guys[mutant].H;		//used to decide if we accept the new guy or not.
		
		if (dH < 0 || exp(-dH / pop.newguys[mutant].t) > uniform_double(&rng, 0, 1)) {
			pop.newguys[mutant].born = pop.generation;
			pop.copy(pop.guys[mutant], pop.newguys[mutant]);

		}
		else { 	
			//Revert changes in newguys,  i.e sync them for the next round
			pop.copy(pop.newguys[mutant], pop.guys[mutant]);
		}
	}


}
void mutation_elite(population &pop, inData& in) {
	//int mutantGenes;	//Number of points to mutate
	int mutant;			//which guy to mutate
	int elite_mutant;	//Which elite guy to receive wisdom from
	double dH;
	//#pragma omp parallel for private( mutant, elite_mutant,dH)
	for (int i = 0; i < N; i++) {
		//mutantGenes = 1;// uniform_integer(&rng, 1, genomeLength - 1);
		int loci = uniform_integer(&rng, 1, genomeLength - 1);
		mutant = i; //uniform_integer(&rng, 0, N - 1);
		//Fill loci with mutantGenes genome points to be mutated
		//rndChoice(loci.data(), mutantGenes, genomeLength);	//Choose locus to mutate
		
		//Copy an elite guy to a new guy to be a guinnea pig
		elite_mutant = uniform_integer(&rng, 0, N_best - 1);
		pop.copy(pop.newguys[mutant].genome, pop.bestguys[elite_mutant].genome);	//Copy DNA only!
		pop.newguys[mutant].genome.flip_loci(loci);									//Flip bits	
		pop.getFitness(pop.newguys[mutant],in);									//Get Fitness score
		
		//Perform metropolis
		dH = pop.newguys[mutant].H - pop.guys[mutant].H;		//used to decide if we accept the new guy or not.
		if (dH < 0 || exp(-dH / pop.newguys[mutant].t) > uniform_double(&rng, 0, 1)) {
			pop.newguys[mutant].born = pop.generation;
			pop.copy(pop.guys[mutant], pop.newguys[mutant]);
		}
		else {
			//Revert changes in newguys of latest events... i.e sync them for the next round
			pop.copy(pop.newguys[mutant], pop.guys[mutant]);
		}

	}
}
void crossover(population &pop, inData& in) {
	int  matings;
	int nMatings = (int)(0.2*N);
	ArrayXi selected(2);
	ArrayXd expX(2);
	ArrayXd expY(2);
	int crossoverPoint;
	double s = 0.5;// pop.guys[(int)(0.9*N)].t;
	double rc, dHt0, dHt1;
	double PXX, PYY;			//Selection probabilities
	double PXY = 1, PYX = 1;	//Generating probabilities (PXY = PYX for this operator)
	double TXY, TYX;			//Transition probability

	long double ZX=0, ZY=0;		//Sum of Boltzmann-weights for current (X) and offspring(Y) populations

	for (matings = 0; matings < nMatings; matings++) {
		roulette_select(pop.guys, selected, ZX, s); //Selected 0 will be good, selected 1 random

		expX(0) = exp(-pop.guys[selected(0)].H / s); //good guy
		expX(1) = exp(-pop.guys[selected(1)].H / s); //bad guy
		PXX = static_cast<double>( 1 / ((N - 1)*ZX)*(expX(0) + expX(1))); //P((xi,xj) | x)
		
		//Now mate the newguys to create offsprings
		crossoverPoint = uniform_integer(&rng, 1, genomeLength-1);
		//cout << selected.transpose() << " " << pop.newguys[selected(0)].H << endl;
		for (int i = crossoverPoint; i < genomeLength; i++) {
			pop.newguys[selected(0)].genome.copy_loci(i, pop.guys[selected(1)].genome(i));
			pop.newguys[selected(1)].genome.copy_loci(i, pop.guys[selected(0)].genome(i));
		}
		for (int i = 0; i < 2; i++) {
			pop.getFitness(pop.newguys[selected(i)],in);
		}

		expY(0) = exp(-pop.newguys[selected(0)].H / s);
		expY(1) = exp(-pop.newguys[selected(1)].H / s);
		for (int i = 0; i < N; i++) {
			ZY += exp(-pop.newguys[i].H / s);
		}
		PYY = static_cast<double>(1 / ((N - 1)*ZY)*(expY(0) + expY(1))); //P((yi,yj) | y) selection probability
		dHt0 = (pop.newguys[selected(0)].H - pop.guys[selected(0)].H) / pop.guys[selected(0)].t; //good 
		dHt1 = (pop.newguys[selected(1)].H - pop.guys[selected(1)].H) / pop.guys[selected(1)].t; //bad 

		TXY = PXX*PXY;
		TYX = PYY*PYX;

		rc = exp(-dHt0 - dHt1)*TXY / TYX;
		//cout << "Selected:	" << selected.transpose()<< "	P: " << rc << " dHt0: " << dHt0 << " dHt1: " << dHt1 << endl;
		//Accept or reject
		if (uniform_double(&rng, 0, 1) < fmin(1, rc)) {
			pop.newguys[selected(0)].born = pop.generation;
			pop.newguys[selected(1)].born = pop.generation;
			pop.copy(pop.guys[selected(0)], pop.newguys[selected(0)]);
			pop.copy(pop.guys[selected(1)], pop.newguys[selected(1)]);
		}
		else {
			//Revert changes on newguys
			pop.copy(pop.newguys[selected(0)], pop.guys[selected(0)]);
			pop.copy(pop.newguys[selected(1)], pop.guys[selected(1)]);
		}
		//Make sure newguys are up to speed on newest events

		

		
	/*	cout << "Parameters 0:	" << pop.guys[selected(0)].genome.parameters.transpose() << endl;
		cout << "Parameters 1:	" << pop.guys[selected(1)].genome.parameters.transpose() << endl;
		getchar();*/
	}
}
void crossover_elite(population &pop, inData& in) {
	//Start roulette selection
	int matings;
	int nMatings = (int)(0.2*N);
	ArrayXi selected(2);
	int crossoverPoint;
	double rc, dHt0, dHt1;
	double PXX, PYY;			//Selection probabilities
	double PXY = 1, PYX = 1;	//Generating probabilities (PXY = PYX for this operator)
	double TXY, TYX;			//Transition probability
	ArrayXd expX(2);
	ArrayXd expY(2);
	long double ZX = 0, ZY = 0;
	double s = 0.5;// pop.guys[(int)(0.9*N)].t;
	int random_bestguy; //Guy to inject into selected[1]
	for (matings = 0; matings < nMatings; matings++) {
		//Selected 0 is a boltzmann "just good" guy. Selected 1 is a random any guy.
		roulette_select(pop.guys, selected, ZX, s);
		//Let the newguy selected[1] impersonate a random bestguy but keep the temperature. Then newguy gets superpowers
		random_bestguy = uniform_integer(&rng, 0, N_best-1);
		pop.copy(pop.newguys[selected(1)], pop.bestguys[random_bestguy]);
		pop.newguys[selected(1)].born = pop.guys[selected(1)].born;		//Keep date of birth count
		pop.newguys[selected(1)].t = pop.guys[selected(1)].t;						//Keep temperature

		//Now selected(1) is some guy high on the temperature ladder with amazing bestguy-genes and fitness
		expX(0) = exp(-pop.newguys[selected(0)].H / s); //good guy
		expX(1) = exp(-pop.newguys[selected(1)].H / s); //gooder guy
		PXX = static_cast<double>(1 / ((N - 1)*ZX)*(expX(0) + expX(1))); //P((xi,xj) | x)
		
		//Now mate the to create offspring. Let a bestGuy inject DNA in this process!
		crossoverPoint = uniform_integer(&rng, 1, genomeLength - 1);
		for (int i = 0; i < genomeLength; i++) {
			if (i < crossoverPoint) {
				pop.newguys[selected(0)].genome.copy_loci(i, pop.guys[selected(0)].genome(i));
				pop.newguys[selected(1)].genome.copy_loci(i, pop.bestguys[random_bestguy].genome(i));
			}
			else {
				pop.newguys[selected(0)].genome.copy_loci(i, pop.bestguys[random_bestguy].genome(i));
				pop.newguys[selected(1)].genome.copy_loci(i, pop.guys[selected(0)].genome(i));
			}
		}

		//From now on the newguys are offsprings. The parents are a good guy selected 1 and a bestguy random_bestguy. Adoptive father is guy selected 1
		
		for (int i = 0; i < 2; i++) {
			pop.getFitness(pop.newguys[selected(i)],in);
		}

		expY(0) = exp(-pop.newguys[selected(0)].H / s);
		expY(1) = exp(-pop.newguys[selected(1)].H / s);
		ZY = 0;
		for (int i = 0; i < N; i++) {
			ZY += exp(-pop.newguys[i].H / s);
		}
		PYY = static_cast<double>(1 / ((N - 1)*ZY)*(expY(0) + expY(1))); //P((xi,xj) | x)
		TXY = PXX*PXY;
		TYX = PYY*PYX;

		dHt0 = (pop.newguys[selected(0)].H - pop.guys[selected(0)].H) / pop.guys[selected(0)].t; //good 
		dHt1 = (pop.newguys[selected(1)].H - pop.guys[selected(1)].H) / pop.guys[selected(1)].t; //bad 

		//cout << fixed << setprecision(6) << "guy[" << selected[0] << "].H = " << pop.guys[selected(0)].H;
		//cout << fixed << setprecision(6) << "	guy[" << selected[1] << "].H = " << pop.guys[selected(1)].H << endl;
		//cout << fixed << setprecision(10) << expX.transpose() << " " << expY.transpose() << endl;
		//cout << fixed << setprecision(10) << "ZX: " << ZX << " ZY: " << ZY << endl;
		//cout << fixed << setprecision(6) << "PXX: " << PXX << " PXY: " << PXY << endl;
		//cout << fixed << setprecision(6) << "PYY: " << PYY << " PYX: " << PYX << endl;
		//cout << p_matrix[0] << " " << p_matrix[1] << " " << p_matrix[2] << " " << p_matrix[3] << endl << endl;



		rc = exp(-dHt0 - dHt1)*TXY / TYX;

		//Accept or reject
		if (uniform_double(&rng, 0, 1) < fmin(1, rc)) {
			//Keep children.
			pop.newguys[selected(0)].born = pop.generation;
			pop.newguys[selected(1)].born = pop.generation;
			pop.copy(pop.guys[selected[0]], pop.newguys[selected[0]]);
			pop.copy(pop.guys[selected[1]], pop.newguys[selected[1]]);
		}
		else {
			//Keep parents.
			pop.copy(pop.newguys[selected[0]], pop.guys[selected[0]]);
			pop.copy(pop.newguys[selected[1]], pop.guys[selected[1]]);
		}

	}

}
void crossover_smartCopy(population &pop, inData& in) {
	//Start roulette selection																							
	int matings;
	int nMatings = (int)(0.2*N);
	ArrayXi selected(2);
	double s = 0.5;// pop.guys[(int)(0.9*N)].t;
	double rc, dHt0, dHt1;
	double PXX, PYY; //Selection probabilities
	double PXY, PYX; //Generating probabilities (PXY != PYX for this operator)
	double TXY, TYX; //Transition probability
	ArrayXd expX(2);
	ArrayXd expY(2);
	ArrayXd n_ab_XY(4); //Exponents for smartCopy probabilities. Note that n_ab.sum() = genomeLength
	ArrayXd n_ab_YX(4); //Exponents for smartCopy probabilities. Note that n_ab.sum() = genomeLength
	long double ZX = 0, ZY = 0;
	for (matings = 0; matings < nMatings; matings++) {
		n_ab_XY.fill(0);
		n_ab_YX.fill(0);
		roulette_select(pop.guys, selected, ZX, s); //Selected 0 will be good, selected 1 random
		expX(0) = exp(-pop.guys[selected(0)].H / s); //good guy
		expX(1) = exp(-pop.guys[selected(1)].H / s); //bad guy
		PXX = static_cast<double>(1 / ((N - 1)*ZX)*(expX(0) + expX(1))); //P((xi,xj) | x)

		//Now mate the newguys to create offsprings		
		bitselector_smartCopy(pop, selected, n_ab_XY,n_ab_YX);
		
		//Now we need generating probabilities PXY and PYX
		PXY = 1;
		PYX = 1;
		for (int i = 0; i < 4; i++) {
			PXY *= pow(p_matrix[i], n_ab_XY(i)/genomeLength); //Convert exponent to frequency
			PYX *= pow(p_matrix[i], n_ab_YX(i)/genomeLength); //Convert exponent to frequency
		}
		for (int i = 0; i < 2; i++) {
			pop.getFitness(pop.newguys[selected(i)],in);
		}
		expY(0) = exp(-pop.newguys[selected(0)].H / s);
		expY(1) = exp(-pop.newguys[selected(1)].H / s);
		for (int i = 0; i < N; i++) {
			ZY += exp(-pop.newguys[i].H / s);
		}
		PYY = static_cast<double>(1 / ((N - 1)*ZY)*(expY(0) + expY(1))); //P((yi,yj) | y) selection probability
	
		dHt0 = (pop.newguys[selected(0)].H - pop.guys[selected(0)].H) / pop.guys[selected(0)].t; //good 
		dHt1 = (pop.newguys[selected(1)].H - pop.guys[selected(1)].H) / pop.guys[selected(1)].t; //bad 
		
		TXY = PXX*PXY;
		TYX = PYY*PYX;
		//cout << fixed << setprecision(6) << "guy[" << selected[0] << "].H = " << pop.guys[selected(0)].H;
		//cout << fixed << setprecision(6) << "	guy[" << selected[1] << "].H = " << pop.guys[selected(1)].H <<  endl;
		//cout << fixed << setprecision(10) << expX.transpose() << " " << expY.transpose() << endl;
		//cout << fixed << setprecision(10) << "ZX: " << ZX << " ZY: " << ZY << endl;
		//cout << "n_ab_XY: " << n_ab_XY.transpose() << endl;
		//cout << "n_ab_YX: " << n_ab_YX.transpose() << endl;
		//cout << fixed << setprecision(6) << "PXX: " << PXX << " PXY: " << PXY << endl;
		//cout << fixed << setprecision(6) << "PYY: " << PYY << " PYX: " << PYX << endl ;
		//cout << p_matrix[0] << " " << p_matrix[1] << " " << p_matrix[2] << " " << p_matrix[3] << endl<< endl;
		rc = exp(-dHt0 - dHt1)*TXY / TYX;
		//Accept or reject
		if (uniform_double(&rng, 0, 1) < fmin(1, rc)) {
			//Keep Children
			pop.newguys[selected(0)].born = pop.generation;
			pop.newguys[selected(1)].born = pop.generation;
			pop.copy(pop.guys[selected[0]], pop.newguys[selected[0]]);
			pop.copy(pop.guys[selected[1]], pop.newguys[selected[1]]);

		}
		else {
			//Keep parents
			pop.copy(pop.newguys[selected[0]], pop.guys[selected[0]]);
			pop.copy(pop.newguys[selected[1]], pop.guys[selected[1]]);
		}

	}
}
void crossover_snooker(population &pop, inData& in) {
	//Perhaps the distribution f(r) can be done with a roulette?
	//f(r) = exp(H(x+re))/integral_-r_min^r_max exp(H(x+re)) typ?

	//Notation from snooker crossover 
	//Choose 3 guys from pop
	
	Array2i selected; //Guy 0 is "x_i", 1 is x_j, "anchor";
	rndChoice(selected.data(), 2, N - 1);
	// ArrayXd &p0 = pop.guys[selected(0)].genome.parameters; //Reference to guy 0, 
	// ArrayXd &p1 = pop.guys[selected(1)].genome.parameters; //Reference to guy 1, who is lower on the ladder (better)
	pop.line.Through(pop.guys[selected(0)].genome.parameters,pop.guys[selected(1)].genome.parameters);
	//Distance in terms of r between guys
	//double distance = pop.line.distance(pop.guys[selected(1)].genome.parameters,pop.guys[selected(0)].genome.parameters);
	//Vary r to find the walls of the parameter domain
    double r_max = pop.line.line_max(bounds.upper_bound);
    double r_min = pop.line.line_min(bounds.lower_bound);
	ArrayXd r_point(r_num);	
	for (int i = 0; i < r_num; i++) {
		//r_point(i) = uniform_double(&rng, fmin(r_min,r_max), fmax(r_min,r_max));
		r_point(i) = gaussian_truncated(&rng,		//Create gaussians points centerd around the good performer
										fmin(r_min, r_max),
										fmax(r_min, r_max),
										1.0,
										0.3);
		pop.snookerGuys[i].genome.set_parameters(pop.line.pointAt(r_point(i)));
		pop.snookerGuys[i].H = fitnessTest(pop.snookerGuys[i], in);
		pop.snookerGuys[i].t = pop.guys[selected(0)].t;
	}
	//Time to make a roulette to see which "r" is chosenpop.

	ArrayXd roulette(r_num);
	double lucky_number = uniform_double(&rng, 0, 1);
	//Make a roulette wheel
	double Z = 0;
	for (int i = 0; i < r_num; i++) {
		Z += exp(-pop.snookerGuys[i].H); //Area on roulette wheel proportional to Boltzmann weights
		roulette(i) = Z; //Cumulative sum - Low fitness gives large area on roulette
	}

	roulette = roulette.array() / Z; //Normalize to a proper distribution
	int chosen = -1;
	for (int i = 0; i < r_num; i++) {
		if (lucky_number <= roulette(i)) {
			chosen = i;
			pop.guys[selected(0)].genome.set_parameters(pop.snookerGuys[i].genome.parameters); //Copy his shit
			pop.guys[selected(0)].H 	= pop.snookerGuys[i].H;
			pop.guys[selected(0)].value = pop.snookerGuys[i].value;
			pop.guys[selected(0)].born 	= pop.generation;
			break;
		}
	}

}
void exchange(population &pop) {
	int i, j, ex; //ex is the exchange iteration
	double re, dt_inv, dH;
	for (ex = 0; ex < N; ex++) {
		i = uniform_integer(&rng, 0, N - 1);										//Chose a guy...
		j = i == 0 ? 1 : (i== N-1? N-2: i + 2 * uniform_integer(&rng, 0, 1) - 1);	//And a neighbor
		dH = pop.guys[i].H - pop.guys[j].H;
		dt_inv = 1 / pop.guys[i].t - 1 / pop.guys[j].t;
		re = exp(dH*dt_inv);
		
		if (uniform_double(&rng, 0, 1) < fmin(1, re)) {
			//Save the date of birth
			pop.newguys[j].born = pop.guys[i].born;
			pop.newguys[i].born = pop.guys[j].born;
			//Swap the guys on the ladder...
			pop.copy(pop.guys[i], pop.newguys[j]);
			pop.copy(pop.guys[j], pop.newguys[i]);
			//But keep the temperature of position i as t_i,and j as t_j
			pop.guys[i].t = pop.newguys[i].t;
			pop.guys[j].t = pop.newguys[j].t;
			//Now the list of newguys should be updated to be identical to guys
			pop.copy(pop.newguys[i], pop.guys[i]);
			pop.copy(pop.newguys[j], pop.guys[j]);
		}
	}
}
void migration(species &sp) {
	//First select randomly m out of M pobulations, subject to migration
	//Then select randomly  n (n=1?) out of N guys to migrate
	//Accept the new population based on the new population fitness score
	int n;										//Guy chosen to migrate
	double dH;									//Probability of migration

	int m = uniform_integer(&rng, 2, M);		//How many populations to migrate
	ArrayXi s_pop(m);							//Array of sender populations
	ArrayXi r_pop(m);							//Array of receiving populations
	rndChoice(s_pop.data(), m, M);				//Fill the array of migrators with m random populations
	
	for (int i = 0; i < m; i++) {
		r_pop(i) = i < m - 1 ? s_pop(i + 1) : s_pop(0) ;						//Define receiving populations cyclically
		n = uniform_integer(&rng, 0, N-1);									//Who will be the lucky one?
		sp.copy(sp.pop[r_pop(i)].newguys[n], sp.pop[s_pop(i)].newguys[n]);	//Inject a guy into another population
		dH = sp.pop[r_pop(i)].newguys[n].H - sp.pop[r_pop(i)].guys[n].H;
		if (dH < 0 || exp(-dH / sp.pop[r_pop(i)].newguys[n].t) > uniform_double(&rng, 0, 1)) {
			sp.pop[r_pop(i)].newguys[n].born = sp.count.generation;
			sp.copy(sp.pop[r_pop(i)].guys[n], sp.pop[r_pop(i)].newguys[n]);
		}
		else {
			//Revert changes in newguys,  i.e sync them for the next round
			sp.copy(sp.pop[r_pop(i)].newguys[n], sp.pop[r_pop(i)].guys[n]);
		}
	}
	
}
void insertguy(population &pop, int from, int to) {
	int i;
	//Push out the worst bestguy and move everybody up until "to"
	for (i = 0; i < to; i++) {
		pop.copy(pop.bestguys[i], pop.bestguys[i + 1]);
	}
	//Insert the exceptional guy into the illustrious group of best guys
	pop.copy(pop.bestguys[to], pop.guys[from]);
}
void find_elite(population &pop) {
	int i;
	int count = 0;
	//int best[N_best];
	int tried_guys[N_best];
	double lowest_H;
	int lowest_i = N;
	int j = 0;
	int success = 1;
	while (success == 1) {
		lowest_H = 1e10;
		//Find the best guy yet among guys
		for (i = 0; i < N; i++) {
			//Check if i is in skip-list
			if (isvalueinarray(i, tried_guys, j) == 1) { continue; }
			if (pop.guys[i].H < lowest_H) {
				lowest_H = pop.guys[i].H;
				lowest_i = i;
			}
		}
		count++;
		//By now we should have a winner, check if he's better than any elite
		//We have already tried to match j guys, so we can start from j?
		for (i = 0; i < N_best; i++) {
			success = 0;
			if (lowest_H == pop.bestguys[N_best - i - 1].H) {
				success = 1;
				tried_guys[j] = lowest_i;
				j++;
				return;//The guy is already here! Do not copy, just check next guy
			}
			if (lowest_H < pop.bestguys[N_best - i - 1].H) {
				insertguy(pop, lowest_i, N_best - i - 1);
				success = 1;
				tried_guys[j] = lowest_i;
				j++;
				break;
			}
		}

	}
}
void evolve(population &pop, inData& in) {
	//This function selects either mutation type operators or crossover type. 
	//Then does an exchange operator, and finally finds the best guys in the population
	pop.generation++;
	int dice;
	//Select mutation or crossover
	if (uniform_double(&rng, 0, 1) < qm) {
		if (uniform_double(&rng, 0, 1) < qma) {
			mutation_elite(pop, in);
		}
		else {
			mutation(pop, in);
		}
	}
	else {
		dice = uniform_integer(&rng, 0, 4);
		if 		(dice == 0)	 { crossover(pop, in); }
		else if (dice == 1)	 { crossover_smartCopy(pop, in); }
		else if (dice == 2)	 { crossover_elite(pop, in); }
		else				 { crossover_snooker(pop, in); }
		
		
	}

	exchange(pop);
	find_elite(pop);
	
}
