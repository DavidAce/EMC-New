#include <bitset>
#include <iostream>
#include <random>
#include "..\Eigen\Dense"
#include "..\Eigen\Core"
#include "constants.h"
#include "randomFunctions.h"
#include "DNA.h"



using namespace std;
using namespace constants;
using namespace Eigen;


std::ostream &operator<<(std::ostream &os, DNA const &genome) {
	for (int i = 0; i < nGenes; i++) {
		os << genome.chromosomes[i] << endl;
	}
	return os;
}

bool DNA::operator==(DNA &target) { //Compare DNA and return true if equal
	bool isequal = true;
	for (int i = 0; i < nGenes; i++) {
		isequal = isequal && (chromosomes[i] == target.chromosomes[i]);
	}
	return isequal;
}

void DNA::flip_loci(const int a) {
	//(loci)
	int gene = a / geneLength;
	int locus = a%geneLength;
	chromosomes[gene].flip(geneLength - locus - 1);
}

void DNA::flip_loci(Ref<ArrayXi> loci) {
	//Flip all bits in array loci
	int gene;// = a / geneLength;
	int locus;// = a%geneLength;
	for (int i = 0; i < loci.size(); i++) {
		gene = loci(i)/geneLength;
		locus = loci(i) % geneLength;
		chromosomes[gene].flip(geneLength - locus - 1);
	}
}


void DNA::copy_loci(const int a, const int bit) {
	//(loci,bit)
	bool b = bit == 1;
	int gene = a / geneLength;
	int locus = a%geneLength;
	chromosomes[gene].set(geneLength - locus - 1, b);
}

double DNA::bin2dec(const int i) {
	auto value = chromosomes[i].to_ulong();
	return value / (pow(2.0, geneLength) - 1) * (bounds.upper_bound(i) - bounds.lower_bound(i)) + bounds.lower_bound(i);
}

bitset<geneLength> DNA::dec2bin(const int i) {
	double value = parameters(i);
	bitset <geneLength> A =  ((value - bounds.lower_bound(i)) / (bounds.upper_bound(i) - bounds.lower_bound(i))*(pow(2.0, geneLength) - 1));
	return A;
}

void DNA::set_parameter(const int i, const double p) {
	//Set one parameter with a double
	parameters(i) = p;
	chromosomes[i] = dec2bin(i);
}

void DNA::update_parameters() {
	for (int i = 0; i < nGenes; i++) {
		parameters(i) = bin2dec(i);
	}
}

void DNA::set_parameters(Ref<ArrayXd> p) {
	//Set all parameters at once with an ArrayXd
	parameters = p;
	for (int i = 0; i < nGenes; i++) {
		chromosomes[i] = dec2bin(i);
	}
}
DNA::DNA(bool) {
	parameters.resize(nGenes);
}

DNA::DNA() {
	parameters.resize(nGenes);
	for (int i = 0; i < nGenes; i++) {
		parameters(i) = uniform_double(&rng, bounds.lower_bound(i),bounds.upper_bound(i));
		chromosomes[i] = dec2bin(i);
	}
}