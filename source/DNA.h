
#ifndef DNA_H   // if x.h hasn't been included yet...
#define DNA_H   //  #define this so the compiler knows it has been included
#include "../Eigen/Dense"
#include "../Eigen/Core"
#include <bitset>
#include <iostream>


using namespace std;
using namespace constants;
using namespace Eigen;
class DNA {
private:
	double bin2dec(const int);
	bitset<geneLength> dec2bin(const int);
public:
	DNA ();
	DNA(bool);
	bitset<geneLength> chromosomes[nGenes]; //Binary representation
	ArrayXd parameters;						//Decimal representation

	bool operator == (DNA&);
	int operator()(int); 
	friend ostream &operator<<(std::ostream &os, DNA const &);
	int length = geneLength*nGenes;
	void flip_loci(const int);
	void flip_loci(Ref<ArrayXi>);
	void copy_loci(const int, const int);
	void set_parameter(const int,const double);
	void set_parameters(Ref<ArrayXd>);
	void update_parameters();
	

}; 

#endif
