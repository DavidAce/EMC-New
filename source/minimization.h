

#ifndef MINIMIZATION_H   // if x.h hasn't been included yet...
#define MINIMIZATION_H   //  #define this so the compiler knows it has been included
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "..\Eigen\Dense"
#include "..\Eigen\Core"
#include "constants.h"
#include "personality.h"
#include "datafiles.h"

using namespace std;
using namespace constants;
using namespace Eigen;

class counters;
class personality;
class inData;

inline void map(ArrayXd &parameter, MatrixXd &xy, MatrixXd &xy_map, int i) {

	//Map the point (x,y) with the parameters
	xy_map(i,0) = parameter(0)
				+ parameter(1) * xy(i,0)
				+ parameter(2) * xy(i,1)
				+ parameter(3) * xy(i,0)* xy(i,0)
				+ parameter(4) * xy(i,1)* xy(i,1)
				+ parameter(5) * xy(i,0)* xy(i,1);
	xy_map(i,1) = parameter(6)
				+ parameter(7) * xy(i,0)
				+ parameter(8) * xy(i,1)
				+ parameter(9) * xy(i,0)* xy(i,0)
				+ parameter(10) * xy(i,1)* xy(i,1)
				+ parameter(11) * xy(i,0)* xy(i,1);
	//xy_map(i, 0) = parameter(0)
	//	+ parameter(1) * xy(i, 0)
	//	+ parameter(2) * xy(i, 1);
	//xy_map(i, 1) = parameter(3)
	//	+ parameter(4) * xy(i, 0)
	//	+ parameter(5) * xy(i, 1);
}

inline double fitnessTest(personality &guy, inData &in) {
	//data[0] has "old data" training-data before mapping
	//data[1] has "new data",training-data after mapping
	double H = 0;
	MatrixXd xy_map(in.data[0].rows(), in.data[0].cols());
	//This part measures the performance, or "fitness" H of the current parameters
	for (int i = 0; i < in.data[0].rows(); i++) {
		map(guy.genome.parameters, in.data[0], xy_map, i);											//Apply the mapping, 
	}
	//H = (xy_map - in.data2).rowwise().norm().sum()/in.data1.rows(); //Add all the differences
	H = (xy_map - in.data[1]).rowwise().squaredNorm().sum() / in.data[0].rows(); //Add all the differences

	guy.value = H; //Record the total distance between mappings


	//(DO NOT CHANGE) I propose to use the function below to make H small overall, and sharp close to H = 0
	H = -1 / log(H + log_param) + log_const + pow(log( 1/(H + 1)), 2);
	return H;
}

//extern void map(ArrayXd &, MatrixXd &, MatrixXd &, int);
//extern double fitnessTest(personality&, inData&);
#endif
