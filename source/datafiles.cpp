#include <iostream>
#include <memory>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include "../Eigen/Dense"
#include "../Eigen/Core"
#include "constants.h"
#include "mymath.h"
#include "DNA.h"
#include "personality.h"
#include "population.h"
#include "species.h"
#include "datafiles.h"
#ifdef __linux__
#define os 0
#elif _WIN32
#define os 1
#else
#define os 2
#endif

namespace patch
{
    template < typename T > std::string to_string( const T& n )
    {
        std::ostringstream stm ;
        stm << n ;
        return stm.str() ;
    }
}

using namespace Eigen;
using namespace std;
using namespace constants;


void inData::ReadDataSize(ifstream &fp, size_t &cols, size_t &rows) {
	string line;
	getline(fp, line);
	istringstream is(line);
	vector<double> v;
	double n;
	int i = 1;
	while (is >> n) {
		v.push_back(n);
	}
	
	cols = v.size();
	while (getline(fp, line)) {
		i++;
	}
	rows = i;
	fp.clear();
	fp.seekg(0);
}
void inData::importData(ifstream &fp, string &fname, MatrixXd &mat) {
	fp.open(fname.c_str());
	streampos current = fp.tellg();
	fp.seekg(0, fp.end);
	
	bool empty = !fp.tellg(); // true if empty file
	fp.seekg(current, fp.beg); //restore stream position
	if (fp.is_open() == false || empty == true) {
		cout << "Error opening file: " << fname <<" | Does it exist? Is Open: "<< fp.is_open()<< endl;
		exit(1);
	}
	else {
		//Read number of columns
		size_t cols, rows;
		ReadDataSize(fp, cols, rows);
		unsigned int i, j;
		mat.resize(rows, cols);
		
		for (i = 0; i < rows; i++) {
			for (j = 0; j < cols; j++) {
				fp >> mat(i, j);
			}
		}
	}
	fp.close();
}
inData::inData(const int argc,const char **argv) :	num_files	(argc-1),
													file		(new ifstream[num_files]),
													filename	(new string[num_files]),
													data		(new MatrixXd[num_files]){


	string mdir;
	switch (os) {
	case 0:
		mdir = "mkdir -p indata";
		folder = "indata";
		break;
	case 1:
		mdir = "mkdir indata";
		folder = "indata";
		break;
	case 2:
		mdir = "mkdir -p indata";
		folder = "indata";
		break;
	}
	if(std::system(mdir.c_str())){}
	int i = 0;
	cout << "Input files:	" << num_files << endl;
	for (i = 0; i < num_files; i++) {
		filename[i] += folder + "/";
		filename[i].append(argv[i+1]);
		cout << "File #" << i+1 << ":	" << filename[i] << " ";
		importData(file[i], filename[i], data[i]);
		cout << "		| Size :	" << data[i].rows() << " x " <<  data[i].cols() << endl;

	}
	i--;
	bounds.all_bounds = data[i];
	bounds.lower_bound =bounds.all_bounds.col(0);
	bounds.upper_bound = bounds.all_bounds.col(1);
}


void outData::print_to_file(population &pop, int &m) {
	int i, j;
	const size_t MAXWIDTH = 15;
	for (j = 0; j < nGenes; j++) {
		guysData[m]		<< right << setw(MAXWIDTH-2) << std::setfill(' ') << "Parameter[" << j << "]";
		eliteGuysData[m]<< right << setw(MAXWIDTH-2) << std::setfill(' ') << "Parameter[" << j << "]";
	}
	guysData[m]	<< right << setw(MAXWIDTH) << std::setfill(' ') << "Value"
				<< right << setw(MAXWIDTH) << std::setfill(' ') << "Fitness"
				<< right << setw(MAXWIDTH) << std::setfill(' ') << "Temperature"
				<< right << setw(MAXWIDTH) << std::setfill(' ') << "Generation" << endl;
	eliteGuysData[m]	<< right << setw(MAXWIDTH) << std::setfill(' ') << "Value"
						<< right << setw(MAXWIDTH) << std::setfill(' ') << "Fitness"
						<< right << setw(MAXWIDTH) << std::setfill(' ') << "Temperature"
						<< right << setw(MAXWIDTH) << std::setfill(' ') << "Generation" << endl;
	for (i = 0; i < N; i++) {
		for (j = 0; j < nGenes; j++) {
			guysData[m] << right << setw(MAXWIDTH) << std::setfill(' ') << pop.guys[i].genome.parameters(j);
		}
		guysData[m] << right << setw(MAXWIDTH) << std::setfill(' ') << pop.guys[i].value
					<< right << setw(MAXWIDTH) << std::setfill(' ')	<< pop.guys[i].H
					<< right << setw(MAXWIDTH) << std::setfill(' ')	<< pop.guys[i].t
					<< right << setw(MAXWIDTH) << std::setfill(' ')	<< pop.guys[i].generation << endl;
	}
	for (i = 0; i < N_best; i++) {
		for (j = 0; j < nGenes; j++) {
			eliteGuysData[m] << right << setw(MAXWIDTH) << std::setfill(' ') << pop.bestguys[i].genome.parameters(j);
		}
		eliteGuysData[m]<< right << setw(MAXWIDTH) << std::setfill(' ') << pop.bestguys[i].value
						<< right << setw(MAXWIDTH) << std::setfill(' ') << pop.bestguys[i].H
						<< right << setw(MAXWIDTH) << std::setfill(' ') << pop.bestguys[i].t
						<< right << setw(MAXWIDTH) << std::setfill(' ') << pop.bestguys[i].generation << endl;
	}
	guysData[m].close();
	eliteGuysData[m].close();
}
void outData::print_to_file(species &sp) {
	for (int m = 0; m < M; m++) {
		print_to_file(sp.pop[m], m);
	}
}

outData::outData() {
	//Create folder for out data storage
	string mdir;
	switch (os) {
	case 0:
		mdir = "mkdir -p outdata";
		folder = "outdata";
		break;
	case 1:
		mdir = "mkdir ..\\outdata";
		folder = "..\\outdata";
		break;
	case 2:
		mdir = "mkdir -p outdata";
		folder = "outdata";
		break;
	}
	if(std::system(mdir.c_str())){}
	for (int m = 0; m < M; m++) {
		filename_guysData[m] = folder + "/parameters" + patch::to_string(m) +  ".dat";
		filename_eliteGuysData[m] = folder + "/elitedata" + patch::to_string(m) + ".dat";
		guysData[m].open(filename_guysData[m].c_str(), ofstream::out | ofstream::trunc);
		guysData[m] << fixed << showpoint << setprecision(10);
		eliteGuysData[m].open(filename_eliteGuysData[m].c_str(), ofstream::out | ofstream::trunc);
		eliteGuysData[m] << fixed << showpoint << setprecision(10);
	}
}
