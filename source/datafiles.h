#ifndef DATAFILES_H
#define DATAFILES_H
#include <iostream>
#include <memory>
#include <string>
#include <fstream>
#include "../Eigen/Dense"
#include "../Eigen/Core"

using namespace std;
using namespace Eigen;
class population; //Forward declaration
class species;

class outData {
private:
	ofstream guysData[M];
	ofstream eliteGuysData[M];
	string filename_guysData[M];
	string filename_eliteGuysData[M];
	string folder;
public:
	outData();
	void print_to_file(species &);
	void print_to_file(population &, int &);

};

class inData {
private:
	void ReadDataSize(ifstream &, size_t &, size_t &);
	void importData(ifstream &, string &, MatrixXd& );


public:
	inData(const int ,const char **);
	const int num_files;
	string folder;
	unique_ptr<ifstream[]>	file;
	unique_ptr<string[]>	filename;
	unique_ptr<MatrixXd[]>	data;

};

#endif
