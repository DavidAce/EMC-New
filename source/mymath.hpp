#ifndef MYMATH_H   
#define MYMATH_H  
#include <Eigen/Dense>
#include <Eigen/Core>
using namespace Eigen;
extern void linspace(ArrayXd &,double, double, int);
extern void logspace(ArrayXd &,double, double, int);
extern void rndChoice(int [], int , int);

extern int get_sign(double);
extern int tri_randint(int, int );
extern int tri_inv_randint(int, int );
extern int mod(int , int );
extern int isvalueinarray(int , int*, int);
extern int heaviside(double );
extern int find_index(int [], int, int);
 
extern double mean(double[], int );
extern double var(double [], int);


#endif