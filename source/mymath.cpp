#include <math.h>
#include <iostream>
#include "randomFunctions.hpp"
#include "mymath.hpp"
#include <Eigen/Dense>
#include <Eigen/Core>
using namespace Eigen;
using namespace std;
void linspace(ArrayXd & array, double xi, double xf, int n) {
	int i;
	double dx, range = xf - xi;
	dx = range / (n - 1);
	array(0) = xi;
	for (i = 1; i < n; i++) {
		array(i) = array(i - 1) + dx;
	}
}
void logspace(ArrayXd & array, double xi, double xf, int n) {
	int i;
	double s = log(xi);
	double f = log(xf);
	linspace(array, s, f, n);
	for (i = 0; i < n; i++) {
		array(i) = exp(array(i));
		//array[i] = pow(10, array[i]);
	}
}
int get_sign(double value)
{
	if (value == 0)
		return 0;
	else if (value > 0)
		return 1;
	else
		return -1;
}
void rndChoice(int randomArray[], int nn, int NN) {
	//Choose nn integers from the set [0,NN-1] without repetition. Used the Knuth algorithm found on the interwebz (stack overflow). Ordo(Npmove).
	int in, im;
	int rn, rm;
	im = 0;
	for (in = 0; in < NN && im < nn; in++) {
		rn = NN - in;
		rm = nn - im;
		if (uniform_double(&rng, 0, 1)*rn < rm) {
			randomArray[im] = in;
			im += 1;
		}
	}
}
int tri_randint(int a, int b) {
	//Returns number between a and b with triangular distribution
	return (int)(a + pow(uniform_double(&rng, 0, 1), 0.5)*(b - a));
}
int tri_inv_randint(int a, int b) {
	//Returns number between a and b with triangular distribution, sloping downwards
	return (int)(b - pow(1 - uniform_double(&rng, 0, 1), 0.5)*(b - a));
}
double mean(double array[], int len) {
	double sum = 0;
	int i;
	for (i = 0; i < len; i++) {
		sum += array[i];
	}
	return sum / len;
}
int mod(int a, int b) {
	if (b < 0)
		return mod(-a, -b);

	int ret = a % b;
	if (ret < 0)
		ret += b;
	return ret;

}
int isvalueinarray(int val, int *arr, int size) {
	int i;
	for (i = 0; i < size; i++) {
		if (arr[i] == val)
			return 1;
	}
	return 0;
}
int heaviside(double x) {
	if (x > 0) {
		return 1;
	}
	else {
		return 0;
	}
}
double var(double a[], int n) {
	if (n == 0) {
		return 0.0;
	}
	double sum = 0;
	for (int i = 0; i < n; i++) {
		sum += a[i];
	}
	double mean = sum / n;
	double sq_diff_sum = 0;
	for (int i = 0; i < n; ++i) {
		double diff = a[i] - mean;
		sq_diff_sum += diff * diff;
	}
	double variance = sq_diff_sum / n;
	return variance;
}
int find_index(int a[], int n, int value) {
	int i;
	for (i = 0; i < n; i++)
	{
		if (a[i] == value)
		{
			return i;  /* it was found */
		}
	}
	return(-1);  /* if it was not found */
}