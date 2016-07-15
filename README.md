# Evolutionary Monte Carlo (EMC)
This program finds the global minimum in a general parameter
space by combining Genetic algorithms and Monte Carlo 
algorithms. 

References:

[On learning strategies for evolutionary Monte Carlo](http://www.people.fas.harvard.edu/~junliu/TechRept/07folder/Goswami%26Liu07.pdf)

[Evolutionary Monte Carlo for protein folding simulations](http://users.phhp.ufl.edu/faliang/papers/2001/JCP2D.pdf)

[Real-parameter evolutionary Monte Carlo with applications to Bayesian mixture models](http://users.phhp.ufl.edu/faliang/papers/2001/RealEMC.pdf)

## Quick Start
From a Linux terminal type

		make
		./EMC xy_old.dat xy_new.dat boundaries.dat
                		

## Usage
The files `source/minimization.hpp` and `source/minimization.cpp`
have the function `fitnessTest()` which is used to map an n-dimensional
point in your parameter space, called a *chromosome*,
to scalar `H`, called the *fitness* (or energy). You will have to re-write this 
function to suit your minimizing needs. Since this function is called **many**
times it is important to write fast code. Do **not** use OpenMP as it is
already used on the outer loop that calls `fitnessTest()`.

Ideally you should aim for H to be small, no larger than 50, say, since differences
in H are present in many Boltzmann weight exponentials, as in `exp((H_new - H_old)/T)`.

If H can become very large, I suggest you leave the function that is already present:

        	H = -1 / log(H + log_param) + log_const + pow(log( 1/(H + 1)), 2);

which basically gives your parameter hypersurface a good shape for the minimization process.

### The default program
If you compile and run the program *as is*, it will try to find
the 12 coefficients of a bivariate polynomial that maps the coordinates
in `indata/xy_old.dat` into `indata/xy_new.dat`. This type of problem is 
found in cartography when one needs a recipe that stretches one
geographical map into another, that may come from a different projection.
All coefficients should approximate the value `2` if the program is working correctly.

## Compilation 
In a linux terminal, enter `make` to generate the binary `EMC`. 
Enter `make clean ` to clean old object files.
The program has been tested with g++ version 5.3.0, 
but should work with anything later than g++ 4.9.


## Excecution  
Enter `./EMC file0 file1 ... fileN` to execute the 
binary EMC with your own input files listed.
Put your input files in a folder "./indata/".
The file run.sh is a bash script for execution.
 
                
All files are optional, except the one containing the
boundaries of your parameter space, which is
mandatory. See format below. This file must be the last
on the argument list.
Optional files may be used for computation of the fitness.
                
##### Example
        ./EMC my_file1.dat my_file2.dat boundaries.dat

                 
                
## File format
Files can contain tab-delimited matrices of size n x m
(n = rows, m = columns).
A matrix from `file0` is accessible from within the 
program as  `in.data[0]`  , which is a pointer to
an Eigen-type matrix, MatrixXd.
    
##### Example
To get element in row 2 and column 5 from 
`file1`, you may write:
                
        double my_number = in.data[1](2,5);
                       
See the Quick Reference to Eigen syntax here:
https://eigen.tuxfamily.org/dox-devel/AsciiQuickReference.txt


The last argument `fileN` must always contain the domain
boundary in two columns. The first column is the lower 
bound, and the second is the upper bound.
There must be as many rows as parameters for the fit!

##### Example 
If parameters `sigma`, `mu`, `rho` are to be minimized
in a 3D cube parameter space with volume `64`, centered at origo,
then boundaries.dat may contain (tab delimited):
                        
                -2  2
                -2  2
                -2  2


## Constants
See the file `constants.h`. In particular, set the following
variables to your liking: 

`int M` sets the number of parallel subpopulations. Preferrably this
number should be the same as the number cpu-threads available for OpenMP.
(Usually 4 to 8, depending on your cpu).



`int generations` controls the maximum duration
of the simulation (default `500000`). 

`double lowest_H` sets threshold for the lowest fitness. If
any fitness is lower than `lowest_H` the program will terminate.


