# Evolutionary Monte Carlo (EMC)
This program finds the global minimum in a general parameter
space by combining Genetic algorithms and Monte Carlo 
algorithms. 

Read more here:

[Evolutionary Monte Carlo for protein folding simulations](http://users.phhp.ufl.edu/faliang/papers/2001/JCP2D.pdf)

[Real-parameter evolutionary Monte Carlo with applications to Bayesian mixture models](http://users.phhp.ufl.edu/faliang/papers/2001/RealEMC.pdf)
## Compilation 
Use `make` to generate the binary `EMC`. 
Use `make clean ` to clean old object files.


## Excecution  
Use `EMC file0 file1 ... fileN` to execute the 
binary EMC with your own input files listed.
Put your input files in a folder "./indata/".
The file run.sh is a bash script for execution.
 
                
All files are optional, except the one containing the
boundaries of your parameter space, which is
mandatory. See format below. This file must be the last
on the argument list.
Optional files may be used for computation of the fitness.
                
##### Example
        EMC my_file1.dat my_file2.dat boundaries.dat

                 
                
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
in a 3D cube parameter space. Then boundaries.dat may 
contain (tab delimited):
                        
                -1  3.5
                10  20
                2   5


## Constants
See the file `constants.h`. In particular, the 
variable "int generations" controls the maximum duration
of the simulation (default 500000).