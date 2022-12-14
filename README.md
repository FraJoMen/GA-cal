# GA-cal
[![DOI](https://zenodo.org/badge/537369162.svg)](https://zenodo.org/badge/latestdoi/537369162)

GA-cal is a Fortran software for automatically calibrating constitutive laws using Genetic Algorithms (GA) optimization. The proposed approach sets the calibration problem as a regression, and the GA optimization is used to adjust the model parameters so that a numerical model matches experimental data. Currently, the code allows the calibration of the Sand Hypoplastic law (SH), proposed by von Wolffersdorff, with the oedometer (OE) and triaxial drained (TD) test data. The implemented subroutines can be easily extended to solve other regression or optimization problems, including different tests and constitutive models.

A copy of the **User manual** is downloadable from the website [https://arxiv.org/abs/2211.13652](https://arxiv.org/abs/2211.13652).

The basic steps in using the code are:
  1. define the experimental data;
  2. specify the initial condition of the tests;
  3. set the optimization parameters;
  4. run the code;
  5. analyze the results. 

In the GitHub repository, you will find the following:
  1. the folder **Examples** contains material helpful in familiarizing and getting started with GA-cal.
  2. the **GA-cal Source code** folder contains the files necessary for compiling the program and creating the executable. The folder contains the project developed with the [Code::Blocks](https://www.codeblocks.org/) environment and the Fortran code. The compilation has been tested and developed using [gfortran](https://gcc.gnu.org/wiki/GFortran). 
  3. the folder **doc** contains a copy of the manual, the .dll libraries of software dependency and some instructions for using them. 

## Reference articles
If you use GA-cal, please cite this reference in your work (books, articles, reports, etc.)
 
  - F. J. Mendez, A. Pasculli, M. A. Mendez, N. Sciarra, *Calibration of a hypoplastic model using genetic algorithms*, Acta Geotechnica 16 (2021) 2031–2047. [https://doi.org/10.1007/s11440-020-01135-z](https://doi.org/10.1007/s11440-020-01135-z) 
  - F. J. Mendez, M. A. Mendez, A. Pasculli, *The GA-cal software for the automatic calibration of soil constitutive laws: a tutorial and a user manual*, arXiv (2022). [https://doi.org/10.48550/arXiv.2211.13652](https://doi.org/10.48550/arXiv.2211.13652) 
  
