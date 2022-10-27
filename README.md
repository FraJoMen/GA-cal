# GA-cal
GA-cal is a Fortran software for automatically calibrating constitutive laws using Genetic Algorithms (GA) optimization. The proposed approach sets the calibration problem as a regression, and the GA optimization is used to adjust the model parameters so that a numerical model matches experimental data. Currently, the code allows the calibration of the Sand Hypoplastic law (SH), proposed by von Wolffersdorff, with the oedometer (OE) and triaxial drained (TD) test data. The implemented subroutines can be easily extended to solve other regression or optimization problems, including different tests and constitutive models.

The basic steps in using the code are:
  1. define the experimental data;
  2. specify the initial condition of the tests;
  3. set the optimization parameters;
  4. run the code;
  5. analyze the results. 

In the GitHub repository, you will find the following:
  1. a copy of the **User manual** can also be downloaded from the website .... If you use GA-cal, please cite this reference in your work (books, articles, reports, etc.).
  2. the **GA-cal Source code** folder contains the files necessary for compiling the program and creating the executable file. In particular, the folder contains the project developed with the Code::Blocks environment and the Fortran code. The compilation has been tested and developed using Gfortran. 
  3. the folder **Examples** contains material helpful in familiarizing and getting started with GA-cal.

## Reference articles.

  - F. J. Mendez, A. Pasculli, M. A. Mendez, N. Sciarra, *Calibration of a hypoplastic model using genetic algorithms*, Acta Geotechnica 16 (2021) 2031–2047. [https://doi.org/10.1007/s11440-020-01135-z](https://doi.org/10.1007/s11440-020-01135-z) 
