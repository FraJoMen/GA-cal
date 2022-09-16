# GA-cal
GA-cal is a Fortran software for automatically calibrating constitutive laws using Genetic Algorithms
(GA) optimization. The proposed approach sets the calibration problem as a regression, and the GA
optimization is used to adjust the model parameters so that a numerical model matches experimental
data. We showcase the software on the calibration of the Sand Hypoplastic law (SH), proposed by
von Wolffersdorff, with the oedometer (OE) and triaxial drained (TD) tests data. The implemented
subroutines can be easily extended to solve other regression or optimization problems, including
different tests and constitutive models

The basic steps in using the code are:
  1. define the experimental data;
  2. specify the initial condition of the tests;
  3. set the optimization parameters;
  4. run the code;
  5. use the Py_Out.py to analyze the results 
