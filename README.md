 
# BHM
Bin histogram method for restoration of smooth functions from noisy integrals

## 1. Directories

    src     BHM headers and source files, example program main.cpp
    tests   BHM unit tests and other tests, header file sput.h required for unit tests

## 2. To run the code

- CMakeLists.txt in src builds main.cpp containing an example of the code  
  parameters can be adjusted in the main.cpp file and there is a choice of several test functions  
  output in build/src  
  
- CMakeLists.txt in tests builds unit tests and other test programs  
  (edit file to change the name of the test program)  
  output in build/tests  
  
## 3. List of files

### src

    basic.\*pp       definition of utility functions, definition of constants and global parameters in basic.hpp  
    slot.\*pp        classes related to histogram bins, including basis projections  
    spline.\*pp      classes related to splines, functions for chi^2 minimization (design matrix initialization, SVD, etc.)  
    histogram.\*pp   heart of the program: class defining the histogram, BHM fit method is implemented there  
    main.cpp         example program  
        
### tests   

    basic_unittest.cpp          unit tests for basic.cpp  
    slot_unittest.cpp           unit tests for slot.cpp  
    spline_unittest.cpp         unit tests for spline.cpp  
    histogram_unittest.cpp      unit tests for histogram.cpp  
        
You can ignore the following tests -- they were used to check the method and produce plots for the paper;  
they have not been commented and might not be well readable  
        
    bootstrap_test.cpp          tests bootstrap error estimate  
    timeevolution_test.cpp      tests fit evolution error estimate  
    error_test.cpp              independent simulation runs to analyse statistics and for robust error estimate  
    threshold_test.cpp          tests how fit acceptance threshold affects the fit result  
    divergence_test.cpp         sampling a divergent function on a semi-infinite domain using appropriate transforms  
    signproblem_test.cpp        sampling function with Monte Carlo Markov chain updates that have severe sign problem  
    signproblem_norm_test.cpp   same as above, but explicit normalization sector  

## 4. To Do

- interface for including into external code
- more safety mechanisms in case user inputs bad parameters
- current output is verbose, this should be optional
- more unit tests and examples?
- unit tests should be part of one test suite
- potential mistakes or not spotted memory leaks
- limiting cases not accounted for that might cause a crash
- improve efficiency where possible
- potentially too many checks that slow things down? checks too strict?
- potential refactoring to improve readability, e.g. some functions have many arguments and are long


