Simulink-Eigen3
===============

Inside this repository, there are different Simulink Models based on 
Eigen3 to accelerate Simulink models by using vector extensions form the
processor;

e.g. for matrix multiplication.
These s-functions can also be used with Simulink Coder; e.g. to accelearte 
huge matrix multiplication by using AVX2, FMA, AVX, SSE3, ... on a x86 platform
running e.g. with RTAI.

The following different s-functions are available:


s-function for matrix multiplication with Eigen3
(just compile the s-function with mex and run Simulink model, see header)


s-function for state-space continuous computation with Eigen3
(just compile the s-function with mex and run Simulink model, see header)


s-function for state-space discrete computation with Eigen3
(just compile the s-function with mex and run Simulink model, see header)


s-function for calculating based on slopes of an SHWFS the Zernike 
coefficients till order 28 with Eigen3
(compile the s-function with mex and run the workspace_tbSlopeToZernike.m 
script)
