/*
* use Eigen3 Library to make state-space computation
* could be used with Simulink Coder and e.g. RTAI
* then, AVX, AVX2 and SSE3 could be used to accelerate state-space
* multiplication (useful for huge controller/systems)
* 
* in contrast to sfun_stateSpace.cpp this is the discrete version
* of the state space implementation
* 
* compile it for normal simulation as follows:
*   mex -O -I/usr/local/include/eigen3/ sfun_stateSpace_discrete.cpp
* 
* [of course eigen3 has to be installed to the mentioned directory]
* 
* based on stspace.c of The MathWorks, Inc.
*
* copyright by:
* Steffen Mauch, (c) 01/2015
* email: steffen.mauch (at) gmail.com
*
* You can redistribute it and/or modify it under the terms of the GNU 
* General Public License as published by the 
* Free Software Foundation, version 2.
* 
* This program is distributed in the hope that it will be useful, but WITHOUT
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
* FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
* details.
* 
* You should have received a copy of the GNU General Public License along with
* this program; if not, write to the Free Software Foundation, Inc., 51
* Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#define S_FUNCTION_NAME sfun_stateSpace_discrete /* Defines and Includes */
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"
#include <Eigen/SVD>
#include <Eigen/Core>
using namespace Eigen;

#include <iostream>
#include <string.h>

/* 
 * define EIGEN3_VARIANT if Eigen3 should be used!
 */
#define EIGEN3_VARIANT

#define NUMBER_OF_PARAMS	(6)

#define MATA_PARAM(S)	ssGetSFcnParam(S,0)
#define MATB_PARAM(S)	ssGetSFcnParam(S,1)
#define MATC_PARAM(S)	ssGetSFcnParam(S,2)
#define MATD_PARAM(S)	ssGetSFcnParam(S,3)
#define X0_PARAM(S)		ssGetSFcnParam(S,4)
#define SAMPLE_TIME_PARAM(S)	ssGetSFcnParam(S,5)

#define SAMPLE_TIME		((real_T) mxGetPr(SAMPLE_TIME_PARAM(S))[0])
#define NSTATES   mxGetM(MATA_PARAM(S))
#define NINPUTS   mxGetN(MATB_PARAM(S))
#define NOUTPUTS  mxGetM(MATC_PARAM(S))

#define IS_PARAM_DOUBLE(pVal) (mxIsNumeric(pVal) && !mxIsLogical(pVal) &&\
!mxIsEmpty(pVal) && !mxIsSparse(pVal) && !mxIsComplex(pVal) && mxIsDouble(pVal))

#define OK_EMPTY_DOUBLE_PARAM(pVal) (mxIsNumeric(pVal) && !mxIsLogical(pVal) &&\
!mxIsSparse(pVal) && !mxIsComplex(pVal) && mxIsDouble(pVal)) 


#define MDL_CHECK_PARAMETERS
#if defined(MDL_CHECK_PARAMETERS) && defined(MATLAB_MEX_FILE)
/* Function: mdlCheckParameters =============================================
 * Abstract:
 *    Validate our parameters to verify they are valid.
 */
static void mdlCheckParameters(SimStruct *S)
{
	/* Check A-Matrix */
	if ( mxGetN(MATA_PARAM(S)) != NSTATES || !IS_PARAM_DOUBLE(MATA_PARAM(S)) ) {
		ssSetErrorStatus(S,"1st parameter to S-function "
							"\"A-Matrix\" must be square and double");
		return;
	}

	/* Check B-Matrix */
	if (mxGetM(MATB_PARAM(S)) != NSTATES || !IS_PARAM_DOUBLE(MATB_PARAM(S)) ) {
		ssSetErrorStatus(S,"2nd parameter to S-function "
							"\"B-Matrix\" is not dimensioned "
							"correctly");
		return;
	}

	/* Check C-Matrix */
	if (mxGetN(MATC_PARAM(S)) != NSTATES || !IS_PARAM_DOUBLE(MATC_PARAM(S)) ) {
		ssSetErrorStatus(S,"3rd parameter to S-function "
						"\"C-Matrix\" is not dimensioned "
						"correctly");
		return;
	}
	
	/* Check D-Matrix */
	if (mxGetM(MATD_PARAM(S)) != NOUTPUTS || 
		mxGetN(MATD_PARAM(S)) != NINPUTS || !IS_PARAM_DOUBLE(MATD_PARAM(S)) ) {
		ssSetErrorStatus(S,"4th parameter to S-function "
							"\"D-Matrix\" is not dimensioned "
							"correctly");
		return;
	}
	
	/* Check X0, initial state */
	if ( ((mxGetM(X0_PARAM(S)) != 0) && 
		(mxGetM(X0_PARAM(S)) != NSTATES)) || !OK_EMPTY_DOUBLE_PARAM(X0_PARAM(S)) ) {
		ssSetErrorStatus(S,"5th parameter to S-function "
							"\"X0-Matrix\" is not dimensioned "
							"correctly");
		return;
	}
	
	/* Check sample time */
	if ( ((mxGetM(SAMPLE_TIME_PARAM(S)) != 0) && 
		(mxGetM(SAMPLE_TIME_PARAM(S)) != 1)) || !OK_EMPTY_DOUBLE_PARAM(SAMPLE_TIME_PARAM(S)) ) {
		ssSetErrorStatus(S,"6th parameter to S-function "
							"sample time must be scalar");
		return;
	}
}
#endif /* MDL_CHECK_PARAMETERS */

static void mdlInitializeSizes(SimStruct *S)
{
	int_T i;
	
	ssSetNumSFcnParams(S, NUMBER_OF_PARAMS);
	if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
		return; /* Parameter mismatch reported by the Simulink engine*/
	}

	for (i = 0; i < NUMBER_OF_PARAMS; i++) {
		ssSetSFcnParamNotTunable(S, i);
	}
	
	ssSetNumContStates(S, 0);
	ssSetNumDiscStates(S, NSTATES);

	if (!ssSetNumInputPorts(S, 1)) return;
	
	ssSetInputPortWidth(S, 0, NINPUTS);
	ssSetInputPortDirectFeedThrough(S, 0, 1);
    ssSetInputPortRequiredContiguous(S, 0, 1);
	ssSetInputPortDataType( S, 0, SS_DOUBLE );

	if (!ssSetNumOutputPorts(S, 1)) return;
	ssSetOutputPortWidth(S, 0, NOUTPUTS);

	ssSetNumSampleTimes(S, 1);
	ssSetNumRWork(S, 0);
	ssSetNumIWork(S, 0);
	ssSetNumPWork(S, 0);
	ssSetNumModes(S, 0);
	ssSetNumNonsampledZCs(S, 0);

	/* Take care when specifying exception free code - see sfuntmpl.doc */
	ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);
}

static void mdlInitializeSampleTimes(SimStruct *S)
{
	ssSetSampleTime(S, 0, SAMPLE_TIME);
	ssSetOffsetTime(S, 0, 0.0);
}

#define MDL_INITIALIZE_CONDITIONS
/* Function: mdlInitializeConditions ========================================
 * Abstract:
 *    If the initial condition parameter (X0) is not empty,
 *    use it as initial conditions, otherwise zero it
 */
static void mdlInitializeConditions(SimStruct *S)
{
	real_T *x0 = ssGetDiscStates(S);
	int_T  i, nStates;
 
	nStates = ssGetNumDiscStates(S);
	if (mxGetM(X0_PARAM(S)) != 0) {
		const real_T *pr = mxGetPr(X0_PARAM(S));

		for (i = 0; i < nStates; i++) {
			*x0++ = *pr++;
		}
	} else {
		for (i = 0; i < nStates; i++) {
			*x0++ = 0.0;
		}
	}
}


/* Function: mdlOutputs =======================================================
 * Abstract:
 *      y = Cx + Du
 */
static void mdlOutputs(SimStruct *S, int_T tid)
{
	real_T            *y       = ssGetOutputPortRealSignal(S,0);
	real_T            *x       = ssGetDiscStates(S);   
    double            *u       = (double*)ssGetInputPortRealSignal(S,0);
	const real_T      *cpr     = mxGetPr(MATC_PARAM(S));
	const real_T      *dpr     = mxGetPr(MATD_PARAM(S));

#ifndef EIGEN3_VARIANT
	int_T             i, j;
	real_T            accum;
#endif

    
	UNUSED_ARG(tid); /* not used in single tasking mode */

#ifdef EIGEN3_VARIANT
	double *xx = (double*)&x[0];
	double *C = (double*)&cpr[0];
	double *D = (double*)&dpr[0];
    
	MatrixXd resMat;
    
	MatrixXd vec_x = Map<MatrixXd>(xx, NSTATES, 1);
	MatrixXd vec_u = Map<MatrixXd>(u, NINPUTS, 1);
    
	MatrixXd mat_C = Map<MatrixXd>(C, NOUTPUTS, NSTATES);
	MatrixXd mat_D = Map<MatrixXd>(D, NOUTPUTS, NINPUTS);

	resMat = mat_C*vec_x + mat_D*vec_u;
	//std::cout << resMat << std::endl;

	memcpy( y,resMat.data(),resMat.size()*sizeof(double));
 
#else
	/* Matrix Multiply: y = Cx + Du */
	for (i = 0; i < (int_T)NOUTPUTS; i++) {
		accum = 0.0;
 
		/* Cx */
		for (j = 0; j < (int_T)NSTATES; j++) {
			accum += cpr[i + (int_T)NOUTPUTS*j] * x[j];
		}
 
		/* Du */
		for (j = 0; j < (int_T)NINPUTS; j++) {
			accum += dpr[i + (int_T)NOUTPUTS*j] * u[j];
		}
 
		y[i] = accum;
	}
#endif

} 
 
#define MDL_UPDATE 
/* Function: mdlDerivatives =================================================
 * Abstract:
 *      x(k+1) = Ax(k) + Bu(k)
 */
static void mdlUpdate(SimStruct *S, int_T tid)
{
	real_T            *x      = ssGetDiscStates(S);
	double            *u      = (double*)ssGetInputPortRealSignal(S,0);
	const real_T      *apr    = mxGetPr(MATA_PARAM(S));
	const real_T      *bpr    = mxGetPr(MATB_PARAM(S));

  #ifndef EIGEN3_VARIANT
	int_T i, j;
	real_T accum;
 #endif
  
 #ifdef EIGEN3_VARIANT
	double *xx = (double*)&x[0];
	double *A = (double*)&apr[0];
	double *B = (double*)&bpr[0];
    
	MatrixXd resMat;
    
	MatrixXd vec_x = Map<MatrixXd>(xx, NSTATES, 1);
	MatrixXd vec_u = Map<MatrixXd>(u, NINPUTS, 1);
    
	MatrixXd mat_A = Map<MatrixXd>(A, NSTATES, NSTATES);
	MatrixXd mat_B = Map<MatrixXd>(B, NSTATES, NINPUTS);

	resMat = mat_A*vec_x + mat_B*vec_u;
	//std::cout << resMat << std::endl;
   	   	
	memcpy( x,resMat.data(),resMat.size()*sizeof(double));
 
 #else
	/* Matrix Multiply: x(k+1) = Ax(k) + Bu(k) */
	real_T tempX[(int_T)NSTATES];
 
	for (i = 0; i < (int_T)NSTATES; i++) {
		accum = 0.0;
 
		/* Ax */
		for (j = 0; j < (int_T)NSTATES; j++) {
			accum += apr[i + (int_T)NSTATES*j] * x[j];
		}
 
		/* Bu */
		for (j = 0; j < (int_T)NINPUTS; j++) {
			accum += bpr[i + (int_T)NSTATES*j] * u[j];
		}
 
		tempX[i] = accum;
	}
	for (i = 0; i < (int_T)NSTATES; i++) { 
		x[i] = tempX[i];
	}
	
#endif
}
 


static void mdlTerminate(SimStruct *S){
	UNUSED_ARG(S); /* unused input argument */
}

#ifdef MATLAB_MEX_FILE /* Is this file being compiled as a MEX-file? */
#include "simulink.c" /* MEX-file interface mechanism */
#else
#include "cg_sfun.h" /* Code generation registration function */
#endif
