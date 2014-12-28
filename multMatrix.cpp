/*
* use Eigen3 Library to multiplicate arbitrary matrices
* could be used with Simulink Coder and e.g. RTAI
* then, AVX, AVX2 and SSE3 could be used to accelerate matrix
* multiplication (useful for huge matrices)
* 
* compile it for normal simulation as follows:
*   mex -O -I/usr/local/include/eigen3/ multMatrix.cpp
* 
* [of course eigen3 has to be installed to the mentioned directory]
*
* copyright by:
* Steffen Mauch, (c) 12/2014
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
 
#define S_FUNCTION_NAME multMatrix /* Defines and Includes */
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"
#include <Eigen/SVD>
#include <Eigen/Core>
using namespace Eigen;

#include <iostream>
#include <string.h>

#define MDL_CHECK_PARAMETERS
#if defined(MDL_CHECK_PARAMETERS) && defined(MATLAB_MEX_FILE)
/* Function: mdlCheckParameters =============================================
 * Abstract:
 *    Validate our parameters to verify they are okay.
 */
static void mdlCheckParameters(SimStruct *S)
{
}
#endif /* MDL_CHECK_PARAMETERS */


static void mdlInitializeSizes(SimStruct *S)
{
	ssSetNumSFcnParams(S, 0);
	if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
		return; /* Parameter mismatch reported by the Simulink engine*/
	}
	
	if( !ssSetNumInputPorts(S, 2) ) return;
	if( !ssSetNumOutputPorts(S,1) ) return;

	if( !ssSetInputPortDimensionInfo( S, 0, DYNAMIC_DIMENSION) ) return;
	if( !ssSetInputPortDimensionInfo( S, 1, DYNAMIC_DIMENSION) ) return;

	if(!ssSetOutputPortDimensionInfo(S, 0, DYNAMIC_DIMENSION)) return;

	ssSetInputPortDataType( S, 0, SS_DOUBLE );
	ssSetInputPortDataType( S, 1, SS_DOUBLE );

	ssSetInputPortDirectFeedThrough(S, 0, 1);
	ssSetInputPortDirectFeedThrough(S, 1, 1);

	ssSetNumSampleTimes(S, 1);

    ssSetNumContStates(S, 0);
    ssSetNumDiscStates(S, 0);

    ssSetNumRWork(S, 0);
    ssSetNumIWork(S, 0);
    ssSetNumPWork(S, 0);
    ssSetNumModes(S, 0);

	/* Take care when specifying exception free code - see sfuntmpl.doc */
	ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);
}
    

#if defined(MATLAB_MEX_FILE)
#define MDL_SET_INPUT_PORT_DIMENSION_INFO
static void mdlSetInputPortDimensionInfo(SimStruct *S,
                                  int_T port,
                                  const DimsInfo_T *dimsInfo)
{       
    if(!ssSetInputPortDimensionInfo(S, port, dimsInfo)) return; 
}

# define MDL_SET_OUTPUT_PORT_DIMENSION_INFO
static void mdlSetOutputPortDimensionInfo(SimStruct        *S, 
                                          int_T            port,
                                          const DimsInfo_T *dimsInfo)
{
	/* This should never occur! */
    ssSetErrorStatus(S, "Error setting output port width.");
}

# define MDL_SET_DEFAULT_PORT_DIMENSION_INFO
static void mdlSetDefaultPortDimensionInfo(SimStruct *S)
{	  
	if (ssGetInputPortWidth(S, 0) == DYNAMICALLY_SIZED) {
		ssSetInputPortMatrixDimensions(S, 0, 1, 1 );
	}
    
	if (ssGetInputPortWidth(S, 1) == DYNAMICALLY_SIZED) {
		ssSetInputPortMatrixDimensions(S, 1, 1, 1 );
	}
	
	int_T port0A = ssGetInputPortDimensionSize(S,  0, 0);
    int_T port0B = ssGetInputPortDimensionSize(S,  0, 1);
    
    int_T port1A = ssGetInputPortDimensionSize(S,  1, 0);
    int_T port1B = ssGetInputPortDimensionSize(S,  1, 1);
    
	if( port0B == port1A )
		ssSetOutputPortMatrixDimensions(S, 0, port0A, port1B);
	else
		ssSetErrorStatus(S,"dimension for matrix multiplication do not match!");

}
#endif


static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, INHERITED_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);
}

static void mdlOutputs(SimStruct *S, int_T tid)
{
    InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S,0);
    InputRealPtrsType uPtrs1 = ssGetInputPortRealSignalPtrs(S,1);
    real_T *y = ssGetOutputPortRealSignal(S,0);
    
    double *u1 = (double*)uPtrs[0];
    double *u2 = (double*)uPtrs1[0];
    
    MatrixXd resMat;
    int_T iRows1 = ssGetInputPortDimensions(S, 0)[0];
    int_T iCols1 = ssGetInputPortDimensions(S, 0)[1];
    int_T iRows2 = ssGetInputPortDimensions(S, 1)[0];
    int_T iCols2 = ssGetInputPortDimensions(S, 1)[1];
    
    MatrixXd mat1 = Map<MatrixXd>(u1, iRows1, iCols1);
    MatrixXd mat2 = Map<MatrixXd>(u2, iRows2, iCols2);  

    resMat = mat1*mat2;

    memcpy( y,resMat.data(),resMat.size()*sizeof(double));
}

static void mdlTerminate(SimStruct *S){}

#ifdef MATLAB_MEX_FILE /* Is this file being compiled as a MEX-file? */
#include "simulink.c" /* MEX-file interface mechanism */
#else
#include "cg_sfun.h" /* Code generation registration function */
#endif
