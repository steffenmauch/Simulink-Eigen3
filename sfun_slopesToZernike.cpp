/*
* 
* compile it for normal simulation as follows:
*   mex -O -I/usr/local/include/eigen3/ sfun_slopesToZernike.cpp
*
* ($Id: sfun_slopesToZernike.cpp 1569 2015-03-04 21:37:21Z smauch $)
* 
* copyright by:
* Steffen Mauch, (c) 03/2015
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

#define S_FUNCTION_NAME sfun_slopesToZernike /* Defines and Includes */
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"
#include <Eigen/SVD>
#include <Eigen/Core>
#include <time.h>
using namespace Eigen;

#include <iostream>
#include <string.h>


#define NUMBER_OF_PARAMS	(3)

#define SLOPES_PARAM        ssGetSFcnParam(S,0)
#define DECIMAL_PARAM   	ssGetSFcnParam(S,1)
#define ZERNIKE_PARAM   	ssGetSFcnParam(S,2)

#define SLOPES_LENGTH       ((uint_T) mxGetPr(SLOPES_PARAM)[0])
#define CNT_OUT             ((uint_T) mxGetPr(DECIMAL_PARAM)[0])
#define ZERNIKE_CNT         ((uint_T) mxGetPr(ZERNIKE_PARAM)[0])

//#define MEASURE_TIME
#ifdef  MEASURE_TIME
    #ifndef MATLAB_MEX_FILE
        extern "C" uint_T rt_get_cpu_time_ns();
    #endif
#endif

template<typename _Matrix_Type_>
bool pseudoInverse(const _Matrix_Type_ &a, _Matrix_Type_ &result, double
epsilon = std::numeric_limits<typename _Matrix_Type_::Scalar>::epsilon())
{
    if(a.rows()<a.cols())
        return false;
    Eigen::JacobiSVD< _Matrix_Type_ > svd = a.jacobiSvd(Eigen::ComputeFullU |
    Eigen::ComputeFullV);
    typename _Matrix_Type_::Scalar tolerance = epsilon * std::max(a.cols(),
    a.rows()) * svd.singularValues().array().abs().maxCoeff();
    result = svd.matrixV() * _Matrix_Type_( (svd.singularValues().array().abs() >
    tolerance).select(svd.singularValues().
    array().inverse(), 0) ).asDiagonal() * svd.matrixU().adjoint();
    
    return true;
}

MatrixXd *dx;
MatrixXd *dy;

MatrixXd *res;

void arrayofslopefunctionsZern( SimStruct *S, double pupildiam, uint_T maxTerm )
{
    dx = new MatrixXd(SLOPES_LENGTH,maxTerm);
	dy = new MatrixXd(SLOPES_LENGTH,maxTerm);
    
    res = new MatrixXd(ZERNIKE_CNT,SLOPES_LENGTH*2);
    
	MatrixXd dx_int = MatrixXd::Zero(SLOPES_LENGTH,maxTerm);
	MatrixXd dy_int = MatrixXd::Zero(SLOPES_LENGTH,maxTerm);
    
    MatrixXd x = MatrixXd::Zero(14,14);
    MatrixXd y = MatrixXd::Zero(14,14);
    
    uint_T k = 0;
    uint_T l = 0;
    
    for( k=0; k<14; k++ ){
        x.row(k) = ArrayXd::LinSpaced(14, -1.0/pupildiam*2, 1.0/pupildiam*2);
        y.col(k) = ArrayXd::LinSpaced(14, -1.0/pupildiam*2, 1.0/pupildiam*2);
    }
    
    MatrixXd vec_x = VectorXd(Map<VectorXd>(x.data(), x.cols()*x.rows()));
    MatrixXd vec_y = VectorXd(Map<VectorXd>(y.data(), y.cols()*y.rows()));
    MatrixXd vec_ones = VectorXd::Ones(y.cols()*y.rows());
    
    MatrixXd x_pow2 = vec_x.cwiseProduct(vec_x);
    MatrixXd y_pow2 = vec_y.cwiseProduct(vec_y);
    
    MatrixXd x_pow3 = x_pow2.cwiseProduct(vec_x);
    MatrixXd y_pow3 = y_pow2.cwiseProduct(vec_y);
    
    MatrixXd x_pow4 = x_pow3.cwiseProduct(vec_x);
    MatrixXd y_pow4 = y_pow3.cwiseProduct(vec_y);
    
    MatrixXd x_pow5 = x_pow4.cwiseProduct(vec_x);
    MatrixXd y_pow5 = y_pow4.cwiseProduct(vec_y);
    
    if( maxTerm >= 1 ){
        dx_int.col(0) = ArrayXd::Zero( SLOPES_LENGTH );
        dy_int.col(0) = ArrayXd::Zero( SLOPES_LENGTH );
    }
    if( maxTerm >= 2 ){
        dy_int.col(1) = ArrayXd::Ones( SLOPES_LENGTH );
    }
    if( maxTerm >= 3 ){
        dx_int.col(2) = ArrayXd::Ones( SLOPES_LENGTH );
    }
    if( maxTerm >= 4 ){
        dx_int.col(3) = 2*vec_y;
        dy_int.col(3) = 2*vec_x;
    }
    if( maxTerm >= 5 ){
        dx_int.col(4) = 4*vec_x;
        dy_int.col(4) = 4*vec_y;
    }
    if( maxTerm >= 6 ){
        dx_int.col(5) = 2*vec_x;
        dy_int.col(5) = -2*vec_y;
    }
    if( maxTerm >= 7 ){
        dx_int.col(6) = 6*vec_x.cwiseProduct(vec_y);
        dy_int.col(6) = 3*x_pow2 - 3*y_pow2;
    }
    if( maxTerm >= 8 ){
        dx_int.col(7) = 6*vec_x.cwiseProduct(vec_y);
        dy_int.col(7) = 3*x_pow2 + 9*y_pow2 - 2*vec_ones;
    }
    if( maxTerm >= 9 ){
        dx_int.col(8) = 9*x_pow2 + 3*y_pow2 - 2*vec_ones;
        dy_int.col(8) = 6*vec_x.cwiseProduct(vec_y);
    }
    if( maxTerm >= 10 ){
        dx_int.col(9) = 3*x_pow2 - 3*y_pow2;
        dy_int.col(9) = -6*vec_x.cwiseProduct(vec_y);
    }
    if( maxTerm >= 11 ){
        dx_int.col(10) = 12*x_pow2.cwiseProduct(vec_y) - 4*y_pow3;
        dy_int.col(10) = 4*x_pow3 - 12*vec_x.cwiseProduct(y_pow2);
    }
    if( maxTerm >= 12 ){
        dx_int.col(11) = 24*x_pow2.cwiseProduct(vec_y) + 8*y_pow3 - 6*vec_y;
        dy_int.col(11) = 8*x_pow3 + 24*vec_x.cwiseProduct(y_pow2) - 6*vec_x;
    }
    if( maxTerm >= 13 ){
        dx_int.col(12) = 24*x_pow3 + 24*vec_x.cwiseProduct(y_pow2) - 12*vec_x;
        dy_int.col(12) = 24*x_pow2.cwiseProduct(vec_y) + 24*y_pow3 - 12*vec_y;
    }
    if( maxTerm >= 14 ){
        dx_int.col(13) = 16*x_pow3 + 8*vec_x.cwiseProduct(y_pow2) - 6*vec_x -8*vec_x.cwiseProduct(y_pow2);
        dy_int.col(13) = 8*x_pow2.cwiseProduct(vec_y) - 8*x_pow2.cwiseProduct(vec_y) - 16*y_pow3 + 6*vec_y;
    }
    if( maxTerm >= 15 ){
        dx_int.col(14) = 4*x_pow3 - 12*vec_x.cwiseProduct(y_pow2);
        dy_int.col(14) = -12*x_pow2.cwiseProduct(vec_y) + 4*y_pow3;
    }
    if( maxTerm >= 16 ){
        dx_int.col(15) = 20*x_pow3.cwiseProduct(vec_y) - 20*vec_x.cwiseProduct(y_pow3);
        dy_int.col(15) = 5*x_pow4 - 30*x_pow2.cwiseProduct(y_pow2) + 5*y_pow4;
    }
    if( maxTerm >= 17 ){
        dx_int.col(16) = 60*x_pow3.cwiseProduct(vec_y) + 30*vec_x.cwiseProduct(y_pow3) - 24*vec_x.cwiseProduct(vec_y) - 10*vec_x.cwiseProduct(y_pow3);
        dy_int.col(16) = 15*x_pow4 + 45*x_pow2.cwiseProduct(y_pow2) - 12*x_pow2 - 15*x_pow2.cwiseProduct(y_pow2) - 25*y_pow4 + 12*y_pow2;
    }
    if( maxTerm >= 18 ){
        dx_int.col(17) = 40*x_pow3.cwiseProduct(vec_y) + 40*vec_x.cwiseProduct(y_pow3) - 24*vec_x.cwiseProduct(vec_y);
        dy_int.col(17) = 10*x_pow4 + 60*x_pow2.cwiseProduct(y_pow2) + 50*y_pow4 - 12*x_pow2 - 36*y_pow2 + 3*vec_ones;
    }
    if( maxTerm >= 19 ){
        dx_int.col(18) = 50*x_pow4 + 60*x_pow2.cwiseProduct(y_pow2) + 10*y_pow4 - 36*x_pow2 - 12*y_pow2 + 3*vec_ones;
        dy_int.col(18) = 40*x_pow3.cwiseProduct(vec_y) + 40*vec_x.cwiseProduct(y_pow3) - 24*vec_x.cwiseProduct(vec_y);
    }
    if( maxTerm >= 20 ){
        dx_int.col(19) = 25*x_pow4 + 15*x_pow2.cwiseProduct(y_pow2) - 12*x_pow2 - 45*x_pow2.cwiseProduct(y_pow2) - 15*y_pow4 + 12*y_pow2;
        dy_int.col(19) = 10*x_pow3.cwiseProduct(vec_y) - 30*x_pow3.cwiseProduct(vec_y) - 60*vec_x.cwiseProduct(y_pow3) + 24*vec_x.cwiseProduct(vec_y);
    }
    if( maxTerm >= 21 ){
        dx_int.col(20) = 5*x_pow4 - 30*x_pow2.cwiseProduct(y_pow2) + 5*y_pow4;
        dy_int.col(20) = -20*x_pow3.cwiseProduct(vec_y) + 20*vec_x.cwiseProduct(y_pow3);
    }
    if( maxTerm >= 22 ){
        dx_int.col(21) = 30*x_pow4.cwiseProduct(vec_y) - 60*x_pow2.cwiseProduct(y_pow3) + 6*y_pow5;
        dy_int.col(21) = 6*x_pow5 - 60*x_pow3.cwiseProduct(y_pow2) + 30*vec_x.cwiseProduct(y_pow4);
    }
    if( maxTerm >= 23 ){
        dx_int.col(22) = 120*x_pow4.cwiseProduct(vec_y) + 72*x_pow2.cwiseProduct(y_pow3) - 60*x_pow2.cwiseProduct(vec_y) - 72*x_pow2.cwiseProduct(y_pow3) - 24*y_pow5 + 20*y_pow3;
        dy_int.col(22) = 24*x_pow5 + 72*x_pow3.cwiseProduct(y_pow2) - 20*x_pow3 - 72*x_pow3.cwiseProduct(y_pow2) - 120*vec_x.cwiseProduct(y_pow4) + 60*vec_x.cwiseProduct(y_pow2);
    }
    if( maxTerm >= 24 ){
        dx_int.col(23) = 150*x_pow4.cwiseProduct(vec_y) + 180*x_pow2.cwiseProduct(y_pow3) + 30*y_pow5 - 120*x_pow2.cwiseProduct(vec_y) - 40*y_pow3 + 12*vec_y;
        dy_int.col(23) = 30*x_pow5 + 180*x_pow3.cwiseProduct(y_pow2) + 150*vec_x.cwiseProduct(y_pow4) - 40*x_pow3 - 120*vec_x.cwiseProduct(y_pow2) + 12*vec_x;
    }
    if( maxTerm >= 25 ){
        dx_int.col(24) = 120*x_pow5 + 240*x_pow3.cwiseProduct(y_pow2) + 120*vec_x.cwiseProduct(y_pow4) - 120*x_pow3 - 120*vec_x.cwiseProduct(y_pow2) + 24*vec_x;
        dy_int.col(24) = 120*x_pow4.cwiseProduct(vec_y) + 240*x_pow2.cwiseProduct(y_pow3) + 120*y_pow5 - 120*x_pow2.cwiseProduct(vec_y) - 120*y_pow3 + 24*vec_y;
    }
    if( maxTerm >= 26 ){
        dx_int.col(25) = 90*x_pow5 + 120*x_pow3.cwiseProduct(y_pow2) + 30*vec_x.cwiseProduct(y_pow4) - 80*x_pow3 - 40*vec_x.cwiseProduct(y_pow2) + 12*vec_x - 60*x_pow3.cwiseProduct(y_pow2) - 60*vec_x.cwiseProduct(y_pow4) + 40*vec_x.cwiseProduct(y_pow2);
        dy_int.col(25) = 60*x_pow4.cwiseProduct(vec_y) + 60*x_pow2.cwiseProduct(y_pow3) - 40*x_pow2.cwiseProduct(vec_y) - 30*x_pow4.cwiseProduct(vec_y) - 120*x_pow2.cwiseProduct(y_pow3) - 90*y_pow5 + 40*x_pow2.cwiseProduct(vec_y) + 80*y_pow3 - 12*vec_y;
    }
    if( maxTerm >= 27 ){
        dx_int.col(26) = 36*x_pow5 + 24*x_pow3.cwiseProduct(y_pow2) - 20*x_pow3 - 144*x_pow3.cwiseProduct(y_pow2) - 72*vec_x.cwiseProduct(y_pow4) + 60*vec_x.cwiseProduct(y_pow2) + 12*vec_x.cwiseProduct(y_pow4);
        dy_int.col(26) = 12*x_pow4.cwiseProduct(vec_y) - 72*x_pow4.cwiseProduct(vec_y) - 144*x_pow2.cwiseProduct(y_pow3) + 60*x_pow2.cwiseProduct(vec_y) + 24*x_pow2.cwiseProduct(y_pow3) + 36*y_pow5 - 20*y_pow3;
    }
    if( maxTerm >= 28 ){
        dx_int.col(27) = 6*x_pow5 - 60*x_pow3.cwiseProduct(y_pow2) + 30*vec_x.cwiseProduct(y_pow4);
        dy_int.col(27) = -30*x_pow4.cwiseProduct(vec_y) + 60*x_pow2.cwiseProduct(y_pow3) - 6*y_pow5;
    }    
    
    for( k=0; k<SLOPES_LENGTH; k++ ){
        for( l=0; l<maxTerm; l++ ){
            (*dx)(k,l) = dx_int(k,l);
            (*dy)(k,l) = dy_int(k,l);
        }
    }
    
    MatrixXd A = MatrixXd::Zero(SLOPES_LENGTH*2,ZERNIKE_CNT);
    A.block(0,0,SLOPES_LENGTH,ZERNIKE_CNT) = *dx;
    A.block(SLOPES_LENGTH,0,SLOPES_LENGTH,ZERNIKE_CNT) = *dy;

	JacobiSVD<MatrixXd> svd(A, ComputeThinU | ComputeThinV);
    
	MatrixXd w = svd.singularValues().asDiagonal();
	MatrixXd pseudoInvW( ZERNIKE_CNT, ZERNIKE_CNT );
	pseudoInverse( w, pseudoInvW, 1e-9);
    
    MatrixXd res_int = svd.matrixV()*pseudoInvW*(svd.matrixU().transpose());
    
    for( k=0; k<ZERNIKE_CNT; k++ ){
        for( l=0; l<SLOPES_LENGTH*2; l++ ){
            (*res)(k,l) = res_int(k,l);
        }
    }
    
    //printf( "voila! \n" );
    //std::cout << dx_int;
    //std::cout << *dx;
}


#define MDL_CHECK_PARAMETERS
#if defined(MDL_CHECK_PARAMETERS) && defined(MATLAB_MEX_FILE)
/* Function: mdlCheckParameters =============================================
 * Abstract:
 *    Validate our parameters to verify they are valid.
 */
static void mdlCheckParameters(SimStruct *S)
{
    uint_T k, i;
    
	/* Check vector */
	if( mxIsInt32(SLOPES_PARAM) ){
		ssSetErrorStatus(S,"parameter to S-function "
							"\"SLOPES_PARAM\" must be int32");
		return;
	}
    
    if( ZERNIKE_CNT > 28 ){
    	ssSetErrorStatus(S,"parameter to S-function "
							"\"ZERNIKE_CNT\" must be smaller or equal 28");
		return;
    }

    /* Check vector */
	if( mxIsInt32(DECIMAL_PARAM) ){
		ssSetErrorStatus(S,"parameter to S-function "
							"\"DECIMAL_PARAM\" must be int32");
		return;
	}

    /* Check vector */
	if( mxIsInt32(ZERNIKE_PARAM) ){
		ssSetErrorStatus(S,"parameter to S-function "
							"\"ZERNIKE_PARAM\" must be int32");
		return;
	}
}
#endif /* MDL_CHECK_PARAMETERS */

static void mdlInitializeSizes(SimStruct *S)
{
	uint_T i;
	
  ssSetNumSFcnParams(S, NUMBER_OF_PARAMS);
#if defined(MATLAB_MEX_FILE)
  if (ssGetNumSFcnParams(S) == ssGetSFcnParamsCount(S)) {
    mdlCheckParameters(S);
    if (ssGetErrorStatus(S) != NULL) {
      return;
    }
  } else {
    return;
  }
#endif

	for (i = 0; i < NUMBER_OF_PARAMS; i++) {
		ssSetSFcnParamNotTunable(S, i);
	}
	
	ssSetNumContStates(S, 0);
	ssSetNumDiscStates(S, 0);

	if (!ssSetNumInputPorts(S, 2)) return;
    
	ssSetInputPortWidth(S, 0, SLOPES_LENGTH);
	ssSetInputPortDirectFeedThrough(S, 0, 1);
    ssSetInputPortRequiredContiguous(S, 0, 1);
	ssSetInputPortDataType( S, 0, SS_DOUBLE );
    
    ssSetInputPortWidth(S, 1, SLOPES_LENGTH);
	ssSetInputPortDirectFeedThrough(S, 1, 1);
    ssSetInputPortRequiredContiguous(S, 1, 1);
	ssSetInputPortDataType( S, 1, SS_DOUBLE );

	if (!ssSetNumOutputPorts(S, 1)) return;
	ssSetOutputPortWidth(S, 0, ZERNIKE_CNT);

	ssSetNumSampleTimes(S, 1);
	ssSetNumRWork(S, 0);
	ssSetNumIWork(S, 0);
	ssSetNumPWork(S, 0);
	ssSetNumModes(S, 0);

	/* Take care when specifying exception free code - see sfuntmpl.doc */
	ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);
}

static void mdlInitializeSampleTimes(SimStruct *S)
{
	ssSetSampleTime(S, 0, -1);
	ssSetOffsetTime(S, 0, 0.0);
	ssSetModelReferenceSampleTimeDefaultInheritance(S);
}

#define MDL_INITIALIZE_CONDITIONS
/* Function: mdlInitializeConditions ========================================
 * Abstract:
 *    
 */
static void mdlInitializeConditions(SimStruct *S)
{
	arrayofslopefunctionsZern( S, 1.0f, ZERNIKE_CNT );	
}

/* Function: mdlOutputs =======================================================
 * Abstract:
 *      y = Cx + Du
 */
static void mdlOutputs(SimStruct *S, int_T tid)
{
    real_T *y       = ssGetOutputPortRealSignal(S,0);
    
    double *slopeX_d = (double*)ssGetInputPortRealSignal(S,0);
    double *slopeY_d = (double*)ssGetInputPortRealSignal(S,1);
    
#ifdef MEASURE_TIME
    int32_t timeS;
    #if defined(MATLAB_MEX_FILE)
        clock_t start, end;
        double cpu_time_used;
        start = clock();
    #else
        int32_t RTTSKinit = 0;
        int32_t RTTSKend = 0;
        RTTSKinit = rt_get_cpu_time_ns();
    #endif
#endif
       
    MatrixXd slopeX = Map<MatrixXd>(slopeX_d, SLOPES_LENGTH, 1);
    MatrixXd slopeY = Map<MatrixXd>(slopeY_d, SLOPES_LENGTH, 1);
    
    MatrixXd slope = MatrixXd::Zero(SLOPES_LENGTH*2,1);
    slope.block(0,0,SLOPES_LENGTH,1) = slopeX;
    slope.block(SLOPES_LENGTH,0,SLOPES_LENGTH,1) = slopeY;

	VectorXd zernikes = (*res)*slope;

#ifdef MEASURE_TIME
    #if defined(MATLAB_MEX_FILE)
        end = clock();
        cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC*1000000;
        printf("time consumption %0.3f (us)\n", cpu_time_used);
    #else
        RTTSKend = rt_get_cpu_time_ns();
        timeS = (RTTSKend-RTTSKinit);
        printf("time consumption %i (ns)\n", timeS);
    #endif
    
#endif
      
    //std::cout << "rows: " << A.rows() << "cols: " << A.cols() << "------------" << std::endl;
    //std::cout << "-------"<< zernikes << std::endl << std::endl;
    memcpy( y,zernikes.data(),ZERNIKE_CNT*sizeof(double));
}

static void mdlTerminate(SimStruct *S){
	UNUSED_ARG(S); /* unused input argument */
    
    delete(dx);
	delete(dy);
    delete(res);
            
}

#ifdef MATLAB_MEX_FILE /* Is this file being compiled as a MEX-file? */
#include "simulink.c" /* MEX-file interface mechanism */
#else
#include "cg_sfun.h" /* Code generation registration function */
#endif
