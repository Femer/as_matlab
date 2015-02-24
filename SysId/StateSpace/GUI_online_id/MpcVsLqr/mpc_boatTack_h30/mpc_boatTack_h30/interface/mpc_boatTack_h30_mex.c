/*
mpc_boatTack_h30 : A fast customized optimization solver.

Copyright (C) 2013-2015 EMBOTECH GMBH [info@embotech.com]. All rights reserved.


This software is intended for simulation and testing purposes only. 
Use of this software for any commercial purpose is prohibited.

This program is distributed in the hope that it will be useful.
EMBOTECH makes NO WARRANTIES with respect to the use of the software 
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
PARTICULAR PURPOSE. 

EMBOTECH shall not have any liability for any damage arising from the use
of the software.

This Agreement shall exclusively be governed by and interpreted in 
accordance with the laws of Switzerland, excluding its principles
of conflict of laws. The Courts of Zurich-City shall have exclusive 
jurisdiction in case of any dispute.

*/

#include "mex.h"
#include "math.h"
#include "../include/mpc_boatTack_h30.h"

/* copy functions */
void copyCArrayToM(float *src, double *dest, int dim) {
    while (dim--) {
        *dest++ = (double)*src++;
    }
}
void copyMArrayToC(double *src, float *dest, int dim) {
    while (dim--) {
        *dest++ = (float) (*src++) ;
    }
}

#if mpc_boatTack_h30_SET_PRINTLEVEL > 0
	FILE *fp = NULL;
#endif

/* Some memory for mex-function */
mpc_boatTack_h30_params params;
mpc_boatTack_h30_output output;
mpc_boatTack_h30_info info;

/* THE mex-function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )  
{
	/* define variables */	
	mxArray *par;
	mxArray *outvar;
	const mxArray *PARAMS = prhs[0];
	double *pvalue;
	int i;
	int exitflag;
	const char *fname;
	const char *outputnames[1] = {"u0"};
	const char *infofields[15] = { "it", "res_eq", "res_ineq",  "pobj",  "dobj",  "dgap", "rdgap",  "mu",  "mu_aff",  "sigma",  "lsit_aff",  "lsit_cc",  "step_aff",   "step_cc",  "solvetime"};
	
	/* Check for proper number of arguments */
    if (nrhs != 1) {
        mexErrMsgTxt("This function requires exactly 1 input: PARAMS struct.\nType 'help mpc_boatTack_h30_mex' for details.");
    }    
	if (nlhs > 3) {
        mexErrMsgTxt("This function returns at most 3 outputs.\nType 'help mpc_boatTack_h30_mex' for details.");
    }

	/* Check whether params is actually a structure */
	if( !mxIsStruct(PARAMS) ) {
		mexErrMsgTxt("PARAMS must be a structure.");
	}

	/* copy parameters into the right location */
	par = mxGetField(PARAMS, 0, "minusAExt_times_x0");
#ifdef MEXARGMUENTCHECKS
    if( par == NULL )	{
        mexErrMsgTxt("PARAMS.minusAExt_times_x0 not found");
    }
    if( !mxIsDouble(par) )
    {
    mexErrMsgTxt("PARAMS.minusAExt_times_x0 must be a double.");
    }
    if( mxGetM(par) != 3 || mxGetN(par) != 1 ) {
    mexErrMsgTxt("PARAMS.minusAExt_times_x0 must be of size [3 x 1]");
    }
#endif	 
    copyMArrayToC(mxGetPr(par), params.minusAExt_times_x0, 3);

	par = mxGetField(PARAMS, 0, "Hessians");
#ifdef MEXARGMUENTCHECKS
    if( par == NULL )	{
        mexErrMsgTxt("PARAMS.Hessians not found");
    }
    if( !mxIsDouble(par) )
    {
    mexErrMsgTxt("PARAMS.Hessians must be a double.");
    }
    if( mxGetM(par) != 4 || mxGetN(par) != 1 ) {
    mexErrMsgTxt("PARAMS.Hessians must be of size [4 x 1]");
    }
#endif	 
    copyMArrayToC(mxGetPr(par), params.Hessians, 4);

	par = mxGetField(PARAMS, 0, "HessiansFinal");
#ifdef MEXARGMUENTCHECKS
    if( par == NULL )	{
        mexErrMsgTxt("PARAMS.HessiansFinal not found");
    }
    if( !mxIsDouble(par) )
    {
    mexErrMsgTxt("PARAMS.HessiansFinal must be a double.");
    }
    if( mxGetM(par) != 4 || mxGetN(par) != 4 ) {
    mexErrMsgTxt("PARAMS.HessiansFinal must be of size [4 x 4]");
    }
#endif	 
    copyMArrayToC(mxGetPr(par), params.HessiansFinal, 16);

	par = mxGetField(PARAMS, 0, "lowerBound");
#ifdef MEXARGMUENTCHECKS
    if( par == NULL )	{
        mexErrMsgTxt("PARAMS.lowerBound not found");
    }
    if( !mxIsDouble(par) )
    {
    mexErrMsgTxt("PARAMS.lowerBound must be a double.");
    }
    if( mxGetM(par) != 2 || mxGetN(par) != 1 ) {
    mexErrMsgTxt("PARAMS.lowerBound must be of size [2 x 1]");
    }
#endif	 
    copyMArrayToC(mxGetPr(par), params.lowerBound, 2);

	par = mxGetField(PARAMS, 0, "upperBound");
#ifdef MEXARGMUENTCHECKS
    if( par == NULL )	{
        mexErrMsgTxt("PARAMS.upperBound not found");
    }
    if( !mxIsDouble(par) )
    {
    mexErrMsgTxt("PARAMS.upperBound must be a double.");
    }
    if( mxGetM(par) != 2 || mxGetN(par) != 1 ) {
    mexErrMsgTxt("PARAMS.upperBound must be of size [2 x 1]");
    }
#endif	 
    copyMArrayToC(mxGetPr(par), params.upperBound, 2);

	par = mxGetField(PARAMS, 0, "C");
#ifdef MEXARGMUENTCHECKS
    if( par == NULL )	{
        mexErrMsgTxt("PARAMS.C not found");
    }
    if( !mxIsDouble(par) )
    {
    mexErrMsgTxt("PARAMS.C must be a double.");
    }
    if( mxGetM(par) != 3 || mxGetN(par) != 4 ) {
    mexErrMsgTxt("PARAMS.C must be of size [3 x 4]");
    }
#endif	 
    copyMArrayToC(mxGetPr(par), params.C, 12);

	par = mxGetField(PARAMS, 0, "D");
#ifdef MEXARGMUENTCHECKS
    if( par == NULL )	{
        mexErrMsgTxt("PARAMS.D not found");
    }
    if( !mxIsDouble(par) )
    {
    mexErrMsgTxt("PARAMS.D must be a double.");
    }
    if( mxGetM(par) != 3 || mxGetN(par) != 4 ) {
    mexErrMsgTxt("PARAMS.D must be of size [3 x 4]");
    }
#endif	 
    copyMArrayToC(mxGetPr(par), params.D, 12);

	#if mpc_boatTack_h30_SET_PRINTLEVEL > 0
		/* Prepare file for printfs */
		fp = freopen("stdout_temp","w+",stdout);
		if( fp == NULL ) {
			mexErrMsgTxt("freopen of stdout did not work.");
		}
		rewind(fp);
	#endif

	/* call solver */
	exitflag = mpc_boatTack_h30_solve(&params, &output, &info);
	
	#if mpc_boatTack_h30_SET_PRINTLEVEL > 0
		/* Read contents of printfs printed to file */
		rewind(fp);
		while( (i = fgetc(fp)) != EOF ) {
			mexPrintf("%c",i);
		}
	#endif

	/* copy output to matlab arrays */
	plhs[0] = mxCreateStructMatrix(1, 1, 1, outputnames);
	outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
	copyCArrayToM( output.u0, mxGetPr(outvar), 1);
	mxSetField(plhs[0], 0, "u0", outvar);	

	/* copy exitflag */
	if( nlhs > 1 )
	{
		plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(plhs[1]) = (double)exitflag;
	}

	/* copy info struct */
	if( nlhs > 2 )
	{
		        plhs[2] = mxCreateStructMatrix(1, 1, 15, infofields);
         
		
		/* iterations */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = (double)info.it;
		mxSetField(plhs[2], 0, "it", outvar);
		
		/* res_eq */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = info.res_eq;
		mxSetField(plhs[2], 0, "res_eq", outvar);

		/* res_ineq */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = info.res_ineq;
		mxSetField(plhs[2], 0, "res_ineq", outvar);

		/* pobj */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = info.pobj;
		mxSetField(plhs[2], 0, "pobj", outvar);

		/* dobj */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = info.dobj;
		mxSetField(plhs[2], 0, "dobj", outvar);

		/* dgap */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = info.dgap;
		mxSetField(plhs[2], 0, "dgap", outvar);

		/* rdgap */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = info.rdgap;
		mxSetField(plhs[2], 0, "rdgap", outvar);

		/* mu */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = info.mu;
		mxSetField(plhs[2], 0, "mu", outvar);

		/* mu_aff */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = info.mu_aff;
		mxSetField(plhs[2], 0, "mu_aff", outvar);

		/* sigma */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = info.sigma;
		mxSetField(plhs[2], 0, "sigma", outvar);

		/* lsit_aff */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = (double)info.lsit_aff;
		mxSetField(plhs[2], 0, "lsit_aff", outvar);

		/* lsit_cc */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = (double)info.lsit_cc;
		mxSetField(plhs[2], 0, "lsit_cc", outvar);

		/* step_aff */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = info.step_aff;
		mxSetField(plhs[2], 0, "step_aff", outvar);

		/* step_cc */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = info.step_cc;
		mxSetField(plhs[2], 0, "step_cc", outvar);

		/* solver time */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = info.solvetime;
		mxSetField(plhs[2], 0, "solvetime", outvar);
	}
}