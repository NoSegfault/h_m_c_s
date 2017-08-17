// output: R(f) = sum_{e in E} w(e) (max_{i in e} f_i - min_{i in e} f_i)

#include <math.h>
#include "mex.h"
#include "matrix.h"

double* ProjAlphaSimplex(double*, int, double);
double* ProjKSimplex(double*, int,int);

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[]) {
  // Test number of parameters.
  if (nrhs != 3 || nlhs != 1) {
    mexWarnMsgTxt("Usage: lambda=mexEvalHypergraphObjVarTransform(W,INC',f)");
    return;
  }
  
  // get important parameters
  int num       = (int)mxGetM(prhs[1]); // INC (number of rows)
  int numEdges  = (int)mxGetN(prhs[1]); // INC (number of columns)
  
  if(!mxIsSparse(prhs[1])) { mexWarnMsgTxt("Matrix is not sparse");}
  
  // Create output array and compute valuesfold
  double* wval = mxGetPr(prhs[0]);   // get edge weights w
  mwIndex* irsTrans = mxGetIr(prhs[1]);   // get rows of INC
  mwIndex* jcsTrans = mxGetJc(prhs[1]);   // get cols of INC
  double* f = mxGetPr(prhs[2]);     // get f
  
  plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL); /* create the output vector */
  double* X = mxGetPr(plhs[0]);
 

  int counter=0,i,j;
  double primalVal=0.0,maxVec,minVec;
  
  //  we have to compute max/min of components of f with indices corresponding
  //  to hyperedges
  
  //  we iterate over the cols of INC
  for(j=0; j<numEdges; j++) 
  {   
	 maxVec=-1E15; minVec=1E15;
	 //mexPrintf("Counter %i\n",counter);
	 for(i=0; i<jcsTrans[j+1]-jcsTrans[j]; i++) // address row of xu; there are jcsTrans[j+1]-jcsTrans[j] elements in col. j
	 {  
		// computation of primal value (maximum and minimum on edge)
		// irsTrans[counter] is one index of an vertex in hyperedge j
        if(f[irsTrans[counter]]>maxVec)
	       maxVec=f[irsTrans[counter]];
		if(f[irsTrans[counter]]<minVec)
		    minVec=f[irsTrans[counter]];
		counter++;
	 }	  
	 primalVal += wval[j]*(maxVec-minVec); // note that f=-vec/norm(vec)
  }
  X[0]=primalVal;
}

