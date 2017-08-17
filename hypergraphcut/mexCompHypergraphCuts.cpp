#include <math.h>
#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[]) {
  // Test number of parameters.
  if (nrhs != 2 || nlhs != 1) {
    mexWarnMsgTxt("Usage: cuts=mexCompHypergraphCuts(W,INC'); expects INC to be sorted in increasing order of f");
    return;
  }
  
  // get important parameters
  int num       = (int)mxGetM(prhs[1]); // transpose of INC (number of rows = number of vertices)
  int numEdges  = (int)mxGetN(prhs[1]); // transpose of INC (number of columns = number of edges)
  
  if(!mxIsSparse(prhs[1])) { mexWarnMsgTxt("Matrix is not sparse"); return;}
  
  // Create output array and compute values
  double* wval = mxGetPr(prhs[0]);   // get edge weights
  mwIndex* irsTrans = mxGetIr(prhs[1]);   // get rows of INC
  mwIndex* jcsTrans = mxGetJc(prhs[1]);   // get cols of INC
  
  plhs[0] = mxCreateDoubleMatrix(num,1,mxREAL); /* create the output vector */
  double* cuts = mxGetPr(plhs[0]);

  int i,j; 

  int* minIndex=new int[numEdges];
  int* maxIndex=new int[numEdges];
  
  for(j=0; j<numEdges; j++) 
  {   
     minIndex[j]=(int)irsTrans[jcsTrans[j]];
	 maxIndex[j]=(int)irsTrans[jcsTrans[j+1]-1];
  }

  for(j=0; j<num; j++)
   cuts[j]=0;

  for(j=0; j<numEdges; j++)
  {
    for(i=minIndex[j]; i<maxIndex[j]; i++)
	 cuts[i]+=wval[j];
  }
  delete minIndex; delete maxIndex;
}

