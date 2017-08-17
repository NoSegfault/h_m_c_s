#include <math.h>
#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[]) {
  // Test number of parameters.
  if (nrhs != 3 || nlhs != 1) {
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
  int bottomK = mxGetScalar(prhs[2]);  
  
  plhs[0] = mxCreateDoubleMatrix(bottomK,1,mxREAL); /* create the output vector */
  double* cuts = mxGetPr(plhs[0]);

  int i,j; 

  int* minIndex=new int[numEdges];
  int* maxIndex=new int[numEdges];
  
  for(j=0; j<numEdges; j++) 
  {   
     minIndex[j]=(int)irsTrans[jcsTrans[j]];
     // minIndex[j] = max(num-bottomK, minIndex[j]); // C indexing from num-bottomK to num-1
     if(minIndex[j] < num-bottomK) minIndex[j] = num-bottomK;
          
	 maxIndex[j]=(int)irsTrans[jcsTrans[j+1]-1];     
  }

  for(j=0; j<bottomK; j++)
   cuts[j]=0;

  for(j=0; j<numEdges; j++)
  {
    for(i=minIndex[j]; i<maxIndex[j]; i++)
	 cuts[i-num+bottomK]+=wval[j];
  }
  delete minIndex; delete maxIndex;
}

