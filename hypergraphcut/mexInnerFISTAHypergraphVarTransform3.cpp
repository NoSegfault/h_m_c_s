// (C)2012-13 Syama Sundar Rangapuram, Matthias Hein and Leonardo Jost 
// Machine Learning Group, Saarland University, Germany
// http://www.ml.uni-saarland.de
 
#include <math.h>
#include <stdio.h>
#include <time.h>

#include "mex.h"
#include "matrix.h"

double* ProjAlphaSimplex(double*, int, double);
double* ProjKSimplex(double*, int,int);
void    ProjectAll(double**, double*,int*, double*, mwIndex*, int,int);

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[]) {
  // Test number of parameters.
  if (nrhs != 11 || nlhs != 5) {
    mexWarnMsgTxt("Usage: [X,u,v,Obj,iter]=mexInnerFISTAHypergraphVarTransform(W,INC,INC',IndexList,Y,u,v,MAXITER,EPS,Lip,verbosity)");
    return;
  }
  
  // get important parameters
  int numEdges = (int)mxGetM(prhs[1]); // INC (number of rows)
  int num      = (int)mxGetN(prhs[1]); // INC (number of columns)
  int lenyu    = (int)mxGetM(prhs[3]); // rval
  
  if(!mxIsSparse(prhs[1])) { mexWarnMsgTxt("Matrix is not sparse");}
  
    
  // Create output array and compute values
  double* wval = mxGetPr(prhs[0]);   // get edge weights
  mwIndex* irs = mxGetIr(prhs[1]);   // get rows of INC
  mwIndex* jcs = mxGetJc(prhs[1]);   // get cols of INC
  mwIndex* irsTrans = mxGetIr(prhs[2]);   // get rows of INC
  mwIndex* jcsTrans = mxGetJc(prhs[2]);   // get cols of INC
  double* IxxList = mxGetPr(prhs[3]);

  int* Ixx = new int[lenyu];
  
  double* Y = mxGetPr(prhs[4]);      // linear term
  double* yu = mxGetPr(prhs[5]); // get old values for u
  double* yv = mxGetPr(prhs[6]); // get old values for v

  int MAXITER = mxGetScalar(prhs[7]); 
  double EPS = mxGetScalar(prhs[8]); 
  double Lip = mxGetScalar(prhs[9]); 
  double verb= mxGetScalar(prhs[10]); 
  //mexPrintf("Lip: %f, lenyu: %i, numEdges: %i, num: %i\n",Lip,lenyu,numEdges,num);

  if(Lip<=0){
	  mexWarnMsgTxt("Lipschitz constant has to be positive");
    return;
  }

  double *X;      /* output matrix */
  plhs[0] = mxCreateDoubleMatrix(num,1,mxREAL); /* create the output vector */
  plhs[1] = mxCreateDoubleMatrix(lenyu,1,mxREAL); /* create the output dual variable */
  plhs[2] = mxCreateDoubleMatrix(lenyu,1,mxREAL); /* create the output dual variable */
  plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL); /* create the output objective value */
  plhs[4] = mxCreateDoubleMatrix(1,1,mxREAL); /* create the final iteration value */
  
  //plhs[1]= (mxArray *)prhs[2];
  //plhs[2]= (mxArray *)prhs[3];

  X = mxGetPr(plhs[0]);
  double* U = mxGetPr(plhs[1]);
  double* V = mxGetPr(plhs[2]);
  double* OutputObj = mxGetPr(plhs[3]);
  double* FinalIter = mxGetPr(plhs[4]);

  int counter=0,i,j;
  int iter=0; int len;
  double tnew=1; double told=1,adummy,bdummy,factor;
  double dummy,dualVal,primalVal,diff,maxVec,minVec;
  double* dummyPointer;

  double* vec     = new double[num];
  double* vecold  = new double[num];
  double* vecBest = new double[num];
  for(i=0; i<num; i++) { vec[i]=0; vecold[i]=0; vecBest[i]=0;}

  double* xu    = new double[lenyu];
  double* xuold = new double[lenyu];
  for(i=0; i<lenyu; i++) { xu[i]=0;  xuold[i]=0;}

  double* xv    = new double[lenyu];
  double* xvold = new double[lenyu];
  for(i=0; i<lenyu; i++) { xv[i]=0;  xvold[i]=0; }

  // copy double list into integer
  for(i=0; i<lenyu; i++) { Ixx[i]=(int)(IxxList[i]);}

  double** IXXDUMMY=new double*[num];
  
  //MaxSumSquaredWeights=max(sum(W.^2,2));
  /*double MaxSumSquaredWeights=0;  
  for(j=0; j<cols; j++) 
  {   
	dummy=0;
	for(i=0; i<jcs[j+1]-jcs[j]; i++) {  dummy+=wval[j]*wval[j]; }
	if(dummy>MaxSumSquaredWeights) { MaxSumSquaredWeights=dummy; }
  }
  MaxSumSquaredWeights=4*MaxSumSquaredWeights;*/

  dualVal=EPS+1; primalVal=EPS; double primalMin=0;

  clock_t proj_start, proj_end, proj_total,primal_start,primal_total,total_start;
  proj_total=0; primal_total=0;
  total_start = clock();

  primalVal=1000;
  
  while(iter<MAXITER && abs(dualVal+primalMin)/dualVal > EPS && dualVal>1E-8)
  {
    // exchange vec and vecold
	dummyPointer=vec; vec=vecold; vecold=dummyPointer;

	//mexPrintf("Exchanged vec \n");

	// exchange xu/xv and xuold/xvold
	dummyPointer=xu; xu=xuold; xuold=dummyPointer;
	dummyPointer=xv; xv=xvold; xvold=dummyPointer;

	//mexPrintf("Exchanged Pointer \n");

	// exchange tnew and told
	told=tnew;
    
	// initialize X with zeros 
    //for(i=0; i<len; i++) { X[i]=0; }

	//mexPrintf("Initialized X, %i \n",X[0]);
  
	//vec = sum((yu - yv),1)' - Y;
    dummy=0; counter=0; dualVal=0;
    for(j=0; j<num; j++) 
    {   
	   dummy=0;
       len = jcs[j+1]-jcs[j];
	   for(i=0; i<len; i++)
	   {  
		 dummy+= yu[counter]-yv[counter];
         counter++;
	   }
	   dummy-=Y[j];
	   vec[j]=dummy;
	   dualVal+=dummy*dummy;
	}
	// split of sums in order to be local in memory
	counter=0;
	for(j=0; j<num; j++) 
    {   
	   dummy=vec[j];
       len = jcs[j+1]-jcs[j];
       for(i=0; i<len; i++)
	   {  
	      // update of xu and xv         
	      xu[counter] = yu[counter] - 1/Lip*dummy;
		  //xv[counter] = yv[counter] + 1/Lip*dummy;	  
		  counter++;
	   }	
	} 
	counter=0;
	for(j=0; j<num; j++) 
    {   
	   dummy=vec[j];
       len = jcs[j+1]-jcs[j];
       for(i=0; i<len; i++)
	   {  
	      // update of xu and xv         
	      //xu[counter] = yu[counter] - 1/Lip*dummy;
		  xv[counter] = yv[counter] + 1/Lip*dummy;	  
		  counter++;
	   }	
	} 
	dualVal=sqrt(dualVal); 

	//mexPrintf("Computed X, %f %f %f %f %f \n",X[0],X[1],X[2],X[3],X[4]);
    // computation of dual value
    /*dualVal=0;
    for(i=0; i<num; i++) 
	{ 
	  dualVal+=vec[i]*vec[i]; 
    }
	dualVal=sqrt(dualVal); */

    
	//xu = yu - 1/L*gradvec_u;  % Nemirovski-update  
    //xv = yv - 1/L*gradvec_v;
	/*counter=0;
    for(j=0; j<num; j++) 
    {   
	   for(i=0; i<jcs[j+1]-jcs[j]; i++)
	   {  
	      // update of xu and xv         
	      xu[counter] = yu[counter] - 1/Lip*vec[j];
		  xv[counter] = yv[counter] + 1/Lip*vec[j];	  
		  counter++;
	   }	  
    }*/
	//compute primal value
	if(iter % 5 ==0)
	{
	primal_start = clock();
    primalVal=0; counter=0;
    for(j=0; j<numEdges; j++) 
    {   
	   maxVec=-1E15; minVec=1E15;
	   len=jcsTrans[j+1]-jcsTrans[j];
	   for(i=0; i<len; i++) // address row of xu
	   {
	     dummy=vec[irsTrans[counter]];
	     if(dummy>maxVec)
	       maxVec=dummy;
		 if(dummy<minVec)
		    minVec=dummy;
		 counter++;
	   }	
	   primalVal -= wval[j]*(minVec-maxVec);
	}
    primal_total+=clock()-primal_start;
	}
    //project each row of xu and xv onto the simplex 
    proj_start=clock();
	ProjectAll(IXXDUMMY,xu,Ixx,wval,jcsTrans,numEdges,num);
	ProjectAll(IXXDUMMY,xv,Ixx,wval,jcsTrans,numEdges,num);
	proj_end = clock();
	proj_total +=proj_end-proj_start;

	
	//yu = xu + i/(i+3)*(xu - xu_old);  factor used below is from FISTA paper
    //yv = xv + i/(i+3)*(xv - xv_old);
	counter=0;
	tnew = (1 + sqrt(1+4*told*told))/2;
    factor = (told-1)/tnew;
    // split of the for loops in order to stay local in memory
	for(j=0; j<num; j++) 
    {   
	   len = jcs[j+1]-jcs[j];
	   for(i=0; i<len; i++)
	   {  
	      // update of xu and xv
	      dummy=xu[counter];
	      yu[counter] = dummy + factor*(dummy-xuold[counter]);
		  counter++;
	   }
	}
	counter=0;
    for(j=0; j<num; j++) 
    {   
	   len = jcs[j+1]-jcs[j];
	   for(i=0; i<len; i++)
	   {  
		  dummy=xv[counter];
		  yv[counter] = dummy + factor*(dummy-xvold[counter]);	  
		  counter++;
	   }	  
    }
    if(iter % 5 ==0)
	{
    primal_start = clock();
	// computation of primal value (continued) - note that f=-vec/norm(vec) !
	for(i=0;i<num;i++)
     primalVal+=Y[i]*vec[i];      // therefore + here
	primalVal=primalVal/dualVal;  // divide by the norm of vec 

	if(primalVal<primalMin)
	{
      primalMin=primalVal;
	  for(i=0; i<num; i++)
       vecBest[i]=vec[i];
	}
	 primal_total+=clock()-primal_start; 
	}
  
    //mexPrintf("Iteration: %i, Fval: %1.15f, RelativeChange %1.15f\n",iter,Fval,relativeChange);
	if(verb>0 && (iter<10 || iter % 100==0))
	 mexPrintf("Iteration: %i, DualValue: %1.15f, PrimalValue %1.15f\n",iter,-dualVal,primalVal);
	iter++;
 }
 double total_time= ((double)(clock() - total_start))/((double)CLOCKS_PER_SEC);
 if(verb>0)
 { 
   mexPrintf("Total time: %f, Time for Projections: %f, Time for primal value: %f\n",total_time,((double)proj_total)/((double)CLOCKS_PER_SEC),((double)primal_total)/((double)CLOCKS_PER_SEC));
   mexPrintf("FINAL: Iterations %i, DualValue: %1.15f, PrimalValue %1.15f\n",iter,-dualVal,primalMin);
 }
 if(primalMin<0)
   for(i=0; i<num; i++) { X[i]=-vecBest[i];} 
 else
   for(i=0; i<num; i++) { X[i]=-vec[i];} 
 for(i=0; i<lenyu; i++) { U[i]=yu[i]; V[i]=yv[i];}
 OutputObj[0]=dualVal;
 FinalIter[0]=iter;
   
 delete vec; delete vecold; delete vecBest; delete xuold; delete xu;
 delete xvold; delete xv;
 delete IXXDUMMY;
}

void ProjectAll(double** IXX, double* xv,int* Indices, double* wval, mwIndex* jcsTrans, int numEdges, int num)
{
  double* ptr; 
  //double* X =new double[num];
  double suma,alpha,nablak,val,t;
  int OUTERcounter=0,i,j,len;
  int counter,counter1,counter2,iteration;
  int check=0;

  for(j=0; j<numEdges; j++) 
  {   
     suma=0;
	 // we go through the rows of xu resp xv
	 len = (jcsTrans[j+1]-jcsTrans[j]);
	 
     // projection onto the wval(j)-simplex that is u_i>=0 and sum_i u_i = wval(j)
     alpha = wval[j];
	 if(len==2)
	 {
	   val = 0.5*(alpha + xv[Indices[OUTERcounter]]- xv[Indices[OUTERcounter+1]]);
	   if(val<0)
        val=0;
	   else if(val>alpha)
	   {
	     val=alpha;
	   }
	   xv[Indices[OUTERcounter]]=val; xv[Indices[OUTERcounter+1]]=alpha-val;
	   OUTERcounter+=2;
	 }
	 else
	 {
	   for(i=0; i<len; i++) // address row of xu
	   {  
	      ptr = &(xv[Indices[OUTERcounter]]);
		  //X[i]=*ptr;//xv[Indices[OUTERcounter]];
		  //ptr = &(X[i]);
	      IXX[i]=ptr;
		  suma +=*ptr;
		  OUTERcounter++;
	   }	  
	   counter = len;  counter1=0; counter2=0;
       
	   iteration=0; 
       t=(suma - alpha)/((double)counter);  
       //mexPrintf("Initial sum: %1.15f, alpha: %1.15f, t: %1.15f\n",suma,alpha,t);

       while(iteration<1000)
       {
         counter1=0; counter2=0;
	     nablak=0; 
         for(i=0;i<counter; i++)
	     {
	       //index=IX[i];
	       //val = copy[index];// - t;
	       ptr= IXX[i];
	       val = *ptr;

	       if(val<=t) // corresponds to x<=0
	        { nablak -= val; counter2++; *ptr=0; }
	       else 
	       {
		    //IX[counter1]=index; counter1++;
		    IXX[counter1]=ptr; counter1++;
	       }
	     } 
	     if(counter2==0)//nablak==0) 
		   break; //iteration=1000; }//mexPrintf("BREAK\n"); 
	     else
	     {
	       counter=counter1;
           t += nablak/((double)counter1) + t*( ((double)counter2) / ((double)counter1));
	     }
	     //mexPrintf("%i %i %1.15f %i %1.15f\n",iteration,counter1,nablak,counter,t);
	     iteration++;
       }  
       suma=0;
       for(i=0; i<counter; i++)
       {
	     *(IXX[i])-=t;
       }
	 }
  }
  //delete IXX;
}


double* ProjAlphaSimplex(double* f,int len, double alpha)
{
  if(len==2)
  {
    double a;
    a = 0.5*(alpha + f[0]-f[1]);
	if(a<0)
     a=0;
	else if(a>alpha)
	{
	  a=alpha;
	}
	f[0]=a; f[1]=alpha-a;
	return f;
  }
  else
  {
  int i,j,index;
  
  double* copy = new double[len]; int* dummyPtr;
  int* IX = new int[len]; double counter=len;
  int* IX1= new int[len]; int counter1=0;
  int* IX2= new int[len]; int counter2=0;
  double ucounter=0,ucounterold=0;
  double suma = 0; double suma1; double suma2;
  
  for(i=0; i<len; i++)
   { IX[i]=i; suma+=f[i]; copy[i]=f[i]; }
  int iteration=0;
  double rk,nablak,deltak,val,aval,t=0;
  
  int MAXITER=10000;
  //return f;
  while(iteration<MAXITER)
  {
    counter1=0; counter2=0; suma1=0; suma2=0;
	nablak=0; deltak=0;
    rk = alpha*(1 - ucounter);
	t = (suma - rk)/counter;
	
    ucounterold=ucounter; 
    for(i=0;i<counter; i++)
	{
      index=IX[i];
	  aval= copy[index];
	  val = aval - t;  
	  
	  if(val<=0) // corresponds to x<=0
	  {
        IX2[counter2]=index; counter2++;
		suma2+=aval;
		//f[index]=0;
		nablak -= val;
	  }
	  else 
	  { 
		  if(val>=alpha)
	      {
             IX1[counter1]=index; counter1++;
			 suma1+=aval;
			 //f[index]=1; 
			 ucounter++;
			 deltak += val-alpha;
	      }
		  else
		  {
            IX2[counter2]=index; counter2++;
			IX1[counter1]=index; counter1++;
			suma1+=aval; suma2+=aval;
			//f[index]=val;
		  }
	  }
	}
	if(nablak==deltak) { iteration=MAXITER;}
	else
	{
	   if(nablak > deltak)
	   { 
		 dummyPtr=IX; IX=IX1; IX1=dummyPtr; 
		 counter=(double)counter1; suma=suma1; ucounter=ucounterold;}
		else
        { 		
		  dummyPtr=IX; IX=IX2; IX2=dummyPtr; 
		  counter=(double)counter2; suma=suma2;}
	}
	//mexPrintf("%i %i %i %f %f %f %f %i\n",iteration,counter1,counter2,nablak,deltak,counter,suma,ucounter);
	iteration++;
  }
  for(i=0;i<len;i++)
  {
	 val=copy[i]-t;
	 if(val<0)
	 {
	   f[i]=0;
	 }
	 else if(val>alpha)
	 {
       f[i]=alpha;
	 
	 }
	 else
		 f[i]=val;
  }
  // check
  //double minf=1E15; double maxf=-10;
  //suma=0; for(j=0;j<len; j++){suma+=f[j]; if(f[j]<minf){minf=f[j]; } if(f[j]>maxf){maxf=f[j]; }}
  //mexPrintf("Computed X, %i, %f %f %f %f %f \n",len,f[0],f[1],f[2],f[3],f[4]);
  if(iteration==MAXITER)
    { suma=0; for(j=0;j<len; j++){suma+=f[j];} mexPrintf("Iter %i, counter %f, ucount %i, suma %f nablak: %1.15f, deltak: %1.15f\n",iteration,counter,ucounter,suma,nablak,deltak);}
  delete IX; delete IX1; delete IX2; delete copy;
  return f;
  }
}
