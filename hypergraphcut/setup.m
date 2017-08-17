% the following need to be compiled before running it:
mex -largeArrayDims mexEvalHypergraphObjVarTransform.cpp
mex -largeArrayDims mexCompHypergraphCuts.cpp
mex -largeArrayDims mexInnerFISTAHypergraphVarTransform3.cpp
mex -largeArrayDims mexCompHypergraphCutsPartial.cpp
mex -largeArrayDims mexCompHypergraphCutsPartialReverse.cpp
disp(['... finished compilation - see ReadMe or start_hyp_mcut how to proceed.']);