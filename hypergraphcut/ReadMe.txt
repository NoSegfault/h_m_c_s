This is an implementation of the hypergraph clustering method
as described in 

M. Hein, S. Setzer, L. Jost, and S. Rangapuram,
The Total Variation on Hypergraphs - Learning on Hypergraphs Revisited
NIPS 2013

All copyrights remain with the authors - you have to cite the above
paper when using this code.

For problems/suggestions/success stories contact
Matthias Hein - hein@cs.uni-sb.de
or
Syama Sundar Rangapuram - srangapu@mpi-inf.mpg.de



First step:
run setup.m to compile the mex files

Usage of code 
You can either just replace the file (zoo.m) in start_hpy_mcut.m with your own file containing
incidence matrix and edge weights (and possible ground-truth labels if available) 

or directly the method for computing the hypergraph cut

[clusters,cuts,cheegers] = computeMultiPartitioningHyperGraph(Method,INC,w,vertex_weights,k,cutType,l2init,nInnerRuns,verbosity);

Method          % 1: normalized hypergraph cut as described in the paper
                % 2: the method by Zhou et al based on hypergraph approximation for comparison
INC             % incidence matrix (numEdges x numVertices)
w               % hypergraph edge weights (numEdges x 1)
vertex_weights  % hypergraph vertex weights (numVertices x 1) - for ratio cut d=ones(numVertices,1); for normalized cut d=INC'*w; (computes hypergraph degrees)
k               % number of clusters
cutType         % 1: Cheeger cut, 2: normalized cut (recommended)
l2init          % true= initialization with second eigenvector (recommended), false=no initialization with second eigenvector
nInnerRuns      % the more the higher quality one can expect (for k=2, we recommend 10 inner iterations and no outer loops, for k>2 we recommend
                % 2 inner iterations and 5 outer loops)
verbosity = 0;  % controls the output level (possible values 0,1,2)

If used as a stand-alone method it makes sense for multi-cut to repeat the method several times (5x recommended) and take the best result
as done in start_hyp_mcut.m
