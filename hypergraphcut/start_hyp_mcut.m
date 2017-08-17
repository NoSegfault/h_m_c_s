% this file runs the hypergraph clustering method based on minimization
% of the normalized hypergraph cut 
%
% M. Hein, S. Setzer, L. Jost, and S. Rangapuram,
% The Total Variation on Hypergraphs - Learning on Hypergraphs Revisited
% NIPS 2013
%
% All copyrights remain with the authors - you have to cite the above
% paper when using this code.
% with the example from the zoo dataset
%
% it runs the file with 2 inner random initializations plus initialization
% using the second eigenvector of the induced graph laplacian found by the
% clique expansion
%
% the whole procedure is repeated 5 times and the best found clustering is
% reported
%
% (C)2012-13 Matthias Hein, Simon Setzer, Leonardo Jost and Syama Sundar Rangapuram

load Flori.mat   % Need three variables: INC (incidence matrix), w (weight vector), and Y (true labels) to compute error

if ~exist('Y', 'var'), Y = zeros(size(INC,2),1); end  
% if labels are available one can compute the clustering error induced
% by majority vote on the clusters - this corresponds to the expected error
% of a semi-supervised learning procedure which classifies each cluster by
% majority vote

Method = 1; % 1: normalized hypergraph cut as described in the paper
            % 2: the method by Zhou et al based on hypergraph approximation
cutType = 2; % 1: Cheeger cut, 2: normalized cut (recommended)
l2init = true; % initialization with second eigenvector (recommended)
nInnerRuns = 1; % the more the higher quality one can expect
verbosity = 1;  % controls the output level (possible values 0,1,2) - the higher the more output
vertex_weights = INC'*w; % this computes the degree of each vertex
nOuterRuns = 1;  % the more the higher quality one can expect
bestCut = inf;
all_cuts = zeros(nOuterRuns,1);
all_errors = zeros(nOuterRuns,1);
k = length(unique(Y));

% for test and understanding
k = 4;

%if k==2, nOuterRuns = 1; nInnerRuns = 10; end
for i=1:nOuterRuns
    [clusters,cuts,cheegers] = computeMultiPartitioningHyperGraph(Method,INC,w,vertex_weights,k,cutType,l2init,nInnerRuns,verbosity);
    all_cuts(i) = cuts(end); 
    %all_errors(i) = cluster_err(clusters(:,end), Y);
    if cuts(end) < bestCut
        bestCut = cuts(end);
        best_i = i;
        bestClusters = clusters(:,end);
        bestAllClusters = clusters;
    end
end
disp(['Best clustering found has normalized hypergraph cut: ',num2str(bestCut,'%1.5f')]);
% bestClusters contains the best found clustering 
% bestAllClusters shows the result of all recursive splits
% bestCut shows the final normalized hypergraph cut
our_balcut = bestCut; assert(bestCut == all_cuts(best_i));
our_clustering = bestClusters; %our_error = all_errors(best_i);
% assert(our_error == cluster_err(our_clustering, Y));

% visualize the recursive splitting!
if length(unique(Y))>1, [a,ix] = sort(Y); imagesc(bestAllClusters(ix,:)); end
disp(bestAllClusters(ix,:));