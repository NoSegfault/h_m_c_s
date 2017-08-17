function [vmin_best, acbest, Bestcut]=hypergraphcut_V2(INC,w,vertex_weights,maxrun,init,p,verbosity)

% Balanced Hypergraph Cuts

% INC: We assume INC exists which is the hypergraph incidence matrix, every
% column represents a hypergraph edge (1 if the correspondent vertex given
% by the row number is in the edge and zero otherwise)

% w: the vector of edge weights

% maxrun: total number of initializations,

% if init is 1 then the first initialization is with the second eigenvector
% of the standard graph laplacian (zhou, sch??lkopf) otherwise all
% initializations are random

% createClustersUniversalHypergraph: computes the optimal clustering based
% on input vector vmin. Last input variable is unimportant (for hard
% balanced cut), second to last is always 1 since we look at l1 norm. Third
% to last is 1 for RCC/NCC and 2 for RCut/NCut. 4th to last is the vertex degree. 
%
% (C)2012-13 Matthias Hein, Simon Setzer, Leonardo Jost and Syama Sundar Rangapuram

[nVertices, nEdges] = size(INC);
%deg_normalized = sum(INC*spdiags(w,0,nEdges,nEdges),2);
Bestcut=inf;
Cuts=zeros(maxrun,1);
for run=1:maxrun

    randinit=randn(nVertices,1);
    if (init==1 && run==1)
        randinit = power_norm_lapl_hypergraph(INC,w,10e-12);
    end
    
    if p==2
        [vmin, RCC_NCC_Clusters, lambda, RCC_NCC]=ComputeEigCleanHypergraph_NCut(w,INC',vertex_weights,randinit,1000,1e-12, false, true,verbosity);
    end
    
    Bestcut_temp=min(RCC_NCC);
    Cuts(run)=Bestcut_temp;
    if Bestcut_temp<Bestcut
        Bestcut=Bestcut_temp;
        acbest=RCC_NCC_Clusters;
        vmin_best = vmin;
    end
    %fprintf('perturbation=> Hyp-IPM NCC: %f \n', RCC_NCC);
end

