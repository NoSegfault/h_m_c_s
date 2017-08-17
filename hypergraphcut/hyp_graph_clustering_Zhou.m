function [v2, clusters_optGraph, balcut_optGraph, clusters0, balcut0, clusters_optHyper, balcut_optHyper] = hyp_graph_clustering_Zhou(INC, w, vertex_weights, efficient)


    if ~exist('efficient', 'var')
        efficient = false;
    end
    
    if ~efficient
        disp('Computing the Laplacian');
        Lz = get_Zhou_Laplacian(INC, w);
        if ~isequal(Lz, Lz')
            Lz = 0.5*(Lz+Lz');
        end
        options.disp = 0;
        [v, l] = eigs(Lz, 2, 'SA', options);

        v2 = v(:,2); l2 = l(2,2);
        %clusters = v2>0;        
    else
        disp('efficient way');
        v2 = power_norm_lapl_hypergraph(INC',w,1e-6);
    end
    
    clusters0 = double(v2>=0);
    balcut0 = balCut_hypGraph(INC, w, vertex_weights, clusters0);
    
    disp(v2);
    [clusters_optHyper, balcut_optHyper] = createClustersUniversalHypergraph(v2,w,INC,vertex_weights,2,1,0);
    balcut_optHyper = balcut_optHyper*2;
    disp(clusters_optHyper);
        
    de = sum(INC,2);
    %W = INC'*diag(w./de)*INC;
    nEdges = length(w);
    W = INC'*spdiags(w./de,0,nEdges,nEdges)*INC;
    n = size(W,1);
    W(1:n+1:n*n) = 0;
    assert(sum(diag(W))==0);
    [clusters_optGraph, cutPart1,cutPart2] = opt_thresh_ncut(v2, W, vertex_weights, 1);
    balcut_optGraph = balCut_hypGraph(INC, w, vertex_weights, clusters_optGraph);
    disp(clusters_optGraph);
end