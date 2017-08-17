function balcut = balCut_hypGraph(INC, w, vertex_weights, clusters)

    balcut = 0;
    classes = unique(clusters);
    for i=1:length(classes)
        unBalCut = mexEvalHypergraphObjVarTransform(w,sparse(INC'),double(clusters==classes(i)));
        balcut = balcut +  unBalCut/ sum(vertex_weights(clusters==classes(i)));
    end
    %disp(balcut);
end