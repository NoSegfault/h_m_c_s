function [clusters,cuts,cheegers] = computeMultiPartitioningHyperGraph(Method, INC,w,vertex_weights,k,cutType,l2init,numTrials,verbosity)
% Computes a multipartitioning of the hypergraph given by the Incidence 
% matrix, INC and weight vector w by recursive biparitioning.
%
% Usage:	
% [clusters,cuts,cheegers] = computeMultiPartitioning(Method,INC,w,vertex_weights,k,cutType,l2init,numTrials,verbosity)
%
% Input:
%   Method: 1 (ours), 2 (Zhou's).
%   INC (numEdges x numVertices) - sparse matrix with entries 0 (edge does not contain the vertex) or 1 (edge contains the vertex), 
%   w: vector of edge weights (numEdges x 1)
%   vertex_weights: deg (=INC'*w) for NCut, ones for RCut
%   k: number of clusters
%   cutType: 1: RCC/NCC, 2: Rcut/Ncut (at the moment only works for Ncut!)
%   l2init: true if you want to perform one initialization with the thresholded second eigenvector of the standard graph Laplacian.
%   numTrials: number of additional runs with different random initializations at each level
%   verbosity: Controls how much information is displayed.
%
% Output:
%   clusters: numVertices x (k-1) matrix containing in each column the computed 
%   clustering for each partitioning step.
%   cuts: (k-1) x 1 vector containing the Ratio/Normalized Cut values after 
%   each partitioning step.
%   cheegers: (k-1) x 1 vector containing the Ratio/Normalized Cheeger Cut 
%   values after each partitioning step.
%
% (C)2012-13 Matthias Hein, Simon Setzer, Leonardo Jost and Syama Sundar Rangapuram
% Machine Learning Group, Saarland University, Germany
% http://www.ml.uni-saarland.de
 
    num=size(INC,2);

    assert(k>=2,'Wrong usage. Number of clusters has to be at least 2.');
    assert(k<=num, 'Wrong usage. Number of clusters is larger than size of the graph.');   
    
    clusters=zeros(num,k-1);
    cuts=zeros(1,k-1);
    cheegers=zeros(1,k-1);
    cutParts=zeros(1,k);
    
    % Check if graph is connected    
    %Wspan = spanningGraph(INC);
    Wclique = INC'*INC;
    [comp,connected,sizes]=connectedComponents(Wclique);
    if(~connected)
        if(verbosity>=1) disp('WARNING! GRAPH IS NOT CONNECTED!');end
        if(verbosity>=2) disp('Optimal Cut achieved by separating connected components.');end
        allClusters = balanceConnectedComponents(comp,sizes,W,normalized);
        [cut,cheeger,cutPart1,cutPart2]=deal(0);
    else
        if(verbosity>=2) disp('Computing partitioning.'); end
        switch Method
            case 1,
                [vmin, allClusters, bestBalCut] = hypergraphcut_V2(INC',w,vertex_weights,numTrials,l2init,cutType,verbosity);
                if cutType==2, bestBalCut = bestBalCut*2; end
            case 2,
                [vmin, clusters_optGraph, balcut_optGraph, clusters0, balcut0, clusters_optHyper, balcut_optHyper] = hyp_graph_clustering_Zhou(INC, w, vertex_weights, false);                
                allClusters = clusters_optGraph; bestBalCut = balcut_optGraph;
        end
        unBalCut = mexEvalHypergraphObjVarTransform(w,sparse(INC'),allClusters);
        cutPart1 = unBalCut/sum(vertex_weights(allClusters==0));   %norm(deg.*fold,1);
        cutPart2 = unBalCut/(sum(vertex_weights)-sum(vertex_weights(allClusters==0)));  %norm(deg.*fold,1);
        cut = cutPart1 + cutPart2; 
        cheeger = max(cutPart1, cutPart2);
        if cutType==1
            assert(bestBalCut == cheeger);
        else
            %assert(abs(bestBalCut-cut) < 1e-8);
        end
    end
           
    allClusters=allClusters+1;
    clusters(:,1)=allClusters;
   
    cuts(:,1)=cut;
    cheegers(:,1)=cheeger;
    cutParts(1)=cutPart1;
    cutParts(2)=cutPart2;

    subCutParts=zeros(k,2);
    subClusters=cell(1,k);

    if(verbosity>=1)
        fprintf('Finished Clustering into 2 parts.\n');
        displayCurrentObjective(cut,cheeger);
        fprintf('\n');
    end

    %Perform the 2nd to (k-1)th partitioning step
    for l=3:k
        bestCut=inf;
        bestCheeger=inf;
        % in each step consider each of the current l-1 clusters
        for m=1:l-1

            index_m=find(allClusters==m);
            
            % if we have already solved this subproblem	
            if (~isempty(subClusters{m}))
				allClustersInCluster_m = subClusters{m};
                cutPart1 = subCutParts(m,1);
                cutPart2 = subCutParts(m,2);
            % if the current cluster has size 1 it cannot be further divided
            elseif(length(index_m)==1)
                allClustersInCluster_m=[];
                cutPart1=inf;
                cutPart2=inf;
                subClusters{m}=allClustersInCluster_m;
                subCutParts(m,1:2)=[cutPart1 cutPart2];
   
            % otherwise we have to compute the partition
            else    
                if(verbosity>=2) fprintf('Computing partitioning of subgraph %d.\n',m);end
                
                % extract subgraph and its connected components
                vwm = vertex_weights(index_m);
                INCm = INC(:, index_m); % this may result in empty hyper-edges.
                dummy_ix = find(sum(INCm,2)>0);
                INCm = INCm(dummy_ix, :);
                wm = w(dummy_ix);
                
                %WmSpan=spanningGraph(INCm);               
                WmClique= INCm'*INCm;
                [comp,connected,sizes]=connectedComponents(WmClique);
                
                if(verbosity>=2 && ~connected) disp('...Subgraph is not connected.'); end
                cutPart1=inf;
                cutPart2=inf;
                
                % if subgraph is not connected
                if (~connected) 
                                                  
                    %check partition which has connected component as one 
                    %cluster and rest as second cluster
                    for m1=1:length(sizes)   
                        %if(true)
                        if(cutPart1+cutPart2>0)
                            if(verbosity>=2) fprintf('...Checking partition found by isolating connected component %d of %d.\n',m1,length(sizes)); end

                            allClustersInCluster_m_temp = double(comp==m1);

                            cluster_m2=zeros(size(allClusters,1),1);
                            cluster_m2(index_m)=allClustersInCluster_m_temp;
                            cutPart2_temp = mexEvalHypergraphObjVarTransform(w,sparse(INC'),cluster_m2)/sum(vertex_weights(cluster_m2>0)); % set C corresponds to 1 in cluster_m2
                            

                            cluster_m1=zeros(size(allClusters,1),1);
                            cluster_m1(index_m)=(allClustersInCluster_m_temp==0);
                            cutPart1_temp = mexEvalHypergraphObjVarTransform(w,sparse(INC'),cluster_m1)/sum(vertex_weights(cluster_m1>0)); % set C corresponds to 1 in cluster_m1
                            
                            % Display current objective
                            if(verbosity>=2) 
                                [cut_temp,cheeger_temp]=computeCutCheeger(cutParts,cutPart1_temp,cutPart2_temp,m,l);
                                displayCurrentObjective(cut_temp,cheeger_temp); 
                            end

                            %Check if we're better						
                            if (cutType==2 && cutPart1_temp+cutPart2_temp<cutPart1+cutPart2 || cutType==1 && max(cutPart1_temp,cutPart2_temp)<max(cutPart1,cutPart2))
                                [cutPart1,cutPart2,allClustersInCluster_m]=deal(cutPart1_temp,cutPart2_temp,allClustersInCluster_m_temp);
                                assert(length(allClustersInCluster_m)==length(index_m));
                            end
                            assert(logical(exist('allClustersInCluster_m','var')));
                        end 
                    end
                end
                %if(true)
                if(cutPart1+cutPart2>0)
                    for m1=1:length(sizes)
                        index_comp=find(comp==m1);
                        % if the size of the current connected component is larger than 1, try to partition it
                        if (length(index_comp)>1)
                            vwm_comp = vwm(index_comp);
                            INCm_comp = INCm(:, index_comp); % this may result in empty hyper-edges.
                            dummy_ix = find(sum(INCm_comp,2)>0);
                            INCm_comp = INCm_comp(dummy_ix, :);
                            wm_comp = wm(dummy_ix);
                            
                            if(~connected && verbosity>=2) fprintf('...Computing partitioning of connected component %d of %d of subgraph %d.\n',m1,length(sizes), m); end
                            switch Method
                                case 1,
                                        [vmin_m_comp, allClusters_dummy, bestBalCut] = hypergraphcut_V2(INCm_comp',wm_comp,vwm_comp,numTrials,l2init,cutType,verbosity);                             
                                case 2,
                                        [vmin_m_comp, allClusters_dummy, betsBalCut] = hyp_graph_clustering_Zhou(INCm_comp, wm_comp, vwm_comp, false);
                            end
                            vmin_m1 = max(vmin_m_comp)*ones(length(index_m), 1) + 1e-6;
                            vmin_m1(index_comp) = vmin_m_comp; 
                            [allClustersInCluster_m_temp, cutPart1_temp, cutPart2_temp] = createSubgraphClustersHypergraph(index_m,vmin_m1,w,INC,vertex_weights(index_m),cutType);
                            %createClustersHypergraphOffset(vmin_m1,w,INC,vertex_weights(index_m),criterion);
                            
                            vmin_m2 = min(vmin_m_comp)*ones(length(index_m), 1) - 1e-6;
                            vmin_m2(index_comp) = vmin_m_comp;
                            [allClustersInCluster_m_temp_dummy, cutPart1_temp_dummy, cutPart2_temp_dummy] = createSubgraphClustersHypergraph(index_m,vmin_m2,w,INC,vertex_weights(index_m),cutType);
                            
                            if cutPart1_temp_dummy + cutPart2_temp_dummy < cutPart1_temp + cutPart2_temp 
                                [allClustersInCluster_m_temp, cutPart1_temp, cutPart2_temp] = deal(allClustersInCluster_m_temp_dummy, cutPart1_temp_dummy, cutPart2_temp_dummy);
                            end
                            

                            % Display current objective
                            if(verbosity>=2) 
                                [cut_temp,cheeger_temp]=computeCutCheeger(cutParts,cutPart1_temp,cutPart2_temp,m,l);
                                displayCurrentObjective(cut_temp,cheeger_temp); 
                            end
                                    
                            if (cutType==2 && cutPart1_temp+cutPart2_temp<cutPart1+cutPart2 || cutType==1 && max(cutPart1_temp,cutPart2_temp)<max(cutPart1,cutPart2))
                                    [cutPart1,cutPart2,allClustersInCluster_m]=deal(cutPart1_temp,cutPart2_temp,allClustersInCluster_m_temp);
                                    assert(length(allClustersInCluster_m)==length(index_m));
                            end
                            
                        end

                    end
                end
                % store current best partition
			    subClusters{m}=allClustersInCluster_m;
				subCutParts(m,1:2)=[cutPart1 cutPart2]; 
   
            end
            
            % print out best cut possible by partitioning of current subgraph
            [cut,cheeger]=computeCutCheeger(cutParts,cutPart1,cutPart2,m,l);
            if(verbosity>=2)
                fprintf('Best result achievable by partitioning of subgraph %d:\n',m);
                displayCurrentObjective(cut,cheeger);
                fprintf('\n');
            end
            
			% check if partitoning of the current subgraph gives better cut
            if (cutType==2 && cut<bestCut) || (cutType==1 && cheeger<bestCheeger)
                [bestCut,bestCheeger,bestCutPart1,bestCutPart2,best_m]= deal(cut,cheeger,cutPart1,cutPart2,m);
                clusters_new=allClusters;
                clusters_new(index_m)=(l-m)*allClustersInCluster_m+clusters_new(index_m);

                % assert(bestCut>=0 && bestCheeger>=0);
            end
            
            % if we have already found a partition with cut 0, we don't
            % need to consider the other subgraphs
            if bestCut==0
                break;
            end
        end
        
        if(bestCut==inf)
             error('OneSpect:cutinf','Clustering initialized with second eigenvector of standard graph Laplacian failed at level %d.',l-1);
        end
        
        % Update
        allClusters=clusters_new;
        cuts(1,l-1)=bestCut;
        cheegers(1,l-1)=bestCheeger;
        clusters(:,l-1)=allClusters;
        
        cutParts(best_m)=bestCutPart1;
        cutParts(l)=bestCutPart2;
        
        % Check that we have the right number of clusters
        assert(length(unique(allClusters))==l);
        
        % Reset subcutparts and subclusters;
        subCutParts(best_m,:)=0;
        subClusters{best_m}= [];
        subCutParts(l,:)=0;
        subClusters{l}= [];
          
        % Print out current objective
        if(verbosity>=1)
            fprintf('Decided to partition subgraph %d. Finished Clustering into %d parts.\n',best_m,l);
            displayCurrentObjective(bestCut,bestCheeger);
            fprintf('\n');
        end

    end
    if(cutType==2)
     disp(['Finished recursive splitting - clustering has normalized/ratio hypergraph cut: ',num2str(cuts(end),'%1.5f')]);
    end
    if(cutType==1)
      disp(['Finished recursive splitting - clustering has normalized/ratio cheeger cut: ',num2str(cuts(end),'%1.5f')]);
    end  
    
end

% Computes Rcut/Ncut and Cheeger Cut values
function [cut,cheeger]=computeCutCheeger(cutParts,cutPart1,cutPart2,m,l)

    cut= sum(cutParts)-cutParts(m)+cutPart1+cutPart2;
    cheeger=max([cutParts((1:l-1)~=m) cutPart1 cutPart2]);
                                
end

% Displays the current objective value
function displayCurrentObjective(cut_temp,cheeger_temp)
    
      fprintf('...Balanced Cut: %.8g   Balanced Cheeger Cut: %.8g\n',cut_temp,cheeger_temp); 

end


% Creates two clusters by thresholding the vector vmin_comp obtained on a
% connected component of a subgraph. Given the two clusters on the 
% connected component, there are two ways of constructing the final clusters 
% on the subgraph, as we can keep each of the clusters on the connected 
% component as cluster and merge the other one with the remaining connected 
% components. The method takes the one yielding the lower Cut/Cheeger.
function [allClustersInClusterM, cutPart1,cutPart2] =  createSubClusters2(vmin_comp,W_comp,normalized,deg,criterion_threshold,index_comp,index_m,cut_rest,size_rest,size_m)
      
        % input parameter deg has to be the degree vector (also in unnormalised case)
        %deg=full(sum(W));
        %Make deg a row vector;
        if (size(deg,1)>1) 
            deg=deg';
        end
        
        [vminM_sorted, index]=sort(vmin_comp);
        [vminU,indexU]=unique(vminM_sorted);
        
        
        W_sorted=W_comp(index,index);

        % calculate cuts
        deg_comp=deg(index_m(index_comp));
        volumes_threshold=cumsum(deg_comp(index));
        triup=triu(W_sorted);
        tempcuts_threshold=volumes_threshold - 2*cumsum(full(sum(triup)));
        tempcuts_threshold2=(volumes_threshold(end)-volumes_threshold) - (sum(sum(W_sorted))-2*cumsum(full(sum(triup,2)))');            

        % it may happen that (due to numerical imprecision) the tempcuts
        % are a small factor of epsilon below zero.
        tempcuts_threshold(tempcuts_threshold<0)=0;
        tempcuts_threshold2(tempcuts_threshold2<0)=0;
        
        tempcuts_threshold=tempcuts_threshold(indexU);
        tempcuts_threshold2=tempcuts_threshold2(indexU);
        volumes_threshold=volumes_threshold(indexU);
        
        
        % divide by size/volume
        if(normalized)
            cutparts1_threshold=(tempcuts_threshold(1:end-1)+cut_rest)./(volumes_threshold(1:end-1)+size_rest);
            cutparts1_threshold(isnan(cutparts1_threshold))=0;
            cutparts2_threshold=tempcuts_threshold2(1:end-1)./(volumes_threshold(end)-volumes_threshold(1:end-1));
            cutparts2_threshold(isnan(cutparts2_threshold))=0;
            
            cutparts1b_threshold=tempcuts_threshold(1:end-1)./volumes_threshold(1:end-1);
            cutparts1b_threshold(isnan(cutparts1b_threshold))=0;
            cutparts2b_threshold=(tempcuts_threshold2(1:end-1)+cut_rest)./((volumes_threshold(end)-volumes_threshold(1:end-1))+size_rest);
            cutparts2b_threshold(isnan(cutparts2b_threshold))=0;
        else
            sizes_threshold=cumsum(ones(1,size(vmin_comp,1)-1));
            sizes_threshold=sizes_threshold(indexU(1:end-1));
            cutparts1_threshold=(tempcuts_threshold(1:end-1)+cut_rest)./(sizes_threshold+size_rest);
            cutparts2_threshold=tempcuts_threshold2(1:end-1)./(size(vmin_comp,1)-sizes_threshold);
            
            cutparts1b_threshold=tempcuts_threshold(1:end-1)./sizes_threshold;
            cutparts2b_threshold=(tempcuts_threshold2(1:end-1)+cut_rest)./((size(vmin_comp,1)-sizes_threshold)+size_rest);
        end

        
        
        %calculate total cuts
        if(criterion_threshold==1)
            cuts_threshold=cutparts1_threshold+cutparts2_threshold;
            [cut1,threshold_index]=min(cuts_threshold);
            
            cuts_threshold_b=cutparts1b_threshold+cutparts2b_threshold;
            [cut1b,threshold_index_b]=min(cuts_threshold_b);
            
            comp_case=1;
            if (cut1b<cut1) 
                comp_case=2;
            end
        else
            cheegers_threshold=max(cutparts1_threshold,cutparts2_threshold);
            [cheeger1,threshold_index]=min(cheegers_threshold);
            
            cheegers_threshold_b=max(cutparts1_threshold,cutparts2_threshold);
            [cheeger1b,threshold_index_b]=min(cheegers_threshold_b);
            
            comp_case=1;
            if (cheeger1b<cheeger1) 
                comp_case=2;
            end
        end

        if(comp_case==1)
            cutPart1=cutparts1_threshold(threshold_index);
            cutPart2=cutparts2_threshold(threshold_index);

            allClustersInClusterM_comp= (vmin_comp>vminU(threshold_index));
        
            allClustersInClusterM= zeros(size_m,1);
            allClustersInClusterM(index_comp)=allClustersInClusterM_comp;
        else
            cutPart1=cutparts1b_threshold(threshold_index_b);
            cutPart2=cutparts2b_threshold(threshold_index_b);

            allClustersInClusterM_comp= (vmin_comp>vminU(threshold_index_b));
        
            allClustersInClusterM= ones(size_m,1);
            allClustersInClusterM(index_comp)=allClustersInClusterM_comp;
        end
        
        

end


% Tries to separate the connected components into two clusters which have
% roughly the same cardinality/volume
function comp2 = balanceConnectedComponents(comp,sizes,W,normalized)

    % for normalized variant, compute the volume for every connected
    % component
    if(normalized)
        deg=sum(W);
        volumes=zeros(length(sizes),1);
        for l=1:length(sizes)
            volumes(l)=sum(deg(comp==l));
        end
        sizes=volumes;
    end
            
    % fill up clusters, trying to balance the size
    [sizes_sort,ind]=sort(sizes,'descend');
    size_a=0;
    size_b=0;
    ind_a=[];
    for l=1:length(sizes_sort)
        if(size_a<=size_b) 
            size_a=size_a+sizes_sort(l);
            ind_a=[ind_a ind(l)];
        else size_b=size_b+sizes_sort(l);
        end
    end
    comp2=double(ismember(comp,ind_a));
end
      
        
