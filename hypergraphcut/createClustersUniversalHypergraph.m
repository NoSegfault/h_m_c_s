function [allClusters, BestBalCut,cutBestBalCut,threshold_index] =  createClustersUniversalHypergraph(vmin,W,INC,deg,criterion,pnorm,K)
% Transforms an eigenvector into a cluster indicator function by thresholding.
% 
% Usage:	[allClusters, BestBalCut,cutBestBalCut,threshold_index] 
%			= createClustersUniversalHypergraph(vmin,W,INC,deg,criterion,pnorm,K)
%
% Input:
%   vmin: The eigenvector.
%   W: Weight matrix (edge weights)
%   INC: Incidence matrix of the hypergraph
%   deg: vertex weights
%   criterion: 1: Cheeger family, 2: Ratio Family, 3: Truncated Cheeger
%   pnorm: parameter for the individual criterion - norm for 1,2 and
%          truncation parameter for 3
%
% Output:
%   allClusters: Obtained clustering after thresholding (indicator vector
%   of the optimal set)
%   BestBalCut: best balanced cut criterion 
%   cutBestBalCut: cut value for the best balanced cut (directed - weight of outgoing edges)
%   threshold_index: The index of the sorted vmin used for thresholding
%
% (C)2012-13 Matthias Hein, Simon Setzer, Leonardo Jost and Syama Sundar Rangapuram
% Machine Learning Group, Saarland University
% http://www.ml.uni-saarland.de

    %Make deg a row vector;
	if (size(deg,1)>1) 
		deg=deg';
    end
 
    [vmin_sorted, index]=sort(vmin,'descend');
    INC_sorted=INC(:,index);
    
    tempcuts_threshold = mexCompHypergraphCuts(W,sparse(INC_sorted'))';
    
   
    % divide by volume/size
    volumes_threshold=cumsum(deg(index));
        
    switch criterion
        case 1,
            if(pnorm==1)
              balcut = tempcuts_threshold(1:end-1)./min(volumes_threshold(1:end-1),volumes_threshold(end)-volumes_threshold(1:end-1));
            else
              balfct = (volumes_threshold(1:end-1).*(volumes_threshold(end)-volumes_threshold(1:end-1))).^(1/pnorm)./(volumes_threshold(1:end-1).^(1/(pnorm-1)) + (volumes_threshold(end)-volumes_threshold(1:end-1)).^(1/(pnorm-1))).^(1-1/pnorm);
              balcut =  tempcuts_threshold(1:end-1)./balfct;
            end
        case 2,
             balfct = (volumes_threshold(1:end-1).*(volumes_threshold(end)-volumes_threshold(1:end-1)).^pnorm + volumes_threshold(1:end-1).^pnorm.*(volumes_threshold(end)-volumes_threshold(1:end-1))).^(1/pnorm)/volumes_threshold(end);
             balcut =  tempcuts_threshold(1:end-1)./balfct;
    end

    [BestBalCut,threshold_index]=min(balcut);
    cutBestBalCut =  tempcuts_threshold(threshold_index);
    allClusters=ones(size(vmin,1),1);
    allClusters(index(threshold_index+1:end))=0;
    threshold=vmin_sorted(threshold_index);
            
