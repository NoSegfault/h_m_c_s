function [allClusters, bestCutPart1, bestCutPart2] =  createSubgraphClustersHypergraph(index_m,vmin_m,W,INC,deg_m,criterion)
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
% Written by Matthias Hein
% Machine Learning Group, Saarland University
% http://www.ml.uni-saarland.de

    %Make deg a row vector;
	if (size(deg_m,1)>1) 
		deg_m=deg_m';
    end
 
    topK = length(vmin_m);
    num = size(INC,2);
    vmin1 = max(vmin_m)*ones(num, 1) + 1e-4;
    vmin1(index_m) = vmin_m;
    
    vmin2 = min(vmin_m)*ones(num, 1) - 1e-4;
    vmin2(index_m) = vmin_m;
                            
    %[vmin_sorted, index]=sort(vmin,'descend');
    [vmin_sorted, index1]=sort(vmin1);
    INC_sorted1=INC(:,index1);
    tempcuts_threshold1 = mexCompHypergraphCutsPartial(W,sparse(INC_sorted1'),topK)';
    
    [vmin_sorted, index2]=sort(vmin2);
    INC_sorted2=INC(:,index2);    
    tempcuts_threshold2 = mexCompHypergraphCutsPartialReverse(W,sparse(INC_sorted2'),topK)';
   
    index = index1(1:topK);
    assert(sum(index ~= index2(end-topK+1:end))==0);
    
    [vmin_m_sorted, index]=sort(vmin_m); % Sorting two times is less efficient; need to change this by modifying above index suitably.
                                         % However this does not cause problem since: If vmin_m has repeated
                                         % elements of equal value, the returned indices preserve the original ordering.
    % divide by volume/size
    volumes_threshold=cumsum(deg_m(index));
        
    switch criterion
        case 1,
%             if(pnorm==1)
%               balcut = tempcuts_threshold(1:end-1)./min(volumes_threshold(1:end-1),volumes_threshold(end)-volumes_threshold(1:end-1));
%             else
%               balfct = (volumes_threshold(1:end-1).*(volumes_threshold(end)-volumes_threshold(1:end-1))).^(1/pnorm)./(volumes_threshold(1:end-1).^(1/(pnorm-1)) + (volumes_threshold(end)-volumes_threshold(1:end-1)).^(1/(pnorm-1))).^(1-1/pnorm);
%               balcut =  tempcuts_threshold(1:end-1)./balfct;
%             end

             balcut1 = tempcuts_threshold1(1:end-1)./volumes_threshold(1:end-1);
             balcut2 = tempcuts_threshold2(1:end-1)./(volumes_threshold(end)-volumes_threshold(1:end-1));  
             balcut = max(balcut1, balcut2);
             
        case 2,
             %balfct = (volumes_threshold(1:end-1).*(volumes_threshold(end)-volumes_threshold(1:end-1)).^pnorm + volumes_threshold(1:end-1).^pnorm.*(volumes_threshold(end)-volumes_threshold(1:end-1))).^(1/pnorm)/volumes_threshold(end);
             balcut1 = tempcuts_threshold1(1:end-1)./volumes_threshold(1:end-1);
             balcut2 = tempcuts_threshold2(1:end-1)./(volumes_threshold(end)-volumes_threshold(1:end-1));  
             balcut = balcut1+balcut2;
             
        case 3, % truncated Cheeger cut - pnorm is the parameter for the truncated cut: number of points from whereon criterion is constant
             balcut = tempcuts_threshold(1:end-1)./min(min(volumes_threshold(1:end-1),volumes_threshold(end)-volumes_threshold(1:end-1)),pnorm);
        case 4, % Hard Cheeger cut
             balcut = tempcuts_threshold(1:end-1)./max( min(volumes_threshold(1:end-1),volumes_threshold(end)-volumes_threshold(1:end-1))-(pnorm-1),0);
        case 6, % truncated hard cheeger cu          
             balcut = tempcuts_threshold(1:end-1)./max(min(min(volumes_threshold(1:end-1),volumes_threshold(end)-volumes_threshold(1:end-1)),pnorm) -(K-1),0);  
    end

    [BestBalCut,threshold_index]=min(balcut);
    allClusters=zeros(size(vmin_m,1),1);
    allClusters(index(threshold_index+1:end))=1;
            
    bestCutPart1 = balcut1(threshold_index);
    bestCutPart2 = balcut2(threshold_index);
    
    if criterion==2, assert(BestBalCut == bestCutPart1+bestCutPart2); end
    if criterion==1, assert(BestBalCut == max(bestCutPart1,bestCutPart2)); end
