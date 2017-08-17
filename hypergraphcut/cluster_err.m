% Computes the error of the clustering. The label of each cluster is
% obtained via a majority vote.
% the first argument is predicted labels and the second one is real labels!
%
% (C)2012-13 Matthias Hein, Simon Setzer, Leonardo Jost and Syama Sundar Rangapuram
function [error, u_reallabels,votes] = cluster_err(u,Y)

    labels=unique(u);
    error=0;
    
    u_reallabels=zeros(length(u));
    reallabels=unique(Y);
    votes=zeros(length(reallabels),length(labels));
    
    
    for k=1:length(labels)
        % extract indices of current cluster
        ukIndex=u==labels(k);

        % extract real labels of current cluster
        Yk=Y(ukIndex);

        % perform majority vote
        currentvotes=zeros(length(reallabels), 1);
        for l=1:length(reallabels)
            currentvotes(l) = sum(Yk==reallabels(l));
        end
        ind=find(currentvotes==max(currentvotes));

        votes(:,k)=currentvotes;
        
        % relabel u
        u_reallabels(ukIndex)=reallabels(ind(1));
        
        % compute error
        currentError= sum(Yk~=reallabels(ind(1)));
        error=error + currentError;
    end
    error=error/length(u);
end