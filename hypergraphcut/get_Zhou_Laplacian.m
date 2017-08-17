function [Lz, Theta] = get_Zhou_Laplacian(INC, W)
%   Computes the Zhou Laplacian defined as Lz = I - Dv^{-1/2} HW De^{-1} H'Dv^{-1/2} 
%   INC: Adjacency matrix: |E|x|V|
%   W: weight vector |E| x 1
%   Lz: Zhou Laplacian
%   Theta: 
	    
    [numEdges, num] = size(INC);
    deltaE=sum(INC,2);          % number of vertices per edge
    dv = sum(spdiags(W,0,numEdges,numEdges)*INC,1); dv=dv';
 
    ReW = spdiags(W./deltaE, 0, numEdges,numEdges); % diagonal weight matrix              
    
    Dv_inv_half=spdiags(1./sqrt(dv),0,num,num);
    M = INC*Dv_inv_half;
    
    Theta = M'*ReW*M;
    clear M; clear ReW;
    Lz = speye(num) -  Theta;
end


% Sanity check on regular graphs:
% [ix, jx, wval] = find(triu(W)); % for normalized Laplacian it does not matter whether take triu(W) or W.
% H = sparse(1:length(ix), ix, 1, length(ix), size(W,1)) + sparse(1:length(ix), jx, 1, length(ix), size(W,1));
% Lz = get_Zhou_Laplacian(H, wval);
% n = size(W,1); d = sum(W,2); Dsqrt =spdiags(1./(d.^0.5), 0, n, n); L = (speye(n)-Dsqrt*W*Dsqrt);
% assert( abs( sum(sum(abs(Lz-0.5*L))) ) <= 1e-8)