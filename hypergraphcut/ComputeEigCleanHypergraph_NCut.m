function [xcur, ac, lambda, balcut, Data]=ComputeEigCleanHypergraph_NCut(W,INC,deg,fold,MAXITER,EPS, thrs_per_iter, Lovasz,verbosity)
% (C)2012-13 Matthias Hein, Simon Setzer, Leonardo Jost and Syama Sundar Rangapuram
  if (~exist('Lovasz','var')), Lovasz = 1; end
  if (~exist('thrs_per_iter','var')), thrs_per_iter = 0; end
  
  L = sqrt(4*max(sum(INC,1).^2)); % Lipschitz constant of the gradient of the objective
  lambdaOld=inf;
  lambda = mexEvalHypergraphObjVarTransform(W,sparse(INC'),fold)/balance_ncut(deg, fold);
  Deg = spdiags(deg,0,size(fold,1),size(fold,1));
  n = length(fold); volV = sum(deg);

    if ~Lovasz  
        disp('ignored');
    else    
        [fsort,sortind]=sort(fold);
        sdeg = deg(sortind);
        
        vec = zeros(n,1);   
        % sum(sdeg) = total vertices weights
        cumvols = sum(sdeg) - [0; cumsum(sdeg(1:n-1))];
        % cumvols is a vector of 17*101 to 17*1 with step 17
        % volV is total vetices weights
        
        % volV - cumvols - [cumvols(2:end); 0] is a vector of 17*-100 to
        % 17*100 with step 34, then each entry * 34 / volV, then sorted
        % corresponds to g
        
        vec(sortind) = 2*sdeg.*(volV - cumvols - [cumvols(2:end); 0])/volV;
    end  
    
  [ac, balcut,cut,threshold_index] =  createClustersUniversalHypergraph(fold,W,INC,deg,2,1,0);
  %disp(ac);
  if(verbosity>=1) disp(['Initial FctValue: ',num2str(lambda,'%1.14f'),' - Initial NCut: ',num2str(balcut,'%1.14f')]); end
  bestf = fold; bestac = ac; best_balcut = balcut; best_lambda = lambda;
  
  Data = [];
  
  [ix,jx]=find(INC);
  SINC = sparse(ix,jx,1:length(ix));
  [ix,jx,Ixx]=find(SINC');
  clear ix; clear jx; clear SINC;
  u = zeros(nnz(INC),1);
  v = zeros(nnz(INC),1);
  
  
  counter=0;
  tic;
  Iterations=[];
  while(abs(lambdaOld-lambda)/lambda>1E-13 && counter<50)
    %PDHG
    [g,u,v,Obj,iter]=mexInnerFISTAHypergraphVarTransform3(W,sparse(INC),sparse(INC'),Ixx-1,lambda*vec,u,v,MAXITER,EPS,L,verbosity);
    %changed because of Lovasz: 
    lambdaNew = mexEvalHypergraphObjVarTransform(W,sparse(INC'),g)/balance_ncut(deg, g);
    if Lovasz
        [gsort,gsortind]=sort(g);
        gsdeg = deg(gsortind);

        gvec = zeros(n,1);
        gcumvols = sum(gsdeg) - [0; cumsum(gsdeg(1:n-1))];
        gvec(gsortind) = 2*gsdeg.*(volV - gcumvols - [gcumvols(2:end); 0])/volV;
        
        % S2 = 0, so the under is S1(f^(k+1)) which is (gvec'*f^(k+1))
        
        % algorithm 1 line 5
        lambdaNew = mexEvalHypergraphObjVarTransform(W,sparse(INC'),g)/(gvec'*g);
    end
    if(verbosity>=1) disp(['New Functional: ',num2str(lambdaNew,'%1.14f')]); end
    if(lambdaNew < lambda)
        fold = g;
        vec = gvec;
        lambdaOld=lambda;
        lambda = lambdaNew;
        if thrs_per_iter
            disp('ignored');
        else
            if(verbosity>=1) disp(['Functional: ',num2str(lambda,'%1.14f')]); end
        end
    else
        if(verbosity>=1) disp(['Iteration stopped - function value not improved']); end
        break;
    end
  counter=counter+1;
  time = toc;
  Data = [Data [lambda;time]];
  end
  
disp('after iterations');
if(verbosity>=1) disp(['Iterations: ',num2str(Iterations)]);end
if thrs_per_iter, 
    disp('ignored');
else [ac, balcut,cut,threshold_index] =  createClustersUniversalHypergraph(fold,W,INC,deg,2,1,0);
    if(verbosity>=1) disp(['Final Result: Functional: ',num2str(lambda,'%1.14f'),' - NCut: ',num2str(balcut,'%1.14f')]); end
end
xcur=fold;

function y = LovaszHyperCut(g,W,INC,numEdges)
y=0;
for i=1:numEdges
 y = y + W(i)*(max(g(INC(i,:)))-min(g(INC(i,:))));
end
  
