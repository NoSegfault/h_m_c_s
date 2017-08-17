function v2 = power_norm_lapl_hypergraph(H,w,eps)
% (C)2012-13 Matthias Hein, Simon Setzer, Leonardo Jost and Syama Sundar Rangapuram

% Zhou:
[n_vert, n_hedg] = size(H);
d = sum(H*spdiags(w,0,n_hedg,n_hedg),2);
d_inv_half = 1./sqrt(d);
card_edg = sum(H,1)';
w_div_card_edg = w./card_edg;

v2 = zeros(n_vert,1);
if mod(n_vert,2) == 0
    v2(1:n_vert/2) = -1;
    v2(n_vert/2+1:end) = 1;
else
    v2(1:(n_vert-1)/2) = -1;
    v2((n_vert+3)/2:end) = 1;
end
v2 = d_inv_half.*v2;

% this is thero: v2'*sqrt(d)

dd = sqrt(d)/norm(sqrt(d));

error = inf;
k = 0;

while error > eps
    v2_old = v2;
    v2 = v2/norm(v2);
    v2 = d_inv_half.*(H*(w_div_card_edg.*(H'*(d_inv_half.*v2))));
    v2 = v2 - (v2'*dd)*dd;
    error = norm(v2-v2_old);
    k = k+1;
end

%k

%T = diag(d_inv_half)*H*diag(w_div_card_edg)*H'*diag(d_inv_half)
%[U S] = eig(T)

%s = sort(diag(S),'descend')
%s(2)
%v2
%T*v2 - s(2)*v2

% [n_vert, n_hedg] = size(H);
% d = sum(H*spdiags(w,0,n_hedg,n_hedg),1);
% d_inv_half = 1./sqrt(d);
% card_edg = sum(H,2);
% w_div_card_edg = w./card_edg;
% 
% v2 = zeros(n_vert,1);
% if mod(n_vert,2) == 0
%     v2(1:n_vert/2) = -1;
%     v2(n_vert/2+1:end) = 1;
% else
%     v2(1:(n_vert-1)/2) = -1;
%     v2((n_vert+1)/2:end) = 1;
% end
% v2 = 
% 
% for i = 1:100
%     v2 = d_inv_half.*(H*(w_div_card_edg.*(H'*d_inv_half.*v2)));
%     %v2 = d_inv_half.*(H*(w_div_card_edg.*(H'*d_inv_half.*v2)));
% end