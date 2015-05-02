function [f, dx] = SC(w,cluster)
% Input
%       w : N x N affinity matrix
% cluster : desired number of clusters
% Output
%      dx : clustering result

%%% compute Laplacian matrix %%%
D=diag(sum(w,1));
L=D-w;
L = (L + L')/2;
L = full(L);
%%% eigen decomposition %%%
% OPTS.disp = 0;
% [f, D_] = eigs((L+L')/2, D, cluster, 'SA', OPTS);%generalized eigenproblem
[eigvec, eigval] = eig(L);
f = eigvec(:, 1:min(cluster, size(eigvec,2)));
f = NormalizeFea(f);
if nargout > 1
    if ~isempty(which('litekmeans.m'))
        dx = litekmeans(f, cluster, 'replicates', 20);
    else
        dx = kmeans(f,cluster,'EmptyAction','drop','Replicates',50);
    end
end
clear D_;
clear D;
clear L;
% clear f;
