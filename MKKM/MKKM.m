function [idx, weight, obj_final]=MKKM(U,cluster,degree,error,kvalue)
% Input
%       U : initial membership matrix
% cluster : desired number of clusters
%  degree : fuzzification degree
%   error : stop threshold
%  kvalue : N x N x k affinity matrices
% 
% Output
% U_final : membership matrix
%  weight : weight assignment to affinity matrices
%
% [1] Multiple Kernel Fuzzy Clustering, IEEE Transactions on Fuzzy Systems,
%
%  Note that multiple kernel c-means is provided by the author at       
%     http://imp.iis.sinica.edu.tw/IVCLab/research/Sean/mkfc/code.rar 
%  The multiple kernel kernel kmeans is modifed according to algorithm 3 in [1], 
%     where U is restricted to be an indicator matrix.
% 

data_n = size(kvalue, 1);
dimension=size(kvalue,3);
if isempty(degree)
    degree = 1;
end
if isempty(error)
    error = 1e-5;
end
% Main loop
U_final=U;
for iter=1:100
    bk=zeros(1,dimension);
    alpha=zeros(data_n,cluster,dimension);
    %calculate normalized memberships
    mf = U_final.^degree;
    mf_tmp = mf*diag(1./sum(mf + eps));
    for k=1:dimension
        alpha(:,:,k)=ones(data_n,1)*diag(mf_tmp'*kvalue(:,:,k)*mf_tmp)'-2*kvalue(:,:,k)*mf_tmp+1;
    end
    %calculate coefficients bk
    for k=1:dimension
        bk(1,k)=sum(sum(mf.*alpha(:,:,k)));
    end
    %update weights w
    w=ones(1,dimension)./bk;
    w=w/sum(w);
    %calculate distances D
    dist=zeros(data_n,cluster);
    wtmp=w.^2;
    for k=1:dimension
        dist=dist+alpha(:,:,k)*wtmp(1,k);
    end
    
    %update memberships U
    
    % tmp = dist.^(-1/(degree-1));
    % U_final = tmp./(sum(tmp,2)*ones(1,cluster));
    U_final = zeros(data_n, cluster);
    [~, idx] = min(dist, [], 2);
    U_final(sub2ind([data_n, cluster], [1:data_n]', idx)) = 1;
    
    % objective function
    obj_fcn(iter) = sum(sum(dist.*mf));
    U_old=U_final;
    % check termination condition
    if iter > 1
        if abs(obj_fcn(iter) - obj_fcn(iter-1))< error
            break;
        end
    end
end
obj_final = obj_fcn(end);
[~, idx] = max(U_final, [], 2);
weight=w;
clear w;
clear mf;
clear mf_tmp;
clear tmp;
clear obj_fcn;
clear dist;
clear bk;
