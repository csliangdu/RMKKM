function [label, kw, center, bCon, sumD, D, objHistory] = RMKKM(Ks, k, varargin)
%   [1] Liang Du, Peng Zhou, Lei Shi, Hanmo Wang, Mingyu Fan, Wenjian Wang, Yi-Dong Shen, 
%   Robust Multiple Kernel Kmeans, IJCAI 2015
%
%**************************************************
%     Author: Liang Du <csliangdu@gmail.com>
%     Version: 1.0
%     Last modified: 2015-04-29 04:37:36
%**************************************************

m = length(Ks);
n = size(Ks{1}, 1);
assert(n == size(Ks{1}, 2), 'The input kernel matrix should be squared');

if ~(isscalar(k) && isnumeric(k) && isreal(k) && k > 0 && (round(k)==k))
    error('RMKKM:Invalid k', 'k must be a positive integer value.');%#ok
elseif n < k
    error('RMKKM:TooManyClusters', ...
        'K must have more rows than the number of clusters.');
end

pnames = {'gamma', 'start'   'maxiter'  'replicates' 'onlinephase' };
dflts =  {'0.7', 'sample'       []        []        'off'                };
[eid,errmsg, gamma, start, maxiter, replicates] = getargs(pnames, dflts, varargin{:});
if ~isempty(eid)
    error(sprintf('RMKKM:%s',eid),errmsg);
end

center = [];
if ischar(start)
    startNames = {'sample'};
    j = find(strncmpi(start,startNames,length(start)));
    if length(j) > 1
        error(message('RMKKM:AmbiguousStart',start));
    elseif isempty(j)
        error(message('RMKKM:UnknownStart', start));
    elseif isempty(k)
        error('RMKKM:MissingK', 'You must specify the number of clusters, K.');
    end
    if j == 2
        if floor(.1*n) < 5*k
            j = 1;
        end
    end
    start = startNames{j};
end

% The maximum iteration number is default 100
if isempty(maxiter)
    maxiter = 100;
end

% Assume one replicate
if isempty(replicates)
    replicates = 1;
end

bestlabel = [];
bestkw = [];
sumD = zeros(1,k);
bCon = false;
bestObjHistory = [];

for t=1:replicates
    objHistory  = [];
    % cluster center initialization, each cluster is represented a linear combination of the data points
    switch start
        case 'sample'
            seed = randsample(n,k);
            center = zeros(n, k);
            for i = 1:k
                center(seed(i),i) = 1;
            end
            center = bsxfun(@rdivide, center, max(sum(center), 1e-10)); % weighted indicator
        case 'numeric'
            center = rand(n, k);
            center = bsxfun(@rdivide, center, max(sum(center), 1e-10)); % weighted indicator
    end
    last = 0;label=1;
    
    % kernel weight initialization
    kw = ones(m,1) / m;
    
    % kernel k-means with L21-norm
    iter = 0;
    while any(label ~= last) && iter < maxiter
        % update Kernel
        Ka = zeros(n);
        for i = 1:m
            Ka = kw(i) * Ks{i} + Ka;
        end
        
        if 1
            [label, center] = L21KKM(Ka, k, 'maxiter', 30, 'Replicates', 10);
        else
            [label, center] = L21KKM_single(Ka, center, label );
        end
        % update kernel weight
        Z = full(sparse(1:n,label,ones(n,1),n,k,n)); % indicator matrix
        A = zeros(n, m);
        for i = 1:m
            bb = sum((Ks{i} * center) .* center);
            ab = Ks{i} * center;
            D = bsxfun(@plus, -2*ab, bb);
            D = bsxfun(@plus, D, diag(Ks{i}));
            A(:, i) = sum(Z .* D, 2);
        end
        A = max(A, eps);
        Aw = bsxfun(@times, A, kw');
        D = 2 * sqrt(sum(Aw,2));
        D = 1 ./ D;
        A = bsxfun(@times, A, D);
        h = sum(A);
        kw = lp_simplex_proj(h, gamma)';
        
        obj = rmkkm_obj(Ks, kw, label, center);
        objHistory = [objHistory; obj];%#ok<AGROW>
        iter = iter + 1;
    end
    
    if iter<maxiter
        bCon = true;
    end
    if isempty(bestlabel)
        bestlabel = label;
        bestkw = kw;
        bestcenter = center;
        bestObjHistory = objHistory;
        if replicates>1
            Ka = zeros(n);
            for i = 1:m
                Ka = kw(i) * Ks{i} + Ka;
            end
            aa = full(diag(Ka));
            bb = sum((Ka * center) .* center);
            ab = Ka * center;
            D = bsxfun(@plus, aa, bb) - 2*ab;
            D(D<0) = 0;
            D = sqrt(D);
            for j = 1:k
                sumD(j) = sum(D(label==j,j));
            end
            bestsumD = sumD;
            bestD = D;
        end
    else
        Ka = zeros(n);
        for i = 1:m
            Ka = kw(i) * Ks{i} + Ka;
        end
        aa = full(diag(Ka));
        bb = sum((Ka * center) .* center);
        ab = Ka * center;
        D = bsxfun(@plus, aa, bb) - 2*ab;
        D(D<0) = 0;
        D = sqrt(D);
        for j = 1:k
            sumD(j) = sum(D(label==j,j));
        end
        if sum(sumD) < sum(bestsumD)
            bestlabel = label;
            bestkw = kw;
            bestcenter = center;
            bestsumD = sumD;
            bestD = D;
            bestObjHistory = objHistory;
        end
    end
end

label = bestlabel;
kw = bestkw;
center = bestcenter;
objHistory = bestObjHistory;
if replicates>1
    sumD = bestsumD;
    D = bestD;
elseif nargout > 3
    Ka = zeros(n);
    for i = 1:m
        Ka = kw(i) * Ks{i} + Ka;
    end
    aa = full(diag(Ka));
    bb = sum((Ka * center) .* center);
    ab = Ka * center;
    D = bsxfun(@plus, aa, bb) - 2*ab;
    D(D<0) = 0;
    D = sqrt(D);
    for j = 1:k
        sumD(j) = sum(D(label==j,j));
    end
end
end

function obj = rmkkm_obj(Ks, kw, label, center)
[n, k] = size(center);
m = length(Ks);
Ka = zeros(n);
for i = 1:m
    Ka = kw(i) * Ks{i} + Ka;
end
bb = sum((Ka * center) .* center);
ab = Ka * center;
D = bsxfun(@plus, -2*ab, bb);
Z = full(sparse(1:n,label,ones(n,1),n,k,n)); % indicator matrix
dist = sum(Z.*D, 2) + diag(Ka);
dist = sqrt(max(dist, eps));
obj = sum(dist);
end

function a = lp_simplex_proj(h, alpha)
% This function solve the following problem
%
%   \min_{a} \quad \sum_i a_i h_i = a^T h
%    s.t.    a_i >=0, \sum_i a_i^alpha = 1
%
% [1]Weighted Feature Subset Non-negative Matrix Factorization and Its Applications to Document Understanding. ICDM 2010
%

assert( (alpha > 0) && (alpha < 1), 'alpha should be (0, 1)');
t1 = 1 / alpha;
t2 = 1 / (alpha - 1);
t3 = alpha / (alpha - 1);

t4 = sum(h .^ t3);
t5 = (1 / t4)^t1;
t6 = h .^ t2;
a = t5 * t6;
end

function [label, center] = L21KKM_single(K, center, label )
n = size(K,1);
k = size(center, 2);
last = 0;
it=0;
aa = full(diag(K));
maxit = 10;
while any(label ~= last) && it<maxit
    last = label;
    
    bb = sum((K * center) .* center);
    ab = K * center;
    D = bsxfun(@plus, -2*ab, bb);
    
    [val,label] = min(D,[],2); % assign samples to the nearest centers
    val = aa + val;
    ll = unique(label);
    if length(ll) < k
        %disp([num2str(k-length(ll)),' clusters dropped at iter ',num2str(it)]);
        missCluster = 1:k;
        missCluster(ll) = [];
        missNum = length(missCluster);
        
        [~,idx] = sort(val,1,'descend');
        label(idx(1:missNum)) = missCluster;
    end
    
    minDist = max(val, eps);
    idx = minDist < 1e-10;
    sw = .5 ./ sqrt(minDist);
    if sum(~idx) > 0
        sw(idx) = mean(sw(~idx));% without this setting, the weight of data point close to cluster center will be infinity!
    end
    sw = sw/max(sw);
    
    center = full(sparse(1:n,label,sw,n,k,n)); % indicator matrix
    center = bsxfun(@rdivide, center, max(sum(center), 1e-10)); % weighted indicator
    it=it+1;
end
end

function [eid,emsg,varargout]=getargs(pnames,dflts,varargin)
%GETARGS Process parameter name/value pairs
%   [EID,EMSG,A,B,...]=GETARGS(PNAMES,DFLTS,'NAME1',VAL1,'NAME2',VAL2,...)
%   accepts a cell array PNAMES of valid parameter names, a cell array
%   DFLTS of default values for the parameters named in PNAMES, and
%   additional parameter name/value pairs.  Returns parameter values A,B,...
%   in the same order as the names in PNAMES.  Outputs corresponding to
%   entries in PNAMES that are not specified in the name/value pairs are
%   set to the corresponding value from DFLTS.  If nargout is equal to
%   length(PNAMES)+1, then unrecognised name/value pairs are an error.  If
%   nargout is equal to length(PNAMES)+2, then all unrecognised name/value
%   pairs are returned in a single cell array following any other outputs.
%
%   EID and EMSG are empty if the arguments are valid.  If an error occurs,
%   EMSG is the text of an error message and EID is the final component
%   of an error message id.  GETARGS does not actually throw any errors,
%   but rather returns EID and EMSG so that the caller may throw the error.
%   Outputs will be partially processed after an error occurs.
%
%   This utility can be used for processing name/value pair arguments.
%
%   Example:
%       pnames = {'color' 'linestyle', 'linewidth'}
%       dflts  = {    'r'         '_'          '1'}
%       varargin = {{'linew' 2 'nonesuch' [1 2 3] 'linestyle'':'}
%       [eid,emsg,c,ls,lw] = statgetargs(pnames,dflts,varargin{:})    % error
%       [eid,emsg,c,ls,lw,ur] = statgetargs(pnames,dflts,varargin{:}) % ok
%
% We always create (nparams+2) outputs:
%    one each for emsg and eid
%    nparams varargs for values corresponding to names in pnames
% If they ask for one more (nargout == nparams+3), it's for unrecognized
% names/values

%   Original Copyright 1993-2008 The MathWorks, Inc.
%   Modified by Deng Cai (dengcai@gmail.com) 2011.11.27

% Initialize some variables
emsg = '';
eid = '';
nparams = length(pnames);
varargout = dflts;
unrecog = {};
nargs = length(varargin);

% Must have name/value pairs
if mod(nargs,2)~=0
    eid = 'WrongNumberArgs';
    emsg = 'Wrong number of arguments.';
else
    % Process name/value pairs
    for j=1:2:nargs
        pname = varargin{j};
        if ~ischar(pname)
            eid = 'BadParamName';
            emsg = 'Parameter name must be text.';
            break;
        end
        i = strcmpi(pname,pnames);
        i = find(i);
        if isempty(i)
            % if they've asked to get back unrecognised names/values, add this
            % one to the list
            if nargout > nparams+2
                unrecog((end+1):(end+2)) = {varargin{j} varargin{j+1}};
                % otherwise, it's an error
            else
                eid = 'BadParamName';
                emsg = sprintf('Invalid parameter name:  %s.', pname);
                break;
            end
        elseif length(i)>1
            eid = 'BadParamName';
            emsg = sprintf('Ambiguous parameter name:  %s.', pname);
            break;
        else
            varargout{i} = varargin{j+1};
        end
    end
end

varargout{nparams+1} = unrecog;
end