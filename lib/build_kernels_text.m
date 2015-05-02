function build_kernels_text(dataset, isTF, prefix_name)
% This function is used to build four text kernels used in [1] and [2]
% where text data is processed as follows:
%	1. word process, 
%		1.1) all terms
%		1.2) no-short
%		1.3) no-stop
%		1.4) stemmed
%	2. tf-idf weight
%	3. data normalization to unit-norm
%	4. kernels (four distance function)
%		4.1) euclidean
%		4.2) cosine, linear kernel
%		4.3) jaccard coefficient
%		4.4) pearson crrelation coefficient
%	
% Input:
%		dataset, each row is a document
%		isTF, if isTF = 1, each feature is TF weighted, (isTF=0, by default)
% Outout:
%	
%
% [1] Multiple Kernel Fuzzy Clustering, IEEE Transactions on Fuzzy Systems, 2012
% [2] Affinity Aggregation for Spectral Clustering, CVPR, 2012
%

if isempty(which('tfidf.m')) || isempty(which('NormalizeFea.m')) || isempty(which('EuDist2.m')) || isempty(which('corr.m'))
	error('Some path for tfidf.m, NormalizeFea.m, EuDist2.m, corr.m is not found');
end

if ~exist('isTF', 'var')
	isTF = 0;
end

if isTF
	beNorm = 1;
	dataset = tfidf(dataset, beNorm);
end

dataset_tfidf_unit = NormalizeFea(dataset);

res_file = fullfile([prefix_name, '_linear.mat']);
if ~exist(res_file, 'file')
	K = dataset_tfidf_unit * dataset_tfidf_unit';
	K = KernelNormalize(K, 'NCW');
	K = (K + K')/2;
	save(res_file, 'K');
	clear K;
end

res_file = fullfile([prefix_name, '_euc.mat']);
if ~exist(res_file, 'file')
	bSqrt = 1;
	D = EuDist2(dataset_tfidf_unit, dataset_tfidf_unit, bSqrt);
	K = max(max(D)) - D;
	K = KernelNormalize(K, 'NCW');
	K = (K + K')/2;
	save(res_file, 'K');
	clear K;
end

res_file = fullfile([prefix_name, '_jaccard.mat']);
if ~exist(res_file, 'file')
	K = dataset_tfidf_unit * dataset_tfidf_unit';
	K = K ./ (2 - K);
	K = KernelNormalize(K, 'NCW');
	K = (K + K')/2;
	save(res_file, 'K');
	clear K;
end

res_file = fullfile([prefix_name, '_pcc.mat']);
if ~exist(res_file, 'file')
	K = corr(dataset_tfidf_unit');
	%K = KernelNormalize(K, 'NCW');
	K = (K + K')/2;
	K = real(K);
	save(res_file, 'K');
	clear K;
end

res_file = fullfile([prefix_name, '_rbm_avg.mat']);
if ~exist(res_file, 'file')
	bSqrt = 0;
	D = EuDist2(dataset_tfidf_unit, dataset_tfidf_unit, bSqrt);
	avgD = mean(D(:));
	K = exp( - D / avgD);
	K = KernelNormalize(K, 'NCW');
	
	K = (K + K')/2;
	K = real(K);
	save(res_file, 'K');
	clear K;
end
