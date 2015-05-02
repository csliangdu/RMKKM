function [res_kkm] = KernelKmeans_single_kernel(K, y, prefix_file_name, nRepeat)

nSmp = size(K, 1);
nClass = length(unique(y));
assert(nSmp == length(y), 'The size of K and y should be the same!');

if ~exist('prefix_file_name', 'var') || isempty('prefix_file_name')
    prefix_file_name = 'some_kernel';
end

disp(['KenelKmeans ', num2str(nRepeat), 'iterations!']);
res_file = strcat(prefix_file_name, '_res_kkm_single_kernel.mat');
if exist(res_file, 'file')
    load(res_file, 'res_kkm');
else
    res_kkm = [];
    rng('default');
    for iRepeat = 1:nRepeat
        t_start = clock;
        disp(['KernelKmeans ', num2str(iRepeat), ' of ' num2str(nRepeat), ' iterations begin ...']);
        label_kkm = KernelKmeans(K, nClass, 'maxiter', 100, 'Replicates', 1);
        res_kkm = [res_kkm; ClusteringMeasure(y, label_kkm)];%#ok<AGROW>
        t_end = clock;
        disp(['KernelKmeans ', num2str(iRepeat), ' of ' num2str(nRepeat), ' iterations done.']);
        disp(['KernelKmeans exe time: ', num2str(etime(t_end, t_start))]);
    end
    save(res_file, 'res_kkm');
end
save(res_file, 'res_kkm');
end
