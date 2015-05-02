function SC_single_data(dataset, kernel_type, nRepeat)

dir1 = pwd;
addpath(fullfile(pwd, '..', 'lib'));
addpath(fullfile(pwd, '..', 'data'));

cd(fullfile(dir1, '..', 'data'));
BuildKernels(dataset);

cd(dir1);
SC_all_affinity(dataset, kernel_type, nRepeat);
SC_equal_weight_multi_affinity(dataset, kernel_type, nRepeat);