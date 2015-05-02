function L21KKM_single_data(dataset, kernel_type, nRepeat)

dir1 = pwd;
addpath(fullfile(pwd, '..', 'lib'));
addpath(fullfile(pwd, '..', 'data'));

cd(fullfile(dir1, '..', 'data'));
BuildKernels(dataset);

cd(dir1);
L21KKM_all_kernel(dataset, kernel_type, nRepeat);
L21KKM_equal_weight_multi_kernel(dataset, kernel_type, nRepeat);