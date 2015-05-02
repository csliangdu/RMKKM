function AASC_single_data(dataset, kernel_type, nRepeat)

dir1 = pwd;
addpath(fullfile(pwd, '..', 'lib'));
addpath(fullfile(pwd, '..', 'data'));

cd(fullfile(dir1, '..', 'data'));
BuildKernels(dataset);

cd(dir1);
AASC_multi_kernel(dataset, kernel_type, nRepeat);