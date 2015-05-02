function [f] = BuildKernels(dataset)
global data_dir
global kernel_dir
global KernelTypes KernelPostProcessTypes PolynomialDegrees PolyPlusDegrees GaussianDegrees 
global kernel_path_part_1
%create a folder named by the name of dataset

data_dir = fullfile(pwd, '..', 'data');

kernel_dir = fullfile(data_dir, [dataset, '_kernel']);
if exist(kernel_dir, 'dir') == 0
    mkdir(kernel_dir);
end

load(fullfile(data_dir, dataset));
if exist('X', 'var')
    % error([dataset, ' do not have variable X!']);
    is_mv = 0;
elseif exist('views', 'var')
        % error([dataset, ' do not have variable views!']);
        is_mv = 1;
        for iView = 1:length(views)
            if eval(sprintf('exist(''%s'', ''var'') == 0', ['X', num2str(iView)]))
                eval(sprintf('error([dataset, '' do not have variable %s!'']);', ['X', num2str(iView)]));
            end
        end
else 
    error([dataset, ' do not have variable X or views!']);
end

KernelTypes = {'Linear', 'PolyPlus', 'Polynomial', 'Gaussian'};
KernelPostProcessTypes = {'Sample-Scale'};
PolynomialDegrees = [2, 4];
PolyPlusDegrees = [2, 4];
GaussianDegrees = [0.01, 0.05, 0.1, 1, 10, 50, 100];


if is_mv
    for iView = 1:length(views)
        kernel_path_part_1 = [dataset, '_kernel_view', num2str(iView) '_'];
        clear X;
        eval(sprintf('X = %s;', ['X', num2str(iView)]));
        BuildSingleKernels(X);
    end
else
    kernel_path_part_1 = [dataset, '_kernel_'];
    BuildSingleKernels(X);
end

f = 1;
end

function f2 = BuildSingleKernels(X)
global kernel_dir
global KernelTypes KernelPostProcessTypes PolynomialDegrees PolyPlusDegrees GaussianDegrees 
global kernel_path_part_1

for kernel_type = KernelTypes
    kernel_option = [];
    switch lower(kernel_type{1})
        case lower('Linear')
            kernel_option.KernelType = 'Linear';
            kernel_path_part_2 = 'linear_';
            for iPost = KernelPostProcessTypes
                kernel_path_part_3 = ['post_', iPost{1}];
                kernel_file = fullfile(kernel_dir, strcat(kernel_path_part_1, kernel_path_part_2, kernel_path_part_3, '.mat'));
                if ~exist(kernel_file, 'file')
                    K0 = constructKernel(X, [], kernel_option);
                    K = KernelNormalize(K0, iPost{1});%#ok
                    save(kernel_file, 'K');
                end
            end
        case lower('Polynomial')
            kernel_option.KernelType = 'Polynomial';
            for iKernelParam = PolynomialDegrees
                kernel_option.d = iKernelParam;
                kernel_path_part_2 = ['polynomial_', num2str(iKernelParam), '_'];
                for iPost = KernelPostProcessTypes
                    kernel_path_part_3 = ['post_', iPost{1}];
                    kernel_file = fullfile(kernel_dir, strcat(kernel_path_part_1, kernel_path_part_2, kernel_path_part_3, '.mat'));
                    if ~exist(kernel_file, 'file')
                        K0 = constructKernel(X, [], kernel_option);
                        K = KernelNormalize(K0, iPost{1});%#ok
                        save(kernel_file, 'K');
                    end
                end
            end
        case lower('PolyPlus')
            kernel_option.KernelType = 'PolyPlus';
            for iKernelParam = PolyPlusDegrees
                kernel_option.d = iKernelParam;
                kernel_path_part_2 = ['polyplus_', num2str(iKernelParam), '_'];
                for iPost = KernelPostProcessTypes
                    kernel_path_part_3 = ['post_', iPost{1}];
                    kernel_file = fullfile(kernel_dir, strcat(kernel_path_part_1, kernel_path_part_2, kernel_path_part_3, '.mat'));
                    if ~exist(kernel_file, 'file')
                        K0 = constructKernel(X, [], kernel_option);
                        K = KernelNormalize(K0, iPost{1});%#ok
                        save(kernel_file, 'K');
                    end
                end
            end
        case lower('Gaussian')
            kernel_option.KernelType = 'Gaussian';
            for iKernelParam = GaussianDegrees
                kernel_option.t = iKernelParam;
                kernel_path_part_2 = ['gaussian_', num2str(iKernelParam), '_'];
                for iPost = KernelPostProcessTypes
                    kernel_path_part_3 = ['post_', iPost{1}];
                    kernel_file = fullfile(kernel_dir, strcat(kernel_path_part_1, kernel_path_part_2, kernel_path_part_3, '.mat'));
                    if ~exist(kernel_file, 'file')
                        D = EuDist2(X, [], 0);
                        max_D = max(D(:));
                        max_D = sqrt(max_D);
                        kernel_option.t = kernel_option.t * max_D;
                        K0 = constructKernel(X, [], kernel_option);
                        % K0 = exp(- D / (2 * kernel_option.t^2) );
                        K = KernelNormalize(K0, iPost{1});%#ok
                        save(kernel_file, 'K');
                    end
                end
            end
        case lower('text')
            kernel_path_part_2 = ['text'];
            t1 = sum(X(:));
            if fix(t1) == t1
                isTF = 1;
            else
                isTF = 0;
            end
            kernel_file = fullfile(kernel_dir, strcat(kernel_path_part_1, kernel_path_part_2));
            build_kernels_text(X, isTF, kernel_file);
        otherwise
            error('KernelType does not exist!');
    end
end
f2 = 1;
end