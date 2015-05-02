function L21KKM_equal_weight_multi_kernel(dataset, kernel_type, nRepeat)

data_dir = fullfile(pwd, '..', 'data');
kernel_dir = fullfile(data_dir, [dataset, '_kernel']);
file_list = dir(kernel_dir);

kernel_list = {};
iKernel = 0;
for iFile = 1:length(file_list)
    sName = file_list(iFile).name;
    if (~strcmp(sName, '.') && ~strcmp(sName, '..'))
        if ~isempty(kernel_type)
            for iType = 1:length(kernel_type)
                if ~isempty(strfind(sName, kernel_type{iType}))
                    iKernel = iKernel + 1;
                    kernel_list{iKernel} = sName; %#ok<AGROW>
                end
            end
        else
            iKernel = iKernel + 1;
            kernel_list{iKernel} = sName; %#ok<AGROW>
        end
    end
end

load(fullfile(data_dir, dataset), 'y');
if ~exist('y', 'var')
    error(['y is not found in ', dataset]);
end

res_dir = fullfile(pwd, [dataset, '_res']);
if ~exist(res_dir, 'dir')
    mkdir(res_dir);
end

disp(['Total number of Kernels: ', num2str(length(kernel_list)) '!']);

Ka = zeros(length(y)); 
nKernel = length(kernel_list);
for iKernel = 1:length(kernel_list)
    iFile = kernel_list{iKernel};
    
    clear K;
    
    load(fullfile(kernel_dir, iFile), 'K');
    Ka = Ka + K / nKernel;
end
clear K; 

disp(['L21KKM on equal weighted multi kernel begin ...']);
ew_l21kkm_res_file = fullfile(res_dir, [dataset, '_res_l21kkm_ew.mat']);
if ~exist(ew_l21kkm_res_file, 'file')
    t_start = clock;
    ew_l21kkm_res = L21KKM_single_kernel(Ka, y, fullfile(res_dir, [dataset, '_ewmk']), nRepeat);
    t_end = clock;
    disp(['L21KKM exe time: ', num2str(etime(t_end, t_start))]);
    save(ew_l21kkm_res_file, 'ew_l21kkm_res');
else
    load(ew_l21kkm_res_file, 'ew_l21kkm_res');
end
disp(['L21KKM on on equal weighted multi kernel done']);
if size(ew_l21kkm_res, 1) > 1
	res_ew_l21kkm_aio = mean(ew_l21kkm_res);
else 
	res_ew_l21kkm_aio = ew_l21kkm_res;
end
	
clear K ew_l21kkm_res;

save(fullfile(res_dir, [dataset, '_res_l21kkm_ew_kernel.mat']), 'res_ew_l21kkm_aio', 'kernel_list');
end