function KernelKmeans_all_kernel(dataset, kernel_type, nRepeat)

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

res_kkm_aio = [];
for iKernel = 1:length(kernel_list)
    iFile = kernel_list{iKernel};
    
    clear K;
    
    load(fullfile(kernel_dir, iFile), 'K');
    
	disp(['KernelKmeans on ',  num2str(iKernel), ' of ' num2str(length(kernel_list)), ' Kernel(',  iFile(1:end-4), ') begin ...']);
    kkm_res_file = fullfile(res_dir, [iFile(1:end-4), '_res_kkm.mat']);
    if ~exist(kkm_res_file, 'file')
        t_start = clock;
        kkm_res = KernelKmeans_single_kernel(K, y, fullfile(res_dir, iFile(1:end-4)), nRepeat);
        t_end = clock;
        disp(['KernelKmeans exe time: ', num2str(etime(t_end, t_start))]);
        save(kkm_res_file, 'kkm_res');
    else
        load(kkm_res_file, 'kkm_res');
    end
    disp(['KernelKmeans on ',  num2str(iKernel), ' of ' num2str(length(kernel_list)), ' Kernel(',  iFile(1:end-4), ') done']);
	if size(kkm_res, 1) > 1
		res_kkm_aio = [res_kkm_aio; mean(kkm_res)]; %#ok<AGROW>
	else
		res_kkm_aio = [res_kkm_aio; kkm_res]; %#ok<AGROW>
	end
    clear K kkm_res;
end

save(fullfile(res_dir, [dataset, '_res_kkm_all_kernel.mat']), 'res_kkm_aio', 'kernel_list');
end