function MKKM_multi_kernel(dataset, kernel_type, nRepeat)
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
nClass = length(unique(y));

res_dir = fullfile(pwd, [dataset, '_res']);
if ~exist(res_dir, 'dir')
    mkdir(res_dir);
end

nKernel = length(kernel_list);
Ks = zeros(length(y), length(y), nKernel);

for iKernel = 1:length(kernel_list)
    iFile = kernel_list{iKernel};
    
    clear K;
    
    load(fullfile(kernel_dir, iFile), 'K');
    Ks(:,:,iKernel) = K;
end
clear K;

mkkm_res = [];
obj_final = [];
kw_aio = cell(nRepeat, 1);
disp(['Total number of Kernels: ', num2str(length(kernel_list)) '!']);

disp(['MKKM on multi kernel begin ...']);
mkkm_res_file = fullfile(res_dir, [dataset, '_res_mkkm.mat']);
if exist(mkkm_res_file, 'file')
    load(mkkm_res_file, 'res_mkkm_aio');
else
    rng('default');
    for iRepeat = 1:nRepeat
        t_start = clock;
        disp(['MKKM ',  num2str(iRepeat), ' of ' num2str(nRepeat), ' iterations begin ...']);
        U = rand(size(Ks,1), nClass);
        [~, uidx] = max(U, [], 2);
        U = zeros(size(U));
        U(sub2ind(size(U), [1:size(U,1)]', uidx)) = 1;
        [label_mkkm, kw, obj] = MKKM(U, nClass, 1, 1e-5, Ks);
        mkkm_res = [mkkm_res; ClusteringMeasure(y, label_mkkm)];%#ok<AGROW>
        obj_final = [obj_final; obj];%#ok<AGROW>
        kw_aio{iRepeat} = kw;
        t_end = clock;
        disp(['MKKM ',  num2str(iRepeat), ' of ' num2str(nRepeat), ' iterations done.']);
        disp(['MKKM exe time: ', num2str(etime(t_end, t_start))]);
    end
    if size(mkkm_res, 1) > 1
        [~, minIdx] = min(obj_final);
        mkkm_res_obj = mkkm_res(minIdx,:);
        kw_obj = kw_aio{iRepeat};
        mkkm_res_mean = mean(mkkm_res);
    else
        mkkm_res_mean = mkkm_res;
    end
    save(mkkm_res_file,  'mkkm_res_mean');
    disp(['MKKM on multi kernel done']);
    
    res_mkkm_aio = mkkm_res_mean;
    
    clear Ks K mkkm_res_mean mkkm_res_obj;
    save(fullfile(res_dir, [dataset, '_res_mkkm_multi_kernel.mat']), 'res_mkkm_aio', 'kernel_list');
end

end