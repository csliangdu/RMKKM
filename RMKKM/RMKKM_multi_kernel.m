function RMKKM_multi_kernel(dataset, kernel_type, nRepeat)
gammaCandidates = (0.1:0.1:0.9);

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
Ks = cell(nKernel,1);

for iKernel = 1:length(kernel_list)
    iFile = kernel_list{iKernel};
    
    clear K;
    
    load(fullfile(kernel_dir, iFile), 'K');
    Ks{iKernel} = K;
end
clear K;

disp(['Total number of Kernels: ', num2str(length(kernel_list)) '!']);

disp('RMKKM on multi kernel begin ...');

res_rmkkm_file = fullfile(res_dir, [dataset, '_res_rmkkm_multi_kernel.mat']);

if exist(res_rmkkm_file, 'file')
    load(res_rmkkm_file, 'res_rmkkm_aio', 'kernel_list', 'best_gamma', 'res_aio_p','best_res_rmkkm_mean');
else
    res_aio_p = [];
    best_val = 0;
    for gammaIdx = 1:length(gammaCandidates)
        gamma = gammaCandidates(gammaIdx);
        
        res_rmkkm = [];
        res_rmkkm_obj = [];
        res_rmkkm_mean = [];
        obj_final = [];
        kw_rmkkm = cell(nRepeat, 1);
        
        res_file_gamma = fullfile(res_dir, strcat(dataset, '_res_rmkkm_gamma=', num2str(gamma), '.mat'));
        if iscell(res_file_gamma); res_file_gamma = res_file_gamma{1}; end
        if exist(res_file_gamma, 'file')
            load(res_file_gamma, 'res_rmkkm', 'obj_final', 'kw_rmkkm');
        else
            rng('default');
            for iRepeat = 1:nRepeat
                t_start = clock;
                disp(['RMKKM with gamma = ', num2str(gamma), ' ', num2str(iRepeat), ' of ' num2str(nRepeat), ' iterations begin ...']);
                
                [label_rmkkm, kw, ~, ~, obj] = RMKKM(Ks, nClass, 'gamma', gamma, 'maxiter', 50, 'replicates', 1);
                res_rmkkm = [res_rmkkm; ClusteringMeasure(y, label_rmkkm)];%#ok<AGROW>
                obj_final = [obj_final; obj];%#ok<AGROW>
                kw_rmkkm{iRepeat} = kw;
                t_end = clock;
                disp(['RMKKM with gamma = ', num2str(gamma), ' ', num2str(iRepeat), ' of ' num2str(nRepeat), ' iterations done']);
                disp(['RMKKM exe time: ', num2str(etime(t_end, t_start))]);
            end
            save(res_file_gamma, 'res_rmkkm', 'obj_final', 'kw_rmkkm', 'gamma');
        end
        if nRepeat > 1
            [~, minIdx] = min(obj_final);
            res_rmkkm_obj = [res_rmkkm_obj; res_rmkkm(minIdx(1),:)];%#ok
            kw_obj = kw_rmkkm{minIdx(1)};
            res_rmkkm_mean = mean(res_rmkkm);
        else
            res_rmkkm_mean = [res_rmkkm_mean; res_rmkkm];%#ok
        end
        res_aio_p = [res_aio_p; res_rmkkm_mean];%#ok
        
        if sum(res_rmkkm_mean(:)) > best_val
            best_gamma = gamma;%#ok
            best_res_rmkkm_mean = res_rmkkm_mean;
            if nRepeat > 1
                best_kw_obj = kw_obj;%#ok
                best_res_rmkkm_obj = res_rmkkm_obj;%#ok
            end
        end
    end
    res_rmkkm_aio = best_res_rmkkm_mean;%#ok
   
    clear Ks K res_rmkkm;
    
    save(res_rmkkm_file, 'res_rmkkm_aio', 'kernel_list', 'best_gamma', 'res_aio_p','best_res_rmkkm_mean');
end

 disp('RMKKM on multi kernel done');
end