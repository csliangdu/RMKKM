function SC_all_affinity(dataset, kernel_type, nRepeat)

data_dir = fullfile(pwd, '..', 'data');
kernel_dir = fullfile(data_dir, [dataset, '_kernel']);
file_list = dir(kernel_dir);

kernel_list = {};
iAffinity = 0;
for iFile = 1:length(file_list)
    sName = file_list(iFile).name;
    if (~strcmp(sName, '.') && ~strcmp(sName, '..'))
        if ~isempty(kernel_type)
            for iType = 1:length(kernel_type)
                if ~isempty(strfind(sName, kernel_type{iType}))
                    iAffinity = iAffinity + 1;
                    kernel_list{iAffinity} = sName; %#ok<AGROW>
                end
            end
        else
            iAffinity = iAffinity + 1;
            kernel_list{iAffinity} = sName; %#ok<AGROW>
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

disp(['Total number of Affinity: ', num2str(length(kernel_list)) '!']);

res_sc_aio = [];
for iAffinity = 1:length(kernel_list)
    iFile = kernel_list{iAffinity};
    
    clear K;
    
    load(fullfile(kernel_dir, iFile), 'K');
    
    disp(['SC on ',  num2str(iAffinity), ' of ' num2str(length(kernel_list)), ' Affinity(',  iFile(1:end-4), ') begin ...']);
    sc_res_file = fullfile(res_dir, [iFile(1:end-4), '_res_sc.mat']);
    if ~exist(sc_res_file, 'file')
        t_start = clock;
        [sc_res_max, sc_res_km, sc_res_uni_km] = SC_single_affinity(K, y, fullfile(res_dir, iFile(1:end-4)), nRepeat);
        t_end = clock;
        disp(['SC exe time: ', num2str(etime(t_end, t_start))]);
        save(sc_res_file, 'sc_res_max', 'sc_res_km', 'sc_res_uni_km');
    else
        load(sc_res_file, 'sc_res_max', 'sc_res_km', 'sc_res_uni_km');
    end
    disp(['SC on ',  num2str(iAffinity), ' of ' num2str(length(kernel_list)), ' Affinity(',  iFile(1:end-4), ') done']);
    if size(sc_res_max, 1) > 1
		res_sc_aio = [res_sc_aio; [mean(sc_res_max), mean(sc_res_km), mean(sc_res_uni_km)]]; %#ok<AGROW>
	else
		res_sc_aio = [res_sc_aio; [sc_res_max, sc_res_km, sc_res_uni_km]]; %#ok<AGROW>
	end
    
    clear K sc_res_max sc_res_km sc_res_uni_km;
end

save(fullfile(res_dir, [dataset, '_res_sc_all_kernel.mat']), 'res_sc_aio', 'kernel_list');
end