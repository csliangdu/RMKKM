function AASC_multi_kernel(dataset, kernel_type, nRepeat)

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

aasc_res_max = [];
aasc_res_km = [];
aasc_res_uni_km = [];
obj_final = [];

disp(['Total number of Kernels: ', num2str(length(kernel_list)) '!']);

disp(['AASC on multi kernel begin ...']);

aasc_res_file = fullfile(res_dir, [dataset, '_res_aasc.mat']);
if exist(aasc_res_file, 'file')
    load(aasc_res_file, 'aasc_res_max', 'aasc_res_km', 'aasc_res_uni_km', 'obj_final');
else
	rng('default');
	for iRepeat = 1:nRepeat
		t_start = clock;
		disp(['AASC ',  num2str(iRepeat), ' of ' num2str(nRepeat), ' iterations begin ...']);
		
		V = AASC(Ks, nClass);
		aasc_res_max = [aasc_res_max; [0,0,0]];%#ok<AGROW>
		label_km = litekmeans(V, nClass, 'maxIter', 1000, 'Replicates', 10);
		aasc_res_km = [aasc_res_km; ClusteringMeasure(y, label_km)];%#ok<AGROW>
		label_uni_km = litekmeans(NormalizeFea(V, 1), nClass, 'maxIter', 1000, 'Replicates', 10);
		aasc_res_uni_km = [aasc_res_uni_km; ClusteringMeasure(y, label_uni_km)];%#ok<AGROW>
		obj_final = [obj_final; 0];%#ok<AGROW>
		t_end = clock;
		disp(['AASC ',  num2str(iRepeat), ' of ' num2str(nRepeat), ' iterations done.']);
		disp(['AASC exe time: ', num2str(etime(t_end, t_start))]);
	end
    save(aasc_res_file, 'aasc_res_max', 'aasc_res_km', 'aasc_res_uni_km', 'obj_final');
end

disp(['AASC on multi kernel done']);
if size(aasc_res_max,1) > 1
	res_aasc_aio = [mean(aasc_res_max), mean(aasc_res_km), mean(aasc_res_uni_km)];
else
	res_aasc_aio = [aasc_res_max, aasc_res_km, aasc_res_uni_km];
end
clear K aasc_res_max aasc_res_km aasc_res_uni_km;

save(fullfile(res_dir, [dataset, '_res_aasc_multi_kernel.mat']), 'res_aasc_aio', 'kernel_list');
end