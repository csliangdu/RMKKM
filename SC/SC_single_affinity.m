function [res_max, res_km, res_uni_km, obj_final] = SC_single_affinity(K, y, prefix_file_name, nRepeat)

nSmp = size(K, 1);
nClass = length(unique(y));
assert(nSmp == length(y), 'The size of K and y should be the same!');

if ~exist('prefix_file_name', 'var') || isempty('prefix_file_name')
    prefix_file_name = 'some_affinity';
end

res_max = [];
res_km = [];
res_uni_km = [];
obj_final = [];

disp(['SC on ', prefix_file_name, ' ', num2str(nRepeat), 'iterations!']);

res_file = strcat(prefix_file_name, '_res_sc_single_affinity.mat');

if exist(res_file, 'file')
    load(res_file, 'res_max', 'res_km', 'res_uni_km', 'obj_final');
else
	rng('default');
	V = SC(K, nClass);
	for iRepeat = 1:nRepeat
		t_start = clock;
		disp(['SC ',  num2str(iRepeat), ' of ' num2str(nRepeat), ' iterations begin ...']);
		res_max = [res_max; [0, 0, 0]];%#ok<AGROW>
		label_km = litekmeans(V, nClass, 'maxIter', 1000, 'Replicates', 10);
		res_km = [res_km; ClusteringMeasure(y, label_km)];%#ok<AGROW>
		label_uni_km = litekmeans(NormalizeFea(V, 1), nClass, 'maxIter', 100, 'Replicates', 10);
		res_uni_km = [res_uni_km; ClusteringMeasure(y, label_uni_km)];%#ok<AGROW>
		obj_final = [obj_final; 0];%#ok<AGROW>
		t_end = clock;
		disp(['SC ',  num2str(iRepeat), ' of ' num2str(nRepeat), ' iterations done.']);
		disp(['SC exe time: ', num2str(etime(t_end, t_start))]);
	end
    save(res_file, 'res_max', 'res_km', 'res_uni_km', 'obj_final');
end