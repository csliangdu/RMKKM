clear;clc;

addpath(fullfile(pwd, '..', 'lib'));

ds = {'YALE_165n_1024d_15c_zscore_uni', 'jaffe_213n_676d_10c_uni', 'ORL_400n_1024d_40c_zscore_uni',...
    'AR_840n_768d_120c_uni', 'COIL20_1440n_1024d_20c', ...
    'tr11_414n_6429d_9c_tfidf_uni',  'tr41_878n_7454d_10c_tfidf_uni', 'tr45_690n_8261d_10c_tfidf_uni',...
    };

for i =1:length(ds)
    run_single_dataset(ds{i});
end

% estimate the best gamma
[best_gamma, best_gamma_idx] = plot_rmkkm_parameter_sensitity_on_gamma(ds);
% best_gamma = 0.4;

res_acc_aio = [];
res_nmi_aio = [];
res_purity_aio = [];
for iData = 1:length(ds)
    clear table_acc table_nmi table_purity
    load(fullfile([ds{iData} '_res'], [ds{iData} '_res_table.mat']));
    % rmkkm with fixed parameter
    if exist('res_rmkkm_p', 'var')
        table_acc(end) = res_rmkkm_p(best_gamma_idx,1);
        table_nmi(end) = res_rmkkm_p(best_gamma_idx,2);
        table_purity(end) = res_rmkkm_p(best_gamma_idx,3);
    end
    res_acc_aio = [res_acc_aio; table_acc];
    res_nmi_aio = [res_nmi_aio; table_nmi];
    res_purity_aio = [res_purity_aio; table_purity];
end

idx = [1,2,4,5,7,8, 10:15];
res_acc_aio = res_acc_aio(:, idx);
res_nmi_aio = res_nmi_aio(:, idx);
res_purity_aio = res_purity_aio(:, idx);

res_ijcai15_aio = [];
for i1 = 1:size(res_purity_aio,1)
    res_ijcai15_aio = [res_ijcai15_aio; res_acc_aio(i1,:)];
    res_ijcai15_aio = [res_ijcai15_aio; res_nmi_aio(i1,:)];
    res_ijcai15_aio = [res_ijcai15_aio; res_purity_aio(i1,:)];
end
res_sort_val = sum(res_ijcai15_aio);
res_ijcai15_aio = [res_ijcai15_aio; mean(res_acc_aio)];
res_ijcai15_aio = [res_ijcai15_aio; mean(res_nmi_aio)];
res_ijcai15_aio = [res_ijcai15_aio; mean(res_purity_aio)];

z = zeros(1, size(res_ijcai15_aio,2));
[~, idx] = sort(res_sort_val, 'descend');
z(idx) = [1:length(idx)];
res_ijcai15_aio = [res_ijcai15_aio; z];

save('rmkkm_ijcai_res_aio.mat', 'res_acc_aio', 'res_nmi_aio', 'res_purity_aio', 'res_ijcai15_aio', 'best_gamma', 'ds');

rowLabels = {'YALE', 'YALE', 'YALE', 'JAFFE', 'JAFFE', 'JAFFE', 'ORL', 'ORL','ORL', ...
    'AR', 'AR', 'AR', 'COIL20', 'COIL20', 'COIL20', 'TR11', 'TR11', 'TR11',...
    'TR41', 'TR41', 'TR41', 'TR45', 'TR45', 'TR45', ...
    'Average', 'Average', 'Average', 'Rank'};
columnLabels = {'KKMb', 'KKMa', 'SCb', 'SCa', 'RKKMb', 'RKKMa', 'KKMew', 'SCew', 'RKKMew', 'MKKM', 'AASC', 'RMKKM'};
matrix2latex(res_ijcai15_aio, 'rmkkm_ijcai_res_aio.tex', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%4.4f','size', 'tiny');
rmpath(fullfile(pwd, '..', 'lib'));