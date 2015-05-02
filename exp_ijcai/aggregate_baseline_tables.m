function aggregate_baseline_tables(dataset, kernel_type)
%
prefix = fullfile(pwd, [dataset, '_res'], [dataset, '_res']);
apps = {'kkm', 'sc', 'l21kkm'};


% load results with single kernel
r_s = [];
for app = apps
    res_file = strcat(prefix, '_', app{1}, '_all_kernel.mat');
    if exist(res_file, 'file');
        eval(sprintf('load(res_file, ''%s'', ''kernel_list'');', ['res_', app{1}, '_aio']));
    end
end
k_idx = [];
if iscell(kernel_type)
    for i = 1:length(kernel_type)
        k_idx = [k_idx; find(~cellfun(@isempty, strfind(kernel_list, kernel_type{i})))];
    end
end
k_idx = unique(k_idx);

r_s = [res_kkm_aio(k_idx,:); res_sc_aio(k_idx,7:9); res_l21kkm_aio(k_idx,:)];

% load results with equal weighted kernel
for app = apps
    res_file = strcat(prefix, '_', app{1}, '_ew_kernel.mat');
    if exist(res_file, 'file');
        eval(sprintf('load(res_file, ''%s'', ''kernel_list'');', ['res_ew_', app{1}, '_aio']));
    end
end
r_ew = [res_ew_kkm_aio; res_ew_sc_aio(7:9); res_ew_l21kkm_aio];

% load results with multiple kernel learning
mkapps = {'mkkm', 'aasc', 'rmkkm'};%
for mkapp = mkapps
    res_file = strcat(prefix, '_', mkapp{1}, '_multi_kernel.mat');
    if exist(res_file, 'file');
        if strcmp(mkapp, 'rmkkm')
            eval(sprintf('load(res_file, ''%s'', ''kernel_list'', ''res_aio_p'');', ['res_', mkapp{1}, '_aio']));
            res_rmkkm_p = res_aio_p;
        else
            eval(sprintf('load(res_file, ''%s'', ''kernel_list'');', ['res_', mkapp{1}, '_aio']));
        end
    end
end
r_m = [res_mkkm_aio; res_aasc_aio(7:9); res_rmkkm_aio];%

% generate acc table
%
% c_measures = {'Accuracy', 'Normalized Mutual Information', 'Purity'};
c_measures = {'acc', 'nmi', 'purity'};
myled = {'KKM-b', 'KKM-a', 'KKM-w', 'SC-b', 'SC-a', 'SC-w', ...
    'L21KKM-b','L21KKM-a','L21KKM-w', 'e-KKM', 'e-SC', 'e-L21KKM',...
    'MKKM', 'AASC', 'RMKKM'%
	};
for idx_m = 1:length(c_measures)
    
    res_s_s = reshape(r_s(:,idx_m), size(r_s,1)/length(apps), length(apps));
    table_single = [];
    for idx_app = 1:length(apps)
        table_single = [table_single, [max(res_s_s(:, idx_app)), mean(res_s_s(:, idx_app)), min(res_s_s(:, idx_app))]];
    end
    table_ew = r_ew(:, idx_m)';
    table_m = r_m(:, idx_m)';
    eval(sprintf('%s = [table_single, table_ew, table_m];', ['table_', c_measures{idx_m}]));
    
    % figure;
    x = (1:length(table_m) + length(table_ew) + length(table_single));
    eval(sprintf('y = %s;', ['table_', c_measures{idx_m}]));
    b=bar(x, y);
    ch = get(b,'children');
    rng('default');
    set(ch,'FaceVertexCData',randperm(length(x))');
    for i = 1:length(x)
        text(x(i),y(i)+0.02,num2str(y(i)));
    end
    title(['Results on ', dataset]);
    ylabel(c_measures{idx_m});
    set(gca, 'XTick', 1:length(x));
    % set(gca, 'XTickLabel', myled);
    % 获取xticklabel的值
    xtl=get(gca,'XTickLabel');
    % 获取xtick的值
    xt=get(gca,'XTick');
    % 获取xtick的值
    yt=get(gca,'YTick');
    % 设置text的x坐标位置们
    xtextp=xt;
    % 设置text的y坐标位置们
    ytextp=(yt(1)-0.2*(yt(2)-yt(1)))*ones(1,length(xt));
    % rotation，正的旋转角度代表逆时针旋转，旋转轴可以由HorizontalAlignment属性来设定，
    % 有3个属性值：left，right，center
    if length(myled) == length(y)
        text(xtextp,ytextp,myled,'HorizontalAlignment','right','rotation',45,'fontsize',9);
    end
    % 取消原始ticklabel
    set(gca,'xticklabel','');
    save2pdf([dataset, '_', c_measures{idx_m}, '_res.pdf'], gcf, 1200);
end
prefix = fullfile(pwd, [dataset, '_res'], [dataset, '_res']);
res_file = [prefix '_table.mat'];
save(res_file, 'table_acc', 'table_nmi', 'table_purity');
if exist('res_rmkkm_p', 'var')
    save(res_file, 'res_rmkkm_p', '-append');
end