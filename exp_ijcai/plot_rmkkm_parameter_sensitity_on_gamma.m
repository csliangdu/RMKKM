function [best_gamma, best_gamma_idx] = plot_rmkkm_parameter_sensitity_on_gamma(ds)
addpath(fullfile(pwd, '..', 'lib'));
addpath(fullfile(pwd, '..', 'data'));

res_aio_p_agg = 0;
for iData = 1:length(ds)
    clear res_aio_p;
    fn = fullfile([ds{iData} '_res'], [ds{iData} '_res_rmkkm_multi_kernel.mat']);
    if exist(fn, 'file')
        load(fn);
        res_aio_p_agg = res_aio_p_agg + res_aio_p;
        z = res_aio_p(1:1:end,:);
        idx = [1:1:size(z,1)];
        z = z(idx,:);
        x = [1:length(idx)]';
        %idx = idx/2;
        yy = z(idx,2);
        plot(x, z(idx,2), '-ro', 'Linewidth', 1.5, 'Markersize', 8);
        hold on;
        
        clear res_ew_kkm_aio;
        load(fullfile([ds{iData} '_res'], [ds{iData} '_res_kkm_ew_kernel.mat']));
        z2 = res_ew_kkm_aio;
        yy = [yy; z2(2)];
        plot(x, z2(2)*ones(length(x),2), '--gd', 'Linewidth', 1.5, 'Markersize', 8);
        hold on;
        
        clear res_kkm_aio;
        load(fullfile([ds{iData} '_res'], [ds{iData} '_res_kkm_all_kernel.mat']));
        z2 = mean(res_kkm_aio);
        yy = [yy; z2(2)];
        p = plot(x, z2(2)*ones(length(x),2), '-.bs', 'Linewidth', 1.5, 'Markersize', 8);
        hold off;
        
        myleds = {'RMKKM', 'KKM-ew', 'KKM-a'};
        legend(myleds, 'location','southeast');
        if length(x) > 10
            x_interval = {'0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95'};
        else
            x_interval = {'0.1', '0.2', '0.3', '0.4','0.5','0.6','0.7','0.8','0.9'};
        end
        xlabel('The effect of the parameter \gamma');
        ylabel('NMI');
        t1 = ds{iData};
        t2 = strfind(t1, '_');
        if strfind(t1, 'webbb'); t1 = 'webkb'; end %...
        title(upper(t1(1:t2(1)-1)));
        yy = [min(yy) - 0.02, max(yy)+0.02];
        ylim(yy);
        xlim([1,length(x)]);
        set(gca, 'XTick', [1:length(x)]);
        set(gca, 'XTickLabel',x_interval);
        
        save2pdf(['gamma_', ds{iData}, '_res.pdf'], gcf, 1200);
    end
end

x = [1:size(res_aio_p_agg,1)];
y = sum(res_aio_p_agg, 2);
b=bar(x, y);
ch = get(b,'children');
rng('default');
set(ch,'FaceVertexCData',randperm(length(x))');
for i = 1:length(x)
    text(x(i),y(i)+0.02,num2str(y(i)));
end
title(['Aggregated Results on All Data Sets']);
ylabel('ACC+NMI+Purity');
% set(gca, 'XTick', 1:length(x));
set(gca,'xticklabel',{'0.1', '0.2', '0.3', '0.4','0.5','0.6','0.7','0.8','0.9'});

[~, best_gamma_idx] = max(y);
gammaCandidates = (0.1:0.1:0.9);
best_gamma = gammaCandidates(best_gamma_idx);