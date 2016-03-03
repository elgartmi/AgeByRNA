function [falling_genes_idx,rising_genes_idx]=find_maternal_zygotic(...
          filtered_data,filtered_indexes,...
          maternal_idx,zygotic_idx,plot_mat_zyg,x,x_label)
fc_points=1.001;
rising_genes_idx=[];
falling_genes_idx=[];

for i=1:length(filtered_data),
    inc=1;
    for j=2:length(filtered_data(i,2:end-1)),
        if (filtered_data(i,j+1)/filtered_data(i,j))<fc_points,
            inc=0;
            break;
        end
    end
    if inc==1,
        %check that it is maternal
        if length(find(zygotic_idx==filtered_indexes(i)))>0,
            rising_genes_idx=[rising_genes_idx;i];
        end
    end
end

for i=1:length(filtered_data),
    inc=1;
    for j=2:length(filtered_data(i,2:end-1)),
        if (filtered_data(i,j)/filtered_data(i,j+1))<fc_points,
            inc=0;
            break;
        end
    end
    if inc==1,
        if length(find(maternal_idx==filtered_indexes(i)))>0,
            falling_genes_idx=[falling_genes_idx;i];
        end
    end
end

if plot_mat_zyg==1,
    cmap1 = redbluecmap(length(x));
    cmap2 = ones(size(cmap1));
    cmap3=(cmap1+cmap2)/2;
    cmap3=(cmap3+cmap2)/2;cmap3=(cmap3+cmap2)/2;cmap3=(cmap3+cmap2)/2;

    figure;
    plot(x(2:end-1),filtered_data(rising_genes_idx,2:end-1),'Color',cmap3(end-1,:));hold on;

    plot(x(2:end-1),filtered_data(falling_genes_idx,2:end-1),'Color',cmap3(2,:));hold on;

    zygotic_mean=mean(filtered_data(rising_genes_idx,2:end-1));
    zygotic_std=std(filtered_data(rising_genes_idx,2:end-1))/sqrt(length(rising_genes_idx));
    zm=plot(x(2:end-1),zygotic_mean,'Color',cmap1(end-1,:), 'LineWidth', 3);hold on;
    errorbar(x(2:end-1),zygotic_mean,zygotic_std,'Color',cmap1(end-1,:), 'LineWidth', 2)
    
    maternal_mean=mean(filtered_data(falling_genes_idx,2:end-1));
    maternal_std=std(filtered_data(falling_genes_idx,2:end-1))/sqrt(length(falling_genes_idx));
    mm=plot(x(2:end-1),maternal_mean,'Color',cmap1(2,:), 'LineWidth', 3);
    errorbar(x(2:end-1),maternal_mean,maternal_std,'Color',cmap1(2,:), 'LineWidth', 2)

    set(gca,'XLim',[2 length(x)-1]);
    set(gca,'XTick',[1:length(x)-1]);
    set(gca,'XTickLabel',x_label);

    title('Selected Maternal & Zygotic genes during stages 10-14D')
    xlabel('Embryonic Development stage')
    ylabel('Relative expression')
    legend([mm,zm],sprintf('Maternal - %d genes',length(falling_genes_idx)),...
        sprintf('Zygotic - %d genes',length(rising_genes_idx)),'Location','North')

    hold off;
end
