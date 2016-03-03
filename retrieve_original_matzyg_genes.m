function [mat_spikes_idx,zyg_spikes_idx]=retrieve_original_matzyg_genes(...
 averaged_normalized_filtered_data,alldata)

[C,ia,ib]=intersect(averaged_normalized_filtered_data(:,1),...
    alldata(:,1));

filt_data=alldata(ib,:);

maternal_idx=find(strcmp('mat',filt_data(:,3)));
[C,mat_spikes_idx,ib]=intersect(averaged_normalized_filtered_data(:,1),...
    filt_data(maternal_idx,1));

zygotic_idx=find(strcmp('zyg',filt_data(:,3)));
[C,zyg_spikes_idx,ib]=intersect(averaged_normalized_filtered_data(:,1),...
    filt_data(zygotic_idx,1));
