function [spikes_idx,mat_spikes_idx,zyg_spikes_idx...
    ]=discover_time_relevant_genes(averaged_normalized_filtered_data)

[spiking_genes_idx]=find_spiking_genes(...
          cell2mat(averaged_normalized_filtered_data(:,2:end)));
      
spikes_idx=[];
for i=1:length(spiking_genes_idx),
    spikes_idx=[spikes_idx;spiking_genes_idx{i}];
end
mat_spikes_idx=spiking_genes_idx{2};
zyg_spikes_idx=spiking_genes_idx{end-1};
