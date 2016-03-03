function [spiking_genes_idx]=find_spiking_genes(filtered_data)
      
spiking_genes_idx=cell(1,length(filtered_data(1,:)));
for i=1:length(filtered_data),
    inc=0;
    peak=0;
    for j=2:length(filtered_data(i,1:end-1)),
        if ((filtered_data(i,j)>filtered_data(i,j+1))&&...
                (filtered_data(i,j)>filtered_data(i,j-1))),
            inc=inc+1;
            peak=j;
        end
    end
    if inc==1,
        spiking_genes_idx{peak}=[spiking_genes_idx{peak};i];
    end
end


