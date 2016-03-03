function [filtered_normalized_data]=normalize_single_efroni(original_data)
%this is preprocessing according to Idan Efroni Arabidopsis article

normalized_data=cell2mat(original_data(2:end,2:end));

for i=1:length(normalized_data(1,:)),
    m=mean(normalized_data(:,i));
    normalized_data(:,i)=normalized_data(:,i)./m*50;
end

%new part
dataNorm=normalized_data';
m=max(dataNorm);
[s,k]=size(dataNorm);
m1=[];
for i=1:s,
    m1=[m1;m];
end
dataNorm=dataNorm./m1;
normalized_data=dataNorm';
%end new

filtered_normalized_data=horzcat(...
    original_data(2:end,1),...
    num2cell(normalized_data));
