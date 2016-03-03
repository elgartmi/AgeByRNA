function [filtered_data,filtered_indexes]=...
    normalize_efroni(original_data)
%this is preprocessing according to Idan Efroni Arabidopsis article
thres=30;

data=original_data(:,2:end);

for i=1:length(data(1,:)),
    m=mean(data(:,i));
    data(:,i)=data(:,i)./m*50;
end

data=[original_data(:,1),data];
normalized_data=data;

maxval=max(data(:,2:end)');
i30=find(maxval>thres);
data30=data(i30,:);

max30=max(data30(:,2:end)');
min30=min(data30(:,2:end)');
min30(find(min30==0))=0.000001;
fc2=max30./min30;
data30fc2=data30(find(fc2>=2),:);

filtered_indexes=data30fc2(:,1);    %these are the indexes into original data

dataNorm=data30fc2(:,2:end)';
m=max(dataNorm);
[s,k]=size(dataNorm);
m1=[];
for i=1:s,
    m1=[m1;m];
end
dataNorm=dataNorm./m1;
dataNorm=dataNorm';

z=zeros(length(dataNorm),1);
filtered_data=[z,dataNorm,z];                %this is the normalized scaled data
