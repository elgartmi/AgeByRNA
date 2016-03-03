function [averaged_normalized_filtered_data,...
          ndata_r, alldata_r]=normalize_timeseries_efroni(...
            timeseriesfile,max_repeat_fc_diff,repeats, x_label,toremove,...
            ndata, alldata)
%%%%%%%%%%%%%%%%%%%%%%%%%%read timeseries data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(ndata)<1,
    [ndata, ~ , alldata] = xlsread(timeseriesfile);
end

ndata_r=ndata;alldata_r=alldata;

%select genes for which male=female and average male and female
male_eq_female=find((ndata(:,3)>(1-max_repeat_fc_diff))&(ndata(:,3)<(1+max_repeat_fc_diff)));
rdata=ndata(male_eq_female,4:end);
x=[1:length(repeats)+2];

% Remove one of the timepoints and use it as unknown to test calculation
if toremove>1,
    repeats=repeats([1:toremove-1,toremove+1:length(repeats)],:);
    x=[1:length(repeats)+2];
    x_label=x_label([1:toremove,toremove+2:length(x_label)],:);
end

%add indexes into original data
idx=1:length(ndata);
idx=idx(male_eq_female);
averaged_normalized_data=[idx',rdata];

%average the repeats
data=[];
for i=1:length(repeats),
    data=[data,mean(rdata(:,repeats(i,:))')'];
end
averaged_normalized_data=[idx',data];

%normalize
[averaged_normalized_filtered_data,filtered_indexes]=...
    normalize_efroni(averaged_normalized_data);

averaged_normalized_filtered_data=horzcat(...
    alldata(filtered_indexes+1,1),...
    num2cell(averaged_normalized_filtered_data));


