function buildmodels(skip_timepoint,matzygall)
if nargin < 2
    skip_timepoint=-1;
    matzygall='mat';
end

timeseriesfile='D:/Work/article2/yoav files/Time course data based on Lott et al 2011 FBgn short2.xlsx';
max_repeat_fc_diff=1.5;
repeats=[1,9;2,10;3,11;4,12;5,13;6,14;7,15;8,16];
x_label=['start';'10   ';'11   ';'12   ';'13   ';...
         '14A  ';'14B  ';'14C  ';'14D  ';'end  '];
Ox_label=['10   ';'11   ';'12   ';'13   ';...
         '14A  ';'14B  ';'14C  ';'14D  '];
data_label=['CanS_GF';'OrR_GF ';'YW_GF  ';'CanS_WT';'OrR_WT ';'YW_WT  '];

if ~exist('ndata','var'),
    ndata=[];alldata=[];
end

if ~exist('expdata','var'),
    expdata=[];allexpdata=[];
end
datafile='D:/Work/article2/yoav files/RNASeq_WT_GF_F1_CanS_OrR_YW_analysis_spikingGenes_RPKM2.xlsx';
if length(expdata)<1,
    [expdata, ~ , allexpdata] = xlsread(datafile);
end

% skip_timepoint=2;   %index of datapoint to remove
showplots=0;

%normalize the timeseries data
[averaged_normalized_filtered_data,...
    ndata, alldata]=normalize_timeseries_efroni(...
 timeseriesfile,max_repeat_fc_diff,repeats, x_label, -1,...
 ndata, alldata);

if skip_timepoint>1,
    allexpdata=horzcat(averaged_normalized_filtered_data(:,1),...
        averaged_normalized_filtered_data(:,skip_timepoint+2));
    allexpdata=vertcat({'ID','NormVal'},allexpdata);
    expdata=cell2mat(averaged_normalized_filtered_data(:,skip_timepoint+2));

    [averaged_normalized_filtered_data,...
        ndata, alldata]=normalize_timeseries_efroni(...
     timeseriesfile,max_repeat_fc_diff,repeats, x_label, skip_timepoint,...
     ndata, alldata);    

    x=1:length(repeats);
    repeats=repeats(find(x~=skip_timepoint),:);
    x_label=x_label([1:skip_timepoint,skip_timepoint+2:length(x_label)],:);
end
%identify spiking genes from the timeseries set
[all_spikes_idx,mat_spikes_idx1,zyg_spikes_idx1]=discover_time_relevant_genes(...
 averaged_normalized_filtered_data);

[mat_spikes_idx,zyg_spikes_idx]=retrieve_original_matzyg_genes(...
 averaged_normalized_filtered_data,alldata);

[C,~,~]=intersect(mat_spikes_idx,mat_spikes_idx1);
mat_spikes_idx=C;
[C,~,~]=intersect(zyg_spikes_idx,zyg_spikes_idx1);
zyg_spikes_idx=C;

%normalize the experiment data
[exp_filtered_normalized_data]=normalize_single_efroni(allexpdata);

%use data
if matzygall=='zyg',
    all_spikes_idx=zyg_spikes_idx;
    length(zyg_spikes_idx)
elseif matzygall=='mat',
    all_spikes_idx=mat_spikes_idx;
    length(mat_spikes_idx)
elseif matzygall=='spk',
    all_spikes_idx=all_spikes_idx;
elseif matzygall=='all',
    all_spikes_idx=1:length(averaged_normalized_filtered_data(:,1));
else
    all_spikes_idx=all_spikes_idx;
end
averaged_normalized_filtered_data_spikes=...
    averaged_normalized_filtered_data(all_spikes_idx,:);
[C,ia,ib]=intersect(averaged_normalized_filtered_data_spikes(:,1),...
    exp_filtered_normalized_data(:,1));
% spM=exp_filtered_normalized_data;
if skip_timepoint>1,
    spM=cell2mat(exp_filtered_normalized_data(ib,2));
else
    spM=[mean(cell2mat(exp_filtered_normalized_data(ib,2:4))')',...
        mean(cell2mat(exp_filtered_normalized_data(ib,5:7))')'];
%     spM=[cell2mat(exp_filtered_normalized_data(ib,4))';...
%         cell2mat(exp_filtered_normalized_data(ib,7))']';
end
sp=cell2mat(averaged_normalized_filtered_data_spikes(ia,3:end-1));

%%%%%%%%%%%%%%%%%%% see where datapoint falls in timeseries
x=1:length(sp(1,:));
x_label1=x_label([2:length(x_label)-1],:);

% a_sp=mean(averaged_normalized_data(:,2:end));
% std_sp=std(averaged_normalized_data(:,2:end))/sqrt(length(averaged_normalized_data(:,1)));

mat_spikes_averaged_normalized_data=cell2mat(averaged_normalized_filtered_data(mat_spikes_idx,3:end-1));
zyg_spikes_averaged_normalized_data=cell2mat(averaged_normalized_filtered_data(zyg_spikes_idx,3:end-1));

a_sp_mat=mean(mat_spikes_averaged_normalized_data(:,1:end));
std_sp_mat=std(mat_spikes_averaged_normalized_data(:,1:end))/sqrt(length(mat_spikes_averaged_normalized_data(:,1)));
a_sp_zyg=mean(zyg_spikes_averaged_normalized_data(:,1:end));
std_sp_zyg=std(zyg_spikes_averaged_normalized_data(:,1:end))/sqrt(length(zyg_spikes_averaged_normalized_data(:,1)));

% avg_spikes=horzcat(a_sp_mat',a_sp_zyg');

cmap1 = redbluecmap(length(x));
figure;
mm=errorbar(x,a_sp_mat,std_sp_mat,'Color',cmap1(1,:), 'LineWidth', 3);hold on;
mz=errorbar(x,a_sp_zyg,std_sp_zyg,'Color',cmap1(end,:), 'LineWidth', 3);hold on;
xlim([x(1) x(end)]);
set(gca,'XTickLabel',x_label1);
ylim([0 1]);
xlabel('Embryonic Development cycle')
ylabel('Relative gene expression')
if skip_timepoint>1,
    sstr1=sprintf('Remove timepoint test: point %s',Ox_label(skip_timepoint,:));
    sstr2=sprintf('Using %s genes as clock',matzygall(:));
    title({sstr1,sstr2});
else
    title({'Development rate of Naive and Germ-Free embryos','collected at same timepoint'});
end
legend([mm,mz],'Maternal genes','Zygotic genes','Location','NorthEast');

barlabel= ['Germ Free (corr)';'Naive (corr)    '];
barlabel1=['Germ Free (fit)'; 'Naive (fit)    '];
barcolor=[cmap1(2,:);cmap1(end-2,:)];

for exp=1:length(spM(1,:)),
%     matzygall='mat';
    if matzygall=='mat',
        datapoint=spM(:,exp);
        bestfit=zeros(length(datapoint),1);
        for l=1:length(datapoint),
            try
                bestfit(l)=find(a_sp_mat<=datapoint(l),1,'first');
            catch
                bestfit(l)=length(a_sp_mat);
            end
        end
        bestfit2=bestfit-1;
        bestfit2(find(bestfit2<1))=1;
        bestfitR=1.0*bestfit2;
        ratio=(datapoint'-a_sp_mat(bestfit))./(a_sp_mat(bestfit2)-a_sp_mat(bestfit));
        ratio(find(isnan(ratio)))=0;
        ratio(find(isinf(ratio)))=0;
        ratio(find(ratio<0))=1;
        ratio(find(isnan(ratio)))=0;
        bestfitR=bestfitR+ratio';
        c_bestfit=mean(bestfitR);
        c_std=std(bestfitR);
%         plot(c_bestfit*ones(length(datapoint),1),  datapoint,'g.');
        c_diff=bestfitR-c_bestfit;
        c_bestfit=mean(bestfitR(find(c_diff<1.5*c_std)));
%         plot(c_bestfit*ones(length(datapoint),1),  datapoint,'Color',barcolor(exp,:));
        y = graph2d.constantline(c_bestfit, 'LineStyle','-.', 'Color',barcolor(exp,:),'LineWidth', 5);
        changedependvar(y,'x');
        x_loc = get(y, 'XData');
        y_height = get(y, 'YData');
        arrayfun(@(x,y) text(x+0.2, y-0.6,barlabel1(exp,:), 'Color', barcolor(exp,:),'fontSize',14,'fontWeight','bold'), x_loc, y_height);
        a=1;
    end
%     matzygall='zyg';
    if matzygall=='zyg',
        datapoint=spM(:,exp);
        bestfit=zeros(length(datapoint),1);
        for l=1:length(datapoint),
            try
                bestfit(l)=find(a_sp_zyg>=datapoint(l),1,'first');
            catch
                bestfit(l)=length(a_sp_mat);
            end
        end
        bestfit2=bestfit-1;
%         bestfit=bestfit(find(bestfit2>0));
%         datapoint=datapoint(find(bestfit2>0));
%         bestfit2=bestfit2(find(bestfit2>0));
        bestfit2(find(bestfit2<1))=1;
        bestfitR=1.0*bestfit2;
        ratio=(datapoint'-a_sp_zyg(bestfit2))./(a_sp_zyg(bestfit)-a_sp_zyg(bestfit2));
        ratio(find(isinf(ratio)))=0;
        ratio(find(ratio<0))=1;
        bestfitR=bestfitR+ratio';
        c_bestfit=mean(bestfitR);
        c_std=std(bestfitR);
%         plot(bestfitR,  datapoint,'g.');
        c_diff=bestfitR-c_bestfit;
        c_bestfit=mean(bestfitR(find(c_diff<1.5*c_std)));
%         plot(c_bestfit*ones(length(datapoint),1),  datapoint,'Color',barcolor(exp,:));
        y = graph2d.constantline(c_bestfit, 'LineStyle',':', 'Color',barcolor(exp,:),'LineWidth', 5);
        changedependvar(y,'x');
        x_loc = get(y, 'XData');
        y_height = get(y, 'YData');
        arrayfun(@(x,y) text(x+0.2, y-0.6,barlabel1(exp,:), 'Color', barcolor(exp,:),'fontSize',14,'fontWeight','bold'), x_loc, y_height);
        a=1;
    end
    
    
    datapoint=spM(:,exp);
    e_dist=zeros(length(sp(1,:))-1,1);
    corr_c=zeros(length(sp(1,:))-1,1);
    for i=2:length(sp(1,:)),
        corr_c(i-1)=corr(sp(:,i),datapoint);
        e_dist(i-1) = norm(sp(:,i) - datapoint);
    end


    %find two most correlated consecutive points
%     horzcat(corr_c,e_dist)
    a=1:length(e_dist)-1;b=a+1;
    kk=sum(corr_c([a(:),b(:)]'));
    [sort_pair_c,ib]=sort(kk,'descend');
    kk=sum(e_dist([a(:),b(:)]'));
    [sort_pair_e,ie]=sort(sum(e_dist([a(:),b(:)]')));
%     ib=ie;
    %decide based on correlation where our point falls
    if (abs(sort_pair_c(1)-sort_pair_c(2))/max(sort_pair_c))<0.05,
        ib=ie;
    end
%     placement=sort_pair_e(ib(1))/(sort_pair_e(ib(1))+sort_pair_e(ib(2)));
    placement=e_dist(ib(1))/(e_dist(ib(1))+e_dist(ib(1)+1));
    ttl=sprintf('Data point %s fall between %s and %s at %.1f\n',...
        barlabel(exp,:),x_label1(ib(1),:),x_label1(ib(1)+1,:),placement);

%     y = bar(ib(1)+placement,1,0.08,barcolor(exp));hold on;
    if ((matzygall=='all')|(matzygall=='spk')),
        y = graph2d.constantline(ib(1)+placement, 'LineStyle',':', 'Color',barcolor(exp,:),'LineWidth', 5);
        changedependvar(y,'x');
        x_loc = get(y, 'XData');
        y_height = get(y, 'YData');
        arrayfun(@(x,y) text(x+0.2, y-0.4,barlabel(exp,:), 'Color', barcolor(exp,:),'fontSize',14,'fontWeight','bold'), x_loc, y_height);
    end
end
figureHandle = gcf;
set(findall(figureHandle,'type','text'),'fontSize',14,'fontWeight','bold')
set(gca,'FontSize',12,'fontWeight','bold')
hold off;


