%% ANOVA for group data
% ZKA Oct 2015
% having run /home/zagnew/data_analysis_code/pitch_pert_stats/get_data_cTBS.m
clear all
close all
set_params_cTBS;

cd /Users/zagnew/cerebellarTBS/raw_data/
load /Users/zagnew/cerebellarTBS/data_analysis/groupdata/GROUPDATA.mat

% for plottin' 
frame_taxis = gp.precereb.patient_dat.frame_taxis;
tlims4anal_req = [-0.2 1.0];
ilims4anal = dsearchn(frame_taxis',tlims4anal_req')';
idxes4anal = ilims4anal(1):ilims4anal(2);
tlims4anal = frame_taxis(ilims4anal);

flip_perts; 
% calculates gp.cereb_diff_eachsub which is POST - PRE. Here a postitive
% value indicates an increase in comp following stimulation. A negative
% value indicates that compensation was smaller following stimulation. 
% also calculates cereb_peakdiff_eachsub which is the peak diff in pre and
% post
% gp.cereb_peakdiff_eachsub(isubj)=max(gp.cereb_diff_eachsub(isubj,:));

calc_stdev;

%% plot all post stim trials
figure
subplot(211)
x=1:size(nanmean(gp.meanpitchpertresp_precereb),2);
shadedErrorBar(x,mean(gp.meanpitchpertresp_precereb), (std(gp.meanpitchpertresp_precereb)/sqrt(13)),{'-','LineWidth', 1.5,'color',[0.8 0 0.2]}, 0.2);
hold
shadedErrorBar(x,mean(gp.meanpitchpertresp_presham), (std(gp.meanpitchpertresp_presham)/sqrt(13)),{'-','LineWidth', 1.5,'color',[0.1 0.1 0.4]}, 0.2);
goodplot
axis([0 500 -20 35])

subplot(212)
x=1:size(nanmean(gp.meanpitchpertresp_precereb),2);
shadedErrorBar(x,mean(gp.meanpitchpertresp_postcereb), (std(gp.meanpitchpertresp_postcereb)/sqrt(13)),{'-','LineWidth', 1.5,'color',[0.8 0 0.2]}, 0.2);
hold
shadedErrorBar(x,mean(gp.meanpitchpertresp_postsham), (std(gp.meanpitchpertresp_presham)/sqrt(13)),{'-','LineWidth', 1.5,'color',[0.1 0.1 0.4]}, 0.2);
goodplot
axis([0 500 -20 35])
print(gcf, '-dpdf', '-r150', '/Users/zagnew/cerebellarTBS/data_analysis/figures/wholetrialcomp.pdf');

% or
figure
subplot(121)
x=1:size(nanmean(gp.meanpitchpertresp_precereb),2);
shadedErrorBar(x,mean(gp.meanpitchpertresp_precereb), (std(gp.meanpitchpertresp_precereb)/sqrt(13)),{'-','LineWidth', 1.5,'color',[0.1 0.1 0.4]}, 0.2);
hold
shadedErrorBar(x,mean(gp.meanpitchpertresp_postcereb), (std(gp.meanpitchpertresp_postcereb)/sqrt(13)),{'-','LineWidth', 1.5,'color',[0.8 0 0.2]}, 0.2);
goodplot
axis([0 500 -20 35])
subplot(122)
x=1:size(nanmean(gp.meanpitchpertresp_precereb),2);
%shadedErrorBar(x,mean(gp.meanpitchpertresp_postcereb), (std(gp.meanpitchpertresp_postcereb)/sqrt(13)),{'-','LineWidth', 1.5,'color',[0.8 0 0.2]}, 0.2);
shadedErrorBar(x,mean(gp.meanpitchpertresp_presham), (std(gp.meanpitchpertresp_presham)/sqrt(13)),{'-','LineWidth', 1.5,'color',[0.1 0.1 0.4]}, 0.2);
hold
shadedErrorBar(x,mean(gp.meanpitchpertresp_postsham), (std(gp.meanpitchpertresp_presham)/sqrt(13)),{'-','LineWidth', 1.5,'color',[0.8 0 0.2]}, 0.2);
goodplot
axis([0 500 -20 35])
legend('moo', 'moo')
print(gcf, '-dpdf', '-r150', '/Users/zagnew/cerebellarTBS/data_analysis/figures/wholetrialcomp2.pdf');





% or 
figure
subplot(311)
plot(frame_taxis,nanmean(gp.meanpitchpertresp_precereb)','r','LineWidth',3);
hold 
plot(frame_taxis,nanmean(gp.meanpitchpertresp_postcereb)','Color',masked_colour,'LineWidth',3);
axis([-0.25 1 -20 30])
title('cerebellar stimulation')
goodplot

subplot(312)
plot(frame_taxis,nanmean(gp.meanpitchpertresp_presham)','r','LineWidth',3);
hold 
plot(frame_taxis,nanmean(gp.meanpitchpertresp_postsham)','Color',masked_colour,'LineWidth',3);
goodplot
title('vertex stimulation')
axis([-0.25 1 -20 30])

subplot(313)
plot(frame_taxis,nanmean(gp.meanpitchpertresp_prebeh)','r','LineWidth',3);
hold 
plot(frame_taxis,nanmean(gp.meanpitchpertresp_postbeh)','Color',masked_colour,'LineWidth',3);
goodplot
axis([-0.25 1 -20 30])
title('no stimulation')
legend('pre stim', 'post stim','Location','SouthEast')
print(gcf, '-dpdf', '-r150', '/Users/zagnew/cerebellarTBS/data_analysis/figures/preandpoststimresponses2.pdf');


%% calc and plot peak difference in compensation
gp.cereb_peakdiff_gpmean=mean(gp.cereb_peakdiff_eachsub)
gp.cereb_peakdiff_gpsem=std(gp.cereb_peakdiff_eachsub)/sqrt(numsubs_1);

gp.sham_peakdiff_gpmean=mean(gp.sham_peakdiff_eachsub)
gp.sham_peakdiff_gpsem=std(gp.sham_peakdiff_eachsub)/sqrt(numsubs_1);

gp.beh_peakdiff_gpmean=mean(gp.beh_peakdiff_eachsub)
gp.beh_peakdiff_gpsem=std(gp.beh_peakdiff_eachsub)/sqrt(numsubs_2);

figure
whitebg('white')
annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Peak Difference in Compensation', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')

y_peakpert=[gp.cereb_peakdiff_gpmean gp.sham_peakdiff_gpmean gp.beh_peakdiff_gpmean];
errY2 = [gp.cereb_peakdiff_gpsem gp.sham_peakdiff_gpsem gp.beh_peakdiff_gpsem];
h = barwitherr(errY2, y_peakpert);% Plot with errorbars

set(gca,'XTickLabel',{'cerebellar','sham', 'beh'})
ylabel('Diff in Peak Compensation to Perturbation (cents) post - pre')
set(h(1),'FaceColor',[1 1 1],'EdgeColor', masked_colour ,'LineWidth',3);
goodplot
print(gcf, '-dpdf', '-r150', '/Users/zagnew/cerebellarTBS/data_analysis/figures/MeanDiffPeakComp.pdf');

%% anova on peak diff data
peakdiffdata=[gp.cereb_peakdiff_eachsub gp.sham_peakdiff_eachsub gp.beh_peakdiff_eachsub]
alldata=length(peakdiffdata);

%create subject group
% test=ones(1,numsubs_1);
% test2=test*2;
% test3=test*3;
% test4=test*4;
% test5=test*5;
% test6=test*6;
% test7=test*7;
% test8=test*8;
% subjectgroup=[test test2 test3 test4 test5 test6 test7 test8 test8+1 test8+2 test8+3 test8+4 test8+5];
% subjectgroup=[subjectgroup subjectgroup];

%create window group
site = cell(1,alldata);
for i=1:numsubs_1
    site{i} = 'cereb';
end
for i=numsubs_1+1:numsubs_1*2
    site{i} = 'sham';
end
for i=(numsubs_1*2)+1:(numsubs_1*2)+numsubs_2
    site{i} = 'behonly';
end
site=site';

group1=site;
[pvals,tbl,stats] = anovan(peakdiffdata, {group1},'model',2, 'random',1,'varnames',{'stim site'});


%% find peaks
gp.gpmean_prebeh=nanmean(gp.meanpitchpertresp_prebeh);
gp.gpmean_postbeh=nanmean(gp.meanpitchpertresp_postbeh);
gp.gpmean_precereb=nanmean(gp.meanpitchpertresp_precereb);
gp.gpmean_postcereb=nanmean(gp.meanpitchpertresp_postcereb);
gp.gpmean_presham=nanmean(gp.meanpitchpertresp_presham);
gp.gpmean_postsham=nanmean(gp.meanpitchpertresp_postsham);

gp.gpmean_prebeh_peak=max(gp.gpmean_prebeh)
gp.gpmean_postbeh_peak=max(gp.gpmean_postbeh)
gp.gpmean_precereb_peak=max(gp.gpmean_precereb)
gp.gpmean_postcereb_peak=max(gp.gpmean_postcereb)
gp.gpmean_presham_peak=max(gp.gpmean_presham)
gp.gpmean_postsham_peak=max(gp.gpmean_postsham)

figure
whitebg('white')
annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'peak difference in compensation', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
y_peakpert=[gp.gpmean_precereb_peak gp.gpmean_postcereb_peak; gp.gpmean_presham_peak gp.gpmean_postsham_peak ; gp.gpmean_prebeh_peak gp.gpmean_postbeh_peak];
errY2 = [2 1; 2.4 4  ; 1.3 2];
h = barwitherr(errY2, y_peakpert);% Plot with errorbars
set(gca,'XTickLabel',{'cerebellar','sham', 'beh'})
ylabel('Peak Compensation to Perturbation (cents)')
set(h(1),'FaceColor',masked_colour,'EdgeColor', masked_colour ,'LineWidth',1.5);
set(h(2),'FaceColor',clear_colour,'EdgeColor', clear_colour ,'LineWidth',1.5);
goodplot
legend('pre stim', 'post stim')
print(gcf, '-dpdf', '-r150', '/Users/zagnew/cerebellarTBS/data_analysis/figures/MeanPeakComp.pdf');


gp.gpmean_precereb_peak
%% anova on whole trial data
anovadata=[...
    gp.cereb_diff_eachsub(1,:), ...
    gp.cereb_diff_eachsub(2,:), ...
    gp.cereb_diff_eachsub(3,:), ...
    gp.cereb_diff_eachsub(4,:), ...
    gp.cereb_diff_eachsub(5,:), ...
    gp.cereb_diff_eachsub(6,:), ...
    gp.cereb_diff_eachsub(7,:), ...
    gp.cereb_diff_eachsub(8,:), ...
    gp.cereb_diff_eachsub(9,:), ...
    gp.cereb_diff_eachsub(10,:), ...
    gp.cereb_diff_eachsub(11,:), ...
    gp.cereb_diff_eachsub(12,:), ...
    gp.cereb_diff_eachsub(13,:), ...   
    gp.sham_diff_eachsub(1,:), ...
    gp.sham_diff_eachsub(2,:), ...
    gp.sham_diff_eachsub(3,:), ...
    gp.sham_diff_eachsub(4,:), ...
    gp.sham_diff_eachsub(5,:), ...
    gp.sham_diff_eachsub(6,:), ...
    gp.sham_diff_eachsub(7,:), ...
    gp.sham_diff_eachsub(8,:), ...
    gp.sham_diff_eachsub(9,:), ...
    gp.sham_diff_eachsub(10,:), ...
    gp.sham_diff_eachsub(11,:), ...
    gp.sham_diff_eachsub(12,:), ...
    gp.sham_diff_eachsub(13,:)];

alldata=10738;
eachsub=413;

%create subject group
test=ones(1,eachsub);
test2=test*2;
test3=test*3;
test4=test*4;
test5=test*5;
test6=test*6;
test7=test*7;
test8=test*8;
subjectgroup=[test test2 test3 test4 test5 test6 test7 test8 test8+1 test8+2 test8+3 test8+4 test8+5];
subjectgroup=[subjectgroup subjectgroup];

%create window group
site = cell(1,alldata);
for i=1:alldata/2
    site{i} = 'cereb';
end
for i=alldata/2+1:alldata
    site{i} = 'sham';
end
site=site';

group1=subjectgroup;
group2=site;
[pvals,tbl,stats] = anovan(anovadata,{group1 group2},'model','interaction','varnames',{'subject' 'stim site'});

% significant effect of stimulation site on pre and post compensation *****
[pvals,tbl,stats] = anovan(anovadata, {group1 group2},'model',2, 'random',1,'varnames',{'subject' 'stim site'});

% plot this
figure
plot(gp.cereb_diff_eachsub_EW_mean,'k','LineWidth', 1.2)
hold 
plot(gp.cereb_diff_eachsub_LW_mean, 'k', 'LineWidth', 1.2)
plot(gp.sham_diff_eachsub_EW_mean, 'r','LineWidth', 1.2)
plot(gp.sham_diff_eachsub_LW_mean', 'r','LineWidth', 1.2)
% plot(gp.beh_diff_eachsub_EW_mean, 'g','LineWidth', 1.2)
% plot(gp.beh_diff_eachsub_LW_mean', 'g','LineWidth', 1.2)
goodplot

figure
whitebg('white')
annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'mean difference in compensation', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')

y_pertcomp=[nanmean(nanmean(gp.cereb_diff_eachsub))  nanmean(nanmean(gp.sham_diff_eachsub)) nanmean(nanmean(gp.beh_diff_eachsub))];
errY2 = [nanstd(nanstd(gp.cereb_diff_eachsub))/sqrt(numsubs_1) nanstd(nanstd(gp.sham_diff_eachsub))/sqrt(numsubs_1) nanstd(nanstd(gp.beh_diff_eachsub))/sqrt(numsubs_2)];
h = barwitherr(errY2, y_pertcomp);% Plot with errorbars
set(gca,'XTickLabel',{'cerebellar','sham', 'beh'})
ylabel('Difference in Compensation to Perturbation (cents)')
set(h(1),'FaceColor',[1 1 1],'EdgeColor', masked_colour ,'LineWidth',1.5);
goodplot

print(gcf, '-dpdf', '-r150', '/Users/zagnew/cerebellarTBS/data_analysis/figures/wholetrial_diffincomp.pdf');

% calculate sems
gp.cereb_diff_eachsub_sem=std(gp.cereb_diff_eachsub)/numsubs_1;
gp.sham_diff_eachsub_sem=std(gp.sham_diff_eachsub)/numsubs_1;
gp.beh_diff_eachsub_sem=std(gp.beh_diff_eachsub)/numsubs_2;

% -------------------------------------------------------------------------
%% early and late windows
EW_LW;








% %% early and late window anova including behavioural condition only
% anova_early_late= ...
%     [gp.cereb_diff_eachsub_EW(1,:) gp.cereb_diff_eachsub_EW(2,:) gp.cereb_diff_eachsub_EW(3,:) gp.cereb_diff_eachsub_EW(4,:) ...
%     gp.cereb_diff_eachsub_EW(5,:) gp.cereb_diff_eachsub_EW(6,:) gp.cereb_diff_eachsub_EW(7,:) gp.cereb_diff_eachsub_EW(8,:)...
%     gp.cereb_diff_eachsub_EW(9,:) gp.cereb_diff_eachsub_EW(10,:) gp.cereb_diff_eachsub_EW(11,:) gp.cereb_diff_eachsub_EW(12,:) ...
%     gp.cereb_diff_eachsub_EW(13,:) 
%     
%     gp.cereb_diff_eachsub_LW(1,:) gp.cereb_diff_eachsub_LW(2,:) gp.cereb_diff_eachsub_LW(3,:) gp.cereb_diff_eachsub_LW(4,:) ...
%     gp.cereb_diff_eachsub_LW(5,:) gp.cereb_diff_eachsub_LW(6,:) gp.cereb_diff_eachsub_LW(7,:) gp.cereb_diff_eachsub_LW(8,:)...
%     gp.cereb_diff_eachsub_LW(9,:) gp.cereb_diff_eachsub_LW(10,:) gp.cereb_diff_eachsub_LW(11,:) gp.cereb_diff_eachsub_LW(12,:) ...
%     gp.cereb_diff_eachsub_LW(13,:) ...
%     
%     gp.sham_diff_eachsub_EW(1,:) gp.sham_diff_eachsub_EW(2,:) gp.sham_diff_eachsub_EW(3,:) gp.sham_diff_eachsub_EW(4,:) ...
%     gp.sham_diff_eachsub_EW(5,:) gp.sham_diff_eachsub_EW(6,:) gp.sham_diff_eachsub_EW(7,:) gp.sham_diff_eachsub_EW(8,:)...
%     gp.sham_diff_eachsub_EW(9,:) gp.sham_diff_eachsub_EW(10,:) gp.sham_diff_eachsub_EW(11,:) gp.sham_diff_eachsub_EW(12,:) ...
%     gp.sham_diff_eachsub_EW(13,:) 
%     
%     gp.sham_diff_eachsub_LW(1,:) gp.sham_diff_eachsub_LW(2,:) gp.sham_diff_eachsub_LW(3,:) gp.sham_diff_eachsub_LW(4,:) ...
%     gp.sham_diff_eachsub_LW(5,:) gp.sham_diff_eachsub_LW(6,:) gp.sham_diff_eachsub_LW(7,:) gp.sham_diff_eachsub_LW(8,:)...
%     gp.sham_diff_eachsub_LW(9,:) gp.sham_diff_eachsub_LW(10,:) gp.sham_diff_eachsub_LW(11,:) gp.sham_diff_eachsub_LW(12,:) ...
%     gp.sham_diff_eachsub_LW(13,:) 
%     
%     gp.beh_diff_eachsub_EW(1,:) gp.beh_diff_eachsub_EW(2,:) gp.beh_diff_eachsub_EW(3,:) gp.beh_diff_eachsub_EW(4,:) ...
%     gp.beh_diff_eachsub_EW(5,:) gp.beh_diff_eachsub_EW(6,:) gp.beh_diff_eachsub_EW(7,:) gp.beh_diff_eachsub_EW(8,:)...
%     gp.beh_diff_eachsub_EW(9,:) 
%     
%     gp.beh_diff_eachsub_LW(1,:) gp.beh_diff_eachsub_LW(2,:) gp.beh_diff_eachsub_LW(3,:) gp.beh_diff_eachsub_LW(4,:) ...
%     gp.beh_diff_eachsub_LW(5,:) gp.beh_diff_eachsub_LW(6,:) gp.beh_diff_eachsub_LW(7,:) gp.beh_diff_eachsub_LW(8,:)...
%     gp.beh_diff_eachsub_LW(9,:)];
% 
% alldata=length(anova_early_late);
% eachsub=151;
% 
% %create subject group
% test=ones(1,eachsub);
% test2=test*2;
% test3=test*3;
% test4=test*4;
% test5=test*5;
% test6=test*6;
% test7=test*7;
% test8=test*8;
% test9=test*9;
% test10=test*10;
% test11=test*11;
% test12=test*12;
% test13=test*13;
% 
% subjectgroup=[test test2 test3 test4 test5 test6 test7 test8 test8+1 test8+2 test8+3 test8+4 test8+5];
% subjectgroup2=[test test2 test3 test4 test5 test6 test7 test8 test8+1];
% subjectgroup=[subjectgroup subjectgroup subjectgroup subjectgroup subjectgroup2 subjectgroup2];
%  
% window_early=ones(1,2678);
% window_late=window_early*2;
% 
% window_early2=ones(1,1854);
% window_late2=window_early2*2;
% window=[window_early window_late window_early window_late window_early2 window_late2];
% 
% cereb=ones(1,206*26);
% sham=cereb*2;
% beh=ones(1,3708)*3;
% stim_site=[cereb sham beh];
% group1=subjectgroup;
% group2=window;
% group3=stim_site;
% 
% % p = anovan(anova_early_late,{group1 group2 group3},'model','interaction');
% % [pvals,tbl,stats]  = anovan(anova_early_late,{group1 group2 group3}, 'full');
% 
% %[pvals,tbl,stats] = anovan(anova_early_late, {group1 group2 group3},'model',2, 'random',1,'varnames',{'subject' 'window' 'stim site'});
% % multcompare(stats)


%% whole trial anova including two stim conditions only
anovadata=[...
    gp.cereb_diff_eachsub(1,:), ...
    gp.cereb_diff_eachsub(2,:), ...
    gp.cereb_diff_eachsub(3,:), ...
    gp.cereb_diff_eachsub(4,:), ...
    gp.cereb_diff_eachsub(5,:), ...
    gp.cereb_diff_eachsub(6,:), ...
    gp.cereb_diff_eachsub(7,:), ...
    gp.cereb_diff_eachsub(8,:), ...
    gp.cereb_diff_eachsub(9,:), ...
    gp.cereb_diff_eachsub(10,:), ...
    gp.cereb_diff_eachsub(11,:), ...
    gp.cereb_diff_eachsub(12,:), ...
    gp.cereb_diff_eachsub(13,:), ...   
    gp.sham_diff_eachsub(1,:), ...
    gp.sham_diff_eachsub(2,:), ...
    gp.sham_diff_eachsub(3,:), ...
    gp.sham_diff_eachsub(4,:), ...
    gp.sham_diff_eachsub(5,:), ...
    gp.sham_diff_eachsub(6,:), ...
    gp.sham_diff_eachsub(7,:), ...
    gp.sham_diff_eachsub(8,:), ...
    gp.sham_diff_eachsub(9,:), ...
    gp.sham_diff_eachsub(10,:), ...
    gp.sham_diff_eachsub(11,:), ...
    gp.sham_diff_eachsub(12,:), ...
    gp.sham_diff_eachsub(13,:)];

alldata=10738;
eachsub=413;

%create subject group
test=ones(1,eachsub);
test2=test*2;
test3=test*3;
test4=test*4;
test5=test*5;
test6=test*6;
test7=test*7;
test8=test*8;
subjectgroup=[test test2 test3 test4 test5 test6 test7 test8 test8+1 test8+2 test8+3 test8+4 test8+5];
subjectgroup=[subjectgroup subjectgroup];

%create window group
site = cell(1,alldata);
for i=1:alldata/2
    site{i} = 'cereb';
end
for i=alldata/2+1:alldata
    site{i} = 'sham';
end
site=site';

group1=subjectgroup;
group2=site;
group3=window
%p = anovan(anovadata,{group1 group2},'model','interaction'
[pvals,tbl,stats] = anovan(anovadata, {group1 group2},'model',1, 'random',1,'varnames',{'subject' 'stim site'});
[pvals,tbl,stats] = anovan(anovadata, {group1 group2},'model',2, 'random',1,'varnames',{'subject' 'stim site'});
[pvals,tbl,stats] = anovan(anovadata, {group1 group2},'model',3,'varnames',{'subject' 'stim site'});

% anovan(mileage, {factory carmod}, 'model',2, 'random',1,'varnames',{'Factory' 'Car Model'});

save GROUPDATA

% do running ttest across all difference trials. 
ttest_compdiff=[]
for iframe=1:413
    ttest_compdiff(:, iframe)=ttest(gp.cereb_diff_eachsub(:,iframe), gp.sham_diff_eachsub(:,iframe));            
end

figure
plot(frame_taxis(idxes4anal),ttest_compdiff, 'k', 'LineWidth',1.3);
hold on
plot(frame_taxis(idxes4anal),ttest(gp.cereb_diff_eachsub, gp.sham_diff_eachsub), 'r','LineWidth',1.3);

%% early and late windows

gp.cereb_diff_eachsub_EW=gp.cereb_diff_eachsub(:,EW);
gp.cereb_diff_eachsub_EW=gp.cereb_diff_eachsub(:,LW);
gp.cereb_diff_eachsub_EW_mean=mean(gp.cereb_diff_eachsub(:,EW));
gp.cereb_diff_eachsub_EW_mean=mean(gp.cereb_diff_eachsub(:,LW));

gp.sham_diff_eachsub_EW=gp.sham_diff_eachsub(:,EW);
gp.sham_diff_eachsub_EW=gp.sham_diff_eachsub(:,LW);
gp.sham_diff_eachsub_EW_mean=mean(gp.sham_diff_eachsub(:,EW));
gp.sham_diff_eachsub_EW_mean=mean(gp.sham_diff_eachsub(:,LW));


%check
subplot(311); plot(gp.cereb_diff_eachsub')
hold
plot(mean(gp.cereb_diff_eachsub)', 'k', 'LineWidth', 1.3)
plot((mean(gp.cereb_diff_eachsub))+gp.cereb_diff_eachsub_sem, 'LineWidth', 1.2)
plot((mean(gp.cereb_diff_eachsub))-gp.cereb_diff_eachsub_sem, 'LineWidth', 1.2)
axis([ 0 450 -25 25])

subplot(312); plot(gp.cereb_diff_eachsub(:,EW)')
hold 
plot(mean(gp.cereb_diff_eachsub(:,EW))', 'k', 'LineWidth', 1.3)
axis([ 0 450 -25 25])

subplot(313); plot(gp.cereb_diff_eachsub(:,LW)')
hold
plot(mean(gp.cereb_diff_eachsub(:,LW))', 'k', 'LineWidth', 1.3)
axis([ 0 450 -25 25])

