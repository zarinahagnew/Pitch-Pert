%% ANOVA for group data
% ZKA Oct 2015

% once you have run get_data_cTBS.m, this will create 4 files that are used here. This script analyses the lot and saves it as GROUPDATA.mat. 


%% the analysis:
% cents4comp(1).pitch_in.dat{1} is the response to the postive shift and looks like a negative vocal response
% here we flip this and add it to the negative pert, so that it all looks like we are seeing responses to the negative perturbation. This looks like an increase in pitch. 
% the flipped data is gp.precereb.patient_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat


% pitch_in.dat{1} - response to down pert
% pitch_in.dat{2} response to up pert
% pitch_in.dat{3} respone to all but *it's not flipped*
% to flip it use: centsdev_dat.(group{igroup}).subj(isubj).absdat = {-the_dat{1} the_dat{2} [-the_dat{1}; the_dat{2}]};

%% interpretation
% all responses are flipped so that they look like responses to the down pert
% diff is post - pre, SO a positive value indiciates increased pitch in the
% post condition compared to the pre stim condition, and a negative value indicates the
% opposite

% TO DO
%for the anova want to run for time point also, and put frame in as 413 levels

% or 
% use mean .peakcomp .comp values post and pre and put into anova
% to see wherepitchpert_fourway this effect is in the trial look at mean(subdata_pre_cereb(isub).pert_resp.tpeak{1})
% /data/bil-mb4/zarinah-data/cerebellar-data/pitch-pert-cTBS/cTBS_data/data_analysis/cerebellarTBS/groupdata
% Also can do an anova on pert_resp.tpeak{1}) to see if the time point of
% the peak comp is different. 
 
% in general, *comp* in name means that the cents response, divded by the cent of the pert and made positive. 

% So look at the peak comp and t peak to get at both amp of comp and time
% peak of comp, and you can check that this fits with this once you have 
% concatednatepitchpert_fourwayd/added pitch_in{1} and {2} of patient_dat.pert_resp(isubj).cents4comp(1).pitch_in
% flip the first one (all negs become pos and all popatient_dat.pert_resps becomes neg)

% = {-the_dat{1} the_dat{2} [-the_dat{1}; the_dat{2}]};

% pert_resp.cents4pert_resp.cents4comp.pitch_in.meancomp.pitch_in.mean
% cd /data/bil-mb4/zarinah-data/cerebellar-data/pitch-pert-cTBS/cTBS_data/data_analysis

%% having run /home/zagnew/data_analysis_code/pitch_pert_stats/get_data_cTBS.m

set(0,'DefaultFigureWindowStyle','docked');
clear all
close all
z_colours;

numsubs_1=13;
numsubs_2=9;

%cd /Users/zagnew/Desktop/cTBSdata
cd /Users/zagnew/cerebellarTBS/raw_data/
load /Users/zagnew/cerebellarTBS/data_analysis/groupdata/GROUPDATA.mat

% for plottin' 
frame_taxis = gp.precereb.patient_dat.frame_taxis;
tlims4anal_req = [-0.2 1.0];
ilims4anal = dsearchn(frame_taxis',tlims4anal_req')';
idxes4anal = ilims4anal(1):ilims4anal(2);
tlims4anal = frame_taxis(ilims4anal);

%% flip pitch and concatenate pitch_in.dat{1} and {2}
for isubj=1:numsubs_1 
    %cerebellar    
    gp.precereb.patient_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat = [-gp.precereb.patient_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{1} ; gp.precereb.patient_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{2}];    
    gp.postcereb.control_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat = [-gp.postcereb.control_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{1} ; gp.postcereb.control_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{2}];   
    
    %vertex
    gp.presham.patient_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat = [-gp.presham.patient_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{1} ; gp.presham.patient_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{2}];
    gp.postsham.control_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat = [-gp.postsham.control_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{1} ; gp.postsham.control_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{2}];
end
for isubj=1:numsubs_2
 % beh only
    gp.prebeh.patient_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat = [-gp.prebeh.patient_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{1} ; gp.prebeh.patient_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{2}];
    gp.postbeh.control_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat = [-gp.postbeh.control_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{1}; gp.postbeh.control_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{2}];
end

% 
% % check
% subplot(311)
% plot(mean(gp.precereb.patient_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{1}))
% axis([0 600 -60 60])
% 
% subplot(312)
% plot(mean(gp.precereb.patient_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{2}))
% axis([0 600 -60 60])
% 
% subplot(313)
% plot(mean(gp.precereb.patient_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat))
% hold
% plot(mean(gp.postcereb.control_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat), 'r')
% axis([0 600 -60 60])



%% calc mean comp
for isubj=1:numsubs_1;
    gp.meanpitchpertresp_precereb(isubj,:)=mean(gp.precereb.patient_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat);
    gp.meanpitchpertresp_postcereb(isubj,:)=mean(gp.postcereb.control_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat);
    gp.meanpitchpertresp_presham(isubj,:)=mean(gp.presham.patient_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat);
    gp.meanpitchpertresp_postsham(isubj,:)=mean(gp.postsham.control_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat);

    gp.cereb_diff_eachsub(isubj,:)=gp.meanpitchpertresp_postcereb(isubj,:)-gp.meanpitchpertresp_precereb(isubj,:);
    gp.sham_diff_eachsub(isubj,:)=gp.meanpitchpertresp_postsham(isubj,:)-gp.meanpitchpertresp_presham(isubj,:);    
end

for isubj=1:numsubs_2;
    gp.meanpitchpertresp_prebeh(isubj,:)=mean(gp.prebeh.patient_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat);
    gp.meanpitchpertresp_postbeh(isubj,:)=mean(gp.postbeh.control_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat);  
    gp.beh_diff_eachsub(isubj,:)=gp.meanpitchpertresp_postbeh(isubj,:)-gp.meanpitchpertresp_prebeh(isubj,:);
end

%% plot all post stim trials
figure
subplot(211)
plot(frame_taxis,nanmean(gp.meanpitchpertresp_precereb)','r','LineWidth',3);
hold 
plot(frame_taxis,nanmean(gp.meanpitchpertresp_presham)','Color',masked_colour,'LineWidth',3);
plot(frame_taxis,nanmean(gp.meanpitchpertresp_prebeh)','Color',clear_colour,'LineWidth',3);
goodplot
legend('cerebellar', 'vertex', 'no stim', 'Location','NorthWest')
subplot(212)
plot(frame_taxis,nanmean(gp.meanpitchpertresp_postcereb)','r','LineWidth',3);
hold 
plot(frame_taxis,nanmean(gp.meanpitchpertresp_postsham)','Color',masked_colour,'LineWidth',3);
plot(frame_taxis,nanmean(gp.meanpitchpertresp_postbeh)','Color',clear_colour,'LineWidth',3);
goodplot

print(gcf, '-dpdf', '-r150', '/Users/zagnew/cerebellarTBS/data_analysis/figures/preandpoststim responses.pdf');

% or 
figure
subplot(311)
plot(frame_taxis,nanmean(gp.meanpitchpertresp_precereb)','r','LineWidth',3);
hold 
plot(frame_taxis,nanmean(gp.meanpitchpertresp_postcereb)','Color',masked_colour,'LineWidth',3);
axis([-0.5 1 -20 30])
title('cerebellar stimulation')
goodplot

subplot(312)
plot(frame_taxis,nanmean(gp.meanpitchpertresp_presham)','r','LineWidth',3);
hold 
plot(frame_taxis,nanmean(gp.meanpitchpertresp_postsham)','Color',masked_colour,'LineWidth',3);
goodplot
title('vertex stimulation')

axis([-0.5 1 -20 30])

subplot(313)
plot(frame_taxis,nanmean(gp.meanpitchpertresp_prebeh)','r','LineWidth',3);
hold 
plot(frame_taxis,nanmean(gp.meanpitchpertresp_postbeh)','Color',masked_colour,'LineWidth',3);
goodplot
axis([-0.5 1 -20 30])
title('no stimulation')
legend('pre stim', 'post stim','Location','SouthWest')


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
errY2 = [5 5; 5 5  ; 5 5];
h = barwitherr(errY2, y_peakpert);% Plot with errorbars
set(gca,'XTickLabel',{'cerebellar','sham', 'beh'})
ylabel('Peak Compensation to Perturbation (cents)')
set(h(1),'FaceColor',masked_colour,'EdgeColor', masked_colour ,'LineWidth',1.5);
set(h(2),'FaceColor',clear_colour,'EdgeColor', clear_colour ,'LineWidth',1.5);
goodplot
legend('pre stim', 'post stim')
print(gcf, '-dpdf', '-r150', '/Users/zagnew/cerebellarTBS/data_analysis/figures/MeanPeakComp.pdf');



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
[pvals,tbl,stats] = anovan(anovadata, {group1 group2},'model',2, 'random',1,'varnames',{'subject' 'stim site'});

% % plot this
% figure
% plot(nanmean(nanmean(gp.cereb_diff_eachsub)),'k','LineWidth', 1.2)
% hold 
% plot(gp.cereb_diff_eachsub_LW_mean, 'k', 'LineWidth', 1.2)
% plot(gp.sham_diff_eachsub_EW_mean, 'r','LineWidth', 1.2)
% plot(gp.sham_diff_eachsub_LW_mean', 'r','LineWidth', 1.2)
% plot(gp.beh_diff_eachsub_EW_mean, 'g','LineWidth', 1.2)
% plot(gp.beh_diff_eachsub_LW_mean', 'g','LineWidth', 1.2)
% goodplot

figure
whitebg('white')
annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'mean difference in compensation, early and late windows', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
y_pertcomp=[nanmean(nanmean(gp.cereb_diff_eachsub))  nanmean(nanmean(gp.sham_diff_eachsub)) nanmean(nanmean(gp.beh_diff_eachsub))];
errY2 = [nanstd(nanstd(gp.cereb_diff_eachsub))/sqrt(numsubs_1) nanstd(nanstd(gp.sham_diff_eachsub))/sqrt(numsubs_1) nanstd(nanstd(gp.beh_diff_eachsub))/sqrt(numsubs_2)];
h = barwitherr(errY2, y_pertcomp);% Plot with errorbars
set(gca,'XTickLabel',{'cerebellar','sham', 'beh'})
ylabel('Difference in Compensation to Perturbation (cents)')
set(h(1),'FaceColor',masked_colour,'EdgeColor', masked_colour ,'LineWidth',1.5);
goodplot

print(gcf, '-dpdf', '-r150', '/Users/zagnew/cerebellarTBS/data_analysis/figures/wholetrial_diffincomp.pdf');

% calculate sems
gp.cereb_diff_eachsub_sem=std(gp.cereb_diff_eachsub)/numsubs_1;
gp.sham_diff_eachsub_sem=std(gp.sham_diff_eachsub)/numsubs_1;
gp.beh_diff_eachsub_sem=std(gp.beh_diff_eachsub)/numsubs_2;

% -------------------------------------------------------------------------
%% early and late windows
EW=100:250;
LW=250:400;

gp.cereb_diff_eachsub_EW=gp.cereb_diff_eachsub(:,EW);
gp.cereb_diff_eachsub_LW=gp.cereb_diff_eachsub(:,LW);
gp.cereb_diff_eachsub_EW_mean=mean(gp.cereb_diff_eachsub(:,EW));
gp.cereb_diff_eachsub_LW_mean=mean(gp.cereb_diff_eachsub(:,LW));

gp.sham_diff_eachsub_EW=gp.sham_diff_eachsub(:,EW);
gp.sham_diff_eachsub_EW=gp.sham_diff_eachsub(:,LW);
gp.sham_diff_eachsub_EW_mean=mean(gp.sham_diff_eachsub(:,EW));
gp.sham_diff_eachsub_EW_mean=mean(gp.sham_diff_eachsub(:,LW));

gp.beh_diff_eachsub_EW=gp.beh_diff_eachsub(:,EW);
gp.beh_diff_eachsub_EW=gp.beh_diff_eachsub(:,LW);
gp.beh_diff_eachsub_EW_mean=mean(gp.beh_diff_eachsub(:,EW));
gp.beh_diff_eachsub_EW_mean=mean(gp.beh_diff_eachsub(:,LW));


%% plot all EW and EW on top of each other
figure
subplot(311)
plot(gp.cereb_diff_eachsub_EW_mean,'k','LineWidth', 3)
hold 
plot(gp.cereb_diff_eachsub_LW_mean, 'k', 'LineWidth', 1)
goodplot
axis([ 0 151 -10 10])
title('cerebellar')

subplot(312)
plot(gp.sham_diff_eachsub_EW_mean,'k','LineWidth', 3)
hold
plot(gp.sham_diff_eachsub_LW_mean','k','LineWidth', 1)
axis([ 0 151 -10 10])
title('vertex')
goodplot

subplot(313)
plot(gp.beh_diff_eachsub_EW_mean, 'k','LineWidth', 3)
hold
plot(gp.beh_diff_eachsub_LW_mean', 'k','LineWidth', 1)
axis([ 0 151 -10 10])
goodplot
title('no stimulation')

print(gcf, '-dpdf', '-r150', '/Users/zagnew/cerebellarTBS/data_analysis/figures/EW_LW_diffincomp.pdf');


figure
whitebg('white')
annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'mean difference in compensation, early and late windows', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
EW_LW_pertcomp=[mean(gp.cereb_diff_eachsub_EW_mean) mean(gp.cereb_diff_eachsub_LW_mean) ; mean(gp.sham_diff_eachsub_EW_mean) mean(gp.sham_diff_eachsub_LW_mean) ; mean(gp.beh_diff_eachsub_EW_mean) mean(gp.beh_diff_eachsub_LW_mean)];
EW_LWerr = [std(gp.cereb_diff_eachsub_EW_mean)/sqrt(length(EW)) std(gp.cereb_diff_eachsub_LW_mean)/sqrt(length(EW)) ; std(gp.sham_diff_eachsub_EW_mean)/sqrt(length(EW)) std(gp.sham_diff_eachsub_LW_mean)/sqrt(length(EW)) ; std(gp.beh_diff_eachsub_EW_mean)/sqrt(length(EW)) std(gp.beh_diff_eachsub_LW_mean)/sqrt(length(EW))];
h = barwitherr(EW_LWerr, EW_LW_pertcomp);% Plot with errorbars

set(gca,'XTickLabel',{'cerebellar','certex', 'no stimulation'})
ylabel('Difference in Compensation to Perturbation (cents)')
set(h(1),'FaceColor',clear_colour,'EdgeColor', clear_colour ,'LineWidth',1.5);
set(h(2),'FaceColor',masked_colour,'EdgeColor', masked_colour ,'LineWidth',1.5);
goodplot
legend('early window', 'late window','Location','SouthEast')
print(gcf, '-dpdf', '-r150', '/Users/zagnew/cerebellarTBS/data_analysis/figures/early_late_winds_diff_comp.pdf');


%% early and late window anova including behavioural condition only
anova_early_late= ...
    [gp.cereb_diff_eachsub_EW(1,:) gp.cereb_diff_eachsub_EW(2,:) gp.cereb_diff_eachsub_EW(3,:) gp.cereb_diff_eachsub_EW(4,:) ...
    gp.cereb_diff_eachsub_EW(5,:) gp.cereb_diff_eachsub_EW(6,:) gp.cereb_diff_eachsub_EW(7,:) gp.cereb_diff_eachsub_EW(8,:)...
    gp.cereb_diff_eachsub_EW(9,:) gp.cereb_diff_eachsub_EW(10,:) gp.cereb_diff_eachsub_EW(11,:) gp.cereb_diff_eachsub_EW(12,:) ...
    gp.cereb_diff_eachsub_EW(13,:) gp.cereb_diff_eachsub_EW(1,:) gp.cereb_diff_eachsub_EW(2,:) gp.cereb_diff_eachsub_EW(3,:) gp.cereb_diff_eachsub_EW(4,:) ...
    gp.cereb_diff_eachsub_EW(5,:) gp.cereb_diff_eachsub_EW(6,:) gp.cereb_diff_eachsub_EW(7,:) gp.cereb_diff_eachsub_EW(8,:)...
    gp.cereb_diff_eachsub_EW(9,:) gp.cereb_diff_eachsub_EW(10,:) gp.cereb_diff_eachsub_EW(11,:) gp.cereb_diff_eachsub_EW(12,:) ...
    gp.cereb_diff_eachsub_EW(13,:) ...
    gp.sham_diff_eachsub_EW(1,:) gp.sham_diff_eachsub_EW(2,:) gp.sham_diff_eachsub_EW(3,:) gp.sham_diff_eachsub_EW(4,:) ...
    gp.sham_diff_eachsub_EW(5,:) gp.sham_diff_eachsub_EW(6,:) gp.sham_diff_eachsub_EW(7,:) gp.sham_diff_eachsub_EW(8,:)...
    gp.sham_diff_eachsub_EW(9,:) gp.sham_diff_eachsub_EW(10,:) gp.sham_diff_eachsub_EW(11,:) gp.sham_diff_eachsub_EW(12,:) ...
    gp.sham_diff_eachsub_EW(13,:) gp.sham_diff_eachsub_EW(1,:) gp.sham_diff_eachsub_EW(2,:) gp.sham_diff_eachsub_EW(3,:) gp.sham_diff_eachsub_EW(4,:) ...
    gp.sham_diff_eachsub_EW(5,:) gp.sham_diff_eachsub_EW(6,:) gp.sham_diff_eachsub_EW(7,:) gp.sham_diff_eachsub_EW(8,:)...
    gp.sham_diff_eachsub_EW(9,:) gp.sham_diff_eachsub_EW(10,:) gp.sham_diff_eachsub_EW(11,:) gp.sham_diff_eachsub_EW(12,:) ...
    gp.sham_diff_eachsub_EW(13,:) gp.beh_diff_eachsub_EW(1,:) gp.beh_diff_eachsub_EW(2,:) gp.beh_diff_eachsub_EW(3,:) gp.beh_diff_eachsub_EW(4,:) ...
    gp.beh_diff_eachsub_EW(5,:) gp.beh_diff_eachsub_EW(6,:) gp.beh_diff_eachsub_EW(7,:) gp.beh_diff_eachsub_EW(8,:)...
    gp.beh_diff_eachsub_EW(9,:) gp.beh_diff_eachsub_EW(1,:) gp.beh_diff_eachsub_EW(2,:) gp.beh_diff_eachsub_EW(3,:) gp.beh_diff_eachsub_EW(4,:) ...
    gp.beh_diff_eachsub_EW(5,:) gp.beh_diff_eachsub_EW(6,:) gp.beh_diff_eachsub_EW(7,:) gp.beh_diff_eachsub_EW(8,:)...
    gp.beh_diff_eachsub_EW(9,:)];



alldata=length(anova_early_late);
eachsub=206;
%create subject group
test=ones(1,eachsub);
test2=test*2;
test3=test*3;
test4=test*4;
test5=test*5;
test6=test*6;
test7=test*7;
test8=test*8;
test9=test*9;
test10=test*10;
test11=test*11;
test12=test*12;
test13=test*13;

subjectgroup=[test test2 test3 test4 test5 test6 test7 test8 test8+1 test8+2 test8+3 test8+4 test8+5];
subjectgroup2=[test test2 test3 test4 test5 test6 test7 test8 test8+1];
subjectgroup=[subjectgroup subjectgroup subjectgroup subjectgroup subjectgroup2 subjectgroup2];
 
window_early=ones(1,2678);
window_late=window_early*2;

window_early2=ones(1,1854);
window_late2=window_early2*2;
window=[window_early window_late window_early window_late window_early2 window_late2];

cereb=ones(1,206*26);
sham=cereb*2;
beh=ones(1,3708)*3;
stim_site=[cereb sham beh];
group1=subjectgroup;
group2=window;
group3=stim_site;

p = anovan(anova_early_late,{group1 group2 group3},'model','interaction');
[pvals,tbl,stats]  = anovan(anova_early_late,{group1 group2 group3}, 'full');

%[pvals,tbl,stats] = anovan(anova_early_late, {group1 group2 group3},'model',2, 'random',1,'varnames',{'subject' 'window' 'stim site'});
multcompare(stats)

save GROUPDATA


%% early and late window anova including two stim conditions only
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
p = anovan(anovadata,{group1 group2},'model','interaction');
% p = anovan(anovadata,{group1 group2 }, 'full')
[pvals,tbl,stats] = anovan(anovadata, {group1 group2},'model',2, 'random',1,'varnames',{'subject' 'stim site'});
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


% xlabel('condition (1:10)')
% ylabel('WT stdev pitch')
% title(sprintf('HC pre-WTstdev pitch across all conditions separately'));
% 
% goodplot
% 
% subplot(212)
% for isubj=1:length(gp.presham.control_dat.pert_resp)
% plot(gp.sham_diff_eachsub(isubj,:))
% hold on
% plot(mean(gp.sham_diff_eachsub), 'm')
% end
% goodplot
% 



% isubj=1;
% figure
% subplot(311)
% for moo=1:26
%     plot(gp.precereb.patient_dat.pert_resp(isubj,1).cents4comp(1).pitch_in.dat{1}(moo,:))
%     hold on
%     plot(mean(gp.precereb.patient_dat.pert_resp(isubj,1).cents4comp(1).pitch_in.dat{1}), 'm')
% end
% 
% subplot(312)
% for moo=1:14
%     plot(gp.precereb.patient_dat.pert_resp(isubj,1).cents4comp(1).pitch_in.dat{2}(moo,:))
%     hold on
%     plot(mean(gp.precereb.patient_dat.pert_resp(isubj,1).cents4comp(1).pitch_in.dat{2}), 'm')
% end
% 
% subplot(313)
% for moo=1:40
%     plot(gp.precereb.patient_dat.pert_resp(isubj,1).cents4comp(1).pitch_in.dat{3}(moo,:))
%     hold on
%     plot(mean(gp.precereb.patient_dat.pert_resp(isubj,1).cents4comp(1).pitch_in.dat{3}), 'm')
% end

%% anova on mean data by time window


% plot this stuffplot(nanmean(gp.postcereb.control_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat)) 
% for i=1:numsubs_1; 
%     plot(gp.precereb.patient_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{1}(i,:))
% hold on
% end
% 
% for i=1:numsubs_1
% plot(gp.precereb.patient_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{2}(i,:), 'm')
% hold on
% end
% title(sprintf('cents4comp(1).pitch_in.dat{1 and 2}'));
% 
% % sham
% figure
% subplot(211)
% plot(nanmean(gp.precereb.patient_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{1}), 'LineWidth', 3);
% hold
% plot(nanmean(gp.precereb.patient_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{2}), 'm', 'LineWidth', 3);
% goodplot
% xlabel('frames')
% title(sprintf('pre cereb - cents4comp(1).pitch_in.dat{1 and 2}'));
% 
% figure
% subplot(311)
% plot(nanmean(gp.precereb.patient_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat),'LineWidth', 3);
% hold
% plot(nanmean(gp.postcereb.control_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat), 'LineWidth', 3);
% 
% subplot(212)
% plot(nanmean(-gp.precereb.patient_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{1}),'LineWidth', 3);
% hold on
% plot(nanmean(gp.precereb.patient_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{2}), 'm', 'LineWidth', 3);
% goodplot
% xlabel('frames')
% title(sprintf('pre cereb - cents4comp(1).pitch_in.dat{1 and 2} with 1 flipped'));
% plot(nanmean(gp.precereb.patient_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat),'c', 'LineWidth', 3)


%

% 
% % plot individual subjects
% figure
% for isubj=1:numsubs_2
%     subplot(numsubs_1/3,numsubs_2/3, isubj)
%     h(1)=plot(frame_taxis(idxes4anal), gp.cereb_diff_eachsub(isubj,:), 'k','LineWidth',3);
%     hold on
%     h(2)=plot(frame_taxis(idxes4anal),gp.sham_diff_eachsub(isubj,:), 'm','LineWidth',3);
%     h(3)=plot(frame_taxis(idxes4anal),gp.beh_diff_eachsub(isubj,:), 'c','LineWidth',3);
% 
%     xlabel('Time (s)')
%     axis([-0.2 1 -100 50])
%     goodplot
% end
% legend('cerebellar stimulation','sham stimulation', 'no stimulation', 'Location','southeast')
% saveas(gcf, sprintf('eachsubject.jpg'))    
% 
% 
%     
% %% plot pre and post comp for each subject
% for isubj=1:numsubs_2% gp.cereb_diff_eachsub_semlength(gp.precereb.patient_dat.pert_resp)
%     figure
%     subplot (3, 1, 1)
%     plot(frame_taxis(idxes4anal),gp.meanpitchpertresp_precereb(isubj,:), 'k','LineWidth',3);
%     hold on
%     plot(frame_taxis(idxes4anal),gp.meanpitchpertresp_postcereb(isubj,:), 'R','LineWidth',3);
%     plot(frame_taxis,ones(1,413),'c','LineWidth',3);
%     axis([-0.2 1 -60 60])
%     title(sprintf('cereb stimulation'));
%     xlabel('Time(s)')
%     ylabel('Mean Compensation (cents)')
%     goodplot
%     
%     subplot (3, 1, 2)
%     plot(frame_taxis(idxes4anal),gp.meanpitchpertresp_presham(isubj,:), 'k','LineWidth',3);
%     hold on
%     plot(frame_taxis(idxes4anal),gp.meanpitchpertresp_postsham(isubj,:), 'R','LineWidth',3);     plot(frame_taxis,ones(1,413),'c','LineWidth',3);
% 
%     axis([-0.2 1 -60 60])
%     title(sprintf('vertex stimulation'));
%     goodplot
%     xlabel('Time(s)')
%     ylabel('Mean Compensation (cents)')
%     saveas(gcf, sprintf('figure%d.jpg', isubj))
%     
%     subplot (3, 1, 3)
%     plot(frame_taxis(idxes4anal),gp.meanpitchpertresp_prebeh(isubj,:), 'k','LineWidth',3);
%     hold on
%     plot(frame_taxis(idxes4anal),gp.meanpitchpertresp_postbeh(isubj,:), 'R','LineWidth',3);     plot(frame_taxis,ones(1,413),'c','LineWidth',3);
% 
%     axis([-0.2 1 -60 60])
%     title(sprintf('no stimulation'));
%     goodplot
%     xlabel('Time(s)')
%     ylabel('Mean Compensation (cents)')
%     saveas(gcf, sprintf('figure%d.jpg', isubj))    
%     
% end

% % plot
% ymax=100;
% ymin=-150;
% figure
% subplot(311)
% for isubj=1:numsubs_1%length(gp.precereb.patient_dat.pert_resp)
%     plot(frame_taxis(idxes4anal), gp.cereb_diff_eachsub(isubj,:), 'k')
%     hold on
%     plot(frame_taxis(idxes4anal),mean(gp.cereb_diff_eachsub), 'r','LineWidth',3);
% end
% title(sprintf('cerebellar stimulation'));
% xlabel('frames')
% ylabel('post stim comp - pre stim comp')
% axis([-0.2 1 ymin ymax])
% goodplot
% 
% subplot(312)
% for isubj=1:numsubs_1 %length(gp.presham.control_dat.pert_resp)
%     plot(frame_taxis(idxes4anal), gp.sham_diff_eachsub(isubj,:), 'k')
%     hold on
%     plot(frame_taxis(idxes4anal), mean(gp.sham_diff_eachsub), 'r','LineWidth',3);
% end
% axis([-0.2 1 ymin ymax])
% xlabel('frames')
% ylabel('post stim comp - pre stim comp')
% title(sprintf('vertex stimulation'));
% goodplot
% 
% subplot(313)
% for isubj=1:numsubs_2 %length(gp.presham.control_dat.pert_resp)
%     plot(frame_taxis(idxes4anal), gp.beh_diff_eachsub(isubj,:), 'k')
%     hold on
%     plot(frame_taxis(idxes4anal), mean(gp.beh_diff_eachsub), 'r','LineWidth',3);
% end
% axis([-0.2 1 ymin ymax])
% xlabel('frames')
% ylabel('post stim compgp.cereb_diff_eachsub_sem - pre stim comp')
% title(sprintf('no stimulation'));
% goodplot


%% group difference for sham and cereb
% 
% % calc sems
% cereb_diff_sem=std(gp.cereb_diff_eachsub)/length(gp.precereb.patient_dat.pert_resp);
% cereb_diff_sem_pos=cereb_diff_sem+mean(gp.cereb_diff_eachsub);
% cereb_diff_sem_neg=mean(gp.cereb_diff_eachsub)-cereb_diff_sem;
% 
% sham_diff_sem=std(gp.sham_diff_eachsub)/length(gp.presham.patient_dat.pert_resp);
% sham_diff_sem_pos=sham_diff_sem+mean(gp.sham_diff_eachsub);
% sham_diff_sem_neg=mean(gp.sham_diff_eachsub)-sham_diff_sem;
% 
% beh_diff_sem=std(gp.beh_diff_eachsub)/length(gp.prebeh.patient_dat.pert_resp);
% beh_diff_sem_pos=beh_diff_sem+mean(gp.beh_diff_eachsub);
% beh_diff_sem_neg=mean(gp.beh_diff_eachsub)-beh_diff_sem;
% 
% figure
% plot(frame_taxis(idxes4anal), mean(gp.cereb_diff_eachsub), 'b','LineWidth', 3);
% hold on
% plot(frame_taxis(idxes4anal), cereb_diff_sem_pos, 'b','LineWidth', 1);
% plot(frame_taxis(idxes4anal), cereb_diff_sem_neg, 'b','LineWidth', 1);
% 
% plot(frame_taxis(idxes4anal), mean(gp.sham_diff_eachsub), 'm','LineWidth', 3);
% plot(frame_taxis(idxes4anal), sham_diff_sem_pos, 'm','LineWidth', 1);
% plot(frame_taxis(idxes4anal), sham_diff_sem_neg, 'm','LineWidth', 1);
% 
% plot(frame_taxis(idxes4anal), mean(gp.beh_diff_eachsub), 'c','LineWidth', 3);
% plot(frame_taxis(idxes4anal),beh_diff_sem_pos, 'c','LineWidth', 1);
% plot(frame_taxis(idxes4anal), beh_diff_sem_neg, 'c','LineWidth', 1);
% legend('cereb', ' ', ' ', 'sham', ' ', ' ','nostim', ' ', ' ')
% goodplot
% title(sprintf('Cerebellar Vs Vertex stimulation'));
% xlabel('Time (ms)')
% ylabel('post stim comp - pre stim comp')
% axis([-0.2 1 -60 30])
% 
% saveas(gcf, sprintf('groupdiff%d.jpg'))


%% check EW and LW trials
%check
% subplot(311); plot(gp.cereb_diff_eachsub')
% hold
% plot(mean(gp.cereb_diff_eachsub)', 'k', 'LineWidth', 1.3)
% plot((mean(gp.cereb_diff_eachsub))+gp.cereb_diff_eachsub_sem, 'LineWidth', 1.2)
% plot((mean(gp.cereb_diff_eachsub))-gp.cereb_diff_eachsub_sem, 'LineWidth', 1.2)
% axis([ 0 450 -25 25])
% 
% subplot(312); plot(gp.cereb_diff_eachsub(:,EW)')
% hold 
% plot(mean(gp.cereb_diff_eachsub(:,EW))', 'k', 'LineWidth', 1.3)
% axis([ 0 450 -25 25])
% 
% subplot(313); plot(gp.cereb_diff_eachsub(:,LW)')
% hold
% plot(mean(gp.cereb_diff_eachsub(:,LW))', 'k', 'LineWidth', 1.3)
% axis([ 0 450 -25 25])



% % do running ttest across all difference trials. 
% ttest_compdiff=[];
% for iframe=1:413
%     ttest_compdiff(:, iframe)=ttest(gp.cereb_diff_eachsub(:,iframe), gp.sham_diff_eachsub(:,iframe));            
% end
% 
% figure
% plot(frame_taxis(idxes4anal),ttest_compdiff, 'k', 'LineWidth',1.3);
% hold on
% plot(frame_taxis(idxes4anal),ttest(gp.cereb_diff_eachsub, gp.sham_diff_eachsub), 'r','LineWidth',1.3);
% 
