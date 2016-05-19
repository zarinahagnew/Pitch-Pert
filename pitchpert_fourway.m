% % ANOVA for group data
% ZKA Oct 2016
% saves the mean data from plot_data
% pert_resp.cents4comp.pitch_in.mean(1,:)
% 
% % Notes on pert.resp
% gp.precereb.patient_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{1}
% pitch_in.dat{1} - response to down pert
% pitch_in.dat{2} response to up pert
% pitch_in.dat{3} respone to all but *it's not flipped*
% to flip it use: 
% centsdev_dat.(group{igroup}).subj(isubj).absdat = {-the_dat{1} the_dat{2} [-the_dat{1}; the_dat{2}]};
% 
% DONE
% for the anova want to run for time point also, and put frame in as 413
% levels
% or 
% use mean .peakcomp .comp values post and pre and put into anova
% to see wherepitchpert_fourway this effect is in the trial look at mean(subdata_pre_cereb(isub).pert_resp.tpeak{1})
% 
% Also can do an anova on pert_resp.tpeak{1}) to see if the time point of
% the peak comp is different. 
%  
% in general, *comp* in name means that the cents response, divded by the cent of the pert and made positive. 
% 
% So look at the peak comp and t peak to get at both amp of comp and time
% peak of comp, and you can check that this fits with this once you have 
% concatednatepitchpert_fourwayd/added pitch_in{1} and {2} of patient_dat.pert_resp(isubj).cents4comp(1).pitch_in
% flip the first one (all negs become pos and all popatient_dat.pert_resps becomes neg)
% 
% = {-the_dat{1} the_dat{2} [-the_dat{1}; the_dat{2}]};
% 
% pert_resp.cents4pert_resp.cents4comp.pitch_in.meancomp.pitch_in.mean
% cd /data/bil-mb4/zarinah-data/cerebellar-data/pitch-pert-cTBS/cTBS_data/data_analysis

set(0,'DefaultFigureWindowStyle','docked');
clear all
close all
z_colours;

cd /Users/zagnew/Cereb_data/cTBS/data_analysis

%% for extra plotz
% pitchpert_fourway_extraplotz;


%% make sure that you have renames the patient.mat patient_pre.mat etc
gp.precereb=load('patient_pre');
gp.postcereb=load('patient_post');
gp.presham=load('control_pre');
gp.postsham=load('control_post');

% for plottin' 
frame_taxis = gp.precereb.patient_dat.frame_taxis;
tlims4anal_req = [-0.2 1.0];
ilims4anal = dsearchn(frame_taxis',tlims4anal_req')';
idxes4anal = ilims4anal(1):ilims4anal(2);
tlims4anal = frame_taxis(ilims4anal);

% flip pitch and concatenate pitch_in.dat{1} and {2}
for isubj=1:13 %length(gp.postcereb.patient_dat.pert_resp)
    gp.precereb.patient_dat.pert_resp(isubj).cents4comp(1).temp.dat=-gp.precereb.patient_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{1};
    gp.precereb.patient_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat = [gp.precereb.patient_dat.pert_resp(isubj).cents4comp(1).temp.dat ; gp.precereb.patient_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{2}];    
    gp.postcereb.patient_dat.pert_resp(isubj).cents4comp(1).temp.dat=-gp.postcereb.patient_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{1};
    gp.postcereb.patient_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat = [gp.postcereb.patient_dat.pert_resp(isubj).cents4comp(1).temp.dat ; gp.postcereb.patient_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{2}];    
    gp.presham.control_dat.pert_resp(isubj).cents4comp(1).temp.dat=-gp.presham.control_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{1};
    gp.presham.control_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat = [gp.presham.control_dat.pert_resp(isubj).cents4comp(1).temp.dat ; gp.presham.control_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{2}];    
    gp.postsham.control_dat.pert_resp(isubj).cents4comp(1).temp.dat=-gp.postsham.control_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{1};
    gp.postsham.control_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat = [gp.postsham.control_dat.pert_resp(isubj).cents4comp(1).temp.dat ; gp.postsham.control_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{2}];    
end

% calc mean comp
for isubj=1:13 %length(gp.precereb.patient_dat.pert_resp)
    gp.meanpitchpertresp_precereb(isubj,:)=mean(gp.precereb.patient_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat);
    gp.meanpitchpertresp_postcereb(isubj,:)=mean(gp.postcereb.patient_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat);
    gp.meanpitchpertresp_presham(isubj,:)=mean(gp.presham.control_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat);
    gp.meanpitchpertresp_postsham(isubj,:)=mean(gp.postsham.control_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat);    
    gp.cereb_diff_eachsub(isubj,:)=gp.meanpitchpertresp_postcereb(isubj,:)-gp.meanpitchpertresp_precereb(isubj,:);
    gp.sham_diff_eachsub(isubj,:)=gp.meanpitchpertresp_postsham(isubj,:)-gp.meanpitchpertresp_presham(isubj,:);       
end

% calc sems
cereb_diff_sem=std(gp.cereb_diff_eachsub)/length(gp.precereb.patient_dat.pert_resp);
cereb_diff_sem_pos=cereb_diff_sem+mean(gp.cereb_diff_eachsub);
cereb_diff_sem_neg=mean(gp.cereb_diff_eachsub)-cereb_diff_sem;

sham_diff_sem=std(gp.sham_diff_eachsub)/length(gp.presham.control_dat.pert_resp);
sham_diff_sem_pos=sham_diff_sem+mean(gp.sham_diff_eachsub);
sham_diff_sem_neg=mean(gp.sham_diff_eachsub)-sham_diff_sem;

%% Final Figures

plot group difference in pre and post pitch pert
figure
plot(frame_taxis(idxes4anal), mean(gp.cereb_diff_eachsub), 'k','LineWidth', 3);
hold on
plot(frame_taxis(idxes4anal), cereb_diff_sem_pos, 'k','LineWidth', 1);
plot(frame_taxis(idxes4anal), cereb_diff_sem_neg, 'k','LineWidth', 1);

plot(frame_taxis(idxes4anal), mean(gp.sham_diff_eachsub), 'r','LineWidth', 3);
plot(frame_taxis(idxes4anal), sham_diff_sem_pos, 'r','LineWidth', 1);
plot(frame_taxis(idxes4anal), sham_diff_sem_neg, 'r','LineWidth', 1);

a=ttest2(gp.cereb_diff_eachsub, gp.sham_diff_eachsub);
a(a == 0) = NaN;
plot(frame_taxis(idxes4anal),a*6, 'c','LineWidth', 3);

%% abs version
figure
plot(frame_taxis(idxes4anal), abs(mean(gp.cereb_diff_eachsub)), 'k','LineWidth', 3);
hold on
plot(frame_taxis(idxes4anal), cereb_diff_sem_pos, 'k','LineWidth', 1);
plot(frame_taxis(idxes4anal), cereb_diff_sem_neg, 'k','LineWidth', 1);

plot(frame_taxis(idxes4anal), abs(mean(gp.sham_diff_eachsub)), 'r','LineWidth', 3);
plot(frame_taxis(idxes4anal), sham_diff_sem_pos, 'r','LineWidth', 1);
plot(frame_taxis(idxes4anal), sham_diff_sem_neg, 'r','LineWidth', 1);
axis([0 1.2 0 15])




goodplot
title(sprintf('Cerebellar Vs Vertex stimulation'));
xlabel('Time (ms)')
ylabel('post stim comp - pre stim comp')
axis([-0.2 1 -6 8])
legend('cerebellar stimulation')
saveas(gcf, sprintf('groupdiff%d.jpg'))

figure
plot(ttest(gp.cereb_diff_eachsub, gp.sham_diff_eachsub), 'k','LineWidth',1.3);
axis([0 450 -0.5 2])

plot(ttest2(gp.meanpitchpertresp_postcereb,gp.meanpitchpertresp_postsham))
moo=ttest2(mean(gp.cereb_diff_eachsub), mean(gp.sham_diff_eachsub));
plot(ttest2(gp.cereb_diff_eachsub, gp.sham_diff_eachsub))

figure
plot(ttest2(gp.cereb_diff_eachsub(:,1:207), gp.sham_diff_eachsub(:,1:207)))
plot(ttest2(gp.cereb_diff_eachsub(:,208:end), gp.sham_diff_eachsub(:,208:end)))
axis([0 250 -2 2])

%% anova
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
p = anovan(anovadata,{group1 group2 }, 'full');
[~,~,stats] = anovan(anovadata, {group1 group2},'model',2, 'random',1,'varnames',{'subject' 'stim site'});

% anovan(mileage, {factory carmod}, 'model',2, 'random',1,'varnames',{'Factory' 'Car Model'});
multcompare(stats)

% do running ttest across all difference trials. 
for iframe=1:413
    ttest_compdiff(:, iframe)=ttest(gp.cereb_diff_eachsub(:,iframe), gp.sham_diff_eachsub(:,iframe));            
end

figure
plot(frame_taxis(idxes4anal),ttest_compdiff, 'k', 'LineWidth',1.3);
hold on
plot(frame_taxis(idxes4anal),ttest(gp.cereb_diff_eachsub, gp.sham_diff_eachsub), 'r','LineWidth',1.3);


% calculate sems
gp.cereb_diff_eachsub_sem=std(gp.cereb_diff_eachsub)/13;
gp.sham_diff_eachsub_sem=std(gp.sham_diff_eachsub)/13;

% early and late windows
gp.cereb_diff_eachsub_EW=gp.cereb_diff_eachsub(:,1:206);
gp.cereb_diff_eachsub_LW=gp.cereb_diff_eachsub(:,208:end);
gp.cereb_diff_eachsub_EW_mean=mean(gp.cereb_diff_eachsub(:,1:206));
gp.cereb_diff_eachsub_LW_mean=mean(gp.cereb_diff_eachsub(:,208:end));

gp.sham_diff_eachsub_EW=gp.sham_diff_eachsub(:,1:206);
gp.sham_diff_eachsub_LW=gp.sham_diff_eachsub(:,208:end);
gp.sham_diff_eachsub_EW_mean=mean(gp.sham_diff_eachsub(:,1:206));
gp.sham_diff_eachsub_LW_mean=mean(gp.sham_diff_eachsub(:,208:end));


%check
subplot(311); plot(gp.cereb_diff_eachsub')
hold
plot(mean(gp.cereb_diff_eachsub)', 'k', 'LineWidth', 1.3)
plot((mean(gp.cereb_diff_eachsub))+gp.cereb_diff_eachsub_sem, 'LineWidth', 1.2)
plot((mean(gp.cereb_diff_eachsub))-gp.cereb_diff_eachsub_sem, 'LineWidth', 1.2)
axis([ 0 450 -25 25])

subplot(312); plot(gp.cereb_diff_eachsub(:,1:206)')
hold 
plot(mean(gp.cereb_diff_eachsub(:,1:206))', 'k', 'LineWidth', 1.3)
axis([ 0 450 -25 25])

subplot(313); plot(gp.cereb_diff_eachsub(:,208:end)')
hold
plot(mean(gp.cereb_diff_eachsub(:,208:end))', 'k', 'LineWidth', 1.3)
axis([ 0 450 -25 25])

figure
plot(gp.cereb_diff_eachsub_EW_mean,'k')
hold 
plot(gp.cereb_diff_eachsub_LW_mean, 'k')
plot(gp.sham_diff_eachsub_EW_mean, 'r')
plot(gp.sham_diff_eachsub_LW_mean', 'r')


figure
whitebg('white')
annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'mean difference in compensation, early and late windows', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
y_pertcomp=[mean(gp.cereb_diff_eachsub_EW_mean) mean(gp.cereb_diff_eachsub_LW_mean) ; mean(gp.sham_diff_eachsub_EW_mean) mean(gp.sham_diff_eachsub_LW_mean) ];
errY2 = [std(gp.cereb_diff_eachsub_EW_mean)/sqrt(206) std(gp.cereb_diff_eachsub_LW_mean)/sqrt(206) ; std(gp.sham_diff_eachsub_EW_mean)/sqrt(206) std(gp.sham_diff_eachsub_LW_mean)/sqrt(206) ];

h = barwitherr(errY2, y_pertcomp);% Plot with errorbars

set(gca,'XTickLabel',{'cerebellar','sham'})
ylabel('Difference in Compensation to Perturbation (cents)')
set(h(1),'FaceColor',clear_colour,'EdgeColor', clear_colour ,'LineWidth',1.5);
set(h(2),'FaceColor',masked_colour,'EdgeColor', masked_colour ,'LineWidth',1.5);
goodplot
legend('early', 'late','Location','SouthEast')
print(gcf, '-dpdf', '-r150', '/Users/zagnew/Cereb_data/cTBS/figures/early_late_winds_diff_comp.pdf');

% early and late window anova

anova_early_late= ...
    [gp.cereb_diff_eachsub_EW(1,:) gp.cereb_diff_eachsub_EW(2,:) gp.cereb_diff_eachsub_EW(3,:) gp.cereb_diff_eachsub_EW(4,:) ...
    gp.cereb_diff_eachsub_EW(5,:) gp.cereb_diff_eachsub_EW(6,:) gp.cereb_diff_eachsub_EW(7,:) gp.cereb_diff_eachsub_EW(8,:)...
    gp.cereb_diff_eachsub_EW(9,:) gp.cereb_diff_eachsub_EW(10,:) gp.cereb_diff_eachsub_EW(11,:) gp.cereb_diff_eachsub_EW(12,:) ...
    gp.cereb_diff_eachsub_EW(13,:) gp.cereb_diff_eachsub_LW(1,:) gp.cereb_diff_eachsub_LW(2,:) gp.cereb_diff_eachsub_LW(3,:) gp.cereb_diff_eachsub_LW(4,:) ...
    gp.cereb_diff_eachsub_LW(5,:) gp.cereb_diff_eachsub_LW(6,:) gp.cereb_diff_eachsub_LW(7,:) gp.cereb_diff_eachsub_LW(8,:)...
    gp.cereb_diff_eachsub_LW(9,:) gp.cereb_diff_eachsub_LW(10,:) gp.cereb_diff_eachsub_LW(11,:) gp.cereb_diff_eachsub_LW(12,:) ...
    gp.cereb_diff_eachsub_LW(13,:) ...
    gp.sham_diff_eachsub_EW(1,:) gp.sham_diff_eachsub_EW(2,:) gp.sham_diff_eachsub_EW(3,:) gp.sham_diff_eachsub_EW(4,:) ...
    gp.sham_diff_eachsub_EW(5,:) gp.sham_diff_eachsub_EW(6,:) gp.sham_diff_eachsub_EW(7,:) gp.sham_diff_eachsub_EW(8,:)...
    gp.sham_diff_eachsub_EW(9,:) gp.sham_diff_eachsub_EW(10,:) gp.sham_diff_eachsub_EW(11,:) gp.sham_diff_eachsub_EW(12,:) ...
    gp.sham_diff_eachsub_EW(13,:) gp.sham_diff_eachsub_LW(1,:) gp.sham_diff_eachsub_LW(2,:) gp.sham_diff_eachsub_LW(3,:) gp.sham_diff_eachsub_LW(4,:) ...
    gp.sham_diff_eachsub_LW(5,:) gp.sham_diff_eachsub_LW(6,:) gp.sham_diff_eachsub_LW(7,:) gp.sham_diff_eachsub_LW(8,:)...
    gp.sham_diff_eachsub_LW(9,:) gp.sham_diff_eachsub_LW(10,:) gp.sham_diff_eachsub_LW(11,:) gp.sham_diff_eachsub_LW(12,:) ...
    gp.sham_diff_eachsub_LW(13,:)];

alldata=10712;
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
subjectgroup=[subjectgroup subjectgroup subjectgroup subjectgroup];

window_early=ones(1,2678);
window_late=window_early*2;
window=[window_early window_late window_early window_late];

cereb=ones(1,206*26);
sham=cereb*2;
stim_site=[cereb sham];

group1=subjectgroup;
group2=window;
group3=stim_site;

p = anovan(anova_early_late,{group1 group2 group3},'model','interaction');
[pvals,tbl,stats]  = anovan(anova_early_late,{group1 group2 group3}, 'full');

%[pvals,tbl,stats] = anovan(anova_early_late, {group1 group2 group3},'model',2, 'random',1,'varnames',{'subject' 'window' 'stim site'});

multcompare(stats)

save GROUPDATA

% 
% cd data-for-r
% for isubj=1:13
% filename = [ 'cereb_prt_diff_EW' num2str(isubj) '.txt' ];
% dlmwrite(filename,gp.cereb_diff_eachsub_EW(isubj,:)')
% 
% filename = [ 'cereb_prt_diff_LW' num2str(isubj) '.txt' ];
% dlmwrite(filename,gp.cereb_diff_eachsub_LW(isubj,:)')
% 
% filename = [ 'sham_prt_diff_EW' num2str(isubj) '.txt' ];
% dlmwrite(filename,gp.sham_diff_eachsub_EW(isubj,:)')
% 
% filename = [ 'sham_prt_diff_LW' num2str(isubj) '.txt' ];
% dlmwrite(filename,gp.sham_diff_eachsub_LW(isubj,:)')
% end
