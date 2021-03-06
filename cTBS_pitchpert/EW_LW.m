% early and late windows

% EW=100:250;
% LW=250:400;


EW=100:199;
MW=200:299;
LW=300:399;

% ewline=[ones(1,EW(1))*-10 ones(1, 150)*30 ones(1,163)*-10]
% lwline=[ones(1,LW(1))*-10 ones(1, 150)*30 ones(1,13)*-10]

gp.cereb_diff_eachsub_EW=gp.cereb_diff_eachsub(:,EW);
gp.cereb_diff_eachsub_MW=gp.cereb_diff_eachsub(:,MW);
gp.cereb_diff_eachsub_LW=gp.cereb_diff_eachsub(:,LW);

gp.cereb_diff_eachsub_EW_mean=mean(gp.cereb_diff_eachsub(:,EW));
gp.cereb_diff_eachsub_MW_mean=mean(gp.cereb_diff_eachsub(:,MW));
gp.cereb_diff_eachsub_LW_mean=mean(gp.cereb_diff_eachsub(:,LW));
gp.cereb_diff_eachsub_LW

gp.sham_diff_eachsub_EW=gp.sham_diff_eachsub(:,EW);
gp.sham_diff_eachsub_MW=gp.sham_diff_eachsub(:,MW);
gp.sham_diff_eachsub_LW=gp.sham_diff_eachsub(:,LW);
gp.sham_diff_eachsub_EW_mean=mean(gp.sham_diff_eachsub(:,EW));
gp.sham_diff_eachsub_MW_mean=mean(gp.sham_diff_eachsub(:,MW));
gp.sham_diff_eachsub_LW_mean=mean(gp.sham_diff_eachsub(:,LW));



% subject means for anova
for isub=1: 13
sub_cereb_diff_EW(isub)=nanmean(gp.cereb_diff_eachsub(isub,EW))
sub_sham_diff_EW(isub)=nanmean(gp.sham_diff_eachsub(isub,EW))

sub_cereb_diff_MW(isub)=nanmean(gp.cereb_diff_eachsub(isub,MW))
sub_sham_diff_MW(isub)=nanmean(gp.sham_diff_eachsub(isub,MW))

sub_cereb_diff_LW(isub)=nanmean(gp.cereb_diff_eachsub(isub,LW))
sub_sham_diff_LW(isub)=nanmean(gp.sham_diff_eachsub(isub,LW))

end




% gp.beh_diff_eachsub_EW=gp.beh_diff_eachsub(:,EW);
% gp.beh_diff_eachsub_LW=gp.beh_diff_eachsub(:,LW);
% gp.beh_diff_eachsub_EW_mean=mean(gp.beh_diff_eachsub(:,EW));
% gp.beh_diff_eachsub_LW_mean=mean(gp.beh_diff_eachsub(:,LW));

% MAIN RESULTS FIGURE
figure
whitebg('white')
annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'mean difference in compensation, early and late windows', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
% EW_LW_pertcomp=[mean(gp.cereb_diff_eachsub_EW_mean) mean(gp.cereb_diff_eachsub_LW_mean) ; mean(gp.sham_diff_eachsub_EW_mean) mean(gp.sham_diff_eachsub_LW_mean) ; mean(gp.beh_diff_eachsub_EW_mean) mean(gp.beh_diff_eachsub_LW_mean)];
% EW_LWerr = [std(gp.cereb_diff_eachsub_EW_mean)/sqrt(length(EW)) std(gp.cereb_diff_eachsub_LW_mean)/sqrt(length(EW)) ; std(gp.sham_diff_eachsub_EW_mean)/sqrt(length(EW)) std(gp.sham_diff_eachsub_LW_mean)/sqrt(length(EW)) ; std(gp.beh_diff_eachsub_EW_mean)/sqrt(length(EW)) std(gp.beh_diff_eachsub_LW_mean)/sqrt(length(EW))];
EW_LW_pertcomp=[mean(gp.cereb_diff_eachsub_EW_mean) mean(gp.cereb_diff_eachsub_MW_mean) mean(gp.cereb_diff_eachsub_LW_mean) ; mean(gp.sham_diff_eachsub_EW_mean) mean(gp.sham_diff_eachsub_MW_mean) mean(gp.sham_diff_eachsub_LW_mean)];
EW_LWerr = [std(gp.cereb_diff_eachsub_EW_mean)/sqrt(length(EW)) std(gp.cereb_diff_eachsub_MW_mean)/sqrt(length(MW)) std(gp.cereb_diff_eachsub_LW_mean)/sqrt(length(LW)) ; std(gp.sham_diff_eachsub_EW_mean)/sqrt(length(EW)) std(gp.sham_diff_eachsub_MW_mean)/sqrt(length(MW)) std(gp.sham_diff_eachsub_LW_mean)/sqrt(length(EW))];

h = barwitherr(EW_LWerr, EW_LW_pertcomp);% Plot with errorbars

set(gca,'XTickLabel',{'cerebellar','vertex'})
ylabel('Difference in Compensation to Perturbation (cents)')
set(h(1),'FaceColor',clear_colour,'EdgeColor', clear_colour ,'LineWidth',1.5);
set(h(2),'FaceColor',masked_colour,'EdgeColor', masked_colour ,'LineWidth',1.5);
set(h(3),'FaceColor',mediumgrey,'EdgeColor', mediumgrey ,'LineWidth',1.5);
goodplot
legend('early window', 'mid window', 'late window','Location','SouthWest')
print(gcf, '-dpdf', '-r150', '/Users/zagnew/cerebellarTBS/data_analysis/figures/early_mid_late_winds_diff_comp.pdf');




% 
% figure
% subplot(211)
% plot(frame_taxis,nanmean(gp.meanpitchpertresp_precereb)','r','LineWidth',3);
% hold 
% plot(frame_taxis, nanmean(gp.meanpitchpertresp_presham)','Color',masked_colour,'LineWidth',3);
% %plot(frame_taxis, nanmean(gp.meanpitchpertresp_prebeh)','Color',clear_colour,'LineWidth',3);
% plot(frame_taxis, ewline,'k','LineWidth',2);
% plot(frame_taxis, lwline,'k','LineWidth',2);
% goodplot
% legend('cerebellar', 'vertex', 'no stim', 'Location','NorthWest')
% 
% subplot(212)
% plot(frame_taxis,nanmean(gp.meanpitchpertresp_postcereb)','r','LineWidth',3);
% hold 
% plot(frame_taxis,nanmean(gp.meanpitchpertresp_postsham)','Color',masked_colour,'LineWidth',3);
% %plot(frame_taxis,nanmean(gp.meanpitchpertresp_postbeh)','Color',clear_colour,'LineWidth',3);
% plot(frame_taxis, ewline,'k','LineWidth',2);
% plot(frame_taxis, lwline,'k','LineWidth',2);
% 
% print(gcf, '-dpdf', '-r150', '/Users/zagnew/cerebellarTBS/data_analysis/figures/preandpoststim responses.pdf');


%% EW MW LW anova
anova_early_mid_late= ...
    [gp.cereb_diff_eachsub_EW(1,:) gp.cereb_diff_eachsub_EW(2,:) gp.cereb_diff_eachsub_EW(3,:) gp.cereb_diff_eachsub_EW(4,:) ...
    gp.cereb_diff_eachsub_EW(5,:) gp.cereb_diff_eachsub_EW(6,:) gp.cereb_diff_eachsub_EW(7,:) gp.cereb_diff_eachsub_EW(8,:)...
    gp.cereb_diff_eachsub_EW(9,:) gp.cereb_diff_eachsub_EW(10,:) gp.cereb_diff_eachsub_EW(11,:) gp.cereb_diff_eachsub_EW(12,:) ...
    gp.cereb_diff_eachsub_EW(13,:) gp.cereb_diff_eachsub_MW(1,:) gp.cereb_diff_eachsub_MW(2,:) gp.cereb_diff_eachsub_MW(3,:) gp.cereb_diff_eachsub_MW(4,:) ...
    gp.cereb_diff_eachsub_MW(5,:) gp.cereb_diff_eachsub_MW(6,:) gp.cereb_diff_eachsub_MW(7,:) gp.cereb_diff_eachsub_MW(8,:)...
    gp.cereb_diff_eachsub_MW(9,:) gp.cereb_diff_eachsub_MW(10,:) gp.cereb_diff_eachsub_MW(11,:) gp.cereb_diff_eachsub_MW(12,:) ...
    gp.cereb_diff_eachsub_MW(13,:) gp.cereb_diff_eachsub_LW(1,:) gp.cereb_diff_eachsub_LW(2,:) gp.cereb_diff_eachsub_LW(3,:) gp.cereb_diff_eachsub_LW(4,:) ...
    gp.cereb_diff_eachsub_LW(5,:) gp.cereb_diff_eachsub_LW(6,:) gp.cereb_diff_eachsub_LW(7,:) gp.cereb_diff_eachsub_LW(8,:)...
    gp.cereb_diff_eachsub_LW(9,:) gp.cereb_diff_eachsub_LW(10,:) gp.cereb_diff_eachsub_LW(11,:) gp.cereb_diff_eachsub_LW(12,:) ...
    gp.cereb_diff_eachsub_LW(13,:) gp.sham_diff_eachsub_EW(1,:) gp.sham_diff_eachsub_EW(2,:) gp.sham_diff_eachsub_EW(3,:) gp.sham_diff_eachsub_EW(4,:) ...
    gp.sham_diff_eachsub_EW(5,:) gp.sham_diff_eachsub_EW(6,:) gp.sham_diff_eachsub_EW(7,:) gp.sham_diff_eachsub_EW(8,:)...
    gp.sham_diff_eachsub_EW(9,:) gp.sham_diff_eachsub_EW(10,:) gp.sham_diff_eachsub_EW(11,:) gp.sham_diff_eachsub_EW(12,:) ...
    gp.sham_diff_eachsub_EW(13,:) gp.sham_diff_eachsub_MW(1,:) gp.sham_diff_eachsub_MW(2,:) gp.sham_diff_eachsub_MW(3,:) gp.sham_diff_eachsub_MW(4,:) ...
    gp.sham_diff_eachsub_MW(5,:) gp.sham_diff_eachsub_MW(6,:) gp.sham_diff_eachsub_MW(7,:) gp.sham_diff_eachsub_MW(8,:)...
    gp.sham_diff_eachsub_MW(9,:) gp.sham_diff_eachsub_MW(10,:) gp.sham_diff_eachsub_MW(11,:) gp.sham_diff_eachsub_MW(12,:) ...
    gp.sham_diff_eachsub_MW(13,:) gp.sham_diff_eachsub_LW(1,:) gp.sham_diff_eachsub_LW(2,:) gp.sham_diff_eachsub_LW(3,:) gp.sham_diff_eachsub_LW(4,:) ...
    gp.sham_diff_eachsub_LW(5,:) gp.sham_diff_eachsub_LW(6,:) gp.sham_diff_eachsub_LW(7,:) gp.sham_diff_eachsub_LW(8,:)...
    gp.sham_diff_eachsub_LW(9,:) gp.sham_diff_eachsub_LW(10,:) gp.sham_diff_eachsub_LW(11,:) gp.sham_diff_eachsub_LW(12,:) ...
    gp.sham_diff_eachsub_LW(13,:)];

alldata=length(anova_early_mid_late);
eachsub=length(EW);
subs=13;
window_1=ones(1,eachsub*subs);
window_2=window_1*2;
window_3=window_1*3;
window=[window_1 window_2 window_3 window_1 window_2 window_3];

cereb=ones(1,eachsub*subs);
sham=cereb*2;
stim_site=[cereb cereb cereb sham sham sham];

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
subjectgroup=[subjectgroup subjectgroup subjectgroup subjectgroup subjectgroup subjectgroup];
 
group1=subjectgroup;
group2=window;
group3=stim_site;

% [pvals,tbl,stats]  = anovan(anova_early_mid_late,{group1 group2 group3},'model','interaction');
% [pvals,tbl,stats]  = anovan(anova_early_mid_late,{group2 group3},'model',1,'varnames',{'window' 'stim site'});
% [pvals,tbl,stats]  = anovan(anova_early_mid_late,{group2 group3},'model',2,'varnames',{'window' 'stim site'});
% [pvals,tbl,stats]  = anovan(anova_early_mid_late,{group2 group3},'model',3,'varnames',{'window' 'stim site'});
% [pvals,tbl,stats]  = anovan(anova_early_mid_late,{group1 group2 group3},'model',3,'varnames',{'subject','window' 'stim site'});
% 
% anovan(anova_early_mid_late, group1)
% p = anovan(anova_early_mid_late,{group1 group2 group3})

% this is the one
p = anovan(anova_early_mid_late,{group2 group3},'model','interaction','varnames',{'window' 'stim site'})


% but why is the
%p = anovan(anova_early_mid_late,{group1 group2 group3},'model','interaction','varnames',{'subject','window' 'stim site'})









%% early and late window anova including two stim conditions only
 anova_early_late= ...
    [gp.cereb_diff_eachsub_EW(1,:) gp.cereb_diff_eachsub_EW(2,:) gp.cereb_diff_eachsub_EW(3,:) gp.cereb_diff_eachsub_EW(4,:) ...
    gp.cereb_diff_eachsub_EW(5,:) gp.cereb_diff_eachsub_EW(6,:) gp.cereb_diff_eachsub_EW(7,:) gp.cereb_diff_eachsub_EW(8,:)...
    gp.cereb_diff_eachsub_EW(9,:) gp.cereb_diff_eachsub_EW(10,:) gp.cereb_diff_eachsub_EW(11,:) gp.cereb_diff_eachsub_EW(12,:) ...
    gp.cereb_diff_eachsub_EW(13,:) gp.cereb_diff_eachsub_LW(1,:) gp.cereb_diff_eachsub_LW(2,:) gp.cereb_diff_eachsub_LW(3,:) gp.cereb_diff_eachsub_LW(4,:) ...
    gp.cereb_diff_eachsub_LW(5,:) gp.cereb_diff_eachsub_LW(6,:) gp.cereb_diff_eachsub_LW(7,:) gp.cereb_diff_eachsub_LW(8,:)...
    gp.cereb_diff_eachsub_LW(9,:) gp.cereb_diff_eachsub_LW(10,:) gp.cereb_diff_eachsub_LW(11,:) gp.cereb_diff_eachsub_LW(12,:) ...
    gp.cereb_diff_eachsub_LW(13,:) gp.sham_diff_eachsub_EW(1,:) gp.sham_diff_eachsub_EW(2,:) gp.sham_diff_eachsub_EW(3,:) gp.sham_diff_eachsub_EW(4,:) ...
    gp.sham_diff_eachsub_EW(5,:) gp.sham_diff_eachsub_EW(6,:) gp.sham_diff_eachsub_EW(7,:) gp.sham_diff_eachsub_EW(8,:)...
    gp.sham_diff_eachsub_EW(9,:) gp.sham_diff_eachsub_EW(10,:) gp.sham_diff_eachsub_EW(11,:) gp.sham_diff_eachsub_EW(12,:) ...
    gp.sham_diff_eachsub_EW(13,:) gp.sham_diff_eachsub_LW(1,:) gp.sham_diff_eachsub_LW(2,:) gp.sham_diff_eachsub_LW(3,:) gp.sham_diff_eachsub_LW(4,:) ...
    gp.sham_diff_eachsub_LW(5,:) gp.sham_diff_eachsub_LW(6,:) gp.sham_diff_eachsub_LW(7,:) gp.sham_diff_eachsub_LW(8,:)...
    gp.sham_diff_eachsub_LW(9,:) gp.sham_diff_eachsub_LW(10,:) gp.sham_diff_eachsub_LW(11,:) gp.sham_diff_eachsub_LW(12,:) ...
    gp.sham_diff_eachsub_LW(13,:)];

alldata=length(anova_early_late);
eachsub=length(EW);
subs=13;

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
 
window_early=ones(1,eachsub*subs);
window_late=window_early*2;

window=[window_early window_late window_early window_late];

cereb=ones(1,eachsub*subs);
sham=cereb*2;
stim_site=[cereb cereb sham sham];
group1=subjectgroup;
group2=window;
group3=stim_site;

% p = anovan(anova_early_late,{group1 group2 group3},'model','interaction');


p = anovan(anova_early_late,{group2 group3},'model','interaction','varnames',{'window' 'stim site'})


[pvals,tbl,stats]  = anovan(anova_early_late,{group2 group3},'model',1,'varnames',{'window' 'stim site'});
[pvals,tbl,stats]  = anovan(anova_early_late,{group2 group3},'model',2,'varnames',{'window' 'stim site'});


[pvals,tbl,stats]  = anovan(anova_early_late,{group1 group2 group3},'model',2,'varnames',{'subject' 'window' 'stim site'});
[pvals,tbl,stats]  = anovan(anova_early_late,{group1 group2 group3},'model',3,'varnames',{'subject' 'window' 'stim site'});

[pvals,tbl,stats] = anovan(anova_early_late, {group1 group2 group3},'model',2, 'random',1,'varnames',{'subject' 'window' 'stim site'});
% multcompare(stats)





%% plot all EW and EW on top of each other
figure
subplot(311)
plot(frame_taxis(EW), gp.cereb_diff_eachsub_EW_mean,'k','LineWidth', 3)
hold 
plot(frame_taxis(EW), gp.cereb_diff_eachsub_LW_mean, 'k', 'LineWidth', 1)
goodplot
axis([ 0 0.6 -3 7])
title('cerebellar')

subplot(312)
plot(frame_taxis(EW),gp.sham_diff_eachsub_EW_mean,'k','LineWidth', 3)
hold
plot(frame_taxis(EW),gp.sham_diff_eachsub_LW_mean','k','LineWidth', 1)
axis([ 0 0.6 -3 7])
title('vertex')
goodplot

subplot(313)
plot(frame_taxis(EW),gp.beh_diff_eachsub_EW_mean, 'k','LineWidth', 3)
hold
plot(frame_taxis(EW),gp.beh_diff_eachsub_LW_mean', 'k','LineWidth', 1)
axis([ 0 0.6 -3 7])
goodplot
title('no stimulation')
print(gcf, '-dpdf', '-r150', '/Users/zagnew/cerebellarTBS/data_analysis/figures/EW_LW_diffincomp.pdf');
legend('early window', 'late window', 'Location','NorthEast')

