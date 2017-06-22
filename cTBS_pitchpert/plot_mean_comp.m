%% loads in the mean compensation from per nad post sessions, takes means, plots and does a ttest
% ZKA Feb 2015


%load /Users/zagnew/Cereb_data/TMS_analysis/PreStim/speak_consolidate_audiodir/TMS1.mat;

load /Users/zagnew/cTBS_data/sham_TMS_pre_data/s06_pre/speak_consolidate_audiodir/TMS1.mat;
TMS1_a=pert_resp.comp{1};
TMS1_b=pert_resp.comp{2};
TMS1=vertcat(TMS1_a, TMS1_b);

%load /Users/zagnew/Cereb_data/TMS_analaysis/PostStim/speak_consolidate_audiodir/TMS2.mat;
load /Users/zagnew/cTBS_data/sham_TMS_post_data/s06_post/speak_consolidate_audiodir/TMS2.mat;
TMS2_a=pert_resp.comp{1};
TMS2_b=pert_resp.comp{2};
TMS2=vertcat(TMS2_a, TMS2_b);

cond=[mean(TMS1(1:length(TMS1)));mean(TMS2)];
errY2=[std(TMS1(1:length(TMS1))/sqrt(length(TMS1)));std(TMS2)/sqrt(length(TMS1))];

figure
subplot(211)
h = barwitherr(errY2, cond);% Plot with errorbars
set(gca,'XTickLabel',{'Pre Stim','Post Stim'})
ylabel('mean compensation (cents)')
set(h(1),'FaceColor','w');
title(sprintf('Compensation to PP before and after cerebellar TBS'));
goodplot

% scatter plot
subplot(212)
scatter(ones(1,length(TMS1)),TMS1,'.','k')
hold
scatter(1, mean(TMS1(1:length(TMS1))), 'filled', 'm')

scatter(2*(ones(1,length(TMS2))),TMS2,'.','k')

scatter(2, mean(TMS2(1:length(TMS2))), 'filled', 'm')
title(sprintf('Compensation to PP before and after cerebellar TBS'));
axis([0 3, -1 1])
goodplot

print(gcf, '-dpdf', '-r150', '/Users/zagnew/cTBS_data/cereb_TMS_post_data/s05_post/PP_TMS.pdf');

ttest2(TMS1, TMS2)

