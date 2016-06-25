%% loads in the mean compensation from per and post sessions, takes means, plots and does a ttest
% ZKA Feb 2016


%%

pre_mean_comp_cereb_WithinSub= [0.2185	0.2465	0.2972	0.232	0.328	0.4943	-0.1104	0.1546	0.46	0.427	0.12	0.227	0.417];
post_mean_comp_cereb_WithinSub= [0.2853	0.3561	0.4374	0.1193	0.11	0.508	0.0739	0.2332	0.397	0.27	-0.009	0.25	0.212];
pre_mean_comp_sham_WithinSub= [0.2763	0.3823	0.319	0.1503	0.2059	0.2443	0.1193	0.2821	-0.0407	-0.0407	0.47	0.1158	0.2536];
post_mean_comp_sham_WithinSub= [0.1662	0.3393	0.3078	0.1491	0.1794	0.4354	0.1643	0.2456	0.1454	0.43	0.389	0.252	0.258];


pre_std_comp_cereb_WithinSub= [0.0448	0.038	0.0434	0.0415	0.0735	0.0302	0.1227	0.0341	0.04	0.038	0.03	0.057	0.043];
post_std_comp_cereb_WithinSub= [0.0602	0.0442	0.0338	0.0401	0.04	0.03	0.0353	0.0379	0.035	0.03	0.024	0.035	0.065];
pre_std_comp_sham_WithinSub= [0.0531	0.0407	0.0329	0.0551	0.0666	0.033	0.035	0.0399	0.0224	0.0224	0.03	0.064	0.057];
post_std_comp_sham_WithinSub= [0.0612	0.0415	0.0461	0.0845	0.0546	0.0855	0.0391	0.0353	0.0515	0.039	0.033	0.035	0.046];

anovandata_meancomp_WithinSub=[pre_mean_comp_cereb_WithinSub post_mean_comp_cereb_WithinSub pre_mean_comp_sham_WithinSub post_mean_comp_sham_WithinSub];
pre=ones(1,13);
post=pre*2;

meancomp_anova_WithinSub.group1_session=[pre post pre post];
meancomp_anova_WithinSub.group2_site=[ones(1,26) (ones(1,26)*2)];
meancomp_anova_WithinSub.group3_subject=[1:13 1:13 1:13 1:13];

[meancomp_anova_WithinSub.p_interaction,meancomp_anova.table,meancomp_anova.stats,meancomp_anova.terms]=...
    anovan(anovandata_meancomp_WithinSub,{meancomp_anova_WithinSub.group1_session meancomp_anova_WithinSub.group2_site meancomp_anova_WithinSub.group3_subject},'model','interaction','varnames',{'session','site','subject'})

[meancomp_anova_WithinSub.p_interaction,meancomp_anova.table,meancomp_anova.stats,meancomp_anova.terms]=...
    anovan(anovandata_meancomp_WithinSub,{meancomp_anova_WithinSub.group1_session meancomp_anova_WithinSub.group2_site},'model','full','varnames',{'session','site'})


% anovan(y,{g1 g2 g3},'model','interaction','varnames',{'g1','g2','g3'})

[meancomp_anova_WithinSub.p_full,meancomp_anova.table,meancomp_anova.stats,meancomp_anova.terms]=...
    anovan(anovandata_meancomp_WithinSub,{meancomp_anova_WithinSub.group1_session meancomp_anova_WithinSub.group2_site meancomp_anova_WithinSub.group3_subject},'full')


%% Mean whole trial compensation
pre_mean_comp_cereb= [0.1317	0.2853	0.2708	0.3561	0.4374	0.1193	0.11	0.508	0.0739	0.2332	0.397	0.27	-0.009	0.25	0.212];
post_mean_comp_cereb=[0.1114	0.2185	0.2511	0.2465	0.2972	0.232	0.328	0.4943	-0.1104	0.1546	0.46	0.427	0.12	0.227	0.417];
pre_mean_comp_sham= [0.1662	0.3393 0.3078 0.1491 0.1794	0.4354 0.1643 0.2456 0.1454	0.43 0.389 0.252 0.258];
pre_mean_comp_sham= [0.1662 0.3393 0.3078 0.1491 0.1794 0.4354 0.1643 0.2456 0.1454 0.43 0.389 0.252 0.258];

post_mean_comp_sham= [0.2763	0.3823	0.319	0.1503	0.2059	0.2443	0.1193	0.2821	-0.0407	-0.0407	0.47	0.1158	0.2536];

anovandata_meancomp=[pre_mean_comp_cereb post_mean_comp_cereb pre_mean_comp_sham post_mean_comp_sham]

pre_cereb=ones(1,15)
post_cereb=pre_cereb*2
pre_sham=ones(1,13)
post_sham=pre_sham*2

meancomp_anova.group1_session=[pre_cereb post_cereb pre_sham post_sham]
meancomp_anova.group2_site=[ones(1,30) (ones(1,26)*2)]
meancomp_anova.group3_subject=[1:15 1:15 1:13 1:13]

[meancomp_anova.p_interaction,meancomp_anova.table,meancomp_anova.stats,meancomp_anova.terms]=...
    anovan(anovandata_meancomp,{meancomp_anova.group1_session meancomp_anova.group2_site meancomp_anova.group3_subject},'model','interaction')

display 'first anova is the interaction'

save /Users/zagnew/Cereb_data/cTBS/data_analysis/stats/meancomp_anova meancomp_anova

[meancomp_anova_full.p_full,meancomp_anova.table,meancomp_anova_full.stats,meancomp_anova_full.terms]...
    = anovan(anovandata_meancomp,{meancomp_anova.group1_session meancomp_anova.group2_site meancomp_anova.group3_subject}, 'full')

display 'second anova is the full model'

save /Users/zagnew/Cereb_data/cTBS/data_analysis/stats/meancomp_anova_full meancomp_anova_full

% plot this
figure
meanmeancomp=[mean(pre_mean_comp_cereb) mean(post_mean_comp_cereb) ;mean(pre_mean_comp_sham) mean(post_mean_comp_sham)]
semmeancomp=[std(pre_mean_comp_cereb)/sqrt(15) std(post_mean_comp_cereb)/sqrt(15) ;std(pre_mean_comp_sham)/sqrt(13) std(post_mean_comp_sham)/sqrt(13)]

h = barwitherr(semmeancomp, meanmeancomp);% Plot with errorbars
set(gca,'XTickLabel',{'Cereb','Sham'})
ylabel('mean compensation (cents)')
set(h(1),'FaceColor',clear_colour,'EdgeColor', clear_colour ,'LineWidth',1.5);
set(h(2),'FaceColor',masked_colour,'EdgeColor', masked_colour ,'LineWidth',1.5);
title(sprintf('Compensation to PP before and after TBS'));
goodplot
legend('pre stim', 'post stim')
print(gcf, '-dpdf', '-r150', '/Users/zagnew/Cereb_data/cTBS/figures/meancomp.pdf');

%% 

load /Users/zagnew/cTBS_data/sham_TMS_pre_data/s06_pre/speak_consolidate_audiodir/TMS1_sham.mat;
TMS1_sham_a=pert_resp.comp{1};
TMS1_sham_b=pert_resp.comp{2};
TMS1_sham=vertcat(TMS1_sham_a, TMS1_sham_b);

load /Users/zagnew/cTBS_data/cereb_TMS_pre_data/s06_pre/speak_consolidate_audiodir/TMS1.mat;
TMS1_cereb_a=pert_resp.comp{1};
TMS1_cereb_b=pert_resp.comp{2};
TMS1_cereb=vertcat(TMS1_cereb_a, TMS1_cereb_b);

load /Users/zagnew/cTBS_data/sham_TMS_post_data/s06_post/speak_consolidate_audiodir/TMS2_sham.mat;
TMS2_sham_a=pert_resp.comp{1};
TMS2_sham_b=pert_resp.comp{2};
TMS2_sham=vertcat(TMS2_a, TMS2_b);

load /Users/zagnew/cTBS_data/cereb_TMS_post_data/s06_post/speak_consolidate_audiodir/TMS2.mat;
TMS2_cereb_a=pert_resp.comp{1};
TMS2_cereb_b=pert_resp.comp{2};
TMS2_cereb=vertcat(TMS2_cereb_a, TMS2_cereb_b);

% anova

anovadata= [TMS1_sham' TMS2_sham' TMS1_cereb' TMS2_cereb'];

% 67
% 74
% 67
% 60

test=ones(1,141);
test2=ones(1,127);
test2=test2*2;
site=[test test2];

%create conditions
condition= cell(1,268);
for i=1:67
    condition{i} = 'pre'; 
end
for i=68:141
    condition{i} = 'post';
end
for i=142:208
    condition{i} = 'pre';
end
for i=209:268
    condition{i} = 'post';
end

condition=condition';

group1=[site];
group2=[condition];
p = anovan(anovadata,{group1 group2 },'model','interaction')


cond_sham=[mean(TMS1_sham(1:length(TMS1_sham)));mean(TMS2_sham)];
errY2_sham=[std(TMS1_sham(1:length(TMS1_sham))/sqrt(length(TMS1_sham)));std(TMS2_sham)/sqrt(length(TMS1_sham))];
cond_cereb=[mean(TMS1_cereb(1:length(TMS1_cereb)));mean(TMS2_cereb)];
errY2_cereb=[std(TMS1_cereb(1:length(TMS1_cereb))/sqrt(length(TMS1_cereb)));std(TMS2_cereb)/sqrt(length(TMS1_cereb))];

figure
subplot(221)
h = barwitherr(errY2_sham, cond_sham);% Plot with errorbars
set(gca,'XTickLabel',{'Pre Stim','Post Stim'})
ylabel('mean compensation (cents)')
set(h(1),'FaceColor','w');
title(sprintf('Compensation to PP before and after sham TBS'));
goodplot

% scatter plot
subplot(222)
scatter(ones(1,length(TMS1_sham)),TMS1_sham,'.','k')
hold
scatter(1, mean(TMS1_sham(1:length(TMS1_sham))), 'filled', 'm')
scatter(2*(ones(1,length(TMS2_sham))),TMS2_sham,'.','k')
scatter(2, mean(TMS2_sham(1:length(TMS2_sham))), 'filled', 'm')
title(sprintf('Compensation to PP before and after sham TBS'));
axis([0 3, -1 1])
goodplot

subplot(223)
h = barwitherr(errY2_cereb, cond_cereb);% Plot with errorbars
set(gca,'XTickLabel',{'Pre Stim','Post Stim'})
ylabel('mean compensation (cents)')
set(h(1),'FaceColor','w');
title(sprintf('Compensation to PP before and after cerebellar TBS'));
goodplot

% scatter plot
subplot(224)
scatter(ones(1,length(TMS1_cereb)),TMS1_cereb,'.','k')
hold
scatter(1, mean(TMS1_cereb(1:length(TMS1_cereb))), 'filled', 'm')
scatter(2*(ones(1,length(TMS2_cereb))),TMS2_cereb,'.','k')
scatter(2, mean(TMS2_cereb(1:length(TMS2_cereb))), 'filled', 'm')
title(sprintf('Compensation to PP before and after cerebellar TBS'));
axis([0 3, -1 1])
goodplot

print(gcf, '-dpdf', '-r150', '/Users/zagnew/cTBS_data/cereb_TMS_post_data/s05_post/PP_TMS.pdf');

ttest2(TMS1_sham, TMS2_sham)

ttest2(TMS1_cereb, TMS2_cereb)


figure
