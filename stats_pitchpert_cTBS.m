%% stats for mean comp for cTBS study
% 15 cereb subs, 13 sham subs, 10 beh subs
% 38 subjects in total, pre and post so 76 data sets in all. 
% =========================================================================
clear all
z_params;
z_colours;

cd /Users/zagnew/cerebellarTBS
%load pitch_pert_cTBS_data.mat; < this is the raw data but has an extra
%data point in that needs removing at line 47. manually removed it and
%saved it as data.mat: 

load data.mat

data.subject{1:30}
data.subject{31:56}
data.subject{58:76}

data.subject{1:2:30};
data.subject{2:2:30};
data.subject{31:2:56};
data.subject(32:2:56);
data.subject(57:2:76);
data.subject(57:2:76);

pre_peakcomp_cereb=data.peakcomp(1:2:30);
post_peakcomp_cereb=data.peakcomp(2:2:30);
pre_peakcomp_sham=data.peakcomp(31:2:56);
post_peakcomp_sham=data.peakcomp(32:2:56);
pre_peakcomp_beh=data.peakcomp(57:2:76);
post_peakcomp_beh=data.peakcomp(57:2:76);

pre_tpeak_cereb=data.tpeak(1:2:30);
post_tpeak_cereb=data.tpeak(2:2:30);
pre_tpeak_sham=data.tpeak(31:2:56);
post_tpeak_sham=data.tpeak(32:2:56);
pre_tpeak_beh=data.tpeak(57:2:76);
post_tpeak_beh=data.tpeak(57:2:76);

%% exploratory ttests
[h,p,ci,stats] = ttest2(pre_peakcomp_cereb, post_peakcomp_cereb)
[h,p,ci,stats] = ttest2(pre_peakcomp_sham, post_peakcomp_sham)
[h,p,ci,stats] = ttest2(pre_peakcomp_sham, post_peakcomp_sham)
[h,p,ci,stats] = ttest2(pre_tpeak_cereb, pre_tpeak_sham)
[h,p,ci,stats] = ttest2(pre_tpeak_cereb, post_tpeak_cereb)
[h,p,ci,stats] = ttest2(pre_tpeak_sham, post_tpeak_sham)


% %%
% 
% x(:,1)=data.peakcomp'
% x(1:30,2)=1;
% 
% x(31:56,2)=2;
% x(57:76,2)=3;
% x(1:2:76,3)=1;
% x(2:2:76,3)=2;
% 
% w = 2;
% n = 76/2;
% v = repmat(1:n,[w 1]);
% x(:,4) = v(:)'
% 
% [SSQs, DFs, MSQs, Fs, Ps]=mixed_between_within_anova(x,1)
%% plot mean peak comp
figure
whitebg('white')
annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'mean comp', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
y_pitch2=[mean(pre_peakcomp_cereb) mean(post_peakcomp_cereb); ...
    mean(pre_peakcomp_sham) mean(post_peakcomp_sham);...
    mean(pre_peakcomp_beh) mean(post_peakcomp_beh)];

errY2 = [std(pre_peakcomp_cereb)/sqrt(length(pre_peakcomp_cereb)) std(post_peakcomp_cereb)/sqrt(length(post_peakcomp_cereb));...
    std(pre_peakcomp_sham)/sqrt(length(pre_peakcomp_sham)) std(post_peakcomp_sham)/sqrt(length(post_peakcomp_sham));...
    std(pre_peakcomp_beh)/sqrt(length(pre_peakcomp_beh)) std(post_peakcomp_beh)/sqrt(length(post_peakcomp_beh))];

h = barwitherr(errY2, y_pitch2);

set(gca,'XTickLabel',{'cereb','sham', 'beh'})
ylabel('Mean Comp')
set(h(1),'FaceColor',clear_colour,'EdgeColor', clear_colour ,'LineWidth',1.5);
set(h(2),'FaceColor',masked_colour,'EdgeColor', masked_colour ,'LineWidth',1.5);
goodplot
print(gcf, '-dpdf', '-r150', '/Users/zagnew/cerebellarTBS/figures/MeanPeakComp_bar.pdf');

%% plot tpeak
figure
whitebg('white')
annotation('textbox', [0 0.9 1 0.1], ...
    'String', 't peak', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
y_pitch2=[mean(pre_tpeak_cereb) mean(post_tpeak_cereb); ...
    mean(pre_tpeak_sham) mean(post_tpeak_sham);...
    mean(pre_tpeak_beh) mean(post_tpeak_beh)];

errY2 = [std(pre_tpeak_cereb)/sqrt(length(pre_tpeak_cereb)) std(post_tpeak_cereb)/sqrt(length(post_tpeak_cereb));...
    std(pre_tpeak_sham)/sqrt(length(pre_tpeak_sham)) std(post_tpeak_sham)/sqrt(length(post_tpeak_sham));...
    std(pre_tpeak_beh)/sqrt(length(pre_tpeak_beh)) std(post_tpeak_beh)/sqrt(length(post_tpeak_beh))];

h = barwitherr(errY2, y_pitch2);% Plot with errorbars

set(gca,'XTickLabel',{'cereb','sham', 'beh'})
ylabel('T peak')
set(h(1),'FaceColor',clear_colour,'EdgeColor', clear_colour ,'LineWidth',1.5);
set(h(2),'FaceColor',masked_colour,'EdgeColor', masked_colour ,'LineWidth',1.5);
goodplot
print(gcf, '-dpdf', '-r150', '/Users/zagnew/cerebellarTBS/figures/Meantpeak_bar.pdf');

%% plot against combined baseline
allpre=[pre_peakcomp_cereb pre_peakcomp_sham pre_peakcomp_beh]
allpremean=mean(allpre)
allpresem=allpremean/sqrt(length(allpre))

figure
whitebg('white')
annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'mean comp', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
y_pitch2=[mean(allpre) mean(post_peakcomp_cereb) mean(post_peakcomp_sham) mean(post_peakcomp_beh)];

errY2 = [std(allpre)/sqrt(length(allpre)) ...
    std(post_peakcomp_cereb)/sqrt(length(post_peakcomp_cereb))...     
    std(post_peakcomp_sham)/sqrt(length(post_peakcomp_sham))...
     std(post_peakcomp_beh)/sqrt(length(post_peakcomp_beh))];

h = barwitherr(errY2, y_pitch2);

set(gca,'XTickLabel',{'pre' 'cereb','sham', 'beh'})
ylabel('Mean Comp')
set(h(1),'FaceColor',clear_colour,'EdgeColor', clear_colour ,'LineWidth',1.5);
set(h(2),'FaceColor',masked_colour,'EdgeColor', masked_colour ,'LineWidth',1.5);
goodplot
print(gcf, '-dpdf', '-r150', '/Users/zagnew/cerebellarTBS/figures/MeanPeakComp_combined baseline_bar.pdf');



%% anova on EW LW

  
%   X: design matrix with four columns (future versions may allow different input configurations)
%       - first column  (i.e., X(:,1)) : all dependent variable values
%       - second column (i.e., X(:,2)) : between-subjects factor (e.g., subject group) level codes (ranging from 1:L where 
%           L is the # of levels for the between-subjects factor)
%       - third column  (i.e., X(:,3)) : within-subjects factor (e.g., condition/task) level codes (ranging from 1:L where 
%           L is the # of levels for the within-subjects factor)
%       - fourth column (i.e., X(:,4)) : subject codes (ranging from 1:N where N is the total number of subjects)
%   

xx=[mean(gp.meanpitchpertresp_precereb_EW')...
mean(gp.meanpitchpertresp_precereb_LW') ...
mean(gp.meanpitchpertresp_postcereb_EW')...
mean(gp.meanpitchpertresp_postcereb_LW') ...
mean(gp.meanpitchpertresp_presham_EW') ...
mean(gp.meanpitchpertresp_presham_LW') ...
mean(gp.meanpitchpertresp_postsham_EW') ...
mean(gp.meanpitchpertresp_postsham_LW') ...
mean(gp.meanpitchpertresp_prebeh_EW') ...
mean(gp.meanpitchpertresp_prebeh_LW') ...
mean(gp.meanpitchpertresp_postbeh_EW') ...
mean(gp.meanpitchpertresp_postbeh_LW')]

% DV data
mixedanova(:,1)=xx' 

% between subjects SITE
mixedanova(1:(13*4),2)=1; 
mixedanova(53:104,2)=2
mixedanova(105:144,2)=3

% within subject PRE POST
test=ones(1,26)
test2=test*2
test3=ones(1,20)
test4=test3*2
test5=[test test2 test test2 test3 test4]
mixedanova(:,3)=test5; 

test=1:13
test2=[14:26]
test3=27:36
mixedanova(:,4)=[test test test test test2 test2 test2 test2 test3 test3 test3 test3]
[SSQs, DFs, MSQs, Fs, Ps]=mixed_between_within_anova(mixedanova,1)

