%% loads in the mean compensation from pre and post sessions, takes means, plots and does a ttest
% ZKA Feb 2015

clear all
close all

workingdir='/Users/zagnew/data_analysis/';
cd(workingdir);

%subjects
nsub=0;
nsub=nsub+1;
subinfo{nsub}.expdir='cerebellarTBS/s01_pre_cereb/';

subdata_pre_cereb(1)=load('/Users/zagnew/data_analysis/cerebellarTBS/s01_pre_cereb/speak_consolidate_audiodir/pert_resp.mat');
subdata_pre_cereb(2)=load('/Users/zagnew/data_analysis/cerebellarTBS/s03_pre_cereb/speak_consolidate_audiodir/pert_resp.mat');
subdata_pre_cereb(3)=load('/Users/zagnew/data_analysis/cerebellarTBS/s04_pre_cereb/speak_consolidate_audiodir/pert_resp.mat');
subdata_pre_cereb(4)=load('/Users/zagnew/data_analysis/cerebellarTBS/s05_pre_cereb/speak_consolidate_audiodir/pert_resp.mat');
subdata_pre_cereb(5)=load('/Users/zagnew/data_analysis/cerebellarTBS/s06_pre_cereb/speak_consolidate_audiodir/pert_resp.mat');
subdata_pre_cereb(6)=load('/Users/zagnew/data_analysis/cerebellarTBS/s07_pre_cereb/speak_consolidate_audiodir/pert_resp.mat');
subdata_pre_cereb(7)=load('/Users/zagnew/data_analysis/cerebellarTBS/s08_pre_cereb/speak_consolidate_audiodir/pert_resp.mat');
subdata_pre_cereb(8)=load('/Users/zagnew/data_analysis/cerebellarTBS/s10_pre_cereb/speak_consolidate_audiodir/pert_resp.mat');

subdata_post_cereb(1)=load('/Users/zagnew/data_analysis/cerebellarTBS/s01_post_cereb/speak_consolidate_audiodir/pert_resp.mat');
subdata_post_cereb(2)=load('/Users/zagnew/data_analysis/cerebellarTBS/s03_post_cereb/speak_consolidate_audiodir/pert_resp.mat');
subdata_post_cereb(3)=load('/Users/zagnew/data_analysis/cerebellarTBS/s04_post_cereb/speak_consolidate_audiodir/pert_resp.mat');
subdata_post_cereb(4)=load('/Users/zagnew/data_analysis/cerebellarTBS/s05_post_cereb/speak_consolidate_audiodir/pert_resp.mat');
subdata_post_cereb(5)=load('/Users/zagnew/data_analysis/cerebellarTBS/s06_post_cereb/speak_consolidate_audiodir/pert_resp.mat');
subdata_post_cereb(6)=load('/Users/zagnew/data_analysis/cerebellarTBS/s07_post_cereb/speak_consolidate_audiodir/pert_resp.mat');
subdata_post_cereb(7)=load('/Users/zagnew/data_analysis/cerebellarTBS/s08_post_cereb/speak_consolidate_audiodir/pert_resp.mat');
subdata_post_cereb(8)=load('/Users/zagnew/data_analysis/cerebellarTBS/s10_post_cereb/speak_consolidate_audiodir/pert_resp.mat');

for isub=1:8
data_pre_cereb{isub}=vertcat(subdata_pre_cereb(isub).pert_resp.comp{1}, subdata_pre_cereb(isub).pert_resp.comp{2})';
data_post_cereb{isub}=vertcat(subdata_post_cereb(isub).pert_resp.comp{1}, subdata_post_cereb(isub).pert_resp.comp{2})';
subject_mean_pre_cereb(isub)=mean(data_pre_cereb{isub});
subject_mean_post_cereb(isub)=mean(data_post_cereb{isub});
subject_std_pre_cereb(isub)=std(data_pre_cereb{isub});
subject_std_post_cereb(isub)=std(data_post_cereb{isub});
singlesub_ttest(isub)=ttest2(data_pre_cereb{isub},data_post_cereb{isub});
end

% plot indivdual means
for isub=1:7
figure
subplot(1,2,1)
cond=[subject_mean_pre_cereb(isub);subject_mean_post_cereb(isub)];
errY2=[subject_std_pre_cereb(isub)/sqrt(length(data_pre_cereb{isub}));subject_std_pre_cereb(isub)/sqrt(length(data_post_cereb{isub}))];

h = barwitherr(errY2, cond);% Plot with errorbars
set(gca,'XTickLabel',{'Pre Stim','Post Stim'})
ylabel('mean compensation (cents)')
set(h(1),'FaceColor','w');
title(sprintf('Compensation to PP before and after cerebellar TBS'));
goodplot
axis([0 3, -1 1])

if singlesub_ttest(isub)==0
   stats='insignificant';
else
   stats='significant';
end

text(1, 1, stats)

% scatter plot
subplot(1,2,2)
scatter(ones(1,length(data_pre_cereb{isub})),data_pre_cereb{isub},'.','k')
hold
scatter(1, mean(data_pre_cereb{isub}(1:length(data_pre_cereb{isub}))), 'filled', 'm')
scatter(2*(ones(1,length(data_post_cereb{isub}))), data_post_cereb{isub}','.','k')
scatter(2, mean(data_post_cereb{isub}(1:length(data_post_cereb{isub}))), 'filled', 'm')
axis([0 3, -1 1])
goodplot

cd figures
filename=sprintf('singlesub_cereb%d.pdf',isub);
%saveas(gcf,filename,'png')
print(gcf, '-dpdf', '-r150', filename);
cd ..
end



%% sham data
subdata_pre_sham(1)=load('/Users/zagnew/data_analysis/shamTBS/s03_pre_sham/speak_consolidate_audiodir/pert_resp.mat');
subdata_pre_sham(2)=load('/Users/zagnew/data_analysis/shamTBS/s05_pre_sham/speak_consolidate_audiodir/pert_resp.mat');
subdata_pre_sham(3)=load('/Users/zagnew/data_analysis/shamTBS/s06_pre_sham/speak_consolidate_audiodir/pert_resp.mat');
subdata_pre_sham(4)=load('/Users/zagnew/data_analysis/shamTBS/s07_pre_sham/speak_consolidate_audiodir/pert_resp.mat');
subdata_pre_sham(5)=load('/Users/zagnew/data_analysis/shamTBS/s09_pre_sham/speak_consolidate_audiodir/pert_resp.mat');
subdata_pre_sham(6)=load('/Users/zagnew/data_analysis/shamTBS/s10_pre_sham/speak_consolidate_audiodir/pert_resp.mat');
subdata_pre_sham(7)=load('/Users/zagnew/data_analysis/shamTBS/s11_pre_sham/speak_consolidate_audiodir/pert_resp.mat');

subdata_post_sham(1)=load('/Users/zagnew/data_analysis/shamTBS/s03_post_sham/speak_consolidate_audiodir/pert_resp.mat');
subdata_post_sham(2)=load('/Users/zagnew/data_analysis/shamTBS/s05_post_sham/speak_consolidate_audiodir/pert_resp.mat');
subdata_post_sham(3)=load('/Users/zagnew/data_analysis/shamTBS/s06_post_sham/speak_consolidate_audiodir/pert_resp.mat');
subdata_post_sham(4)=load('/Users/zagnew/data_analysis/shamTBS/s07_post_sham/speak_consolidate_audiodir/pert_resp.mat');
subdata_post_sham(5)=load('/Users/zagnew/data_analysis/shamTBS/s09_post_sham/speak_consolidate_audiodir/pert_resp.mat');
subdata_post_sham(6)=load('/Users/zagnew/data_analysis/shamTBS/s10_post_sham/speak_consolidate_audiodir/pert_resp.mat');
subdata_post_sham(7)=load('/Users/zagnew/data_analysis/shamTBS/s11_post_sham/speak_consolidate_audiodir/pert_resp.mat');


for isub=1:7
data_pre_sham{isub}=vertcat(subdata_pre_sham(isub).pert_resp.comp{1}, subdata_pre_sham(isub).pert_resp.comp{2})';
data_post_sham{isub}=vertcat(subdata_post_sham(isub).pert_resp.comp{1}, subdata_post_sham(isub).pert_resp.comp{2})';
subject_mean_pre_sham(isub)=mean(data_pre_sham{isub});
subject_mean_post_sham(isub)=mean(data_post_sham{isub});
subject_std_pre_sham(isub)=std(data_pre_sham{isub});
subject_std_post_sham(isub)=std(data_post_sham{isub});
singlesub_ttest_sham(isub)=ttest2(data_pre_sham{isub},data_post_sham{isub});
end

subject_mean_pre_sham(isub)=mean(data_pre_sham{isub});
subject_mean_post_sham(isub)=mean(data_post_sham{isub});
subject_std_pre_sham(isub)=std(data_pre_sham{isub});
subject_std_post_sham(isub)=std(data_post_sham{isub});

% sham figure
figure
subplot(1,2,1)
cond=[subject_mean_pre_sham(isub);subject_mean_post_cereb(isub)];
errY2=[subject_std_pre_cereb(isub)/sqrt(length(data_pre_cereb{isub}));subject_std_pre_cereb(isub)/sqrt(length(data_post_cereb{isub}))];

h = barwitherr(errY2, cond);% Plot with errorbars
set(gca,'XTickLabel',{'Pre Stim','Post Stim'})
ylabel('mean compensation (cents)')
set(h(1),'FaceColor','w');
title(sprintf('Compensation to PP before and after cerebellar TBS'));
goodplot
axis([0 3, -1 1])

if singlesub_ttest_sham(isub)==0
   stats='insignificant';
else
   stats='significant';
end

text(1, 1, stats)

% scatter plot
subplot(1,2,2)
scatter(ones(1,length(data_pre_cereb{isub})),data_pre_cereb{isub},'.','k')
hold
scatter(1, mean(data_pre_cereb{isub}(1:length(data_pre_cereb{isub}))), 'filled', 'm')
scatter(2*(ones(1,length(data_post_cereb{isub}))), data_post_cereb{isub}','.','k')
scatter(2, mean(data_post_cereb{isub}(1:length(data_post_cereb{isub}))), 'filled', 'm')
axis([0 3, -1 1])
goodplot

cd figures
filename=sprintf('singlesub_sham%d.png',isub);
%saveas(gcf,filename,'png')
print(gcf, '-dpdf', '-r150', filename);
cd ..


%% single subject anova

anovasub=isub;
anovadata= [data_pre_sham{anovasub} data_post_sham{anovasub} data_pre_cereb{anovasub} data_post_cereb{anovasub}];

a=length(data_pre_sham{anovasub});
b=length(data_post_sham{anovasub});
c=length(data_pre_cereb{anovasub});
d=length(data_post_cereb{anovasub});
e=length(anovadata);

test=ones(1,a+b);
test2=ones(1,c+d);
test2=test2*2;
site=[test test2];

%create conditions
condition= cell(1,e);
for i=1:a
    condition{i} = 'pre'; 
end
for i=a+1:a+b
    condition{i} = 'post';
end
for i=a+b+1:a+b+c
    condition{i} = 'pre';
end
for i=a+b+c+1:e
    condition{i} = 'post';
end

condition=condition';

group1=[site];
group2=[condition];
p = anovan(anovadata,{group1 group2 },'model','interaction')

%plot anova

cond_sham=[mean(data_pre_sham{anovasub}(1:length(data_pre_sham{anovasub})));mean(data_post_sham{anovasub})];
errY2_sham=[std(data_pre_sham{anovasub}(1:length(data_pre_sham{anovasub}))/sqrt(length(data_pre_sham{anovasub})));std(data_post_sham{anovasub})/sqrt(length(data_pre_sham{anovasub}))];
cond_cereb=[mean(data_pre_cereb{anovasub}(1:length(data_pre_cereb{anovasub})));mean(data_post_cereb{anovasub})];
errY2_cereb=[std(data_pre_cereb{anovasub}(1:length(data_pre_cereb{anovasub}))/sqrt(length(data_pre_cereb{anovasub})));std(data_post_cereb{anovasub})/sqrt(length(data_post_cereb{anovasub}))];

figure
subplot(221)
h = barwitherr(errY2_sham, cond_sham);% Plot with errorbars
set(gca,'XTickLabel',{'Pre Stim','Post Stim'})
ylabel('mean compensation (cents)')
set(h(1),'FaceColor','w');
title(sprintf('                                      Compensation to PP before and after TBS'));
goodplot

% scatter plot
subplot(222)
scatter(ones(1,length(data_pre_sham{anovasub})),data_pre_sham{anovasub},'.','k')
hold
scatter(1, mean(data_pre_sham{anovasub}(1:length(data_pre_sham{anovasub}))), 'filled', 'm')
scatter(2*(ones(1,length(data_post_sham{anovasub}))),data_post_sham{anovasub},'.','k')
scatter(2, mean(data_post_sham{anovasub}(1:length(data_post_sham{anovasub}))), 'filled', 'm')
axis([0 3, -1 1])
goodplot

subplot(223)
h = barwitherr(errY2_cereb, cond_cereb);% Plot with errorbars
set(gca,'XTickLabel',{'Pre Stim','Post Stim'})
ylabel('mean compensation (cents)')
set(h(1),'FaceColor','w');
goodplot

% scatter plot
subplot(224)
scatter(ones(1,length(data_pre_cereb{anovasub})),data_pre_cereb{anovasub},'.','k')
hold
scatter(1, mean(data_pre_cereb{anovasub}(1:length(data_pre_cereb{anovasub}))), 'filled', 'm')
scatter(2*(ones(1,length(data_post_cereb{anovasub}))),data_post_cereb{anovasub},'.','k')
scatter(2, mean(data_post_cereb{anovasub}(1:length(data_post_cereb{anovasub}))), 'filled', 'm')
axis([0 3, -1 1])
goodplot

filename=sprintf('singlesub_anova%d.png',isub);
print(gcf, '-dpdf', '-r150', filename);

%end



%% group anova
gp_anovadata= [subject_mean_pre_cereb subject_mean_post_cereb subject_mean_pre_sham subject_mean_post_sham];

test=ones(1,length(subject_mean_pre_cereb)*2);
test2=ones(1,length(subject_mean_pre_sham)*2);
test2=test2*2;
site=[test test2];

%create conditions
clear condition
for i=1:length(subject_mean_pre_cereb)
    condition1{i} = 'pre'; 
end
for i=1:length(subject_mean_pre_cereb)
    condition2{i} = 'post'; 
end

for i=1:length(subject_mean_pre_sham)
    condition3{i} = 'pre'; 
end

for i=1:length(subject_mean_post_sham)
    condition4{i} = 'post'; 
end

condition=[condition1 condition2 condition3 condition4];
condition=condition';

group1=[site];
group2=[condition];
p = anovan(gp_anovadata,{group1 group2 },'model','interaction')

% plot group anova data
cond_cereb=[mean(subject_mean_pre_cereb) mean(subject_mean_post_cereb)];
cond_sham=[mean(subject_mean_pre_sham) mean(subject_mean_post_sham)];

errY2_cereb=[std(subject_mean_pre_cereb)/sqrt(length(subject_mean_pre_cereb)) std(subject_mean_post_cereb)/sqrt(length(subject_mean_post_cereb))]
errY2_sham=[std(subject_mean_pre_sham)/sqrt(length(subject_mean_pre_sham)) std(subject_mean_post_sham)/sqrt(length(subject_mean_post_sham))]    

figure
subplot(221)
h = barwitherr(errY2_sham, cond_sham);% Plot with errorbars
set(gca,'XTickLabel',{'Pre Sham','Post Sham'})
ylabel('mean compensation (cents)')
set(h(1),'FaceColor','w');
title(sprintf('                                      Compensation to PP before and after sham TBS'));
axis([0 3, 0 0.4])
goodplot

% scatter plot
subplot(222)
scatter(ones(1,8),subject_mean_pre_cereb,'.','k')
hold 
scatter(1, mean(subject_mean_pre_cereb), 'filled', 'm')
scatter(2*(ones(1,8)),subject_mean_post_cereb,'.','k')
scatter(2, mean(subject_mean_post_cereb), 'filled', 'm')
axis([0 3, -0.2 0.5])
goodplot

subplot(223)
h = barwitherr(errY2_cereb, cond_cereb);% Plot with errorbars
set(gca,'XTickLabel',{'Pre cereb','Post cereb'})
ylabel('mean compensation (cents)')
set(h(1),'FaceColor','w');
title(sprintf('cerebellar TBS'));
axis([0 3, 0 0.4])

goodplot

ttest(subject_mean_pre_cereb, subject_mean_post_cereb)
ttest2(subject_mean_post_sham, subject_mean_post_cereb)
% scatter plot
subplot(224)
scatter(ones(1,7),subject_mean_pre_sham,'.','k')
hold 
scatter(1, mean(subject_mean_pre_sham), 'filled', 'm')
scatter(2*(ones(1,7)),subject_mean_post_sham,'.','k')
scatter(2, mean(subject_mean_post_cereb), 'filled', 'm')
axis([0 3, -0.2 0.5])
goodplot

filename=sprintf('singlesub_anova%d.png',isub);
print(gcf, '-dpdf', '-r150', filename);

