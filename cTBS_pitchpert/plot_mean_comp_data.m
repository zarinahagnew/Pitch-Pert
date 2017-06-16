%% plots mean subject data

figure
subplot(221);
plot(gp.meanpitchpertresp_precereb')
axis([ 0 500 -20 60])
title('pre cerebellar')
goodplot

subplot(222);
plot(gp.meanpitchpertresp_postcereb')
axis([ 0 500 -20 60])
title('post cerebellar')
goodplot

subplot(223);
plot(gp.meanpitchpertresp_presham')
axis([ 0 500 -20 60])
title('pre sham')
goodplot

subplot(224);
plot(gp.meanpitchpertresp_postsham')
axis([ 0 500 -20 60])
title('post sham')
goodplot
print(gcf, '-dpdf', '-r150', '/Users/zagnew/cerebellarTBS/data_analysis/figures/subject_mean_comp.pdf');





subplot(224);
plot(gp.meanpitchpertresp_postsham')


gp.meanpitchpertresp_pre_all=[gp.meanpitchpertresp_precereb; gp.meanpitchpertresp_presham]


figure
subplot(311);
plot(gp.meanpitchpertresp_pre_all')
axis([ 0 500 -50 100])
title('all pre stim')
goodplot

subplot(312);
plot(gp.meanpitchpertresp_postcereb')
axis([ 0 500 -50 100])
title('post cereb')
goodplot

subplot(313);
plot(gp.meanpitchpertresp_postsham')
axis([ 0 500 -50 100])
title('post sham')
goodplot


x=1:size(nanmean(gp.meanpitchpertresp_precereb),2);
shadedErrorBar(x,mean(gp.meanpitchpertresp_postcereb), (std(gp.meanpitchpertresp_postcereb)/sqrt(13)),{'-','LineWidth', 1.5,'color',[0.8 0 0.2]}, 0.2);

figure
subplot(311);
plot(mean(gp.meanpitchpertresp_pre_all))
plot(mean(gp.meanpitchpertresp_pre_all))




