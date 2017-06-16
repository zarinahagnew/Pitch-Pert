%% STDEV
% figure
% subplot(411);plot(gp.stdpitchpertresp_presham');goodplot
% subplot(412);plot(gp.stdpitchpertresp_postsham');goodplot
% subplot(413);plot(gp.stdpitchpertresp_precereb');goodplot
% subplot(414);plot(gp.stdpitchpertresp_postcereb');goodplot

figure
x=1:size(std(gp.stdpitchpertresp_presham),2);
subplot(411); shadedErrorBar(x,std(gp.stdpitchpertresp_presham), (std(gp.stdpitchpertresp_presham)/sqrt(13)),{'-','LineWidth', 1.5,'color',[0.8 0 0.2]}, 0.2); goodplot; title('stdev presham')
subplot(412); shadedErrorBar(x,std(gp.stdpitchpertresp_postsham), (std(gp.stdpitchpertresp_postsham)/sqrt(13)),{'-','LineWidth', 1.5,'color',[0.8 0 0.2]}, 0.2); goodplot;title('stdev postsham')
subplot(413); shadedErrorBar(x,std(gp.stdpitchpertresp_precereb), (std(gp.stdpitchpertresp_precereb)/sqrt(13)),{'-','LineWidth', 1.5,'color',[0.8 0 0.2]}, 0.2); goodplot;title('stdev precereb')
subplot(414); shadedErrorBar(x,std(gp.stdpitchpertresp_postcereb), (std(gp.stdpitchpertresp_postcereb)/sqrt(13)),{'-','LineWidth', 1.5,'color',[0.8 0 0.2]}, 0.2); goodplot;title('stdev postcereb')

b=60;
x=1:b,2;
figure
subplot(411); shadedErrorBar(x,std(gp.stdpitchpertresp_presham(:,1:b)), (std(gp.stdpitchpertresp_presham(:,1:b))/sqrt(13)),{'-','LineWidth', 1.5,'color',[0.8 0 0.2]}, 0.2); goodplot; title('stdev presham')
subplot(412); shadedErrorBar(x,std(gp.stdpitchpertresp_postsham(:,1:b)), (std(gp.stdpitchpertresp_postsham(:,1:b))/sqrt(13)),{'-','LineWidth', 1.5,'color',[0.2 0 0.2]}, 0.2); goodplot; title('stdev postsham')
subplot(413); shadedErrorBar(x,std(gp.stdpitchpertresp_precereb(:,1:b)), (std(gp.stdpitchpertresp_precereb(:,1:b))/sqrt(13)),{'-','LineWidth', 1.5,'color',[0.8 0 0.2]}, 0.2); goodplot; title('stdev precereb')
subplot(414); shadedErrorBar(x,std(gp.stdpitchpertresp_postcereb(:,1:b)), (std(gp.stdpitchpertresp_postcereb(:,1:b))/sqrt(13)),{'-','LineWidth', 1.5,'color',[0.2 0 0.2]}, 0.2); goodplot; title('stdev postcereb')

figure
x=1:size(std(gp.stdpitchpertresp_presham),2);
subplot(411); shadedErrorBar(x,std(gp.stdpitchpertresp_presham), (std(gp.stdpitchpertresp_presham)/sqrt(13)),{'-','LineWidth', 1.5,'color',[0.8 0 0.2]}, 0.2); goodplot; title('stdev presham')
subplot(412); shadedErrorBar(x,std(gp.stdpitchpertresp_postsham), (std(gp.stdpitchpertresp_postsham)/sqrt(13)),{'-','LineWidth', 1.5,'color',[0.8 0 0.2]}, 0.2); goodplot;title('stdev postsham')
subplot(413); shadedErrorBar(x,std(gp.stdpitchpertresp_precereb), (std(gp.stdpitchpertresp_precereb)/sqrt(13)),{'-','LineWidth', 1.5,'color',[0.8 0 0.2]}, 0.2); goodplot;title('stdev precereb')
subplot(414); shadedErrorBar(x,std(gp.stdpitchpertresp_postcereb), (std(gp.stdpitchpertresp_postcereb)/sqrt(13)),{'-','LineWidth', 1.5,'color',[0.8 0 0.2]}, 0.2); goodplot;title('stdev postcereb')


%% WT STDEV
figure
location=ones(1, length(gp.wt_stdpitchpertresp_precereb));
scatter(ones(1, length(gp.wt_stdpitchpertresp_precereb)),gp.wt_stdpitchpertresp_precereb)
hold
scatter(2*location,gp.wt_stdpitchpertresp_postcereb)
scatter(3*location,gp.wt_stdpitchpertresp_presham)
scatter(4*location,gp.wt_stdpitchpertresp_postsham)
c = 70
scatter(1, mean(gp.wt_stdpitchpertresp_precereb),c,'filled')
scatter(2, mean(gp.wt_stdpitchpertresp_postcereb),c,'filled')
scatter(3, mean(gp.wt_stdpitchpertresp_presham),c,'filled')
scatter(4, mean(gp.wt_stdpitchpertresp_postsham),c,'filled')
legend('pre cereb', 'post cereb', 'pre sham', 'post sham')
axis([0 5 0 30])
goodplot

pause
close all

