% extra plots

%% plot one subject
    fig1=figure
    isubj=1
    annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Pert Resps for Four Sessions', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')

    subplot(221)
    title('pre cereb stim')
    plot(gp.precereb.patient_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat', 'k')
    hold on
    plot(mean(gp.precereb.patient_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat), 'LineWidth', 1.3, 'Color','r');
    goodplot    
    
    subplot(222)
    title('post cereb stim')
    plot(gp.postcereb.patient_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat', 'k')
    hold on
    plot(mean(gp.postcereb.patient_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat), 'LineWidth', 1.3, 'Color','r');
    goodplot    
    
    subplot(223)
    title('pre sham stim')
    plot(gp.presham.control_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat', 'k')
    hold on
    plot(mean(gp.presham.control_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat), 'LineWidth', 1.3, 'Color','r');
    goodplot    
    
    subplot(224)
    title('post sham stim')
    plot(gp.postsham.control_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat', 'k')
    hold on
    plot(mean(gp.postsham.control_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat), 'LineWidth', 1.3, 'Color','r');
    goodplot
    
%% plot each subject
ymax=100;
ymin=-100;

for isubk=1:13
fig2=figure
subplot(211)
plot(frame_taxis,gp.meanpitchpertresp_precereb(isubj,:), 'r','LineWidth',2);
hold on    
plot(frame_taxis,gp.meanpitchpertresp_postcereb(isubj,:),'--r','LineWidth',2);   
goodplot
axis([0 1.2 ymin ymax])
title('pert resp - cerebellar')
legend('pre cereb pert resp', 'post cereb pert resp')

subplot(212)
plot(frame_taxis,gp.meanpitchpertresp_presham(isubj,:), 'k','LineWidth',2);
hold on    
plot(frame_taxis,gp.meanpitchpertresp_postsham(isubj,:),'--k','LineWidth',2);
goodplot
axis([0 1.2 ymin ymax])
title('pert resp - sham')
legend('pre sham pert resp', 'post sham pert resp')
cd
saveas(fig2, (['figures/SingleSub_pertresps_final_',num2str(isubj) '.jpg']))
end



%% 
figure
% plot individual subjects  - diff pre and post
ymax=100;
ymin=-100;
for isubj=1:13
    subplot(5, 3, isubj)
    h(1)=plot(frame_taxis(idxes4anal), gp.cereb_diff_eachsub(isubj,:), 'k','LineWidth',3);
    hold on
    h(2)=plot(frame_taxis(idxes4anal),gp.sham_diff_eachsub(isubj,:), 'r','LineWidth',3);
    xlabel('Time (s)')
    axis([-0.2 1 -100 50])
    goodplot
end
legend('cerebellar stimulation','sham stimulation', 'Location','northeast')

%%
ymax=50;
ymin=-50;
figure
subplot(211)
for isubj=1:13 %length(gp.precereb.patient_dat.pert_resp)
    plot(frame_taxis(idxes4anal), gp.cereb_diff_eachsub(isubj,:), 'k')
    hold on
    plot(frame_taxis(idxes4anal),mean(gp.cereb_diff_eachsub), 'r','LineWidth',3);
end
title(sprintf('cerebellar stimulation'));
xlabel('frames')
ylabel('post stim comp - pre stim comp')
axis([-0.2 1 ymin ymax])
goodplot

subplot(212)
for isubj=1:13 %length(gp.presham.control_dat.pert_resp)
    plot(frame_taxis(idxes4anal), gp.sham_diff_eachsub(isubj,:), 'k')
    hold on
    plot(frame_taxis(idxes4anal), mean(gp.sham_diff_eachsub), 'r','LineWidth',3);
end
axis([-0.2 1 ymin ymax])
xlabel('frames')
ylabel('post stim comp - pre stim comp')
title(sprintf('vertex stimulation'));
goodplot

%% Plot pre and post compensation for each subject
for isubj=1:13 % length(gp.precereb.patient_dat.pert_resp)
    figure
    subplot (2, 1, 1)
    plot(frame_taxis(idxes4anal),gp.meanpitchpertresp_precereb(isubj,:), 'k','LineWidth',3);
    hold on
    plot(frame_taxis(idxes4anal),gp.meanpitchpertresp_postcereb(isubj,:), 'R','LineWidth',3);
    axis([-0.2 1 -60 60])
    title(sprintf('cereb stimulation'));
    xlabel('Time(s)')
    ylabel('Mean Compensation (cents)')
    goodplot
    
    subplot (2, 1, 2)
    plot(frame_taxis(idxes4anal),gp.meanpitchpertresp_presham(isubj,:), 'k','LineWidth',3);
    hold on
    plot(frame_taxis(idxes4anal),gp.meanpitchpertresp_postsham(isubj,:), 'R','LineWidth',3);
    axis([-0.2 1 -60 60])
    title(sprintf('vertex stimulation'));
    
    xlabel('Time(s)')
    ylabel('Mean Compensation (cents)')
    saveas(gcf, sprintf('figure%d.jpg', isubj))
    goodplot
end

%% plot mean pitch pert in general
figure
subplot(311)
plot(frame_taxis(idxes4anal),nanmean(gp.meanpitchpertresp_precereb), 'k','LineWidth',3);
hold on
plot(frame_taxis(idxes4anal),nanmean(gp.meanpitchpertresp_postcereb),'r','LineWidth',3);
axis([-0.2 1 -60 60])
goodplot
xlabel('Time(s)')
ylabel('Mean Compensation (cents)')
title(sprintf('cerebellar stimulation'));
legend('prestim','post stim')

subplot(312)  
plot(frame_taxis(idxes4anal),nanmean(gp.meanpitchpertresp_presham), 'k','LineStyle','--','LineWidth',3);
hold on
plot(frame_taxis(idxes4anal),nanmean(gp.meanpitchpertresp_postsham),'r','LineStyle','--','LineWidth',3);
axis([-0.2 1 -60 60])
goodplot
xlabel('Time(s)')
ylabel('Mean Compensation (cents)')
legend('prestim','post stim')
title(sprintf('vertex stimulation'));

subplot(313)  
plot(frame_taxis(idxes4anal),nanmean(gp.meanpitchpertresp_prebeh), 'k','LineStyle','--','LineWidth',3);
hold on
plot(frame_taxis(idxes4anal),nanmean(gp.meanpitchpertresp_postbeh),'r','LineStyle','--','LineWidth',3);
axis([-0.2 1 -60 60])
goodplot
xlabel('Time(s)')
ylabel('Mean Compensation (cents)')
legend('prestim','post stim')
title(sprintf('no stimulation'));




%%
figure
plot(frame_taxis,-nanmean(gp.meanpitchpertresp_postcereb)','r','LineWidth',3);
hold 
plot(frame_taxis,-nanmean(gp.meanpitchpertresp_postsham)','Color',masked_colour,'LineWidth',3);
plot(frame_taxis,-nanmean(gp.meanpitchpertresp_postbeh)','Color',clear_colour,'LineWidth',3);
goodplot
legend('cerebellar stimulation', 'vertex', 'no stim', 'Location','SouthEast')
print(gcf, '-dpdf', '-r150', '/Users/zagnew/cerebellarTBS/data_analysis/figures/poststim responses.pdf');


%% meanpitchpertresp_postcereb

% plot post stim responses only
figure
subplot(211)
plot(gp.meanpitchpertresp_postcereb')
hold 
plot(nanmean(gp.meanpitchpertresp_postcereb), 'k','LineWidth',1.3);
goodplot
axis([0 450 -60 60])

subplot(212)
plot(gp.meanpitchpertresp_postsham')
hold 
plot(nanmean(gp.meanpitchpertresp_postsham), 'k','LineWidth',1.3);
goodplot
axis([0 450 -60 60])


