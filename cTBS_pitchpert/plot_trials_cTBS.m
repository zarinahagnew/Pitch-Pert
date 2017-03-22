

addpath('c:/matlab/myfiles')

figure
subplot(311)
plot(mean(gp.precereb.patient_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{1}), 'Color',clear_colour);
hold 
plot(mean(gp.precereb.patient_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{2}), 'm')
plot(mean(gp.precereb.patient_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat),'Color',masked_colour,'LineWidth',3);
axis([0 600 -60 60])
goodplot
title('cereb')

subplot(312)
plot(mean(gp.presham.patient_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{1}), 'Color',clear_colour);
hold 
plot(mean(gp.presham.patient_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{2}), 'm')
plot(mean(gp.presham.patient_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat),'Color',masked_colour,'LineWidth',3);
axis([0 600 -60 60])
goodplot
title('vertex')

subplot(313)
plot(mean(gp.prebeh.patient_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{1}), 'Color',clear_colour);
hold 
plot(mean(gp.prebeh.patient_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{2}), 'm')
plot(mean(gp.prebeh.patient_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat),'Color',masked_colour,'LineWidth',3);
axis([0 600 -60 60])
goodplot
title('no stim')
legend('response to down pert', 'response to up pert', 'flipped and averaged','Location','SouthWest')
