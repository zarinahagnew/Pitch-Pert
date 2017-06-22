clear all
close all
cd /Users/zagnew/Cereb_data
set(0,'DefaultFigureWindowStyle','docked');


%% bring in the predictable data
% full_pert_resps.patients.subj(1,1).pert_resp.cents4comp.pitch_in.dat{1,1}
% flip this to match the unpredicable

load /Users/zagnew/Cereb_data/predpitchpert/full_pert_resps.mat;
npats_p=full_pert_resps.patients.nsubj;
nhcs_p=full_pert_resps.controls.nsubj;

for isubj=1:npats_p 
compdata.patients_pred(isubj,:)=mean(-full_pert_resps.patients.subj(1,isubj).pert_resp.cents4comp.pitch_in.dat{1,1}); 
compdata.patients_pred_sem(isubj,:)=std(-full_pert_resps.patients.subj(1,isubj).pert_resp.cents4comp.pitch_in.dat{1,1})/sqrt(npats_p);
end

% Sanity check:- 
% figure
% for isubj=1:9
% plot(mean(full_pert_resps.patients.subj(1,isubj).pert_resp.cents4comp.pitch_in.dat{1,1}))
% hold on
% plot(mean(compdata.patients_pred),'m')
% pause
% end
% 

for isubj=1:nhcs_p
compdata.controls_pred(isubj,:)=mean(-full_pert_resps.controls.subj(1,isubj).pert_resp.cents4comp.pitch_in.dat{1,1}); 
compdata.controls_pred_sem(isubj,:)=std(-full_pert_resps.controls.subj(1,isubj).pert_resp.cents4comp.pitch_in.dat{1,1})/sqrt(nhcs_p);
end

%% bring in the unpredictable data
load /Users/zagnew/Cereb_data/unpredpitchpert/patient;
load /Users/zagnew/Cereb_data/unpredpitchpert/control.mat;
npats_up=13;
npats_up=11;
for isubj=1:npats_up   
unpred_data_pats{isubj}=[-patient_dat.pert_resp(1,isubj).cents4comp.pitch_in.dat{1,1};patient_dat.pert_resp(1,isubj).cents4comp.pitch_in.dat{1,2}];
compdata.patients_unpred(isubj,:)=mean(unpred_data_pats{isubj});
compdata.patients_unpred_sem(isubj,:)=std(unpred_data_pats{isubj})/sqrt(npats_up);
end

for isubj=1:npats_up   
unpred_data_hcs{isubj}=[-control_dat.pert_resp(1,isubj).cents4comp.pitch_in.dat{1,1};control_dat.pert_resp(1,isubj).cents4comp.pitch_in.dat{1,2}];
compdata.controls_unpred(isubj,:)=mean(unpred_data_hcs{isubj});
compdata.controls_unpred_sem(isubj,:)=std(unpred_data_hcs{isubj})/sqrt(npats_up);
end

%% calc group means and sems
compdata.patients_pred_mean=(mean(compdata.patients_pred))
compdata.controls_pred_mean=(mean(compdata.controls_pred))
compdata.patients_unpred_mean=(mean(compdata.patients_unpred))
compdata.controls_unpred_mean=(mean(compdata.controls_unpred))

compdata.patients_pred_sem=(std(compdata.patients_pred)/sqrt(npats_p))
compdata.controls_pred_sem=(std(compdata.controls_pred)/sqrt(nhcs_p))
compdata.patients_unpred_sem=(std(compdata.patients_unpred)/sqrt(npats_up))
compdata.controls_unpred_sem=(std(compdata.controls_unpred)/sqrt(nhcs_p))

%% now we have
% compdata.patients_unpred
% compdata.controls_unpred
% compdata.patients_pred
% compdata.controls_pred

figure
subplot(211)
plot(compdata.patients_unpred_mean, 'r')
hold
plot(compdata.controls_unpred_mean)
axis([0 400 -30 40])
goodplot
title('unpredictable')

subplot(212)
plot(compdata.patients_pred_mean, 'r')
hold
plot(compdata.controls_pred_mean)
axis([0 400 -30 40])
goodplot
title('predictable')

figure
subplot(211)
x=1:size(compdata.patients_pred_mean,2);
shadedErrorBar(x,compdata.patients_pred_mean, compdata.patients_pred_sem,{'-','LineWidth', 1.5,'color',[0.8 0 0.2]});
hold 
shadedErrorBar(x,compdata.patients_unpred_mean, compdata.patients_unpred_sem,{'-','LineWidth', 1.5,'color',[0.1 0.1 0.4]});
axis([0 400 -30 40])
goodplot
title('patients')

subplot(212)
shadedErrorBar(x,compdata.controls_pred_mean, compdata.controls_pred_sem,{'-','LineWidth', 1.5,'color',[0.8 0 0.2]});
hold 
shadedErrorBar(x,compdata.controls_unpred_mean, compdata.controls_unpred_sem,{'-','LineWidth', 1.5,'color',[0.1 0.1 0.4]});
axis([0 400 -30 40])
goodplot
title('controls')
legend('','','','pred','','','' ,'unpred')


shadedErrorBar(x,compdata.controls_pred_mean, compdata.controls_pred_sem,{'-','LineWidth', 1.5,'color',[0.1 0.1 0.4]});
goodplot







% figure
% subplot(121)
% 
% x=1:size(dafdata.patient_utterance_mean,2);
% shadedErrorBar(x,dafdata.patient_utterance_mean, dafdata.patient_utterance_sem,{'-','LineWidth', 1.5,'color',[0.1 0.1 0.4]});
% hold on
% shadedErrorBar(x,dafdata.control_utterance_mean, dafdata.control_utterance_sem,{'-','LineWidth', 1.5,'color',[0.8 0 0.2]});
% axis([1 6 -0.1 0.7])
% goodplot_wide





% 
% 
% moo=[-patient_dat.pert_resp(1,isubj).cents4comp.pitch_in.dat{1} ; patient_dat.pert_resp(1,isubj).cents4comp.pitch_in.dat{2}];    
% moo2=patient_dat.pert_resp(1,isubj).cents4comp.pitch_in.dat{1,1}
% 
% figure
% subplot(311)
% plot(mean(patient_dat.pert_resp(1,isubj).cents4comp.pitch_in.dat{1,1}))
% subplot(312)
% plot(mean(patient_dat.pert_resp(1,isubj).cents4comp.pitch_in.dat{1,2}))
% 
% subplot(313)
% plot(mean(patient_dat.pert_resp(1,isubj).cents4comp.pitch_in.dat{1,3}))
