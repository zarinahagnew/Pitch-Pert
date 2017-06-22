
cd /Users/zagnew/Cereb_data/predpitchpert;
set(0,'DefaultFigureWindowStyle','docked');

%% unpred data
% scp zagnew@carvaka.radiology.ucsf.edu:/data/bil-mb4/zarinah-data/cerebellar-data/pitch-pert-ataxia/controls/control.mat .


% load full_pert_resps
%  full_pert_resps.patients.subj.pert_resp
% full_pert_resps.patients.subj(1,1).pert_resp.pert_types
% full_pert_resps.patients.subj(1,1).pert_resp.peakcomp{1,1}
% cents4comp(1).abspitch_in.dat
% frame_taxis = full_pert_resps.patients.subj.pert_resp.frame_taxis;
% tlims4anal_req = [-0.2 1.0];
% ilims4anal = dsearchn(frame_taxis',tlims4anal_req')';
% idxes4anal = ilims4anal(1):ilims4anal(2);
% tlims4anal = frame_taxis(ilims4anal);

numpats=9
numhcs=10
load('centsdev_dat.mat')
for isubj=1:numpats 
patients(isubj,:)=mean(centsdev_dat.patients.subj(1,isubj).absdat{1,1})
end
for isubj=1:numhcs 
controls(isubj,:)=mean(centsdev_dat.controls.subj(1,isubj).absdat{1,1})
end

figure
subplot(121)
plot(patients')
axis([0 400 -100 200])
goodplot
title('patients')

subplot(122)
plot(controls')
axis([0 400 -100 200])
goodplot
title('controls')

figure
plot(mean(patients), 'r')
hold
plot(mean(controls))




