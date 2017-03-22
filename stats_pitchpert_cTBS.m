%% stats for mean comp

cd /Users/zagnew/cerebellarTBS
load pitch_pert_cTBS_data.mat;

% peak comp
pre_peakcomp_cereb=data.peakcomp(1:2:12);
post_peakcomp_cereb=data.peakcomp(2:3:12);

pre_peakcomp_sham=data.peakcomp(13:2:22);

patients_peakcomp=data.peakcomp(12:end);
controls_meanpeakcomp=mean(controls_peakcomp):
patients_meanpeakcomp=mean(patients_peakcomp);
controls_sempeakcomp=std(controls_peakcomp)/sqrt(length(controls_peakcomp));
patients_sempeakcomp=std(patients_peakcomp)/sqrt(length(patients_peakcomp));

[h,p,ci,stats] = ttest2(controls_peakcomp, patients_peakcomp)

% t peak
controls_tpeak=data.tpeak(1:11);
patients_tpeak=data.tpeak(12:end);
controls_meantpeak=mean(controls_tpeak);
patients_meantpeak=mean(patients_tpeak);
controls_semtpeak=std(controls_tpeak)/sqrt(length(controls_tpeak));
patients_semtpeak=std(patients_tpeak)/sqrt(length(patients_tpeak));

[h,p,ci,stats] = ttest2(controls_tpeak, patients_tpeak)








data.subject{:}

controls=data.subject{1:11}
patients=data.subject{12:end}