%% pitch pert cerebellar cTBS data analysis 2015
% Notes on pert.resp
% gp.precereb.patient_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{1}
% pitch_in.dat{1} - response to down pert
% pitch_in.dat{2} response to up pert
% pitch_in.dat{3} respone to all but *it's not flipped*
% to flip it use: 
% centsdev_dat.(group{igroup}).subj(isubj).absdat = {-the_dat{1} the_dat{2} [-the_dat{1}; the_dat{2}]};


TODO
plot peaks 




% anova on mean comp values
plot_mean_comp_anova_cTBS


lab code is saved locally in 
/home/zagnew/matlab/runpitchpertanalysis.m
needs vline

data is in 
/data/bil-mb4/zarinah-data/cerebellar-data/pitch-pert-cTBS/cTBS_data/data_analysis


%my code:
/home/zagnew/data_analysis_code/
/home/zagnew/data_analysis_code/pitch_pert_stats/edit plot_mean_comp.m
/home/zagnew/data_analysis_code/pitch_pert_stats/pitchpert_fourway.m

%% now need to be here
cd /data/bil-mb4/zarinah-data/cerebellar-data/pitch-pert-cTBS/cTBS_data/data_analysis/cerebellarTBS
% needs an excell file here with all the data in:
/data/bil-mb4/zarinah-data/cerebellar-data/pitch-pert-cTBS/cTBS_data/data_analysis/cerebellarTBS/cTBS_info_files
post_cTBS.xls

[summary data is all in cTBSinfo_2016.xls]

cd /data/bil-mb4/zarinah-data/cerebellar-data/pitch-pert-cTBS/cTBS_data/data_analysis/cerebellarTBS
% run
get_data
%creates 
/data/bil-mb4/zarinah-data/cerebellar-data/pitch-pert-cTBS/cTBS_data/data_analysis/
patient_dat.mat and control_dat.mat


put all data together
/data/bil-mb4/zarinah-data/cerebellar-data/pitch-pert-cTBS/cTBS_data/data_analysis/patient.m

%% rename patient and control files from patient to patient_post.m etc


pre_cereb_data=load('patient_pre.mat')
pre_sham_data=load('control_pre.mat')
save pre_cereb_data pre_cereb_data
save pre_sham_data pre_sham_data

post_cereb_data=load('patient_post.mat')
post_sham_data=load('control_post.mat')
save post_cereb_data post_cereb_data
save post_sham_data post_sham_data




% now do get_data
load /data/bil-mb4/zarinah-data/cerebellar-data/pitch-pert-cTBS/cTBS_data/data_analysis/pre_cereb_data

pre_cereb_data.patient_dat=load('patient_pre')
pre_sham_data.patient_dat=load('control_pre')

save pre_cereb_data pre_cereb_data
save pre_sham_data pre_sham_data

%% for group analysis
run pitchpert_fourway.m
