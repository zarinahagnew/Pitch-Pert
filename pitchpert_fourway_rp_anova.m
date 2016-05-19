
clear all
close all
cd /Users/zagnew/Cereb_data/cTBS/data_analysis
load /Users/zagnew/Cereb_data/cTBS/data_analysis/GROUPDATA.mat

% individual data
gp.precereb.patient_dat.pert_resp(1).cents4comp(1).abspitch_in.dat

% mean data
gp.meanpitchpertresp_precereb(1,:)
gp.meanpitchpertresp_postcereb(1,:)

gp.meanpitchpertresp_postcereb(1,:)-gp.meanpitchpertresp_precereb(1,:)

% mean pre post difference
gp.cereb_diff_eachsub(1,:)

cd data-for-r
for isubj=1:13
filename = [ 'cereb_prt_diff_' num2str(isubj) '.txt' ];
dlmwrite(filename,gp.cereb_diff_eachsub(isubj,:))
filename_sham = [ 'sham_prt_diff_' num2str(isubj) '.txt' ];
dlmwrite(filename,gp.sham_diff_eachsub(isubj,:))
end
