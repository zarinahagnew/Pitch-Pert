
%% 1. flip pitch and concatenate pitch_in.dat{1} and {2}
for isubj=1:numsubs_1 
    %cerebellar    
    gp.precereb.patient_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat = [-gp.precereb.patient_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{1} ; gp.precereb.patient_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{2}];    
    gp.postcereb.control_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat = [-gp.postcereb.control_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{1} ; gp.postcereb.control_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{2}];   
    %vertex
    gp.presham.patient_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat = [-gp.presham.patient_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{1} ; gp.presham.patient_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{2}];
    gp.postsham.control_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat = [-gp.postsham.control_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{1} ; gp.postsham.control_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{2}];
end
for isubj=1:numsubs_2
 % beh only
    gp.prebeh.patient_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat = [-gp.prebeh.patient_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{1} ; gp.prebeh.patient_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{2}];
    gp.postbeh.control_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat = [-gp.postbeh.control_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{1}; gp.postbeh.control_dat.pert_resp(isubj).cents4comp(1).pitch_in.dat{2}];
end

% to check the up and down trials are all good and correct:
% calc mean compensation, and the difference between post and pre compensation for each subject

%% 2. calculate mean comp and difference in compensation
for isubj=1:numsubs_1;
    %mean comp
    gp.meanpitchpertresp_precereb(isubj,:)=mean(gp.precereb.patient_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat);
    gp.meanpitchpertresp_postcereb(isubj,:)=mean(gp.postcereb.control_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat);
    gp.meanpitchpertresp_presham(isubj,:)=mean(gp.presham.patient_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat);
    gp.meanpitchpertresp_postsham(isubj,:)=mean(gp.postsham.control_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat);
    
    %std comp
    gp.stdpitchpertresp_precereb(isubj,:)=std(gp.precereb.patient_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat);
    gp.stdpitchpertresp_postcereb(isubj,:)=std(gp.postcereb.control_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat);
    gp.stdpitchpertresp_presham(isubj,:)=std(gp.presham.patient_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat);
    gp.stdpitchpertresp_postsham(isubj,:)=std(gp.postsham.control_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat);
    
    %std within trial comp
    gp.wt_stdpitchpertresp_precereb(isubj,:)=std(mean(gp.precereb.patient_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat));
    gp.wt_stdpitchpertresp_postcereb(isubj,:)=std(mean(gp.postcereb.control_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat));
    gp.wt_stdpitchpertresp_presham(isubj,:)=std(mean(gp.presham.patient_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat));
    gp.wt_stdpitchpertresp_postsham(isubj,:)=std(mean(gp.postsham.control_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat));
    
    %mean difference between post and pre
    gp.cereb_diff_eachsub(isubj,:)=gp.meanpitchpertresp_postcereb(isubj,:)-gp.meanpitchpertresp_precereb(isubj,:);
    gp.sham_diff_eachsub(isubj,:)=gp.meanpitchpertresp_postsham(isubj,:)-gp.meanpitchpertresp_presham(isubj,:);       
    
    gp.cereb_peakdiff_eachsub(isubj)=max(gp.cereb_diff_eachsub(isubj,:));
    gp.sham_peakdiff_eachsub(isubj)=max(gp.sham_diff_eachsub(isubj,:));
    
end

% beh
for isubj=1:numsubs_2;
    gp.meanpitchpertresp_prebeh(isubj,:)=mean(gp.prebeh.patient_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat);
    gp.meanpitchpertresp_postbeh(isubj,:)=mean(gp.postbeh.control_dat.pert_resp(isubj).cents4comp(1).abspitch_in.dat);  
    
    gp.beh_diff_eachsub(isubj,:)=gp.meanpitchpertresp_postbeh(isubj,:)-gp.meanpitchpertresp_prebeh(isubj,:);
    gp.beh_peakdiff_eachsub(isubj)=max(gp.beh_diff_eachsub(isubj,:));
end
