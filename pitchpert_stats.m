pert_resp.mean_signed_response{1}
pert_resp.mean_signed_response{2}

sub04_pre_negative=pert_resp.mean_signed_response{1}
sub04_pre_positive=pert_resp.mean_signed_response{2}

% is there a difference between the up and down trials
ttest(sub04_pre_negative, sub04_pre_positive)

%make all values positive
sub04_pre_negative=abs(sub04_pre_negative) 

% append them
sub04_pre=[sub04_pre_negative ; sub04_pre_positive] 

meancomp_sun04=mean(sub04_pre)
stdevcomp_sun04=std(sub04_pre)
