for j= 1:79
    temp_bs_score_smi_np=bs_score_smi_np(:,:,j);
    temp_bs_score_smi_np=temp_bs_score_smi_np';
    temp_bs_score_smi_np=temp_bs_score_smi_np(:);
    temp_bs_score_smi_np=reshape(temp_bs_score_smi_np(12:end),325,1);
    bs_score_smi_np_trans(:,j)=temp_bs_score_smi_np;
end
for i=1:79
    [acf(:,i),lags(:,i),bounds(:,i)] = autocorr(bs_score_smi_np_trans(:,i),12);
end
