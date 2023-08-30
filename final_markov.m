clc;
clearvars -except  best_cop_indx best_cop_values index_cat_ci_bsn index_cat point_MK_v_bsn...
    point_sp_v_bsn point_adf_v_bsn point_kpss_v_bsn point_lmc_v_bsn point_pp_const_v_bsn...
    point_adf_const_v_bsn point_pp_v_bsn point_pp_tr_v_bsn index_cat_ci_bsn_vec index_cat_ci_ds_calib;
%% Reading Input Data File and Data Preparation
ind_xls=xlsread('index_basin.xlsx','ind','D2:F80');
ind_basin=ind_xls(:,1);
n_c=4;%can be diffrent
n_ci=3;%can be diffrent
nm=n_ci*n_c;
mm=size(index_cat,1);
nn=size(index_cat,2);
m=ceil(mm*(2/3));%217+108=325
l_q=(mm-m)+1;
%% Constants - Defining Variables In Advance to Make it Faster
e_train_drought_temp=zeros(n_c,l_q-1,nn);
e_train_temp=zeros(nm,l_q-1,nn);
e_train_drought_temp_bin=zeros(n_c,l_q-1,nn);
e_train_temp_bin=zeros(nm,l_q-1,nn);
RPS_train=zeros(nn,1);
BS_train=zeros(nn,1);
Rel_train=zeros(nn,1);
Res_train=zeros(nn,1);
RPS_train_bin=zeros(nn,1);
BS_train_bin=zeros(nn,1);
Rel_train_bin=zeros(nn,1);
Res_train_bin=zeros(nn,1);
mks_train=zeros(nn,1);
lss_train=zeros(nn,1);
gds_train=zeros(nn,1);
% GMSS_train=zeros(nn,1);
% MIGN_train=zeros(nn,1);
RPS_train_drought=zeros(nn,1);
BS_train_drought=zeros(nn,1);
Rel_train_drought=zeros(nn,1);
Res_train_drought=zeros(nn,1);
RPS_train_drought_bin=zeros(nn,1);
BS_train_drought_bin=zeros(nn,1);
Rel_train_drought_bin=zeros(nn,1);
Res_train_drought_bin=zeros(nn,1);
mks_train_drought=zeros(nn,1);
lss_train_drought=zeros(nn,1);
gds_train_drought=zeros(nn,1);
% GMSS_train_drought=zeros(nn,1);
% MIGN_train_drought=zeros(nn,1);
e_train_drought=eye(n_c,n_c);
e_train=eye(nm,nm);
e_train_drought_tot=zeros(n_c,n_c,nn);
e_train_tot=zeros(nm,nm,nn);
Relative_state=zeros(nm,nn);
Relative_state_drought=zeros(n_c,nn);
m_e_train_temp_t=zeros(nm,nn);
m_e_train_drought_temp_t=zeros(n_c,nn);
m_e_train_temp_t_bin=zeros(nm,nn);
m_e_train_drought_temp_t_bin=zeros(n_c,nn);
for ij=1:nn%number of location points within the basin
   ind_state=((index_cat(m:end,ij)-1)*n_ci)+index_cat_ci_bsn_vec;
   index_cat_ds_calib=index_cat(1:m,ij);
   countIndex_cat_ds_calib=histc(index_cat_ds_calib,1:n_c);
   probIndex_cat_ds_calib=countIndex_cat_ds_calib/sum(countIndex_cat_ds_calib);
   ind_state_calib=((index_cat_ds_calib-1)*n_ci)+index_cat_ci_ds_calib;
   countInd_state_calib=histc(ind_state_calib,1:nm);
   probInd_state_calib=countInd_state_calib/sum(countInd_state_calib);
%%  calculating transition probabilities (e_train) using whole training data
train_drought =  sparse(index_cat_ds_calib(1:end-1),index_cat_ds_calib(2:end),1);
train =  sparse(ind_state_calib(1:end-1),ind_state_calib(2:end),1);
train_drought=full(train_drought);
train=full(train);
sum_train2_drought=sum(train_drought,2);
sum_train2=sum(train,2);
e1_train_drought= train_drought./ repmat(sum_train2_drought,1,size(train_drought,2));
e1_train= train./ repmat(sum_train2,1,size(train,2));
[d1_etrain1_drought,d2_etrain1_drought]=size(e1_train_drought);
[d1_etrain1,d2_etrain1]=size(e1_train);
e_train_drought(1:d1_etrain1_drought,1:d2_etrain1_drought)=e1_train_drought;
e_train(1:d1_etrain1,1:d2_etrain1)=e1_train;
e_train_drought(isnan(e_train_drought))=0;
e_train(isnan(e_train))=0;
I_etrain=find(~(sum(e_train,2)));
nz_I_etrain=~(isempty(I_etrain));
if nz_I_etrain
        eye_I_etrain1=zeros(length(I_etrain),nm);
        as=[1:length(I_etrain);I_etrain'];
        as_prd=((as(2,:)-1)*length(I_etrain))+(as(1,:));
        eye_I_etrain1(as_prd)=1;
for i_nz=1:length(I_etrain)
    e_train(I_etrain(i_nz),:)=eye_I_etrain1(i_nz,:);
end
end
I_etrain_drought=find(~(sum(e_train_drought,2)));
nz_I_etrain_drought=~(isempty(I_etrain_drought));
if nz_I_etrain_drought
        eye_I_etrain1_drought=zeros(length(I_etrain_drought),n_c);
        as_drought=[1:length(I_etrain_drought);I_etrain_drought'];
        as_prd_drought=((as_drought(2,:)-1)*length(I_etrain_drought))+(as_drought(1,:));
        eye_I_etrain1_drought(as_prd_drought)=1;
for i_nz=1:length(I_etrain_drought)
    e_train_drought(I_etrain_drought(i_nz),:)=eye_I_etrain1_drought(i_nz,:);
end
end
e_train_tot(:,:,ij)=e_train;
e_train_drought_tot(:,:,ij)=e_train_drought;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%  Check transition probabilities to reach steady-state condition
diff_e=10;
diff_e_drought=10;
n_e_train=e_train*e_train;
n_e_train_drought=e_train_drought*e_train_drought;
n_e=2;n_e_drought=2;
while and(diff_e>0.5,n_e<1000)
   diff_e= norm((n_e_train*e_train),'fro');
   n_e=n_e+1;
   n_e_train=n_e_train*e_train;
end
while and(diff_e_drought>0.5,n_e_drought<1000)
   diff_e_drought= norm((n_e_train_drought*e_train_drought),'fro');
   n_e_drought=n_e_drought+1;
   n_e_train_drought=n_e_train_drought*e_train_drought;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Predicting
for q=m+1:(mm)
    initial_drought=index_cat(q-1,ij);
    initial=ind_state(q-m);
    cur_TC=index_cat_ci_bsn_vec(q-m+1);
    e_train_temp2_bin=zeros(nm,1);
    e_train_drought_temp2_bin=zeros(n_c,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %for 2D case(with SST information)
    e_train_temp1=e_train(initial,:);
    ind_nm=1:nm;
    ind_cur_TC=cur_TC:n_ci:nm;
    ind_nm(ind_cur_TC)=[];
    e_train_temp1(ind_nm)=0;
    I_ef_train=(~(sum(e_train_temp1)));
    if I_ef_train
            eye_I_ef_train=zeros(1,nm);
            eye_I_ef_train(ind_cur_TC)=(1/n_c);
            e_train_temp1=eye_I_ef_train;
    end
    e_train_temp2=e_train_temp1./sum(e_train_temp1);%Correcting not being unity of summation 
    e_train_temp(:,q-m,ij)=e_train_temp2;
    [~,I_bin]=max(e_train_temp2);
    e_train_temp2_bin(I_bin)=1;
    e_train_temp_bin(:,q-m,ij)=e_train_temp2_bin;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %for 1D case(without SST information)
    e_train_drought_temp1=e_train_drought(initial_drought,:);
%      if (~(sum(e_f_drought_temp1)))
%         e_f_drought_temp1=zeros(1,n_c);
%         e_f_drought_temp1(initial_drought)=1;
%     end
     if (~(sum(e_train_drought_temp1)))
        e_train_drought_temp1=zeros(1,n_c);
        e_train_drought_temp1(initial_drought)=1;
     end
    e_train_drought_temp2=e_train_drought_temp1./sum(e_train_drought_temp1);%Correcting not being unity of summation 
    e_train_drought_temp(:,q-m,ij)=e_train_drought_temp2;
    [~,I_drought_bin]=max(e_train_drought_temp2);
    e_train_drought_temp2_bin(I_drought_bin)=1;
    e_train_drought_temp_bin(:,q-m,ij)=e_train_drought_temp2_bin;
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  control_stats=ind_state(2:end);
  [N_state,~] = histcounts(control_stats,'BinMethod','integers','BinLimits',[1,nm]);
  Relative_state(:,ij)=N_state/(l_q-1);
  control_stats_drought=index_cat(m+1:end,ij);
  [N_state_drought,~] = histcounts(control_stats_drought,'BinMethod','integers','BinLimits',[1,n_c]);
  Relative_state_drought(:,ij)=N_state_drought/(l_q-1);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  e_train_temp_t=e_train_temp(:,:,ij);
  m_e_train_temp_t(:,ij)=mean(e_train_temp_t,2);
  e_train_drought_temp_t=e_train_drought_temp(:,:,ij);
  m_e_train_drought_temp_t(:,ij)=mean(e_train_drought_temp_t,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For 0 and 1 Forecasts
  e_train_temp_t_bin=e_train_temp_bin(:,:,ij);
  m_e_train_temp_t_bin(:,ij)=mean(e_train_temp_t_bin,2);
  e_train_drought_temp_t_bin=e_train_drought_temp_bin(:,:,ij);
  m_e_train_drought_temp_t_bin(:,ij)=mean(e_train_drought_temp_t_bin,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [~,MRPS1_train]=rps(control_stats,e_train_temp_t);
  [~,MRPS1_train_drought]=rps(control_stats_drought,e_train_drought_temp_t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [BS1_train,Rel1_train,Res1_train]=bscore(control_stats,e_train_temp_t);
  [BS1_train_drought,Rel1_train_drought,Res1_train_drought]=bscore(control_stats_drought,e_train_drought_temp_t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [mks1_train]=ksharp(e_train_temp_t);
  [mks1_train_drought]=ksharp(e_train_drought_temp_t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [Likelihood_train]=lss(control_stats,e_train_temp_t);
  [Likelihood_train_drought]=lss(control_stats_drought,e_train_drought_temp_t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [GDS1_train]=gds(control_stats,e_train_temp_t,ng=);
  [GDS1_train_drought]=gds(control_stats_drought,e_train_drought_temp_t,ng=);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   [GMSS1_train]=gmss(control_stats,e_train_temp_t);
%   [GMSS1_train_drought]=gmss(control_stats_drought,e_train_drought_temp_t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   [IGN1_train,MIGN1_train]=ign(control_stats,e_train_temp_t);
%   [IGN1_train_drought,MIGN1_train_drought]=ign(control_stats_drought,e_train_drought_temp_t);
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % For 0 and 1 Forecasts
  [~,MRPS1_train_bin]=rps(control_stats,e_train_temp_t_bin);
  [~,MRPS1_train_drought_bin]=rps(control_stats_drought,e_train_drought_temp_t_bin);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [BS1_train_bin,Rel1_train_bin,Res1_train_bin]=bscore(control_stats,e_train_temp_t_bin);
  [BS1_train_drought_bin,Rel1_train_drought_bin,Res1_train_drought_bin]=bscore(control_stats_drought,e_train_drought_temp_t_bin);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%	classic MC model for whole basin (2D)
    RPS_train_bin(ij)=MRPS1_train_bin;
    RPS_train(ij)=MRPS1_train;
    BS_train_bin(ij)=BS1_train_bin;
    BS_train(ij)=BS1_train;
    Rel_train_bin(ij)=Rel1_train_bin;
    Rel_train(ij)=Rel1_train;
    Res_train_bin(ij)=Res1_train_bin;
    Res_train(ij)=Res1_train;
%   GMSS_train(ij)=GMSS1_train;
%   MIGN_train(ij)=MIGN1_train;
    lss_train(ij)=Likelihood_train;
    gds_train(ij)=GDS1_train;
    mks_train(ij)=mks1_train;
%%	classic MC model for whole basin only dealing with drought index (1D)
    RPS_train_drought_bin(ij)=MRPS1_train_drought_bin;
    RPS_train_drought(ij)=MRPS1_train_drought;
    BS_train_drought_bin(ij)=BS1_train_drought_bin;
    BS_train_drought(ij)=BS1_train_drought;
    Rel_train_drought_bin(ij)=Rel1_train_drought_bin;
    Rel_train_drought(ij)=Rel1_train_drought;
    Res_train_drought_bin(ij)=Res1_train_drought_bin;
    Res_train_drought(ij)=Res1_train_drought;
%   GMSS_train_drought(ij)=GMSS1_train_drought;
%   MIGN_train_drought(ij)=MIGN1_train_drought;
    lss_train_drought(ij)=Likelihood_train_drought;
    gds_train_drought(ij)=GDS1_train_drought;
    mks_train_drought(ij)=mks1_train_drought;
end
ave_RPS=mean(RPS_train);
ave_BS=mean(BS_train);
ave_Rel_train= mean(Rel_train);
ave_Res_train=mean(Res_train);
ave_RPS_drought=mean(RPS_train_drought);
ave_BS_drought=mean(BS_train_drought);
ave_Rel_train_drought= mean(Rel_train_drought);
ave_Res_train_drought=mean(Res_train_drought);
ave_mks_drought=mean(mks_train_drought);
ave_mks=mean(mks_train);
ave_lss_drought=mean(lss_train_drought);
ave_lss=mean(lss_train);
ave_gds_drought=mean(gds_train_drought);
ave_gds=mean(gds_train);
ave_m_e_train_temp_t=mean(m_e_train_temp_t,2);
ave_m_e_train_drought_temp_t=mean(m_e_train_drought_temp_t,2);
ave_Relative_state=mean(Relative_state,2);
ave_Relative_state_drought=mean(Relative_state_drought,2);
% ave_GMSSS=mean(GMSS_train);
% ave_GMSS_drought=mean(GMSS_train_drought); 
ave_RPS_bin=mean(RPS_train_bin);
ave_BS_bin=mean(BS_train_bin);
ave_Rel_train_bin= mean(Rel_train_bin);
ave_Res_train_bin=mean(Res_train_bin);
ave_RPS_drought_bin=mean(RPS_train_drought_bin);
ave_BS_drought_bin=mean(BS_train_drought_bin);
ave_Rel_train_drought_bin= mean(Rel_train_drought_bin);
ave_Res_train_drought_bin=mean(Res_train_drought_bin);
ave_m_e_train_temp_t_bin=mean(m_e_train_temp_t_bin,2);
ave_m_e_train_drought_temp_t_bin=mean(m_e_train_drought_temp_t_bin,2);
m_e_train_drought_tot=mean(e_train_drought_tot,3);
m_e_train_tot=mean(e_train_tot,3);
save('F:\Program Files\MATLAB\MATLAB Production Server\R2015a\bin\work_thesis\Res_MAT\final.mat',...
    'ind_basin','ave_RPS','ave_BS','ave_RPS_drought','ave_BS_drought',...
    'ave_RPS_drought','ave_BS_drought','ave_RPS','ave_BS',...
    'n_e','n_e_drought','RPS_train',...
    'Res_train_drought','Rel_train_drought','Rel_train','Res_train',...
    'e_train_drought','e_train','ave_mks','ave_mks_drought','mks_train',...
    'mks_train_drought','RPS_train_drought','BS_train_drought','BS_train',...
    'ave_Rel_train','ave_Res_train','ave_Rel_train_drought',...
    'ave_Res_train_drought','e_train_drought_temp','e_train_temp','Relative_state',...
    'Relative_state_drought','ave_Relative_state', 'ave_Relative_state_drought', ...
    'm_e_train_temp_t','m_e_train_drought_temp_t','ave_m_e_train_temp_t',...
    'ave_m_e_train_drought_temp_t','ave_RPS_bin','ave_BS_bin',...
    'ave_Rel_train_bin','ave_Res_train_bin','ave_RPS_drought_bin',...
    'ave_BS_drought_bin','ave_Rel_train_drought_bin','ave_Res_train_drought_bin',...
    'ave_m_e_train_temp_t_bin','ave_m_e_train_drought_temp_t_bin','m_e_train_tot',...
    'm_e_train_drought_tot','ave_lss_drought','lss_train_drought','ave_lss',...
    'lss_train','ave_gds_drought','gds_train_drought','ave_gds','gds_train')

















