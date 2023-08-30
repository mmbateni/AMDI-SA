clc;
clearvars -except  best_cop_indx best_cop_values H_V h_adf h_kpss h_pp index_cat...
    point_MK_bsn point_sp_bsn point_adf_bsn point_adf_v_bsn point_kpss_bsn point_lmc_v_bsn...
    point_pp_const_v_bsn  point_adf_const_v_bsn point_pp_v_bsn point_pp_tr_v_bsn;
ind_xls=xlsread('index_basin.xlsx','ind','D2:F80');
ind_xls_t_bsn=ind_xls(:,1);
n_bs=length(ind_xls_t_bsn);
scale=1:12;
n_scale=length(scale);
overlap=n_scale-1;
bc=size(best_cop_indx,2);
bd=size(best_cop_indx,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reading XLS files for Information of Mediterranean Sea
%CI_tot has dimentions of (number of point within meditranean*(28*12))
%SST_from1_1983_to_12_2010_1to1_tot_Medi
CI_tot1=xlsread('SST_from1_1983to12_2010_1to1_tot_Medi.xlsx','sst');
CI_tot2=CI_tot1';
nn=16*43;
BB = reshape(CI_tot2,nn,[]);
CI_tot3=BB';
t_reg=(1:bd)';
t1_reg=[ones(bd,1) t_reg];
Ci_ds_bsn=[];
lag_TC_bn=[];
index_cat_ci_bsn=[];
index_cat_ci_or_bsn=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reading XLS files for Domain of Mediterranean Sea
ind_xls_md=xlsread('index_basin.xlsx','ind_Medi','D2:F130');
ind_xls_md_t=ind_xls_md(:,1)';
zr_mult_ind=zeros(size(BB,2),nn);
zr_mult_ind(:,ind_xls_md_t)=ones(size(BB,2),size(ind_xls_md_t,2));
CI_tot_zr=zr_mult_ind.*CI_tot3;
CI_tot_zr( :, all( ~any( CI_tot_zr ), 1 ) ) = [];
CI_tot=CI_tot_zr;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Defining Some Constants
%Classification Could be replaced with Below Values
%y_quantile=[0    0.0909    0.1818    0.2727    0.3636    0.4545    0.5454    0.6363 0.7272    0.8181    0.9090    1];
%[0 0.2 0.4 0.6 0.8 1];%[0 0.33 0.667 1];%[0 0.159 0.841 1];
%[0 0.023 0.067 0.159 0.841 0.933 0.977 1];
%[0    0.1429    0.2858    0.4287    0.5716    0.7145    0.8574 1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y_quantile=[0 0.333 0.667 1];
hk_limit=12;
n_sig=size(CI_tot,2);
m=ceil(bd*(2/3));
Cov_calib=zeros(n_sig,n_bs,(hk_limit-1));
norm_lags=zeros(1,hk_limit);
Ci_ds=zeros(bd,bd-m+1);
index_cat_ci_or=zeros(bd,bd-m+1);
index_cat_ci=zeros(bd,bd-m+1);
nseas=12;
dim_t_ind=(floor(325/nseas)+1);
n_tot_ind=(dim_t_ind*nseas);
t_ind2=reshape(1:n_tot_ind,nseas,dim_t_ind);
t_ind1=t_ind2;
t_ind3=t_ind2';
t_ind1((325+1):(dim_t_ind*nseas))=0;
t_ind= t_ind1';
Repros_tr=[];
Ci_tot=zeros(bd,bd-m+1);
Ci=zeros(bd,1);
H_MK_ci1 = false(1, bc);
H_SP_ci1= false(1, bc);
CI_calib=zeros(n_sig,m,(hk_limit-1));
Ci_calib_f_tot=zeros(n_tot_ind,1);
Ci_s_calib=zeros(m,1);
index_cat_ci_bsn_vec=zeros(bd-m,1);
n_singul=1;%  n_singul can be diffrent
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%  Selecting Optimum Lag For whole Basin
 best_cop_indx_calib=best_cop_indx(1:m,ind_xls_t_bsn);
     for hk=1:(hk_limit-1)
        CI_calib(:,:,hk)=CI_tot(n_scale-hk:((m+(n_scale-1)-hk)),:)';
        Cov_calib(:,:,hk)=(CI_calib(:,:,hk)*best_cop_indx_calib)./m;
        norm_lags(hk)=norm(Cov_calib(:,:,hk),1);
     end
[~,In] = max(norm_lags);
lag_TC=In;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Calculating TC series for Calibration Period
    ind_max_calib=m+(n_scale-1)-lag_TC;
    ind_min_calib=(n_scale-lag_TC);
    CI1_calib=CI_tot(ind_min_calib:ind_max_calib,:)';
    Cov1__calib=(CI1_calib*best_cop_indx_calib)./m;
    [U_calib,S_calib,V_calib]=svd(Cov1__calib);
    Ci_calib_f=CI1_calib'*U_calib(:,n_singul);
    for ko=1:12;
        ind_CI_t_calib=t_ind3(:,ko);
        ind_CI_t_calib(ind_CI_t_calib>ind_max_calib)=[];
        ind_CI_t_calib(ind_CI_t_calib<ind_min_calib)=[];
        ind_CI_temp_calib=ind_CI_t_calib-ind_min_calib+1;
        m_Ci_calib=mean(Ci_calib_f(ind_CI_temp_calib));
        std_Ci_calib=std(Ci_calib_f(ind_CI_temp_calib));
        Ci_s_calib(ind_CI_temp_calib)=((Ci_calib_f(ind_CI_temp_calib))-m_Ci_calib)./std_Ci_calib;
    end
    edges=quantile(Ci_s_calib,y_quantile);
    edges(1)=-inf;
    edges(end)=inf;
    index_cat_ci_ds_calib=discretize(Ci_s_calib,edges);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%For each time step in Forecasting Period
for hj=m:bd
    %%  Calculating Covariance for each Forecast(Updating Continiusly))
    ind_max=(hj+(n_scale-1)-lag_TC);
    ind_min=(n_scale-lag_TC);
    CI1=CI_tot(ind_min:ind_max,:)';
    if hj==m+22;
        xyz=1;
    end
    best_cop_indx_t=best_cop_indx(1:hj,ind_xls_t_bsn);
    Cov_t=(CI1*best_cop_indx_t)./hj;
%   Cov_best=Cov(:,In(hj),hj);
%   CI_best=CI_tot(In(hj):(end-(12-In(hj))),:);
%   lag_TC(hj)=12-In(hj);
    %%   Compressing Data, Carried in whole Domain of Mediterranean, into a Single Series
    [U,S,V]=svd(Cov_t);
    Ci(1:hj)=CI1'*U(:,n_singul);
    if sum(isnan(Ci))
        xyz=1;
    end
    Ci_temp=Ci;
    Ci_temp(hj+1:end)=[];
    Ci_tot(1:hj,hj-m+1)=Ci_temp;
    ini_month=mod(hj,nseas);
    for jo=1:12;
        ind_CI_t=t_ind3(:,jo);
        ind_CI_t(ind_CI_t>ind_max)=[];
        ind_CI_t(ind_CI_t<ind_min)=[];
        ind_CI_temp=ind_CI_t-ind_min+1;
        m_Ci=mean(Ci_temp(ind_CI_temp));
        std_Ci=std(Ci_temp(ind_CI_temp));
        Ci_ds(ind_CI_temp,hj-m+1)=((Ci_temp(ind_CI_temp))-m_Ci)./std_Ci;
    end
%%  Discretizing Tele Connection Information
% edges=quantile(Ci_ds(1:hj,hj-m),y_quantile);
index_cat_ci_or(1:hj,hj-m+1)=ordinal(Ci_ds(1:hj,hj-m+1),{'L1','N','H1'},[],edges);
index_cat_ci(1:hj,hj-m+1)=discretize(Ci_ds(1:hj,hj-m+1),edges,'IncludedEdge','right');
Ci_ds_bsn=[Ci_ds_bsn,Ci_ds(:,hj-m+1)];
index_cat_ci_bsn=[index_cat_ci_bsn,index_cat_ci(:,hj-m+1)];
index_cat_ci_or_bsn=[index_cat_ci_or_bsn,index_cat_ci_or(:,hj-m+1)];  
index_cat_ci_bsn_vec(hj-m+1)=index_cat_ci_bsn(hj,hj-m+1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('F:\Program Files\MATLAB\MATLAB Production Server\R2015a\bin\work_thesis\Res_MAT\climate_TC_bestLag.mat', ...
    'index_cat_ci_bsn','index_cat_ci_or_bsn','index_cat_ci','Ci_ds','Ci_ds_bsn','lag_TC',...
    'edges','Ci_s_calib','index_cat_ci_ds_calib','index_cat_ci_bsn_vec');