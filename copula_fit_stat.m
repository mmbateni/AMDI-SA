clc;
clearvars -except var_spi_p scor_spi_p var_spi_np scor_spi_np var_smi_np ...
    score_smi_np var_smi_w_pca score_smi_w_pca score_smi_wo_pca;
ind_xls=xlsread('index_basin.xlsx','ind','D2:F80');
ind_xls_t=ind_xls(:,1);
sz_bsn=size(ind_xls_t,1);
scale=1:12;
n_scale=length(scale);
overlap=n_scale-1;
mm=size(scor_spi_np,3);
nn=(size(scor_spi_np,1)*size(scor_spi_np,2))-overlap;
cdf_cop_name={'Gaussian','t','Clayton','Frank','Gumbel','Empirical'};
cdf_cop_nameLD={'Clayton LD','Frank LD','Gumbel LD'};
index_cat=[];
y_edges_p=[0	0.023	0.055	0.097	0.212	0.309	0.691	0.788	0.903	0.945	0.977	1];
y_edges=norminv(y_edges_p,0,1);%for 11 classes
% %for 11 classes it can be[-inf,-2,-1.6,-1.3,-0.8,-0.5,0.5,0.8,1.3,1.6,2,inf];
% flag_ind=0;
% flag_bcop=0;
% fs=2*(2^8+1);%choosed for number of data in this case
% alpha=0.01;% for MK test: could be 0.05
pd= makedist('Uniform','lower',0,'upper',1);
lag=0;%set the lag between soil moisture and percipitationit,can be diffrent
y_archem_K=zeros(nn,3,mm);
K_c_theo_MSE=zeros(nn,3,mm);
Rho_copula_all=zeros(nn,5);
cop_empri=zeros(nn,mm);
cop_stat2=ones(5,1,mm);
cop_stat1=zeros(5,1,mm);
cop_stat=[cop_stat1 cop_stat2];
cop_values=zeros(nn,6,mm);
w=zeros(nn,2,mm);
best_cop_values_1=zeros(nn,mm);
aic=zeros(mm,5);
best_cop_indx=zeros(nn,mm);
P_ks=zeros(1,5);P_cm=zeros(1,5);
CvMSTAT=zeros(1,5);H_cm=ones(1,5); H_ks=ones(1,5);CvMSTAT_ks=zeros(1,5);
% ind_seas=zeros(nn,mm);
% H_V=ones(1,mm);
% p_value_V=zeros(1,mm);
% H_V_sp=zeros(1,mm);
% hurst_ex=zeros(1,mm);
% p_value_V_sp=zeros(1,mm);
% h_adf=ones(1,mm);
% h_pp=ones(1,mm);
% h_adf_const=ones(1,mm);
% h_pp_const=ones(1,mm);
% h_adf_tr=ones(1,mm);
% h_pp_tr=ones(1,mm);
% h_kpss=zeros(1,mm);
% % h_kpss=zeros(1,14,mm);
% h_lmc=zeros(1,14,mm);
% lags=(-1)*ones(1,mm);
% lagspp=(-1)*ones(1,mm);
% lags_const=(-1)*ones(1,mm);
% lagspp_tr=(-1)*ones(1,mm);
% lags_tr=(-1)*ones(1,mm);
% lagspp_const=(-1)*ones(1,mm);
% point_tr_MK1=zeros(1,mm);
% point_tr_sp1=zeros(1,mm);
% point_adf1=zeros(1,mm);
% point_pp1=zeros(1,mm);
% point_adf1_const=zeros(1,mm);
% point_pp1_const=zeros(1,mm);
% point_adf1_tr=zeros(1,mm);
% point_pp1_tr=zeros(1,mm);
% % point_kpss1=zeros(14,mm);
% point_kpss1=zeros(1,mm);
% point_lmc1=zeros(14,mm);
% LOC=ones(nn,nn);
% LA=zeros(nn,nn);
% CSum=zeros(nn,nn);
Prob_CM=zeros(mm,5);
Prob_KS=zeros(mm,5);
Res_CM=zeros(mm,5);
Res_KS=zeros(mm,5);
Kendall_gauss=zeros(1,mm);Kendall_frank=zeros(1,mm);
Kendall_clayton=zeros(1,mm);Kendall_gumbel=zeros(1,mm);
Kendall_t=zeros(1,mm);
C_tau1=zeros(1,mm);
% C_tau2=zeros(1,mm);
% S=zeros(1,(nn*nn));
K_c_bestcop=zeros(nn,mm);
K_c=zeros(nn,mm);
K_c_t=zeros(nn,mm);
K_c_gauss=zeros(nn,mm);
K_c_ne=zeros(nn,2,mm);
K_c_ro=zeros(nn,mm);
K_c_theo=zeros(nn,3,mm);
K_c_theo_KS=zeros(nn,5,mm);
% point_tr_MK_bsn=zeros(1,sz_bsn);
% point_tr_sp_bsn=zeros(1,sz_bsn);
% point_lp_hrst1=zeros(1,sz_bsn);
% point_lp_hrst_bsn=zeros(1,sz_bsn);
% point_adf_bsn=ones(1,sz_bsn);
% point_pp_bsn=ones(1,sz_bsn);
% point_adf_const_bsn=ones(1,sz_bsn);
% point_pp_const_bsn=ones(1,sz_bsn);
% point_adf_tr_bsn=ones(1,sz_bsn);
% point_pp_tr_bsn=ones(1,sz_bsn);
% % point_kpss_bsn=zeros(14,sz_bsn);
% point_kpss_bsn=zeros(1,sz_bsn);
% point_lmc_bsn=zeros(14,sz_bsn);
% cop_v=zeros(nn,sz_bsn);
% h_unif_cop=zeros(sz_bsn,2);
% cop_empri_nelson=zeros(nn,mm);
% z_war=zeros(1,nn);
% t_emp_war=zeros(nn,mm);
% ind_H_KC_bs=zeros(sz_bsn,1);
% best_cop= cell(1,sz_bsn);
% best_cop_KC=cell(1,sz_bsn);
K_c_best=zeros(nn,mm);
cdf_cop_param=zeros(nn,6,mm);
h_uniform_ks=ones(1,mm);
h_uniform_chi2=ones(1,mm);
ind_H_KC=zeros(mm,1);
P_KC_ks=zeros(mm,3);
H_KC_ks=ones(mm,3);
P_KC_ks_MSE=zeros(mm,3);
H_KC_ks_MSE=ones(mm,3);
Rho_K_c_all=zeros(mm,3);
rmse_all=zeros(mm,3);
good_cop = cell(5,mm);
best_cop_all= cell(1,mm);
best_cop_all_KC= cell(1,mm);
I_rmse=zeros(mm,3);
nbins = 50; % number of bin for uniform dist.
edges = linspace(0,1,nbins+1); % edges of the bins
N=nn; % sample size
E = N/nbins*ones(nbins,1); % expected value (equal for uniform dist)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main Loop of The Program Which Is Repeated For Each Spatial Point
for ui=1:mm
    percip_ind1=scor_spi_np(:,:,ui);
    percip_ind2=percip_ind1';
    u2=percip_ind2(n_scale:end);
    u1=u2(1:(end-lag));
%%  Working with ranks instead of raw data
	[ranking_SP1] = tiedrank(u1);
    ranking_SP=ranking_SP1';
    %[ranking_SP] =ceil(tiedrank(u1));
    AS_SP=(ranking_SP)/(length(u1)+1);
    smoist_ind1=score_smi_np(:,:,ui);
    smoist_ind2=smoist_ind1';
    v2=smoist_ind2(n_scale:end);
    v1=v2(lag+1:end);
    [ranking_SM1] = tiedrank(v1);
    ranking_SM=ranking_SM1';
    %[ranking_SM] =ceil(tiedrank(v1));
    AS_SM=(ranking_SM)/(length(v1)+1);
    w(:,:,ui)=[AS_SP AS_SM];
    w_ui=w(:,:,ui);
%%  Fitting Parametric copulas
    Rho1 = copulafit('Gaussian',w_ui);
    Rho=mean(diag(fliplr(Rho1)));%Make it Scalar 
    F1 = negloglike_gs(Rho,w_ui);
    aic(ui,1)=aicbic((-1*F1),1);
    Kendall_gauss(ui)=copulastat('Gaussian',Rho);
    [rhohat1,nuhat] = copulafit('t',w_ui); %also returns an approximate 95 confidence interval, nuci, for the degrees of freedom estimated in nuhat.
    rhohat=mean(diag(fliplr(rhohat1)));%Make it Scalar 
    F2=negloglike_t(nuhat,chol(rhohat),w_ui);
    aic(ui,2)=aicbic((-1*F2),2);
    Kendall_t(ui)=copulastat('t',rhohat,nuhat);
    [alpha_clayton] = copulafit('Clayton',w_ui);%also returns an approximate 95% confidence interval, paramci, for the copula parameter estimated in paramhat.
    F3=negloglike_clayton(alpha_clayton,w_ui);
    aic(ui,3)=aicbic((-1*F3),1);
    Kendall_clayton(ui)=copulastat('Clayton',alpha_clayton);
    [alpha_frank] = copulafit('Frank',w_ui);%also returns an approximate 95% confidence interval, paramci, for the copula parameter estimated in paramhat.
    F4=negloglike_frank(alpha_frank,w_ui);
    aic(ui,4)=aicbic((-1*F4),1);
    Kendall_frank(ui)=copulastat('Frank',alpha_frank);
    [alpha_gumbel] = copulafit('Gumbel',w_ui);%also returns an approximate 95% confidence interval, paramci, for the copula parameter estimated in paramhat.
    F5=negloglike_gumbel(alpha_gumbel,w_ui);
    aic(ui,5)=aicbic((-1*F5),1);
    Kendall_gumbel(ui)=copulastat('Gumbel',alpha_gumbel);
    y_gauss = copulacdf('Gaussian',w_ui,Rho);
    y_t = copulacdf('t',w_ui,rhohat,nuhat);
    y_clayton = copulacdf('Clayton',w_ui,alpha_clayton);
    y_frank = copulacdf('Frank',w_ui,alpha_frank);
    y_gumbel = copulacdf('Gumbel',w_ui,alpha_gumbel);
    Rho_copula_all(ui,:)=[Rho rhohat alpha_clayton alpha_frank alpha_gumbel];
    %for best fitted based on Kendall fcn.
    [Rho_K_c,rmse]=copulafit2(w_ui);
    Rho_K_c_all(ui,:)=Rho_K_c;
    rmse_all(ui,:)=rmse;
    y_clayton_K = copulacdf('Clayton',w_ui,Rho_K_c(1));
    y_frank_K = copulacdf('Frank',w_ui,Rho_K_c(2));
    y_gumbel_K = copulacdf('Gumbel',w_ui,Rho_K_c(3));
    y_archem_K(:,:,ui)=[y_clayton_K y_frank_K y_gumbel_K];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Calculating Copula Values (Empirical).
    cop_empri1=empcop(w_ui);
    cop_empri(:,ui)=cop_empri1;
    cop_values(:,6,ui)=cop_empri1;
    cdf_cop_param(:,:,ui)=[y_gauss y_t y_clayton y_frank y_gumbel cop_empri1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculating Kendall's coefficient (Empirical)
    C_tau1(ui)= corr(w_ui(:,1),w_ui(:,2),'type','kendall');%By MATLAB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Using AIC and GOF Tests to select the best
    aic_ui=aic(ui,:);
    [~,I_cand_best]=min(aic_ui);
    A_aic_ui =sort(aic_ui(:));
    I_cand_2nd = A_aic_ui(2);
    I_cand_3rd = A_aic_ui(3);
    I_cand_4th = A_aic_ui(4);
    I_cand_5th = A_aic_ui(5);
    [I_best]=6;
    for hj=1:5
       y_cop_param=cdf_cop_param(:,hj,ui);
      [H_cm(hj),P_cm(hj),CvMSTAT(hj)]=cmtest2(y_cop_param,cop_empri(:,ui),0.05);
      [H_ks(hj),P_ks(hj),CvMSTAT_ks(hj)]=kstest2(y_cop_param,cop_empri(:,ui),0.05);
       if all([~H_cm(hj),~H_ks(hj),eq(hj,I_cand_best)])
        good_cop(hj,ui)={cdf_cop_name(hj)};
        cop_stat(hj,:,ui)=[P_cm(hj) P_ks(hj)];
        cdf_cop_param_temp=cdf_cop_param(:,hj,ui);
        cop_values(:,hj,ui)=cdf_cop_param_temp;
        [I_best]=I_cand_best;
        elseif  all([~H_cm(hj),~H_ks(hj),eq(hj,I_cand_2nd)])
       good_cop(hj,ui)={cdf_cop_name(hj)};
        cop_stat(hj,:,ui)=[P_cm(hj) P_ks(hj)];
        cdf_cop_param_temp=cdf_cop_param(:,hj,ui);
        cop_values(:,hj,ui)=cdf_cop_param_temp;
        [I_best]=I_cand_2nd;
       elseif  all([~H_cm(hj),~H_ks(hj),eq(hj,I_cand_3rd)])
       good_cop(hj,ui)={cdf_cop_name(hj)};
        cop_stat(hj,:,ui)=[P_cm(hj) P_ks(hj)];
        cdf_cop_param_temp=cdf_cop_param(:,hj,ui);
        cop_values(:,hj,ui)=cdf_cop_param_temp;
        [I_best]=I_cand_3rd;
       elseif  all([~H_cm(hj),~H_ks(hj),eq(hj,I_cand_4th)])
       good_cop(hj,ui)={cdf_cop_name(hj)};
        cop_stat(hj,:,ui)=[P_cm(hj) P_ks(hj)];
        cdf_cop_param_temp=cdf_cop_param(:,hj,ui);
        cop_values(:,hj,ui)=cdf_cop_param_temp;
        [I_best]=I_cand_4th;
         elseif  all([~H_cm(hj),~H_ks(hj),eq(hj,I_cand_5th)])
       good_cop(hj,ui)={cdf_cop_name(hj)};
        cop_stat(hj,:,ui)=[P_cm(hj) P_ks(hj)];
        cdf_cop_param_temp=cdf_cop_param(:,hj,ui);
        cop_values(:,hj,ui)=cdf_cop_param_temp;
        [I_best]=I_cand_5th;      
       end   
    end
    Res_CM(ui,:)=H_cm;
    Res_KS(ui,:)=H_ks;
    Prob_CM(ui,:)=P_cm;
    Prob_KS(ui,:)=P_ks;
    aic_ui=aic(ui,:);
    best_cop_all(ui)={cdf_cop_name(I_best)};
    best_cop_values_1(:,ui)=cop_values(:,I_best,ui);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculating Empirical Kendall Function Based on Nelson's Work
    K_c_bestcop(:,ui)=empkend(w_ui,best_cop_values_1(:,ui));
    K_c(:,ui)=empkend(w_ui,cop_empri1);
    K_c_t(:,ui)=empkend(w_ui,y_t);
    K_c_gauss(:,ui)=empkend(w_ui,y_gauss);
    K_c_ne(:,1:2,ui)=[K_c_gauss(:,ui) K_c_t(:,ui)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculating Theoritical Kendall Function 
    K_c_clayton=y_clayton.*((1+alpha_clayton-(y_clayton.^alpha_clayton))/alpha_clayton);
    K_c_gumbel=y_gumbel-((y_gumbel.*log(y_gumbel))/(alpha_gumbel+1));
    ex_K_frank=exp(-alpha_frank*y_frank);
    K_c_frank=y_frank+(((1-ex_K_frank)./(alpha_frank*ex_K_frank)).*(log((1-exp(-alpha_frank))./(1-ex_K_frank))));
    K_c_theo(:,:,ui)=[K_c_clayton K_c_frank K_c_gumbel];
    K_c_theo1=K_c_theo(:,:,ui);
    %for best fitted based on Cop fcn.
    K_c_clayton_KS=cop_empri1.*((1+alpha_clayton-(cop_empri1.^alpha_clayton))/alpha_clayton);
    K_c_gumbel_KS=cop_empri1-((cop_empri1.*log(cop_empri1))/(alpha_gumbel+1));
    ex_K_frank_KS=exp(-alpha_frank*cop_empri1);
    K_c_frank_KS=cop_empri1+(((1-ex_K_frank_KS)./(alpha_frank*ex_K_frank_KS)).*(log((1-exp(-alpha_frank))./(1-ex_K_frank_KS))));
    K_c_theo_KS(:,:,ui)=[K_c_gauss(:,ui) K_c_t(:,ui) K_c_clayton_KS K_c_frank_KS K_c_gumbel_KS];
    K_c_theo1_KS=K_c_theo_KS(:,:,ui);
    %for best fitted based on Kendall fcn.
    K_c_clayton_MSE=cop_empri1.*((1+Rho_K_c(1)-(cop_empri1.^Rho_K_c(1)))/Rho_K_c(1));
    ex_K_frank_MSE=exp(-Rho_K_c(2)*cop_empri1);
    K_c_frank_MSE=cop_empri1+(((1-ex_K_frank_MSE)./(Rho_K_c(2)*ex_K_frank_MSE)).*(log((1-exp(-Rho_K_c(2)))./(1-ex_K_frank_MSE))));
    K_c_gumbel_MSE=cop_empri1-((cop_empri1.*log(cop_empri1))/(Rho_K_c(3)+1));
    K_c_theo_MSE(:,:,ui)=[K_c_clayton_MSE K_c_frank_MSE K_c_gumbel_MSE];
    K_c_theo1_MSE=K_c_theo_MSE(:,:,ui);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Comparing Kendall Function with Empirical ones
    for ks=1:5
   [H_KC_ks(ui,ks),P_KC_ks(ui,ks),~]=kstest2(K_c_theo1_KS(:,ks),K_c(:,ui),0.05);%1 means not similar distributions
    end
    for kse=1:3
    [H_KC_ks_MSE(ui,kse),P_KC_ks_MSE(ui,kse),~]=kstest2(K_c_theo1_MSE(:,kse),K_c(:,ui),0.05);%1 means not similar distributions
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Selecting the Best based on Closness of Kendall Function to Empirical Kc
    H_KC=~(H_KC_ks(ui,:));
    [~,I_min]=min(rmse);
    [Y,I_rmse(ui,:)]=sort(rmse);
    I_min_2=find(I_rmse(ui,:)==2);
    I_min_3=find(I_rmse(ui,:)==3);
    H_KC_MSE=~(H_KC_ks_MSE(ui,:));
    K_c_best(:,ui)=K_c(:,ui);
    best_cop_all_KC(ui)={cdf_cop_name(6)};
    ind_H_KC(ui)=H_KC(I_best);
    if (I_best==6)
       K_c_best(:,ui)=K_c(:,ui);
       best_cop_all_KC(ui)={cdf_cop_name(6)};
    elseif  H_KC(I_best)
        if I_best>2
        K_c_best(:,ui)=K_c_theo1(:,I_best-2);
        best_cop_all_KC(ui)={cdf_cop_name(I_best)};
        else
        K_c_best(:,ui)=K_c_ne(:,I_best);
        best_cop_all_KC(ui)={cdf_cop_name(I_best)};
        end
    end
    ind_K_c_best=~K_c_best(:,ui);
    if all(ind_K_c_best)
    if H_KC_MSE(I_min_3)
        K_c_best(:,ui)=K_c_theo1_MSE(:,I_min_3);
        best_cop_all_KC(ui)={cdf_cop_nameLD(I_min_3)};
    elseif H_KC_MSE(I_min_2)
        K_c_best(:,ui)=K_c_theo1_MSE(:,I_min_2);
        best_cop_all_KC(ui)={cdf_cop_nameLD(I_min_2)};
    elseif H_KC_MSE(I_min)
        K_c_best(:,ui)=K_c_theo1_MSE(:,I_min);
        best_cop_all_KC(ui)={cdf_cop_nameLD(I_min)};
    end
    end
    ind_K_c_best=~K_c_best(:,ui);
    if all(ind_K_c_best)
    K_c_best(:,ui)=K_c(:,ui);
    best_cop_all_KC(ui)={cdf_cop_name(6)};   
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Transformation of data to Normal distribution
% Data rearrangement and testing uniformity
    input_Cop=K_c_best(:,ui);
    h_uniform_ks(ui) = kstest(input_Cop,'CDF',pd,'Alpha',0.05);%1 means NOT uniform
    [h_uniform_chi2(ui),~,~] = chi2gof(input_Cop,'Expected',E,'Edges',edges);%1 means NOT uniform
    V=norminv(input_Cop,0,1);
    V(V==-inf)=-4;
    V(V==inf)=4;
    best_cop_indx(:,ui)=V;
    index_cat_1=discretize(best_cop_indx(:,ui),y_edges,'IncludedEdge','right');
    index_cat=[index_cat,index_cat_1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Testing for Stationarity (Within the Basin)
%  Mann-Kendal Test
%      if sum(ind_xls_t==ui)>0 
%             flag_ind=flag_ind+1;
%             best_cop_KC(flag_ind)=best_cop_all_KC(ui);
%             best_cop(flag_ind)=best_cop_all(ui);
%             cop_v(:,flag_ind)=best_cop_indx(:,ui);
%             ind_H_KC_bs(flag_ind)=ind_H_KC(ui);
%             h_unif_cop(flag_ind,:)=[h_uniform_chi2(ui) h_uniform_ks(ui)];
%             [H_V(ui),p_value_V(ui)]=Mann_Kendall(V,alpha,1);
%             point_tr_MK1(ui)=abs(H_V(ui)-1)*ui;
%             point_tr_MK_bsn(flag_ind)=abs(H_V(ui)-1)*ui;
%             [H_V_sp(ui),p_value_V_sp(ui)]=SpearmanRho(V,alpha,1);
%             point_tr_sp1(ui)=(abs(H_V_sp(ui)))*ui;
%             point_tr_sp_bsn(flag_ind)=(abs(H_V_sp(ui)))*ui;
%             [hurst_ex(ui)]=hurst(V);
%             point_lp_hrst1(ui)=round(hurst_ex(ui)-0.0000001)*ui;
%             point_lp_hrst_bsn(flag_ind)=(round(hurst_ex(ui)-0.0000001))*ui;
            %   Investigating the Spectrum
%             [pxxv,vw,pxxcv] = periodogram(V,[],[],fs,'ConfidenceLevel', 0.99);
%             ind_opv=(or (pxxv>pxxcv(:,2),pxxv<pxxcv(:,1) ));
%             if (sum(ind_opv)>0)
%                ind_seas(:,ui)=find(ind_opv);
%             end
        %   Unit Root(Stationary) Tests
%             [h_adf(ui),~,~,~,lags(ui),~]=augdfautolag(V,0);
%             [h_pp(ui),~,~,~,lagspp(ui),~]=ppautolag(V,0);
%             [h_adf_const(ui),~,~,~,lags_const(ui),~]=augdfautolag(V,1);
%             [h_pp_const(ui),~,~,~,lagspp_const(ui),~]=ppautolag(V,1);
%             [h_adf_tr(ui),~,~,~,lags_tr(ui),~]=augdfautolag(V,2);
%             [h_pp_tr(ui),~,~,~,lagspp_tr(ui),~]=ppautolag(V,2);
%             point_pp1_tr(ui)=abs(h_pp_tr(ui)-1)*ui;
%             point_pp_tr_bsn(flag_ind)=abs(h_pp_tr(ui)-1)*ui;
%             point_adf1_tr(ui)=abs(h_adf_tr(ui)-1)*ui;
%             point_adf_tr_bsn(flag_ind)=abs(h_adf_tr(ui)-1)*ui;
%             point_pp1_const(ui)=abs(h_pp_const(ui)-1)*ui;
%             point_pp_const_bsn(flag_ind)=abs(h_pp_const(ui)-1)*ui;
%             point_adf1_const(ui)=abs(h_adf_const(ui)-1)*ui;
%             point_adf_const_bsn(flag_ind)=abs(h_adf_const(ui)-1)*ui;
%             point_pp1(ui)=abs(h_pp(ui)-1)*ui;
%             point_pp_bsn(flag_ind)=abs(h_pp(ui)-1)*ui;
%             point_adf1(ui)=abs(h_adf(ui)-1)*ui;
%             point_adf_bsn(flag_ind)=abs(h_adf(ui)-1)*ui;%non stationary points
%             [~,h_kpss(:,ui)]=wtest(V,1,0,[1;floor((min(nn/3,12)*((nn/100)^0.25)))],0.99);
% %%%             h_kpss(:,:,ui)=kpsstest(V,'trend',false,'lags',(floor(4*((nn/100)^(0.25))):floor(sqrt(nn))));
%             h_lmc(:,:,ui)=lmctest(V,'trend',false,'lags',(floor(4*((nn/100)^(0.25))):floor(sqrt(nn))));
%             point_kpss1(ui)=(h_kpss(:,ui)).*ui;
%             point_kpss_bsn(flag_ind)=(h_kpss(:,ui)).*ui;%non stationary points
% %%%             point_kpss1(:,ui)=(h_kpss(:,:,ui)).*ui;
% %%%             point_kpss_bsn(:,flag_ind)=(h_kpss(:,:,ui)).*ui;
%             point_lmc1(:,ui)=(h_lmc(:,:,ui)).*ui;
%             point_lmc_bsn(:,flag_ind)=(h_lmc(:,:,ui)).*ui;
%     end
end
mu=length(y_edges)/2;
index_cat=((-1)*(index_cat-mu))+mu;%Now, 1 means wettest and 11(end) means driest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [~,~,point_MK_v] = find(point_tr_MK1);
% [~,~,point_sp_v] = find(point_tr_sp1);
% [~,~,point_hr_v] = find(point_lp_hrst1);
% [~,~,point_adf_v] = find(point_adf1);
% [~,~,point_pp_v] = find(point_pp1);
% [~,~,point_adf_tr_v] = find(point_adf1_tr);
% [~,~,point_pp_tr_v] = find(point_pp1_tr);
% [~,~,point_adf_const_v] = find(point_adf1_const);
% [~,~,point_pp_const_v] = find(point_pp1_const);
% [~,~,point_kpss_v] = find(point_kpss1);
% [~,~,point_lmc_v] = find(point_lmc1);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [~,~,point_MK_v_bsn] = find(point_tr_MK_bsn);
% [~,~,point_sp_v_bsn] = find(point_tr_sp_bsn);
% [~,~,point_hrp_v_bsn] = find(point_lp_hrst_bsn);
% [~,~,point_adf_v_bsn] = find(point_adf_bsn);
% [~,~,point_pp_v_bsn] = find(point_pp_bsn);
% [~,~,point_adf_tr_v_bsn] = find(point_adf_tr_bsn);
% [~,~,point_pp_tr_v_bsn] = find(point_pp_tr_bsn);
% [~,~,point_adf_const_v_bsn] = find(point_adf_const_bsn);
% [~,~,point_pp_const_v_bsn] = find(point_pp_const_bsn);
% [~,~,point_kpss_v_bsn] = find(point_kpss_bsn);
% [~,~,point_lmc_v_bsn] = find(point_lmc_bsn);
  save('F:\Program Files\MATLAB\MATLAB Production Server\R2015a\bin\work_thesis\Res_MAT\coupla_derived.mat','cdf_cop_param','K_c_theo',...
      'K_c_theo_KS','K_c_theo_MSE','best_cop_indx','index_cat','Rho_K_c_all',...
      'cop_empri','H_KC_ks','P_KC_ks','best_cop_all_KC','rmse_all',...
      'H_KC_ks_MSE','P_KC_ks_MSE','y_archem_K','Rho_copula_all','best_cop_all',...
      'ind_H_KC');