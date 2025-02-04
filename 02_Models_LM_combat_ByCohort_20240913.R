# 
# load("~/Desktop/Lund/20_TauPET_plasma/Data_BF2/BF2_PlasmaPET_20230907_CSF.RData")
# dataBF1=read.csv("~/Desktop/Lund/20_TauPET_plasma/Data_BF1/DataLong_BF1_plasmaTauPET_20230907.csv")
# 
# dataBF2=dataFinal
# 
# load("~/Desktop/Lund/20_TauPET_plasma/Data_AIBL/AIBL_PlasmaPET_20230907.RData")
# dataAIBL=dataFinal
# dataAIBL$sex=as.factor(ifelse(dataAIBL$sex=="Female",1,0))
# 
# dataTRIAD=read.csv("~/Desktop/Lund/20_TauPET_plasma/Data_TRIAD/TRIAD_PlasmaPET_20230907.csv")
# dataTRIAD$sex=as.factor(ifelse(dataTRIAD$sex=="F",1,0))
# 
# dataPREVENT=read.csv("~/Desktop/Lund/20_TauPET_plasma/Data_PREVENTAD/PREVENTAD_PlasmaPET_20230907.csv")
# dataPREVENT$sex=as.factor(ifelse(dataPREVENT$Sex=="Female",1,0))
# 
# dataBF2=dataBF2[,c("sid","age","sex","education_level_years_baseline_variable","timeDiff",
#                    "mmse_score_FU","mpacc","z_MTL","z_NeoT","z_ptau217","Ab_pos")]
# dataBF1=dataBF1[,c("sid","Age","sex","Education","timeDiff","mmse_score","mpacc","z_MTL","z_NeoT","z_ptau217","ab_pos")]
# dataAIBL=dataAIBL[,c("AIBL ID","Age","sex","edu","timeDiff","Neuropsych.MMSE_FU","mpacc","z_MTL","z_NeoT","z_ptau217","Ab_pos")]
# dataTRIAD=dataTRIAD[,c("ID","age","sex","edu","timeDiff","MMSE","mpacc","z_MTL","z_NeoT","z_ptau217","Ab_pos")]
# dataPREVENT=dataPREVENT[,c("PSCID.x","age_PET","sex","Education_years","time_cogn_PET","MMSE","mpacc","z_MTL","z_NeoT","z_ptau217","Ab_pos")]
# 
# colnames(dataBF2)=c("sid","age","sex","edu","timeDiff","mmse_score","mpacc","z_MTL","z_NeoT","z_ptau217","ab_pos")
# colnames(dataBF1)=c("sid","age","sex","edu","timeDiff","mmse_score","mpacc","z_MTL","z_NeoT","z_ptau217","ab_pos")
# colnames(dataAIBL)=c("sid","age","sex","edu","timeDiff","mmse_score","mpacc","z_MTL","z_NeoT","z_ptau217","ab_pos")
# colnames(dataTRIAD)=c("sid","age","sex","edu","timeDiff","mmse_score","mpacc","z_MTL","z_NeoT","z_ptau217","ab_pos")
# colnames(dataPREVENT)=c("sid","age","sex","edu","timeDiff","mmse_score","mpacc","z_MTL","z_NeoT","z_ptau217","ab_pos")
# 
# dataBF2$cohort="BF2"
# dataBF1$cohort="BF1"
# dataAIBL$cohort="AIBL"
# dataTRIAD$cohort="TRIAD"
# dataPREVENT$cohort="PREVENT"
# 
# dataBF1$tracer="FTP"
# dataBF2$tracer="RO"
# dataAIBL$tracer="MK"
# dataTRIAD$tracer="MK"
# dataPREVENT$tracer="FTP"
# 
# dataAll=rbind(dataBF2,dataBF1,dataAIBL,dataTRIAD,dataPREVENT)
# dataAll$tracer=as.factor(dataAll$tracer)
# 
# dataAll=dataAll[with(dataAll, order(sid, -timeDiff)),]
# dataAll=dataAll %>% filter(age>50)
# dataBl=dataAll %>% filter(!duplicated(sid))
# 
# 
# # c_var=c("z_MTL","z_NeoT","z_ptau217")
# # 
# # for(i_ROIs in 1:3){#:length(c_var)
# #   Q <- quantile(as.numeric(unlist(dataBl[,c_var[i_ROIs]])), probs=c(.25, .75), na.rm = TRUE)
# #   iqr <- IQR(as.numeric(unlist(dataBl[,c_var[i_ROIs]])),na.rm = TRUE)
# #   a=subset(dataBl, dataBl[,c_var[i_ROIs]] < (Q[1] - 5*iqr))
# #   b=subset(dataBl,dataBl[,c_var[i_ROIs]] > (Q[2]+5*iqr))
# #   c=c(as.numeric(unlist(a[,c_var[i_ROIs]])),as.numeric(unlist(b[,c_var[i_ROIs]])))
# #   
# #   d=which(as.numeric(unlist(dataBl[,c_var[i_ROIs]])) %in% c)
# #   
# #   list_ID=dataBl$sid[d]
# #   dataAll[which(dataAll$sid %in% list_ID),c_var[i_ROIs]]=NA
# # }
# 
# 
# dataAll=dataAll %>% filter(!is.na(z_ptau217) & !is.na(z_MTL) & !is.na(z_NeoT))
# dataBl=dataAll %>% filter(!duplicated(sid))
# 
# table1::table1(~age+sex+edu+as.factor(ab_pos)+z_ptau217+z_MTL+z_NeoT+mmse_score+mpacc+timeDiff|cohort, dataBl)
# 
# dataBl %>%
#   ggplot( aes(x=z_NeoT, fill=as.factor(cohort))) +
#   geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
#   # scale_fill_manual(values=c("#69b3a2", "#404080")) +
#   theme_classic() +
#   labs(fill="")
# 
# dataBl %>%
#   ggplot( aes(x=z_ptau217, fill=as.factor(cohort))) +
#   geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
#   # scale_fill_manual(values=c("#69b3a2", "#404080")) +
#   theme_classic() +
#   labs(fill="")
# 
# dataBl %>%
#   ggplot( aes(x=z_MTL, fill=as.factor(cohort))) +
#   geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
#   # scale_fill_manual(values=c("#69b3a2", "#404080")) +
#   theme_classic() +
#   labs(fill="")
# 
# ####
dir_res="~/Desktop/Lund/20_TauPET_plasma/Results_Merged/20240917_Rev/Cogn_decline_LM"
dir_fig="~/Desktop/Lund/20_TauPET_plasma/Figures_Merged/20240917_Rev/Cogn_decline_LM"
load("~/Desktop/Lund/20_TauPET_plasma/Data_all/Data_all_cogn_comb_20240917.RData")
############## Analyses
# dataAll=dataAll %>% filter(cohort!="Gen")
# save(dataAll,file = "~/Desktop/Lund/20_TauPET_plasma/Data_all/Data_all_cogn_comb_20240206_2.RData")


dataBl=dataAll %>% filter(!duplicated(sid))
# median(dataBl$z_ptau217_comb)

# dataAll$plasma_pos=as.factor(ifelse(dataAll$z_ptau217_comb>median(dataBl$z_ptau217_comb),1,0))#1.96,1,0))
# dataAll$MTL_pos=as.factor(ifelse(dataAll$z_MTL_comb>median(dataBl$z_MTL_comb),1,0))#1.96,1,0))
# dataAll$NeoT_pos=as.factor(ifelse(dataAll$z_NeoT_comb>median(dataBl$z_NeoT_comb),1,0))#1.96,1,0))
# dataAll$cohort=as.factor(dataAll$cohort)

f_R2=function(data,indices,l_cohortR){

  d=data[indices,]
  
  if (l_cohortR==1) {
    fml_0=as.formula(paste("slope_cog~scale(age)+sex+scale(edu)+cohort+APOE_e4",sep = ""))
    fml_00=as.formula(paste("slope_cog~scale(age)+sex+scale(edu)+cohort",sep = ""))
    fml_plasma=as.formula(paste("slope_cog~scale(age)+sex+scale(edu)+APOE_e4+cohort+z_ptau217",sep = ""))
    fml_MTL=as.formula(paste("slope_cog~scale(age)+sex+scale(edu)+APOE_e4+cohort+z_MTL",sep = ""))
    fml_NeoT=as.formula(paste("slope_cog~scale(age)+sex+scale(edu)+APOE_e4+cohort+z_NeoT",sep = ""))
    fml_mix_MTL=as.formula(paste("slope_cog~scale(age)+sex+scale(edu)+APOE_e4+cohort+z_ptau217+z_MTL",sep = ""))
    fml_mix_NeoT=as.formula(paste("slope_cog~scale(age)+sex+scale(edu)+APOE_e4+cohort+z_ptau217+z_NeoT",sep = ""))
  }else{
    fml_0=as.formula(paste("slope_cog~scale(age)+sex+scale(edu)+APOE_e4",sep = ""))
    fml_00=as.formula(paste("slope_cog~scale(age)+sex+scale(edu)",sep = ""))
    fml_plasma=as.formula(paste("slope_cog~scale(age)+sex+scale(edu)+APOE_e4+z_ptau217",sep = ""))
    fml_MTL=as.formula(paste("slope_cog~scale(age)+sex+scale(edu)+APOE_e4+z_MTL",sep = ""))
    fml_NeoT=as.formula(paste("slope_cog~scale(age)+sex+scale(edu)+APOE_e4+z_NeoT",sep = ""))
    fml_mix_MTL=as.formula(paste("slope_cog~scale(age)+sex+scale(edu)+APOE_e4+z_ptau217+z_MTL",sep = ""))
    fml_mix_NeoT=as.formula(paste("slope_cog~scale(age)+sex+scale(edu)+APOE_e4+z_ptau217+z_NeoT",sep = ""))
  }
    
  mdl_0=lm(fml_0,d)
  mdl_00=lm(fml_00,d)
  mdl_plasma=lm(fml_plasma,d)
  mdl_MTL=lm(fml_MTL,d)
  mdl_NeoT=lm(fml_NeoT,d)
  mdl_mix_MTL=lm(fml_mix_MTL,d)
  mdl_mix_NeoT=lm(fml_mix_NeoT,d)
  
  mdl_sum_0=summary(mdl_0)
  mdl_sum_00=summary(mdl_00)
  mdl_sum_plasma=summary(mdl_plasma)
  mdl_sum_MTL=summary(mdl_MTL)
  mdl_sum_NeoT=summary(mdl_NeoT)
  mdl_sum_mix_MTL=summary(mdl_mix_MTL)
  mdl_sum_mix_NeoT=summary(mdl_mix_NeoT)
  
  l_mdl=list(mdl_sum_0,mdl_sum_00,mdl_sum_plasma,mdl_sum_NeoT,mdl_sum_MTL,mdl_sum_mix_NeoT,mdl_sum_mix_MTL)
  v_diff_R2=NA
  v_diff_B=NA
  i_count=1
  # i_count2=1
  for (i in 1:7) {
    
    for (j in (i):7) {
      
      v_diff_R2[i_count]=l_mdl[[i]]$adj.r.squared-l_mdl[[j]]$adj.r.squared
      if (i>2 & j>2) {
        v_diff_B[i_count]=l_mdl[[i]]$coefficients[nrow(l_mdl[[i]]$coefficients),1]-l_mdl[[j]]$coefficients[nrow(l_mdl[[j]]$coefficients),1]
        # i_count2=i_count2+1
      }else{v_diff_B[i_count]=NA}
      i_count=i_count+1
    }
    
  }
    
  return(c(mdl_sum_0$adj.r.squared,mdl_sum_00$adj.r.squared,mdl_sum_plasma$adj.r.squared,mdl_sum_NeoT$adj.r.squared,mdl_sum_MTL$adj.r.squared,
           mdl_sum_mix_NeoT$adj.r.squared,mdl_sum_mix_MTL$adj.r.squared,v_diff_R2,v_diff_B))
}

c_colors_PET=c("#3c5488ff","#00a087ff")#c("#ffa600","#bc5090")#c("#31a354","#a1d99b")
c_colors_mix=c("#8491b4","#91d1c2ff")#c("#ff6361","#58508d")

c_age="age"
c_sex="sex"
c_edu="edu"
c_apoe="APOE_e4"
c_ID="sid"
c_cohort="cohort"
c_tracer="tracer"
n_boot=1000
# v_out=c(19,6)#c(238,235)#251:256,
# c_out=names(dataAll)[v_out]
c_out_name=c("mPACC")#,"MMSE")
n_out=length(c_out_name)

# v_PET=c(18,17)#c(244:247)
# c_PET=names(dataAll)[v_PET]
# n_PET=length(v_PET)
c_name=c("NeoT","ERC/Amygd")
# v_PET_pos=c(22,21)
# c_PET_pos=names(dataAll)[v_PET_pos]

l_ab="Abpos"#"all"#
l_DX="allEtiol"
l_out="_cov_Int_time_APOE_cohort_Q4"#"noOutlier"
l_adj="raw"#"combat"#

l_cohort=c(1,1,1,1,1,1,1,1,1)#,1)#c(0,0,1,0,0)# 1: BF2, 2: BF1, 3: AIBL, 4: TRIAD, 5: PREVENT
name_cohort=c("BF2","BF1","AIBL","TRIAD","PREVENT","WRAP","Ams","WU","MCSA")#,"Gen")

############

for (i_out in 1:n_out) {#1:n_out
  
    
  for (i_cohort in (length(l_cohort)+1)) {#1:1:(length(l_cohort))
    
    dataLong=dataAll[complete.cases(dataAll[,c("APOE_e4","sex","age","edu","z_MTL","z_NeoT","z_ptau217")]),]
    if (l_adj=="combat") {
      v_out=c(19,6)#c(238,235)#251:256,
      v_PET=c(18,17)#c(244:247)
      v_pl=16
      dataLong$plasma_pos=as.factor(ifelse(dataLong$z_ptau217_comb>quantile(dataBl$z_ptau217_comb,3/4),1,0))#1.96,1,0))
      dataLong$MTL_pos=as.factor(ifelse(dataLong$z_MTL_comb>quantile(dataBl$z_MTL_comb,3/4),1,0))#1.96,1,0))
      dataLong$NeoT_pos=as.factor(ifelse(dataLong$z_NeoT_comb>quantile(dataBl$z_NeoT_comb,3/4),1,0))#1.96,1,0))
      
    }else{
      v_out=c(7,6)#c(238,235)#251:256,
      v_PET=c(9,8)#c(244:247)
      v_pl=10
      dataLong$plasma_pos=as.factor(ifelse(dataLong$z_ptau217>quantile(dataBl$z_ptau217,3/4),1,0))#1.96,1,0))
      dataLong$MTL_pos=as.factor(ifelse(dataLong$z_MTL>quantile(dataBl$z_MTL,3/4),1,0))#1.96,1,0))
      dataLong$NeoT_pos=as.factor(ifelse(dataLong$z_NeoT>quantile(dataBl$z_NeoT,3/4),1,0))#1.96,1,0))
    }
    c_pl=names(dataLong)[v_pl]
    c_out=names(dataLong)[v_out]
    n_out=length(v_out)
    c_PET=names(dataLong)[v_PET]
    n_PET=length(v_PET)
    v_PET_pos=c(21,20)
    c_PET_pos=names(dataLong)[v_PET_pos]
    
    dataLong=dataLong[complete.cases(dataLong[,c_out[i_out]]),]
    
    s_cohort=""
    
    l_cohort=rep(0,length(name_cohort))
    
    if (i_cohort<(length(name_cohort)+1)) {
      s_cohort=name_cohort[i_cohort]
      dataLong=dataLong %>% filter(cohort==name_cohort[i_cohort])
      l_cohort[i_cohort]=1
      l_cohortRef=0
    }else{
      s_cohort="allCohorts"
      l_cohort=rep(1,length(name_cohort))
      l_cohortRef=1
    }
    
  
  
  if (l_ab=="Abpos") {
    dataLong=dataLong %>% filter(ab_pos==1)
  }
  
  # Select repetition data
  n_occur <- data.frame(table(dataLong$sid))
  list_bl=n_occur[n_occur$Freq >1,1]
  
  dataLong=filter(dataLong, sid %in% list_bl)
  
  
  fml_cog=as.formula(paste(c_out[i_out],"~timeDiff+(1+timeDiff|sid)",sep=""))
  mdl_cog=lmer(fml_cog,dataLong)
  v_cog=coef(mdl_cog)$sid
  v_cog$sid=rownames(v_cog)
  colnames(v_cog)[2]="slope_cog"
  
  df=dataLong %>% filter(!duplicated(sid))
  df=merge(df,v_cog[,c("sid","slope_cog")])
  
  results_R2=data.frame(matrix(nrow = (2*n_PET)+3,ncol = 14))
  colnames(results_R2)=c("Biomarker","name","B_plasma","p_plasma","B_PET","p_PET","R2","AICc","Ftest_plasma","pComp_plasma","Ftest_PET","pComp_PET","Ftest_0","pComp_0")
  
  results_plot=data.frame(matrix(nrow = (2*n_PET)+3,ncol = 7))
  colnames(results_plot)=c("Biomarker","name","R2","R2_l","R2_h","AICc","N")

  results_partialR2=data.frame(matrix(nrow = (2*n_PET)+3,ncol = 6))
  colnames(results_partialR2)=c("Model","R2","R2_covs","R2_plasma","R2_PET","R2_shared")
  
  results_B=data.frame(matrix(nrow = (2*n_PET)+3,ncol = 6))
  colnames(results_B)=c("Biomarker","name","B","B_h","B_l","N")
  
  results_p=data.frame(matrix(nrow = (2*n_PET)+2,ncol = (2*n_PET)+2))
  colnames(results_p)=c("Basic_A","plasma p-tau217",c_name[1],c_name[2],paste("Plasma+",c_name[1],sep = ""),paste("Plasma+",c_name[2],sep = ""))
  rownames(results_p)=c("Basic_NA","Basic_A","plasma p-tau217",c_name[1],c_name[2],paste("Plasma+",c_name[1],sep = ""))
  
  set.seed(1234)
  results_boot=boot(df,f_R2,n_boot,l_cohortR=l_cohortRef)
  
  # l_mdl=list(mdl_sum_0,mdl_sum_00,mdl_sum_plasma,mdl_sum_MTL,mdl_sum_NeoT,mdl_sum_mix_MTL,mdl_sum_mix_NeoT)
  v_diff_R2=data.frame(matrix(nrow = 7,ncol = 7))
  colnames(v_diff_R2)=c("A","NA","Plasma","NeoT","MTL","P-NeoT","P-MTL")
  rownames(v_diff_R2)=c("A","NA","Plasma","NeoT","MTL","P-NeoT","P-MTL")
  v_diff_B=data.frame(matrix(nrow = 7,ncol = 7))
  colnames(v_diff_B)=c("A","NA","Plasma","NeoT","MTL","P-NeoT","P-MTL")
  rownames(v_diff_B)=c("A","NA","Plasma","NeoT","MTL","P-NeoT","P-MTL")
  i_count=8
  # i_count2=8
  for (i1 in 1:7) {
    
    for (j1 in (i1):7) {
      
      # boot_ci_norm <- boot.ci(results_boot, type = "norm",seed=12345,index = i_count)
      results_under_H0 <- results_boot$t[,i_count] - mean(results_boot$t[,i_count])
      boot_pvalue <- mean(abs(results_under_H0) >= abs(results_boot$t0[i_count]))
      v_diff_R2[i1,j1]=boot_pvalue
      
      if (!is.na(results_boot$t0[i_count+28])) {
        results_under_H0 <- results_boot$t[,i_count+28] - mean(results_boot$t[,i_count+28])
        boot_pvalue <- mean(abs(results_under_H0) >= abs(results_boot$t0[i_count+28]))
        v_diff_B[i1,j1]=boot_pvalue
      }
      i_count=i_count+1
    }
  }
  
  # fml_0=as.formula(paste("slope_cog~scale(",c_age,")+",c_sex,"+scale(",c_edu,")+APOE_e4+",c_cohort,sep = ""))
  if (sum(l_cohort)>1) {
    fml_0=as.formula(paste("slope_cog~scale(",c_age,")+",c_sex,"+scale(",c_edu,")+APOE_e4+",c_cohort,sep = ""))
    fml_00=as.formula(paste("slope_cog~scale(",c_age,")+",c_sex,"+scale(",c_edu,")+",c_cohort,sep = ""))
    
    }else{
      fml_0=as.formula(paste("slope_cog~scale(",c_age,")+",c_sex,"+scale(",c_edu,")+APOE_e4",sep = ""))
      fml_00=as.formula(paste("slope_cog~scale(",c_age,")+",c_sex,"+scale(",c_edu,")",sep = ""))
  }
  mdl_0=lm(fml_0,df)
  mdl_sum_0=summary(mdl_0)
  ci_R2=boot.ci(results_boot,type="norm",index=1)
  c=rsq::rsq.partial(mdl_0,adj = T)    
  
  results_R2$Biomarker[(2*n_PET)+2]="Basic_APOE"
  results_R2$name[(2*n_PET)+2]="Basic_A"
  results_R2$R2[(2*n_PET)+2]=paste(round(mdl_sum_0$adj.r.squared,2),"[",round(ci_R2$normal[2],2),", ",round(ci_R2$normal[3],2),"]",sep = "")
  results_R2$AICc[(2*n_PET)+2]=paste(round(MuMIn::AICc(mdl_0),1))
  
  results_plot$Biomarker[(2*n_PET)+2]="Basic_APOE"
  results_plot$name[(2*n_PET)+2]="Basic_A"
  results_plot$R2[(2*n_PET)+2]=mdl_sum_0$adj.r.squared
  results_plot$R2_l[(2*n_PET)+2]=ci_R2$normal[2]
  results_plot$R2_h[(2*n_PET)+2]=ci_R2$normal[3]
  results_plot$AICc[(2*n_PET)+2]=MuMIn::AICc(mdl_0)
  
  results_partialR2$Model[(2*n_PET)+2]="Basic_APOE"
  results_partialR2$R2_covs[(2*n_PET)+2]=sum(c$partial.rsq[1:(length(c$partial.rsq))])
  results_partialR2$R2[(2*n_PET)+2]=mdl_sum_0$adj.r.squared
  results_plot$N=nrow(df)
  
  
  mdl_00=lm(fml_00,df)
  mdl_sum_00=summary(mdl_00)
  ci_R2=boot.ci(results_boot,type="norm",index=2)
  c=rsq::rsq.partial(mdl_00,adj = T)    
  
  results_R2$Biomarker[(2*n_PET)+3]="Basic_noAPOE"
  results_R2$name[(2*n_PET)+3]="Basic_NA"
  results_R2$R2[(2*n_PET)+3]=paste(round(mdl_sum_00$adj.r.squared,2),"[",round(ci_R2$normal[2],2),", ",round(ci_R2$normal[3],2),"]",sep = "")
  results_R2$AICc[(2*n_PET)+3]=paste(round(MuMIn::AICc(mdl_00),1))

  results_plot$Biomarker[(2*n_PET)+3]="Basic"
  results_plot$name[(2*n_PET)+3]="Basic_NA"
  results_plot$R2[(2*n_PET)+3]=mdl_sum_00$adj.r.squared
  results_plot$R2_l[(2*n_PET)+3]=ci_R2$normal[2]
  results_plot$R2_h[(2*n_PET)+3]=ci_R2$normal[3]
  results_plot$AICc[(2*n_PET)+3]=MuMIn::AICc(mdl_00)
  
  results_partialR2$Model[(2*n_PET)+3]="Basic_noAPOE"
  results_partialR2$R2_covs[(2*n_PET)+3]=sum(c$partial.rsq[1:(length(c$partial.rsq))])
  results_partialR2$R2[(2*n_PET)+3]=mdl_sum_00$adj.r.squared
  
  
  if (sum(l_cohort)>1) {
    fml_plasma=as.formula(paste("slope_cog~scale(",c_age,")+",c_sex,"+scale(",c_edu,")+APOE_e4+",c_cohort,"+",c_pl,sep = ""))
  }else{
    fml_plasma=as.formula(paste("slope_cog~scale(",c_age,")+",c_sex,"+scale(",c_edu,")+APOE_e4+",c_pl,sep = ""))
  }
  mdl_plasma=lm(fml_plasma,df)
  mdl_sum_plasma=summary(mdl_plasma)
  testPlasma0=anova(mdl_0,mdl_plasma,test="F")
  ci_R2=boot.ci(results_boot,type="norm",index=3)
  ci_B=confint(mdl_plasma)
  c=rsq::rsq.partial(mdl_plasma,adj = T)    
  
  results_R2$Biomarker[1]="Plasma"
  results_R2$name[1]="plasma p-tau217"
  results_R2$R2[1]=paste(round(mdl_sum_plasma$adj.r.squared,2),"[",round(ci_R2$normal[2],2),", ",round(ci_R2$normal[3],2),"]",sep = "")
  results_R2$AICc[1]=paste(round(MuMIn::AICc(mdl_plasma),1))
  results_R2$Ftest_0[1]=paste(round(testPlasma0$F[2],2))
  results_R2$pComp_0[1]=paste(round(testPlasma0$`Pr(>F)`[2],3))
  results_R2$B_plasma[1]=paste(format(round(mdl_sum_plasma$coefficients[nrow(mdl_sum_plasma$coefficients),1],2),nsmall=2)," [",format(round(ci_B[nrow(mdl_sum_plasma$coefficients),1],2),nsmall=2),", ",format(round(ci_B[nrow(mdl_sum_plasma$coefficients),2],2),nsmall=2),"]",sep = "")
  results_R2$p_plasma[1]=paste(format(round(mdl_sum_plasma$coefficients[nrow(mdl_sum_plasma$coefficients),4],3),nsmall=3),sep = "")
  
  results_plot$Biomarker[1]="Plasma"
  results_plot$name[1]="plasma p-tau217"
  results_plot$R2[1]=mdl_sum_plasma$adj.r.squared
  results_plot$R2_l[1]=ci_R2$normal[2]
  results_plot$R2_h[1]=ci_R2$normal[3]
  results_plot$AICc[1]=MuMIn::AICc(mdl_plasma)
  # results_plot$pComp_0[1]=testPlasma0$`Pr(>F)`[2]

  results_B$Biomarker[1]="Plasma"
  results_B$name[1]="plasma p-tau217"
  results_B$B[1]=mdl_sum_plasma$coefficients[nrow(mdl_sum_plasma$coefficients),1]
  results_B$B_l[1]=ci_B[nrow(mdl_sum_plasma$coefficients),1]
  results_B$B_h[1]=ci_B[nrow(mdl_sum_plasma$coefficients),2]
  results_B$N=nrow(df)
  
  results_partialR2$Model[1]="Plasma"
  results_partialR2$R2_plasma[1]=c$partial.rsq[length(c$partial.rsq)]
  results_partialR2$R2_covs[1]=sum(c$partial.rsq[1:(length(c$partial.rsq)-1)])
  results_partialR2$R2[1]=mdl_sum_plasma$adj.r.squared
  results_partialR2$R2_shared[1]=mdl_sum_plasma$adj.r.squared-sum(c$partial.rsq[1:(length(c$partial.rsq))])
  
    
  for (i_PET in 1:n_PET) {
    
    if (sum(l_cohort)>1) {
      fml_PET=as.formula(paste("slope_cog~scale(",c_age,")+",c_sex,"+scale(",c_edu,")+APOE_e4+",c_cohort,"+",c_PET[i_PET],"",sep = ""))
    }else{
      fml_PET=as.formula(paste("slope_cog~scale(",c_age,")+",c_sex,"+scale(",c_edu,")+APOE_e4+",c_PET[i_PET],"",sep = ""))
    }
    mdl_PET=lm(fml_PET,df)
    mdl_sum_PET=summary(mdl_PET)
    testPET0=anova(mdl_0,mdl_PET,test="F")
    testplasma=nonnest2::vuongtest(mdl_plasma,mdl_PET)
    ci_R2=boot.ci(results_boot,type="norm",index=(3+i_PET))
    ci_B=confint(mdl_PET)
    c=rsq::rsq.partial(mdl_PET,adj = T)    
    
    results_R2$Biomarker[1+i_PET]=c_PET[i_PET]
    results_R2$name[1+i_PET]=c_name[i_PET]
    results_R2$R2[1+i_PET]=paste(round(mdl_sum_PET$adj.r.squared,2),"[",round(ci_R2$normal[2],2),", ",round(ci_R2$normal[3],2),"]",sep = "")
    results_R2$AICc[1+i_PET]=paste(round(MuMIn::AICc(mdl_PET),1))
    results_R2$Ftest_0[1+i_PET]=paste(round(testPET0$F[2],2))
    results_R2$pComp_0[1+i_PET]=paste(round(testPET0$`Pr(>F)`[2],3))
    results_R2$B_PET[1+i_PET]=paste(format(round(mdl_sum_PET$coefficients[nrow(mdl_sum_PET$coefficients),1],2),nsmall=2)," [",format(round(ci_B[nrow(mdl_sum_PET$coefficients),1],2),nsmall=2),", ",format(round(ci_B[nrow(mdl_sum_PET$coefficients),2],2),nsmall=2),"]",sep = "")
    results_R2$p_PET[1+i_PET]=paste(format(round(mdl_sum_PET$coefficients[nrow(mdl_sum_PET$coefficients),4],3),nsmall=3),sep = "")
    
    results_plot$Biomarker[1+i_PET]=c_PET[i_PET]
    results_plot$name[1+i_PET]=c_name[i_PET]
    results_plot$R2[1+i_PET]=mdl_sum_PET$adj.r.squared
    results_plot$R2_l[1+i_PET]=ci_R2$normal[2]
    results_plot$R2_h[1+i_PET]=ci_R2$normal[3]
    results_plot$AICc[1+i_PET]=MuMIn::AICc(mdl_PET)
    # results_plot$pComp_0[1+i_PET]=testPET0$`Pr(>F)`[2]
    
    results_B$Biomarker[1+i_PET]=c_PET[i_PET]
    results_B$name[1+i_PET]=c_name[i_PET]
    results_B$B[1+i_PET]=mdl_sum_PET$coefficients[nrow(mdl_sum_PET$coefficients),1]
    results_B$B_l[1+i_PET]=ci_B[nrow(mdl_sum_PET$coefficients),1]
    results_B$B_h[1+i_PET]=ci_B[nrow(mdl_sum_PET$coefficients),2]
    
    results_partialR2$Model[1+i_PET]=c_PET[i_PET]
    results_partialR2$R2_PET[1+i_PET]=c$partial.rsq[length(c$partial.rsq)]
    results_partialR2$R2_covs[1+i_PET]=sum(c$partial.rsq[1:(length(c$partial.rsq)-1)])
    results_partialR2$R2[1+i_PET]=mdl_sum_PET$adj.r.squared
    results_partialR2$R2_shared[1+i_PET]=mdl_sum_PET$adj.r.squared-sum(c$partial.rsq[1:(length(c$partial.rsq))])
    
    
    if (sum(l_cohort)>1) {
      fml_mix=as.formula(paste("slope_cog~scale(",c_age,")+",c_sex,"+scale(",c_edu,")+APOE_e4+",c_cohort,"+",c_pl,"+",c_PET[i_PET],"",sep = ""))
    }else{
      fml_mix=as.formula(paste("slope_cog~scale(",c_age,")+",c_sex,"+scale(",c_edu,")+APOE_e4+",c_pl,"+",c_PET[i_PET],"",sep = ""))
      }
    mdl_mix=lm(fml_mix,df)
    mdl_sum_mix=summary(mdl_mix)
    ci_R2=boot.ci(results_boot,type="norm",index=(5+i_PET))
    ci_B=confint(mdl_mix)
    c=rsq::rsq.partial(mdl_mix,adj = T)    
    
    results_R2$Biomarker[1+n_PET+i_PET]=paste("Plasma+",c_PET[i_PET],sep = "")
    results_R2$name[1+n_PET+i_PET]=paste("Plasma+",c_name[i_PET],sep = "")
    results_R2$R2[1+n_PET+i_PET]=paste(round(mdl_sum_mix$adj.r.squared,2),"[",round(ci_R2$normal[2],2),", ",round(ci_R2$normal[3],2),"]",sep = "")
    results_R2$AICc[1+n_PET+i_PET]=paste(round(MuMIn::AICc(mdl_mix),1))
    results_R2$B_PET[1+n_PET+i_PET]=paste(format(round(mdl_sum_mix$coefficients[nrow(mdl_sum_mix$coefficients),1],2),nsmall=2)," [",format(round(ci_B[nrow(mdl_sum_mix$coefficients),1],2),nsmall=2),", ",format(round(ci_B[nrow(mdl_sum_mix$coefficients),2],2),nsmall=2),"]",sep = "")
    results_R2$B_plasma[1+n_PET+i_PET]=paste(format(round(mdl_sum_mix$coefficients[nrow(mdl_sum_mix$coefficients)-1,1],2),nsmall=2)," [",format(round(ci_B[nrow(mdl_sum_mix$coefficients)-1,1],2),nsmall=2),", ",format(round(ci_B[nrow(mdl_sum_mix$coefficients)-1,2],2),nsmall=2),"]",sep = "")
    results_R2$p_plasma[1+n_PET+i_PET]=paste(format(round(mdl_sum_mix$coefficients[nrow(mdl_sum_mix$coefficients)-1,4],3),nsmall=3),sep = "")
    results_R2$p_PET[1+n_PET+i_PET]=paste(format(round(mdl_sum_mix$coefficients[nrow(mdl_sum_mix$coefficients),4],3),nsmall=3),sep = "")
    results_R2$p_PET[1+n_PET+i_PET]=paste(format(round(mdl_sum_mix$coefficients[nrow(mdl_sum_mix$coefficients),4],3),nsmall=3),sep = "")
    
    results_plot$Biomarker[1+n_PET+i_PET]=paste("Plasma+",c_PET[i_PET],sep = "")
    results_plot$name[1+n_PET+i_PET]=paste("Plasma+",c_name[i_PET],sep = "")
    results_plot$R2[1+n_PET+i_PET]=mdl_sum_mix$adj.r.squared
    results_plot$R2_l[1+n_PET+i_PET]=ci_R2$normal[2]
    results_plot$R2_h[1+n_PET+i_PET]=ci_R2$normal[3]
    results_plot$AICc[1+n_PET+i_PET]=MuMIn::AICc(mdl_mix)
    
    results_B$Biomarker[1+n_PET+i_PET+2]=c_PET[i_PET]
    results_B$name[1+n_PET+i_PET+2]=paste("Plasma+",c_name[i_PET],sep = "")
    results_B$B[1+n_PET+i_PET+2]=mdl_sum_mix$coefficients[nrow(mdl_sum_mix$coefficients),1]
    results_B$B_l[1+n_PET+i_PET+2]=ci_B[nrow(mdl_sum_mix$coefficients),1]
    results_B$B_h[1+n_PET+i_PET+2]=ci_B[nrow(mdl_sum_mix$coefficients),2]
    
    results_B$Biomarker[1+n_PET+i_PET]="Plasma"
    results_B$name[1+n_PET+i_PET]=paste("Plasma+",c_name[i_PET],sep = "")
    results_B$B[1+n_PET+i_PET]=mdl_sum_mix$coefficients[nrow(mdl_sum_mix$coefficients)-1,1]
    results_B$B_l[1+n_PET+i_PET]=ci_B[nrow(mdl_sum_mix$coefficients)-1,1]
    results_B$B_h[1+n_PET+i_PET]=ci_B[nrow(mdl_sum_mix$coefficients)-1,2]
    
    results_partialR2$Model[1+n_PET+i_PET]=paste("Plasma+",c_PET[i_PET],sep = "")
    results_partialR2$R2_plasma[1+n_PET+i_PET]=c$partial.rsq[length(c$partial.rsq)-1]
    results_partialR2$R2_PET[1+n_PET+i_PET]=c$partial.rsq[length(c$partial.rsq)]
    results_partialR2$R2_covs[1+n_PET+i_PET]=sum(c$partial.rsq[1:(length(c$partial.rsq)-2)])
    results_partialR2$R2[1+n_PET+i_PET]=mdl_sum_mix$adj.r.squared
    results_partialR2$R2_shared[1+n_PET+i_PET]=mdl_sum_mix$adj.r.squared-sum(c$partial.rsq[1:(length(c$partial.rsq))])
    
    testmix0=anova(mdl_0,mdl_mix,test="F")
    testPlasma=anova(mdl_plasma,mdl_mix,test="F")
    testPET=anova(mdl_PET,mdl_mix,test="F")
    results_R2$Ftest_plasma[1+n_PET+i_PET]=paste(round(testPlasma$F[2],2))
    results_R2$pComp_plasma[1+n_PET+i_PET]=paste(round(testPlasma$`Pr(>F)`[2],3))
    results_R2$Ftest_PET[1+n_PET+i_PET]=paste(round(testPET$F[2],2))
    results_R2$pComp_PET[1+n_PET+i_PET]=paste(round(testPET$`Pr(>F)`[2],3))
    results_R2$Ftest_0[1+n_PET+i_PET]=paste(round(testmix0$F[2],2))
    results_R2$pComp_0[1+n_PET+i_PET]=paste(round(testmix0$`Pr(>F)`[2],3))
    
    # results_plot$pComp_plasma[1+n_PET+i_PET]=testPlasma$`Pr(>F)`[2]
    # results_plot$pComp_PET[1+n_PET+i_PET]=testPET$`Pr(>F)`[2]
    # results_plot$pComp_0[1+n_PET+i_PET]=testmix0$`Pr(>F)`[2]
    
    # df$mix_group=interaction(df$plasma_pos,df[,c_PET_pos[i_PET]])
    # levels(df$mix_group)=c("p-PET-","p+PET-","p-PET+","p+PET+")
    
    
    if (i_PET==1) {
      v_R2=c(mdl_sum_PET$adj.r.squared,mdl_sum_mix$adj.r.squared)
      v_AIC=c(MuMIn::AICc(mdl_PET),MuMIn::AICc(mdl_mix))
      
    }
  }
  
  results_R2=results_R2[c(7,6,1,3,2,5,4),]
  
  writexl::write_xlsx(results_R2,path = paste(dir_res,"/Results_LM_",s_cohort,"_",c_out[i_out],"_DX_",l_DX,"_Ab_",l_ab,"_adj_",l_adj,"_",l_out,".xlsx",sep = ""))
  writexl::write_xlsx(results_plot,path = paste(dir_res,"/Results_plot_LM_",s_cohort,"_",c_out[i_out],"_DX_",l_DX,"_Ab_",l_ab,"_adj_",l_adj,"_",l_out,".xlsx",sep = ""))
  writexl::write_xlsx(results_B,path = paste(dir_res,"/Results_B_LM_",s_cohort,"_",c_out[i_out],"_DX_",l_DX,"_Ab_",l_ab,"_adj_",l_adj,"_",l_out,".xlsx",sep = ""))
  openxlsx::write.xlsx(v_diff_R2,file = paste(dir_res,"/pVal_LM_",s_cohort,"_",c_out[i_out],"_DX_",l_DX,"_Ab_",l_ab,"_adj_",l_adj,"_",l_out,".xlsx",sep = ""),rowNames=T)
  writexl::write_xlsx(results_partialR2,path = paste(dir_res,"/Results_partialR2_LM_",s_cohort,"_",c_out[i_out],"_DX_",l_DX,"_Ab_",l_ab,"_adj_",l_adj,"_",l_out,".xlsx",sep = ""))
  
  if (sum(l_cohort)>1) {
  c_colors_bar=c("Basic_NA"="#d9d9d9",
              "Basic_A"="#737373",
             "plasma p-tau217"="#e64b35",
             "NeoT"=c_colors_PET[1],
             "ERC/Amygd"=c_colors_PET[2],
             "Plasma+NeoT"=c_colors_mix[1],
             "Plasma+ERC/Amygd"=c_colors_mix[2])

  c_colors_bar2=c(
                 "Plasma"="#e64b35",
                 "z_NeoT"=c_colors_PET[1],
                 "z_MTL"=c_colors_PET[2])
  
  plot_bar=results_plot%>%
    mutate(index = reorder(name, R2)) %>%
    ggplot( aes(x=index, y=R2))+
    geom_bar(stat="identity", aes(fill = name),width = 1) +
    # labs(title=paste("Prediction of ", c_out_name[i_out],sep=""))+
    ylab(label = bquote(R^2))+#superscript
    xlab(label="Models")+
    geom_errorbar(aes(ymin=R2_l, ymax=R2_h), width=.2,
                  position=position_dodge(.9))+
    theme_classic()+
    scale_y_continuous(expand=c(0,0), limits=c(0,0.5)) +
    geom_text(aes(label=paste(round(AICc,0)), y=R2+0.1), size=2)+
    geom_text(aes(label=paste(format(round(R2,2),nsmall=2)), y=(R2)/2),angle=90 ,size=4, color="black")+
    # geom_text(aes(label = ifelse(pComp_PET<0.05 & pComp_plasma<0.05, "*", ""),x=name,y=R2+0.05), stat="identity", vjust=0.1,size=4,face="bold")+
    theme(legend.position = "none",axis.text=element_text(size=8),
          axis.title=element_text(size=10,face="bold"),
          # theme(legend.position = "none",axis.text=element_text(size=8),
          # axis.title=element_text(size=10,face="bold"),legend.text = element_text(size=6),
          axis.line.x = element_line(size = 0.3),axis.line.y = element_line(size = 0.3),
          axis.ticks.x=element_line(size = 0.3),axis.ticks.y = element_line(size = 0.3))+
          # legend.key.size = unit(0.2, 'cm'),legend.title = element_blank())+
    scale_fill_manual(values=c_colors_bar)
  
  results_B$name=factor(results_B$name,levels = rev(c("plasma p-tau217","ERC/Amygd","NeoT",
                                                  "Plasma+ERC/Amygd","Plasma+NeoT")))
  
  results_B$Model=as.factor(c(1,1,1,3,2,3,2))
  plot_lines=ggplot(results_B,aes(y = Biomarker, x = B,group=Biomarker))+
    geom_linerange(aes(xmax = B_h, xmin = B_l,y=Biomarker), size = 0.4, color ="black", position=position_dodge(0.5))+
    geom_point(aes(color=Biomarker),position=position_dodge(0.5),size=2.5)+
    # geom_point(aes(color=name),position=position_dodge(0.5),size=0.7)+
    scale_color_manual(values=c_colors_bar2)+
    # geom_errorbar(aes(xmax = B_h, xmin = B_l,y=name, color =name), size = 0.5, width = 0.2,position=position_dodge(0.5))+
    facet_grid(Model~.,scales = "free")+
    # geom_text(aes(x=as.numeric(B_h)+0.015,label=paste(format(round(B,3), nsmall =3))), size=3, color="black",hjust=1)+
    # geom_text(aes(label = ifelse(p_npv <0.05, "*", ""),x=as.numeric(B)), stat="identity", vjust=0.1,size=3)+
    xlab("B") + ylab("")+
    geom_vline(xintercept = 0,linetype="dashed")+
    theme_classic()+theme(legend.position = "none",axis.text=element_text(size=8),
                          axis.title=element_text(size=10,face="bold"),axis.text.y = element_blank(),
                          axis.line.x = element_line(size = 0.3),axis.line.y = element_line(size = 0.3),
                          axis.ticks.x=element_line(size = 0.3),axis.ticks.y = element_line(size = 0.3))#+
    # scale_x_continuous(limits=100*v_lim)
  
  # pf2=ggarrange(plot_bar+rremove("ylab"),
  #               # labels = c("A", "B", "C","D","E"),
  #               ncol = 5, nrow = 3)
  # 
  # ggexport(pf2,filename=paste(dir_fig,"/Barplot_LM_mergedCohorts_",s_cohort,c_out[i_out],"_DX_",l_DX,"_Ab_",l_ab,"_adj_",l_adj,"_allPlots_",l_out,"_5.pdf",sep = ""))

  pf2=ggarrange(plot_bar+rremove("ylab"),plot_lines+rremove("ylab"),
                # labels = c("A", "B", "C","D","E"),
                ncol = 3, nrow = 3)
  
  ggexport(pf2,filename=paste(dir_fig,"/Barplot_LM_mergedCohorts_",s_cohort,c_out[i_out],"_DX_",l_DX,"_Ab_",l_ab,"_adj_",l_adj,"_allPlots_",l_out,"_3.pdf",sep = ""))

  results_plot2=reshape2::melt(results_partialR2[which(results_partialR2$Model=="Plasma+z_NeoT" | results_partialR2$Model=="Plasma+z_MTL"),c(1,3:6)],id="Model")
  results_plot2=results_plot2 %>% filter(!is.na("value"))
  results_plot2$variable=as.character(results_plot2$variable)
  results_plot2$variable[which(results_plot2$variable=="R2_PET" & results_plot2$Model=="Plasma+z_NeoT")]="R2_PET_N"
  results_plot2$variable=factor(results_plot2$variable,levels = rev(c("R2_covs","R2_plasma","R2_PET","R2_PET_N","R2_shared")))
  
  plot=ggplot(results_plot2, aes(fill=variable, y=value, x=Model)) + 
    geom_bar(position="stack", stat="identity",width=0.75)+
    scale_fill_manual(values=c("R2_plasma"="#e64b35",
                               "R2_PET"="#00a087ff",
                               "R2_PET_N"="#3c5488ff",
                               "R2_covs"="#d3d3d3",
                               "R2_shared"="#4dbbd5"))+
    ylab(label = bquote(R^2))+#superscript
    xlab(label="Models")+
    theme_classic()+
    # ylim(c(0,max(results_plot$R2)))+
    theme(#aspect.ratio = 3,
          legend.position = "none",
          axis.text.y=element_text(size=8),
          axis.title=element_text(size=10,face="bold"),
          axis.ticks = element_line(size=0.25),
          axis.line = element_line(size=0.25))+
    scale_y_continuous(expand = expansion(mult = c(0)))
  
  pf2=ggarrange(plot+rremove("ylab"),
                # labels = c("A", "B", "C","D","E"),
                ncol = 3, nrow = 3)
  
  ggexport(pf2,filename=paste(dir_fig,"/Partial_R2_",s_cohort,c_out[i_out],"_DX_",l_DX,"_Ab_",l_ab,"_adj_",l_adj,"_",l_out,"_CovsOnlyCohort.pdf",sep = ""))
  
  }
  }
}

