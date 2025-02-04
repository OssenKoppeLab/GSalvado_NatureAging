# 
# dataBF2=read.csv("~/Desktop/Lund/20_TauPET_plasma/Data_BF2/BF2_PlasmaPET_progression_20230907.csv")
# dataBF1=read.csv("~/Desktop/Lund/20_TauPET_plasma/Data_BF1/DataConv_BF1_plasmaTauPET_20230907.csv")
# dataAIBL=read.csv("~/Desktop/Lund/20_TauPET_plasma/Data_AIBL/AIBL_conv_PlasmaPET_20230907.csv")
# dataTRIAD=read.csv("~/Desktop/Lund/20_TauPET_plasma/Data_TRIAD/TRIAD_conv_PlasmaPET_20240318.csv")
# dataPREVENT=read.csv("~/Desktop/Lund/20_TauPET_plasma/Data_PREVENTAD/PREVENTAD_conv_PlasmaPET_20240318.csv")
# dataWRAP=read.csv("~/Desktop/Lund/20_TauPET_plasma/Data_WRAP/WRAP_PlasmaPET_progression_20230912_CSF.csv")
# dataAms=read.csv("~/Desktop/Lund/20_TauPET_plasma/Data_Amsterdam/Amsterdam_conv_PlasmaPET_20240318.csv")
# dataWU=read.csv("~/Desktop/Lund/20_TauPET_plasma/Data_WU/WU_conv_PlasmaPET_20230929.csv")
# load("~/Desktop/Lund/20_TauPET_plasma/Data_MCSA/MCSA_PlasmaPET_20240117_CSF.RData")
# dataMCSA=dataMCSA %>% filter(!duplicated(ID))
# 
# dataBF1$APOE_e4=as.factor(ifelse(dataBF1$APOE.genotyp==24 | dataBF1$APOE.genotyp==34 | dataBF1$APOE.genotyp==44,1,0))
# dataBF1$sid=paste0("BF1_",dataBF1$sid)
# 
# dataBF2$conv_AD_dich=ifelse(dataBF2$converted_dementia_dich_baseline_variable==1 & dataBF2$underlying_etiology_text_baseline_variable=="AD",1,0)
# dataBF2$conv_AD_dich[which(is.na(dataBF2$converted_dementia_dich_baseline_variable))]=NA
# dataBF2$APOE_e4=as.factor(ifelse(dataBF2$apoe==24 | dataBF2$apoe==34 | dataBF2$apoe==44,1,0))
# dataBF2$sid=paste0("BF2_",dataBF2$sid)
# 
# dataAIBL$timeRisk=dataAIBL$timeRisk*365.25
# dataAIBL$sex=as.factor(ifelse(dataAIBL$sex=="Female",1,0))
# dataAIBL$APOE_e4=as.factor(ifelse(dataAIBL$ApoE=="E4",1,0))
# dataAIBL$AIBL.ID=paste0("AIBL_",dataAIBL$AIBL.ID)
# 
# 
# dataTRIAD$timeRisk=dataTRIAD$timeRisk*365.25
# dataTRIAD$sex=as.factor(ifelse(dataTRIAD$sex=="F",1,0))
# dataTRIAD$APOE_e4=as.factor(ifelse(dataTRIAD$apoe==24 | dataTRIAD$apoe==34 | dataTRIAD$apoe==44,1,0))
# dataTRIAD$ID=paste0("TRIAD_",dataTRIAD$ID)
# 
# dataPREVENT$timeRisk=dataPREVENT$conv_time*365.25
# dataPREVENT$sex=as.factor(ifelse(dataPREVENT$Sex=="Female",1,0))
# dataPREVENT$APOE_e4=as.factor(dataPREVENT$APOE_e4)
# dataPREVENT$PSCID.x=paste0("PREV_",dataPREVENT$PSCID.x)
# 
# dataWRAP$sex=as.factor(ifelse(dataWRAP$sex=="F",1,0))
# dataWRAP$apoe=paste(dataWRAP$apoe_e1,dataWRAP$apoe_e2,sep = "")
# dataWRAP$APOE_e4=as.factor(ifelse(dataWRAP$apoe==24 | dataWRAP$apoe==34 | dataWRAP$apoe==44,1,0))
# dataWRAP$subject_id=paste0("WRAP_",dataWRAP$subject_id)
# 
# dataAms$sex=as.factor(ifelse(dataAms$Sex=="female",1,0))
# dataAms$APOE_e4=as.factor(dataAms$APOE_status)
# dataAms$timeRisk=dataAms$timeRisk*365.25
# dataAms$I_ID=paste0("Ams_",dataAms$I_ID)
# 
# dataWU$APOE_e4=as.factor(ifelse(dataWU$apoe==24 | dataWU$apoe==34 | dataWU$apoe==44,1,0))
# dataWU$sex=as.factor(ifelse(dataWU$sex=="F",1,0))
# dataWU$age_pl=as.numeric(as.Date(dataWU$Collection.Date,format="%Y-%m-%d")-as.Date(dataWU$BIRTH,format="%Y-%m-%d"))/365.25
# dataWU$age_tau=as.numeric(as.Date(dataWU$PET_Date.x,format="%Y-%m-%d")-as.Date(dataWU$BIRTH,format="%Y-%m-%d"))/365.25
# dataWU$age=rowMeans(dataWU[,c("age_pl","age_tau")])
# dataWU$id=paste0("WU_",dataWU$id)
# 
# dataMCSA$APOE=as.factor(ifelse(dataMCSA$Any_E4=="1=Yes",1,0))
# dataMCSA$timeRisk=dataMCSA$timeRisk*365.25
# dataMCSA$ID=paste0("MCSA",dataMCSA$ID)
# 
# dataBF1=dataBF1[,c("sid","Age","sex","Education","ab_pos","Conversion_MCI",#"Conversion_Dementia",#"Conversion_AD","Time_days_AD",
#                    "Time_days","z_ptau217","z_MTL","z_NeoT","APOE_e4")]
# 
# dataBF2=dataBF2[,c("sid","age","sex","education_level_years_baseline_variable","Ab_pos",
#                    "converted_MCI_dich_baseline_variable",#"converted_dementia_dich_baseline_variable","conv_AD_dich","time_Dem",
#                    "time_MCI","z_ptau217","z_MTL","z_NeoT","APOE_e4")]
# 
# dataAIBL=dataAIBL[,c("AIBL.ID","Age","sex","edu","Ab_pos","conv_MCI_AD","timeRisk","z_ptau217","z_MTL","z_NeoT","APOE_e4")]#"conv_AD","conv_AD","timeRisk",
# 
# dataTRIAD=dataTRIAD[,c("ID","age","sex","edu","Ab_pos","conv_MCI","timeRisk","z_ptau217","z_MTL","z_NeoT","APOE_e4")]
# dataPREVENT=dataPREVENT[,c("PSCID.x","age_PET","sex","Education_years","Ab_pos","conv_mci","timeRisk","z_ptau217","z_MTL","z_NeoT","APOE_e4")]
# dataWRAP=dataWRAP[,c("subject_id","age_tau_plasma","sex","ed_years","Ab_pos","conv_MCIdem","timeRisk","z_ptau217","z_MTL","z_NeoT","APOE_e4")]
# dataAms=dataAms[,c("I_ID","Age_at_tau_PET","sex","Education_years","Amyloid_status","conv_MCI","timeRisk","z_ptau217","z_MTL","z_NeoT","APOE_e4")]
# dataWU=dataWU[,c("id","age","sex","EDUC","CL20","conv_MCI","timeRisk","z_ptau217","z_MTL","z_NeoT","APOE_e4")]
# dataMCSA=dataMCSA[,c("ID","age_bl","sex","EDUC","Ab_pos","incidentmci_before90","timeRisk","z_ptau217","z_MTL","z_NeoT","APOE")]
# 
# colnames(dataBF1)=c("sid","age","sex","edu","ab_pos","conv_MCI","time_MCI","z_ptau217","z_MTL","z_NeoT","APOE_e4")#"conv_dem","conv_AD","time_AD",
# colnames(dataBF2)=c("sid","age","sex","edu","ab_pos","conv_MCI","time_MCI","z_ptau217","z_MTL","z_NeoT","APOE_e4")#"conv_dem","conv_AD","time_AD",
# colnames(dataAIBL)=c("sid","age","sex","edu","ab_pos","conv_MCI","time_MCI","z_ptau217","z_MTL","z_NeoT","APOE_e4")#"conv_dem","conv_AD","time_AD",
# colnames(dataTRIAD)=c("sid","age","sex","edu","ab_pos","conv_MCI","time_MCI","z_ptau217","z_MTL","z_NeoT","APOE_e4")#"conv_dem","conv_AD","time_AD",
# colnames(dataPREVENT)=c("sid","age","sex","edu","ab_pos","conv_MCI","time_MCI","z_ptau217","z_MTL","z_NeoT","APOE_e4")#"conv_dem","conv_AD","time_AD",
# colnames(dataWRAP)=c("sid","age","sex","edu","ab_pos","conv_MCI","time_MCI","z_ptau217","z_MTL","z_NeoT","APOE_e4")#"conv_dem","conv_AD","time_AD",
# colnames(dataAms)=c("sid","age","sex","edu","ab_pos","conv_MCI","time_MCI","z_ptau217","z_MTL","z_NeoT","APOE_e4")#"conv_dem","conv_AD","time_AD",
# colnames(dataWU)=c("sid","age","sex","edu","ab_pos","conv_MCI","time_MCI","z_ptau217","z_MTL","z_NeoT","APOE_e4")#"conv_dem","conv_AD","time_AD",
# colnames(dataMCSA)=c("sid","age","sex","edu","ab_pos","conv_MCI","time_MCI","z_ptau217","z_MTL","z_NeoT","APOE_e4")#"conv_dem","conv_AD","time_AD",
# 
# dataBF1$cohort="BF1"
# dataBF2$cohort="BF2"
# dataAIBL$cohort="AIBL"
# dataTRIAD$cohort="TRIAD"
# dataPREVENT$cohort="PREVENT"
# dataWRAP$cohort="WRAP"
# dataAms$cohort="Ams"
# dataWU$cohort="WU"
# dataMCSA$cohort="MCSA"
# 
# dataAll=rbind(dataBF2,dataBF1,dataAIBL,dataTRIAD,dataPREVENT,dataWRAP,dataAms,dataWU,dataMCSA)
# # dataAll=rbind(dataAll)
# 
# dataAll=dataAll %>% filter(age>=40)
# 
# # c_var=c("z_MTL","z_NeoT","z_ptau217")
# #
# # for(i_ROIs in 1:3){#:length(c_var)
# #   Q <- quantile(as.numeric(unlist(dataAll[,c_var[i_ROIs]])), probs=c(.25, .75), na.rm = TRUE)
# #   iqr <- IQR(as.numeric(unlist(dataAll[,c_var[i_ROIs]])),na.rm = TRUE)
# #   a=subset(dataAll, dataAll[,c_var[i_ROIs]] < (Q[1] - 5*iqr))
# #   b=subset(dataAll,dataAll[,c_var[i_ROIs]] > (Q[2]+5*iqr))
# #   c=c(as.numeric(unlist(a[,c_var[i_ROIs]])),as.numeric(unlist(b[,c_var[i_ROIs]])))
# #
# #   d=which(as.numeric(unlist(dataAll[,c_var[i_ROIs]])) %in% c)
# #
# #   list_ID=dataAll$sid[d]
# #   dataAll[which(dataAll$sid %in% list_ID),c_var[i_ROIs]]=NA
# # }
# 
# 
# dataAll=dataAll %>% filter(!is.na(z_ptau217) & !is.na(z_MTL) & !is.na(z_NeoT))
# dataAll=dataAll %>% filter(!duplicated(sid))
# 
# dataAll$time_MCI=dataAll$time_MCI/365.25
# # table1::table1(~age+sex+edu+as.factor(ab_pos)+z_ptau217+z_MTL+z_NeoT+mmse_score+mpacc+timeDiff|cohort, dataAll)
# dataAll=dataAll %>% filter(!is.na(conv_MCI))
# 
# write.csv(dataAll,"~/Desktop/Lund/20_TauPET_plasma/Data_all/Data_all_conv_20240117.csv")
# save(dataAll,file="~/Desktop/Lund/20_TauPET_plasma/Data_all/Data_all_conv_20240117.RData")
# ##### Kaplan-Meier
# load("~/Desktop/Lund/20_TauPET_plasma/Data_all/Data_all_cogn_comb_20240117.Rdata")
# data_comb=dataAll
# data_comb=data_comb %>% filter(!duplicated(sid))
# load("~/Desktop/Lund/20_TauPET_plasma/Data_all/Data_all_conv_20240117.RData")
# # dataAll=dataAll %>% filter(cohort!="WU")
# 
# dataAll=merge(dataAll,data_comb[,c("sid","z_ptau217_comb","z_MTL_comb","z_NeoT_comb")],by="sid")
# save(dataAll,file = "~/Desktop/Lund/20_TauPET_plasma/Data_all/Data_all_conv_20240117.RData")

library(survival)
library(survminer)
library(nonnestcox)
library(boot)

dir_res="~/Desktop/Lund/20_TauPET_plasma/Results_Merged/20240917_Rev/Conversion"
dir_fig="~/Desktop/Lund/20_TauPET_plasma/Figures_Merged/20240917_Rev/Conversion"
load("~/Desktop/Lund/20_TauPET_plasma/Data_all/Data_all_cogn_comb_20240917.RData")


f_R2=function(data,indices,s_cohort0){
  
  d=data[indices,]
 
  if (s_cohort0!="allCohorts") {
    fml_0=as.formula(paste("Surv(time_MCI, conv_MCI)~ scale(age)+sex+APOE_e4+scale(edu)"))
    fml_00=as.formula(paste("Surv(time_MCI, conv_MCI)~ scale(age)+sex+scale(edu)"))
    fml_plasma=as.formula(paste("Surv(time_MCI, conv_MCI) ~ scale(age)+sex+APOE_e4+scale(edu)+z_ptau217"))
    fml_MTL=as.formula(paste("Surv(time_MCI, conv_MCI)~ scale(age)+sex+APOE_e4+scale(edu)+z_MTL"))
    fml_mix_MTL=as.formula(paste("Surv(time_MCI, conv_MCI)~ scale(age)+sex+APOE_e4+scale(edu)+z_ptau217+z_MTL"))
    fml_NeoT=as.formula(paste("Surv(time_MCI, conv_MCI)~ scale(age)+sex+APOE_e4+scale(edu)+z_NeoT"))
    fml_mix_NeoT=as.formula(paste("Surv(time_MCI, conv_MCI)~ scale(age)+sex+APOE_e4+scale(edu)+z_ptau217+z_NeoT"))
  }else{
  fml_0=as.formula(paste("Surv(time_MCI, conv_MCI)~ scale(age)+sex+APOE_e4+scale(edu)+cohort"))
  fml_00=as.formula(paste("Surv(time_MCI, conv_MCI)~ scale(age)+sex+scale(edu)+cohort"))
  fml_plasma=as.formula(paste("Surv(time_MCI, conv_MCI) ~ scale(age)+sex+APOE_e4+scale(edu)+cohort+z_ptau217"))
  fml_MTL=as.formula(paste("Surv(time_MCI, conv_MCI)~ scale(age)+sex+APOE_e4+scale(edu)+cohort+z_MTL"))
  fml_mix_MTL=as.formula(paste("Surv(time_MCI, conv_MCI)~ scale(age)+sex+APOE_e4+scale(edu)+cohort+z_ptau217+z_MTL"))
  fml_NeoT=as.formula(paste("Surv(time_MCI, conv_MCI)~ scale(age)+sex+APOE_e4+scale(edu)+cohort+z_NeoT"))
  fml_mix_NeoT=as.formula(paste("Surv(time_MCI, conv_MCI)~ scale(age)+sex+APOE_e4+scale(edu)+cohort+z_ptau217+z_NeoT"))
  }
  
  mdl_0=coxph(fml_0,d)
  mdl_00=coxph(fml_00,d)
  mdl_plasma=coxph(fml_plasma,d)
  mdl_MTL=coxph(fml_MTL,d)
  mdl_NeoT=coxph(fml_NeoT,d)
  mdl_mix_MTL=coxph(fml_mix_MTL,d)
  mdl_mix_NeoT=coxph(fml_mix_NeoT,d)
  
  mdl_sum_0=summary(mdl_0)
  mdl_sum_00=summary(mdl_00)
  mdl_sum_plasma=summary(mdl_plasma)
  mdl_sum_MTL=summary(mdl_MTL)
  mdl_sum_NeoT=summary(mdl_NeoT)
  mdl_sum_mix_MTL=summary(mdl_mix_MTL)
  mdl_sum_mix_NeoT=summary(mdl_mix_NeoT)
  
  l_mdl=list(mdl_0,mdl_00,mdl_plasma,mdl_NeoT,mdl_MTL,mdl_mix_NeoT,mdl_mix_MTL)
  v_diff_R2=NA
  i_count=1
  for (i in 1:7) {
    
    for (j in (i):7) {
      
      v_diff_R2[i_count]=-(AIC(l_mdl[[i]])-AIC(l_mdl[[j]])) #Change sign because lower is better
      i_count=i_count+1
    }
    
  }
  
  return(c(AIC(mdl_0),AIC(mdl_00),AIC(mdl_plasma),AIC(mdl_NeoT),AIC(mdl_MTL),
           AIC(mdl_mix_NeoT),AIC(mdl_mix_MTL),
           mdl_sum_MTL$conf.int[nrow(mdl_sum_MTL$conf.int),1]-mdl_sum_plasma$conf.int[nrow(mdl_sum_plasma$conf.int),1],
           mdl_sum_plasma$conf.int[nrow(mdl_sum_plasma$conf.int),1]-mdl_sum_NeoT$conf.int[nrow(mdl_sum_NeoT$conf.int),1],
           mdl_sum_MTL$conf.int[nrow(mdl_sum_MTL$conf.int),1]-mdl_sum_NeoT$conf.int[nrow(mdl_sum_NeoT$conf.int),1],v_diff_R2))
}


dataAll=dataAll %>% filter(!duplicated(sid))
dataBl=dataAll %>% filter(!duplicated(sid))
n_boot=1000

dataAll$time_MCI[which(dataAll$cohort=="WU" | dataAll$cohort=="WRAP")]=dataAll$time_MCI[which(dataAll$cohort=="WU" | dataAll$cohort=="WRAP")]*365.25
dataAll=dataAll %>% filter(!is.na(time_MCI))

v_time_min=NA#4#5# #change Max time of progression
dataBl=dataAll %>% filter(!duplicated(sid))
if (!is.na(v_time_min)) {
  dataAll=dataAll %>% filter(timeDiff>=v_time_min)
}

# median(dataBl$z_ptau217_comb)

# dataAll$plasma_pos=as.factor(ifelse(dataAll$z_ptau217_comb>median(dataBl$z_ptau217_comb),1,0))#1.96,1,0))
# dataAll$MTL_pos=as.factor(ifelse(dataAll$z_MTL_comb>median(dataBl$z_MTL_comb),1,0))#1.96,1,0))
# dataAll$NeoT_pos=as.factor(ifelse(dataAll$z_NeoT_comb>median(dataBl$z_NeoT_comb),1,0))#1.96,1,0))
# dataAll$cohort=as.factor(dataAll$cohort)

c_colors_PET=c("#3c5488ff","#00a087ff")#c("#ffa600","#bc5090")#c("#31a354","#a1d99b")
c_colors_mix=c("#8491b4","#91d1c2ff")#c("#ff6361","#58508d")

l_ab="Abpos"##"all"#
l_DX="allEtiol"
l_out="Q4"#"noOutlier"
l_adj="raw"#"combat"#

# l_cohort=c(1,1,1,1,1)#c(0,0,1,0,0)# 1: BF2, 2: BF1, 3: AIBL, 4: TRIAD, 5: PREVENT
name_cohort=c("BF2","BF1","AIBL","TRIAD","PREVENT","WRAP","Ams","WU","MCSA")
n_cohort=length(name_cohort)

# v_PET=c(15,14)#c(244:247)
# c_PET=names(dataAll)[v_PET]
# n_PET=length(v_PET)
# c_name=c("NeoT","ERC/Amygd")
# v_PET_pos=c(18,17)
# c_PET_pos=names(dataAll)[v_PET_pos]


for (i_cohort in (n_cohort+1)) {##1:(n_cohort)
  
  results_plot=data.frame(matrix(nrow = 7,ncol = 19))
  colnames(results_plot)=c("Biomarker","Model","N","Nnp","Np","HR_plasma","HR_plasma_l","HR_plasma_h","p_plasma","HR_PET","HR_PET_l","HR_PET_h","p_PET","c-index","AIC","AIC_l","AIC_h","p_modelPET","p_modelPlasma")

  results_plot2=data.frame(matrix(nrow = 7,ncol = 5))
  colnames(results_plot2)=c("Biomarker","Model","HR","HR_l","HR_h")
  
  results_B=data.frame(matrix(nrow = 7,ncol = 13))
  colnames(results_B)=c("Biomarker","Model","N","Nnp","Np","HR_plasma","p_plasma","HR_PET","p_PET","c-index","AIC","p_modelPET","p_modelPlasma")
  
  v_diff_R2=data.frame(matrix(nrow = 7,ncol = 7))
  colnames(v_diff_R2)=c("A","NA","Plasma","NeoT","MTL","P-NeoT","P-MTL")
  rownames(v_diff_R2)=c("A","NA","Plasma","NeoT","MTL","P-NeoT","P-MTL")
  
    
  df=dataAll[complete.cases(dataAll[,c("conv_MCI","z_ptau217","z_MTL","z_NeoT","APOE_e4","sex","age")]),]
  
  if (l_ab=="Abpos") {
    df=df %>% filter(ab_pos==1)
  }
  
  if (l_adj=="combat") {
    v_PET=c(15,14)+3#c(244:247)
    v_pl=13+3
    df$plasma_pos=as.factor(ifelse(df$z_ptau217_comb>quantile(dataBl$z_ptau217_comb,3/4),1,0))#1.96,1,0))#1.96,1,0))
    df$MTL_pos=as.factor(ifelse(df$z_MTL_comb>quantile(dataBl$z_MTL_comb,3/4),1,0))#1.96,1,0))#1.96,1,0))
    df$NeoT_pos=as.factor(ifelse(df$z_NeoT_comb>quantile(dataBl$z_NeoT_comb,3/4),1,0))#1.96,1,0))#1.96,1,0))
    
  }else{
    v_PET=c(10,9)-1#c(244:247)
    v_pl=8+2
    df$plasma_pos=as.factor(ifelse(df$z_ptau217>quantile(dataBl$z_ptau217,3/4),1,0))#1.96,1,0))#
    df$MTL_pos=as.factor(ifelse(df$z_MTL>quantile(dataBl$z_MTL,3/4),1,0))#1.96,1,0))#1.96,1,0))
    df$NeoT_pos=as.factor(ifelse(df$z_NeoT>quantile(dataBl$z_NeoT,3/4),1,0))#1.96,1,0))#1.96,1,0))
  }
  c_pl=names(df)[v_pl]
  c_PET=names(df)[v_PET]
  n_PET=length(v_PET)
  v_PET_pos=c(18,17)+7
  c_PET_pos=c("NeoT_pos","MTL_pos")#names(df)[v_PET_pos]
  
  
  if (i_cohort<(n_cohort+1)) {
    l_cohort=rep(0,5)
    l_cohort[i_cohort]=1
    s_cohort=name_cohort[i_cohort]
    df=df %>% filter(cohort==name_cohort[i_cohort])
    fml_ptau=as.formula(paste("Surv(time_MCI, conv_MCI) ~ scale(age)+sex+APOE_e4+scale(edu)+",c_pl))
    
  }else{
    l_cohort=rep(1,5)
    s_cohort="allCohorts"
    fml_ptau=as.formula(paste("Surv(time_MCI, conv_MCI) ~ scale(age)+sex+cohort+APOE_e4+scale(edu)+",c_pl))
  }
  
  if (nrow(df)>=50) {
    
  set.seed(1234)
  results_boot=boot(df,f_R2,n_boot,strata = df$conv_MCI,s_cohort0=s_cohort)
  
  
  v_diff_R2=data.frame(matrix(nrow = 7,ncol = 7))
  colnames(v_diff_R2)=c("A","NA","Plasma","NeoT","MTL","P-NeoT","P-MTL")
  rownames(v_diff_R2)=c("A","NA","Plasma","NeoT","MTL","P-NeoT","P-MTL")
  i_count=8
  
  v_diff_HR=data.frame(matrix(nrow = 3,ncol = 2))
  colnames(v_diff_HR)=c("Difference","p_val")
  v_diff_HR$Difference=c("MTL-Plasma","Plasma-NeoT","MTL-NeoT")
  for (i2 in 1:3) {
    results_under_H0 <- results_boot$t[,i_count] - mean(results_boot$t[,i_count])
    boot_pvalue <- mean(abs(results_under_H0) >= abs(results_boot$t0[i_count]))
    v_diff_HR[i2,2]=boot_pvalue
    i_count=i_count+1
  }
  
  for (i1 in 1:7) {
    
    for (j1 in (i1):7) {
      
      # boot_ci_norm <- boot.ci(results_boot, type = "norm",seed=12345,index = i_count)
      results_under_H0 <- results_boot$t[,i_count] - mean(results_boot$t[,i_count])
      boot_pvalue <- mean(abs(results_under_H0) >= abs(results_boot$t0[i_count]))
      v_diff_R2[i1,j1]=boot_pvalue
      i_count=i_count+1
    }
  }
  }
   
  mdl_ptau=coxph(fml_ptau, data = df)
  mdl_sum=summary(mdl_ptau)
  ci_R2=boot.ci(results_boot,type="norm",index=3)
  
  results_B$N=nrow(df)
  results_B$Nnp=sum(df$conv_MCI==0)
  results_B$Np=sum(df$conv_MCI==1)
  results_plot$N=nrow(df)
  results_plot$Nnp=sum(df$conv_MCI==0)
  results_plot$Np=sum(df$conv_MCI==1)
  

  results_B$Biomarker[1]="Plasma"
  results_B$Model[1]="Plasma"
  results_B$AIC[1]=paste(round(AIC(mdl_ptau),0))#," [",round(ci_R2$normal[2],0),", ",round(ci_R2$normal[3],0),"]",sep = "")
  results_B$p_plasma[1]=paste(format(round(mdl_sum$coefficients[nrow(mdl_sum$conf.int),5],3),nsmall=3))
  results_B$`c-index`[1]=paste(format(round(mdl_sum[["concordance"]][["C"]],2),nsmall=2))
  if (nrow(df)<50) {
    results_B$HR_plasma[1]=paste(format(round(mdl_sum$conf.int[nrow(mdl_sum$conf.int),1],2),nsmall=2),sep = "")
  }else{
    results_B$HR_plasma[1]=paste(format(round(mdl_sum$conf.int[nrow(mdl_sum$conf.int),1],2),nsmall=2)," [",format(round(mdl_sum$conf.int[nrow(mdl_sum$conf.int),3],2),nsmall=2),", ",format(round(mdl_sum$conf.int[nrow(mdl_sum$conf.int),4],2),nsmall=2),"]",sep = "")
  }
  
  results_plot$Biomarker[1]="Plasma"
  results_plot$Model[1]="Plasma"
  results_plot$HR_plasma[1]=mdl_sum$conf.int[nrow(mdl_sum$conf.int),1]
  if (nrow(df)>=50) {
  results_plot$HR_plasma_l[1]=mdl_sum$conf.int[nrow(mdl_sum$conf.int),3]
  results_plot$HR_plasma_h[1]=mdl_sum$conf.int[nrow(mdl_sum$conf.int),4]
  results_plot$p_plasma[1]=mdl_sum$coefficients[nrow(mdl_sum$conf.int),5]
  results_plot$AIC_l[1]=ci_R2$normal[2]
  results_plot$AIC_h[1]=ci_R2$normal[3]
  }
  results_plot$`c-index`[1]=mdl_sum[["concordance"]][["C"]]
  results_plot$AIC[1]=AIC(mdl_ptau)
  
  
  results_plot2$Biomarker[1]="Plasma"
  results_plot2$Model[1]="Simple"
  results_plot2$HR[1]=mdl_sum$conf.int[nrow(mdl_sum$conf.int),1]
  if (nrow(df)>=50) {
  results_plot2$HR_l[1]=mdl_sum$conf.int[nrow(mdl_sum$conf.int),3]
  results_plot2$HR_h[1]=mdl_sum$conf.int[nrow(mdl_sum$conf.int),4]
  }
  
  if (i_cohort==(length(name_cohort)+1)) {
    
  f_ptau <- survfit(Surv(time_MCI, conv_MCI) ~ plasma_pos, data = df,na.action = na.exclude)
  
  fig_ptau=ggsurvplot(
    fit = f_ptau, 
    xlab = "Years", 
    ylab = "MCI survival probability",
    size = 1,                 # change line size
    palette =  c("lightgrey","#e64b35ff"),# custom color palettes
    conf.int = TRUE,          # Add confidence interval
    pval = F,              # Add p-value
    risk.table = TRUE,        # Add risk table
    surv.median.line = "hv",
    # risk.table.col = "strata",# Risk table color by groups
    # legend.labs =
    #   c("< mean p-tau/ab42", "> mean p-tau/ab42"),    # Change legend labels
    risk.table.y.text = FALSE,# show bars instead of names in text annotations
    risk.table.height = 0.20, # Useful to change when you have multiple groups
    ggtheme = theme_classic(),xlim=c(0,6)  )
  
  # ggexport(fig_ptau$plot,filename = paste(dir_fig,"/Kaplan_MCI_",s_cohort,"_ptau_Ab_",l_ab,"_20240318.pdf",sep = ""),width = 8.27, height = 8.27,useDingbats = FALSE)
  
  }
  
  if (i_cohort<(n_cohort+1)) {
  fml_0=as.formula(paste("Surv(time_MCI, conv_MCI)~ scale(age)+sex+APOE_e4+scale(edu)"))
  }else{
    fml_0=as.formula(paste("Surv(time_MCI, conv_MCI)~ scale(age)+sex+APOE_e4+scale(edu)+cohort"))
  }
  mdl_0=coxph(fml_0, data = df)
  mdl_sum=summary(mdl_0)
  ci_R2=boot.ci(results_boot,type="norm",index=1)
  
  results_B$Biomarker[6]="Basic_APOE"
  results_B$Model[6]="Basic_APOE"
  results_B$`c-index`[6]=paste(format(round(mdl_sum[["concordance"]][["C"]],2),nsmall=2))
  results_B$AIC[6]=paste(round(AIC(mdl_0),0))#," [",round(ci_R2$normal[2],0),", ",round(ci_R2$normal[3],0),"]",sep = "")
  
  results_plot$Biomarker[6]="Basic_APOE"
  results_plot$Model[6]="Basic_A"
  results_plot$`c-index`[6]=mdl_sum[["concordance"]][["C"]]
  results_plot$AIC[6]=AIC(mdl_0)
  if (nrow(df)>=50) {
  results_plot$AIC_l[6]=ci_R2$normal[2]
  results_plot$AIC_h[6]=ci_R2$normal[3]
  }
  
  if (i_cohort<(n_cohort+1)) {
  fml_00=as.formula(paste("Surv(time_MCI, conv_MCI)~ scale(age)+sex+scale(edu)"))
  }else{
    fml_00=as.formula(paste("Surv(time_MCI, conv_MCI)~ scale(age)+sex+scale(edu)+cohort"))
  }
  mdl_00=coxph(fml_00, data = df)
  mdl_sum=summary(mdl_00)
  ci_R2=boot.ci(results_boot,type="norm",index=2)
  
  results_B$Biomarker[7]="Basic_noAPOE"
  results_B$Model[7]="Basic_noAPOE"
  results_B$`c-index`[7]=paste(format(round(mdl_sum[["concordance"]][["C"]],2),nsmall=2))
  results_B$AIC[7]=paste(round(AIC(mdl_00),0))#," [",round(ci_R2$normal[2],0),", ",round(ci_R2$normal[3],0),"]",sep = "")
  
  results_plot$Biomarker[7]="Basic_noAPOE"
  results_plot$Model[7]="Basic_NA"
  results_plot$`c-index`[7]=mdl_sum[["concordance"]][["C"]]
  results_plot$AIC[7]=AIC(mdl_00)
  if (nrow(df)>=50) {
  results_plot$AIC_l[7]=ci_R2$normal[2]
  results_plot$AIC_h[7]=ci_R2$normal[3]
  }
  
  i_pos=2
  for (i_PET in 1:n_PET) {
    
    if (i_cohort<(n_cohort+1)) {
      fml_PET=as.formula(paste("Surv(time_MCI, conv_MCI)~ scale(age)+sex+APOE_e4+scale(edu)+",c_PET[i_PET]))
      fml_mix=as.formula(paste("Surv(time_MCI, conv_MCI)~ scale(age)+sex+APOE_e4+scale(edu)+",c_pl,"+",c_PET[i_PET]))
    }else{
      fml_PET=as.formula(paste("Surv(time_MCI, conv_MCI)~ scale(age)+sex+APOE_e4+scale(edu)+ cohort+",c_PET[i_PET]))
      fml_mix=as.formula(paste("Surv(time_MCI, conv_MCI)~ scale(age)+sex+APOE_e4+scale(edu)+cohort+",c_pl,"+",c_PET[i_PET]))
    }
    
    mdl_PET=coxph(fml_PET, data = df)
    mdl_sum=summary(mdl_PET)
    ci_R2=boot.ci(results_boot,type="norm",index=(3+i_PET))
    
    results_B$Biomarker[1+i_PET]=c_PET[i_PET]
    results_B$Model[1+i_PET]=c_PET[i_PET]
    results_B$p_PET[1+i_PET]=paste(format(round(mdl_sum$coefficients[nrow(mdl_sum$conf.int),5],3),nsmall=3))
    results_B$`c-index`[1+i_PET]=paste(format(round(mdl_sum[["concordance"]][["C"]],2),nsmall=2))
    results_B$AIC[1+i_PET]=paste(round(AIC(mdl_PET),0))#," [",round(ci_R2$normal[2],0),", ",round(ci_R2$normal[3],0),"]",sep = "")
    if (nrow(df)<50) {
      results_B$HR_PET[1+i_PET]=paste(format(round(mdl_sum$conf.int[nrow(mdl_sum$conf.int),1],2),nsmall=2),sep = "")
    }else{
      results_B$HR_PET[1+i_PET]=paste(format(round(mdl_sum$conf.int[nrow(mdl_sum$conf.int),1],2),nsmall=2)," [",format(round(mdl_sum$conf.int[nrow(mdl_sum$conf.int),3],2),nsmall=2),", ",format(round(mdl_sum$conf.int[nrow(mdl_sum$conf.int),4],2),nsmall=2),"]",sep = "")
    }
    
    results_plot$Biomarker[1+i_PET]=c_PET[i_PET]
    results_plot$Model[1+i_PET]=c_PET[i_PET]
    results_plot$HR_PET[1+i_PET]=mdl_sum$conf.int[nrow(mdl_sum$conf.int),1]
    if (nrow(df)>=50) {
    results_plot$HR_PET_l[1+i_PET]=mdl_sum$conf.int[nrow(mdl_sum$conf.int),3]
    results_plot$HR_PET_h[1+i_PET]=mdl_sum$conf.int[nrow(mdl_sum$conf.int),4]
    results_plot$AIC_l[1+i_PET]=ci_R2$normal[2]
    results_plot$AIC_h[1+i_PET]=ci_R2$normal[3]
    }
    results_plot$p_PET[1+i_PET]=mdl_sum$coefficients[nrow(mdl_sum$conf.int),5]
    results_plot$`c-index`[1+i_PET]=mdl_sum[["concordance"]][["C"]]
    results_plot$AIC[1+i_PET]=AIC(mdl_PET)
    
    
    results_plot2$Biomarker[i_pos]=c_PET[i_PET]
    results_plot2$Model[i_pos]="Simple"
    results_plot2$HR[i_pos]=mdl_sum$conf.int[nrow(mdl_sum$conf.int),1]
    if (nrow(df)>=50) {
    results_plot2$HR_l[i_pos]=mdl_sum$conf.int[nrow(mdl_sum$conf.int),3]
    results_plot2$HR_h[i_pos]=mdl_sum$conf.int[nrow(mdl_sum$conf.int),4]
    }
    
    i_pos=i_pos+1
    
    if (i_cohort==(length(name_cohort)+1)) {
      
    fml_PET_pos=as.formula(paste("Surv(time_MCI, conv_MCI)~ ",c_PET_pos[i_PET]))
    if (i_PET==1) {
      f_PET <- survfit(Surv(time_MCI, conv_MCI) ~ NeoT_pos, data = df,na.action = na.exclude)
    }else{
      f_PET <- survfit(Surv(time_MCI, conv_MCI) ~ MTL_pos, data = df,na.action = na.exclude)
    }

    #survfit(fml_PET_pos, data = df,na.action = na.exclude)
    
    fig_PET=ggsurvplot(
      fit = f_PET, 
      xlab = "Years", 
      ylab = "MCI survival probability",
      size = 1,                 # change line size
      palette =  c("lightgrey",c_colors_PET[i_PET]),# custom color palettes
      conf.int = TRUE,          # Add confidence interval
      pval = F,              # Add p-value
      risk.table = TRUE,        # Add risk table
      surv.median.line = "hv",
      # risk.table.col = "strata",# Risk table color by groups
      # legend.labs =
      #   c("< mean p-tau/ab42", "> mean p-tau/ab42"),    # Change legend labels
      risk.table.y.text = FALSE,# show bars instead of names in text annotations
      risk.table.height = 0.20, # Useful to change when you have multiple groups
      ggtheme = theme_classic(),xlim=c(0,6)   )
    
    # ggexport(fig_PET$plot,filename = paste(dir_fig,"/Kaplan_MCI_",s_cohort,"_",c_PET[i_PET],"_Ab_",l_ab,"_20240318.pdf",sep = ""),width = 8.27, height = 8.27,useDingbats = FALSE)
    }
    
    mdl_mix=coxph(fml_mix, data = df)
    mdl_sum=summary(mdl_mix)
    ci_R2=boot.ci(results_boot,type="norm",index=(5+i_PET))
    
    results_B$Biomarker[1+n_PET+i_PET]=c_PET[i_PET]
    results_B$Model[1+n_PET+i_PET]=paste("Plasma_",c_PET[i_PET],sep = "")
    results_B$p_plasma[1+n_PET+i_PET]=paste(format(round(mdl_sum$coefficients[nrow(mdl_sum$conf.int)-1,5],3),nsmall=3))
    results_B$p_PET[1+n_PET+i_PET]=paste(format(round(mdl_sum$coefficients[nrow(mdl_sum$conf.int),5],3),nsmall=3))
    results_B$`c-index`[1+n_PET+i_PET]=paste(format(round(mdl_sum[["concordance"]][["C"]],2),nsmall=2))
    results_B$AIC[1+n_PET+i_PET]=paste(round(AIC(mdl_mix),0))#," [",round(ci_R2$normal[2],0),", ",round(ci_R2$normal[3],0),"]",sep = "")
    if (nrow(df)<50) {
      results_B$HR_plasma[1+n_PET+i_PET]=paste(format(round(mdl_sum$conf.int[nrow(mdl_sum$conf.int)-1,1],2),nsmall=2),sep = "")
      results_B$HR_PET[1+n_PET+i_PET]=paste(format(round(mdl_sum$conf.int[nrow(mdl_sum$conf.int),1],2),nsmall=2),sep = "")
    }else{
      results_B$HR_plasma[1+n_PET+i_PET]=paste(format(round(mdl_sum$conf.int[nrow(mdl_sum$conf.int)-1,1],2),nsmall=2)," [",format(round(mdl_sum$conf.int[nrow(mdl_sum$conf.int)-1,3],2),nsmall=2),", ",format(round(mdl_sum$conf.int[nrow(mdl_sum$conf.int)-1,4],2),nsmall=2),"]",sep = "")
      results_B$HR_PET[1+n_PET+i_PET]=paste(format(round(mdl_sum$conf.int[nrow(mdl_sum$conf.int),1],2),nsmall=2)," [",format(round(mdl_sum$conf.int[nrow(mdl_sum$conf.int),3],2),nsmall=2),", ",format(round(mdl_sum$conf.int[nrow(mdl_sum$conf.int),4],2),nsmall=2),"]",sep = "")
    }
    results_plot$Biomarker[1+n_PET+i_PET]="Plasma"
    results_plot$Model[1+n_PET+i_PET]=paste("Plasma_",c_PET[i_PET],sep = "")
    results_plot$HR_plasma[1+n_PET+i_PET]=mdl_sum$conf.int[nrow(mdl_sum$conf.int)-1,1]
    if (nrow(df)>=50) {
    results_plot$HR_plasma_l[1+n_PET+i_PET]=mdl_sum$conf.int[nrow(mdl_sum$conf.int)-1,3]
    results_plot$HR_plasma_h[1+n_PET+i_PET]=mdl_sum$conf.int[nrow(mdl_sum$conf.int)-1,4]
    results_plot$HR_PET_l[1+n_PET+i_PET]=mdl_sum$conf.int[nrow(mdl_sum$conf.int),3]
    results_plot$HR_PET_h[1+n_PET+i_PET]=mdl_sum$conf.int[nrow(mdl_sum$conf.int),4]
    results_plot$AIC_l[1+n_PET+i_PET]=ci_R2$normal[2]
    results_plot$AIC_h[1+n_PET+i_PET]=ci_R2$normal[3]
    }
    results_plot$p_plasma[1+n_PET+i_PET]=mdl_sum$coefficients[nrow(mdl_sum$conf.int)-1,5]
    results_plot$HR_PET[1+n_PET+i_PET]=mdl_sum$conf.int[nrow(mdl_sum$conf.int),1]
    results_plot$p_PET[1+n_PET+i_PET]=mdl_sum$coefficients[nrow(mdl_sum$conf.int),5]
    results_plot$`c-index`[1+n_PET+i_PET]=mdl_sum[["concordance"]][["C"]]
    results_plot$AIC[1+n_PET+i_PET]=AIC(mdl_mix)
    
    
    results_plot2$Biomarker[i_pos]=c_PET[i_PET]
    results_plot2$Model[i_pos]=paste("Plasma_",c_PET[i_PET],sep = "")
    results_plot2$HR[i_pos]=mdl_sum$conf.int[nrow(mdl_sum$conf.int),1]
    if (nrow(df)>=50) {
    results_plot2$HR_l[i_pos]=mdl_sum$conf.int[nrow(mdl_sum$conf.int),3]
    results_plot2$HR_h[i_pos]=mdl_sum$conf.int[nrow(mdl_sum$conf.int),4]
    }
    
    i_pos=i_pos+1
    
    results_plot2$Biomarker[i_pos]="Plasma"
    results_plot2$Model[i_pos]=paste("Plasma_",c_PET[i_PET],sep = "")
    results_plot2$HR[i_pos]=mdl_sum$conf.int[nrow(mdl_sum$conf.int)-1,1]
    if (nrow(df)>=50) {
    results_plot2$HR_l[i_pos]=mdl_sum$conf.int[nrow(mdl_sum$conf.int)-1,3]
    results_plot2$HR_h[i_pos]=mdl_sum$conf.int[nrow(mdl_sum$conf.int)-1,4]
    }
    i_pos=i_pos+1
    
    # df$group=interaction(df$plasma_pos,df[,c_PET_pos[i_PET]])
    # 
    # fml_mix_pos=as.formula(paste("Surv(time_MCI, conv_MCI)~ group"))
    # f_mix <- survfit(Surv(time_MCI, conv_MCI)~ group, data = df,na.action = na.exclude)
    # 
    # fig_mix=ggsurvplot(
    #   fit = f_mix, 
    #   xlab = "Years", 
    #   ylab = "MCI survival probability",
    #   size = 1,                 # change line size
    #   palette =  c("lightgrey","#3182bd",c_colors_PET[i_PET],c_colors_mix[i_PET]),# custom color palettes
    #   conf.int = TRUE,          # Add confidence interval
    #   pval = F,              # Add p-value
    #   risk.table = TRUE,        # Add risk table
    #   surv.median.line = "hv",
    #   # risk.table.col = "strata",# Risk table color by groups
    #   # legend.labs =
    #   #   c("< mean p-tau/ab42", "> mean p-tau/ab42"),    # Change legend labels
    #   risk.table.y.text = FALSE,# show bars instead of names in text annotations
    #   risk.table.height = 0.20, # Useful to change when you have multiple groups
    #   ggtheme = theme_classic()   )
    # 
    # ggexport(fig_mix$plot,filename = paste(dir_fig,"/Kaplan_MCI_",s_cohort,"_",c_PET[i_PET],"_plasma_Ab_",l_ab,"_20240318.pdf",sep = ""),width = 8.27, height = 8.27,useDingbats = FALSE)
    
    # mdl_ptau=coxph(fml_ptau, data = df,x=T)
    # mdl_PET=coxph(fml_PET, data = df,x=T)
    # mdl_NeoT=coxph(Surv(time_MCI, conv_MCI) ~ scale(age)+sex+ scale(z_NeoT)+cohort, data = dataAll,x=T)
    # mdl_MTL=coxph(Surv(time_MCI, conv_MCI) ~ scale(age)+sex+ scale(z_MTL)+cohort, data = dataAll,x=T)
    # mdl_mix=coxph(fml_mix, data = df,x=T)
    # mdl_plasma_NeoT=coxph(Surv(time_MCI, conv_MCI) ~ scale(age)+sex+ scale(z_ptau217) + scale(z_NeoT)+cohort, data = dataAll,x=T)
    # mdl_plasma_MTL=coxph(Surv(time_MCI, conv_MCI) ~ scale(age)+sex+ scale(z_ptau217) + scale(z_MTL)+cohort, data = dataAll,x=T)
    # 
    # p_comp_plasma_mix=nonnestcox::plrtest(mdl_ptau,mdl_mix,nested = T)# H1: Full model fits better than reduced model LR = 16.501,   p = 0.000137
    # p_comp_PET_mix=nonnestcox::plrtest(mdl_PET,mdl_mix,nested = T)#H1: Full model fits better than reduced model LR = 10.911,   p = 0.00272
    # # p_comp_plasma_MTL=nonnestcox::plrtest(mdl_MTL,mdl_plasma_MTL,nested = T)# H1: Full model fits better than reduced model LR = 13.455,   p = 0.000154
    # # p_comp_plasma_NeoT=nonnestcox::plrtest(mdl_NeoT,mdl_plasma_NeoT,nested = T)#H1: Full model fits better than reduced model LR = 25.065,   p = 1.99e-06
    # p_comp_mix_NeoT_MTL=nonnestcox::plrtest(mdl_plasma_NeoT,mdl_plasma_MTL,nested = F)# Are different but none better than the other!
    # nonnestcox::plrtest(mdl_NeoT,mdl_MTL,nested = F)# Are different but none better than the other!
    # 
    # results_plot$p_modelPlasma[1+n_PET+i_PET]=p_comp_plasma_mix[["pLRTAB"]]
    # results_plot$p_modelPET[1+n_PET+i_PET]=p_comp_PET_mix[["pLRTAB"]]
    # results_B$p_modelPlasma[1+n_PET+i_PET]=paste(format(round(p_comp_plasma_mix[["pLRTAB"]],3),nsmall=3))
    # results_B$p_modelPET[1+n_PET+i_PET]=paste(format(round(p_comp_PET_mix[["pLRTAB"]],3),nsmall=3))

    if (i_cohort==(length(name_cohort)+1)) {
      if (i_PET==1) {
        fig_NeoT=fig_PET
        # fig_NeoT_plasma=fig_mix
      }
    }
  }
  
  # mdl_0=coxph(Surv(time_MCI, conv_MCI) ~ scale(age)+sex+ APOE_e4+cohort+scale(edu), data = df,x=T)
  # mdl_00=coxph(Surv(time_MCI, conv_MCI) ~ scale(age)+sex+cohort+scale(edu), data = df,x=T)
  # mdl_ptau=coxph(Surv(time_MCI, conv_MCI) ~ scale(age)+sex+ APOE_e4+scale(edu)+scale(z_ptau217)+cohort, data = df,x=T)
  # mdl_NeoT=coxph(Surv(time_MCI, conv_MCI) ~ scale(age)+sex+ APOE_e4+scale(edu)+scale(z_NeoT)+cohort, data = df,x=T)
  # mdl_MTL=coxph(Surv(time_MCI, conv_MCI) ~ scale(age)+sex+ APOE_e4+scale(edu)+scale(z_MTL)+cohort, data = df,x=T)
  # mdl_plasma_NeoT=coxph(Surv(time_MCI, conv_MCI) ~ scale(age)+sex+ APOE_e4+scale(edu)+scale(z_ptau217) + scale(z_NeoT)+cohort, data = df,x=T)
  # mdl_plasma_MTL=coxph(Surv(time_MCI, conv_MCI) ~ scale(age)+sex+ APOE_e4+scale(edu)+scale(z_ptau217) + scale(z_MTL)+cohort, data = df,x=T)
  

  # v_diff_stats=data.frame(matrix(nrow = 7,ncol = 7))
  # colnames(v_diff_stats)=c("A","NA","Plasma","NeoT","MTL","P-NeoT","P-MTL")
  # rownames(v_diff_stats)=c("A","NA","Plasma","NeoT","MTL","P-NeoT","P-MTL")
  # p_comp=nonnestcox::plrtest(mdl_00,mdl_0,nested = T)
  # v_diff_stats[1,2]=p_comp[["pLRTAB"]]
  # p_comp=nonnestcox::plrtest(mdl_ptau,mdl_0,nested = T)
  # v_diff_stats[1,3]=p_comp[["pLRTAB"]]
  # p_comp=nonnestcox::plrtest(mdl_NeoT,mdl_0,nested = T)
  # v_diff_stats[1,4]=p_comp[["pLRTAB"]]
  # p_comp=nonnestcox::plrtest(mdl_MTL,mdl_0,nested = T)
  # v_diff_stats[1,5]=p_comp[["pLRTAB"]]
  # p_comp=nonnestcox::plrtest(mdl_plasma_NeoT,mdl_0,nested = T)
  # v_diff_stats[1,6]=p_comp[["pLRTAB"]]
  # p_comp=nonnestcox::plrtest(mdl_plasma_MTL,mdl_0,nested = T)
  # v_diff_stats[1,7]=p_comp[["pLRTAB"]]
  # 
  # p_comp=nonnestcox::plrtest(mdl_ptau,mdl_00,nested = T)
  # v_diff_stats[2,3]=p_comp[["pLRTAB"]]
  # p_comp=nonnestcox::plrtest(mdl_NeoT,mdl_00,nested = T)
  # v_diff_stats[2,4]=p_comp[["pLRTAB"]]
  # p_comp=nonnestcox::plrtest(mdl_MTL,mdl_00,nested = T)
  # v_diff_stats[2,5]=p_comp[["pLRTAB"]]
  # p_comp=nonnestcox::plrtest(mdl_plasma_NeoT,mdl_00,nested = T)
  # v_diff_stats[2,6]=p_comp[["pLRTAB"]]
  # p_comp=nonnestcox::plrtest(mdl_plasma_MTL,mdl_00,nested = T)
  # v_diff_stats[2,7]=p_comp[["pLRTAB"]]
  # 
  # p_comp=nonnestcox::plrtest(mdl_NeoT,mdl_ptau,nested = F)
  # v_diff_stats[3,4]=min(p_comp[["pLRTB"]],p_comp[["pLRTA"]])
  # p_comp=nonnestcox::plrtest(mdl_MTL,mdl_ptau,nested = F)
  # v_diff_stats[3,5]=min(p_comp[["pLRTB"]],p_comp[["pLRTA"]])
  # p_comp=nonnestcox::plrtest(mdl_plasma_NeoT,mdl_ptau,nested = T)
  # v_diff_stats[3,6]=p_comp[["pLRTAB"]]
  # p_comp=nonnestcox::plrtest(mdl_plasma_MTL,mdl_ptau,nested = T)
  # v_diff_stats[3,7]=p_comp[["pLRTAB"]]
  # 
  # p_comp=nonnestcox::plrtest(mdl_MTL,mdl_NeoT,nested = F)
  # v_diff_stats[4,5]=min(p_comp[["pLRTB"]],p_comp[["pLRTA"]])
  # p_comp=nonnestcox::plrtest(mdl_plasma_NeoT,mdl_NeoT,nested = T)
  # v_diff_stats[4,6]=p_comp[["pLRTAB"]]
  # p_comp=nonnestcox::plrtest(mdl_plasma_MTL,mdl_NeoT,nested = F)
  # v_diff_stats[4,7]=min(p_comp[["pLRTB"]],p_comp[["pLRTA"]])
  # 
  # p_comp=nonnestcox::plrtest(mdl_plasma_NeoT,mdl_MTL,nested = F)
  # v_diff_stats[5,6]=min(p_comp[["pLRTB"]],p_comp[["pLRTA"]])
  # p_comp=nonnestcox::plrtest(mdl_plasma_MTL,mdl_MTL,nested = T)
  # v_diff_stats[5,7]=p_comp[["pLRTAB"]]
  # 
  # p_comp=nonnestcox::plrtest(mdl_plasma_NeoT,mdl_NeoT,nested = T)
  # v_diff_stats[6,7]=p_comp[["pLRTAB"]]
  
  # openxlsx::write.xlsx(v_diff_stats,file  = paste(dir_res,"/p_stats_",s_cohort,"_Ab_",l_ab,"_adj_",l_adj,"_20240126.xlsx",sep = ""),rowNames=T)
  openxlsx::write.xlsx(v_diff_R2,file  = paste(dir_res,"/p_boot_",s_cohort,"_Ab_",l_ab,"_adj_",l_adj,"_SensTime_",v_time_min,"_20241211.xlsx",sep = ""),rowNames=T)
  openxlsx::write.xlsx(v_diff_HR,file  = paste(dir_res,"/p_boot_HR_",s_cohort,"_Ab_",l_ab,"_adj_",l_adj,"_SensTime_",v_time_min,"_20241211.xlsx",sep = ""),rowNames=T)
  
  
  if (l_adj=="combat") {
    results_plot$Model=factor(results_plot$Model,levels=c("Basic_NA","Basic_A","Plasma_z_NeoT_comb","Plasma_z_MTL_comb","z_NeoT_comb","z_MTL_comb","Plasma"))
  }else{
    results_plot$Model=factor(results_plot$Model,levels=c("Basic_NA","Basic_A","Plasma_z_NeoT","Plasma_z_MTL","z_NeoT","z_MTL","Plasma"))
  }
  
  if (i_cohort==(length(name_cohort)+1)) {
    
    c_colors=c("Plasma"="#e64b35",
               "z_NeoT"="#3c5488ff",
               "z_MTL"="#00a087ff")
    
    results_plot2$Model=factor(results_plot2$Model,levels = c("Simple","Plasma_z_MTL","Plasma_z_NeoT"))
    
    plot_plasma=results_plot2 %>% filter(!is.na(HR)) %>% #group_by(Model) %>%
      # mutate(index = reorder(Cohort_excl, -abs(HR))) %>%
      ggplot( aes(y=Biomarker, x=HR, xmin=HR_l, xmax=HR_h)) +
      facet_grid(Model~.,scales = "free")+
      # geom_vline(xintercept = results_plot$HR_plasma[results_plot$Model=="Plasma"],linetype="dashed")+
      geom_vline(xintercept = 1,linetype="dotted")+
      geom_errorbarh(height=.3,color="black",size=0.5)+
      geom_point(shape=22,aes(color=Biomarker,fill=Biomarker),size=4) +
      labs( x=paste("HR"), y = ' ') +
      theme_classic()+#ggtitle("MCI")+
      theme(axis.text.x =element_text(size=8),
            axis.text.y =element_text(size=4),
            axis.title = element_text(size=10),
            legend.position = "none",plot.title = element_text(hjust = 0.5),
            strip.text.y = element_text(size = 10))+
      # scale_color_manual(values=c(rep("black",length(name_cohort))))+
      scale_fill_manual(values=c_colors)+#(rep(length(name_cohort))))+
      scale_color_manual(values=c_colors)+#(rep(length(name_cohort))))+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            strip.background = element_rect(size=0.25, linetype="solid"))+
      coord_cartesian(xlim = c(1,1.8))
    
    
  #   results_plot2=results_plot[which(results_plot$Model=="Plasma"),c(1,2,6:8)]
  #   results_plot2$Biomarker=as.character(results_plot2$Biomarker)
  #   results_plot2$Model=as.character(results_plot2$Model)
  #   results_plot2[2,]=c("MTL","MTL",results_plot$HR_PET[which(results_plot$Model=="z_MTL")],results_plot$HR_PET_l[which(results_plot$Model=="z_MTL")],results_plot$HR_PET_h[which(results_plot$Model=="z_MTL")])
  #   results_plot2[3,]=c("NeoT","NeoT",results_plot$HR_PET[which(results_plot$Model=="z_NeoT")],results_plot$HR_PET_l[which(results_plot$Model=="z_NeoT")],results_plot$HR_PET_h[which(results_plot$Model=="z_NeoT")])
  #   results_plot2$HR_plasma=as.numeric(results_plot2$HR_plasma)
  #   results_plot2$HR_plasma_l=as.numeric(results_plot2$HR_plasma_l)
  #   results_plot2$HR_plasma_h=as.numeric(results_plot2$HR_plasma_h)
  #   results_plot2$Model=factor(results_plot2$Model,levels = c("NeoT","MTL","Plasma"))
  #   
  #   
  # plot_plasma=results_plot2 %>% filter(!is.na(HR_plasma)) %>% #group_by(Model) %>%
  #   # mutate(index = reorder(Cohort_excl, -abs(HR))) %>%
  #   ggplot( aes(y=Model, x=HR_plasma, xmin=HR_plasma_l, xmax=HR_plasma_h)) +
  #   # geom_vline(xintercept = results_plot$HR_plasma[results_plot$Model=="Plasma"],linetype="dashed")+
  #   geom_vline(xintercept = 1,linetype="dotted")+
  #   geom_errorbarh(height=.3,color="black",size=0.5)+
  #   geom_point(shape=22,aes(fill=Model,color=Model),size=4) +
  #   labs( x=paste("HR"), y = ' ') +
  #   theme_classic()+#ggtitle("MCI")+
  #   theme(axis.text.x =element_text(size=8),
  #         axis.text.y =element_text(size=4),
  #         axis.title = element_text(size=10),
  #         legend.position = "none",plot.title = element_text(hjust = 0.5),
  #         strip.text.y = element_text(size = 10))+
  #   # scale_color_manual(values=c(rep("black",length(name_cohort))))+
  #   scale_fill_manual(values=c_colors)+#(rep(length(name_cohort))))+
  #   scale_color_manual(values=c_colors)+#(rep(length(name_cohort))))+
  #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+    
  #   coord_cartesian(xlim = c(1,1.8))
  # 
  # results_plot2=results_plot[which(results_plot$Model=="Plasma_z_MTL"),c(1,2,6:8)]
  # results_plot2$Biomarker=as.character(results_plot2$Biomarker)
  # results_plot2$Model=as.character(results_plot2$Model)
  # results_plot2[2,]=c("MTL","MTL",results_plot$HR_PET[which(results_plot$Model=="Plasma_z_MTL")],results_plot$HR_PET_l[which(results_plot$Model=="Plasma_z_MTL")],results_plot$HR_PET_h[which(results_plot$Model=="Plasma_z_MTL")])
  # results_plot2$HR_plasma=as.numeric(results_plot2$HR_plasma)
  # results_plot2$HR_plasma_l=as.numeric(results_plot2$HR_plasma_l)
  # results_plot2$HR_plasma_h=as.numeric(results_plot2$HR_plasma_h)
  # results_plot2$Model[1]="Plasma"
  # results_plot2$Model=factor(results_plot2$Model,levels = c("MTL","Plasma"))
  # 
  # 
  # plot_MTL=results_plot2 %>%  #group_by(Model) %>%
  #   # mutate(index = reorder(Cohort_excl, -abs(HR))) %>%
  #   ggplot( aes(y=Model, x=HR_plasma, xmin=HR_plasma_l, xmax=HR_plasma_h)) +
  #   geom_vline(xintercept = results_plot$HR_PET[results_plot$Model=="z_MTL"],linetype="dashed",color="#00a087ff")+
  #   geom_vline(xintercept = results_plot$HR_plasma[results_plot$Model=="Plasma"],linetype="dashed",color="#e64b35")+
  #   geom_vline(xintercept = 1,linetype="dotted")+
  #   geom_errorbarh(height=.3,color="black",size=0.5)+
  #   geom_point(shape=22,aes(fill=Model,color=Model),size=4) +
  #   labs( x=paste("HR"), y = ' ') +
  #   theme_classic()+#ggtitle("MCI")+
  #   theme(axis.text.x =element_text(size=8),
  #         axis.text.y =element_text(size=4),
  #         axis.title = element_text(size=10),
  #         legend.position = "none",plot.title = element_text(hjust = 0.5),
  #         strip.text.y = element_text(size = 10))+
  #   # scale_color_manual(values=c(rep("black",length(name_cohort))))+
  #   scale_fill_manual(values=c_colors)+#rep("#bc5090",length(name_cohort))))+
  #   scale_color_manual(values=c_colors)+#rep("#bc5090",length(name_cohort))))+
  #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+    
  #   coord_cartesian(xlim = c(1,1.8))
  # 
  # results_plot2=results_plot[which(results_plot$Model=="Plasma_z_NeoT"),c(1,2,6:8)]
  # results_plot2$Biomarker=as.character(results_plot2$Biomarker)
  # results_plot2$Model=as.character(results_plot2$Model)
  # results_plot2[2,]=c("NeoT","NeoT",results_plot$HR_PET[which(results_plot$Model=="Plasma_z_NeoT")],results_plot$HR_PET_l[which(results_plot$Model=="Plasma_z_NeoT")],results_plot$HR_PET_h[which(results_plot$Model=="Plasma_z_NeoT")])
  # results_plot2$HR_plasma=as.numeric(results_plot2$HR_plasma)
  # results_plot2$HR_plasma_l=as.numeric(results_plot2$HR_plasma_l)
  # results_plot2$HR_plasma_h=as.numeric(results_plot2$HR_plasma_h)
  # results_plot2$Model[1]="Plasma"
  # results_plot2$Model=factor(results_plot2$Model,levels = c("NeoT","Plasma"))
  # 
  # 
  # plot_NeoT=results_plot2 %>%  #group_by(Model) %>%
  #   # mutate(index = reorder(Cohort_excl, -abs(HR))) %>%
  #   ggplot( aes(y=Model, x=HR_plasma, xmin=HR_plasma_l, xmax=HR_plasma_h)) +
  #   geom_vline(xintercept = results_plot$HR_PET[results_plot$Model=="z_NeoT"],linetype="dashed",color="#3c5488ff")+
  #   geom_vline(xintercept = results_plot$HR_plasma[results_plot$Model=="Plasma"],linetype="dashed",color="#e64b35")+
  #   geom_vline(xintercept = 1,linetype="dotted")+
  #   geom_errorbarh(height=.3,color="black",size=0.5)+
  #   geom_point(shape=22,aes(fill=Model,color=Model),size=4) +
  #   labs( x=paste("HR"), y = ' ') +
  #   theme_classic()+#ggtitle("MCI")+
  #   theme(axis.text.x =element_text(size=8),
  #         axis.text.y =element_text(size=4),
  #         axis.title = element_text(size=10),
  #         legend.position = "none",plot.title = element_text(hjust = 0.5),
  #         strip.text.y = element_text(size = 10))+
  #   # scale_color_manual(values=c(rep("black",length(name_cohort))))+
  #   scale_fill_manual(values=c_colors)+#rep("#bc5090",length(name_cohort))))+
  #   scale_color_manual(values=c_colors)+#rep("#bc5090",length(name_cohort))))+
  #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+    
  #   coord_cartesian(xlim = c(1,1.8))
  # 
  pf2=ggarrange(
                plot_plasma,
                # labels = c("A", "B", "C","D","E"),
                ncol = 2, nrow = 3)
  
  
  ggexport(pf2,filename=paste(dir_fig,"/Kaplan_MCI_",s_cohort,"_",c_PET[i_PET],"_plasma_allPlots_Ab_",l_ab,"_adj_",l_adj,"_SensTime_",v_time_min,"_20241211.pdf",sep = ""))
  
  if (l_adj=="combat") {
    c_colors_bar=c(
      "Basic_NA"="#d9d9d9",
      "Basic_A"="#737373",
      "Plasma"="#e64b35ff",
      "z_NeoT_comb"=c_colors_PET[1],
      "z_MTL_comb"=c_colors_PET[2],
      "Plasma_z_NeoT_comb"=c_colors_mix[1],
      "Plasma_z_MTL_comb"=c_colors_mix[2])
  }else{c_colors_bar=c(
    "Basic_NA"="#d9d9d9",
    "Basic_A"="#737373",
    "Plasma"="#e64b35ff",
    "z_NeoT"=c_colors_PET[1],
    "z_MTL"=c_colors_PET[2],
    "Plasma_z_NeoT"=c_colors_mix[1],
    "Plasma_z_MTL"=c_colors_mix[2])
  }
  
  # plot_bar=results_plot%>%
  #   mutate(index = reorder(Model, `c-index`)) %>%
  #   ggplot( aes(x=index, y=`c-index`))+
  #   geom_bar(stat="identity", aes(fill = Model)) +
  #   # labs(title=paste("Prediction of ", c_out_name[i_out],sep=""))+
  #   ylab(label = "C-index")+#superscript
  #   xlab(label="Models")+
  #   theme_classic()+
  #   scale_y_continuous(expand=c(0,0), limits=c(0,1)) +
  #   geom_text(aes(label=paste(format(round(`c-index`,2),nsmall=2)), y=`c-index`-0.025), size=8, color="black")+
  #   geom_text(aes(label = ifelse(p_modelPET<0.05 & p_modelPlasma<0.05, "*", ""),x=Model,y=`c-index`+0.05), stat="identity", vjust=0.1,size=12,face="bold")+
  #   theme(legend.position = "none",axis.text=element_text(size=10),
  #         axis.title=element_text(size=26,face="bold"))+
  #   # theme(legend.position = "none",axis.text=element_text(size=8),
  #   # axis.title=element_text(size=10,face="bold"),legend.text = element_text(size=6),
  #   # axis.line.x = element_line(size = 0.3),axis.line.y = element_line(size = 0.3),
  #   # axis.ticks.x=element_line(size = 0.3),axis.ticks.y = element_line(size = 0.3),
  #   # legend.key.size = unit(0.2, 'cm'),legend.title = element_blank())+
  #   scale_fill_manual(values=c_colors_bar)
  # 
  # ggexport(plot_bar,filename=paste(dir_fig,"/Barplot_Cstat_mergedCohorts_",s_cohort,"_PET_plasma_Ab_",l_ab,"_adj_",l_adj,"_20240318.pdf",sep = ""))

  plot_bar=results_plot%>%
    mutate(index = reorder(Model, -AIC)) %>%
    ggplot( aes(x=index, y=AIC))+
    geom_bar(stat="identity", aes(fill = Model)) +
    # labs(title=paste("Prediction of ", c_out_name[i_out],sep=""))+
    ylab(label = "AIC")+#superscript
    xlab(label="Models")+
    # geom_errorbar(aes(ymin=AIC_l, ymax=AIC_h), width=.2,
    #               position=position_dodge(.9))+
    theme_classic()+
    # scale_y_continuous(expand=c(0,0), limits=c(0,1)) +
    geom_text(aes(label=paste(format(round(AIC,0))), y=AIC/2), size=4, color="black",angle=90)+
    # geom_text(aes(label = ifelse(p_modelPET<0.05 & p_modelPlasma<0.05, "*", ""),x=Model,y=AIC+0.05), stat="identity", vjust=0.1,size=12,face="bold")+
    theme(legend.position = "none",axis.text=element_text(size=8),
          axis.title=element_text(size=10,face="bold"))+
    # theme(legend.position = "none",axis.text=element_text(size=8),
    # axis.title=element_text(size=10,face="bold"),legend.text = element_text(size=6),
    # axis.line.x = element_line(size = 0.3),axis.line.y = element_line(size = 0.3),
    # axis.ticks.x=element_line(size = 0.3),axis.ticks.y = element_line(size = 0.3),
    # legend.key.size = unit(0.2, 'cm'),legend.title = element_blank())+
    scale_fill_manual(values=c_colors_bar)
  
  pf3=ggarrange(plot_bar+rremove("ylab"),
                # labels = c("A", "B", "C","D","E"),
                ncol = 3, nrow = 3)
  
  ggexport(pf3,filename=paste(dir_fig,"/Barplot_AIC_mergedCohorts_",s_cohort,"_PET_plasma_Ab_",l_ab,"_adj_",l_adj,"_SensTime_",v_time_min,"_20241211.pdf",sep = ""))
  
    }
  
  results_B=results_B[c(7,6,1,3,2,5,4),]
  writexl::write_xlsx(results_B,path = paste(dir_res,"/Results_Kaplan_MCI_",s_cohort,"_Ab_",l_ab,"_adj_",l_adj,"_SensTime_",v_time_min,"_20241211.xlsx",sep = ""))
  writexl::write_xlsx(results_plot,path = paste(dir_res,"/Plot_Kaplan_MCI_",s_cohort,"_Ab_",l_ab,"_adj_",l_adj,"_SensTime_",v_time_min,"_20241211.xlsx",sep = ""))
  writexl::write_xlsx(results_plot2,path = paste(dir_res,"/Plot_HR_Kaplan_MCI_",s_cohort,"_Ab_",l_ab,"_adj_",l_adj,"_SensTime_",v_time_min,"_20241211.xlsx",sep = ""))
  
  
}

pf2=ggarrange(fig_ptau$plot,  fig_PET$plot , fig_NeoT$plot, 
              fig_ptau$table,  fig_PET$table , fig_NeoT$table,
              # plot_plasma,plot_MTL,plot_NeoT,
              # labels = c("A", "B", "C","D","E"),
              ncol = 3, nrow = 4,heights = c(1,0.4,0.4))


ggexport(pf2,filename=paste(dir_fig,"/Kaplan_MCI_",l_out,"_",s_cohort,"_",c_PET[i_PET],"_plasma_allPlots_Ab_",l_ab,"_adj_",l_adj,"_SensTime_",v_time_min,"_20241211.pdf",sep = ""))

