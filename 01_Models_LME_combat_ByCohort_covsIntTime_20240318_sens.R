# 
# load("~/Desktop/Lund/20_TauPET_plasma/Data_BF2/BF2_PlasmaPET_20230907_CSF.RData")
# dataBF1=read.csv("~/Desktop/Lund/20_TauPET_plasma/Data_BF1/DataLong_BF1_plasmaTauPET_20230907.csv")
# dataBF1$APOE_e4=as.factor(ifelse(dataBF1$APOE.genotyp==24 | dataBF1$APOE.genotyp==34 | dataBF1$APOE.genotyp==44,1,0))
# dataBF1$sid=paste0("BF1_",dataBF1$sid)
# 
# dataBF2=dataFinal
# dataBF2$APOE_e4=as.factor(ifelse(dataBF2$apoe==24 | dataBF2$apoe==34 | dataBF2$apoe==44,1,0))
# dataBF2$sid=paste0("BF2_",dataBF2$sid)
# 
# load("~/Desktop/Lund/20_TauPET_plasma/Data_AIBL/AIBL_PlasmaPET_20230907.RData")
# dataAIBL=dataFinal
# dataAIBL$sex=as.factor(ifelse(dataAIBL$sex=="Female",1,0))
# dataAIBL$APOE_e4=as.factor(ifelse(dataAIBL$ApoE=="E4",1,0))
# dataAIBL$time_tau_blood=abs(as.numeric(as.Date(dataAIBL$Tau.Acquisition.Date,format="%Y-%m-%d")-as.Date(dataAIBL$`Blood date`,format="%Y-%m-%d")))/365.25
# dataAIBL$`AIBL ID`=paste0("AIBL_",dataAIBL$`AIBL ID`)
# 
# dataTRIAD=read.csv("~/Desktop/Lund/20_TauPET_plasma/Data_TRIAD/TRIAD_PlasmaPET_20230927.csv")
# dataTRIAD$sex=as.factor(ifelse(dataTRIAD$sex=="F",1,0))
# dataTRIAD$APOE_e4=as.factor(ifelse(dataTRIAD$apoe==24 | dataTRIAD$apoe==34 | dataTRIAD$apoe==44,1,0))
# dataTRIAD$ID=paste0("TRIAD_",dataTRIAD$ID)
# 
# dataPREVENT=read.csv("~/Desktop/Lund/20_TauPET_plasma/Data_PREVENTAD/PREVENTAD_PlasmaPET_20230918.csv")
# dataPREVENT$sex=as.factor(ifelse(dataPREVENT$Sex=="Female",1,0))
# dataPREVENT$APOE_e4=as.factor(dataPREVENT$APOE_e4)
# dataPREVENT$time_tau_blood=abs(as.numeric(as.Date(dataPREVENT$Date_PET,format="%Y-%m-%d")-as.Date(dataPREVENT$Date_cog,format="%Y-%m-%d")))/365.25
# dataPREVENT$PSCID.x=paste0("PREV_",dataPREVENT$PSCID.x)
# 
# load("~/Desktop/Lund/20_TauPET_plasma/Data_WRAP/WRAP_PlasmaPET_20230912_CSF.RData")
# dataWRAP=data4
# dataWRAP$sex=as.factor(ifelse(dataWRAP$sex=="F",1,0))
# dataWRAP$apoe=paste(dataWRAP$apoe_e1,dataWRAP$apoe_e2,sep = "")
# dataWRAP$APOE_e4=as.factor(ifelse(dataWRAP$apoe==24 | dataWRAP$apoe==34 | dataWRAP$apoe==44,1,0))
# dataWRAP$subject_id=paste0("WRAP_",dataWRAP$subject_id)
# 
# load("~/Desktop/Lund/20_TauPET_plasma/Data_Amsterdam/Amsterdam_PlasmaPET_20230927.RData")
# dataAms=data_long
# dataAms$sex=as.factor(ifelse(dataAms$Sex=="female",1,0))
# dataAms$time_tau_plasma=NA
# dataAms$I_ID=paste0("Ams_",dataAms$I_ID)
# 
# dataWU=read.csv("~/Desktop/Lund/20_TauPET_plasma/Data_WU/DataLong_WU_plasmaTauPET_20230929.csv")
# dataWU$APOE_e4=as.factor(ifelse(dataWU$apoe==24 | dataWU$apoe==34 | dataWU$apoe==44,1,0))
# dataWU$sex=as.factor(ifelse(dataWU$sex=="F",1,0))
# dataWU$age_pl=as.numeric(as.Date(dataWU$Collection.Date,format="%Y-%m-%d")-as.Date(dataWU$BIRTH,format="%Y-%m-%d"))/365.25
# dataWU$age_tau=as.numeric(as.Date(dataWU$PET_Date.x,format="%Y-%m-%d")-as.Date(dataWU$BIRTH,format="%Y-%m-%d"))/365.25
# dataWU$age=rowMeans(dataWU[,c("age_pl","age_tau")])
# dataWU$id=paste0("WU_",dataWU$id)
# 
# load("~/Desktop/Lund/20_TauPET_plasma/Data_MCSA/MCSA_PlasmaPET_20240117_CSF.RData")
# dataMCSA$APOE=as.factor(ifelse(dataMCSA$Any_E4=="1=Yes",1,0))
# dataMCSA$ID=paste0("MCSA",dataMCSA$ID)
# 
# dataBF2=dataBF2[,c("sid","age","sex","education_level_years_baseline_variable","timeDiff",
#                    "mmse_score_FU","mpacc","z_MTL","z_NeoT","z_ptau217","Ab_pos","APOE_e4","time_tau_blood")]
# dataBF1=dataBF1[,c("sid","Age","sex","Education","timeDiff","mmse_score","mpacc","z_MTL","z_NeoT","z_ptau217","ab_pos","APOE_e4","tau_blood_diff")]
# dataAIBL=dataAIBL[,c("AIBL ID","Age","sex","edu","timeDiff","Neuropsych.MMSE_FU","mpacc","z_MTL","z_NeoT","z_ptau217","Ab_pos","APOE_e4","time_tau_blood")]
# dataTRIAD=dataTRIAD[,c("ID","age","sex","edu","timeDiff","MMSE","mpacc","z_MTL","z_NeoT","z_ptau217","Ab_pos","APOE_e4","time_tau_blood")]
# dataPREVENT=dataPREVENT[,c("PSCID.x","age_PET","sex","Education_years","time_cogn_PET","MMSE","mpacc","z_MTL","z_NeoT","z_ptau217","Ab_pos","APOE_e4","time_tau_blood")]
# dataWRAP=dataWRAP[,c("subject_id","age_tau_plasma","sex","ed_years","time_age_cogn","mmseTot","mpacc","z_MTL","z_NeoT","z_ptau217","Ab_pos","APOE_e4","time_tau_plasma")]
# dataAms=dataAms[,c("I_ID","Age_at_tau_PET","sex","Education_years","timeDiff","V_MMSE","mpacc","z_MTL","z_NeoT","z_ptau217","Amyloid_status","APOE_status","time_tau_plasma")]
# dataWU=dataWU[,c("id","age","sex","EDUC","timeDiff","MMSE","mpacc","z_MTL","z_NeoT","z_ptau217","CL20","APOE_e4","timeDiff_Tau_plasma")]
# dataMCSA=dataMCSA[,c("ID","agevis","sex","EDUC","timeDiff","MMSEcalc","mpacc","z_MTL","z_NeoT","z_ptau217","Ab_pos","APOE","duration_plasma_tau")]
# 
# colnames(dataBF2)=c("sid","age","sex","edu","timeDiff","mmse_score","mpacc","z_MTL","z_NeoT","z_ptau217","ab_pos","APOE_e4","time_tau_blood")
# colnames(dataBF1)=c("sid","age","sex","edu","timeDiff","mmse_score","mpacc","z_MTL","z_NeoT","z_ptau217","ab_pos","APOE_e4","time_tau_blood")
# colnames(dataAIBL)=c("sid","age","sex","edu","timeDiff","mmse_score","mpacc","z_MTL","z_NeoT","z_ptau217","ab_pos","APOE_e4","time_tau_blood")
# colnames(dataTRIAD)=c("sid","age","sex","edu","timeDiff","mmse_score","mpacc","z_MTL","z_NeoT","z_ptau217","ab_pos","APOE_e4","time_tau_blood")
# colnames(dataPREVENT)=c("sid","age","sex","edu","timeDiff","mmse_score","mpacc","z_MTL","z_NeoT","z_ptau217","ab_pos","APOE_e4","time_tau_blood")
# colnames(dataWRAP)=c("sid","age","sex","edu","timeDiff","mmse_score","mpacc","z_MTL","z_NeoT","z_ptau217","ab_pos","APOE_e4","time_tau_blood")
# colnames(dataAms)=c("sid","age","sex","edu","timeDiff","mmse_score","mpacc","z_MTL","z_NeoT","z_ptau217","ab_pos","APOE_e4","time_tau_blood")
# colnames(dataWU)=c("sid","age","sex","edu","timeDiff","mmse_score","mpacc","z_MTL","z_NeoT","z_ptau217","ab_pos","APOE_e4","time_tau_blood")
# colnames(dataMCSA)=c("sid","age","sex","edu","timeDiff","mmse_score","mpacc","z_MTL","z_NeoT","z_ptau217","ab_pos","APOE_e4","time_tau_blood")
# 
# dataBF2$cohort="BF2"
# dataBF1$cohort="BF1"
# dataAIBL$cohort="AIBL"
# dataTRIAD$cohort="TRIAD"
# dataPREVENT$cohort="PREVENT"
# dataWRAP$cohort="WRAP"
# dataAms$cohort="Ams"
# dataWU$cohort="WU"
# dataMCSA$cohort="MCSA"
# 
# dataBF1$tracer="FTP"
# dataBF2$tracer="RO"
# dataAIBL$tracer="MK"
# dataTRIAD$tracer="MK"
# dataPREVENT$tracer="FTP"
# dataWRAP$tracer="MK"
# dataAms$tracer="FTP" #???
# dataWU$tracer="FTP" #???
# dataMCSA$tracer="FTP" #???
# 
# dataAll=rbind(dataBF2,dataBF1,dataAIBL,dataTRIAD,dataPREVENT,dataWRAP,dataAms,dataWU,dataMCSA)
# dataAll$tracer=as.factor(dataAll$tracer)
# 
# dataAll=dataAll[with(dataAll, order(sid, -timeDiff)),]
# dataAll=dataAll %>% filter(age>=40)
# dataAll=dataAll[complete.cases(dataAll[,c("z_MTL","z_NeoT","z_ptau217")]),]
# 
# ##Neurocombat
# # library(devtools)
# # install_github("jfortin1/neuroCombatData")
# # install_github("jfortin1/neuroCombat_Rpackage")
# 
# library(neuroCombat)
# v_dat=rbind(t(dataAll$z_ptau217))
# dataNew=neuroCombat(dat=v_dat,batch = dataAll$cohort,eb = F)
# 
# v_dat_tau=rbind(t(dataAll[,c("z_MTL","z_NeoT")]))
# dataNew_tau=neuroCombat(dat=v_dat_tau,batch = dataAll$cohort,eb = F)
# 
# v_dat_pacc=rbind(t(dataAll$mpacc))
# dataNew_pacc=neuroCombat(dat=v_dat_pacc,batch = dataAll$cohort,eb = F)
# 
# dataAll$z_ptau217_comb=as.numeric(dataNew$dat.combat)
# dataAll$z_MTL_comb=as.numeric(dataNew_tau$dat.combat[1,])
# dataAll$z_NeoT_comb=as.numeric(dataNew_tau$dat.combat[2,])
# dataAll$mpacc_comb=as.numeric(dataNew_pacc$dat.combat)
# 
# dataBl=dataAll %>% filter(!duplicated(sid))
# # 
# # 
# # # c_var=c("z_MTL","z_NeoT","z_ptau217")
# # #   
# # # for(i_ROIs in 1:3){#:length(c_var)
# # #   Q <- quantile(as.numeric(unlist(dataBl[,c_var[i_ROIs]])), probs=c(.25, .75), na.rm = TRUE)
# # #   iqr <- IQR(as.numeric(unlist(dataBl[,c_var[i_ROIs]])),na.rm = TRUE)
# # #   a=subset(dataBl, dataBl[,c_var[i_ROIs]] < (Q[1] - 5*iqr))
# # #   b=subset(dataBl,dataBl[,c_var[i_ROIs]] > (Q[2]+5*iqr))
# # #   c=c(as.numeric(unlist(a[,c_var[i_ROIs]])),as.numeric(unlist(b[,c_var[i_ROIs]])))
# # #   
# # #   d=which(as.numeric(unlist(dataBl[,c_var[i_ROIs]])) %in% c)
# # #   
# # #   list_ID=dataBl$sid[d]
# # #   dataAll[which(dataAll$sid %in% list_ID),c_var[i_ROIs]]=NA
# # # }
# # 
# # 
# # dataAll=dataAll %>% filter(!is.na(z_ptau217) & !is.na(z_MTL) & !is.na(z_NeoT))
# # dataBl=dataAll %>% filter(!duplicated(sid))
# # 
# table1::table1(~age+sex+edu+APOE_e4+as.factor(ab_pos)+z_ptau217+z_ptau217_comb+z_MTL+z_MTL_comb+z_NeoT+z_NeoT_comb+mmse_score+mpacc+mpacc_comb+timeDiff|cohort, dataBl)#as.factor(conv_MCI)+
# table1::table1(~age+sex+edu+APOE_e4+as.factor(ab_pos)+z_ptau217_comb+z_MTL_comb+z_NeoT_comb+mmse_score+mpacc_comb+timeDiff|cohort, dataBl)#as.factor(conv_MCI)+
# # 
# dataBl %>%
#   ggplot( aes(x=z_NeoT, fill=as.factor(cohort))) +
#   geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
#   # scale_fill_manual(values=c("#69b3a2", "#404080")) +
#   theme_classic() +
#   labs(fill="")
# 
# dataBl %>%
#   ggplot( aes(x=z_NeoT_comb, fill=as.factor(cohort))) +
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
#   ggplot( aes(x=z_ptau217_comb, fill=as.factor(cohort))) +
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
# dataBl %>%
#   ggplot( aes(x=z_MTL_comb, fill=as.factor(cohort))) +
#   geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
#   # scale_fill_manual(values=c("#69b3a2", "#404080")) +
#   theme_classic() +
#   labs(fill="")
# 
# dataBl %>%
#   ggplot( aes(x=mpacc, fill=as.factor(cohort))) +
#   geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
#   # scale_fill_manual(values=c("#69b3a2", "#404080")) +
#   theme_classic() +
#   labs(fill="")
# 
# dataBl %>%
#   ggplot( aes(x=mpacc_comb, fill=as.factor(cohort))) +
#   geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
#   # scale_fill_manual(values=c("#69b3a2", "#404080")) +
#   theme_classic() +
#   labs(fill="")
# 
# # 
# # dataBl %>%
# #   ggplot( aes(x=mmse_score, fill=as.factor(cohort))) +
# #   geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
# #   # scale_fill_manual(values=c("#69b3a2", "#404080")) +
# #   theme_classic() +
# #   labs(fill="")
# # 
# save(dataAll,file="~/Desktop/Lund/20_TauPET_plasma/Data_all/Data_all_cogn_comb_20240117.RData")

####
dir_res="~/Desktop/Lund/20_TauPET_plasma/Results_Merged/20240917_Rev/Cogn_decline"
dir_fig="~/Desktop/Lund/20_TauPET_plasma/Figures_Merged/20240917_Rev/Cogn_decline"
load("~/Desktop/Lund/20_TauPET_plasma/Data_all/Data_all_cogn_comb_20240917.RData")
############# Analyses
library(lme4)
library(lmerTest)
library(ggeffects)
library(boot)

# f_mdl=function(fml,data,indices,v_row){
#   d=data[indices,]
#   mdl_boot=lmer(fml,d,REML=F)
#   mdl_sum_boot=summary(mdl_boot)
#   return(mdl_sum_boot$coefficients[v_row,1])
# }

dataAll=dataAll[with(dataAll, order(sid, -timeDiff)),]

# median(dataBl$z_ptau217_comb)
n_rep=500#10#

v_time_min=4#NA#4#Change maximum time!
dataBl=dataAll %>% filter(!duplicated(sid))
if (!is.na(v_time_min)) {
  dataBl=dataBl %>% filter(timeDiff>v_time_min)
}

dataAll=dataAll[with(dataAll, order(sid, timeDiff)),]
dataAll=filter(dataAll,sid %in%  dataBl$sid)

dataAll$cohort=as.factor(dataAll$cohort)

c_colors_PET=c("#3c5488ff","#00a087ff")#c("#ffa600","#bc5090")#c("#31a354","#a1d99b")
c_colors_mix=c("#ff6361","#58508d")

c_age="age"
c_sex="sex"
c_edu="edu"
c_apoe="APOE_e4"
c_ID="sid"
c_cohort="cohort"
c_tracer="tracer"

# v_out=c(19,6)#c(238,235)#251:256,
# c_out=names(dataAll)[v_out]
c_out_name=c("mPACC")#,"MMSE")
n_out=length(c_out_name)

# v_PET=c(18,17)#c(244:247)
# c_PET=names(dataAll)[v_PET]
# n_PET=length(v_PET)
c_name=c("NeoT","ERC/Amygd")

l_ab="all"#"Abpos"#
l_DX="allEtiol"
l_adj="raw"#"combat"#
l_out="_cov_Int_time_APOE"#"noOutlier"

l_cohort=c(1,1,1,1,1,1,1,1,1)#,1)#c(0,0,1,0,0)# 1: BF2, 2: BF1, 3: AIBL, 4: TRIAD, 5: PREVENT
name_cohort=c("BF2","BF1","AIBL","TRIAD","PREVENT","WRAP","Ams","WU","MCSA")#,"A4")

for (i_cohort2 in (length(name_cohort)+1)) {##(length(name_cohort)+1)#
  
  df2=dataAll[complete.cases(dataAll[,c("APOE_e4","sex","age","edu","z_MTL","z_NeoT","z_ptau217")]),]
  s_cohort=NA
  l_cohort=rep(0,length(name_cohort))
  
  if (i_cohort2<(length(name_cohort)+1)) {
    s_cohort=paste(name_cohort[i_cohort2],"_",sep = "")
    df2=df2 %>% filter(cohort==name_cohort[i_cohort2])
    l_cohort[i_cohort2]=1
  }else{
    s_cohort="allCohorts"
    l_cohort=rep(1,length(name_cohort))
  }
  
  
  if (l_adj=="combat") {
    v_out=c(19,6)#c(238,235)#251:256,
    v_PET=c(18,17)#c(244:247)
    v_pl=16
    df2$plasma_pos=as.factor(ifelse(df2$z_ptau217_comb>1.96,1,0))#median(dataBl$z_ptau217_comb),1,0))#
    df2$MTL_pos=as.factor(ifelse(df2$z_MTL_comb>1.96,1,0))#median(dataBl$z_MTL_comb),1,0))#
    df2$NeoT_pos=as.factor(ifelse(df2$z_NeoT_comb>1.96,1,0))#median(dataBl$z_NeoT_comb),1,0))#
    
  }else{
    v_out=c(7,6)#c(238,235)#251:256,
    v_PET=c(9,8)#c(244:247)
    v_pl=10
    df2$plasma_pos=as.factor(ifelse(df2$z_ptau217>quantile(dataBl$z_ptau217,3/4),1,0))#1.96,1,0))
    df2$MTL_pos=as.factor(ifelse(df2$z_MTL>quantile(dataBl$z_MTL,3/4),1,0))#1.96,1,0))
    df2$NeoT_pos=as.factor(ifelse(df2$z_NeoT>quantile(dataBl$z_NeoT,3/4),1,0))#1.96,1,0))
    # v_PET_pos=c(20,19)
    c_PET_pos=c("NeoT_pos","MTL_pos")#names(df2)[v_PET_pos]
  }
  c_pl=names(df2)[v_pl]
  c_out=names(df2)[v_out]
  n_out=length(v_out)
  c_PET=names(df2)[v_PET]
  n_PET=length(v_PET)
  
  
  for (i_out in 1) {#:n_out
    
    results_B=data.frame(matrix(nrow = (2*n_PET)+3,ncol = 12))
    colnames(results_B)=c("Biomarker","name","N","B_plasma","p_plasma","B_PET","p_PET","R2","AICc","p_modelPET","p_modelPlasma","p_modelBasic")
  
    results_plot=data.frame(matrix(nrow = (2*n_PET)+3,ncol = 12))
    colnames(results_plot)=c("Biomarker","name","N","B_plasma","B_plasma_l","B_plasma_h","p_plasma","B_PET","B_PET_l","B_PET_h","p_PET","R2")
    
    df=df2[complete.cases(df2[,c(c_out[i_out])]),]
    
    
    if (l_ab=="Abpos") {
      df=df %>% filter(ab_pos==1)
    }
  
    # Select repetition data
    n_occur <- data.frame(table(df$sid))
    list_bl=n_occur[n_occur$Freq >1,1]
    
    df=filter(df, sid %in% list_bl)
    
    if (sum(!duplicated(df$sid))>10) {
      
    
    if (sum(l_cohort)>1) {
      fml_plasma=as.formula(paste(c_out[i_out],"~scale(",c_age,")*scale(timeDiff)+",c_sex,"*scale(timeDiff)+scale(",c_edu,")*scale(timeDiff)+",c_apoe,"*scale(timeDiff)+",c_cohort,"*scale(timeDiff)+",c_pl,"*scale(timeDiff)+(1+scale(timeDiff)|sid)"))
      fml_0=as.formula(paste(c_out[i_out],"~scale(",c_age,")*scale(timeDiff)+",c_sex,"*scale(timeDiff)+scale(",c_edu,")*scale(timeDiff)+",c_apoe,"*scale(timeDiff)+(1+scale(timeDiff)|sid)"))
      fml_00=as.formula(paste(c_out[i_out],"~scale(",c_age,")*scale(timeDiff)+",c_sex,"*scale(timeDiff)+scale(",c_edu,")*scale(timeDiff)+(1+scale(timeDiff)|sid)"))
    }else{
      fml_plasma=as.formula(paste(c_out[i_out],"~scale(",c_age,")*scale(timeDiff)+",c_sex,"*scale(timeDiff)+scale(",c_edu,")*scale(timeDiff)+",c_apoe,"*scale(timeDiff)+",c_pl,"*scale(timeDiff)+(1+scale(timeDiff)|sid)"))
      fml_0=as.formula(paste(c_out[i_out],"~scale(",c_age,")*scale(timeDiff)+",c_sex,"*scale(timeDiff)+scale(",c_edu,")*scale(timeDiff)+",c_apoe,"*scale(timeDiff)+(1+scale(timeDiff)|sid)"))
      fml_00=as.formula(paste(c_out[i_out],"~scale(",c_age,")*scale(timeDiff)+",c_sex,"*scale(timeDiff)+scale(",c_edu,")*scale(timeDiff)+(1+scale(timeDiff)|sid)"))
    }
    
    mdl_0=lmer(fml_0,df,REML = F)
    mdl_sum_0=summary(mdl_0)
    
    results_B$N=sum(!duplicated(df$sid))
    results_B$Biomarker[(2*n_PET)+2]="Basic_APOE"
    results_B$name[(2*n_PET)+2]="Basic_A"
    results_B$R2[(2*n_PET)+2]=paste(format(round(MuMIn::r.squaredGLMM(mdl_0)[1],2),nsmall=2))
    results_B$AICc[(2*n_PET)+2]=paste(format(round(MuMIn::AICc(mdl_0),1),nsmall=1))
    
    results_plot$N[(2*n_PET)+2]=sum(!duplicated(df$sid))
    results_plot$Biomarker[(2*n_PET)+2]="Basic_APOE"
    results_plot$name[(2*n_PET)+2]="Basic_A"
    results_plot$R2[(2*n_PET)+2]=MuMIn::r.squaredGLMM(mdl_0)[1]
    
    mdl_00=lmer(fml_00,df,REML = F)
    mdl_sum_00=summary(mdl_00)
    
    results_B$N=sum(!duplicated(df$sid))
    results_B$Biomarker[(2*n_PET)+3]="Basic_noAPOE"
    results_B$name[(2*n_PET)+3]="Basic_NA"
    results_B$R2[(2*n_PET)+3]=paste(format(round(MuMIn::r.squaredGLMM(mdl_00)[1],2),nsmall=2))
    results_B$AICc[(2*n_PET)+3]=paste(format(round(MuMIn::AICc(mdl_00),1),nsmall=1))
    
    results_plot$N[(2*n_PET)+3]=sum(!duplicated(df$sid))
    results_plot$Biomarker[(2*n_PET)+3]="Basic_noAPOE"
    results_plot$name[(2*n_PET)+3]="Basic_NA"
    results_plot$R2[(2*n_PET)+3]=MuMIn::r.squaredGLMM(mdl_00)[1]
    
    
    mdl_plasma=lmer(fml_plasma,df,REML = F)
    mdl_sum_plasma=summary(mdl_plasma)
    c_plasma=confint(mdl_plasma)
    # boot_test=boot(data=df,f_mdl,n_rep,fml=fml_plasma,v_row=c(nrow(mdl_sum_plasma$coefficients)))
    # boot_plasma=boot.ci(boot_test,type="norm",index=1)
    
    d_p0=anova(mdl_0,mdl_plasma)
    
    results_B$N=sum(!duplicated(df$sid))
    results_B$Biomarker[1]="Plasma"
    results_B$name[1]="Plasma p-tau217"
    results_B$B_plasma[1]=paste(format(round(mdl_sum_plasma$coefficients[nrow(mdl_sum_plasma$coefficients),1],2),nsmall=2)," [",format(round(c_plasma[nrow(mdl_sum_plasma$coefficients)+4,1],2),nsmall=2),", ",format(round(c_plasma[nrow(mdl_sum_plasma$coefficients)+4,2],2),nsmall=2),"]",sep = "")
    results_B$p_plasma[1]=paste(format(round(mdl_sum_plasma$coefficients[nrow(mdl_sum_plasma$coefficients),5],3),nsmall=3))
    results_B$R2[1]=paste(format(round(MuMIn::r.squaredGLMM(mdl_plasma)[1],2),nsmall=2))
    results_B$AICc[1]=paste(format(round(MuMIn::AICc(mdl_plasma),1),nsmall=1))
    results_B$p_modelBasic[1]=paste(format(round(d_p0$`Pr(>Chisq)`[2],3),nsmall=3))
    
    results_plot$N[1]=sum(!duplicated(df$sid))
    results_plot$Biomarker[1]="Plasma"
    results_plot$name[1]="Plasma p-tau217"
    results_plot$B_plasma[1]=mdl_sum_plasma$coefficients[nrow(mdl_sum_plasma$coefficients),1]
    results_plot$B_plasma_l[1]=c_plasma[nrow(mdl_sum_plasma$coefficients)+4,1]#boot_plasma$normal[2]
    results_plot$B_plasma_h[1]=c_plasma[nrow(mdl_sum_plasma$coefficients)+4,2]#boot_plasma$normal[3]
    results_plot$p_plasma[1]=mdl_sum_plasma$coefficients[nrow(mdl_sum_plasma$coefficients),5]
    results_plot$R2[1]=MuMIn::r.squaredGLMM(mdl_plasma)[1]
    
    if (i_cohort2==(length(name_cohort)+1)) {

    c_sign=""
    if(mdl_sum_plasma$coefficients[nrow(mdl_sum_plasma$coefficients),5]<0.05){c_sign="*"}
    if(mdl_sum_plasma$coefficients[nrow(mdl_sum_plasma$coefficients),5]<0.01){c_sign="**"}
    if(mdl_sum_plasma$coefficients[nrow(mdl_sum_plasma$coefficients),5]<0.001){c_sign="***"}
    
    if (sum(l_cohort)>1) {
      fml_line_plasma=as.formula(paste(c_out[i_out],"~",c_age,"*timeDiff+",c_sex,"*timeDiff+",c_edu,"*timeDiff+",c_apoe,"*timeDiff+",c_cohort,"*timeDiff+plasma_pos*timeDiff+(1+timeDiff|sid)"))
    }else{
      fml_line_plasma=as.formula(paste(c_out[i_out],"~",c_age,"*timeDiff+",c_sex,"*timeDiff+",c_edu,"*timeDiff+",c_apoe,"*timeDiff+plasma_pos*timeDiff+(1+timeDiff|sid)"))
    }
    mdl_line_plasma=lmer(fml_line_plasma,data = df,REML = F)
    
    d_line=ggeffect(mdl_line_plasma,terms = c("timeDiff","plasma_pos")) #It takes into account the proportion of factors as compared to ggpredict
    
    fig1=ggplot(data = df, aes(y =get(c_out[i_out]), x = timeDiff)) +
      geom_line(data = df, aes(x = timeDiff, y = get(c_out[i_out]),group = sid,color=plasma_pos), alpha=.15,size=0.5) +
      geom_point(data = df, aes(y = get(c_out[i_out]), x = timeDiff,color=plasma_pos, alpha=plasma_pos), size=0.5)+
      labs(x = "Time (years)", y =paste(c_out_name[i_out])) + theme_classic()+ #xlim(c(v_lim[i_count],v_lim[i_count+1]))+#geom_smooth(method = "lm",na.rm = T,color="black")#+
      scale_alpha_manual(values = c(0.3,0.7))+ylim(c(-4,2))+xlim(c(-1,8))+
      geom_line(data = d_line,aes(x=x,y=predicted,color=group),size=1)+
      geom_ribbon(data = d_line,aes(x= x, y= predicted, ymin = conf.low, ymax = conf.high,color=group), alpha = .1, linetype = 0)+
      scale_color_manual(values=c("lightgrey","#e64b35ff"))+#"#3182bd"))+
      geom_text(x=-Inf,y=-Inf,label=paste("R2=",round(MuMIn::r.squaredGLMM(mdl_plasma)[1],2),"   ",sep = ""),size=3,color="black",hjust   = -0.1,
                vjust   = -2)+
      geom_text(x=-Inf,y=-Inf,label=paste("B=",round(mdl_sum_plasma$coefficients[nrow(mdl_sum_plasma$coefficients),1],2),c_sign,sep = ""),size=3,color="black",hjust   = -0.1,#size=10
                vjust   = -0.6)+
      theme(legend.position = "none",axis.text=element_text(size=22),
                                      axis.title=element_text(size=26,face="bold"))+
      theme(legend.position = "none",axis.text=element_text(size=8),
            axis.title=element_text(size=10,face="bold"),legend.text = element_text(size=6),
            axis.line.x = element_line(size = 0.3),axis.line.y = element_line(size = 0.3),
            axis.ticks.x=element_line(size = 0.3),axis.ticks.y = element_line(size = 0.3),
            legend.key.size = unit(0.2, 'cm'),legend.title = element_blank())
      
    }

    # ggexport(fig1,filename=paste(dir_fig,"/LME_mergedCohorts_",s_cohort,c_out[i_out],"_DX_",l_DX,"_Ab_",l_ab,"_adj_",l_adj,"_plasma",l_out,"_SensTime_",v_time_min,".pdf",sep = ""))
    
    print("Plasma")
    
    for (i_PET in 1:n_PET) {
      
      if (sum(l_cohort)>1) {
        fml_PET=as.formula(paste(c_out[i_out],"~scale(",c_age,")*scale(timeDiff)+",c_sex,"*scale(timeDiff)+scale(",c_edu,")*scale(timeDiff)+",c_apoe,"*scale(timeDiff)+",c_cohort,"*scale(timeDiff)+",c_PET[i_PET],"*scale(timeDiff)+(1+scale(timeDiff)|sid)"))
      }else{
        fml_PET=as.formula(paste(c_out[i_out],"~scale(",c_age,")*scale(timeDiff)+",c_sex,"*scale(timeDiff)+scale(",c_edu,")*scale(timeDiff)+",c_apoe,"*scale(timeDiff)+",c_PET[i_PET],"*scale(timeDiff)+(1+scale(timeDiff)|sid)"))
      }
      mdl_PET=lmer(fml_PET,df,REML = F)
      mdl_sum_PET=summary(mdl_PET)
      c=confint(mdl_PET)
      # boot_test=boot(data=df,f_mdl,n_rep,fml=fml_PET,v_row=c(nrow(mdl_sum_PET$coefficients)))
      # boot_PET=boot.ci(boot_test,type="norm",index=1)
      
      d_PET0=anova(mdl_0,mdl_PET)
      
      results_B$Biomarker[1+i_PET]=c_PET[i_PET]
      results_B$name[1+i_PET]=c_name[i_PET]
      results_B$B_PET[1+i_PET]=paste(format(round(mdl_sum_PET$coefficients[nrow(mdl_sum_PET$coefficients),1],2),nsmall=2)," [",format(round(c[nrow(mdl_sum_PET$coefficients)+4,1],2),nsmall=2),", ",format(round(c[nrow(mdl_sum_PET$coefficients)+4,2],2),nsmall=2),"]",sep = "")
      results_B$p_PET[1+i_PET]=paste(format(round(mdl_sum_PET$coefficients[nrow(mdl_sum_PET$coefficients),5],3),nsmall=3))
      results_B$R2[1+i_PET]=paste(format(round(MuMIn::r.squaredGLMM(mdl_PET)[1],2),nsmall=2))
      results_B$AICc[1+i_PET]=paste(format(round(MuMIn::AICc(mdl_PET),1),nsmall=1))
      results_B$p_modelBasic[1+i_PET]=paste(format(round(d_PET0$`Pr(>Chisq)`[2],3),nsmall=3))
      
      results_plot$Biomarker[1+i_PET]=c_PET[i_PET]
      results_plot$name[1+i_PET]=c_name[i_PET]
      results_plot$B_PET[1+i_PET]=mdl_sum_PET$coefficients[nrow(mdl_sum_PET$coefficients),1]
      results_plot$B_PET_l[1+i_PET]=c[nrow(mdl_sum_PET$coefficients)+4,1]#boot_PET$normal[2]
      results_plot$B_PET_h[1+i_PET]=c[nrow(mdl_sum_PET$coefficients)+4,1]#boot_PET$normal[3]
      results_plot$p_PET[1+i_PET]=mdl_sum_PET$coefficients[nrow(mdl_sum_PET$coefficients),5]
      results_plot$R2[1+i_PET]=MuMIn::r.squaredGLMM(mdl_PET)[1]
      
      if (i_cohort2==(length(name_cohort)+1)) {
        
      c_sign=""
      if(mdl_sum_PET$coefficients[nrow(mdl_sum_PET$coefficients),5]<0.05){c_sign="*"}
      if(mdl_sum_PET$coefficients[nrow(mdl_sum_PET$coefficients),5]<0.01){c_sign="**"}
      if(mdl_sum_PET$coefficients[nrow(mdl_sum_PET$coefficients),5]<0.001){c_sign="***"}
      
      
      if (sum(l_cohort)>1) {
        fml_line_PET=as.formula(paste(c_out[i_out],"~",c_age,"*timeDiff+",c_sex,"*timeDiff+",c_edu,"*timeDiff+",c_apoe,"*timeDiff+",c_cohort,"*timeDiff+",c_PET_pos[i_PET],"*timeDiff+(1+timeDiff|sid)"))
      }else{
        fml_line_PET=as.formula(paste(c_out[i_out],"~",c_age,"*timeDiff+",c_sex,"*timeDiff+",c_edu,"*timeDiff+",c_apoe,"*timeDiff+",c_PET_pos[i_PET],"*timeDiff+(1+timeDiff|sid)"))
      }
      mdl_line_PET=lmer(fml_line_PET,data = df,REML = F)
      
      d_line=ggeffect(mdl_line_PET,terms = c("timeDiff",c_PET_pos[i_PET])) #It takes into account the proportion of factors as compared to ggpredict
      
      fig2=ggplot(data = df, aes(y =get(c_out[i_out]), x = timeDiff)) +
        geom_line(data = df, aes(x = timeDiff, y = get(c_out[i_out]),group = sid,color=get(c_PET_pos[i_PET])), alpha=.15,size=0.5) +
        geom_point(data = df, aes(y = get(c_out[i_out]), x = timeDiff,color=get(c_PET_pos[i_PET]),alpha=get(c_PET_pos[i_PET])), size=0.5)+
        labs(x = "Time (years)", y =paste(c_out_name[i_out])) + theme_classic()+ #xlim(c(v_lim[i_count],v_lim[i_count+1]))+#geom_smooth(method = "lm",na.rm = T,color="black")#+
        scale_alpha_manual(values = c(0.3,0.7))+ylim(c(-4,2))+xlim(c(-1,8))+
        geom_line(data = d_line,aes(x=x,y=predicted,color=group),size=1)+
        geom_ribbon(data = d_line,aes(x= x, y= predicted, ymin = conf.low, ymax = conf.high,color=group), alpha = .1, linetype = 0)+
        scale_color_manual(values=c("lightgrey",c_colors_PET[i_PET]))+
        geom_text(x=-Inf,y=-Inf,label=paste("R2=",round(MuMIn::r.squaredGLMM(mdl_PET)[1],2),"   ",sep = ""),size=3,color="black",hjust   = -0.1,
                  vjust   = -2)+
        geom_text(x=-Inf,y=-Inf,label=paste("B=",round(mdl_sum_PET$coefficients[nrow(mdl_sum_PET$coefficients),1],2),c_sign,sep = ""),size=3,color="black",hjust   = -0.1,
                  vjust   = -0.6)+
        theme(legend.position = "none",axis.text=element_text(size=22),
              axis.title=element_text(size=26,face="bold"))+
        theme(legend.position = "none",axis.text=element_text(size=8),
              axis.title=element_text(size=10,face="bold"),legend.text = element_text(size=6),
              axis.line.x = element_line(size = 0.3),axis.line.y = element_line(size = 0.3),
              axis.ticks.x=element_line(size = 0.3),axis.ticks.y = element_line(size = 0.3),
              legend.key.size = unit(0.2, 'cm'),legend.title = element_blank())
        
      }
      # ggexport(fig2,filename=paste(dir_fig,"/LME_mergedCohorts_",s_cohort,c_out[i_out],"_DX_",l_DX,"_Ab_",l_ab,"_adj_",l_adj,"_PET_",c_PET[i_PET],l_out,"_SensTime_",v_time_min,".pdf",sep = ""))
      
      if (sum(l_cohort)>1) {
        fml_mix=as.formula(paste(c_out[i_out],"~scale(",c_age,")*scale(timeDiff)+",c_sex,"*scale(timeDiff)+
                               scale(",c_edu,")*scale(timeDiff)+",c_apoe,"*scale(timeDiff)+",c_cohort,"*scale(timeDiff)+
                               ",c_pl,"*scale(timeDiff)+",
                               c_PET[i_PET],"*scale(timeDiff)+(1+scale(timeDiff)|sid)"))
        
      }else{
        fml_mix=as.formula(paste(c_out[i_out],"~scale(",c_age,")*scale(timeDiff)+",c_sex,"*scale(timeDiff)+
                               scale(",c_edu,")*scale(timeDiff)+",c_apoe,"*scale(timeDiff)+
                               ",c_pl,"*scale(timeDiff)+",
                               c_PET[i_PET],"*scale(timeDiff)+(1+scale(timeDiff)|sid)"))
      }
      
      mdl_mix=lmer(fml_mix,df,REML = F)
      mdl_sum_mix=summary(mdl_mix)
      set.seed(1234)
      c_mix=confint(mdl_mix,devtol=Inf)
      
      # boot_test=boot(data=df,f_mdl,n_rep,fml=fml_mix,v_row=c(nrow(mdl_sum_mix$coefficients)-1,nrow(mdl_sum_mix$coefficients)))
      # boot_plasma=boot.ci(boot_test,type="norm",index=1)
      # boot_PET=boot.ci(boot_test,type="norm",index=2)
      
      d_PET=anova(mdl_PET,mdl_mix)
      d_plasma=anova(mdl_plasma,mdl_mix)
      d_mix0=anova(mdl_0,mdl_mix)
      
      results_B$Biomarker[1+n_PET+i_PET]=paste("Plasma+",c_PET[i_PET],sep = "")
      results_B$name[1+n_PET+i_PET]=paste("Plasma+",c_name[i_PET],sep = "")
      results_B$B_plasma[1+n_PET+i_PET]=paste(format(round(mdl_sum_mix$coefficients[nrow(mdl_sum_mix$coefficients)-1,1],2),nsmall=2)," [",format(round(c_mix[nrow(mdl_sum_mix$coefficients)-1+4,1],2),nsmall=2),", ",format(round(c_mix[nrow(mdl_sum_mix$coefficients)-1+4,2],2),nsmall=2),"]",sep = "")#c[nrow(c)-1,1],2),nsmall=2),", ",format(round(c[nrow(c)-1,2],2),nsmall=2),"]",sep = "")
      results_B$p_plasma[1+n_PET+i_PET]=paste(format(round(mdl_sum_mix$coefficients[nrow(mdl_sum_mix$coefficients)-1,5],3),nsmall=3))
      results_B$B_PET[1+n_PET+i_PET]=paste(format(round(mdl_sum_mix$coefficients[nrow(mdl_sum_mix$coefficients),1],2),nsmall=2)," [",format(round(c_mix[nrow(mdl_sum_mix$coefficients)+4,1],2),nsmall=2),", ",format(round(c_mix[nrow(mdl_sum_mix$coefficients)+4,2],2),nsmall=2),"]",sep = "")
      results_B$p_PET[1+n_PET+i_PET]=paste(format(round(mdl_sum_mix$coefficients[nrow(mdl_sum_mix$coefficients),5],3),nsmall=3))
      results_B$R2[1+n_PET+i_PET]=paste(format(round(MuMIn::r.squaredGLMM(mdl_mix)[1],2),nsmall=2))
      results_B$AICc[1+n_PET+i_PET]=paste(format(round(MuMIn::AICc(mdl_mix),1),nsmall=1))
      results_B$p_modelPlasma[1+n_PET+i_PET]=paste(format(round(d_plasma$`Pr(>Chisq)`[2],3),nsmall=3))
      results_B$p_modelPET[1+n_PET+i_PET]=paste(format(round(d_PET$`Pr(>Chisq)`[2],3),nsmall=3))
      results_B$p_modelBasic[1+n_PET+i_PET]=paste(format(round(d_mix0$`Pr(>Chisq)`[2],3),nsmall=3))
      
      results_plot$Biomarker[1+n_PET+i_PET]=paste("Plasma+",c_PET[i_PET],sep = "")
      results_plot$name[1+n_PET+i_PET]=paste("Plasma+",c_name[i_PET],sep = "")
      results_plot$B_plasma[1+n_PET+i_PET]=mdl_sum_mix$coefficients[nrow(mdl_sum_mix$coefficients)-1,1]
      results_plot$B_plasma_l[1+n_PET+i_PET]=c_mix[nrow(mdl_sum_mix$coefficients)-1+4,1]#boot_plasma$normal[2]
      results_plot$B_plasma_h[1+n_PET+i_PET]=c_mix[nrow(mdl_sum_mix$coefficients)-1+4,2]#boot_plasma$normal[3]
      results_plot$p_plasma[1+n_PET+i_PET]=mdl_sum_mix$coefficients[nrow(mdl_sum_mix$coefficients)-1,5]
      results_plot$B_PET[1+n_PET+i_PET]=mdl_sum_mix$coefficients[nrow(mdl_sum_mix$coefficients),1]
      results_plot$B_PET_l[1+n_PET+i_PET]=c_mix[nrow(mdl_sum_mix$coefficients)+4,1]#boot_PET$normal[2]
      results_plot$B_PET_h[1+n_PET+i_PET]=c_mix[nrow(mdl_sum_mix$coefficients)+4,2]#boot_PET$normal[3]
      results_plot$p_PET[1+n_PET+i_PET]=mdl_sum_mix$coefficients[nrow(mdl_sum_mix$coefficients),5]
      results_plot$R2[1+n_PET+i_PET]=MuMIn::r.squaredGLMM(mdl_mix)[1]
      #
      c_sign1=""
      if(mdl_sum_mix$coefficients[nrow(mdl_sum_mix$coefficients)-1,5]<0.05){c_sign1="*"}
      if(mdl_sum_mix$coefficients[nrow(mdl_sum_mix$coefficients)-1,5]<0.01){c_sign1="**"}
      if(mdl_sum_mix$coefficients[nrow(mdl_sum_mix$coefficients)-1,5]<0.001){c_sign1="***"}

      c_sign2=""
      if(mdl_sum_mix$coefficients[nrow(mdl_sum_mix$coefficients),5]<0.05){c_sign2="*"}
      if(mdl_sum_mix$coefficients[nrow(mdl_sum_mix$coefficients),5]<0.01){c_sign2="**"}
      if(mdl_sum_mix$coefficients[nrow(mdl_sum_mix$coefficients),5]<0.001){c_sign2="***"}


      df$mix_group=interaction(df$plasma_pos,df[,c_PET_pos[i_PET]])
      levels(df$mix_group)=c("p-PET-","p+PET-","p-PET+","p+PET+")

      if (sum(l_cohort)>1) {
        fml_line_mix=as.formula(paste(c_out[i_out],"~",c_age,"+",c_sex,"+",c_edu,"+",c_apoe,"+",c_cohort,"*timeDiff+plasma_pos*timeDiff+",c_PET_pos[i_PET],"*timeDiff+(1+timeDiff|sid)"))
      }else{
        fml_line_mix=as.formula(paste(c_out[i_out],"~",c_age,"+",c_sex,"+",c_edu,"+",c_apoe,"+plasma_pos*timeDiff+",c_PET_pos[i_PET],"*timeDiff+(1+timeDiff|sid)"))
      }
      mdl_line_mix=lmer(fml_line_mix,data = df,REML = F)

      d_line_plasma=ggeffect(mdl_line_mix,terms = c("timeDiff","plasma_pos")) #It takes into account the proportion of factors as compared to ggpredict
      d_line_PET=ggeffect(mdl_line_mix,terms = c("timeDiff",(c_PET_pos[i_PET]))) #It takes into account the proportion of factors as compared to ggpredict

      d_line_plasma$group=ifelse(d_line_plasma$group==0,"p-","p+")
      d_line_PET$group=ifelse(d_line_PET$group==0,"T-","T+")

      d_line_neg=d_line_plasma[d_line_plasma$group=="p-",]
      d_line_neg=merge(d_line_neg,d_line_PET[d_line_PET$group=="T-",],by="x")
      d_line_neg$predicted=rowMeans(d_line_neg[,c("predicted.x","predicted.y")])
      d_line_neg$conf.low=rowMeans(d_line_neg[,c("conf.low.x","conf.low.y")])
      d_line_neg$conf.high=rowMeans(d_line_neg[,c("conf.high.x","conf.high.y")])

      c_colors=c("p-PET-"="lightgrey",
                 "p-"="lightgrey",
                 "p+PET-"="#3182bd",
                 "p+"="#3182bd",
                 "p-PET+"=c_colors_PET[i_PET],
                 "T+"=c_colors_PET[i_PET],
                 "p+PET+"=c_colors_mix[i_PET])

      fig3=ggplot(data = df, aes(y =get(c_out[i_out]), x = timeDiff)) +
        geom_line(data = df, aes(x = timeDiff, y = get(c_out[i_out]),group = sid,color=mix_group), alpha=.15,size=0.5) +
        geom_point(data = df, aes(y = get(c_out[i_out]), x = timeDiff,color=mix_group,alpha=mix_group), size=0.5)+
        labs(x = "Time (years)", y =paste(c_out_name[i_out])) + theme_classic()+ #xlim(c(v_lim[i_count],v_lim[i_count+1]))+#geom_smooth(method = "lm",na.rm = T,color="black")#+
        scale_alpha_manual(values = c(0.3,0.7,0.7,0.7))+
        geom_line(data = d_line_plasma[d_line_plasma$group=="p+",],aes(x=x,y=predicted,color=group),size=1)+
        geom_ribbon(data = d_line_plasma[d_line_plasma$group=="p+",],aes(x= x, y= predicted, ymin = conf.low, ymax = conf.high,color=group), alpha = .1, linetype = 0)+
        geom_line(data = d_line_PET[d_line_PET$group=="T+",],aes(x=x,y=predicted,color=group),size=1)+
        geom_ribbon(data = d_line_PET[d_line_PET$group=="T+",],aes(x= x, y= predicted, ymin = conf.low, ymax = conf.high,color=group), alpha = .1, linetype = 0)+
        geom_line(data = d_line_neg,aes(x=x,y=predicted,color=group.x),size=1)+
        geom_ribbon(data = d_line_neg,aes(x= x, y= predicted, ymin = conf.low, ymax = conf.high,color=group.x), alpha = .1, linetype = 0)+
        scale_color_manual(values=c_colors)+#("lightgrey","#3182bd","#31a354","#e6550d"))+
        geom_text(x=-Inf,y=-Inf,label=paste("R2=",round(MuMIn::r.squaredGLMM(mdl_mix)[1],2),"   ",sep = ""),size=3,color="black",hjust   = -0.1,
                  vjust   = -3.6)+
        geom_text(x=-Inf,y=-Inf,label=paste("Bplasma=",round(mdl_sum_mix$coefficients[nrow(mdl_sum_mix$coefficients)-1,1],2),c_sign1," ",sep = ""),size=3,color="black",hjust   = -0.1,
                  vjust   = -2)+
        geom_text(x=-Inf,y=-Inf,label=paste("Bpet=",round(mdl_sum_mix$coefficients[nrow(mdl_sum_mix$coefficients),1],2),c_sign2,sep = ""),size=3,color="black",hjust   = -0.1,
                  vjust   = -0.4)+
        theme(legend.position = "none",axis.text=element_text(size=22),
              axis.title=element_text(size=26,face="bold"))+
        theme(legend.position = "none",axis.text=element_text(size=8),
              axis.title=element_text(size=10,face="bold"),legend.text = element_text(size=6),
              axis.line.x = element_line(size = 0.3),axis.line.y = element_line(size = 0.3),
              axis.ticks.x=element_line(size = 0.3),axis.ticks.y = element_line(size = 0.3),
              legend.key.size = unit(0.2, 'cm'),legend.title = element_blank())


      # ggexport(fig3,filename=paste(dir_fig,"/LME_mergedCohorts_",s_cohort,c_out[i_out],"_DX_",l_DX,"_Ab_",l_ab,"_plasmaPET_",c_PET[i_PET],l_out,"_SensTime_",v_time_min,".pdf",sep = ""))
      
      if (i_cohort2==(length(name_cohort)+1)) {
        if (i_PET==1) {
          fig_NeoT=fig2
          fig_NeoT_plasma=fig3
        }
      }
      
      print(c_PET[i_PET])
      
    }
    
    if (i_cohort2==(length(name_cohort)+1)) {
      
    # pf2=ggarrange(fig1+rremove("ylab"),  fig2 +rremove("ylab"), fig_NeoT+rremove("ylab"), #fig3, fig_NeoT_plasma,
    #               # labels = c("A", "B", "C","D","E"),
    #               ncol = 6, nrow = 6)
    
    # ggexport(pf2,filename=paste(dir_fig,"/LME_mergedCohorts_",s_cohort,c_out[i_out],"_DX_",l_DX,"_Ab_",l_ab,"_adj_",l_adj,"_allPlots_",l_out,"_SensTime_",v_time_min,".pdf",sep = ""))
    
    pf5=ggarrange(fig1+rremove("ylab"),  fig2 +rremove("ylab"), fig_NeoT+rremove("ylab"), #fig3, fig_NeoT_plasma,
                  # labels = c("A", "B", "C","D","E"),
                  ncol = 3, nrow = 3)

    ggexport(pf5,filename=paste(dir_fig,"/LME_mergedCohorts_",s_cohort,c_out[i_out],"_DX_",l_DX,"_Ab_",l_ab,"_adj_",l_adj,"_allPlots_",l_out,"_SensTime_",v_time_min,".pdf",sep = ""))

   }
    
    writexl::write_xlsx(results_B,path = paste(dir_res,"/Results_LME_mergedCohorts_",s_cohort,c_out[i_out],"_DX_",l_DX,"_Ab_",l_ab,"_adj_",l_adj,l_out,"_SensTime_",v_time_min,".xlsx",sep = ""))
    writexl::write_xlsx(results_plot,path = paste(dir_res,"/Results_plot_LME_mergedCohorts_",s_cohort,c_out[i_out],"_DX_",l_DX,"_Ab_",l_ab,"_adj_",l_adj,l_out,"_SensTime_",v_time_min,".xlsx",sep = ""))
    
    }
  }
}



