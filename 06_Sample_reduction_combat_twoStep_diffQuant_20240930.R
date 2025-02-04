dir_res="~/Desktop/Lund/20_TauPET_plasma/Results_Merged/20240917_Rev/Two_step_sample_red"
dir_fig="~/Desktop/Lund/20_TauPET_plasma/Figures_Merged/20240917_Rev/Two_step_sample_red"
load("~/Desktop/Lund/20_TauPET_plasma/Data_all/Data_all_cogn_comb_20240917.RData")

library(longpower)
# 
# dataBl=dataAll %>% filter(!duplicated(sid))
# median(dataBl$z_ptau217_comb)

# dataAll$plasma_pos=as.factor(ifelse(dataAll$z_ptau217_comb>median(dataBl$z_ptau217_comb),1,0))#1.96,1,0))
# dataAll$MTL_pos=as.factor(ifelse(dataAll$z_MTL_comb>median(dataBl$z_MTL_comb),1,0))#1.96,1,0))
# dataAll$NeoT_pos=as.factor(ifelse(dataAll$z_NeoT_comb>median(dataBl$z_NeoT_comb),1,0))#1.96,1,0))
# dataAll$cohort=as.factor(dataAll$cohort)


###Analyses
# Function to use bootstrap for calculating CI and p-value
f_N <- function(data, indices,c_out0) {
  d <- data[indices,] # allows boot to select sample
  
  fml=as.formula(paste(c_out[i_out],"~timeDiff+(1+timeDiff|sid)"))
  
  mdl=lmer(fml,d)
  
  N_or=lmmpower(mdl, pct.change = v_change, t = seq(0,4,1), power = v_power)
  
  mdl_plasma=lmer(fml,d[which(d$plasma_pos==1),])
  N_plasma=lmmpower(mdl_plasma, pct.change = v_change, t = seq(0,4,1), power = v_power)
  
  # mdl_MTL=lmer(fml,d[which(d$MTL_pos==1),],REML = F)
  # N_MTL=lmmpower(mdl_MTL, pct.change = v_change, t = seq(0,4,1), power = v_power)
  # 
  # mdl_NeoT=lmer(fml,d[which(d$NeoT_pos==1),],REML = F)
  # N_NeoT=lmmpower(mdl_NeoT, pct.change = v_change, t = seq(0,4,1), power = v_power)
  
  mdl_MTL_plasma=lmer(fml,d[which(d$MTL_pos==1 & d$plasma_pos==1),],REML = F)
  N_MTL_plasma=lmmpower(mdl_MTL_plasma, pct.change = v_change, t = seq(0,4,1), power = v_power)
  
  mdl_NeoT_plasma=lmer(fml,d[which(d$NeoT_pos==1 & d$plasma_pos==1),],REML = F)
  N_NeoT_plasma=lmmpower(mdl_NeoT_plasma, pct.change = v_change, t = seq(0,4,1), power = v_power)
  
  p_plasma0=(N_or$N-N_plasma$N)
  # p_MTL0=(N_or$N-N_MTL$N)
  # p_NeoT0=(N_or$N-N_NeoT$N)
  p_pMTL0=(N_or$N-N_MTL_plasma$N)
  p_pNeoT0=(N_or$N-N_NeoT_plasma$N)
  # p_MTL_plasma=-(N_MTL$N-N_plasma$N)
  # p_NeoT_plasma=-(N_NeoT$N-N_plasma$N)
  # p_MTL_NeoT=(N_NeoT$N-N_MTL$N)
  # p_MTL_MTLplasma=(N_MTL$N-N_MTL_plasma$N)
  p_plasma_MTLplasma=(N_plasma$N-N_MTL_plasma$N)
  p_plasma_NeoTplasma=(N_plasma$N-N_NeoT_plasma$N)
  # p_NeoT_NeoTplasma=(N_NeoT$N-N_NeoT_plasma$N)
  p_pMTL_pNeoT=-(N_MTL_plasma$N-N_NeoT_plasma$N)#ULL!
  
  
  return(c(N_or$N,N_plasma$N,
           # N_MTL$N,N_NeoT$N,
           N_MTL_plasma$N,N_NeoT_plasma$N,
           p_plasma0,
           # p_MTL0,p_NeoT0,
           p_pMTL0,p_pNeoT0,
           # p_MTL_plasma,p_NeoT_plasma,p_MTL_NeoT,
           p_plasma_MTLplasma,#p_MTL_MTLplasma,
           p_plasma_NeoTplasma,#p_NeoT_NeoTplasma,
           100*N_plasma$N/N_or$N,
           #100*N_MTL$N/N_or$N,100*N_NeoT$N/N_or$N,
           100*N_MTL_plasma$N/N_or$N,100*N_NeoT_plasma$N/N_or$N,p_pMTL_pNeoT,
           100*N_MTL_plasma$N/N_plasma$N,100*N_NeoT_plasma$N/N_plasma$N))
}
dataAll$cohort=as.factor(dataAll$cohort)
dataBl=dataAll %>% filter(!duplicated(sid))

c_colors=c("Basic_NA"="#d9d9d9",
           "Basic_A"="#737373",#"Original"="lightgrey",
           "Plasma"="#e64b35ff",
           # "MTL"="#a1d99b",
           # "NeoT"="#31a354",
           "p-MTL"="#00a087ff",
           "p-NeoT"="#3c5488ff")


# c_age="age"
# c_sex="sex"
# c_edu="edu"
# c_apoe="APOE_e4"
# c_ID="sid"
# c_cohort="cohort"
# c_tracer="tracer"

# v_out=c(19,6)#c(238,235)#251:256,
# c_out=names(dataAll)[v_out]
c_out_name=c("mPACC","MMSE")
n_out=length(c_out_name)

# v_PET=c(18,17)#c(244:247)
# c_PET=names(dataAll)[v_PET]
# n_PET=length(v_PET)
# c_name=c("NeoT","ERC/Amygd")
# v_PET_pos=c(22,21)
# c_PET_pos=names(dataAll)[v_PET_pos]

n_boot=1000
v_power=0.8#0.6
v_change=0.3
v_time_final=4
v_time_int=1#2#1

l_ab="all"#"Abpos"#
l_DX="allEtiol"
l_adj="raw"#"combat"
l_out="_noCovs"#"noOutlier"
c_method_plasma=c("Q2_Q4","Q3_Q4","Q4")#"zscore2"
c_method_PET=c("Q2_Q4","Q3_Q4","Q4")#"zscore2"

for (i_method in 1:length(c_method_plasma)) {#
  
  for (i_method2 in 1:length(c_method_PET)) {#
    
    if (c_method_plasma[i_method]=="Q4") {
      v_perc=0.75
    }
    if (c_method_plasma[i_method]=="Q3_Q4") {
      v_perc=0.5
    }
    if (c_method_plasma[i_method]=="Q2_Q4") {
      v_perc=0.25
    }
    
    if (c_method_PET[i_method2]=="Q4") {
      v_perc_PET=0.75
    }
    if (c_method_PET[i_method2]=="Q3_Q4") {
      v_perc_PET=0.5
    }
    if (c_method_PET[i_method2]=="Q2_Q4") {
      v_perc_PET=0.25
    }
  
  
  results_print=data.frame(matrix(nrow = n_out,ncol = 16))
  colnames(results_print)=c("Outcome",
                            "N_or","N_plasma",#"N_MTL","N_NeoT",
                            "N_plasmaMTL","N_plasmaNeoT",
                            "p_plasma_or",#"p_MTL_or","p_NeoT_or",
                            "p_pMTL_or","p_pNeoT_or",
                            #"p_MTL_plasma","p_NeoT_plasma","p_MTL_NeoT",
                            "p_plasma_pMTL",#"p_MTL_pMTL",
                            "p_plasma_pNeoT",#"p_NeoT_pNeoT",
                            "p_pMTL_pNeoT",
                            "perc_plasma",#"perc_MTL","perc_NeoT",
                            "perc_plasmaMTL","perc_plasmaNeoT",
                            "perc_plasmaMTL_plasma","perc_plasmaNeoT_plasma")
  
  
  ### Change n-boot for final calculations
  
  for (i_out in 1) {#:n_out
    
    results_plot=data.frame(matrix(nrow = 5,ncol = 4))
    colnames(results_plot)=c("Model",
                             "Perc","Perc_l","Perc_h")
    df=dataAll
    
    if (l_adj=="combat") {
      v_out=c(19,6)#c(238,235)#251:256,
      v_PET=c(18,17)#c(244:247)
      v_pl=16
      # df$plasma_pos=as.factor(ifelse(df$z_ptau217_comb>median(dataBl$z_ptau217_comb),1,0))#1.96,1,0))
      # df$MTL_pos=as.factor(ifelse(df$z_MTL_comb>median(dataBl$z_MTL_comb),1,0))#1.96,1,0))
      # df$NeoT_pos=as.factor(ifelse(df$z_NeoT_comb>median(dataBl$z_NeoT_comb),1,0))#1.96,1,0))
      
    }else{
      v_out=c(7,6)#c(238,235)#251:256,
      v_PET=c(9,8)#c(244:247)
      v_pl=10
      # df$plasma_pos=as.factor(ifelse(df$z_ptau217>median(dataBl$z_ptau217),1,0))#1.96,1,0))
      # df$MTL_pos=as.factor(ifelse(df$z_MTL>median(dataBl$z_MTL),1,0))#1.96,1,0))
      # df$NeoT_pos=as.factor(ifelse(df$z_NeoT>median(dataBl$z_NeoT),1,0))#1.96,1,0))
    }
    c_pl=names(df)[v_pl]
    c_out=names(df)[v_out]
    c_PET=names(df)[v_PET]
    n_PET=length(v_PET)
      
    df=df[complete.cases(df[,c_out[i_out]]),]
    
    if (l_ab=="Abpos") {
      df=df %>% filter(ab_pos==1)
    }
    
    # Select repetition data
    n_occur <- data.frame(table(df$sid))
    list_bl=n_occur[n_occur$Freq >1,1]
    
    df=filter(df, sid %in% list_bl)
    
    df$plasma_pos=as.factor(ifelse(df[,c_pl]>quantile(df[!duplicated(df$sid),c_pl],v_perc),1,0))#as.factor(ifelse(dataAll$z_ptau217_comb>1.96,1,0))
    
    df2=df %>% filter(!duplicated(sid))
    df2=df2 %>% filter(plasma_pos==1)
    
    #We select the quartiles of those participants already plasma positive to assess the additive reduction of performing a PET after being plasma+
    df$MTL_pos=as.factor(ifelse(df[,c_PET[2]]>quantile(df2[,c_PET[2]],v_perc_PET),1,0))
    df$NeoT_pos=as.factor(ifelse(df[,c_PET[1]]>quantile(df2[,c_PET[1]],v_perc_PET),1,0))
    
    # v_PET_pos=c(22,21)
    # c_PET_pos=names(df)[v_PET_pos]
    
    fml=as.formula(paste(c_out[i_out],"~timeDiff+(1+timeDiff|sid)"))
    
    mdl=lmer(fml,df)
    
    N_or=lmmpower(mdl, pct.change = v_change, t = seq(0,v_time_final,v_time_int), power = v_power)
    
    mdl_plasma=lmer(fml,df[which(df$plasma_pos==1),])
    N_plasma=lmmpower(mdl_plasma, pct.change = v_change, t = seq(0,v_time_final,v_time_int), power = v_power)
    
    # mdl_MTL=lmer(fml,df[which(df$MTL_pos==1),],REML = F)
    # N_MTL=lmmpower(mdl_MTL, pct.change = v_change, t = seq(0,4,1), power = v_power)
    # 
    # mdl_NeoT=lmer(fml,df[which(df$NeoT_pos==1),],REML = F)
    # N_NeoT=lmmpower(mdl_NeoT, pct.change = v_change, t = seq(0,4,1), power = v_power)
    
    mdl_MTL_plasma=lmer(fml,df[which(df$MTL_pos==1 & df$plasma_pos==1),],REML = F)
    N_MTL_plasma=lmmpower(mdl_MTL_plasma, pct.change = v_change, t = seq(0,v_time_final,v_time_int), power = v_power)
    
    mdl_NeoT_plasma=lmer(fml,df[which(df$NeoT_pos==1 & df$plasma_pos==1),],REML = F)
    N_NeoT_plasma=lmmpower(mdl_NeoT_plasma, pct.change = v_change, t = seq(0,v_time_final,v_time_int), power = v_power)
    
  
    #Bootstrapping for obtaining p-values and confidence intervals
    set.seed(12345)
    boot_results= boot(data = df, statistic = f_N, R = n_boot,
                       c_out0=c_out[i_out],strata = df$plasma_pos)#,seed=12345
    
  
    # ci_N_or=boot.ci(boot_results, type = "norm",seed=12345,index = 1)
    # ci_N_p=boot.ci(boot_results, type = "norm",seed=12345,index = 2)
    # # ci_N_MTL=boot.ci(boot_results, type = "norm",seed=12345,index = 3)
    # # ci_N_NeoT=boot.ci(boot_results, type = "norm",seed=12345,index = 4)
    # ci_N_pMTL=boot.ci(boot_results, type = "norm",seed=12345,index = 3)
    # ci_N_pNeoT=boot.ci(boot_results, type = "norm",seed=12345,index = 4)
    # 
    
  
    
    # results_under_H0 <- boot_results$t[,7] - mean(boot_results$t[,7])
    # p_perc_plasma0 <- mean(abs(results_under_H0) >= abs(boot_results$t0[7]))
    
    p_perc_plasma0 =(1 + sum(boot_results$t[,5] <= 0)) / (nrow(boot_results$t)+1) 
  
    # # results_under_H0 <- boot_results$t[,8] - mean(boot_results$t[,8])
    # # p_perc_MTL0 <- mean(abs(results_under_H0) >= abs(boot_results$t0[8]))
    # 
    # p_perc_MTL0 =(1 + sum(boot_results$t[,8] <= 0)) / (nrow(boot_results$t)+1) 
    # 
    # 
    # # results_under_H0 <- boot_results$t[,9] - mean(boot_results$t[,9])
    # # p_perc_NeoT0 <- mean(abs(results_under_H0) >= abs(boot_results$t0[9]))
    # p_perc_NeoT0 =(1 + sum(boot_results$t[,9] <= 0)) / (nrow(boot_results$t)+1) 
    
    # results_under_H0 <- boot_results$t[,10] - mean(boot_results$t[,10])
    # p_perc_pMTL0 <- mean(abs(results_under_H0) >= abs(boot_results$t0[10]))
    p_perc_pMTL0 =(1 + sum(boot_results$t[,6] <= 0)) / (nrow(boot_results$t)+1) 
    
    # results_under_H0 <- boot_results$t[,11] - mean(boot_results$t[,11])
    # p_perc_pNeoT0 <- mean(abs(results_under_H0) >= abs(boot_results$t0[11]))
    p_perc_pNeoT0 =(1 + sum(boot_results$t[,7] <= 0)) / (nrow(boot_results$t)+1) 
    
    # # results_under_H0 <- boot_results$t[,12] - mean(boot_results$t[,12])
    # # p_perc_MTL_plasma <- mean(abs(results_under_H0) >= abs(boot_results$t0[12]))
    # p_perc_MTL_plasma =(1 + sum(boot_results$t[,8] <= 0)) / (nrow(boot_results$t)+1) 
    # 
    # # results_under_H0 <- boot_results$t[,13] - mean(boot_results$t[,13])
    # # p_perc_NeoT_plasma <- mean(abs(results_under_H0) >= abs(boot_results$t0[13]))
    # p_perc_NeoT_plasma =(1 + sum(boot_results$t[,9] <= 0)) / (nrow(boot_results$t)+1) 
    # 
    # # results_under_H0 <- boot_results$t[,14] - mean(boot_results$t[,14])
    # # p_perc_MTL_NeoT <- mean(abs(results_under_H0) >= abs(boot_results$t0[14]))
    # p_perc_MTL_NeoT =(1 + sum(boot_results$t[,10] <= 0)) / (nrow(boot_results$t)+1) 
    
    # results_under_H0 <- boot_results$t[,15] - mean(boot_results$t[,15])
    # p_perc_plasma_MTLplasma <- mean(abs(results_under_H0) >= abs(boot_results$t0[15]))
    p_perc_plasma_MTLplasma =(1 + sum(boot_results$t[,8] <= 0)) / (nrow(boot_results$t)+1) 
    
    # # results_under_H0 <- boot_results$t[,16] - mean(boot_results$t[,16])
    # # p_perc_MTL_MTLplasma <- mean(abs(results_under_H0) >= abs(boot_results$t0[16]))
    # p_perc_MTL_MTLplasma =(1 + sum(boot_results$t[,16] <= 0)) / (nrow(boot_results$t)+1) 
    
    # results_under_H0 <- boot_results$t[,17] - mean(boot_results$t[,17])
    # p_perc_plasma_NeoTplasma <- mean(abs(results_under_H0) >= abs(boot_results$t0[17]))
    p_perc_plasma_NeoTplasma =(1 + sum(boot_results$t[,9] <= 0)) / (nrow(boot_results$t)+1) 
    
    # # results_under_H0 <- boot_results$t[,18] - mean(boot_results$t[,18])
    # # p_perc_NeoT_NeoTplasma <- mean(abs(results_under_H0) >= abs(boot_results$t0[18]))
    # p_perc_NeoT_NeoTplasma =(1 + sum(boot_results$t[,18] <= 0)) / (nrow(boot_results$t)+1) 
    
    # results_under_H0 <- boot_results$t[,24] - mean(boot_results$t[,24])
    # p_perc_pMTL_pNeoT <- mean(abs(results_under_H0) >= abs(boot_results$t0[24]))
    p_perc_pMTL_pNeoT =(1 + sum(boot_results$t[,13] <= 0)) / (nrow(boot_results$t)+1) 
    
    ci_perc_plasma0=boot.ci(boot_results, type = "norm",seed=12345,index = 10)
    # ci_perc_MTL0=boot.ci(boot_results, type = "norm",seed=12345,index = 20)
    # ci_perc_NeoT0=boot.ci(boot_results, type = "norm",seed=12345,index = 21)
    ci_perc_pMTL0=boot.ci(boot_results, type = "norm",seed=12345,index = 11)
    ci_perc_pNeoT0=boot.ci(boot_results, type = "norm",seed=12345,index = 12)
    ci_perc_pMTL_plasma=boot.ci(boot_results, type = "norm",seed=12345,index = 14)
    ci_perc_pNeoT_plasma=boot.ci(boot_results, type = "norm",seed=12345,index = 15)
    
   results_print$Outcome[i_out]=c_out[i_out]
   results_print$N_or[i_out]=paste(round(N_or$N,0))#,"[",round(ci_N_or[["normal"]][2]),", ",round(ci_N_or[["normal"]][3]),"]",sep = "")
   results_print$N_plasma[i_out]=paste(round(N_plasma$N,0))#,"[",round(ci_N_p[["normal"]][2]),", ",round(ci_N_p[["normal"]][3]),"]",sep = "")
   # results_print$N_MTL[i_out]=paste(round(N_MTL$N,0),"[",round(ci_N_MTL[["normal"]][2]),", ",round(ci_N_MTL[["normal"]][3]),"]",sep = "")
   # results_print$N_NeoT[i_out]=paste(round(N_NeoT$N,0),"[",round(ci_N_NeoT[["normal"]][2]),", ",round(ci_N_NeoT[["normal"]][3]),"]",sep = "")
   results_print$N_plasmaMTL[i_out]=paste(round(N_MTL_plasma$N,0))#,"[",round(ci_N_pMTL[["normal"]][2]),", ",round(ci_N_pMTL[["normal"]][3]),"]",sep = "")
   results_print$N_plasmaNeoT[i_out]=paste(round(N_NeoT_plasma$N,0))#,"[",round(ci_N_pNeoT[["normal"]][2]),", ",round(ci_N_pNeoT[["normal"]][3]),"]",sep = "")
   
   results_print$p_plasma_or[i_out]=paste(format(round(p_perc_plasma0,3),nsmall=3))
   # results_print$p_MTL_or[i_out]=paste(format(round(p_perc_MTL0,3),nsmall=3))
   # results_print$p_NeoT_or[i_out]=paste(format(round(p_perc_MTL0,3),nsmall=3))
   results_print$p_pMTL_or[i_out]=paste(format(round(p_perc_pMTL0,3),nsmall=3))
   results_print$p_pNeoT_or[i_out]=paste(format(round(p_perc_pNeoT0,3),nsmall=3))
   
   # results_print$p_MTL_plasma[i_out]=paste(format(round(p_perc_MTL_plasma,3),nsmall=3))
   # results_print$p_NeoT_plasma[i_out]=paste(format(round(p_perc_NeoT_plasma,3),nsmall=3))
   # results_print$p_MTL_NeoT[i_out]=paste(format(round(p_perc_MTL_NeoT,3),nsmall=3))
   
   
   results_print$p_plasma_pMTL[i_out]=paste(format(round(p_perc_plasma_MTLplasma,3),nsmall=3))
   results_print$p_plasma_pNeoT[i_out]=paste(format(round(p_perc_plasma_NeoTplasma,3),nsmall=3))
   
   # results_print$p_MTL_pMTL[i_out]=paste(format(round(p_perc_MTL_MTLplasma,3),nsmall=3))
   # results_print$p_NeoT_pNeoT[i_out]=paste(format(round(p_perc_NeoT_NeoTplasma,3),nsmall=3))
   
   results_print$p_pMTL_pNeoT[i_out]=paste(format(round(p_perc_pMTL_pNeoT,3),nsmall=3))
   
   c_sign=NA
   v_p=c(#p_perc_MTL_plasma,p_perc_MTL_NeoT,p_perc_NeoT_plasma,
         p_perc_pMTL_pNeoT,#p_perc_MTL_MTLplasma,
         p_perc_plasma_MTLplasma,#p_perc_NeoT_NeoTplasma,
         p_perc_plasma_NeoTplasma)
   for (i_sign in 1:length(v_p)) {
     if (v_p[i_sign]<0.05) {
       c_sign[i_sign]="*"
     }
     if (v_p[i_sign]<0.01) {
       c_sign[i_sign]="**"
     }
     if (v_p[i_sign]<0.05) {
       c_sign[i_sign]="*"
     }
     if (v_p[i_sign]<0.001) {
       c_sign[i_sign]="***"
     }
     if (v_p[i_sign]>0.05) {
       c_sign[i_sign]=NA
     }
   }
   
   results_print$perc_plasma[i_out]=paste(round(100*N_plasma$N/N_or$N,0),"[",round(ci_perc_plasma0[["normal"]][2]),", ",round(ci_perc_plasma0[["normal"]][3]),"]",sep = "")
   # results_print$perc_MTL[i_out]=paste(round(100*N_MTL$N/N_or$N,0),"[",round(ci_perc_MTL0[["normal"]][2]),", ",round(ci_perc_MTL0[["normal"]][3]),"]",sep = "")
   # results_print$perc_NeoT[i_out]=paste(round(100*N_NeoT$N/N_or$N,0),"[",round(ci_perc_NeoT0[["normal"]][2]),", ",round(ci_perc_NeoT0[["normal"]][3]),"]",sep = "")
   results_print$perc_plasmaMTL[i_out]=paste(round(100*N_MTL_plasma$N/N_or$N,0),"[",round(ci_perc_pMTL0[["normal"]][2]),", ",round(ci_perc_pMTL0[["normal"]][3]),"]",sep = "")
   results_print$perc_plasmaNeoT[i_out]=paste(round(100*N_NeoT_plasma$N/N_or$N,0),"[",round(ci_perc_pNeoT0[["normal"]][2]),", ",round(ci_perc_pNeoT0[["normal"]][3]),"]",sep = "")
   results_print$perc_plasmaMTL_plasma[i_out]=paste(round(100*N_MTL_plasma$N/N_plasma$N,0),"[",round(ci_perc_pMTL_plasma[["normal"]][2]),", ",round(ci_perc_pMTL_plasma[["normal"]][3]),"]",sep = "")
   results_print$perc_plasmaNeoT_plasma[i_out]=paste(round(100*N_NeoT_plasma$N/N_plasma$N,0),"[",round(ci_perc_pNeoT_plasma[["normal"]][2]),", ",round(ci_perc_pNeoT_plasma[["normal"]][3]),"]",sep = "")
   
   
   # results_plot$Model[1]="Original"
   results_plot$Model[1]="Plasma"
   # results_plot$Model[2]="MTL"
   # results_plot$Model[3]="NeoT"
   results_plot$Model[2]="p-MTL"
   results_plot$Model[3]="p-NeoT"
   results_plot$Model[4]="p-MTL_plasma"
   results_plot$Model[5]="p-NeoT_plasma"
   
   # results_plot$Perc[1]=100
   # results_plot$Perc_l[1]=NA
   # results_plot$Perc_h[1]=NA
   
   results_plot$Perc[1]=100*(N_plasma$N)/N_or$N#-100*(N_plasma$N-N_or$N)/N_or$N
   results_plot$Perc_l[1]=ci_perc_plasma0[["normal"]][2]
   results_plot$Perc_h[1]=ci_perc_plasma0[["normal"]][3]
   
   # results_plot$Perc[2]=100*(N_MTL$N)/N_or$N#-100*(N_MTL$N-N_or$N)/N_or$N
   # results_plot$Perc_l[2]=ci_perc_MTL0[["normal"]][2]
   # results_plot$Perc_h[2]=ci_perc_MTL0[["normal"]][3]
   # 
   # results_plot$Perc[3]=100*(N_NeoT$N)/N_or$N#-100*(N_NeoT$N-N_or$N)/N_or$N
   # results_plot$Perc_l[3]=ci_perc_NeoT0[["normal"]][2]
   # results_plot$Perc_h[3]=ci_perc_NeoT0[["normal"]][3]
   
   results_plot$Perc[2]=100*(N_MTL_plasma$N)/N_or$N#-100*(N_MTL_plasma$N-N_or$N)/N_or$N
   results_plot$Perc_l[2]=ci_perc_pMTL0[["normal"]][2]
   results_plot$Perc_h[2]=ci_perc_pMTL0[["normal"]][3]
   
   results_plot$Perc[3]=100*(N_NeoT_plasma$N)/N_or$N#-100*(N_NeoT_plasma$N-N_or$N)/N_or$N
   results_plot$Perc_l[3]=ci_perc_pNeoT0[["normal"]][2]
   results_plot$Perc_h[3]=ci_perc_pNeoT0[["normal"]][3]
   
   results_plot$Perc[4]=100*(N_MTL_plasma$N)/N_plasma$N#-100*(N_MTL_plasma$N-N_or$N)/N_or$N
   results_plot$Perc_l[4]=ci_perc_pMTL_plasma[["normal"]][2]
   results_plot$Perc_h[4]=ci_perc_pMTL_plasma[["normal"]][3]
   
   results_plot$Perc[5]=100*(N_NeoT_plasma$N)/N_plasma$N#-100*(N_NeoT_plasma$N-N_or$N)/N_or$N
   results_plot$Perc_l[5]=ci_perc_pNeoT_plasma[["normal"]][2]
   results_plot$Perc_h[5]=ci_perc_pNeoT_plasma[["normal"]][3]
   
   results_plot=results_plot %>% filter(!is.na(Model))
   
   results_line=data.frame(matrix(nrow = 6,ncol = 5))
   colnames(results_line)=c("Model","step",
                            "N","N_l","N_h")
   
   results_line$Model=c(rep("MTL",3),rep("NeoT",3))
   results_line$step=c(0,1,2,0,1,2)
   results_line$N=c(100,100*N_plasma$N/N_or$N,100*N_MTL_plasma$N/N_or$N,100,100*N_plasma$N/N_or$N,100*N_NeoT_plasma$N/N_or$N)
   results_line$N_l=c(NA,ci_perc_plasma0[["normal"]][2],ci_perc_pMTL0[["normal"]][2],
                      NA,ci_perc_plasma0[["normal"]][2],ci_perc_pNeoT0[["normal"]][2])
   results_line$N_h=c(NA,ci_perc_plasma0[["normal"]][3],ci_perc_pMTL0[["normal"]][3],
                      NA,ci_perc_plasma0[["normal"]][3],ci_perc_pNeoT0[["normal"]][3])
   results_line$Model=factor(results_line$Model,levels=c("MTL","NeoT"))
   results_line$step=as.factor(results_line$step)
   
   # plot_line=ggplot(results_line,aes(x=step,y=N,group=Model,fill=Model,color=Model))+
   #   geom_line()+   geom_errorbar(aes(ymin=N_l, ymax=N_h), width=.2,alpha=0.5)+
   #   geom_point(size = 2, shape = 21) +
   #   geom_text(aes(label = paste(round(N,0)),x=step,y=N+5), stat="identity", vjust=0.1,size=4,color="black")+
   #   theme_classic()+
   #   scale_fill_manual(values = c("MTL"="#00a087ff",
   #                                "NeoT"="#3c5488ff"))+
   #   scale_color_manual(values = c("MTL"="#00a087ff",
   #                                 "NeoT"="#3c5488ff"))+
   #   theme(legend.position = "none",axis.text=element_text(size=8),
   #                          axis.title=element_text(size=10,face="bold"),legend.text = element_text(size=6),
   #                          axis.line.x = element_line(size = 0.3),axis.line.y = element_line(size = 0.3),
   #                          axis.ticks.x=element_line(size = 0.3),axis.ticks.y = element_line(size = 0.3),
   #                          legend.key.size = unit(0.2, 'cm'),legend.title = element_blank())+
   #                              geom_hline(yintercept=100,linetype="dotted")+
   #                          scale_y_continuous(limits=c(0,100))
   # 
   # pf2=ggarrange(plot_line+rremove("ylab"),
   #               labels = c("c"),
   #               ncol = 3, nrow = 3)
   # ggexport(pf2,filename=paste(dir_fig,"/SampleRed_line_",c_out[i_out],"_Plasma_",c_method_plasma[i_method],"_PET_",c_method_PET[i_method2],"_Ab_",l_ab,"_adj_",l_adj,"_power_",v_power,"_change_",v_change,".pdf",sep = ""))
   
   # p6=ggplot(results_plot, aes(x = factor(Model,levels = c("Plasma","p-MTL","p-NeoT")), y = Perc,fill=Model)) +
   #   geom_bar(stat="identity", position=position_dodge()) +
   #   geom_signif(comparisons = list(#c("Plasma","MTL"),c("MTL","NeoT"),c("Plasma","NeoT"),
   #                                  c("p-MTL","p-NeoT"),#c("MTL","p-MTL"),
   #                                  c("Plasma","p-MTL"),
   #                                  #c("NeoT","p-NeoT"),
   #                                  c("Plasma","p-NeoT")),annotations = c_sign,step_increase = 0.1,y_position = max(results_plot$Perc))+
   #   geom_errorbar(aes(ymin=Perc_l, ymax=Perc_h), width=.2,
   #                 position=position_dodge(.9))+
   #    scale_fill_manual(values=c_colors)+ylab("Final sample (%)") + xlab("Model")+
   #   # geom_text(aes(label = paste(round(Perc,0),"%",sep = ""),x=Model,y=Perc_h+25), stat="identity", vjust=0.1,size=4)+
   #   theme_classic()+theme(legend.position = "none",axis.text=element_text(size=8),
   #                         axis.title=element_text(size=10,face="bold"),legend.text = element_text(size=6),
   #                         axis.line.x = element_line(size = 0.3),axis.line.y = element_line(size = 0.3),
   #                         axis.ticks.x=element_line(size = 0.3),axis.ticks.y = element_line(size = 0.3),
   #                         legend.key.size = unit(0.2, 'cm'),legend.title = element_blank())+
   #      geom_hline(yintercept=100,linetype="dotted")#+
   #   # scale_y_continuous(limits=c(0,130))
   # 
   # 
   # ggexport(p6,filename=paste(dir_fig,"/SampleRed_",c_out[i_out],"_",c_method[i_method],"_Ab_",l_ab,"_adj_",l_adj,"_power_",v_power,"_change_",v_change,".pdf",sep = ""))
   
   writexl::write_xlsx(results_plot,path = paste(dir_res,"/Plot_sampleRed_",c_out[i_out],"_Plasma_",c_method_plasma[i_method],"_PET_",c_method_PET[i_method2],"_Ab_",l_ab,"_adj_",l_adj,"_power_",v_power,"_change_",v_change,"_proportion_time_",v_time_final,"_interval_",v_time_int,"_20241010.xlsx",sep = ""))
   
  }
  
  writexl::write_xlsx(results_print,path = paste(dir_res,"/SampleRed_Plasma_",c_method_plasma[i_method],"_PET_",c_method_PET[i_method2],"_Ab_",l_ab,"_adj_",l_adj,"_power_",v_power,"_change_",v_change,"_proportion_time_",v_time_final,"_interval_",v_time_int,"_20241010.xlsx",sep = ""))
  }
}
