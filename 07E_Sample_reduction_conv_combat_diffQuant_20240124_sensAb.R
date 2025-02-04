  library(powerSurvEpi)
  library(boot)
  dir_res="~/Desktop/Lund/20_TauPET_plasma/Results_Merged/20240917_Rev/Two_step_sample_red_conv/"
  dir_fig="~/Desktop/Lund/20_TauPET_plasma/Figures_Merged/20240917_Rev/Two_step_sample_red_conv/"
  load("~/Desktop/Lund/20_TauPET_plasma/Data_all/Data_all_cogn_comb_20240917.RData")
  
  dataAll$time_MCI[which(dataAll$cohort=="WU" | dataAll$cohort=="WRAP")]=dataAll$time_MCI[which(dataAll$cohort=="WU" | dataAll$cohort=="WRAP")]*365.25
  dataAll=dataAll %>% filter(!is.na(time_MCI))
  dataAll=dataAll %>% filter(!duplicated(sid))
  
  f_N <- function(data, indices) {
    d <- data[indices,] # allows boot to select sample
    
    v_perc_conv0=sum(d$conv_MCI==1)/sum(!is.na(d$conv_MCI))
    v_perc_conv_plasma=sum(d$conv_MCI[d$plasma_pos==1]==1)/sum(!is.na(d$conv_MCI[d$plasma_pos==1]))
    # Added:
        v_perc_conv_MTL=sum(d$conv_MCI[d$plasma_pos==1 & d$ab_pos==1]==1)/sum(!is.na(d$conv_MCI[d$plasma_pos==1 & d$ab_pos==1]))
    # v_perc_conv_NeoT=sum(d$conv_MCI[d$NeoT_pos==1]==1)/sum(!is.na(d$conv_MCI[d$NeoT_pos==1]))
    v_perc_conv_pMTL=sum(d$conv_MCI[d$MTL_pos==1 & d$plasma_pos==1 & d$ab_pos==1]==1)/sum(!is.na(d$conv_MCI[d$MTL_pos==1 & d$plasma_pos==1 & d$ab_pos==1]))
    v_perc_conv_pNeoT=sum(d$conv_MCI[d$NeoT_pos==1 & d$plasma_pos==1 & d$ab_pos==1]==1)/sum(!is.na(d$conv_MCI[d$NeoT_pos==1 & d$plasma_pos==1 & d$ab_pos==1]))
    
    v_N0=ssizeCT.default(power=v_power0,
                         k=1,
                         pE=v_redPerc*v_perc_conv0, #this should be the percentage of converters in the ”treated group”, put that at 70% of the observed in the whole data set
                         pC=v_perc_conv0, #should be the percentage of converters in the placebo group, putt hat at the observed % in your whole data set  
                         RR = exp ( log ( log (v_redPerc*v_perc_conv0)/log(v_perc_conv0) ) ),#exp ( log ( log (pE)/log(pC) ) ), # see below
                         alpha = 0.05)
    
    v_Nplasma=ssizeCT.default(power=v_power0,
                              k=1,
                              pE=v_redPerc*v_perc_conv_plasma, #this should be the percentage of converters in the ”treated group”, put that at 70% of the observed in the whole data set
                              pC=v_perc_conv_plasma, #should be the percentage of converters in the placebo group, putt hat at the observed % in your whole data set  
                              RR = exp ( log ( log (v_redPerc*v_perc_conv_plasma)/log(v_perc_conv_plasma) ) ),#exp ( log ( log (pE)/log(pC) ) ), # see below
                              alpha = 0.05)
    
    v_NMTL=ssizeCT.default(power=v_power0,
                           k=1,
                           pE=v_redPerc*v_perc_conv_MTL, #this should be the percentage of converters in the ”treated group”, put that at 70% of the observed in the whole data set
                           pC=v_perc_conv_MTL, #should be the percentage of converters in the placebo group, putt hat at the observed % in your whole data set
                           RR = exp ( log ( log (v_redPerc*v_perc_conv_MTL)/log(v_perc_conv_MTL) ) ),#exp ( log ( log (pE)/log(pC) ) ), # see below
                           alpha = 0.05)
    # 
    # v_NNeoT=ssizeCT.default(power=v_power0,
    #                         k=1,
    #                         pE=v_redPerc*v_perc_conv_NeoT, #this should be the percentage of converters in the ”treated group”, put that at 70% of the observed in the whole data set
    #                         pC=v_perc_conv_NeoT, #should be the percentage of converters in the placebo group, putt hat at the observed % in your whole data set  
    #                         RR = exp ( log ( log (v_redPerc*v_perc_conv_NeoT)/log(v_perc_conv_NeoT) ) ),#exp ( log ( log (pE)/log(pC) ) ), # see below
    #                         alpha = 0.05)
    
    v_NpMTL=ssizeCT.default(power=v_power0,
                            k=1,
                            pE=v_redPerc*v_perc_conv_pMTL, #this should be the percentage of converters in the ”treated group”, put that at 70% of the observed in the whole data set
                            pC=v_perc_conv_pMTL, #should be the percentage of converters in the placebo group, putt hat at the observed % in your whole data set  
                            RR = exp ( log ( log (v_redPerc*v_perc_conv_pMTL)/log(v_perc_conv_pMTL) ) ),#exp ( log ( log (pE)/log(pC) ) ), # see below
                            alpha = 0.05)
    
    v_NpNeoT=ssizeCT.default(power=v_power0,
                             k=1,
                             pE=v_redPerc*v_perc_conv_pNeoT, #this should be the percentage of converters in the ”treated group”, put that at 70% of the observed in the whole data set
                             pC=v_perc_conv_pNeoT, #should be the percentage of converters in the placebo group, putt hat at the observed % in your whole data set  
                             RR = exp ( log ( log (v_redPerc*v_perc_conv_pNeoT)/log(v_perc_conv_pNeoT) ) ),#exp ( log ( log (pE)/log(pC) ) ), # see below
                             alpha = 0.05)
    
    
    p_plasma0=(v_N0[[1]]-v_Nplasma[[1]])
    p_MTL0=(v_N0[[1]]-v_NMTL[[1]])#Step with Abpet after plasma and before Tau-PET
    # p_NeoT0=(v_N0[[1]]-v_NNeoT[[1]])
    p_pMTL0=(v_N0[[1]]-v_NpMTL[[1]])
    p_pNeoT0=(v_N0[[1]]-v_NpNeoT[[1]])
    p_MTL_plasma=-(v_NMTL[[1]]-v_Nplasma[[1]])#Step with Abpet after plasma and before Tau-PET
    # p_NeoT_plasma=-(v_NNeoT[[1]]-v_Nplasma[[1]])
    # p_MTL_NeoT=(v_NNeoT[[1]]-v_NMTL[[1]])
    p_MTL_MTLplasma=(v_NMTL[[1]]-v_NpMTL[[1]])#Step with Abpet after plasma and before Tau-PET
    p_plasma_MTLplasma=(v_Nplasma[[1]]-v_NpMTL[[1]])
    p_plasma_NeoTplasma=(v_Nplasma[[1]]-v_NpNeoT[[1]])
    # p_NeoT_NeoTplasma=(v_NNeoT[[1]]-v_NpNeoT[[1]])
    p_pMTL_pNeoT=-(v_NpMTL[[1]]-v_NpNeoT[[1]])#ULL!
    
    return(c(100*v_Nplasma[[1]]/v_N0[[1]],#100*v_NMTL[[1]]/v_N0[[1]],100*v_NNeoT[[1]]/v_N0[[1]],
             100*v_NpMTL[[1]]/v_N0[[1]],100*v_NpNeoT[[1]]/v_N0[[1]],
             p_plasma0,#p_MTL0,p_NeoT0,
             p_pMTL0,p_pNeoT0,
             #p_MTL_plasma,p_NeoT_plasma,p_MTL_NeoT,
             p_plasma_MTLplasma,#p_MTL_MTLplasma,
             p_plasma_NeoTplasma,#p_NeoT_NeoTplasma,
             100*v_Nplasma[[1]]/v_N0[[1]],
             100*v_NpMTL[[1]]/v_N0[[1]],100*v_NpNeoT[[1]]/v_N0[[1]],p_pMTL_pNeoT,
             #Mod:
             100*v_NpMTL[[1]]/v_NMTL[[1]],100*v_NpNeoT[[1]]/v_NMTL[[1]],
             #New
             v_NMTL[[1]],p_pMTL0,p_MTL_plasma,p_plasma_MTLplasma,100*v_NMTL[[1]]/v_Nplasma[[1]],100*v_NMTL[[1]]/v_N0[[1]]))
    }
  
  dataAll$cohort=as.factor(dataAll$cohort)
  n_boot=1000
  
  c_colors=c("Basic_NA"="#d9d9d9",
             "Basic_A"="#737373",#"Original"="lightgrey",
             "Plasma"="#e64b35ff",
             # "MTL"="#a1d99b",
             # "NeoT"="#31a354",
             "p-Ab"="#4DBBD5FF",
             "p-MTL"="#00a087ff",
             "p-NeoT"="#3c5488ff")
  
  
  v_power0=0.8
  v_redPerc=0.7
  l_ab="all"#"Abpos"#
  l_adj="raw"#"combat"#
  c_method_plasma=c("Q2_Q4","Q3_Q4","Q4")#"zscore2"
  c_method_PET=c("Q2_Q4","Q3_Q4","Q4")#"zscore2"
  
  n_method=length(c_method_plasma)
  v_method_plasma=c(0.25,0.5,0.75)
  v_method_PET=c(0.25,0.5,0.75)
  
  results_print=data.frame(matrix(nrow = n_method,ncol = 16+6))
  colnames(results_print)=c("Outcome",
                            "N_or","N_plasma","N_MTL",#"N_NeoT",
                            "N_plasmaMTL","N_plasmaNeoT",
                            "p_plasma_or","p_MTL_or",#"p_NeoT_or",
                            "p_pMTL_or","p_pNeoT_or",
                            "p_MTL_plasma",#"p_NeoT_plasma","p_MTL_NeoT",
                            "p_plasma_pMTL","p_MTL_pMTL",
                            "p_plasma_pNeoT",#"p_NeoT_pNeoT",
                            "p_pMTL_pNeoT",
                            "perc_plasma","perc_MTL0","perc_MTL",#"perc_NeoT",
                            "perc_plasmaMTL","perc_plasmaNeoT",
                            "perc_plasmaMTL_plasma","perc_plasmaNeoT_plasma")
  
  
  for (i_method in 1:length(c_method_plasma)) {
    
    for (i_method2 in 1:length(c_method_PET)) {
      
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
      
      df=dataAll[complete.cases(dataAll[,c("conv_MCI","z_ptau217","z_MTL","z_NeoT")]),]
      
      if (l_ab=="Abpos") {
        df=df %>% filter(ab_pos==1)
      }
      
      if (l_adj=="combat") {
        v_PET=c(15,14)#c(244:247)
        v_pl=13
        
      }else{
        v_PET=c(10,9)#c(244:247)
        v_pl=8
          }
      c_pl=names(df)[v_pl]
      c_PET=names(df)[v_PET]
      n_PET=length(v_PET)
      
      
      df$plasma_pos=as.factor(ifelse(df[,c_pl]>quantile(df[!duplicated(df$sid),c_pl],v_method_plasma[i_method]),1,0))#as.factor(ifelse(df$z_ptau217_comb>1.96,1,0))
      
      df2=df %>% filter(!duplicated(sid))
      df2=df2 %>% filter(plasma_pos==1)
      
      #We select the quartiles of those participants already plasma positive to assess the additive reduction of performing a PET after being plasma+
      df$MTL_pos=as.factor(ifelse(df[,c_PET[2]]>quantile(df2[,c_PET[2]],v_method_PET[i_method2]),1,0))
      df$NeoT_pos=as.factor(ifelse(df[,c_PET[1]]>quantile(df2[,c_PET[1]],v_method_PET[i_method2]),1,0))
      
      v_perc_conv0=sum(df$conv_MCI==1)/sum(!is.na(df$conv_MCI))
      v_perc_conv_plasma=sum(df$conv_MCI[df$plasma_pos==1]==1)/sum(!is.na(df$conv_MCI[df$plasma_pos==1]))
      v_perc_conv_MTL=sum(df$conv_MCI[df$plasma_pos==1 & df$ab_pos==1]==1)/sum(!is.na(df$conv_MCI[df$plasma_pos==1 & df$ab_pos==1]))
      # v_perc_conv_NeoT=sum(df$conv_MCI[df$NeoT_pos==1]==1)/sum(!is.na(df$conv_MCI[df$NeoT_pos==1]))
      v_perc_conv_pMTL=sum(df$conv_MCI[df$MTL_pos==1 & df$plasma_pos==1 & df$ab_pos==1]==1)/sum(!is.na(df$conv_MCI[df$MTL_pos==1 & df$plasma_pos==1 & df$ab_pos==1]))
      v_perc_conv_pNeoT=sum(df$conv_MCI[df$NeoT_pos==1 & df$plasma_pos==1 & df$ab_pos==1]==1)/sum(!is.na(df$conv_MCI[df$NeoT_pos==1 & df$plasma_pos==1 & df$ab_pos==1]))
      
      v_N0=ssizeCT.default(power=v_power0,
                           k=1,
                           pE=v_redPerc*v_perc_conv0, #this should be the percentage of converters in the ”treated group”, put that at 70% of the observed in the whole data set
                           pC=v_perc_conv0, #should be the percentage of converters in the placebo group, putt hat at the observed % in your whole data set  
                           RR = exp ( log ( log (v_redPerc*v_perc_conv0)/log(v_perc_conv0) ) ),#exp ( log ( log (pE)/log(pC) ) ), # see below
                           alpha = 0.05)
      
      v_Nplasma=ssizeCT.default(power=v_power0,
                                k=1,
                                pE=v_redPerc*v_perc_conv_plasma, #this should be the percentage of converters in the ”treated group”, put that at 70% of the observed in the whole data set
                                pC=v_perc_conv_plasma, #should be the percentage of converters in the placebo group, putt hat at the observed % in your whole data set  
                                RR = exp ( log ( log (v_redPerc*v_perc_conv_plasma)/log(v_perc_conv_plasma) ) ),#exp ( log ( log (pE)/log(pC) ) ), # see below
                                alpha = 0.05)
      
      v_NMTL=ssizeCT.default(power=v_power0,
                             k=1,
                             pE=v_redPerc*v_perc_conv_MTL, #this should be the percentage of converters in the ”treated group”, put that at 70% of the observed in the whole data set
                             pC=v_perc_conv_MTL, #should be the percentage of converters in the placebo group, putt hat at the observed % in your whole data set
                             RR = exp ( log ( log (v_redPerc*v_perc_conv_MTL)/log(v_perc_conv_MTL) ) ),#exp ( log ( log (pE)/log(pC) ) ), # see below
                             alpha = 0.05)
      # 
      # v_NNeoT=ssizeCT.default(power=v_power0,
      #                         k=1,
      #                         pE=v_redPerc*v_perc_conv_NeoT, #this should be the percentage of converters in the ”treated group”, put that at 70% of the observed in the whole data set
      #                         pC=v_perc_conv_NeoT, #should be the percentage of converters in the placebo group, putt hat at the observed % in your whole data set  
      #                         RR = exp ( log ( log (v_redPerc*v_perc_conv_NeoT)/log(v_perc_conv_NeoT) ) ),#exp ( log ( log (pE)/log(pC) ) ), # see below
      #                         alpha = 0.05)
      
      v_NpMTL=ssizeCT.default(power=v_power0,
                              k=1,
                              pE=v_redPerc*v_perc_conv_pMTL, #this should be the percentage of converters in the ”treated group”, put that at 70% of the observed in the whole data set
                              pC=v_perc_conv_pMTL, #should be the percentage of converters in the placebo group, putt hat at the observed % in your whole data set  
                              RR = exp ( log ( log (v_redPerc*v_perc_conv_pMTL)/log(v_perc_conv_pMTL) ) ),#exp ( log ( log (pE)/log(pC) ) ), # see below
                              alpha = 0.05)
      
      v_NpNeoT=ssizeCT.default(power=v_power0,
                               k=1,
                               pE=v_redPerc*v_perc_conv_pNeoT, #this should be the percentage of converters in the ”treated group”, put that at 70% of the observed in the whole data set
                               pC=v_perc_conv_pNeoT, #should be the percentage of converters in the placebo group, putt hat at the observed % in your whole data set  
                               RR = exp ( log ( log (v_redPerc*v_perc_conv_pNeoT)/log(v_perc_conv_pNeoT) ) ),#exp ( log ( log (pE)/log(pC) ) ), # see below
                               alpha = 0.05)
      
     
      
      set.seed(12345)
      boot_results= boot(data = df, statistic = f_N, R = n_boot)
      
      # 1 100*v_Nplasma[[1]]/v_N0[[1]],#100*v_NMTL[[1]]/v_N0[[1]],100*v_NNeoT[[1]]/v_N0[[1]],
      # 2-3 100*v_NpMTL[[1]]/v_N0[[1]],100*v_NpNeoT[[1]]/v_N0[[1]],
      # 4 p_plasma0,#p_MTL0,p_NeoT0,
      # 5-6 p_pMTL0,p_pNeoT0,
      # #p_MTL_plasma,p_NeoT_plasma,p_MTL_NeoT,
      # 7 p_plasma_MTLplasma,#p_MTL_MTLplasma,
      # 8 p_plasma_NeoTplasma,#p_NeoT_NeoTplasma,
      # 9 100*v_Nplasma[[1]]/v_N0[[1]],
      # 10-12 100*v_NpMTL[[1]]/v_N0[[1]],100*v_NpNeoT[[1]]/v_N0[[1]],p_pMTL_pNeoT,
      # #Mod:
      # 13-14 100*v_NpMTL[[1]]/v_NMTL[[1]],100*v_NpNeoT[[1]]/v_NMTL[[1]],
      # #New
      # 15 v_NMTL[[1]],
      # 16-18 p_pMTL0,p_MTL_plasma,p_plasma_MTLplasma,
      # 19-20 100*v_NMTL[[1]]/v_Nplasma[[1]],100*v_NMTL[[1]]/v_N0[[1]]))
  
      
      #New
      # N_MTL$N,p_MTL0,p_MTL_plasma,p_MTL_MTLplasma,100*N_MTL$N/N_or$N
      p_perc_MTL0 =(1 + sum(boot_results$t[,16] <= 0)) / (nrow(boot_results$t)+1)
      p_perc_MTL_plasma =(1 + sum(boot_results$t[,17] <= 0)) / (nrow(boot_results$t)+1)
      p_perc_MTL_MTLplasma =(1 + sum(boot_results$t[,18] <= 0)) / (nrow(boot_results$t)+1)
      ci_perc_Ab_plasma=boot.ci(boot_results, type = "norm",seed=12345,index = 19)#Compared to plasma+
      ci_perc_Ab0=boot.ci(boot_results, type = "norm",seed=12345,index = 20)#Compared to plasma+
      
      
      p_perc_plasma0 =(1 + sum(!is.na(boot_results$t[,4]) <= 0)) / (nrow(boot_results$t)+1) 
      # p_perc_MTL0 =(1 + sum(boot_results$t[,7] <= 0)) / (nrow(boot_results$t)+1) 
      # p_perc_NeoT0 =(1 + sum(boot_results$t[,8] <= 0)) / (nrow(boot_results$t)+1) 
      p_perc_pMTL0 =(1 + sum(!is.na(boot_results$t[,5]) <= 0)) / (nrow(boot_results$t)+1) 
      p_perc_pNeoT0 =(1 + sum(!is.na(boot_results$t[,6]) <= 0)) / (nrow(boot_results$t)+1) 
      # p_perc_MTL_plasma =(1 + sum(boot_results$t[,11] <= 0)) / (nrow(boot_results$t)+1) 
      # p_perc_NeoT_plasma =(1 + sum(boot_results$t[,12] <= 0)) / (nrow(boot_results$t)+1) 
      # p_perc_MTL_NeoT =(1 + sum(boot_results$t[,13] <= 0)) / (nrow(boot_results$t)+1) 
      p_perc_plasma_MTLplasma =(1 + sum(!is.na(boot_results$t[,7]) <= 0)) / (nrow(boot_results$t)+1) 
      # p_perc_MTL_MTLplasma =(1 + sum(boot_results$t[,15] <= 0)) / (nrow(boot_results$t)+1)
      p_perc_plasma_NeoTplasma =(1 + sum(!is.na(boot_results$t[,8]) <= 0)) / (nrow(boot_results$t)+1) 
      # p_perc_NeoT_NeoTplasma =(1 + sum(boot_results$t[,17] <= 0)) / (nrow(boot_results$t)+1) 
      p_perc_pMTL_pNeoT =(1 + sum(!is.na(boot_results$t[,12]) <= 0)) / (nrow(boot_results$t)+1) 
      
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
      
      
      ci_perc_plasma0=boot.ci(boot_results, type = "norm",seed=12345,index = 9)
      # ci_perc_MTL0=boot.ci(boot_results, type = "norm",seed=12345,index = 2)
      # ci_perc_NeoT0=boot.ci(boot_results, type = "norm",seed=12345,index = 3)
      ci_perc_pMTL0=boot.ci(boot_results, type = "norm",seed=12345,index = 10)
      ci_perc_pNeoT0=boot.ci(boot_results, type = "norm",seed=12345,index = 11)  
  
      ci_perc_pMTL_plasma=boot.ci(boot_results, type = "norm",seed=12345,index = 13)
      ci_perc_pNeoT_plasma=boot.ci(boot_results, type = "norm",seed=12345,index = 14)  
      
      results_print$Outcome[i_method2]=c_method_PET[i_method2]
      results_print$N_or[i_method2]=paste(round(v_N0[[1]],0))
      results_print$N_plasma[i_method2]=paste(round(v_Nplasma[[1]],0))
      results_print$N_MTL[i_method]=paste(round(v_NMTL[[1]],0))#Mod
      # results_print$N_NeoT[i_method]=paste(round(v_NNeoT[[1]],0))
      results_print$N_plasmaMTL[i_method2]=paste(round(v_NpMTL[[1]],0))
      results_print$N_plasmaNeoT[i_method2]=paste(round(v_NpNeoT[[1]],0))
      
      results_print$p_plasma_or[i_method2]=paste(format(round(p_perc_plasma0,3),nsmall=3))
      results_print$p_MTL_or[i_method]=paste(format(round(p_perc_MTL0,3),nsmall=3))
      # results_print$p_NeoT_or[i_method]=paste(format(round(p_perc_MTL0,3),nsmall=3))
      results_print$p_pMTL_or[i_method2]=paste(format(round(p_perc_pMTL0,3),nsmall=3))
      results_print$p_pNeoT_or[i_method2]=paste(format(round(p_perc_pNeoT0,3),nsmall=3))
      # results_print$p_MTL_plasma[i_method]=paste(format(round(p_perc_MTL_plasma,3),nsmall=3))
      # results_print$p_NeoT_plasma[i_method]=paste(format(round(p_perc_NeoT_plasma,3),nsmall=3))
      # results_print$p_MTL_NeoT[i_method]=paste(format(round(p_perc_MTL_NeoT,3),nsmall=3))
      results_print$p_plasma_pMTL[i_method2]=paste(format(round(p_perc_plasma_MTLplasma,3),nsmall=3))
      results_print$p_plasma_pNeoT[i_method2]=paste(format(round(p_perc_plasma_NeoTplasma,3),nsmall=3))
      # results_print$p_MTL_pMTL[i_method]=paste(format(round(p_perc_MTL_MTLplasma,3),nsmall=3))
      # results_print$p_NeoT_pNeoT[i_method]=paste(format(round(p_perc_NeoT_NeoTplasma,3),nsmall=3))
      results_print$p_pMTL_pNeoT[i_method2]=paste(format(round(p_perc_pMTL_pNeoT,3),nsmall=3))
      
      results_print$perc_plasma[i_method2]=paste(round(100*v_Nplasma[[1]]/v_N0[[1]],0),"[",round(ci_perc_plasma0[["normal"]][2]),", ",round(ci_perc_plasma0[["normal"]][3]),"]",sep = "")
      results_print$perc_MTL[i_method]=paste(round(100*v_NMTL[[1]]/v_Nplasma[[1]],0),"[",round(ci_perc_Ab_plasma[["normal"]][2]),", ",round(ci_perc_Ab_plasma[["normal"]][3]),"]",sep = "")#Mod:
      # results_print$perc_NeoT[i_method]=paste(round(100*v_NNeoT[[1]]/v_N0[[1]],0),"[",round(ci_perc_NeoT0[["normal"]][2]),", ",round(ci_perc_NeoT0[["normal"]][3]),"]",sep = "")
      results_print$perc_plasmaMTL[i_method2]=paste(round(100*v_NpMTL[[1]]/v_N0[[1]],0),"[",round(ci_perc_pMTL0[["normal"]][2]),", ",round(ci_perc_pMTL0[["normal"]][3]),"]",sep = "")#Compared to initial participants
      results_print$perc_plasmaNeoT[i_method2]=paste(round(100*v_NpNeoT[[1]]/v_N0[[1]],0),"[",round(ci_perc_pNeoT0[["normal"]][2]),", ",round(ci_perc_pNeoT0[["normal"]][3]),"]",sep = "")#Compared to initial participants
      results_print$perc_plasmaMTL_plasma[i_method2]=paste(round(100*v_NpMTL[[1]]/v_NMTL[[1]],0),"[",round(ci_perc_pMTL_plasma[["normal"]][2]),", ",round(ci_perc_pMTL_plasma[["normal"]][3]),"]",sep = "")#Compared to Ab-PET positive
      results_print$perc_plasmaNeoT_plasma[i_method2]=paste(round(100*v_NpNeoT[[1]]/v_NMTL[[1]],0),"[",round(ci_perc_pNeoT_plasma[["normal"]][2]),", ",round(ci_perc_pNeoT_plasma[["normal"]][3]),"]",sep = "")#Compared to Ab-PET positive
      
      
      
      results_plot=data.frame(matrix(nrow = 5+2,ncol = 4))
      colnames(results_plot)=c("Model",
                               "Perc","Perc_l","Perc_h")
      
      # results_plot$Model[1]="Original"
      results_plot$Model[1]="Plasma"
      results_plot$Model[2]="Ab-PET_plasma"#Mod
      results_plot$Model[2+1]="Ab-PET"
      # results_plot$Model[3]="NeoT"
      results_plot$Model[2+2]="p-MTL"
      results_plot$Model[3+2]="p-NeoT"
      results_plot$Model[4+2]="p-MTL_plasma"
      results_plot$Model[5+2]="p-NeoT_plasma"
      
      # results_plot$Perc[1]=100
      # results_plot$Perc_l[1]=NA
      # results_plot$Perc_h[1]=NA
      
      results_plot$Perc[1]=100*(v_Nplasma[[1]])/v_N0[[1]]
      results_plot$Perc_l[1]=ci_perc_plasma0[["normal"]][2]
      results_plot$Perc_h[1]=ci_perc_plasma0[["normal"]][3]

      #Mod      
      results_plot$Perc[2]=100*(v_NMTL[[1]])/v_Nplasma[[1]]
      results_plot$Perc_l[2]=ci_perc_Ab_plasma[["normal"]][2]#Compared to plasma+
      results_plot$Perc_h[2]=ci_perc_Ab_plasma[["normal"]][3]#Compared to plasma+

      results_plot$Perc[2+1]=100*(v_NMTL[[1]])/v_N0[[1]]
      results_plot$Perc_l[2+1]=ci_perc_Ab0[["normal"]][2]#Compared to all
      results_plot$Perc_h[2+1]=ci_perc_Ab0[["normal"]][3]#Compared to all
      # 
      # results_plot$Perc[3]=100*(v_NNeoT[[1]])/v_N0[[1]]
      # results_plot$Perc_l[3]=ci_perc_NeoT0[["normal"]][2]
      # results_plot$Perc_h[3]=ci_perc_NeoT0[["normal"]][3]
      
      results_plot$Perc[2+2]=100*(v_NpMTL[[1]])/v_N0[[1]]
      results_plot$Perc_l[2+2]=ci_perc_pMTL0[["normal"]][2]
      results_plot$Perc_h[2+2]=ci_perc_pMTL0[["normal"]][3]
      
      results_plot$Perc[3+2]=100*(v_NpNeoT[[1]])/v_N0[[1]]
      results_plot$Perc_l[3+2]=ci_perc_pNeoT0[["normal"]][2]
      results_plot$Perc_h[3+2]=ci_perc_pNeoT0[["normal"]][3]
  
      #Mod
      results_plot$Perc[4+2]=100*(v_NpMTL[[1]])/v_NMTL[[1]]
      results_plot$Perc_l[4+2]=ci_perc_pMTL_plasma[["normal"]][2]
      results_plot$Perc_h[4+2]=ci_perc_pMTL_plasma[["normal"]][3]
      
      #Mod
      results_plot$Perc[5+2]=100*(v_NpNeoT[[1]])/v_NMTL[[1]]
      results_plot$Perc_l[5+2]=ci_perc_pNeoT_plasma[["normal"]][2]
      results_plot$Perc_h[5+2]=ci_perc_pNeoT_plasma[["normal"]][3]
      
      #Mod:
      results_line=data.frame(matrix(nrow = 8,ncol = 5))
      colnames(results_line)=c("Model","step",
                               "N","N_l","N_h")
      
      results_line$Model=c(rep("MTL",4),rep("NeoT",4))
      results_line$step=c(0,1,2,3,0,1,2,3)
      results_line$N=c(100,100*(v_Nplasma[[1]])/v_N0[[1]],100*(v_NMTL[[1]])/v_N0[[1]],100*(v_NpMTL[[1]])/v_N0[[1]],
                       100,100*(v_Nplasma[[1]])/v_N0[[1]],100*(v_NMTL[[1]])/v_N0[[1]],100*(v_NpNeoT[[1]])/v_N0[[1]])
      results_line$N_l=c(NA,ci_perc_plasma0[["normal"]][2],ci_perc_Ab0[["normal"]][2],ci_perc_pMTL0[["normal"]][2],
                         NA,ci_perc_plasma0[["normal"]][2],ci_perc_Ab0[["normal"]][2],ci_perc_pNeoT0[["normal"]][2])
      results_line$N_h=c(NA,ci_perc_plasma0[["normal"]][3],ci_perc_Ab0[["normal"]][3],ci_perc_pMTL0[["normal"]][3],
                         NA,ci_perc_plasma0[["normal"]][3],ci_perc_Ab0[["normal"]][3],ci_perc_pNeoT0[["normal"]][3])
      results_line$Model=factor(results_line$Model,levels=c("MTL","NeoT"))
      results_line$step=as.factor(results_line$step)
      
      # p6=ggplot(results_plot, aes(x = factor(Model,levels = c("Plasma","p-MTL","p-NeoT","p-MTL_plasma","p-NeoT_plasma")), y = Perc,fill=Model)) +
      #   geom_bar(stat="identity", position=position_dodge()) +
      #   geom_errorbar(aes(ymin=Perc_l, ymax=Perc_h), width=.2,
      #                 position=position_dodge(.9))+
      #   geom_signif(comparisons = list(#c("Plasma","MTL"),c("MTL","NeoT"),c("Plasma","NeoT"),
      #                                  c("p-MTL","p-NeoT"),#c("MTL","p-MTL"),
      #                                  c("Plasma","p-MTL"),
      #                                  #c("NeoT","p-NeoT"),
      #                                  c("Plasma","p-NeoT")),annotations = c_sign,step_increase = 0.1,y_position = max(results_plot$Perc))+
      #   scale_fill_manual(values=c_colors)+ylab("Final sample (%)") + xlab("Model")+
      #   # geom_text(aes(label = paste(round(Perc,0),"%",sep = ""),x=Model,y=Perc_h+25), stat="identity", vjust=0.1,size=4)+
      #   theme_classic()+theme(legend.position = "none",axis.text=element_text(size=8),
      #                         axis.title=element_text(size=10,face="bold"),legend.text = element_text(size=6),
      #                         axis.line.x = element_line(size = 0.3),axis.line.y = element_line(size = 0.3),
      #                         axis.ticks.x=element_line(size = 0.3),axis.ticks.y = element_line(size = 0.3),
      #                         legend.key.size = unit(0.2, 'cm'),legend.title = element_blank())+
      #   geom_hline(yintercept=100,linetype="dotted")#+
      # 
      #   # scale_y_continuous(limits=c(0,100))
      # 
      # 
      # ggexport(p6,filename=paste(dir_fig,"/SampleRed_convMCI_power_",v_power0,"_red_",v_redPerc,"_",c_method[i_method],"_Ab_",l_ab,"_adj_",l_adj,"_20240116.pdf",sep = ""))
  
      
      writexl::write_xlsx(results_plot,path = paste(dir_res,"/Plot_sampleRed_convMCI_power_",v_power0,"_red_",v_redPerc,"_Plasma_",c_method_plasma[i_method],"_PET_",c_method_PET[i_method2],"_Ab_",l_ab,"_adj_",l_adj,"_20241009.xlsx",sep = ""))
      
    }
    
    writexl::write_xlsx(results_print,path = paste(dir_res,"/SampleRed_convMCI_power_",v_power0,"_red_",v_redPerc,"_Plasma_",c_method_plasma[i_method],"_PET_Ab_",l_ab,"_adj_",l_adj,"_20241009.xlsx",sep = ""))
    
  }
  
