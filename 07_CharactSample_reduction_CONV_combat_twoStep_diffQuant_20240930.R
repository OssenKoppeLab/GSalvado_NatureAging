dir_res="~/Desktop/Lund/20_TauPET_plasma/Results_Merged/20240917_Rev/Two_step_sample_red_conv/"
dir_fig="~/Desktop/Lund/20_TauPET_plasma/Figures_Merged/20240917_Rev/Two_step_sample_red_conv/"
load("~/Desktop/Lund/20_TauPET_plasma/Data_all/Data_all_cogn_comb_20240917.RData")

dataAll$time_MCI[which(dataAll$cohort=="WU" | dataAll$cohort=="WRAP")]=dataAll$time_MCI[which(dataAll$cohort=="WU" | dataAll$cohort=="WRAP")]*365.25
dataAll=dataAll %>% filter(!is.na(time_MCI))
dataAll=dataAll %>% filter(!duplicated(sid))

c_out_name=c("convMCI")#,"MMSE")
l_PET="NeoT"#"MTL"#
c_PET=c("z_NeoT","z_MTL")
n_out=length(c_out_name)
l_ab="Abpos"#"all"#
l_DX="allEtiol"
l_adj="raw"#"combat"
l_out="_noCovs"#"noOutlier"
c_method_plasma=c("Q2_Q4","Q3_Q4","Q4")#"zscore2"
c_method_PET=c("Q2_Q4","Q3_Q4","Q4")#"zscore2"

v_cost_plasma=300
v_cost_PET=4000

i_out=1

results_print_MTL=data.frame(matrix(nrow = 12,ncol = 16))
colnames(results_print_MTL)=c("Group_plasma","Group_PET","N_neg","N_pos","Ab_neg","Ab_pos","conv_neg","conv_pos",
                              "APOE_neg","APOE_pos","Fem_neg","Fem_pos","Age_neg","Age_pos","Edu_neg","Edu_pos")

results_plot_MTL=data.frame(matrix(nrow = 2*12,ncol = 12))
colnames(results_plot_MTL)=c("Group_plasmaPET","Group","N","Ab",
                             "APOE","Fem","Conv","Age","Age_SD","Edu","Edu_SD","USD")

i_count=1
i_count2=1


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
    
    df=dataAll
    
    df=df[complete.cases(df[,"conv_MCI"]),]
    
    if (l_ab=="Abpos") {
      df=df %>% filter(ab_pos==1)
    }
    
    df$plasma_pos=as.factor(ifelse(df[,"z_ptau217"]>quantile(df[!duplicated(df$sid),"z_ptau217"],v_perc),1,0))#as.factor(ifelse(dataAll$z_ptau217_comb>1.96,1,0))
    
   
    #We select the quartiles of those participants already plasma positive to assess the additive reduction of performing a PET after being plasma+
    if (l_PET=="MTL") {
      df$MTL_pos=as.factor(ifelse(df[,c_PET[2]]>quantile(df[,c_PET[2]],v_perc_PET),1,0))
    }else{
      df$MTL_pos=as.factor(ifelse(df[,c_PET[1]]>quantile(df[,c_PET[1]],v_perc_PET),1,0))
    }
    
    df3=df %>% filter(!duplicated(sid))

    if (i_method2==1) {
      results_print_MTL$Group_plasma[i_count]=c_method_plasma[i_method]
      results_print_MTL$Group_PET[i_count]="all"
      results_print_MTL$N_neg[i_count]=sum((df3$plasma_pos==0)==T)
      results_print_MTL$N_pos[i_count]=sum((df3$plasma_pos==1)==T)
      results_print_MTL$Ab_neg[i_count]=paste(sum((df3[which(df3$plasma_pos==0),"ab_pos"]==1)==T,na.rm = T)," (",format(round(100*sum((df3[which(df3$plasma_pos==0),"ab_pos"]==1)==T,na.rm = T)/sum((df3$plasma_pos==0)==T),1),nsmall=1),"%)",sep = "")
      results_print_MTL$Ab_pos[i_count]=paste(sum((df3[which(df3$plasma_pos==1),"ab_pos"]==1)==T,na.rm = T)," (",format(round(100*sum((df3[which(df3$plasma_pos==1),"ab_pos"]==1)==T,na.rm = T)/sum((df3$plasma_pos==1)==T),1),nsmall=1),"%)",sep = "")
      results_print_MTL$Fem_neg[i_count]=paste(sum((df3[which(df3$plasma_pos==0),"sex"]==1)==T,na.rm = T)," (",format(round(100*sum((df3[which(df3$plasma_pos==0),"sex"]==1)==T,na.rm = T)/sum((df3$plasma_pos==0)==T),1),nsmall=1),"%)",sep = "")
      results_print_MTL$Fem_pos[i_count]=paste(sum((df3[which(df3$plasma_pos==1),"sex"]==1)==T,na.rm = T)," (",format(round(100*sum((df3[which(df3$plasma_pos==1),"sex"]==1)==T,na.rm = T)/sum((df3$plasma_pos==1)==T),1),nsmall=1),"%)",sep = "")
      results_print_MTL$APOE_neg[i_count]=paste(sum((df3[which(df3$plasma_pos==0),"APOE_e4"]==1)==T,na.rm = T)," (",format(round(100*sum((df3[which(df3$plasma_pos==0),"APOE_e4"]==1)==T,na.rm = T)/sum((df3$plasma_pos==0)==T),1),nsmall=1),"%)",sep = "")
      results_print_MTL$APOE_pos[i_count]=paste(sum((df3[which(df3$plasma_pos==1),"APOE_e4"]==1)==T,na.rm = T)," (",format(round(100*sum((df3[which(df3$plasma_pos==1),"APOE_e4"]==1)==T,na.rm = T)/sum((df3$plasma_pos==1)==T),1),nsmall=1),"%)",sep = "")
      results_print_MTL$conv_neg[i_count]=paste(sum((df3[which(df3$plasma_pos==0),"conv_MCI"]==1)==T,na.rm = T)," (",format(round(100*sum((df3[which(df3$plasma_pos==0),"conv_MCI"]==1)==T,na.rm = T)/sum((df3$plasma_pos==0)==T),1),nsmall=1),"%)",sep = "")
      results_print_MTL$conv_pos[i_count]=paste(sum((df3[which(df3$plasma_pos==1),"conv_MCI"]==1)==T,na.rm = T)," (",format(round(100*sum((df3[which(df3$plasma_pos==1),"conv_MCI"]==1)==T,na.rm = T)/sum((df3$plasma_pos==1)==T),1),nsmall=1),"%)",sep = "")
      results_print_MTL$Age_neg[i_count]=paste(format(round(mean(df3[which(df3$plasma_pos==0),"age"],na.rm = T),1),nsmall=1)," (",format(round(sd(df3[which(df3$plasma_pos==0),"age"],na.rm = T),1),nsmall=1),")",sep = "")
      results_print_MTL$Age_pos[i_count]=paste(format(round(mean(df3[which(df3$plasma_pos==1),"age"],na.rm = T),1),nsmall=1)," (",format(round(sd(df3[which(df3$plasma_pos==1),"age"],na.rm = T),1),nsmall=1),")",sep = "")
      results_print_MTL$Edu_neg[i_count]=paste(format(round(mean(df3[which(df3$plasma_pos==0),"edu"],na.rm = T),1),nsmall=1)," (",format(round(sd(df3[which(df3$plasma_pos==0),"edu"],na.rm = T),1),nsmall=1),")",sep = "")
      results_print_MTL$Edu_pos[i_count]=paste(format(round(mean(df3[which(df3$plasma_pos==1),"edu"],na.rm = T),1),nsmall=1)," (",format(round(sd(df3[which(df3$plasma_pos==1),"edu"],na.rm = T),1),nsmall=1),")",sep = "")
      
      i_count=i_count+1
      
      results_plot_MTL$Group_plasmaPET[i_count2]=paste(c_method_plasma[i_method],"_all",sep="")
      results_plot_MTL$Group[i_count2]="Neg"
      results_plot_MTL$N[i_count2]=sum(df3$plasma_pos==0)
      results_plot_MTL$Ab[i_count2]=sum((df3[which(df3$plasma_pos==0),"ab_pos"]==1)==T,na.rm = T)
      results_plot_MTL$APOE[i_count2]=sum((df3[which(df3$plasma_pos==0),"APOE_e4"]==1)==T,na.rm = T)
      results_plot_MTL$Fem[i_count2]=sum((df3[which(df3$plasma_pos==0),"sex"]==1)==T,na.rm = T)
      results_plot_MTL$Conv[i_count2]=sum((df3[which(df3$plasma_pos==0),"conv_MCI"]==1)==T,na.rm = T)
      results_plot_MTL$Age[i_count2]=mean(df3$age[which(df3$plasma_pos==0)])
      results_plot_MTL$Age_SD[i_count2]=sd(df3$age[which(df3$plasma_pos==0)])
      results_plot_MTL$Edu[i_count2]=mean(df3$edu[which(df3$plasma_pos==0)])
      results_plot_MTL$Edu_SD[i_count2]=sd(df3$edu[which(df3$plasma_pos==0)])
      results_plot_MTL$USD[i_count2]=nrow(df3)*v_cost_plasma#+1492*(1-v_perc_PET)*(1-v_perc)*v_cost_PET ###We only count the cost of plasma
        
      i_count2=i_count2+1
      
      results_plot_MTL$Group_plasmaPET[i_count2]=paste(c_method_plasma[i_method],"_all",sep="")
      results_plot_MTL$Group[i_count2]="Pos"
      results_plot_MTL$N[i_count2]=sum(df3$plasma_pos==1)
      results_plot_MTL$Ab[i_count2]=sum((df3[which(df3$plasma_pos==1),"ab_pos"]==1)==T,na.rm = T)
      results_plot_MTL$APOE[i_count2]=sum((df3[which(df3$plasma_pos==1),"APOE_e4"]==1)==T,na.rm = T)
      results_plot_MTL$Fem[i_count2]=sum((df3[which(df3$plasma_pos==1),"sex"]==1)==T,na.rm = T)
      results_plot_MTL$Conv[i_count2]=sum((df3[which(df3$plasma_pos==1),"conv_MCI"]==1)==T,na.rm = T)
      results_plot_MTL$Age[i_count2]=mean(df3$age[which(df3$plasma_pos==1)])
      results_plot_MTL$Age_SD[i_count2]=sd(df3$age[which(df3$plasma_pos==1)])
      results_plot_MTL$Edu[i_count2]=mean(df3$edu[which(df3$plasma_pos==1)])
      results_plot_MTL$Edu_SD[i_count2]=sd(df3$edu[which(df3$plasma_pos==1)])
      results_plot_MTL$USD[i_count2]=nrow(df3)*v_cost_plasma#+1492*(1-v_perc_PET)*(1-v_perc)*v_cost_PET ###We only count the cost of plasma
      
      i_count2=i_count2+1
      
    }
    
    results_print_MTL$Group_plasma[i_count]=c_method_plasma[i_method]
    results_print_MTL$Group_PET[i_count]=c_method_PET[i_method2]
    results_print_MTL$N_neg[i_count]=sum((df3$plasma_pos==0 | df3$MTL_pos==0)==T)
    results_print_MTL$N_pos[i_count]=sum((df3$plasma_pos==1 & df3$MTL_pos==1)==T)
    results_print_MTL$Ab_neg[i_count]=paste(sum((df3[which(df3$plasma_pos==0 | df3$MTL_pos==0),"ab_pos"]==1)==T,na.rm = T)," (",format(round(100*sum((df3[which(df3$plasma_pos==0 | df3$MTL_pos==0),"ab_pos"]==1)==T,na.rm = T)/sum((df3$plasma_pos==0 | df3$MTL_pos==0)==T),1),nsmall=1),"%)",sep = "")
    results_print_MTL$Ab_pos[i_count]=paste(sum((df3[which(df3$plasma_pos==1 & df3$MTL_pos==1),"ab_pos"]==1)==T,na.rm = T)," (",format(round(100*sum((df3[which(df3$plasma_pos==1 & df3$MTL_pos==1),"ab_pos"]==1)==T,na.rm = T)/sum((df3$plasma_pos==1 & df3$MTL_pos==1)==T),1),nsmall=1),"%)",sep = "")
    results_print_MTL$Fem_neg[i_count]=paste(sum((df3[which(df3$plasma_pos==0 | df3$MTL_pos==0),"sex"]==1)==T,na.rm = T)," (",format(round(100*sum((df3[which(df3$plasma_pos==0 | df3$MTL_pos==0),"sex"]==1)==T,na.rm = T)/sum((df3$plasma_pos==0 | df3$MTL_pos==0)==T),1),nsmall=1),"%)",sep = "")
    results_print_MTL$Fem_pos[i_count]=paste(sum((df3[which(df3$plasma_pos==1 & df3$MTL_pos==1),"sex"]==1)==T,na.rm = T)," (",format(round(100*sum((df3[which(df3$plasma_pos==1 & df3$MTL_pos==1),"sex"]==1)==T,na.rm = T)/sum((df3$plasma_pos==1 & df3$MTL_pos==1)==T),1),nsmall=1),"%)",sep = "")
    results_print_MTL$APOE_neg[i_count]=paste(sum((df3[which(df3$plasma_pos==0 | df3$MTL_pos==0),"APOE_e4"]==1)==T,na.rm = T)," (",format(round(100*sum((df3[which(df3$plasma_pos==0 | df3$MTL_pos==0),"APOE_e4"]==1)==T,na.rm = T)/sum((df3$plasma_pos==0 | df3$MTL_pos==0)==T),1),nsmall=1),"%)",sep = "")
    results_print_MTL$APOE_pos[i_count]=paste(sum((df3[which(df3$plasma_pos==1 & df3$MTL_pos==1),"APOE_e4"]==1)==T,na.rm = T)," (",format(round(100*sum((df3[which(df3$plasma_pos==1 & df3$MTL_pos==1),"APOE_e4"]==1)==T,na.rm = T)/sum((df3$plasma_pos==1 & df3$MTL_pos==1)==T),1),nsmall=1),"%)",sep = "")
    results_print_MTL$conv_neg[i_count]=paste(sum((df3[which(df3$plasma_pos==0 | df3$MTL_pos==0),"conv_MCI"]==1)==T,na.rm = T)," (",format(round(100*sum((df3[which(df3$plasma_pos==0 | df3$MTL_pos==0),"conv_MCI"]==1)==T,na.rm = T)/sum((df3$plasma_pos==0 | df3$MTL_pos==0)==T),1),nsmall=1),"%)",sep = "")
    results_print_MTL$conv_pos[i_count]=paste(sum((df3[which(df3$plasma_pos==1 & df3$MTL_pos==1),"conv_MCI"]==1)==T,na.rm = T)," (",format(round(100*sum((df3[which(df3$plasma_pos==1 & df3$MTL_pos==1),"conv_MCI"]==1)==T,na.rm = T)/sum((df3$plasma_pos==1 & df3$MTL_pos==1)==T),1),nsmall=1),"%)",sep = "")
    results_print_MTL$Age_neg[i_count]=paste(format(round(mean(df3[which(df3$plasma_pos==0 | df3$MTL_pos==0),"age"],na.rm = T),1),nsmall=1)," (",format(round(sd(df3[which(df3$plasma_pos==0 | df3$MTL_pos==0),"age"],na.rm = T),1),nsmall=1),")",sep = "")
    results_print_MTL$Age_pos[i_count]=paste(format(round(mean(df3[which(df3$plasma_pos==1 & df3$MTL_pos==1),"age"],na.rm = T),1),nsmall=1)," (",format(round(sd(df3[which(df3$plasma_pos==1 & df3$MTL_pos==1),"age"],na.rm = T),1),nsmall=1),")",sep = "")
    results_print_MTL$Edu_neg[i_count]=paste(format(round(mean(df3[which(df3$plasma_pos==0 | df3$MTL_pos==0),"edu"],na.rm = T),1),nsmall=1)," (",format(round(sd(df3[which(df3$plasma_pos==0 | df3$MTL_pos==0),"edu"],na.rm = T),1),nsmall=1),")",sep = "")
    results_print_MTL$Edu_pos[i_count]=paste(format(round(mean(df3[which(df3$plasma_pos==1 & df3$MTL_pos==1),"edu"],na.rm = T),1),nsmall=1)," (",format(round(sd(df3[which(df3$plasma_pos==1 & df3$MTL_pos==1),"edu"],na.rm = T),1),nsmall=1),")",sep = "")
    
    i_count=i_count+1
    
    results_plot_MTL$Group_plasmaPET[i_count2]=paste(c_method_plasma[i_method],"_",c_method_PET[i_method2],sep="")
    results_plot_MTL$Group[i_count2]="Neg"
    results_plot_MTL$N[i_count2]=sum(df3$plasma_pos==0 | df3$MTL_pos==0)
    results_plot_MTL$Ab[i_count2]=sum((df3[which(df3$plasma_pos==0 | df3$MTL_pos==0),"ab_pos"]==1)==T,na.rm = T)
    results_plot_MTL$APOE[i_count2]=sum((df3[which(df3$plasma_pos==0 | df3$MTL_pos==0),"APOE_e4"]==1)==T,na.rm = T)
    results_plot_MTL$Fem[i_count2]=sum((df3[which(df3$plasma_pos==0 | df3$MTL_pos==0),"sex"]==1)==T,na.rm = T)
    results_plot_MTL$Conv[i_count2]=sum((df3[which(df3$plasma_pos==0 | df3$MTL_pos==0),"conv_MCI"]==1)==T,na.rm = T)
    results_plot_MTL$Age[i_count2]=mean(df3$age[which(df3$plasma_pos==0 | df3$MTL_pos==0)])
    results_plot_MTL$Age_SD[i_count2]=sd(df3$age[which(df3$plasma_pos==0 | df3$MTL_pos==0)])
    results_plot_MTL$Edu[i_count2]=mean(df3$edu[which(df3$plasma_pos==0 | df3$MTL_pos==0)])
    results_plot_MTL$Edu_SD[i_count2]=sd(df3$edu[which(df3$plasma_pos==0 | df3$MTL_pos==0)])
    results_plot_MTL$USD[i_count2]=nrow(df3)*v_cost_plasma+nrow(df3)*(1-v_perc_PET)*(1-v_perc)*v_cost_PET
    
    i_count2=i_count2+1
    
    results_plot_MTL$Group_plasmaPET[i_count2]=paste(c_method_plasma[i_method],"_",c_method_PET[i_method2],sep="")
    results_plot_MTL$Group[i_count2]="Pos"
    results_plot_MTL$N[i_count2]=sum(df3$plasma_pos==1 & df3$MTL_pos==1)
    results_plot_MTL$Ab[i_count2]=sum((df3[which(df3$plasma_pos==1 & df3$MTL_pos==1),"ab_pos"]==1)==T,na.rm = T)
    results_plot_MTL$APOE[i_count2]=sum((df3[which(df3$plasma_pos==1 & df3$MTL_pos==1),"APOE_e4"]==1)==T,na.rm = T)
    results_plot_MTL$Fem[i_count2]=sum((df3[which(df3$plasma_pos==1 & df3$MTL_pos==1),"sex"]==1)==T,na.rm = T)
    results_plot_MTL$Conv[i_count2]=sum((df3[which(df3$plasma_pos==1 & df3$MTL_pos==1),"conv_MCI"]==1)==T,na.rm = T)
    results_plot_MTL$Age[i_count2]=mean(df3$age[which(df3$plasma_pos==1 & df3$MTL_pos==1)])
    results_plot_MTL$Age_SD[i_count2]=sd(df3$age[which(df3$plasma_pos==1 & df3$MTL_pos==1)])
    results_plot_MTL$Edu[i_count2]=mean(df3$edu[which(df3$plasma_pos==1 & df3$MTL_pos==1)])
    results_plot_MTL$Edu_SD[i_count2]=sd(df3$edu[which(df3$plasma_pos==1 & df3$MTL_pos==1)])
    results_plot_MTL$USD[i_count2]=nrow(df3)*v_cost_plasma+nrow(df3)*(1-v_perc_PET)*(1-v_perc)*v_cost_PET
    
    i_count2=i_count2+1
    
    
  }
}

# results_plot_MTL_pos=results_plot_MTL[which(results_plot_MTL$Group=="Pos"),]


writexl::write_xlsx(results_print_MTL,path = paste(dir_res,"/Charact_sample_red_convMCI_DX_",l_DX,"_Ab_",l_ab,"_",l_out,"_",l_PET,".xlsx",sep = ""))
writexl::write_xlsx(results_plot_MTL,path = paste(dir_res,"/Plot_charact_sample_red_convMCI_DX_",l_DX,"_Ab_",l_ab,"_",l_out,"_",l_PET,".xlsx",sep = ""))

# results_plot_MTL=readxl::read_xlsx(paste(dir_res,"/Plot_charact_sample_red_convMCI_DX_",l_DX,"_Ab_",l_ab,"_",l_out,"_",l_PET,".xlsx",sep = ""))

results_plot_MTL2=results_plot_MTL[which(results_plot_MTL$Group_plasmaPET=="Q2_Q4_Q2_Q4" |
                                           results_plot_MTL$Group_plasmaPET=="Q3_Q4_Q3_Q4" |
                                           results_plot_MTL$Group_plasmaPET=="Q4_Q4" | 
                                           results_plot_MTL$Group_plasmaPET=="Q4_all"),]

# p0=ggplot(data=results_plot_MTL2,aes(x = Group_plasmaPET, y = 100*N/nrow(df),fill=Group)) +
#   geom_bar(stat="identity",position="dodge")+theme_classic() +#scale_fill_manual(values=c())+
#   # geom_errorbar(aes(ymin=l_Pos, ymax=h_Pos), width=.2,
#   #               position=position_dodge(.9))+
#   ylab("N (%)") + xlab("Sample reduction groups")+
#   scale_fill_manual(values = c("#d3d3d3","#33a02c"))+
#   geom_text(aes(label=paste(round(100*N/nrow(df),1),"%",sep = "")),position = "dodge",size=6)+
#   theme_classic()+theme(legend.position = "none",axis.text=element_text(size=18),
#                         axis.title=element_text(size=24,face="bold"),legend.text = element_text(size=6),
#                         axis.line.x = element_line(size = 0.3),axis.line.y = element_line(size = 0.3),
#                         axis.ticks.x=element_line(size = 0.3),axis.ticks.y = element_line(size = 0.3),
#                         legend.key.size = unit(0.2, 'cm'),legend.title = element_blank())+
#   scale_y_continuous(limits=c(0,100))
# 
# ggexport(p0,filename=paste(dir_fig,"/Barplots_charact_sample_red_conv_Npos_DX_",l_DX,"_Ab_",l_ab,"_",l_out,".pdf",sep = "")) #Do not change
# 
# p1=ggplot(data=results_plot_MTL2,aes(x = Group_plasmaPET, y = 100*Ab/N,fill=Group)) +
#   geom_bar(stat="identity",position="dodge")+theme_classic() +#scale_fill_manual(values=c())+
#   # geom_errorbar(aes(ymin=l_Pos, ymax=h_Pos), width=.2,
#   #               position=position_dodge(.9))+
#   ylab("Ab+ (%)") + xlab("Sample reduction groups")+
#   scale_fill_manual(values = c("#d3d3d3","#33a02c"))+
#   geom_text(aes(label=paste(round(100*Ab/N,1),"%",sep = "")),position = "dodge",size=6)+
#   theme_classic()+theme(legend.position = "none",axis.text=element_text(size=18),
#                         axis.title=element_text(size=24,face="bold"),legend.text = element_text(size=6),
#                         axis.line.x = element_line(size = 0.3),axis.line.y = element_line(size = 0.3),
#                         axis.ticks.x=element_line(size = 0.3),axis.ticks.y = element_line(size = 0.3),
#                         legend.key.size = unit(0.2, 'cm'),legend.title = element_blank())+
#   scale_y_continuous(limits=c(0,100))
# 
# ggexport(p1,filename=paste(dir_fig,"/Barplots_charact_sample_red_conv_Abpos_DX_",l_DX,"_Ab_",l_ab,"_",l_out,".pdf",sep = "")) #Do not change
# 
# p2=ggplot(data=results_plot_MTL2,aes(x = Group_plasmaPET, y = 100*APOE/N,fill=Group)) +
#   geom_bar(stat="identity",position="dodge")+theme_classic() +#scale_fill_manual(values=c())+
#   # geom_errorbar(aes(ymin=l_Pos, ymax=h_Pos), width=.2,
#   #               position=position_dodge(.9))+
#   ylab("APOE+ (%)") + xlab("Sample reduction groups")+
#   geom_text(aes(label=paste(round(100*APOE/N,1),"%",sep = "")),position = "dodge",size=12)+
#   theme_classic()+theme(legend.position = "none",axis.text=element_text(size=18),
#                         axis.title=element_text(size=24,face="bold"),legend.text = element_text(size=6),
#                         axis.line.x = element_line(size = 0.3),axis.line.y = element_line(size = 0.3),
#                         axis.ticks.x=element_line(size = 0.3),axis.ticks.y = element_line(size = 0.3),
#                         legend.key.size = unit(0.2, 'cm'),legend.title = element_blank())+
#   scale_y_continuous(limits=c(0,100))
# 
# ggexport(p2,filename=paste(dir_fig,"/Barplots_charact_sample_red_conv_APOE_DX_",l_DX,"_Ab_",l_ab,"_",l_out,".pdf",sep = "")) #Do not change
# 
# p3=ggplot(data=results_plot_MTL2,aes(x = Group_plasmaPET, y = 100*Fem/N,fill=Group)) +
#   geom_bar(stat="identity",position="dodge")+theme_classic() +#scale_fill_manual(values=c())+
#   # geom_errorbar(aes(ymin=l_Pos, ymax=h_Pos), width=.2,
#   #               position=position_dodge(.9))+
#   ylab("Women (%)") + xlab("Sample reduction groups")+
#   geom_text(aes(label=paste(round(100*Fem/N,1),"%",sep = "")),position = "dodge",size=12)+
#   theme_classic()+theme(legend.position = "none",axis.text=element_text(size=18),
#                         axis.title=element_text(size=24,face="bold"),legend.text = element_text(size=6),
#                         axis.line.x = element_line(size = 0.3),axis.line.y = element_line(size = 0.3),
#                         axis.ticks.x=element_line(size = 0.3),axis.ticks.y = element_line(size = 0.3),
#                         legend.key.size = unit(0.2, 'cm'),legend.title = element_blank())+
#   scale_y_continuous(limits=c(0,100))
# 
# ggexport(p3,filename=paste(dir_fig,"/Barplots_charact_sample_red_conv_sex_DX_",l_DX,"_Ab_",l_ab,"_",l_out,".pdf",sep = "")) #Do not change
# 
# 
# p4=ggplot(data=results_plot_MTL2,aes(x = Group_plasmaPET, y = 100*Conv/N,fill=Group)) +
#   geom_bar(stat="identity",position="dodge")+theme_classic() +#scale_fill_manual(values=c())+
#   # geom_errorbar(aes(ymin=l_Pos, ymax=h_Pos), width=.2,
#   #               position=position_dodge(.9))+
#   ylab("Progressed MCI (%)") + xlab("Sample reduction groups")+
#   scale_fill_manual(values = c("#d3d3d3","#33a02c"))+
#   geom_text(aes(label=paste(round(100*Conv/N,1),"%",sep = "")),position = "dodge",size=6)+
#   theme_classic()+theme(legend.position = "none",axis.text=element_text(size=18),
#                         axis.title=element_text(size=24,face="bold"),legend.text = element_text(size=6),
#                         axis.line.x = element_line(size = 0.3),axis.line.y = element_line(size = 0.3),
#                         axis.ticks.x=element_line(size = 0.3),axis.ticks.y = element_line(size = 0.3),
#                         legend.key.size = unit(0.2, 'cm'),legend.title = element_blank())+
#   scale_y_continuous(limits=c(0,100))
# 
# ggexport(p4,filename=paste(dir_fig,"/Barplots_charact_sample_red_conv_ConvMCI_DX_",l_DX,"_Ab_",l_ab,"_",l_out,".pdf",sep = "")) #Do not change
# 
# pcost=ggplot(data=results_plot_MTL2[which(results_plot_MTL2$Group=="Pos"),],aes(x = Group_plasmaPET, y = USD/(1000*Conv),fill=Group)) +
#   geom_bar(stat="identity",position="dodge")+theme_classic() +#scale_fill_manual(values=c())+
#   # geom_errorbar(aes(ymin=l_Pos, ymax=h_Pos), width=.2,
#   #               position=position_dodge(.9))+
#   ylab("Cost (kUSD/conv)") + xlab("Sample reduction groups")+
#   scale_fill_manual(values = c("#33a02c"))+
#   geom_text(aes(label=paste(round(USD/(1000*Conv),1),sep = "")),position = "dodge",size=6)+
#   geom_hline(yintercept = nrow(df)*v_cost_PET/(1000*sum(df3$conv_MCI)),linetype="dashed")+#172 number of conversions in the full sample
#   theme_classic()+theme(legend.position = "none",axis.text=element_text(size=18),
#                         axis.title=element_text(size=24,face="bold"),legend.text = element_text(size=6),
#                         axis.line.x = element_line(size = 0.3),axis.line.y = element_line(size = 0.3),
#                         axis.ticks.x=element_line(size = 0.3),axis.ticks.y = element_line(size = 0.3),
#                         legend.key.size = unit(0.2, 'cm'),legend.title = element_blank())
# 
# 
# ggexport(pcost,filename=paste(dir_fig,"/Barplots_charact_sample_red_conv_Cost_DX_",l_DX,"_Ab_",l_ab,"_",l_out,".pdf",sep = "")) #Do not change
# 
# pl2=ggarrange(p4,p1,p0,pcost,#plines_part_CSF,
#               # labels = c("A", "B", "C","D","E"),
#               ncol = 2, nrow = 2)
# ggexport(pl2,filename=paste(dir_fig,"/Barplots_charact_sample_red_conv_summary_",c_out[i_out],"_DX_",l_DX,"_Ab_",l_ab,"_",l_out,".pdf",sep = "")) #Do not change
# 
# 
# p5=ggplot(data=results_plot_MTL2,aes(x = Group_plasmaPET, y = Age,fill=Group)) +
#   geom_bar(stat="identity",position="dodge")+theme_classic() +#scale_fill_manual(values=c())+
#   geom_errorbar(aes(ymin=Age-Age_SD/sqrt(N), ymax=Age+Age_SD/sqrt(N)), width=.2,
#                 position=position_dodge(.9))+
#   ylab("Age") + xlab("Sample reduction groups")+
#   #  geom_text(aes(label=paste(round(100*APOE/N,1),"%",sep = "")),position = "dodge",size=12)+
#   theme_classic()+theme(legend.position = "none",axis.text=element_text(size=18),
#                         axis.title=element_text(size=24,face="bold"),legend.text = element_text(size=6),
#                         axis.line.x = element_line(size = 0.3),axis.line.y = element_line(size = 0.3),
#                         axis.ticks.x=element_line(size = 0.3),axis.ticks.y = element_line(size = 0.3),
#                         legend.key.size = unit(0.2, 'cm'),legend.title = element_blank())
# # scale_y_continuous(limits=c(0,100))
# 
# ggexport(p5,filename=paste(dir_fig,"/Barplots_charact_sample_red_conv_Age_DX_",l_DX,"_Ab_",l_ab,"_",l_out,".pdf",sep = "")) #Do not change
# 
# p6=ggplot(data=results_plot_MTL2,aes(x = Group_plasmaPET, y = Edu,fill=Group)) +
#   geom_bar(stat="identity",position="dodge")+theme_classic() +#scale_fill_manual(values=c())+
#   geom_errorbar(aes(ymin=Edu-Edu_SD/sqrt(N), ymax=Edu+Edu_SD/sqrt(N)), width=.2,
#                 position=position_dodge(.9))+
#   ylab("Education") + xlab("Sample reduction groups")+
#   #  geom_text(aes(label=paste(round(100*APOE/N,1),"%",sep = "")),position = "dodge",size=12)+
#   theme_classic()+theme(legend.position = "none",axis.text=element_text(size=18),
#                         axis.title=element_text(size=24,face="bold"),legend.text = element_text(size=6),
#                         axis.line.x = element_line(size = 0.3),axis.line.y = element_line(size = 0.3),
#                         axis.ticks.x=element_line(size = 0.3),axis.ticks.y = element_line(size = 0.3),
#                         legend.key.size = unit(0.2, 'cm'),legend.title = element_blank())
# # scale_y_continuous(limits=c(0,100))
# 
# ggexport(p6,filename=paste(dir_fig,"/Barplots_charact_sample_red_conv_edu_DX_",l_DX,"_Ab_",l_ab,"_",l_out,".pdf",sep = "")) #Do not change
# 

#####
if (l_PET=="MTL") {
  v_color="#00a087ff"
}else{v_color="#1f78b4"}

p0=ggplot(data=results_plot_MTL2,aes(x = Group_plasmaPET, y = 100*N/nrow(df),fill=Group)) +
  geom_bar(stat="identity",position="dodge")+theme_classic() +#scale_fill_manual(values=c())+
  # geom_errorbar(aes(ymin=l_Pos, ymax=h_Pos), width=.2,
  #               position=position_dodge(.9))+
  ylab("N (%)") + xlab("Sample reduction groups")+
  scale_fill_manual(values = c("#d3d3d3",v_color))+
  geom_text(aes(label=paste(round(100*N/nrow(df),0),"%",sep = "")),position = position_dodge(.9),size=3)+
  theme_classic()+theme(legend.position = "none",axis.text=element_text(size=8),#axis.text.x = element_text(angle=45),
                        axis.title=element_text(size=10,face="bold"),legend.text = element_text(size=6),
                        axis.line.x = element_line(size = 0.1),axis.line.y = element_line(size = 0.1),
                        axis.ticks.x=element_line(size = 0.1),axis.ticks.y = element_line(size = 0.1),
                        legend.key.size = unit(0.2, 'cm'),legend.title = element_blank())+scale_y_continuous(limits=c(0,100))

p1=ggplot(data=results_plot_MTL2,aes(x = Group_plasmaPET, y = 100*Ab/N,fill=Group)) +
  geom_bar(stat="identity",position="dodge")+theme_classic() +#scale_fill_manual(values=c())+
  # geom_errorbar(aes(ymin=l_Pos, ymax=h_Pos), width=.2,
  #               position=position_dodge(.9))+
  ylab("Ab+ (%)") + xlab("Sample reduction groups")+
  scale_fill_manual(values = c("#d3d3d3",v_color))+
  geom_text(aes(label=paste(round(100*Ab/N,0),"%",sep = "")),position = position_dodge(.9),size=3)+
  theme_classic()+theme(legend.position = "none",axis.text=element_text(size=8),#axis.text.x = element_text(angle=45),
                        axis.title=element_text(size=10,face="bold"),legend.text = element_text(size=6),
                        axis.line.x = element_line(size = 0.1),axis.line.y = element_line(size = 0.1),
                        axis.ticks.x=element_line(size = 0.1),axis.ticks.y = element_line(size = 0.1),
                        legend.key.size = unit(0.2, 'cm'),legend.title = element_blank())+scale_y_continuous(limits=c(0,100))




p4=ggplot(data=results_plot_MTL2,aes(x = Group_plasmaPET, y = 100*Conv/N,fill=Group)) +
  geom_bar(stat="identity",position="dodge")+theme_classic() +#scale_fill_manual(values=c())+
  # geom_errorbar(aes(ymin=l_Pos, ymax=h_Pos), width=.2,
  #               position=position_dodge(.9))+
  ylab("Progressed MCI (%)") + xlab("Sample reduction groups")+
  scale_fill_manual(values = c("#d3d3d3",v_color))+
  geom_text(aes(label=paste(round(100*Conv/N,0),"%",sep = "")),position = position_dodge(.9),size=3)+
  theme_classic()+theme(legend.position = "none",axis.text=element_text(size=8),#axis.text.x = element_text(angle=45),
                        axis.title=element_text(size=10,face="bold"),legend.text = element_text(size=6),
                        axis.line.x = element_line(size = 0.1),axis.line.y = element_line(size = 0.1),
                        axis.ticks.x=element_line(size = 0.1),axis.ticks.y = element_line(size = 0.1),
                        legend.key.size = unit(0.2, 'cm'),legend.title = element_blank())+scale_y_continuous(limits=c(0,100))


pcost=ggplot(data=results_plot_MTL2[which(results_plot_MTL2$Group=="Pos"),],aes(x = Group_plasmaPET, y = USD/(1000*Conv),fill=Group)) +
  geom_bar(stat="identity",position="dodge")+theme_classic() +#scale_fill_manual(values=c())+
  # geom_errorbar(aes(ymin=l_Pos, ymax=h_Pos), width=.2,
  #               position=position_dodge(.9))+
  ylab("Cost (kUSD/conv)") + xlab("Sample reduction groups")+
  scale_fill_manual(values = c(v_color))+
  geom_text(aes(label=paste(round(USD/(1000*Conv),1),sep = "")),position = position_dodge(.9),size=3)+
  geom_hline(yintercept = nrow(df)*v_cost_PET/(1000*sum(df$conv_MCI)),linetype="dashed")+#172 number of conversions in the full sample
  theme_classic()+theme(legend.position = "none",axis.text=element_text(size=8),#axis.text.x = element_text(angle=45),
                        axis.title=element_text(size=10,face="bold"),legend.text = element_text(size=6),
                        axis.line.x = element_line(size = 0.1),axis.line.y = element_line(size = 0.1),
                        axis.ticks.x=element_line(size = 0.1),axis.ticks.y = element_line(size = 0.1),
                        legend.key.size = unit(0.2, 'cm'),legend.title = element_blank())



pl2=ggarrange(p4,p1,p0,pcost,#plines_part_CSF,
              # labels = c("A", "B", "C","D","E"),
              ncol = 3, nrow = 3)
c_out=c_out_name
ggexport(pl2,filename=paste(dir_fig,"/Barplots_charact_sample_red_conv_summary_",c_out[i_out],"_DX_",l_DX,"_Ab_",l_ab,"_",l_out,"_",l_PET,".pdf",sep = "")) #Do not change
