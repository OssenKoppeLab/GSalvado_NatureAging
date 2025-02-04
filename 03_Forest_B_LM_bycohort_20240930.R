####
dir_res="~/Desktop/Lund/20_TauPET_plasma/Results_Merged/20240917_Rev/Cogn_decline_LM"
dir_fig="~/Desktop/Lund/20_TauPET_plasma/Figures_Merged/20240917_Rev/Cogn_decline_LM"
####
l_ab="all"#"Abpos"#
l_DX="allEtiol"
l_out="__cov_Int_time_APOE_cohort"#"_cov_Int_time_APOE"#"noOutlier"
l_adj="raw"

c_out=c("mpacc","mmse_score")
list_cohort=c("dataBF2","dataBF1","dataAIBL","dataTRIAD","dataPREVENT","dataWRAP","dataAms","dataKA","dataMCSA")#,"dataGen")
name_cohort=c("B2","B1","Ai","T","PA","W","Am","KA","MCSA")#,"Gen")

for (i_out in 1) {#:length(c_out)
  
  if (i_out==1) {
    v_lim=c(-0.065,0.018)
    x_pos=NA# 0.1
  }else{
    v_lim=c(-0.25,0.27)
    x_pos=NA#0.2
  }
  
  # dataBF1=readxl::read_xlsx(paste(dir_res,"/Results_plot_LME_BF2_AIBL_TRIAD_PREVENT_",c_out[i_out],"_DX_allEtiol_Ab_all",l_out,".xlsx",sep = ""))
  # dataBF2=readxl::read_xlsx(paste(dir_res,"/Results_plot_LME_BF1_AIBL_TRIAD_PREVENT_",c_out[i_out],"_DX_allEtiol_Ab_all",l_out,".xlsx",sep = ""))
  # dataAIBL=readxl::read_xlsx(paste(dir_res,"/Results_plot_LME_BF2_BF1_TRIAD_PREVENT_",c_out[i_out],"_DX_allEtiol_Ab_all",l_out,".xlsx",sep = ""))
  # dataTRIAD=readxl::read_xlsx(paste(dir_res,"/Results_plot_LME_BF2_BF1_AIBL_PREVENT_",c_out[i_out],"_DX_allEtiol_Ab_all",l_out,".xlsx",sep = ""))
  # dataPREVENT=readxl::read_xlsx(paste(dir_res,"/Results_plot_LME_BF2_BF1_AIBL_TRIAD_",c_out[i_out],"_DX_allEtiol_Ab_all",l_out,".xlsx",sep = ""))

  dataBF1=readxl::read_xlsx(paste(dir_res,"/Results_B_LM_BF1_",c_out[i_out],"_DX_allEtiol_Ab_",l_ab,"_adj_",l_adj,l_out,"_Q4.xlsx",sep = ""))
  dataBF2=readxl::read_xlsx(paste(dir_res,"/Results_B_LM_BF2_",c_out[i_out],"_DX_allEtiol_Ab_",l_ab,"_adj_",l_adj,l_out,"_Q4.xlsx",sep = ""))
  dataAIBL=readxl::read_xlsx(paste(dir_res,"/Results_B_LM_AIBL_",c_out[i_out],"_DX_allEtiol_Ab_",l_ab,"_adj_",l_adj,l_out,"_Q4.xlsx",sep = ""))
  dataTRIAD=readxl::read_xlsx(paste(dir_res,"/Results_B_LM_TRIAD_",c_out[i_out],"_DX_allEtiol_Ab_",l_ab,"_adj_",l_adj,l_out,"_Q4.xlsx",sep = ""))
  dataPREVENT=readxl::read_xlsx(paste(dir_res,"/Results_B_LM_PREVENT_",c_out[i_out],"_DX_allEtiol_Ab_",l_ab,"_adj_",l_adj,l_out,"_Q4.xlsx",sep = ""))
  dataWRAP=readxl::read_xlsx(paste(dir_res,"/Results_B_LM_WRAP_",c_out[i_out],"_DX_allEtiol_Ab_",l_ab,"_adj_",l_adj,l_out,"_Q4.xlsx",sep = ""))
  dataAms=readxl::read_xlsx(paste(dir_res,"/Results_B_LM_Ams_",c_out[i_out],"_DX_allEtiol_Ab_",l_ab,"_adj_",l_adj,l_out,"_Q4.xlsx",sep = ""))
  dataKA=readxl::read_xlsx(paste(dir_res,"/Results_B_LM_WU_",c_out[i_out],"_DX_allEtiol_Ab_",l_ab,"_adj_",l_adj,l_out,"_Q4.xlsx",sep = ""))
  dataMCSA=readxl::read_xlsx(paste(dir_res,"/Results_B_LM_MCSA_",c_out[i_out],"_DX_allEtiol_Ab_",l_ab,"_adj_",l_adj,l_out,"_Q4.xlsx",sep = ""))
  # dataGen=readxl::read_xlsx(paste(dir_res,"/Results_B_LM_Gen_",c_out[i_out],"_DX_allEtiol_Ab_all_adj_",l_adj,l_out,".xlsx",sep = ""))
  
  dataAll=readxl::read_xlsx(paste(dir_res,"/Results_B_LM_allCohorts_",c_out[i_out],"_DX_allEtiol_Ab_",l_ab,"_adj_",l_adj,l_out,"_Q4.xlsx",sep = ""))
  
  table_plasma=data.frame(matrix(nrow = length(name_cohort)*3,ncol = 6))
  colnames(table_plasma)=c("N","B","B_l","B_h","Model","Cohort_excl")#data.frame(beta=NA,se=NA,group=NA,pval=NA,mdl=NA)
  
  i_count=1
  for (i_cohort in 1:length(name_cohort)) {
    df=get(list_cohort[i_cohort])
    table_plasma$N[i_count]=df$N[1]
    table_plasma$B[i_count]=df$B[which(df$name=="plasma p-tau217")]
    table_plasma$B_l[i_count]=df$B_l[which(df$name=="plasma p-tau217")]
    table_plasma$B_h[i_count]=df$B_h[which(df$name=="plasma p-tau217")]
    table_plasma$Model[i_count]="Single"
    table_plasma$Cohort_excl[i_count]=name_cohort[i_cohort]
    
    table_plasma$N[i_count+1]=df$N[1]
    table_plasma$B[i_count+1]=df$B[which(df$name=="Plasma+ERC/Amygd")]
    table_plasma$B_l[i_count+1]=df$B_l[which(df$name=="Plasma+ERC/Amygd")]
    table_plasma$B_h[i_count+1]=df$B_h[which(df$name=="Plasma+ERC/Amygd")]
    table_plasma$Model[i_count+1]="With MTL"
    table_plasma$Cohort_excl[i_count+1]=name_cohort[i_cohort]

    table_plasma$N[i_count+2]=df$N[1]
    table_plasma$B[i_count+2]=df$B[which(df$name=="Plasma+NeoT")]
    table_plasma$B_l[i_count+2]=df$B_l[which(df$name=="Plasma+NeoT")]
    table_plasma$B_h[i_count+2]=df$B_h[which(df$name=="Plasma+NeoT")]
    table_plasma$Model[i_count+2]="With NeoT"
    table_plasma$Cohort_excl[i_count+2]=name_cohort[i_cohort]
    
    i_count=i_count+3
  }
  
  # table_plasma$Cohort_excl=factor(table_plasma$Cohort_excl,levels=rev(c("AIBL","BF1","BF2","PREVENT-AD","TRIAD")))
  
  plot_plasma=table_plasma %>% filter(Model=="Single") %>% #group_by(Model) %>%
    mutate(index = reorder(Cohort_excl, N)) %>%
    ggplot( aes(y=index, x=B, xmin=B_l, xmax=B_h)) +
    geom_rect(aes(ymin = -Inf, xmin = dataAll$B_h[which(dataAll$name=="plasma p-tau217")],
                  ymax = Inf, xmax = dataAll$B_l[which(dataAll$name=="plasma p-tau217")]),alpha=0.2,
              fill = "lightgrey")+
    geom_point(aes(size=N),color="#e64b35ff",shape=23,fill="#e64b35ff") +
    geom_errorbarh(height=.3,color="black",size=0.5)+
    labs( x=paste("Beta"), y = ' ') +xlim(v_lim)+
    geom_vline(xintercept = dataAll$B[which(dataAll$name=="plasma p-tau217")],linetype="dashed")+
    geom_vline(xintercept = 0,linetype="dotted")+
    geom_text(aes(label=paste("N=",N,sep = ""),y=Cohort_excl),x=x_pos,hjust=0,size=4)+
    theme_classic()+ggtitle(paste(c_out[i_out]))+
    theme(axis.text=element_text(size=14),axis.title = element_text(size=14),
          legend.position = "none",plot.title = element_text(hjust = 0.5),
          strip.text.y = element_text(size = 10))+
    scale_color_manual(values=c(rep("black",length(name_cohort))))+
    scale_fill_manual(values=c(rep("#e64b35ff",length(name_cohort))))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.ticks = element_line(size=0.1),axis.line = element_line(size=0.1))#+scale_x_continuous(v_lim)
  #   scale_x_reverse(limits=v_lim)
  
  plot_plasma_MTL=table_plasma %>% filter(Model=="With MTL") %>% #group_by(Model) %>%
    mutate(index = reorder(Cohort_excl, N)) %>%
    ggplot( aes(y=index, x=B, xmin=B_l, xmax=B_h)) +
    geom_rect(aes(ymin = -Inf, xmin = dataAll$B_h[which(dataAll$name=="Plasma+ERC/Amygd" & dataAll$Biomarker=="Plasma")],
                  ymax = Inf, xmax = dataAll$B_l[which(dataAll$name=="Plasma+ERC/Amygd" & dataAll$Biomarker=="Plasma")]),alpha=0.2,
              fill = "lightgrey")+
    geom_point(aes(size=N),color="#f39b7fff",shape=23,fill="#f39b7fff") +
    geom_errorbarh(height=.3,color="black",size=0.5)+
    labs( x=paste("Beta"), y = ' ') +xlim(v_lim)+
    geom_vline(xintercept = dataAll$B[which(dataAll$name=="Plasma+ERC/Amygd" & dataAll$Biomarker=="Plasma")],linetype="dashed")+
    geom_vline(xintercept = 0,linetype="dotted")+
    geom_text(aes(label=paste("N=",N,sep = ""),y=Cohort_excl),x=x_pos,hjust=0,size=4)+
    theme_classic()+ggtitle(paste(c_out[i_out]))+
    theme(axis.text=element_text(size=14),axis.title = element_text(size=14),
          legend.position = "none",plot.title = element_text(hjust = 0.5),
          strip.text.y = element_text(size = 10))+
    scale_color_manual(values=c(rep("black",length(name_cohort))))+
    scale_fill_manual(values=c(rep("#f39b7fff",length(name_cohort))))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.ticks = element_line(size=0.1),axis.line = element_line(size=0.1))#+scale_x_continuous(v_lim)
  # +
  #   scale_x_reverse(limits=v_lim)
  
  plot_plasma_NeoT=table_plasma %>% filter(Model=="With NeoT") %>% #group_by(Model) %>%
    mutate(index = reorder(Cohort_excl, N)) %>%
    ggplot( aes(y=index, x=B, xmin=B_l, xmax=B_h)) +
    geom_rect(aes(ymin = -Inf, xmin = dataAll$B_h[which(dataAll$name=="Plasma+NeoT" & dataAll$Biomarker=="Plasma")],
                  ymax = Inf, xmax = dataAll$B_l[which(dataAll$name=="Plasma+NeoT" & dataAll$Biomarker=="Plasma")]),alpha=0.2,
              fill = "lightgrey")+
    geom_point(aes(size=N),color="#f39b7fff",shape=23,fill="#f39b7fff") +
    geom_errorbarh(height=.3,aes(color=Cohort_excl),size=0.5)+
    labs( x=paste("Beta"), y = ' ') +xlim(v_lim)+
    geom_vline(xintercept = dataAll$B[which(dataAll$name=="Plasma+NeoT" & dataAll$Biomarker=="Plasma")],linetype="dashed")+
    geom_vline(xintercept = 0,linetype="dotted")+
    geom_text(aes(label=paste("N=",N,sep = ""),y=Cohort_excl),x=x_pos,hjust=0,size=4)+
    theme_classic()+ggtitle(paste(c_out[i_out]))+
    theme(axis.text=element_text(size=14),axis.title = element_text(size=14),
          legend.position = "none",plot.title = element_text(hjust = 0.5),
          strip.text.y = element_text(size = 10))+
    scale_color_manual(values=c(rep("black",length(name_cohort))))+
    scale_fill_manual(values=c(rep("#f39b7fff",length(name_cohort))))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.ticks = element_line(size=0.1),axis.line = element_line(size=0.1))#+scale_x_continuous(v_lim)
  # +
  #   scale_x_reverse(limits=v_lim)
  
  table_PET=data.frame(matrix(nrow = length(name_cohort)*2,ncol = 6))
  colnames(table_PET)=c("N","B","B_l","B_h","Model","Cohort_excl")#data.frame(beta=NA,se=NA,group=NA,pval=NA,mdl=NA)
  
  i_count=1
  for (i_cohort in 1:length(name_cohort)) {
    df=get(list_cohort[i_cohort])
    table_PET$N[i_count]=df$N[1]
    table_PET$B[i_count]=df$B[which(df$Biomarker=="z_NeoT" & df$name=="NeoT")]
    table_PET$B_l[i_count]=df$B_l[which(df$Biomarker=="z_NeoT" & df$name=="NeoT")]
    table_PET$B_h[i_count]=df$B_h[which(df$Biomarker=="z_NeoT" & df$name=="NeoT")]
    table_PET$Model[i_count]="Single"
    table_PET$Cohort_excl[i_count]=name_cohort[i_cohort]
    
    table_PET$N[i_count+1]=df$N[1]
    table_PET$B[i_count+1]=df$B[which(df$Biomarker=="z_NeoT" & df$name=="Plasma+NeoT")]
    table_PET$B_l[i_count+1]=df$B_l[which(df$Biomarker=="z_NeoT" & df$name=="Plasma+NeoT")]
    table_PET$B_h[i_count+1]=df$B_h[which(df$Biomarker=="z_NeoT" & df$name=="Plasma+NeoT")]
    table_PET$Model[i_count+1]="With plasma"
    table_PET$Cohort_excl[i_count+1]=name_cohort[i_cohort]
    
    i_count=i_count+2
  }
  

  plot_NeoT=table_PET %>% filter(Model=="Single") %>% #group_by(Model) %>%
    mutate(index = reorder(Cohort_excl, N)) %>%
    ggplot( aes(y=index, x=B, xmin=B_l, xmax=B_h)) +
    geom_rect(aes(ymin = -Inf, xmin = dataAll$B_h[which(dataAll$name=="NeoT")],
                  ymax = Inf, xmax = dataAll$B_l[which(dataAll$name=="NeoT")]),alpha=0.2,
              fill = "lightgrey")+
    geom_point(aes(size=N),color="#3c5488ff",shape=23,fill="#3c5488ff") +
    geom_errorbarh(height=.3,color="black",size=0.5)+
    labs( x=paste("Beta"), y = ' ') +xlim(v_lim)+
    geom_vline(xintercept = dataAll$B[which(dataAll$name=="NeoT")],linetype="dashed")+
    geom_vline(xintercept = 0,linetype="dotted")+
    geom_text(aes(label=paste("N=",N,sep = ""),y=Cohort_excl),x=x_pos,hjust=0,size=4)+
    theme_classic()+ggtitle(paste(c_out[i_out]))+
    theme(axis.text=element_text(size=14),axis.title = element_text(size=14),
          legend.position = "none",plot.title = element_text(hjust = 0.5),
          strip.text.y = element_text(size = 10))+
    scale_color_manual(values=c(rep("black",length(name_cohort))))+
    scale_fill_manual(values=c(rep("#3c5488ff",length(name_cohort))))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.ticks = element_line(size=0.1),axis.line = element_line(size=0.1))#+scale_x_continuous(v_lim)
  # +
  #   scale_x_reverse(limits=v_lim)
  
  plot_NeoT2=table_PET %>% filter(Model=="With plasma") %>% #group_by(Model) %>%
    mutate(index = reorder(Cohort_excl, N)) %>%
    ggplot( aes(y=index, x=B, xmin=B_l, xmax=B_h)) +
    geom_rect(aes(ymin = -Inf, xmin = dataAll$B_h[which(dataAll$name=="Plasma+NeoT" & dataAll$Biomarker=="z_NeoT")],
                  ymax = Inf, xmax = dataAll$B_l[which(dataAll$name=="Plasma+NeoT" & dataAll$Biomarker=="z_NeoT")]),alpha=0.2,
              fill = "lightgrey")+
    geom_point(aes(size=N),color="#8491b4",shape=23,fill="#8491b4") +
    geom_errorbarh(height=.3,aes(color=Cohort_excl),size=0.5)+
    labs( x=paste("Beta"), y = ' ') +xlim(v_lim)+
    geom_vline(xintercept = dataAll$B[which(dataAll$name=="Plasma+NeoT" & dataAll$Biomarker=="z_NeoT")],linetype="dashed")+
    geom_vline(xintercept = 0,linetype="dotted")+
    geom_text(aes(label=paste("N=",N,sep = ""),y=Cohort_excl),x=x_pos,hjust=0,size=4)+
    theme_classic()+ggtitle(paste(c_out[i_out]))+
    theme(axis.text=element_text(size=14),axis.title = element_text(size=14),
          legend.position = "none",plot.title = element_text(hjust = 0.5),
          strip.text.y = element_text(size = 10))+
    scale_color_manual(values=c(rep("black",length(name_cohort))))+
    scale_fill_manual(values=c(rep("#8491b4",length(name_cohort))))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.ticks = element_line(size=0.1),axis.line = element_line(size=0.1))#+scale_x_continuous(v_lim)
  # +
  #   scale_x_reverse(limits=v_lim)
  
  ##MTL
  i_count=1
  for (i_cohort in 1:length(name_cohort)) {
    df=get(list_cohort[i_cohort])
    table_PET$N[i_count]=df$N[1]
    table_PET$B[i_count]=df$B[which(df$Biomarker=="z_MTL" & df$name=="ERC/Amygd")]
    table_PET$B_l[i_count]=df$B_l[which(df$Biomarker=="z_MTL" & df$name=="ERC/Amygd")]
    table_PET$B_h[i_count]=df$B_h[which(df$Biomarker=="z_MTL" & df$name=="ERC/Amygd")]
    table_PET$Model[i_count]="Single"
    table_PET$Cohort_excl[i_count]=name_cohort[i_cohort]
    
    table_PET$N[i_count+1]=df$N[1]
    table_PET$B[i_count+1]=df$B[which(df$Biomarker=="z_MTL" & df$name=="Plasma+ERC/Amygd")]
    table_PET$B_l[i_count+1]=df$B_l[which(df$Biomarker=="z_MTL" & df$name=="Plasma+ERC/Amygd")]
    table_PET$B_h[i_count+1]=df$B_h[which(df$Biomarker=="z_MTL" & df$name=="Plasma+ERC/Amygd")]
    table_PET$Model[i_count+1]="With plasma"
    table_PET$Cohort_excl[i_count+1]=name_cohort[i_cohort]
    
    i_count=i_count+2
  }
  

  plot_MTL=table_PET %>% filter(Model=="Single") %>% #group_by(Model) %>%
    mutate(index = reorder(Cohort_excl, N)) %>%
    ggplot( aes(y=index, x=B, xmin=B_l, xmax=B_h)) +
    geom_rect(aes(ymin = -Inf, xmin = dataAll$B_h[which(dataAll$name=="ERC/Amygd")],
                  ymax = Inf, xmax = dataAll$B_l[which(dataAll$name=="ERC/Amygd")]),alpha=0.2,
              fill = "lightgrey")+
    geom_point(aes(size=N),color="#00a087ff",shape=23,fill="#00a087ff") +
    geom_errorbarh(height=.3,color="black",size=0.5)+
    labs( x=paste("Beta"), y = ' ') +xlim(v_lim)+
    geom_vline(xintercept = dataAll$B[which(dataAll$name=="ERC/Amygd")],linetype="dashed")+
    geom_vline(xintercept = 0,linetype="dotted")+
    geom_text(aes(label=paste("N=",N,sep = ""),y=Cohort_excl),x=x_pos,hjust=0,size=4)+
    theme_classic()+ggtitle(paste(c_out[i_out]))+
    theme(axis.text=element_text(size=14),axis.title = element_text(size=14),
          legend.position = "none",plot.title = element_text(hjust = 0.5),
          strip.text.y = element_text(size = 10))+
    scale_color_manual(values=c(rep("black",length(name_cohort))))+
    scale_fill_manual(values=c(rep("#00a087ff",length(name_cohort))))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.ticks = element_line(size=0.1),axis.line = element_line(size=0.1))#+scale_x_continuous(v_lim)
  # +
  #   scale_x_reverse(limits=v_lim)
  
    plot_MTL2=table_PET %>% filter(Model=="With plasma") %>% #group_by(Model) %>%
      mutate(index = reorder(Cohort_excl, N)) %>%
      ggplot( aes(y=index, x=B, xmin=B_l, xmax=B_h)) +
    geom_rect(aes(ymin = -Inf, xmin = dataAll$B_h[which(dataAll$name=="Plasma+ERC/Amygd" & dataAll$Biomarker=="z_MTL")],
                  ymax = Inf, xmax = dataAll$B_l[which(dataAll$name=="Plasma+ERC/Amygd" & dataAll$Biomarker=="z_MTL")]),alpha=0.2,
              fill = "lightgrey")+
    geom_point(aes(size=N),color="#91d1c2ff",shape=23,fill="#91d1c2ff") +
    geom_errorbarh(height=.3,color="black",size=0.5)+
    labs( x=paste("Beta"), y = ' ') +xlim(v_lim)+
    geom_text(aes(label=paste("N=",N,sep = ""),y=Cohort_excl),x=x_pos,hjust=0,size=4)+
    geom_vline(xintercept = dataAll$B[which(dataAll$name=="Plasma+ERC/Amygd" & dataAll$Biomarker=="z_MTL")],linetype="dashed")+
    geom_vline(xintercept = 0,linetype="dotted")+
      theme_classic()+ggtitle(paste(c_out[i_out]))+
    theme(axis.text=element_text(size=14),axis.title = element_text(size=14),
          legend.position = "none",plot.title = element_text(hjust = 0.5),
          strip.text.y = element_text(size = 10))+
    scale_color_manual(values=c(rep("black",length(name_cohort))))+
    scale_fill_manual(values=c(rep("#91d1c2ff",length(name_cohort))))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.ticks = element_line(size=0.1),axis.line = element_line(size=0.1))#+scale_x_continuous(v_lim)
    # +
    # scale_x_reverse(limits=v_lim)
  
  
    pf2=ggarrange(plot_plasma+rremove("y.text"),  plot_MTL +rremove("y.text"), plot_NeoT+rremove("y.text"),
                plot_plasma_MTL+rremove("y.text"), plot_MTL2+rremove("y.text"),"",
                plot_plasma_NeoT+rremove("y.text"),"",plot_NeoT2+rremove("y.text"),
                # labels = c("A", "B", "C","D","E"),
                ncol = 3, nrow = 3)
  
  
  ggexport(pf2,filename=paste(dir_fig,"/Compare_LM_B_",c_out[i_out],"_DX_",l_DX,"_Ab_",l_ab,"_adj_",l_adj,"_allPlots",l_out,"_onlyOneCohort.pdf",sep = ""))
  
  
}



##### R2
for (i_out in 1:length(c_out)) {
  
  if (i_out==1) {
    v_lim=c(-0.5,0.8)
    x_pos=NA#0.22
  }else{
    v_lim=c(-0.5,0.8)
    x_pos=NA#0.55
  }
  
  # dataBF1=readxl::read_xlsx(paste(dir_res,"/Results_plot_LME_BF2_AIBL_TRIAD_PREVENT_",c_out[i_out],"_DX_allEtiol_Ab_all",l_out,".xlsx",sep = ""))
  # dataBF2=readxl::read_xlsx(paste(dir_res,"/Results_plot_LME_BF1_AIBL_TRIAD_PREVENT_",c_out[i_out],"_DX_allEtiol_Ab_all",l_out,".xlsx",sep = ""))
  # dataAIBL=readxl::read_xlsx(paste(dir_res,"/Results_plot_LME_BF2_BF1_TRIAD_PREVENT_",c_out[i_out],"_DX_allEtiol_Ab_all",l_out,".xlsx",sep = ""))
  # dataTRIAD=readxl::read_xlsx(paste(dir_res,"/Results_plot_LME_BF2_BF1_AIBL_PREVENT_",c_out[i_out],"_DX_allEtiol_Ab_all",l_out,".xlsx",sep = ""))
  # dataPREVENT=readxl::read_xlsx(paste(dir_res,"/Results_plot_LME_BF2_BF1_AIBL_TRIAD_",c_out[i_out],"_DX_allEtiol_Ab_all",l_out,".xlsx",sep = ""))
  
  dataBF1=readxl::read_xlsx(paste(dir_res,"/Results_plot_LM_BF1_",c_out[i_out],"_DX_allEtiol_Ab_all_adj_",l_adj,l_out,".xlsx",sep = ""))
  dataBF2=readxl::read_xlsx(paste(dir_res,"/Results_plot_LM_BF2_",c_out[i_out],"_DX_allEtiol_Ab_all_adj_",l_adj,l_out,".xlsx",sep = ""))
  dataAIBL=readxl::read_xlsx(paste(dir_res,"/Results_plot_LM_AIBL_",c_out[i_out],"_DX_allEtiol_Ab_all_adj_",l_adj,l_out,".xlsx",sep = ""))
  dataTRIAD=readxl::read_xlsx(paste(dir_res,"/Results_plot_LM_TRIAD_",c_out[i_out],"_DX_allEtiol_Ab_all_adj_",l_adj,l_out,".xlsx",sep = ""))
  dataPREVENT=readxl::read_xlsx(paste(dir_res,"/Results_plot_LM_PREVENT_",c_out[i_out],"_DX_allEtiol_Ab_all_adj_",l_adj,l_out,".xlsx",sep = ""))
  dataWRAP=readxl::read_xlsx(paste(dir_res,"/Results_plot_LM_WRAP_",c_out[i_out],"_DX_allEtiol_Ab_all_adj_",l_adj,l_out,".xlsx",sep = ""))
  dataAms=readxl::read_xlsx(paste(dir_res,"/Results_plot_LM_Ams_",c_out[i_out],"_DX_allEtiol_Ab_all_adj_",l_adj,l_out,".xlsx",sep = ""))
  dataKA=readxl::read_xlsx(paste(dir_res,"/Results_plot_LM_WU_",c_out[i_out],"_DX_allEtiol_Ab_all_adj_",l_adj,l_out,".xlsx",sep = ""))
  dataMCSA=readxl::read_xlsx(paste(dir_res,"/Results_plot_LM_MCSA_",c_out[i_out],"_DX_allEtiol_Ab_all_adj_",l_adj,l_out,".xlsx",sep = ""))
  # dataGen=readxl::read_xlsx(paste(dir_res,"/Results_plot_LM_Gen_",c_out[i_out],"_DX_allEtiol_Ab_all_adj_",l_adj,l_out,".xlsx",sep = ""))
  
  dataAll=readxl::read_xlsx(paste(dir_res,"/Results_plot_LM_allCohorts_",c_out[i_out],"_DX_allEtiol_Ab_all_adj_",l_adj,l_out,"_Q4.xlsx",sep = ""))
  
  table_plasma=data.frame(matrix(nrow = length(name_cohort),ncol = 6))
  colnames(table_plasma)=c("N","R2","R2_l","R2_h","Model","Cohort_excl")#data.frame(beta=NA,se=NA,group=NA,pval=NA,mdl=NA)
  
  i_count=1
  for (i_cohort in 1:length(name_cohort)) {
    df=get(list_cohort[i_cohort])
    table_plasma$N[i_count]=df$N[1]
    table_plasma$R2[i_count]=df$R2[which(df$name=="plasma p-tau217")]
    table_plasma$R2_l[i_count]=df$R2_l[which(df$name=="plasma p-tau217")]
    table_plasma$R2_h[i_count]=df$R2_h[which(df$name=="plasma p-tau217")]
    table_plasma$Model[i_count]="Single"
    table_plasma$Cohort_excl[i_count]=name_cohort[i_cohort]
    
    # table_plasma$N[i_count+1]=df$N[1]
    # table_plasma$R2[i_count+1]=df$R2[which(df$name=="Plasma+ERC/Amygd")]
    # table_plasma$R2_l[i_count+1]=df$R2_l[which(df$name=="Plasma+ERC/Amygd")]
    # table_plasma$R2_h[i_count+1]=df$R2_h[which(df$name=="Plasma+ERC/Amygd")]
    # table_plasma$Model[i_count+1]="With MTL"
    # table_plasma$Cohort_excl[i_count+1]=name_cohort[i_cohort]
    # 
    # table_plasma$N[i_count+2]=df$N[1]
    # table_plasma$R2[i_count+2]=df$R2[which(df$name=="Plasma+NeoT")]
    # table_plasma$R2_l[i_count+2]=df$R2_l[which(df$name=="Plasma+NeoT")]
    # table_plasma$R2_h[i_count+2]=df$R2_h[which(df$name=="Plasma+NeoT")]
    # table_plasma$Model[i_count+2]="With NeoT"
    # table_plasma$Cohort_excl[i_count+2]=name_cohort[i_cohort]
    
    i_count=i_count+1#3
  }
  
  # table_plasma$Cohort_excl=factor(table_plasma$Cohort_excl,levels=rev(c("AIBL","BF1","BF2","PREVENT-AD","TRIAD")))
  
  plot_plasma=table_plasma %>% filter(Model=="Single") %>% #group_by(Model) %>%
    mutate(index = reorder(Cohort_excl, N)) %>%
    ggplot( aes(y=index, x=R2, xmin=R2_l, xmax=R2_h)) +
    geom_rect(aes(ymin = -Inf, xmin = dataAll$R2_h[which(dataAll$name=="plasma p-tau217")],
                  ymax = Inf, xmax = dataAll$R2_l[which(dataAll$name=="plasma p-tau217")]),alpha=0.2,
              fill = "lightgrey")+
    geom_point(aes(size=N),color="#e64b35ff",shape=23,fill="#e64b35ff") +
    geom_errorbarh(height=.3,color="black",size=0.5)+
    labs( x=paste("R2"), y = ' ') +xlim(v_lim)+
    geom_vline(xintercept = dataAll$R2[which(dataAll$name=="plasma p-tau217")],linetype="dashed")+
    geom_vline(xintercept = 0,linetype="dotted")+
    # geom_text(aes(label=paste("N=",N,sep = ""),y=Cohort_excl),x=x_pos,hjust=0,size=4)+
    theme_classic()+ggtitle(paste(c_out[i_out]))+
    theme(axis.text=element_text(size=14),axis.title = element_text(size=14),
          legend.position = "none",plot.title = element_text(hjust = 0.5),
          strip.text.y = element_text(size = 10))+
    # scale_color_manual(values=c(rep("black",length(name_cohort))))+
    scale_fill_manual(values=c(rep("#e64b35ff",length(name_cohort))))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.ticks = element_line(size=0.1),axis.line = element_line(size=0.1))#+scale_x_continuous(v_lim)
  #   scale_x_reverse(limits=v_lim)
  
  # plot_plasma_MTL=table_plasma %>% filter(Model=="With MTL") %>% #group_by(Model) %>%
  #   # mutate(index = reorder(Cohort_excl, -abs(R2))) %>%
  #   ggplot( aes(y=Cohort_excl, x=R2, xmin=R2_l, xmax=R2_h)) +
  #   geom_rect(aes(ymin = -Inf, xmin = dataAll$R2_h[5],
  #                 ymax = Inf, xmax = dataAll$R2_l[5]),alpha=0.2,
  #             fill = "lightgrey")+
  #   geom_point(aes(size=N),color="#f39b7fff",shape=23,fill="#f39b7fff") +
  #   geom_errorbarh(height=.3,aes(color=Cohort_excl),size=0.5)+
  #   labs( x=paste("Beta"), y = ' ') +
  #   geom_vline(xintercept = dataAll$R2[5],linetype="dashed")+
  #   geom_vline(xintercept = 0,linetype="dotted")+
  #   geom_text(aes(label=paste("N=",N,sep = ""),y=Cohort_excl),x=x_pos,hjust=0,size=4)+
  #   theme_classic()+ggtitle(paste(c_out[i_out]))+
  #   theme(axis.text=element_text(size=14),axis.title = element_text(size=14),
  #         legend.position = "none",plot.title = element_text(hjust = 0.5),
  #         strip.text.y = element_text(size = 10))+
  #   scale_color_manual(values=c(rep("black",length(name_cohort))))+
  #   scale_fill_manual(values=c(rep("#f39b7fff",length(name_cohort))))+
  #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+scale_x_continuous(v_lim)
  # # +
  # #   scale_x_reverse(limits=v_lim)
  # 
  # plot_plasma_NeoT=table_plasma %>% filter(Model=="With NeoT") %>% #group_by(Model) %>%
  #   # mutate(index = reorder(Cohort_excl, -abs(R2))) %>%
  #   ggplot( aes(y=Cohort_excl, x=R2, xmin=R2_l, xmax=R2_h)) +
  #   geom_rect(aes(ymin = -Inf, xmin = dataAll$R2_h[4],
  #                 ymax = Inf, xmax = dataAll$R2_l[4]),alpha=0.2,
  #             fill = "lightgrey")+
  #   geom_point(aes(size=N),color="#f39b7fff",shape=23,fill="#f39b7fff") +
  #   geom_errorbarh(height=.3,aes(color=Cohort_excl),size=0.5)+
  #   labs( x=paste("Beta"), y = ' ') +
  #   geom_vline(xintercept = dataAll$R2[4],linetype="dashed")+
  #   geom_vline(xintercept = 0,linetype="dotted")+
  #   geom_text(aes(label=paste("N=",N,sep = ""),y=Cohort_excl),x=x_pos,hjust=0,size=4)+
  #   theme_classic()+ggtitle(paste(c_out[i_out]))+
  #   theme(axis.text=element_text(size=14),axis.title = element_text(size=14),
  #         legend.position = "none",plot.title = element_text(hjust = 0.5),
  #         strip.text.y = element_text(size = 10))+
  #   scale_color_manual(values=c(rep("black",length(name_cohort))))+
  #   scale_fill_manual(values=c(rep("#f39b7fff",length(name_cohort))))+
  #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+scale_x_continuous(v_lim)
  # +
  #   scale_x_reverse(limits=v_lim)
  
  table_PET=data.frame(matrix(nrow = length(name_cohort)*2,ncol = 6))
  colnames(table_PET)=c("N","R2","R2_l","R2_h","Model","Cohort_excl")#data.frame(beta=NA,se=NA,group=NA,pval=NA,mdl=NA)
  
  i_count=1
  for (i_cohort in 1:length(name_cohort)) {
    df=get(list_cohort[i_cohort])
    table_PET$N[i_count]=df$N[1]
    table_PET$R2[i_count]=df$R2[which(df$Biomarker=="z_NeoT" & df$name=="NeoT")]
    table_PET$R2_l[i_count]=df$R2_l[which(df$Biomarker=="z_NeoT" & df$name=="NeoT")]
    table_PET$R2_h[i_count]=df$R2_h[which(df$Biomarker=="z_NeoT" & df$name=="NeoT")]
    table_PET$Model[i_count]="Single"
    table_PET$Cohort_excl[i_count]=name_cohort[i_cohort]
    
    table_PET$N[i_count+1]=df$N[1]
    table_PET$R2[i_count+1]=df$R2[which(df$name=="Plasma+NeoT")]
    table_PET$R2_l[i_count+1]=df$R2_l[which(df$name=="Plasma+NeoT")]
    table_PET$R2_h[i_count+1]=df$R2_h[which(df$name=="Plasma+NeoT")]
    table_PET$Model[i_count+1]="With plasma"
    table_PET$Cohort_excl[i_count+1]=name_cohort[i_cohort]
    
    i_count=i_count+2
  }
  
  
  plot_NeoT=table_PET %>% filter(Model=="Single") %>% #group_by(Model) %>%
    mutate(index = reorder(Cohort_excl, N)) %>%
    ggplot( aes(y=index, x=R2, xmin=R2_l, xmax=R2_h)) +
    geom_rect(aes(ymin = -Inf, xmin = dataAll$R2_h[which(dataAll$name=="NeoT")],
                  ymax = Inf, xmax = dataAll$R2_l[which(dataAll$name=="NeoT")]),alpha=0.2,
              fill = "lightgrey")+
    geom_point(aes(size=N),color="#3c5488ff",shape=23,fill="#3c5488ff") +
    geom_errorbarh(height=.3,aes(color=Cohort_excl),size=0.5)+
    labs( x=paste("Beta"), y = ' ') +
    geom_vline(xintercept = dataAll$R2[which(dataAll$name=="NeoT")],linetype="dashed")+
    geom_vline(xintercept = 0,linetype="dotted")+
    geom_text(aes(label=paste("N=",N,sep = ""),y=Cohort_excl),x=x_pos,hjust=0,size=4)+
    theme_classic()+ggtitle(paste(c_out[i_out]))+
    theme(axis.text=element_text(size=14),axis.title = element_text(size=14),
          legend.position = "none",plot.title = element_text(hjust = 0.5),
          strip.text.y = element_text(size = 10))+
    scale_color_manual(values=c(rep("black",length(name_cohort))))+
    scale_fill_manual(values=c(rep("#3c5488ff",length(name_cohort))))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.ticks = element_line(size=0.1),axis.line = element_line(size=0.1))#+scale_x_continuous(v_lim)
  # +
  #   scale_x_reverse(limits=v_lim)
  
  plot_NeoT2=table_PET %>% filter(Model=="With plasma") %>% #group_by(Model) %>%
    mutate(index = reorder(Cohort_excl, N)) %>%
    ggplot( aes(y=index, x=R2, xmin=R2_l, xmax=R2_h)) +
    geom_rect(aes(ymin = -Inf, xmin = dataAll$R2_h[which(dataAll$name=="Plasma+NeoT")],
                  ymax = Inf, xmax = dataAll$R2_l[which(dataAll$name=="Plasma+NeoT")]),alpha=0.2,
              fill = "lightgrey")+
    geom_point(aes(size=N),color="#8491b4",shape=23,fill="#8491b4") +
    geom_errorbarh(height=.3,aes(color=Cohort_excl),size=0.5)+
    labs( x=paste("Beta"), y = ' ') +
    geom_vline(xintercept = dataAll$R2[which(dataAll$name=="Plasma+NeoT")],linetype="dashed")+
    geom_vline(xintercept = 0,linetype="dotted")+
    geom_text(aes(label=paste("N=",N,sep = ""),y=Cohort_excl),x=x_pos,hjust=0,size=4)+
    theme_classic()+ggtitle(paste(c_out[i_out]))+
    theme(axis.text=element_text(size=14),axis.title = element_text(size=14),
          legend.position = "none",plot.title = element_text(hjust = 0.5),
          strip.text.y = element_text(size = 10))+
    scale_color_manual(values=c(rep("black",length(name_cohort))))+
    scale_fill_manual(values=c(rep("#8491b4",length(name_cohort))))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.ticks = element_line(size=0.1),axis.line = element_line(size=0.1))#+scale_x_continuous(v_lim)
  #   scale_x_reverse(limits=v_lim)
  
  ##MTL
  i_count=1
  for (i_cohort in 1:length(name_cohort)) {
    df=get(list_cohort[i_cohort])
    table_PET$N[i_count]=df$N[1]
    table_PET$R2[i_count]=df$R2[which(df$Biomarker=="z_MTL" & df$name=="ERC/Amygd")]
    table_PET$R2_l[i_count]=df$R2_l[which(df$Biomarker=="z_MTL" & df$name=="ERC/Amygd")]
    table_PET$R2_h[i_count]=df$R2_h[which(df$Biomarker=="z_MTL" & df$name=="ERC/Amygd")]
    table_PET$Model[i_count]="Single"
    table_PET$Cohort_excl[i_count]=name_cohort[i_cohort]
    
    table_PET$N[i_count+1]=df$N[1]
    table_PET$R2[i_count+1]=df$R2[which(df$name=="Plasma+ERC/Amygd")]
    table_PET$R2_l[i_count+1]=df$R2_l[which(df$name=="Plasma+ERC/Amygd")]
    table_PET$R2_h[i_count+1]=df$R2_h[which(df$name=="Plasma+ERC/Amygd")]
    table_PET$Model[i_count+1]="With plasma"
    table_PET$Cohort_excl[i_count+1]=name_cohort[i_cohort]
    
    i_count=i_count+2
  }
  
  
  plot_MTL=table_PET %>% filter(Model=="Single") %>% #group_by(Model) %>%
    mutate(index = reorder(Cohort_excl, N)) %>%
    ggplot( aes(y=index, x=R2, xmin=R2_l, xmax=R2_h)) +
    geom_rect(aes(ymin = -Inf, xmin = dataAll$R2_h[which(dataAll$name=="ERC/Amygd")],
                  ymax = Inf, xmax = dataAll$R2_l[which(dataAll$name=="ERC/Amygd")]),alpha=0.2,
              fill = "lightgrey")+
    geom_point(aes(size=N),color="#00a087ff",shape=23,fill="#00a087ff") +
    geom_errorbarh(height=.3,aes(color=Cohort_excl),size=0.5)+
    labs( x=paste("Beta"), y = ' ') +
    geom_vline(xintercept = dataAll$R2[which(dataAll$name=="ERC/Amygd")],linetype="dashed")+
    geom_vline(xintercept = 0,linetype="dotted")+
    geom_text(aes(label=paste("N=",N,sep = ""),y=Cohort_excl),x=x_pos,hjust=0,size=4)+
    theme_classic()+ggtitle(paste(c_out[i_out]))+
    theme(axis.text=element_text(size=14),axis.title = element_text(size=14),
          legend.position = "none",plot.title = element_text(hjust = 0.5),
          strip.text.y = element_text(size = 10))+
    scale_color_manual(values=c(rep("black",length(name_cohort))))+
    scale_fill_manual(values=c(rep("#00a087ff",length(name_cohort))))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.ticks = element_line(size=0.1),axis.line = element_line(size=0.1))#+scale_x_continuous(v_lim)
  # +
  #   scale_x_reverse(limits=v_lim)
  
  plot_MTL2=table_PET %>% filter(Model=="With plasma") %>% #group_by(Model) %>%
    mutate(index = reorder(Cohort_excl, N)) %>%
    ggplot( aes(y=index, x=R2, xmin=R2_l, xmax=R2_h)) +
    geom_rect(aes(ymin = -Inf, xmin = dataAll$R2_h[which(dataAll$name=="Plasma+ERC/Amygd")],
                  ymax = Inf, xmax = dataAll$R2_l[which(dataAll$name=="Plasma+ERC/Amygd")]),alpha=0.2,
              fill = "lightgrey")+
    geom_point(aes(size=N),color="#91d1c2ff",shape=23,fill="#91d1c2ff") +
    geom_errorbarh(height=.3,aes(color=Cohort_excl),size=0.5)+
    labs( x=paste("Beta"), y = ' ') +
    geom_text(aes(label=paste("N=",N,sep = ""),y=Cohort_excl),x=x_pos,hjust=0,size=4)+
    geom_vline(xintercept = dataAll$R2[which(dataAll$name=="Plasma+ERC/Amygd")],linetype="dashed")+
    geom_vline(xintercept = 0,linetype="dotted")+
    theme_classic()+ggtitle(paste(c_out[i_out]))+
    theme(axis.text=element_text(size=14),axis.title = element_text(size=14),
          legend.position = "none",plot.title = element_text(hjust = 0.5),
          strip.text.y = element_text(size = 10))+
    scale_color_manual(values=c(rep("black",length(name_cohort))))+
    scale_fill_manual(values=c(rep("#91d1c2ff",length(name_cohort))))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.ticks = element_line(size=0.1),axis.line = element_line(size=0.1))#+scale_x_continuous(v_lim)
  # +
  # scale_x_reverse(limits=v_lim)
  
  
  pf2=ggarrange(plot_plasma+rremove("y.text"),  plot_MTL+rremove("y.text") , plot_NeoT+rremove("y.text"),
                #plot_plasma_MTL, 
                "",plot_MTL2+rremove("y.text"),plot_NeoT2+rremove("y.text"),
                #plot_plasma_NeoT,
                
                # labels = c("A", "B", "C","D","E"),
                ncol = 3, nrow = 3)
  
  ggexport(pf2,filename=paste(dir_fig,"/Compare_LM_R2_",c_out[i_out],"_DX_",l_DX,"_Ab_",l_ab,"_adj_",l_adj,"_allPlots",l_out,"_onlyOneCohort.pdf",sep = ""))
  
  
}
