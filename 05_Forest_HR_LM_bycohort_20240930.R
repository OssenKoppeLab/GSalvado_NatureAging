####
dir_res="~/Desktop/Lund/20_TauPET_plasma/Results_Merged/20240917_Rev/Conversion"
dir_fig="~/Desktop/Lund/20_TauPET_plasma/Figures_Merged/20240917_Rev/Conversion"
####
l_ab="Abpos"#"all"#
l_adj="raw"

list_cohort=c("dataBF2","dataBF1","dataAIBL","dataTRIAD","dataPREVENT","dataWRAP","dataAms","dataKA","dataMCSA")
name_cohort=c("B2","B1","Ai","T","PA","W","Am","KA","MCSA")

    v_lim=c(0.2,4.5)#6.7)#c(0.15,-0.15)

  dataBF1=readxl::read_xlsx(paste(dir_res,"/Plot_Kaplan_MCI_BF1_Ab_",l_ab,"_adj_",l_adj,"_SensTime_NA.xlsx",sep = ""))
  dataBF2=readxl::read_xlsx(paste(dir_res,"/Plot_Kaplan_MCI_BF2_Ab_",l_ab,"_adj_",l_adj,"_SensTime_NA.xlsx",sep = ""))
  dataAIBL=readxl::read_xlsx(paste(dir_res,"/Plot_Kaplan_MCI_AIBL_Ab_",l_ab,"_adj_",l_adj,"_SensTime_NA.xlsx",sep = ""))
  dataTRIAD=readxl::read_xlsx(paste(dir_res,"/Plot_Kaplan_MCI_TRIAD_Ab_",l_ab,"_adj_",l_adj,"_SensTime_NA.xlsx",sep = ""))
  dataPREVENT=readxl::read_xlsx(paste(dir_res,"/Plot_Kaplan_MCI_PREVENT_Ab_",l_ab,"_adj_",l_adj,"_SensTime_NA.xlsx",sep = ""))
  dataWRAP=readxl::read_xlsx(paste(dir_res,"/Plot_Kaplan_MCI_WRAP_Ab_",l_ab,"_adj_",l_adj,"_SensTime_NA.xlsx",sep = ""))
  dataAms=readxl::read_xlsx(paste(dir_res,"/Plot_Kaplan_MCI_Ams_Ab_",l_ab,"_adj_",l_adj,"_SensTime_NA.xlsx",sep = ""))
  dataKA=readxl::read_xlsx(paste(dir_res,"/Plot_Kaplan_MCI_WU_Ab_",l_ab,"_adj_",l_adj,"_SensTime_NA.xlsx",sep = ""))
  dataMCSA=readxl::read_xlsx(paste(dir_res,"/Plot_Kaplan_MCI_MCSA_Ab_",l_ab,"_adj_",l_adj,"_SensTime_NA.xlsx",sep = ""))
  
  dataAll=readxl::read_xlsx(paste(dir_res,"/Plot_Kaplan_MCI_allCohorts_Ab_",l_ab,"_adj_",l_adj,"_SensTime_NA.xlsx",sep = ""))
  
  table_plasma=data.frame(matrix(nrow = length(name_cohort)*3,ncol = 6))
  colnames(table_plasma)=c("N","HR","HR_l","HR_h","Model","Cohort_excl")#data.frame(beta=NA,se=NA,group=NA,pval=NA,mdl=NA)
  
  i_count=1
  for (i_cohort in 1:length(name_cohort)) {
    df=get(list_cohort[i_cohort])
    table_plasma$N[i_count]=df$N[1]
    table_plasma$HR[i_count]=df$HR_plasma[1]
    table_plasma$HR_l[i_count]=df$HR_plasma_l[1]
    table_plasma$HR_h[i_count]=df$HR_plasma_h[1]
    table_plasma$Model[i_count]="Single"
    table_plasma$Cohort_excl[i_count]=name_cohort[i_cohort]
    
    table_plasma$N[i_count+1]=df$N[5]
    table_plasma$HR[i_count+1]=df$HR_plasma[5]
    table_plasma$HR_l[i_count+1]=df$HR_plasma_l[5]
    table_plasma$HR_h[i_count+1]=df$HR_plasma_h[5]
    table_plasma$Model[i_count+1]="With MTL"
    table_plasma$Cohort_excl[i_count+1]=name_cohort[i_cohort]

    table_plasma$N[i_count+2]=df$N[4]
    table_plasma$HR[i_count+2]=df$HR_plasma[4]
    table_plasma$HR_l[i_count+2]=df$HR_plasma_l[4]
    table_plasma$HR_h[i_count+2]=df$HR_plasma_h[4]
    table_plasma$Model[i_count+2]="With NeoT"
    table_plasma$Cohort_excl[i_count+2]=name_cohort[i_cohort]
    
    i_count=i_count+3
  }
  
  # table_plasma$Cohort_excl=factor(table_plasma$Cohort_excl,levels=rev(c("AIBL","BF1","BF2","PREVENT-AD","TRIAD")))
  table_plasma$HR_h[which((table_plasma$HR_h)=="Inf")]=99
  table_plasma$HR_h=as.numeric(table_plasma$HR_h)
  
  plot_plasma=table_plasma %>% filter(Model=="Single") %>% #group_by(Model) %>%
    mutate(index = reorder(Cohort_excl, N)) %>%
    ggplot( aes(y=index, x=HR, xmin=HR_l, xmax=HR_h)) +
    geom_rect(aes(ymin = -Inf, xmin = dataAll$HR_plasma_h[1],
                  ymax = Inf, xmax = dataAll$HR_plasma_l[1]),alpha=0.2,
              fill = "lightgrey")+
    geom_point(aes(size=N),shape=23,fill="#e64b35ff",color="#e64b35ff") +
    # facet_grid(index~., scales= "free", space="free") +
    geom_errorbarh(height=.3,aes(color=Cohort_excl),size=0.5)+
    labs( x=paste("HR"), y = ' ') +
    geom_vline(xintercept = dataAll$HR_plasma[1],linetype="dashed")+
    geom_vline(xintercept = 1,linetype="dotted")+
    geom_text(aes(label=paste("N=",N,sep = ""),y=Cohort_excl),x=5,hjust=0,size=4)+
    theme_classic()+ggtitle("MCI")+
    theme(axis.text=element_text(size=14),axis.title = element_text(size=14),
          legend.position = "none",plot.title = element_text(hjust = 0.5),
          strip.text.y = element_text(size = 10))+
    scale_color_manual(values=c(rep("black",length(name_cohort))))+
    scale_fill_manual(values=c(rep("#e64b35ff",length(name_cohort))))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.ticks = element_line(size=0.1),axis.line = element_line(size=0.1))+    
    coord_cartesian(xlim = v_lim)
#+
    # scale_x_reverse(limits=v_lim)
  
  plot_plasma_MTL=table_plasma %>% filter(Model=="With MTL") %>% #group_by(Model) %>%
    mutate(index = reorder(Cohort_excl, N)) %>%
    ggplot( aes(y=index, x=HR, xmin=HR_l, xmax=HR_h)) +
    geom_rect(aes(ymin = -Inf, xmin = dataAll$HR_plasma_h[5],
                  ymax = Inf, xmax = dataAll$HR_plasma_l[5]),alpha=0.2,
              fill = "lightgrey")+
    geom_point(aes(size=N),shape=23,fill="#f39b7fff",color="#f39b7fff") +
    geom_errorbarh(height=.3,aes(color=Cohort_excl),size=0.5)+
    labs( x=paste("HR"), y = ' ') +
    geom_vline(xintercept = dataAll$HR_plasma[5],linetype="dashed")+
    geom_vline(xintercept = 1,linetype="dotted")+
    geom_text(aes(label=paste("N=",N,sep = ""),y=Cohort_excl),x=5,hjust=0,size=4)+
    theme_classic()+ggtitle("MCI")+
    theme(axis.text=element_text(size=14),axis.title = element_text(size=14),
          legend.position = "none",plot.title = element_text(hjust = 0.5),
          strip.text.y = element_text(size = 10))+
    scale_color_manual(values=c(rep("black",length(name_cohort))))+
    scale_fill_manual(values=c(rep("#f39b7fff",length(name_cohort))))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.ticks = element_line(size=0.1),axis.line = element_line(size=0.1))+
    coord_cartesian(xlim = v_lim)
  #+
    # scale_x_reverse(limits=v_lim)
  
  plot_plasma_NeoT=table_plasma %>% filter(Model=="With NeoT") %>% #group_by(Model) %>%
    mutate(index = reorder(Cohort_excl, N)) %>%
    ggplot( aes(y=index, x=HR, xmin=HR_l, xmax=HR_h)) +
    geom_rect(aes(ymin = -Inf, xmin = dataAll$HR_plasma_h[4],
                  ymax = Inf, xmax = dataAll$HR_plasma_l[4]),alpha=0.2,
              fill = "lightgrey")+
    geom_point(aes(size=N),shape=23,fill="#f39b7fff",color="#f39b7fff") +
    geom_errorbarh(height=.3,aes(color=Cohort_excl),size=0.5)+
    labs( x=paste("HR"), y = ' ') +
    geom_vline(xintercept = dataAll$HR_plasma[4],linetype="dashed")+
    geom_vline(xintercept = 1,linetype="dotted")+
    geom_text(aes(label=paste("N=",N,sep = ""),y=Cohort_excl),x=5,hjust=0,size=4)+
    theme_classic()+ggtitle("MCI")+
    theme(axis.text=element_text(size=14),axis.title = element_text(size=14),
          legend.position = "none",plot.title = element_text(hjust = 0.5),
          strip.text.y = element_text(size = 10))+
    scale_color_manual(values=c(rep("black",length(name_cohort))))+
    scale_fill_manual(values=c(rep("#f39b7fff",length(name_cohort))))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.ticks = element_line(size=0.1),axis.line = element_line(size=0.1))+
    coord_cartesian(xlim = v_lim)
  #+
    # scale_x_reverse(limits=v_lim)
  
  table_PET=data.frame(matrix(nrow = length(name_cohort)*2,ncol = 6))
  colnames(table_PET)=c("N","HR","HR_l","HR_h","Model","Cohort_excl")#data.frame(beta=NA,se=NA,group=NA,pval=NA,mdl=NA)
  
  i_count=1
  for (i_cohort in 1:length(name_cohort)) {
    df=get(list_cohort[i_cohort])
    table_PET$N[i_count]=df$N[2]
    table_PET$HR[i_count]=df$HR_PET[2]
    table_PET$HR_l[i_count]=df$HR_PET_l[2]
    table_PET$HR_h[i_count]=df$HR_PET_h[2]
    table_PET$Model[i_count]="Single"
    table_PET$Cohort_excl[i_count]=name_cohort[i_cohort]
    
    table_PET$N[i_count+1]=df$N[4]
    table_PET$HR[i_count+1]=df$HR_PET[4]
    table_PET$HR_l[i_count+1]=df$HR_PET_l[4]
    table_PET$HR_h[i_count+1]=df$HR_PET_h[4]
    table_PET$Model[i_count+1]="With plasma"
    table_PET$Cohort_excl[i_count+1]=name_cohort[i_cohort]
    
    i_count=i_count+2
  }
  
  table_PET$HR_h[which((table_PET$HR_h)=="Inf")]=99
  table_PET$HR_h=as.numeric(table_PET$HR_h)
  
  plot_NeoT=table_PET %>% filter(Model=="Single") %>% #group_by(Model) %>%
    mutate(index = reorder(Cohort_excl, N)) %>%
    ggplot( aes(y=index, x=HR, xmin=HR_l, xmax=HR_h)) +
    geom_rect(aes(ymin = -Inf, xmin = dataAll$HR_PET_h[2],
                  ymax = Inf, xmax = dataAll$HR_PET_l[2]),alpha=0.2,
              fill = "lightgrey")+
    geom_point(aes(size=N),shape=23,fill="#3c5488ff",color="#3c5488ff") +
    geom_errorbarh(height=.3,aes(color=Cohort_excl),size=0.5)+
    labs( x=paste("HR"), y = ' ') +
    geom_vline(xintercept = dataAll$HR_PET[2],linetype="dashed")+
    geom_vline(xintercept = 1,linetype="dotted")+
    geom_text(aes(label=paste("N=",N,sep = ""),y=Cohort_excl),x=3.5,hjust=0,size=4)+
    theme_classic()+ggtitle("MCI")+
    theme(axis.text=element_text(size=14),axis.title = element_text(size=14),
          legend.position = "none",plot.title = element_text(hjust = 0.5),
          strip.text.y = element_text(size = 10))+
    scale_color_manual(values=c(rep("black",length(name_cohort))))+
    scale_fill_manual(values=c(rep("#3c5488ff",length(name_cohort))))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.ticks = element_line(size=0.1),axis.line = element_line(size=0.1))+
    coord_cartesian(xlim = v_lim)
  #+
    # scale_x_reverse(limits=v_lim)#
  
  plot_NeoT2=table_PET %>% filter(Model=="With plasma") %>% #group_by(Model) %>%
    mutate(index = reorder(Cohort_excl, N)) %>%
    ggplot( aes(y=index, x=HR, xmin=HR_l, xmax=HR_h)) +
    geom_rect(aes(ymin = -Inf, xmin = dataAll$HR_PET_h[4],
                  ymax = Inf, xmax = dataAll$HR_PET_l[4]),alpha=0.2,
              fill = "lightgrey")+
    geom_errorbarh(height=.3,aes(color=Cohort_excl),size=0.5)+
    geom_point(aes(size=N),shape=23,fill="#8491b4",color="#8491b4") +
    geom_errorbarh(height=.3,aes(color=Cohort_excl),size=0.5)+
    labs( x=paste("HR"), y = ' ') +
    geom_vline(xintercept = dataAll$HR_PET[4],linetype="dashed")+
    geom_vline(xintercept = 1,linetype="dotted")+
    geom_text(aes(label=paste("N=",N,sep = ""),y=Cohort_excl),x=3.5,hjust=0,size=4)+
    theme_classic()+ggtitle("MCI")+
    theme(axis.text=element_text(size=14),axis.title = element_text(size=14),
          legend.position = "none",plot.title = element_text(hjust = 0.5),
          strip.text.y = element_text(size = 10))+
    scale_color_manual(values=c(rep("black",length(name_cohort))))+
    scale_fill_manual(values=c(rep("#8491b4",length(name_cohort))))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.ticks = element_line(size=0.1),axis.line = element_line(size=0.1))+
    coord_cartesian(xlim = v_lim)
  #+
    # scale_x_reverse(limits=v_lim)
  
  ##MTL
  i_count=1
  for (i_cohort in 1:length(name_cohort)) {
    df=get(list_cohort[i_cohort])
    table_PET$N[i_count]=df$N[3]
    table_PET$HR[i_count]=df$HR_PET[3]
    table_PET$HR_l[i_count]=df$HR_PET_l[3]
    table_PET$HR_h[i_count]=df$HR_PET_h[3]
    table_PET$Model[i_count]="Single"
    table_PET$Cohort_excl[i_count]=name_cohort[i_cohort]
    
    table_PET$N[i_count+1]=df$N[5]
    table_PET$HR[i_count+1]=df$HR_PET[5]
    table_PET$HR_l[i_count+1]=df$HR_PET_l[5]
    table_PET$HR_h[i_count+1]=df$HR_PET_h[5]
    table_PET$Model[i_count+1]="With plasma"
    table_PET$Cohort_excl[i_count+1]=name_cohort[i_cohort]
    
    i_count=i_count+2
  }
  
 
  table_PET$HR_h[which((table_PET$HR_h)=="Inf")]=99
  table_PET$HR_h=as.numeric(table_PET$HR_h)
  
  plot_MTL=table_PET %>% filter(Model=="Single") %>% #group_by(Model) %>%
    mutate(index = reorder(Cohort_excl, N)) %>%
    ggplot( aes(y=index, x=HR, xmin=HR_l, xmax=HR_h)) +
    geom_rect(aes(ymin = -Inf, xmin = dataAll$HR_PET_h[3],
                  ymax = Inf, xmax = dataAll$HR_PET_l[3]),alpha=0.2,
              fill = "lightgrey")+
    geom_errorbarh(height=.3,aes(color=Cohort_excl),size=0.5)+
    geom_point(aes(size=N),shape=23,fill="#00a087ff",color="#00a087ff") +
    geom_errorbarh(height=.3,aes(color=Cohort_excl),size=0.5)+
    labs( x=paste("HR"), y = ' ') +
    geom_vline(xintercept = dataAll$HR_PET[3],linetype="dashed")+
    geom_vline(xintercept = 1,linetype="dotted")+
    geom_text(aes(label=paste("N=",N,sep = ""),y=Cohort_excl),x=5,hjust=0,size=4)+
    theme_classic()+ggtitle("MCI")+
    theme(axis.text=element_text(size=14),axis.title = element_text(size=14),
          legend.position = "none",plot.title = element_text(hjust = 0.5),
          strip.text.y = element_text(size = 10))+
    scale_color_manual(values=c(rep("black",length(name_cohort))))+
    scale_fill_manual(values=c(rep("#00a087ff",length(name_cohort))))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.ticks = element_line(size=0.1),axis.line = element_line(size=0.1))+
    coord_cartesian(xlim = v_lim)
    # scale_x_reverse(limits=v_lim)
  
  plot_MTL2=table_PET %>% filter(Model=="With plasma") %>% #group_by(Model) %>%
    mutate(index = reorder(Cohort_excl, N)) %>%
    ggplot( aes(y=index, x=HR, xmin=HR_l, xmax=HR_h)) +
    geom_rect(aes(ymin = -Inf, xmin = dataAll$HR_PET_h[5],
                  ymax = Inf, xmax = dataAll$HR_PET_l[5]),alpha=0.2,
              fill = "lightgrey")+
    geom_point(aes(size=N),shape=23,fill="#91d1c2ff",color="#91d1c2ff") +
    geom_errorbarh(height=.3,aes(color=Cohort_excl),size=0.5)+
    labs( x=paste("HR"), y = ' ') +
    geom_vline(xintercept = dataAll$HR_PET[5],linetype="dashed")+
    geom_vline(xintercept = 1,linetype="dotted")+
    geom_text(aes(label=paste("N=",N,sep = ""),y=Cohort_excl),x=3.5,hjust=0,size=4)+
    theme_classic()+ggtitle("MCI")+
    theme(axis.text=element_text(size=14),axis.title = element_text(size=14),
          legend.position = "none",plot.title = element_text(hjust = 0.5),
          strip.text.y = element_text(size = 10))+
    scale_color_manual(values=c(rep("black",length(name_cohort))))+
    scale_fill_manual(values=c(rep("#91d1c2ff",length(name_cohort))))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.ticks = element_line(size=0.1),axis.line = element_line(size=0.1))+
    coord_cartesian(xlim = v_lim)
  
  
  # pf2=ggarrange(plot_plasma+ rremove("y.text"),  plot_MTL + rremove("y.text"), plot_NeoT+ rremove("y.text"), 
  #               plot_plasma_MTL+ rremove("y.text"), plot_MTL2+ rremove("y.text"),"",
  #               plot_plasma_NeoT+ rremove("y.text"),"",plot_NeoT2+ rremove("y.text"),
  #               # labels = c("A", "HR", "C","D","E"),
  #               ncol = 3, nrow = 3)
  pf2=ggarrange(plot_plasma+rremove("y.text"),  plot_MTL +rremove("y.text"), plot_NeoT+rremove("y.text"),
                plot_plasma_MTL+rremove("y.text"), plot_MTL2+rremove("y.text"),"",
                plot_plasma_NeoT+rremove("y.text"),"",plot_NeoT2+rremove("y.text"),
                # labels = c("A", "B", "C","D","E"),
                ncol = 3, nrow = 3)
  
  
  ggexport(pf2,filename=paste(dir_fig,"/Compare_Kaplan_MCI_Ab_",l_ab,"_adj_",l_adj,"SensTime_NA.pdf",sep = ""))
  


