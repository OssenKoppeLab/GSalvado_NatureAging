dir_res="~/Desktop/Lund/20_TauPET_plasma/Results_Merged/20240917_Rev/Two_step_sample_red_conv"
dir_fig="~/Desktop/Lund/20_TauPET_plasma/Figures_Merged/20240917_Rev/Two_step_sample_red_conv"

c_colors=c(
           "Plasma"="#e64b35ff",
           "p-MTL_plasma"="#00a087ff",
           "p-NeoT_plasma"="#3c5488ff")

c_out=c("convMCI")#"mmse_score",
v_power=0.8
v_change=0.7
l_ab="all"#"Abpos"#
l_DX="allEtiol"
l_adj="raw"#"combat"
l_out="_noCovs"#"noOutlier"
c_method_plasma=c("Q2_Q4","Q3_Q4","Q4")#"zscore2"
c_method_PET=c("Q2_Q4","Q3_Q4","Q4")#"zscore2"

for (i_out in 1:length(c_out)) {
  
  for (i_method in 1:3) {
    
    results_Q2_Q4=readxl::read_xlsx(paste(dir_res,"/Plot_sampleRed_",c_out[i_out],"_power_",v_power,"_red_",v_change,"_Plasma_",c_method_plasma[i_method],"_PET_Q2_Q4_Ab_",l_ab,"_adj_",l_adj,"_2steps_20241015.xlsx",sep = ""))
    results_Q3_Q4=readxl::read_xlsx(paste(dir_res,"/Plot_sampleRed_",c_out[i_out],"_power_",v_power,"_red_",v_change,"_Plasma_",c_method_plasma[i_method],"_PET_Q3_Q4_Ab_",l_ab,"_adj_",l_adj,"_2steps_20241015.xlsx",sep = ""))
    results_Q4=readxl::read_xlsx(paste(dir_res,"/Plot_sampleRed_",c_out[i_out],"_power_",v_power,"_red_",v_change,"_Plasma_",c_method_plasma[i_method],"_PET_Q4_Ab_",l_ab,"_adj_",l_adj,"_2steps_20241015.xlsx",sep = ""))
    
    results_PET=results_Q2_Q4[which(results_Q2_Q4$Model=="p-MTL_plasma" | results_Q2_Q4$Model=="p-NeoT_plasma"),]
    results_PET=rbind(results_PET,results_Q3_Q4[which(results_Q3_Q4=="p-MTL_plasma" | results_Q3_Q4=="p-NeoT_plasma"),],results_Q4[which(results_Q4=="p-MTL_plasma" | results_Q4=="p-NeoT_plasma"),])
    results_PET$Quant=c(rep("Q2-Q4",2),rep("Q3-Q4",2),rep("Q4",2))
    results_PET$Quant=as.factor(results_PET$Quant)
    
    p_PET1=ggplot(results_PET[which(results_PET$Model=="p-MTL_plasma"),], aes(x = Quant, y = Perc,fill=Model)) +
      geom_bar(stat="identity", position=position_dodge()) +
      geom_errorbar(aes(ymin=Perc_l, ymax=Perc_h), width=.2,
                    position=position_dodge(.9))+
      scale_fill_manual(values=c_colors)+ylab("Final sample (%)") + xlab("Quantiles")+
      # geom_text(aes(label = paste(round(Perc,0),"%",sep = ""),x=Model,y=Perc_h+25), stat="identity", vjust=0.1,size=4)+
      theme_classic()+theme(legend.position = "none",axis.text=element_text(size=8),
                            axis.title=element_text(size=10,face="bold"),legend.text = element_text(size=6),
                            axis.line.x = element_line(size = 0.3),axis.line.y = element_line(size = 0.3),
                            axis.ticks.x=element_line(size = 0.3),axis.ticks.y = element_line(size = 0.3),
                            legend.key.size = unit(0.2, 'cm'),legend.title = element_blank())+
      geom_hline(yintercept=100,linetype="dotted")+coord_cartesian(ylim = c(0, 100))

    p_PET2=ggplot(results_PET[which(results_PET$Model=="p-NeoT_plasma"),], aes(x = Quant, y = Perc,fill=Model)) +
      geom_bar(stat="identity", position=position_dodge()) +
      geom_errorbar(aes(ymin=Perc_l, ymax=Perc_h), width=.2,
                    position=position_dodge(.9))+
      scale_fill_manual(values=c_colors)+ylab("Final sample (%)") + xlab("Quantiles")+
      # geom_text(aes(label = paste(round(Perc,0),"%",sep = ""),x=Model,y=Perc_h+25), stat="identity", vjust=0.1,size=4)+
      theme_classic()+theme(legend.position = "none",axis.text=element_text(size=8),
                            axis.title=element_text(size=10,face="bold"),legend.text = element_text(size=6),
                            axis.line.x = element_line(size = 0.3),axis.line.y = element_line(size = 0.3),
                            axis.ticks.x=element_line(size = 0.3),axis.ticks.y = element_line(size = 0.3),
                            legend.key.size = unit(0.2, 'cm'),legend.title = element_blank())+
      geom_hline(yintercept=100,linetype="dotted")+coord_cartesian(ylim = c(0, 100))
    
    # results_Q2_Q4=readxl::read_xlsx(paste(dir_res,"/SampleRed_Plasma_",c_method_plasma[i_method],"_PET_Q2_Q4_Ab_",l_ab,"_adj_",l_adj,"_power_",v_power,"_change_",v_change,"_proportion.xlsx",sep = ""))
    # results_Q3_Q4=readxl::read_xlsx(paste(dir_res,"/SampleRed_Plasma_",c_method_plasma[i_method],"_PET_Q3_Q4_Ab_",l_ab,"_adj_",l_adj,"_power_",v_power,"_change_",v_change,"_proportion.xlsx",sep = ""))
    # results_Q4=readxl::read_xlsx(paste(dir_res,"/SampleRed_Plasma_",c_method_plasma[i_method],"_PET_Q4_Ab_",l_ab,"_adj_",l_adj,"_power_",v_power,"_change_",v_change,"_proportion.xlsx",sep = ""))
    
    # results_plot=results_Q2_Q4[which(results_Q2_Q4$Outcome==c_out[i_out]),1:5]
    # results_plot=rbind(results_plot,results_Q3_Q4[which(results_Q3_Q4$Outcome==c_out[i_out]),1:5],results_Q4[which(results_Q4$Outcome==c_out[i_out]),1:5])
    # results_plot$Quant=c("Q2-Q4","Q3-Q4","Q4")
    # results_plot$Quant=as.factor(results_plot$Quant)
    
    results_line=data.frame(matrix(nrow = 18,ncol = 4))
    colnames(results_line)=c("Quant","step",
                             "N","col")
    results_line$N=c(100,100,100,results_Q2_Q4$Perc[1],results_Q2_Q4$Perc[1],results_Q2_Q4$Perc[1],
                     results_Q2_Q4$Perc[2],results_Q3_Q4$Perc[2],results_Q4$Perc[2],
                     100,100,100,results_Q2_Q4$Perc[1],results_Q2_Q4$Perc[1],results_Q2_Q4$Perc[1],
                     results_Q2_Q4$Perc[3],results_Q3_Q4$Perc[3],results_Q4$Perc[3])
    results_line$N=as.numeric(results_line$N)
    results_line$step=rep(c(rep(0,3),rep(1,3),rep(2,3)),2)
    results_line$step=as.factor(results_line$step)
    results_line$Quant=c(rep(c("Q2-Q4","Q3-Q4","Q4"),3),rep(c("Q2-Q4_","Q3-Q4_","Q4_"),3))
    results_line$col=c(rep("Original",3),rep("Plasma",3),rep("MTL",3),rep("Original",3),rep("Plasma",3),rep("NeoT",3))#
    
    plot_line=ggplot(results_line[which(results_line$col!="NeoT"),],aes(x=step,y=N,group=Quant,fill=step,color=interaction(step,col)))+
      geom_line()+   #geom_errorbar(aes(ymin=N_l/100, ymax=N_h/100), width=.2,alpha=0.5)+
      geom_point(size = 2, shape = 16) +
      geom_text(aes(label = paste(round(N,0),"%",sep = ""),x=step,y=N-7), stat="identity", vjust=0.1,size=4,color="black")+
      theme_classic()+ylab("Participants (%)")+xlab("Step")+
      scale_color_manual(values = c("0.Original"="#d3d3d3",
                                    "1.Plasma"="#e64b35ff",
                                    "2.MTL"="#00a087ff",
                                    "2.NeoT"="#3c5488ff"))+
      theme(legend.position = "none",axis.text=element_text(size=8),
            axis.title=element_text(size=10,face="bold"),legend.text = element_text(size=6),
            axis.line.x = element_line(size = 0.3),axis.line.y = element_line(size = 0.3),
            axis.ticks.x=element_line(size = 0.3),axis.ticks.y = element_line(size = 0.3),
            legend.key.size = unit(0.2, 'cm'),legend.title = element_blank())+
      geom_hline(yintercept=results_line$N[1],linetype="dotted")+
    scale_y_continuous(limits=c(0,100))
    
    plot_line2=ggplot(results_line[which(results_line$col!="MTL"),],aes(x=step,y=N,group=Quant,fill=step,color=interaction(step,col)))+
      geom_line()+   #geom_errorbar(aes(ymin=N_l/100, ymax=N_h/100), width=.2,alpha=0.5)+
      geom_point(size = 2, shape = 16) +
      geom_text(aes(label = paste(round(N,0),"%",sep = ""),x=step,y=N), stat="identity", vjust=0.1,size=4,color="black")+
      theme_classic()+ylab("Participants (%)")+xlab("Step")+
      scale_color_manual(values = c("0.Original"="#d3d3d3",
                                    "1.Plasma"="#e64b35ff",
                                    "2.MTL"="#00a087ff",
                                    "2.NeoT"="#3c5488ff"))+
      theme(legend.position = "none",axis.text=element_text(size=8),
            axis.title=element_text(size=10,face="bold"),legend.text = element_text(size=6),
            axis.line.x = element_line(size = 0.3),axis.line.y = element_line(size = 0.3),
            axis.ticks.x=element_line(size = 0.3),axis.ticks.y = element_line(size = 0.3),
            legend.key.size = unit(0.2, 'cm'),legend.title = element_blank())+
      geom_hline(yintercept=results_line$N[1],linetype="dotted")+
      scale_y_continuous(limits=c(0,100))
    
    pf2=ggarrange(plot_line+rremove("ylab"),plot_line2+rremove("ylab"),p_PET1+rremove("ylab"),p_PET2+rremove("ylab"),
                  # labels = c("a","","b"),
                  ncol = 4, nrow = 4)
    ggexport(pf2,filename=paste(dir_fig,"/SampleRed_step2_conv_perc_",c_out[i_out],"_Plasma_",c_method_plasma[i_method],"_Ab_",l_ab,"_adj_",l_adj,"_power_",v_power,"_change_",v_change,"_2steps_20241015.pdf",sep = ""))
    
  }
}
  


 for (i_out in 1:length(c_out)) {
  
  results_Q2_Q4=readxl::read_xlsx(paste(dir_res,"/Plot_sampleRed_",c_out[i_out],"_power_",v_power,"_red_",v_change,"_Plasma_Q2_Q4_PET_Q2_Q4_Ab_",l_ab,"_adj_",l_adj,"_2steps_20241015.xlsx",sep = ""))
  results_Q3_Q4=readxl::read_xlsx(paste(dir_res,"/Plot_sampleRed_",c_out[i_out],"_power_",v_power,"_red_",v_change,"_Plasma_Q3_Q4_PET_Q2_Q4_Ab_",l_ab,"_adj_",l_adj,"_2steps_20241015.xlsx",sep = ""))
  results_Q4=readxl::read_xlsx(paste(dir_res,"/Plot_sampleRed_",c_out[i_out],"_power_",v_power,"_red_",v_change,"_Plasma_Q4_PET_Q2_Q4_Ab_",l_ab,"_adj_",l_adj,"_2steps_20241015.xlsx",sep = ""))
  
  results_plasma=results_Q2_Q4[which(results_Q2_Q4$Model=="Plasma"),]
  results_plasma=rbind(results_plasma,results_Q3_Q4[which(results_Q3_Q4$Model=="Plasma"),],results_Q4[which(results_Q4$Model=="Plasma"),])
  results_plasma$Quant=c("Q2-Q4","Q3-Q4","Q4")
  results_plasma$Quant=as.factor(results_plasma$Quant)
  
  p_plasma=ggplot(results_plasma, aes(x = Quant, y = Perc,fill=Model)) +
    geom_bar(stat="identity", position=position_dodge()) +
    geom_errorbar(aes(ymin=Perc_l, ymax=Perc_h), width=.2,
                  position=position_dodge(.9))+
    scale_fill_manual(values=c_colors)+ylab("Final sample (%)") + xlab("Quantiles")+
    # geom_text(aes(label = paste(round(Perc,0),"%",sep = ""),x=Model,y=Perc_h+25), stat="identity", vjust=0.1,size=4)+
    theme_classic()+theme(legend.position = "none",axis.text=element_text(size=8),
                          axis.title=element_text(size=10,face="bold"),legend.text = element_text(size=6),
                          axis.line.x = element_line(size = 0.3),axis.line.y = element_line(size = 0.3),
                          axis.ticks.x=element_line(size = 0.3),axis.ticks.y = element_line(size = 0.3),
                          legend.key.size = unit(0.2, 'cm'),legend.title = element_blank())+
       geom_hline(yintercept=100,linetype="dotted")#+
    # scale_y_continuous(limits=c(0,130))
  
   pf2=ggarrange(p_plasma+rremove("ylab"),
                labels = c("a"),
                ncol = 4, nrow = 4)
  ggexport(pf2,filename=paste(dir_fig,"/SampleRed_step1_conv_",c_out[i_out],"_Ab_",l_ab,"_adj_",l_adj,"_power_",v_power,"_change_",v_change,"_2steps_20241015.pdf",sep = ""))
  
}

####lines
# 
# c_colors=c(
#   "Plasma"="#e64b35ff",
#   "p-MTL_plasma"="#00a087ff",
#   "p-NeoT_plasma"="#3c5488ff")
# 
# c_out=c("mmse_score","mpacc")
# v_power=0.8
# v_change=0.3
# l_ab="all"#"Abpos"#
# l_DX="allEtiol"
# l_adj="raw"#"combat"
# l_out="_noCovs"#"noOutlier"
# c_method=c("Q2_Q4","Q3_Q4","Q4")#"zscore2"
# 
# 
#   results_Q2_Q4=readxl::read_xlsx(paste(dir_res,"/SampleRed_Q2_Q4_Ab_",l_ab,"_adj_",l_adj,"_power_",v_power,"_change_",v_change,"_proportion.xlsx",sep = ""))
#   results_Q3_Q4=readxl::read_xlsx(paste(dir_res,"/SampleRed_Q3_Q4_Ab_",l_ab,"_adj_",l_adj,"_power_",v_power,"_change_",v_change,"_proportion.xlsx",sep = ""))
#   results_Q4=readxl::read_xlsx(paste(dir_res,"/SampleRed_Q4_Ab_",l_ab,"_adj_",l_adj,"_power_",v_power,"_change_",v_change,"_proportion.xlsx",sep = ""))
#   
#   for (i_out in 1:length(c_out)) {
#     
#   results_plot=results_Q2_Q4[which(results_Q2_Q4$Outcome==c_out[i_out]),1:5]
#   results_plot=rbind(results_plot,results_Q3_Q4[which(results_Q3_Q4$Outcome==c_out[i_out]),1:5],results_Q4[which(results_Q4$Outcome==c_out[i_out]),1:5])
#   results_plot$Quant=c("Q2-Q4","Q3-Q4","Q4")
#   results_plot$Quant=as.factor(results_plot$Quant)
#   
#   results_line=data.frame(matrix(nrow = 18,ncol = 4))
#   colnames(results_line)=c("Quant","step",
#                            "N","col")
#   results_line$N=c(results_plot$N_or,results_plot$N_plasma,results_plot$N_plasmaMTL,results_plot$N_or,results_plot$N_plasma,results_plot$N_plasmaNeoT)
#   results_line$N=as.numeric(results_line$N)
#   results_line$step=rep(c(rep(0,3),rep(1,3),rep(2,3)),2)
#   results_line$step=as.factor(results_line$step)
#   results_line$Quant=c(rep(c("Q2-Q4","Q3-Q4","Q4"),3),rep(c("Q2-Q4_","Q3-Q4_","Q4_"),3))
#   results_line$col=c(rep("Original",3),rep("Plasma",3),rep("MTL",3),rep("Original",3),rep("Plasma",3),rep("NeoT",3))#
#   
#   plot_line=ggplot(results_line,aes(x=step,y=N,group=Quant,fill=step,color=interaction(step,col)))+
#     geom_line()+   #geom_errorbar(aes(ymin=N_l/100, ymax=N_h/100), width=.2,alpha=0.5)+
#     geom_point(size = 2, shape = 16) +
#     geom_text(aes(label = paste(round(N,0)),x=step,y=N+25), stat="identity", vjust=0.1,size=4,color="black")+
#     theme_classic()+
#     scale_color_manual(values = c("0.Original"="#d3d3d3",
#                                  "1.Plasma"="#e64b35ff",
#                                   "2.MTL"="#00a087ff",
#                                  "2.NeoT"="#3c5488ff"))+
#     theme(legend.position = "none",axis.text=element_text(size=8),
#           axis.title=element_text(size=10,face="bold"),legend.text = element_text(size=6),
#           axis.line.x = element_line(size = 0.3),axis.line.y = element_line(size = 0.3),
#           axis.ticks.x=element_line(size = 0.3),axis.ticks.y = element_line(size = 0.3),
#           legend.key.size = unit(0.2, 'cm'),legend.title = element_blank())+
#     geom_hline(yintercept=results_line$N[1],linetype="dotted")#+
#   # scale_y_continuous(limits=c(0,130))
#   
#   pf2=ggarrange(plot_line+rremove("ylab"),
#                 labels = c("c"),
#                 ncol = 3, nrow = 3)
#   ggexport(pf2,filename=paste(dir_fig,"/SampleRed_line__Plasma_",c_method_plasma[i_method],"_",c_out[i_out],"_Ab_",l_ab,"_adj_",l_adj,"_power_",v_power,"_change_",v_change,".pdf",sep = ""))
#   
#   }
# }
# 
