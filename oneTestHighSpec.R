#One screen algorithm
#High specificity screen (98%)
#Varied sensitivity


# Packages and base set-up ------------------------------------------------

library(tidyverse)
options(scipen = 999)
rm(list=ls())
N=100000#screening population size

vals=c(0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.80,0.85,0.9,0.95,0.98) #sensitivity range for novel screen
startedPrevs=c(0.001,0.0025,0.005,0.01) #considered prevalences

spec1=0.98 #Specificity of novel screen
confSens=0.9 #Sensitivity of confirmatory test
confSpec=0.96 #Specificity of confirmatory test


# Simulate decision tree ---------------------------------------------------
for (p in 1:4){
  prev=startedPrevs[p]
  
  #sens1
  for (i in 1: length(vals)){
    sens1=vals[i]
    
    TBP=prev*N
    TBN=(1-prev)*N
    
    FP1=TBN*(1-spec1)
    FP2=FP1*(1-confSpec)
    TN1=TBN*(spec1)
    TN2=FP1*confSpec
    
    TP1=TBP*sens1
    TP2=TP1*confSens
    FN1=TBP*(1-sens1)
    FN2=TP1*(1-confSens)
    
    totalScreenPos=TP1+FP1
    preTestPrev=(TP1/totalScreenPos)
    percCorTB=TP2/TBP
    
    outIntermediate=t(c(prev,sens1,spec1,TBP, FP1,TP1,FN1,TN1,FP2,TP2,totalScreenPos,preTestPrev,percCorTB))
    outIntermediate=data.frame((outIntermediate))
    
    
    if (i==1& p==1){
      out=outIntermediate
    }else{
      out=rbind(out,outIntermediate)
    }
    
  }
  
}

colnames(out)=c("prev","sens1","spec1","numTB" ,"FP1","TP1","FN1","TN1","FP2","TP2","totalScreenPos","preTestPrev","percCorTB")
outFinal=out
outFinal$prev=factor(outFinal$prev*100)
outFinal$screenRatio=outFinal$totalScreenPos/outFinal$TP1
outFinal$numDiagnostic=outFinal$TP1+outFinal$FP1
outFinal$yield=100*(outFinal$TP2/outFinal$numDiagnostic)
outFinal=outFinal %>% group_by(prev)


# Figure development ------------------------------------------------------

#tested vs found
testedFound=outFinal %>% select(prev,sens1,numDiagnostic,TP2,percCorTB,yield)
testedFound=testedFound %>% gather("group","value",3:6)

TBPos=testedFound%>% filter(group=="TP2")
diags=testedFound%>% filter(group=="numDiagnostic")
corTB=testedFound%>% filter(group=="percCorTB")
yield=testedFound%>% filter(group=="yield")


#######Manuscript Figures
#Pre-screen Prevalence
fig3a=ggplot(outFinal,aes(sens1*100,preTestPrev*100,group=prev,colour=prev))+geom_point()+geom_line()+
  theme_bw()+xlab("Sensitivity(%)")+ylab("Post-screen prevalence (%)")+geom_hline(aes(yintercept = 5, linetype = "5% post-screen prevalence"))+
  labs(subtitle = paste("Specificity of screen at ",spec1*100,"%",sep=""))+
  scale_y_continuous(breaks=seq(0,35,2))+scale_x_continuous(breaks=seq(40,100,5))+theme(text=element_text(size=16))+
  scale_linetype_manual(name="",values = c("5% post-screen prevalence" = "solid"))+
  scale_color_discrete(name="Prevalence (%)")

fig3a

#substitute path
ggsave(file="/Users/adenooy/Library/CloudStorage/OneDrive-Personal/AMC/TB/Screening Paper/Revision 1 - LGH/figures/Screening Figure 3A_R1.svg", plot=fig3a, width=10, height=5)

diags$prevlab=paste("Prev: ",diags$prev,"%",sep="")
corTB$prevlab=paste("Prev: ",corTB$prev,"%",sep="")
TBPos$prevlab=paste("Prev: ",corTB$prev,"%",sep="")


fig3b=ggplot(diags,aes(sens1*100,value,fill="Confirmatory Tests"))+geom_bar(stat="identity")+facet_wrap(.~prevlab)+
  theme_bw()+xlab("Screen Sensitivity(%)")+
  labs(subtitle = paste("Specificity of screen at ",spec1*100,"%",sep=""))+
  scale_x_continuous(breaks=seq(40,100,5))+ scale_color_manual(name="",values=c("black","black"))+ scale_fill_manual(name="",values=c("#B4CF66"))+ 
  theme(legend.position = "bottom",text=element_text(size=16))+ylab("Number of Confirmatory Tests")


fig3b
ggsave(file="/Users/adenooy/Library/CloudStorage/OneDrive-Personal/AMC/TB/Screening Paper/Revision 1 - LGH/figures/Screening Figure 3B_R1.svg", plot=fig3b, width=12, height=6)


coeff=20

fig3c=ggplot(data=TBPos,aes(sens1*100,value,fill="Cases detected"))+geom_bar(stat="identity")+facet_wrap(.~prevlab)+
  geom_line(data=corTB,aes(x=sens1*100,y=value*100*coeff,colour="TB Case Detection (%)"))+
  geom_point(data=corTB,aes(x=sens1*100,y=value*100*coeff,colour="TB Case Detection (%)"))+
  geom_text(data=corTB,aes(x=sens1*100,y=(value)*100*coeff+120,label=value*100,colour="TB Case Detection (%)"))+
  theme_bw()+xlab("Screen Sensitivity(%)")+
  labs(subtitle = paste("Specificity of screen at ",spec1*100,"%",sep=""))+
  scale_x_continuous(breaks=seq(40,100,5))+ scale_color_manual(name="",values=c("black","black"))+ scale_fill_manual(name="",values=c("#146152"))+ 
  theme(legend.position = "bottom",text=element_text(size=16))+ylab("Number of Individuals")+scale_y_continuous(
    name = "Number of positive TB tests",
    sec.axis = sec_axis(~ . /coeff , name = "TB Case Detection (%)")  # transformation
  )

fig3c

ggsave(file="/Users/adenooy/Library/CloudStorage/OneDrive-Personal/AMC/TB/Screening Paper/Revision 1 - LGH/figures/Screening Figure 3C_R1.svg", plot=fig3c, width=14, height=6)


# Additional Costing ------------------------------------------------------
costDiag=8
cost_screen=2

outFinal$diagCost=costDiag*outFinal$numDiagnostic
outFinal$total_cost_screen_2 =outFinal$diagCost+N*cost_screen

outFinal$cost_correct=outFinal$total_cost_screen_2/outFinal$TP2

costing_data=outFinal %>% select(prev,sens1,spec1,cost_correct)

library(RColorBrewer)

costing_data$prev=factor(costing_data$prev)
costing_data$sens1=factor(costing_data$sens1*100)

figS3b=ggplot(costing_data,aes(sens1,prev,fill=cost_correct))+geom_tile(color = "grey50",
                                                         lwd = 0.5,
                                                         linetype = 1)+theme_bw()+xlab("Sensitivity (%)")+ylab("Prevalence (%)")+
  theme(text=element_text(size=16))+
  scale_fill_distiller(
    palette = "RdYlGn",    # reversed green-to-red, use "YlGnBu" etc. for others
    direction = -1,        # -1 to flip colors (green = low, red = high)
    limits = c(0, 6100),    # set fixed min and max
    name = "Cost per correct diagnosis"
  )+geom_text(aes(label=round(cost_correct,0)))+labs(title="One-screen - High-specificity Test")

figS3b

write.csv(costing_data,"/Users/adenooy/Downloads/costing_one_high_spec.csv")




