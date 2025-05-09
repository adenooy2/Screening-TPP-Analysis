#Two screen algorithm
#High sensitivity screen (90%)
#Varied specificity


# Packages and base set-up ------------------------------------------------
rm(list=ls())
library(tidyverse)
N=100000 #screening population size
vals=c(0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,0.98) #specificity range
startedPrevs=c(0.001,0.0025,0.005,0.01) #considered prevalences

sens1=0.9 #sensitivity of nove screen
confSens=0.9 #Sensitivity of confirmatory test
confSpec=0.96 #Specificity of confirmatory test


cxrSpec=0.70 #specificity of second screen (CXR)
cxrSens=0.85 #sensitivity of second screen (CXR)

# Simulate decision tree ---------------------------------------------------
for (p in 1:4){
  prev=startedPrevs[p]
  
  #sens1
  for (i in 1: length(vals)){
    spec1=vals[i]
    
    TBP=prev*N
    TBN=(1-prev)*N
    
    FP1=TBN*(1-spec1)
    FP2=FP1*(1-cxrSpec)
    TN1=TBN*(spec1)
    TN2=FP1*cxrSpec
    
    TP1=TBP*sens1
    TP2=TP1*cxrSens
    FN1=TBP*(1-sens1)
    FN2=TP1*(1-cxrSens)
    
    FP3=FP2*(1-confSpec)
    TB3=FP2*confSpec
    TP3=TP2*confSens
    FN3=TP2*(1-confSens)
    
    screen1Pos=TP1+FP1
    cxrScreenPos=TP2+FP2
    preTestPrev=(TP2/cxrScreenPos)
    percCorTB=TP3/TBP
    
    
    outIntermediate=t(c(prev,sens1,spec1,TBP, FP1,TP1,FN1,TN1,FP2,TP2,FP3,TP3,screen1Pos,cxrScreenPos,preTestPrev,percCorTB))
    outIntermediate=data.frame((outIntermediate))
    
    
    if (i==1& p==1){
      out=outIntermediate
    }else{
      out=rbind(out,outIntermediate)
    }
    
    
  }
  
}



colnames(out)=c("prev","sens1","spec1","TBP", "FP1","TP1","FN1","TN1","FP2","TP2","FP3","TP3","screen1Pos","cxrScreenPos","preTestPrev","percCorTB")

outFinal=out

outFinal$prev=factor(outFinal$prev*100)
outFinal$numDiagnostic=outFinal$TP2+outFinal$FP2
outFinal$yield=100*(outFinal$TP3/outFinal$numDiagnostic)
outFinal=outFinal %>% group_by(prev)
outFinal$algType="TwoStep"
outFinal$testConstraint="HighSens"
outFinal$totalCohort=N
outFinal$confirmatoryPos=outFinal$TP3+outFinal$FP3

# Figure development ------------------------------------------------------

#tested vs found
testedFound=outFinal %>% select(prev,spec1,numDiagnostic,TP3,screen1Pos,percCorTB,yield)
testedFound=testedFound %>% gather("group","value",3:7)

TBPos=testedFound%>% filter(group=="TP3")
diags=testedFound%>% filter(group=="numDiagnostic")
corTB=testedFound%>% filter(group=="percCorTB")
yield=testedFound%>% filter(group=="yield")
cxrs=testedFound %>% filter(group=="screen1Pos")

diags$prevlab=paste("Prev: ",diags$prev,"%",sep="")
corTB$prevlab=paste("Prev: ",corTB$prev,"%",sep="")
TBPos$prevlab=paste("Prev: ",corTB$prev,"%",sep="")
cxrs$prevlab=paste("Prev: ",cxrs$prev,"%",sep="")

##Manuscript visuals
#Pre-diagnostic Prevalence
figS1a=ggplot(outFinal,aes(spec1*100,preTestPrev*100,group=prev,colour=prev))+geom_point()+geom_line()+
  theme_bw()+xlab("Specificity(%)")+ylab("Post-screen prevalence (%)")+geom_hline(aes(yintercept = 5, linetype = "5% post-screen prevalence"))+
  labs(subtitle = paste("Sensitivity of screen at: ",sens1*100,"%",sep=""))+
  scale_y_continuous(breaks=seq(0,90,10))+scale_x_continuous(breaks=seq(40,100,5))+theme(legend.position = "right")+
  theme(text=element_text(size=16))+scale_linetype_manual(name="",values = c("5% post-screen prevalence" = "solid"))+
  scale_color_discrete(name="Prevalence (%)")

figS1a

ggsave(file="/Users/adenooy/Library/CloudStorage/OneDrive-Personal/AMC/TB/Screening Paper/Revision 1 - LGH/figures/Screening Figure S1A_R1.svg", plot=figS1a, width=10, height=5)




figS1b=ggplot(data=cxrs ,aes(spec1*100,value,fill="# Chest Xrays"))+geom_bar(stat="identity")+facet_wrap(.~prevlab)+
  geom_bar(data=diags,aes(spec1*100,value,fill="# Confirmatory Tests"),stat="identity")+
  theme_bw()+xlab("Screen Specificity(%)")+
  labs(subtitle = "Sensitivity of screen at 90% ")+
  scale_color_manual(name="",values=c("black","black"))+ scale_fill_manual(name="",values=c("#45C4B0","#B4CF66"))+ 
  theme(legend.position = "bottom",text=element_text(size=16))+ylab("Number of Tests")+scale_x_continuous(breaks=seq(40,100,5))

figS1b
ggsave(file="/Users/adenooy/Library/CloudStorage/OneDrive-Personal/AMC/TB/Screening Paper/Revision 1 - LGH/figures/Screening Figure S1B_R1.svg", plot=figS1b, width=10, height=5)


coeff=10
figS1c=ggplot(data=TBPos,aes(spec1*100,value,fill="Cases detected"))+geom_bar(stat="identity")+facet_wrap(.~prevlab)+
  geom_line(data=corTB,aes(x=spec1*100,y=value*100*coeff,colour="TB Case Detection (%)"))+
  geom_point(data=corTB,aes(x=spec1*100,y=value*100*coeff,colour="TB Case Detection (%)"))+
  geom_text(data=corTB,aes(x=spec1*100,y=(value)*100*coeff+100,label=round(value*100,0)))+
  theme_bw()+xlab("Screen Specificity(%)")+
  labs(subtitle = "Sensitivity of screen at 90% ")+
  scale_color_manual(name="",values=c("black","black"))+ scale_fill_manual(name="",values=c("#146152"))+ 
  theme(legend.position = "bottom",text=element_text(size=16))+ylab("Number of Individuals")+scale_x_continuous(breaks=seq(40,100,5))+
  scale_y_continuous(
    name = "Number of positive TB tests",
    sec.axis = sec_axis(~ . /coeff , name = "TB Case Detection (%)")  # transformation
  )

figS1c


ggsave(file="/Users/adenooy/Library/CloudStorage/OneDrive-Personal/AMC/TB/Screening Paper/Revision 1 - LGH/figures/Screening Figure S1C_R1.svg", plot=figS1c, width=10, height=5)


# Additional Costing ------------------------------------------------------
costDiag=8
cost_screen=2
cost_cxr=3.8

outFinal$diagCost=costDiag*outFinal$numDiagnostic
outFinal$cxrCost=cost_cxr*outFinal$screen1Pos
outFinal$total_cost_screen =outFinal$cxrCost+outFinal$diagCost+N*cost_screen

outFinal$cost_correct=outFinal$total_cost_screen/outFinal$TP3

costing_data=outFinal %>% select(prev,sens1,spec1,cost_correct)

library(RColorBrewer)

costing_data$prev=factor(costing_data$prev)
costing_data$spec1=factor(costing_data$spec1*100)

figS3c=ggplot(costing_data,aes(spec1,prev,fill=cost_correct))+geom_tile(color = "grey50",
                                                                        lwd = 0.5,
                                                                        linetype = 1)+theme_bw()+xlab("Specificity (%)")+ylab("Prevalence (%)")+
  theme(text=element_text(size=16))+
  scale_fill_distiller(
    palette = "RdYlGn",    # reversed green-to-red, use "YlGnBu" etc. for others
    direction = -1,        # -1 to flip colors (green = low, red = high)
    limits = c(0, 8500),    # set fixed min and max
    name = "Cost per correct diagnosis"
  )+geom_text(aes(label=round(cost_correct,0)))+labs(title="Two-screen - High-sensitivity Test")

figS3c

write.csv(costing_data,"/Users/adenooy/Downloads/costing_two_high_sens.csv")




