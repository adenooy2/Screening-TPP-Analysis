#One screen algorithm
#High sensitivity screen (90%)
#Varied specificity


# Packages and base set-up ------------------------------------------------
library(tidyverse)
rm(list=ls())
N=100000 #screening population size

vals=c(0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.80,0.85,0.90,0.95,0.98) #specificity range for novel screen
startedPrevs=c(0.001,0.0025,0.005,0.01) #considered prevalences

sens1=0.9 #sensitivity of novel screen
confSens=0.9 #Sensitivity of confirmatory test
confSpec=0.96 #Specificity of confirmatory test


# Simulate decision tree ---------------------------------------------------
for (p in 1:4){
  prev=startedPrevs[p]
  
  #spec1
  for (i in 1: length(vals)){
    spec1=vals[i]
    
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

outFinal$algType="OneStep"
outFinal$testConstraint="HighSens"
outFinal$totalCohort=N
outFinal$confirmatoryPos=outFinal$TP2+outFinal$FP2
outFinal$cxrScreenPos="-"
outFinal$screen1Pos=outFinal$totalScreenPos
outFinal$TBP=outFinal$numTB
outTable=outFinal %>% select(algType,testConstraint,prev,sens1,spec1,totalCohort,screen1Pos,cxrScreenPos,confirmatoryPos,TP2,TBP,percCorTB,preTestPrev)


# Figure development ------------------------------------------------------
#tested vs found
testedFound=outFinal %>% select(prev,spec1,numDiagnostic,TP2,percCorTB,yield)
testedFound=testedFound %>% gather("group","value",3:6)

TBPos=testedFound%>% filter(group=="TP2")
diags=testedFound%>% filter(group=="numDiagnostic")
corTB=testedFound%>% filter(group=="percCorTB")
yield=testedFound%>% filter(group=="yield")


##############Manuscript figures
fig2a=ggplot(outFinal,aes(spec1*100,preTestPrev*100,group=prev,colour=prev))+geom_point()+geom_line()+
  theme_bw()+xlab("Specificity(%)")+ylab("Post-screen prevalence (%)")+geom_hline(aes(yintercept = 5, linetype = "5% post-screen prevalence"))+
  labs(subtitle = paste("Sensitivity of screen at ",sens1*100,"%",sep=""))+
  scale_y_continuous(breaks=seq(0,35,2))+scale_x_continuous(breaks=seq(40,100,5))+theme(text = element_text(size = 16))+
  scale_linetype_manual(name="",values = c("5% post-screen prevalence" = "solid"))+
  scale_color_discrete(name="Prevalence (%)")
  

fig2a

ggsave(file="/Users/adenooy/Library/CloudStorage/OneDrive-Personal/AMC/TB/Screening Paper/Revision 1 - LGH/figures/Screening Figure 2A_R1.svg", plot=fig2a, width=10, height=5)

diags$prevlab=paste("Prev: ",diags$prev,"%",sep="")
corTB$prevlab=paste("Prev: ",corTB$prev,"%",sep="")
TBPos$prevlab=paste("Prev: ",TBPos$prev,"%",sep="")



fig2b=ggplot(diags ,aes(spec1*100,value,fill="Confirmatory Tests"))+geom_bar(stat="identity")+facet_wrap(.~prevlab)+
  theme_bw()+xlab("Screen Specificity(%)")+
  labs(subtitle = paste("Sensitivity of screen at ",sens1*100,"%",sep=""))+
  scale_x_continuous(breaks=seq(40,95,5))+ scale_color_manual(name="",values=c("black"))+ scale_fill_manual(name="",values=c("#B4CF66"))+ 
  theme(legend.position = "bottom",text=element_text(size=14))+ylab("Number of Confirmatory Tests")+theme(text = element_text(size = 16))

fig2b


ggsave(file="/Users/adenooy/Library/CloudStorage/OneDrive-Personal/AMC/TB/Screening Paper/Revision 1 - LGH/figures/Screening Figure 2B_R1.svg", plot=fig2b, width=10, height=5)

coeff=10

fig2c=ggplot(data=TBPos,aes(spec1*100,value,fill="Cases detected"),stat="identity")+geom_bar(stat="identity")+facet_wrap(.~prevlab)+
  geom_line(data=corTB,aes(x=spec1*100,y=value*100*coeff,colour="TB Case Detection (%)"))+
  geom_point(data=corTB,aes(x=spec1*100,y=value*100*coeff,colour="TB Case Detection (%)"))+
  geom_text(data=corTB,aes(x=spec1*100,y=(value)*100*coeff+100,label=round(value*100,0),colour="TB Case Detection (%)"))+
  theme_bw()+xlab("Screen Specificity(%)")+
  labs(subtitle = paste("Sensitivity of screen at ",sens1*100,"%",sep=""))+
  scale_x_continuous(breaks=seq(40,95,5))+ scale_color_manual(name="",values=c("black","black"))+ scale_fill_manual(name="",values=c("#146152"))+ 
  theme(legend.position = "bottom",text=element_text(size=14))+ylab("Number of Individuals")+theme(text = element_text(size = 16))+
  theme(legend.position = "bottom",text=element_text(size=16))+ylab("Number of Individuals")+scale_y_continuous(
    name = "Number of positive TB tests",
    sec.axis = sec_axis(~ . /coeff , name = "TB Case Detection (%)")  # transformation
  )


fig2c

ggsave(file="/Users/adenooy/Library/CloudStorage/OneDrive-Personal/AMC/TB/Screening Paper/Revision 1 - LGH/figures/Screening Figure 2C_R1.svg", plot=fig2c, width=10, height=6)



# Additional Costing ------------------------------------------------------
costDiag=8
cost_screen=3

outFinal$diagCost=costDiag*outFinal$numDiagnostic
outFinal$total_cost_screen_2 =outFinal$diagCost+N*cost_screen

outFinal$cost_correct=outFinal$total_cost_screen_2/outFinal$TP2

costing_data=outFinal %>% select(prev,sens1,spec1,cost_correct)

library(RColorBrewer)

costing_data$prev=factor(costing_data$prev)
costing_data$spec1=factor(costing_data$spec1*100)

figS3a=ggplot(costing_data,aes(spec1,prev,fill=cost_correct))+geom_tile(color = "grey50",
                                                                        lwd = 0.5,
                                                                        linetype = 1)+theme_bw()+xlab("Specificity (%)")+ylab("Prevalence (%)")+
  theme(text=element_text(size=16))+
  scale_fill_distiller(
    palette = "RdYlGn",    # reversed green-to-red, use "YlGnBu" etc. for others
    direction = -1,     # set fixed min and max
    name = "Cost per correct diagnosis"
  )+geom_text(aes(label=round(cost_correct,0)))+labs(title="One-screen - High-sensitivity Test")

figS3a

write.csv(costing_data,"/Users/adenooy/Downloads/costing_one_high_sens.csv")





