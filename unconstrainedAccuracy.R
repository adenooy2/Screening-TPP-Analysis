#Unconstrained accuracy analysis
#Both sensitivity and specificity vary simultaneously


# Packages and base set-up ------------------------------------------------
library(tidyverse)
library(svglite)
options(scipen = 999)
rm(list=ls())
N=100000 #screening population size

vals=c(0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.80,0.85,0.90,0.95)
startedPrevs=c(0.001,0.0025,0.005,0.01)

confSens=0.9
confSpec=0.96


# One Screen Algorithm ----------------------------------------------------

for (p in 1:4){
  prev=startedPrevs[p]
  
  #sens1
  for (i in 1: length(vals)){
    for (j in 1:length(vals)){
      
    
    sens1=vals[i]
    spec1=vals[j]
    
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
}



colnames(out)=c("prev","sens1","spec1","numTB" ,"FP1","TP1","FN1","TN1","FP2","TP2","totalScreenPos","preTestPrev","percCorTB")

outFinal=out
outFinal$prev=factor(outFinal$prev*100)
outFinal$screenRatio=outFinal$totalScreenPos/outFinal$TP1
outFinal$numDiagnostic=outFinal$TP1+outFinal$FP1
outFinal$yield=100*(outFinal$TP2/outFinal$numDiagnostic)

oneScreen=outFinal
oneScreen$alg="one_screen"
###acceptable combinations
accCombOneSc=oneScreen %>% select(prev, alg,sens1,spec1,preTestPrev,percCorTB)


# Two Screen Algorithm ----------------------------------------------------

cxrSpec=0.70
cxrSens=0.85

screenCost=7

for (p in 1:4){
  prev=startedPrevs[p]
  
  #sens1
  for (i in 1: length(vals)){
    for (j in 1:length(vals)){
    spec1=vals[i]
    sens1=vals[j]
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
  
}



colnames(out)=c("prev","sens1","spec1","TBP", "FP1","TP1","FN1","TN1","FP2","TP2","FP3","TP3","screen1Pos","cxrScreenPos","preTestPrev","percCorTB")

outFinal=out

outFinal$prev=factor(outFinal$prev*100)
outFinal$numDiagnostic=outFinal$TP2+outFinal$FP2
outFinal$yield=100*(outFinal$TP3/outFinal$numDiagnostic)
outFinal=outFinal %>% group_by(prev)

twoScreen=outFinal
twoScreen$alg="two_screen"
accCombTwoSc=twoScreen  %>% select(prev, alg,sens1,spec1,preTestPrev,percCorTB)


###Final Acceptable Combinations for both algorithms
final=rbind(accCombOneSc,accCombTwoSc)
final=final %>% filter(sens1>0.68)
final=final %>% filter(spec1>=0.5)
final$prev=paste("Prev: ",final$prev,"%",sep="")

final$alg_labs="One Screen"
final$alg_labs[final$alg=="two_screen"]="Two Screen"

final2=final 
final2$sens1=final2$sens1*100
final2$spec1=final2$spec1*100 
final2$prod=-final2$sens1*final2$spec1
final2$sens1=factor(final2$sens1)
final2$spec1=factor(final2$spec1)


final2$include=0
final2$include[final2$percCorTB>0.6&final2$preTestPrev>0.05]=1
final2$CDR=final2$percCorTB*100
final2$CDR[final2$include==0]=NA

final2$`Post-screen prevalence`=final2$preTestPrev*100
final2$`Post-screen prevalence`[final2$include==0]=NA


library(RColorBrewer)

fig4a=ggplot(final2,aes(sens1,spec1,fill=CDR))+geom_tile(color = "grey50",
                                                   lwd = 0.5,
                                                   linetype = 1)+facet_grid(cols=vars(prev),rows=vars(alg_labs))+theme_bw()+xlab("Sensitivity (%)")+ylab("Specificity (%)")+
  theme(text=element_text(size=16))+
  scale_fill_distiller(palette="Greens",na.value = "white",direction=1,limits=c(60, max(final2$CDR)+5))+coord_fixed()

fig4a

ggsave(file="/Users/adenooy/Library/CloudStorage/OneDrive-Personal/AMC/TB/Screening Paper/Revision 1 - LGH/figures/Screening Figure 4a.svg", plot=fig4a, width=10, height=6)


fig4b=ggplot(final2,aes(sens1,spec1,fill=`Post-screen prevalence`))+geom_tile(color = "grey50",
                                                   lwd = 0.5,
                                                   linetype = 1)+facet_grid(cols=vars(prev),rows=vars(alg_labs))+theme_bw()+xlab("Sensitivity (%)")+ylab("Specificity (%)")+
  theme(text=element_text(size=16))+
  scale_fill_distiller(palette=7,na.value = "white",direction=1,limits=c(5, 38),breaks = seq(5, 40, by = 5))+coord_fixed()

fig4b
ggsave(file="/Users/adenooy/Library/CloudStorage/OneDrive-Personal/AMC/TB/Screening Paper/Revision 1 - LGH/figures/Screening Figure4b.svg", plot=fig4b, width=10, height=6)



