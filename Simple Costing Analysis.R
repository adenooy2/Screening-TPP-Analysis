
library(tidyverse)
library(openxlsx)
rm(list=ls())

#Load data
one_high_spec=read.csv("/Users/adenooy/Downloads/costing_one_high_spec.csv")
one_high_spec$algorithm="One-screen"
one_high_spec$test="High-specificity"

one_high_sens=read.csv("/Users/adenooy/Downloads/costing_one_high_sens.csv")
one_high_sens$algorithm="One-screen"
one_high_sens$test="High-sensitivity"


two_high_spec=read.csv("/Users/adenooy/Downloads/costing_two_high_spec.csv")
two_high_spec$algorithm="Two-screen"
two_high_spec$test="High-specificity"

two_high_sens=read.csv("/Users/adenooy/Downloads/costing_two_high_sens.csv")
two_high_sens$algorithm="Two-screen"
two_high_sens$test="High-sensitivity"

cost_data=rbind(one_high_spec,one_high_sens,two_high_spec,two_high_sens)

cost_data$X=NULL
cost_data$prev=factor(cost_data$prev)
cost_data$sens1=factor(cost_data$sens1)
cost_data$spec1=factor(cost_data$spec1)


high_spec=cost_data %>% filter(test=="High-specificity")
high_spec$test=gsub("High-specificity","High-specificity test (98%)",high_spec$test)

high_sens=cost_data %>% filter(test!="High-specificity")
high_sens$test=gsub("High-sensitivity","High-sensitivity test (90%)",high_sens$test)

fig_high_spec=ggplot(high_spec,aes(sens1,prev,fill=cost_correct))+geom_tile(color = "grey50",
                                                                        lwd = 0.5,
                                                                        linetype = 1)+theme_bw()+xlab("Sensitivity (%)")+ylab("Prevalence (%)")+
  theme(text=element_text(size=16))+facet_grid(rows=vars(algorithm),cols=vars(test))+
  scale_fill_distiller(
    palette = "RdYlGn",    # reversed green-to-red, use "YlGnBu" etc. for others
    direction = -1,        # -1 to flip colors (green = low, red = high)
    limits = c(0, 10000),    # set fixed min and max
    name = "Cost per correct\n diagnosis ($)"
  )+geom_text(aes(label=round(cost_correct,0)))

fig_high_spec


fig_high_sens=ggplot(high_sens,aes(spec1,prev,fill=cost_correct))+geom_tile(color = "grey50",
                                                                            lwd = 0.5,
                                                                            linetype = 1)+theme_bw()+xlab("Specificity (%)")+ylab("Prevalence (%)")+
  theme(text=element_text(size=16))+facet_grid(rows=vars(algorithm),cols=vars(test))+
  scale_fill_distiller(
    palette = "RdYlGn",    # reversed green-to-red, use "YlGnBu" etc. for others
    direction = -1,        # -1 to flip colors (green = low, red = high)
    limits = c(0, 10000),    # set fixed min and max
    name = "Cost per correct\n diagnosis ($)"
  )+geom_text(aes(label=round(cost_correct,0)))

fig_high_sens

library("ggpubr")

figure <- ggarrange(fig_high_sens,fig_high_spec,
                              labels = c("", ""),
                              ncol = 1, nrow = 2,
                    legend = "right")
figure


ggsave(file="/Users/adenooy/Library/CloudStorage/OneDrive-Personal/AMC/TB/Screening Paper/Revision 1 - LGH/figures/Screening Figure S3_R1.svg", plot=figure, width=10, height=15)

