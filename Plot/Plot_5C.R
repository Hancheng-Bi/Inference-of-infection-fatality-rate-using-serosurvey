library(sf)
library(lemon)
library(ggpubr)
library(ggrepel)
library(cowplot)
library(ggspatial)
library(lubridate)
library(tidyverse)
library(readr)
IFR_rate <- read_csv("IFR_rate.csv")
View(IFR_rate)
colsF5 <- c("royalblue4","mediumpurple4","maroon","lightcoral")
brksF5 <- c('1','2','3','4')
legend <- c("0-40" = "royalblue4" ,"40-60"="mediumpurple4","60-80"="maroon","80+"="lightcoral")

nm <- "Ratio between the IFR values of\nthe low and high SES categories"
fig5C <- IFR_rate %>% mutate(age=as_factor(age)) %>%
  ggplot(aes(y=age, col = age)) + theme_classic() +  
  geom_vline(xintercept=1, lty="dashed", col = "gray70") + 
  geom_point(aes(x=ratem), size=2) + 
  geom_errorbarh(aes(x=ratem, xmin=ratel, xmax=rateu), 
                 height=0, size=0.8, alpha = 0.6, show.legend=F) + 
  scale_color_manual('Age group', values=colsF5, labels=brksF5,
                     guide=guide_legend(override.aes=list(size=1.5), reverse=TRUE)) + 
  scale_y_discrete('Age group', 
                   labels = c("1 and 2"="0-60")) +
  scale_x_discrete(name = nm, limits=c(0:12),expand = c(0.1,0.1)) + 
  theme(legend.justification = c(1, 1), 
        legend.position = c(1, 1),
        legend.key.height=unit(0.3,"cm"),
        legend.key.width=unit(0.3,"cm"),
        legend.title = element_text(size=9),
        legend.text = element_text(size=7),
        legend.background = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        axis.line = element_line(size = 0.4),
        axis.text = element_text(size=8, colour= "black"),
        axis.title = element_text(size=9, colour= "black")) 
fig5C
