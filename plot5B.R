library(sf)
library(lemon)
library(ggpubr)
library(ggrepel)
library(cowplot)
library(ggspatial)
library(lubridate)
library(tidyverse)

### B. Estimated IFR per age
colsF5 <- c("royalblue4","mediumpurple4","maroon","lightcoral")
brksF5 <- c('1','2','3','4')
legend <- c("0-40" = "royalblue4" ,"40-60"="mediumpurple4","60-80"="maroon","80+"="lightcoral")


IFR <- read_csv("IFR_nonprefer_0.001.csv")
IFR$Age <-as.factor(IFR$Age)

fig5B <- IFR %>% 
  ggplot(aes(x=Ses,y=IFRm,color=Age,fill=Age)) +
  geom_smooth(method = "lm", formula = y~x, size = 0.5, 
              alpha = 0.3, show.legend = FALSE) +
  geom_linerange(aes(ymin=IFRl,ymax=IFRu), size=0.3, 
                 alpha = 0.5, show.legend=F) +
  geom_point(aes(color=Age), size=1, alpha = 0.6) + 
  scale_y_continuous("Infection fatality rate (IFR)", trans = "log10",
                     labels=scales::label_percent(accuracy = 0.01)) + 
  
  scale_color_manual(breaks= brksF5, values=colsF5) +
  scale_fill_manual(breaks=brksF5, values=colsF5) +
  
  theme_classic() + scale_x_continuous("Socioeconomic status\n",
                                       breaks = seq(0,100,20), limits = c(16,94)) + 
  
  labs(color='Age group',fill='Age group') + 

  theme(legend.position = "none",
        legend.key.height=unit(0.4,"cm"),
        legend.key.width=unit(0.4,"cm"),
        legend.title = element_text(size=8),
        legend.text = element_text(size=7),
        legend.background = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        axis.line = element_line(size = 0.4),
        axis.text = element_text(size=8, colour= "black"),
        axis.title = element_text(size=9, colour= "black")) 

fig5B                  
