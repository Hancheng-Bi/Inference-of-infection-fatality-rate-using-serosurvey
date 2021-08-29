library(readr)
parameter <- read_csv("parameter.csv")
View(parameter)

nm <- "Credible interval of posterior of parameters"
fig_para <- parameter %>% mutate(para=as_factor(para)) %>%
  ggplot(aes(y=para, col = para)) + theme_classic() +  
  
  geom_point(aes(x=median), size=2) + 
  geom_errorbarh(aes(x=median, xmin=lower, xmax=upper), 
                 height=0, size=0.8, alpha = 0.6, show.legend=F) + 
 
  scale_x_discrete(name = nm, limits=c(-5:5),expand = c(0.1,0.1)) + 
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
fig_para