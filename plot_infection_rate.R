### read in data 
IR <- read_csv("IR.csv")
IR %>% as_tibble() %>% filter() %>% ggplot(aes(x=ses,y=IRm,color=factor(age))) +geom_point(position=position_dodge(0.1))+
  geom_errorbar(aes(ymin=IRl,ymax=IRu))+facet_wrap(~age)+geom_smooth(method=lm)+
  stat_cor(aes(label =  paste(..rr.label.., ..p.label.., sep = "~`,`~")))