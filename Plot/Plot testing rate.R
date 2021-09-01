library(lemon)
library(ggpubr)
library(ggrepel)
library(cowplot)
library(ggspatial)
library(lubridate)
library(readr)
### read in data 
data <- read_csv("data_combined.csv")
View(data)

ggplot(data,aes(x=ses,y=Age_group))+
  geom_point(aes(size=100*Sample_size/Pop)) +
  scale_size_area()
