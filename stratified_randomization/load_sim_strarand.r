library(tidyverse)
require(gridExtra)
library(gapminder)
library(ggpubr)
theme_set(theme_pubr())

output_nms <- list.files(path="./output")
for (t in 1:length(output_nms)){
  if (t == 1){
    result <- readRDS(paste0("output/",output_nms[1]))
  } else
  {
    result <- rbind(result,readRDS(paste0("output/",output_nms[t])))
  }
}

pt <- result %>%
  as_tibble %>%
  #arrange(Rand,I,icc,theta,Method,covar_adj) %>%
  group_by(I,tau,eta,theta,covar_adj,Method,Rand) %>%
  summarise(power = mean(Reject, na.rm = TRUE))

temp <- pt %>%
  filter(theta!=0)





s1 <-  result %>%
  filter(theta!=0 & I==8 & tau==0.1 & eta==0)
s2 <-  result %>%
  filter(theta!=0 & I==8 & tau==0.1 & eta==0.1)
s3 <-  result %>%
  filter(theta!=0 & I==8 & tau==0.3 & eta==0)
s4 <-  result %>%
  filter(theta!=0 & I==8 & tau==0.3 & eta==0.1)
s_combined <- rbind(s1,s2,s3,s4)
ymin = min(s_combined$Estimate)
ymax = max(s_combined$Estimate)

plot_s1 <-  s1 %>% 
  ggplot(aes(x=Method,y=Estimate, fill=factor(covar_adj))) +
  geom_boxplot() +
  facet_wrap(~Rand) + 
  ggtitle("Scenario 1") +
  ylim(ymin,ymax) + 
  geom_hline(aes(yintercept=0.25),colour="#990000", linetype="dashed")

plot_s2 <-  s2 %>% 
  ggplot(aes(x=Method,y=Estimate, fill=factor(covar_adj))) +
  geom_boxplot() +
  facet_wrap(~Rand) + 
  ggtitle("Scenario 2") +
  ylim(ymin,ymax) + 
  geom_hline(aes(yintercept=0.25),colour="#990000", linetype="dashed")

plot_s3 <-  s3 %>% 
  ggplot(aes(x=Method,y=Estimate, fill=factor(covar_adj))) +
  geom_boxplot() +
  facet_wrap(~Rand) + 
  ggtitle("Scenario 3") +
  ylim(ymin,ymax) + 
  geom_hline(aes(yintercept=0.25),colour="#990000", linetype="dashed")

plot_s4 <-  s4 %>% 
  ggplot(aes(x=Method,y=Estimate, fill=factor(covar_adj))) +
  geom_boxplot() +
  facet_wrap(~Rand) + 
  ggtitle("Scenario 4") +
  ylim(ymin,ymax) + 
  geom_hline(aes(yintercept=0.25),colour="#990000", linetype="dashed")


figure <- ggarrange(plot_s1,plot_s2,plot_s3,plot_s4,
                ncol = 2, nrow = 2,
                common.legend = TRUE)
figure  
  





