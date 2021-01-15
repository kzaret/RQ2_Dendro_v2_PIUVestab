#Climate data for examining associations with PIUV estab at annual resolution

library(dplyr)
library(tidyverse)
library(patchwork)


##====================SAM (Villalba et al. 2012)===========================
#Southern Annular Mode

SAMall<-read_csv("C:/Users/Marsh Hawk Zaret/Documents/Big_D/Data_Analysis/RQ1v2_PIUVestab/Data/SAM_Villalba2012kz.csv")

View(SAMall)


##or LC data, earliest estab. date of the tree core and interpolated data was 1646.

#LC formatting
SAMall %>% filter(age_AD > 1639) %>%
  ggplot(aes(x=age_AD)) +
  geom_line(aes(y=SAM_Mr))+ geom_line(aes(y=SAM_Mr30), color="blue", size=1) +
  geom_hline(yintercept=0, color="orange")+
  expand_limits(x=c(1600, 2020)) +
  scale_x_continuous(breaks=seq(1600, 2020, 100), minor_breaks=seq(1600, 2020, 20)) +
  xlab("Year") + ylab("Annual SAM Index (DJF) with 30-year spline") +
  theme_bw()



##====================SADA (Morales et al. 2020)===========================
#South American Drought Atlas

sada<-read_csv("C:/Users/Marsh Hawk Zaret/Documents/Big_D/Data_Analysis/RQ1v2_PIUVestab/Data/SADA_Morales2020_LC_RMB.csv")

View(sada)

