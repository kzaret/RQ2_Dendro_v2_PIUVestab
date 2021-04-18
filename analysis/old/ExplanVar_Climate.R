#Climate data for examining associations with PIUV estab at annual resolution

library(feasts)
library(tsibble)
library(tsibbledata)
library(tidyverse)
library(dplyr)
library(lubridate)
library(patchwork)


##====================SAM (Villalba et al. 2012)===========================
#Southern Annular Mode

SAMall<-read_csv("C:/Users/Marsh Hawk Zaret/Documents/Big_D/Data_Analysis/RQ1v2_PIUVestab/Data/Climate/SAM_Villalba2012kz.csv")

View(SAMall)


##For LC data, earliest estab. date of the tree core and interpolated data was 1646.

#LC formatting
SAMall %>% filter(age_AD > 1639) %>%
  ggplot(aes(x=age_AD)) +
  geom_line(aes(y=SAM_Mr))+ geom_line(aes(y=SAM_Mr30), color="blue", size=1) +
  geom_hline(yintercept=0, color="orange")+
  expand_limits(x=c(1600, 2020)) +
  scale_x_continuous(breaks=seq(1600, 2020, 100), minor_breaks=seq(1600, 2020, 20)) +
  xlab("Year") + ylab("Annual SAM Index (DJF) with 30-year spline") +
  theme_bw()

#gen formatting
SAMall %>% #filter(age_AD > 1639) %>%
  ggplot(aes(x=age_AD)) +
  geom_line(aes(y=SAM_Mr))+ geom_line(aes(y=SAM_Mr30), color="blue", size=1) +
  geom_hline(yintercept=0, color="orange")+
  expand_limits(x=c(1400, 2020)) +
  scale_x_continuous(breaks=seq(1400, 2020, 100), minor_breaks=seq(1400, 2020, 20)) +
  xlab("Year") +
  theme_bw()

##====================SADA (Morales et al. 2020)===========================
#South American Drought Atlas
#These data from 0.5 degree cells that surround the LC and RMB sites, respectively

sada<-read_csv("C:/Users/Marsh Hawk Zaret/Documents/Big_D/Data_Analysis/RQ1v2_PIUVestab/Data/Climate/SADA_Morales2020_LC_RMB.csv")

View(sada)

#### Tor & RMB side-by-side from B. Nanavati w/ some edits
library(tidyr)
Kyla.pdsiBySite.long <- sada %>% pivot_longer(-c(time), names_to = 'site', values_to = 'val')

View(Kyla.pdsiBySite.long)

ggplot(Kyla.pdsiBySite.long, aes(time, val, col = val))+
  geom_bar(stat = 'identity', position = 'identity')+
  geom_smooth(method = 'loess', span = 0.05, se = FALSE, col = 1)+
  scale_color_gradient2(low = "red", mid = 'white', high = "blue", midpoint = 0)+ 
  facet_wrap(~site)+
  labs(x = 'Years CA', y = 'scPDSI (DJF)')+
  theme_minimal()

#### individual sites
ggplot(sada, aes(x=time, y=RMB.sada, col=RMB.sada))+
  geom_bar(stat = 'identity', position = 'identity')+
  geom_smooth(method = 'loess', span = 0.05, se = FALSE, col = 1)+
  scale_color_gradient2(low = "red", mid = 'white', high = "blue", midpoint = 0)+ 
  labs(x = 'Years CA', y = 'scPDSI (DJF)')+
  theme_minimal()

colnames(sada)

#gen formatting
sada %>%
  ggplot(aes(x=time, y=RMB.sada)) +
  geom_line()+
  geom_smooth(method = 'loess', span = 0.05, se = FALSE, col = "purple") +
  geom_hline(yintercept=0, color="orange")+
  expand_limits(x=c(1400, 2020)) +
  scale_x_continuous(breaks=seq(1400, 2020, 100), minor_breaks=seq(1400, 2020, 20)) +
  xlab("Year") +
  theme_bw()

##====================Instrumental Precip (DGA)===========================
#Direccion General de Aguas, 2015 - 2018
#Note:  these data from instruments in the lower Palena and lower Baker rivers watersheds, respectively


### Daily
dga.day<-read_csv("C:/Users/Marsh Hawk Zaret/Documents/Big_D/Data_Analysis/RQ1v2_PIUVestab/Data/Climate/DGA_Tortel_RMB_instrumental_daily_climate.csv")

dga.day <- mutate(dga.day, Time = as.Date(Time, format= "%m/%d/%Y")) #'time' variable changed to date class
class(dga.day$Time)
tibble(dga.day)

dga.day <- as_tsibble(dga.day, key = Record_ID, index = Time) #create tsibble

#Note:  to display RMB and LC data in single plot, could use Buzz's code as for SADA dataset to pivot, then assign site as key. . . .


#### Monthly
dga.mon<-read_csv("C:/Users/Marsh Hawk Zaret/Documents/Big_D/Data_Analysis/RQ1v2_PIUVestab/Data/Climate/DGA_Tortel_RMB_instrumental_monthly_climate.csv")

colnames(dga.mon)

dga.mon <- mutate(dga.mon, Time = as.Date(Time, format= "%m/%d/%Y"))

dga.mon <- as_tsibble(dga.mon, key = Key, index = Time)

ggplot(dga.mon, aes(Time, Tor_precip_mean_mm))+
  geom_bar(stat = 'identity', position = 'identity')+
  labs(x = 'Time', y = 'Monthly Precip. (mm)')+
  expand_limits(y=c(0, 800)) +
  theme_minimal()


ggplot(dga.mon, aes(x=Time, y=Tor_precip_mean_mm)) +
  geom_line() + geom_point()+
  expand_limits(y=c(0, 800)) +
  xlab("") +
  theme_minimal()+
  scale_x_date(date_breaks="1 year", date_labels="%Y") +
  theme(axis.text.x=element_text(angle=60, hjust=1)) 

colnames(dga.mon)


#plot both sites
dga.mon<-read_csv("C:/Users/Marsh Hawk Zaret/Documents/Big_D/Data_Analysis/RQ1v2_PIUVestab/Data/Climate/DGA_Tortel_RMB_instrumental_monthly_climate.csv")

dga.mon.long <- dga.mon %>%
  select(-c("Key", "Tor_precip_RecordDays", "RMB_precip_RecordDays")) %>%
  pivot_longer(-c(Time), names_to = 'site', values_to = 'val') %>% 
  mutate(Time = as.Date(Time, format= "%m/%d/%Y")) %>%
  mutate(Key=c(1:120)) %>%
  as_tsibble(key = Key, index = Time)

View(dga.mon.long)  
  

ggplot(dga.mon.long, aes(Time, val, col = site))+
  geom_line() + geom_point()




##====================Modeled/Instrumental Climate (Camels Cr2)================
#Alvarez-Garreton et al. 2018
#Note:  these data for watersheds close to but not including the LC or RMB sites

cr2 <- read_csv("C:/Users/Marsh Hawk Zaret/Documents/Big_D/Data_Analysis/RQ1v2_PIUVestab/Data/Climate/Cr2_Chaiten_RioPascua_climatevars.csv")


