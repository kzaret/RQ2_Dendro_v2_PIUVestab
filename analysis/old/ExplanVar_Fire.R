###Fire data for use in analysis of factors associated with annual tree establishment counts

##Fire scar dates from Holz 2009 (dissertation), LC:  7 trees/fire scar records

library(tidyverse)


fs <- read.csv("C:/Users/Marsh Hawk Zaret/Documents/Big_D/Data_Analysis/RQ1v2_PIUVestab/Data/FireScarDates_LagoCute_TOR_kz_20201207.csv", header = TRUE)

head(fs, 10)
colnames(fs)
sapply(fs, class)


fs %>%  ggplot(aes(x=year, y = FR.29)) + geom_point()


#count of years between fire events by various proportions of fire recording years with fire scars (e.g., fr.29 indicates >30% of trees sampled recorded fire, and those trees had been scarred once previously).  The larger the proportion, the greater the intensity and/or the greater the extent of the fire.

fs.0 <- fs %>%
  filter(FR.0 == 1) %>% #filter for event years
  mutate(tsle.0 = year - lag(year,1)) #calculate the difference between one event year and the previous event's year

fs.0 %>%  ggplot(aes(x=year, y = tsle.0)) + geom_point() + geom_line() + ylim(0, 105)


fs.19 <- fs %>%
  filter(FR.19 == 1) %>% #filter for event years
  mutate(tsle.19 = year - lag(year,1)) #calculate the difference between one event year and the previous event's year

fs.19 %>%  ggplot(aes(x=year, y = tsle.19)) + geom_point() + geom_line() + ylim(0, 105)


fs.29 <- fs %>%
  filter(FR.29 == 1) %>% #filter for event years
  mutate(tsle.29 = year - lag(year,1)) #calculate the difference between one event year and the previous event's year

fs.29 %>%  ggplot(aes(x=year, y = tsle.29)) + geom_point() + geom_line() + ylim(0, 105)


fs.39 <- fs %>%
  filter(FR.39 == 1) %>% #filter for event years
  mutate(tsle.39 = year - lag(year,1)) #calculate the difference between one event year and the previous event's year

fs.39 %>%  ggplot(aes(x=year, y = tsle.39)) + geom_point() + geom_line() + ylim(0, 105)


fs.49 <- fs %>%
  filter(FR.49 == 1) %>% #filter for event years
  mutate(tsle.49 = year - lag(year,1)) #calculate the difference between one event year and the previous event's year

fs.49 %>%  ggplot(aes(x=year, y = tsle.49)) + geom_point() + geom_line() + ylim(0, 105)


fs.59 <- fs %>%
  filter(FR.59 == 1) %>% #filter for event years
  mutate(tsle.59 = year - lag(year, 1)) #calculate the difference between one event year and the previous event's year

fs.59 %>%  ggplot(aes(x=year, y = tsle.59)) + geom_point() + geom_line() + ylim(0, 105)


fs.79 <- fs %>%
  filter(FR.79 == 1) %>% #filter for event years
  mutate(tsle.79 = year - lag(year,1)) #calculate the difference between one event year and the previous event's year

fs.79 %>%  ggplot(aes(x=year, y = tsle.79)) + geom_point() + geom_line() + ylim(0, 105)


#plot proportion of fire recording trees scarred by year
fs %>%  filter(FR_prop_a != 0) %>%
  ggplot(aes(x=year, y = FR_prop_a)) + geom_point() + geom_line()


#add sample depth to above
fs %>% filter(FR_prop_a != 0) %>%
  ggplot(aes(x=year, y = FR_prop_a)) + geom_point(color="#69b3a2") + geom_line(color="#69b3a2", size=0.5) + geom_line(aes(y=FR_cnt_a/7), color="orange") + scale_y_continuous(name="Prop. FR Trees with Fire Scars", sec.axis=sec_axis(trans=~.*7, name="Sample Depth (FR Trees)")) + scale_x_continuous(name="Year")
