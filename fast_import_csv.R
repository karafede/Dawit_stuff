library(data.table)



setwd("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/dawit Data/Hourly Database format CSV/Arranged dates/R files/filtered_4_box")



EAD_data_2013 <- read_csv("database_EAD_ 2014 _hourly_filtered.csv")[-1]
EAD_data_2014 <- read_csv("database_EAD_ 2014 _hourly_filtered.csv")[-1]
EAD_data_2015 <- read_csv("database_EAD_ 2015 _hourly_filtered.csv")[-1]
EAD_data_2016 <- read_csv("database_EAD_ 2016 _hourly_filtered.csv")[-1]

DM_data_2013 <- read_csv("database_DM_ 2013 _hourly_filtered.csv")[-1]
DM_data_2014 <- read_csv("database_DM_ 2014 _hourly_filtered.csv")[-1]
DM_data_2015 <- read_csv("database_DM_ 2015 _hourly_filtered.csv")[-1]
DM_data_2016 <- read_csv("database_DM_ 2016 _hourly_filtered.csv")[-1]

NCMS_data_2013 <- read_csv("database_NCMS_ 2013 _hourly_filtered.csv")[-1]
NCMS_data_2014 <- read_csv("database_NCMS_ 2014 _hourly_filtered.csv")[-1]
NCMS_data_2015 <- read_csv("database_NCMS_ 2015 _hourly_filtered.csv")[-1]
NCMS_data_2016 <- read_csv("database_NCMS_ 2016 _hourly_filtered.csv")[-1]


system.time(fread('../data/2008.csv', header = T, sep = ',')) 
#   user  system elapsed 
#  4.740   0.048   4.785

     
     