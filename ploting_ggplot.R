
jpeg('D:/R_processing/plots/Abu_Dhabi_Asthsma_respiratory.jpg',   
     quality = 100, bg = "white", res = 200, width = 13, height = 7, units = "in")
par(mar=c(4, 10, 9, 2) + 0.3)
oldpar <- par(las=1)


min <- as.Date("2011-06-01") 
max <- as.Date("2013-05-31") 


# count patients in each day (all genders all ages)
health_data_sum <- health_data %>%
  group_by(Date) %>%
  summarise(sum_patients = sum(sum_patients, na.rm = TRUE))

plot <- ggplot(health_data_sum, aes(Date, sum_patients)) + 
  theme_bw() +
  geom_line(aes(y = sum_patients, col = "sum_patients"), alpha=1, col="blue") +
  ggtitle("counts admissions (Abu Dhabi - Asthma, 2013-2015)") + 
  theme(plot.title = element_text(lineheight=.8, face="bold", size = 20, hjust = 0.5)) +
  stat_smooth(method = "loess") +
  theme(legend.position="none") + 
  ylab(expression("Sum Patients (counts) per day")) + 
  theme(axis.title.x=element_blank(),
        axis.text.x  = element_text(angle=0, vjust=0.5, hjust = 0.5, size=22, colour = "black", face="bold")) +
  theme(axis.title.y = element_text(face="bold", colour="black", size=22),
        axis.text.y  = element_text(angle=0, vjust=0.5, size=20, colour = "black")) +
  ylim(0, 600) + 
  xlim(min, max) 
plot

par(oldpar)
dev.off()

dawit<-ChickWeight

# First plot
library(ggplot2)
p1<- qplot( dawi$localR2, binwidth= 0.01, xlim=c(0,1),
           geom="histogram",xlab = paste0(month.name[01]), fill=I("blue" ),alpha=I(.4),ylab ="")

p1 <- ggplot(ChickWeight, aes(x=Time, y=weight, colour=Diet, group=Chick)) +
  geom_line() +
  ggtitle("Growth curve for individual chicks")


# Second plot
p2 <- ggplot(ChickWeight, aes(x=Time, y=weight, colour=Diet)) +
  geom_point(alpha=.3) +
  geom_smooth(alpha=.2, size=1) +
  ggtitle("Fitted growth curve per diet")

# Third plot
p3 <- ggplot(subset(ChickWeight, Time==21), aes(x=weight, colour=Diet)) +
  geom_density() +
  ggtitle("Final weight, by diet")

# Fourth plot
p4 <- ggplot(subset(ChickWeight, Time==21), aes(x=weight, fill=Diet)) +
  geom_histogram(colour="black", binwidth=50) +
  facet_grid(Diet ~ .) +
  ggtitle("Final weight, by diet") +
  theme(legend.position="none") 
