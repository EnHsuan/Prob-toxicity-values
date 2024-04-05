library(reshape2)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)

bio.data <- read.csv("bimonitoring data.csv")
be.blood <- read.csv("BEMI blood.csv")
be.urine <- read.csv("BEMI urine.csv")
chem <- read.csv("regulatory reference dose.csv")[,1:2]

cs.hdmi <- read.csv("HDMI.csv")
expo.dis <- read.csv("Expocast distribution by row.csv")

be.blood.df <- merge(be.blood, chem, by="CAS")
be.urine.df <- merge(be.urine, chem, by="CAS")

#HDMI and Expocast
cs.hdmi$quantile <- apply(cs.hdmi[,3:10002], 1, quantile, probs=0.05) #5th quantile
cs.hdmi.sen <- cs.hdmi %>% group_by(CAS) %>% slice(which.min(quantile))
cs.hdmi.df <- merge(chem, cs.hdmi.sen, by="CAS")
cs.hdmi.m <- melt(cs.hdmi.df, id=c(1:3,10004))
cs.hdmi.m$name <- "Chemical-specific"

expo.df <- expo.dis %>% arrange(CAS)
expo.df$quantile <- apply(expo.df[,3:10002], 1, quantile, probs=0.95) #95th quantile
expo.m <- melt(expo.df, id=c(1,2,10003))
expo.m$name <- "ExpoCast"

hdmi.expo_bind <- rbind(cs.hdmi.m[,c(1,2,4:7)], expo.m)
hdmi.expo_bind$order <- hdmi.expo_bind$quantile
hdmi.expo_bind[190001:380000,7] <- hdmi.expo_bind[1:190000,7]

hdmi.expo_point <- rbind(cs.hdmi.df[,c(1,2,10004)], expo.df[,c(1,2,10003)])
hdmi.expo_point$name <- "Chemical-specific"
hdmi.expo_point$name[20:38] <- "ExpoCast"

#set whisker to [5%,95%]
quantiles_95 <- function(x) {
  r <- quantile(x, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

hdmi.expo.boxplot <- ggplot()+
  stat_summary(fun.data = quantiles_95, data = hdmi.expo_bind, geom="boxplot", aes(x=value, y=reorder(Chemical,order), fill=name), outlier.shape=NA, position = position_dodge(width = 0.75))+
  scale_x_log10(limits=c(1e-013, 1e+02),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme_bw()+
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 12), 
        axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 12), 
        legend.text = element_text(size = 12), legend.title = element_blank(), legend.position = "top", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(5.5, 5.5, 5.5, 0),"pt"))+
  xlab(bquote(HD[M]^I~or~ExpoCast~(mg/kg/day)))+
  scale_fill_manual(#name="legend",
                    breaks=c("Chemical-specific","ExpoCast"), 
                    labels=c(bquote(Chemical-specific~HD[M]^I), "ExpoCast"),
                    values=c(`Chemical-specific`="red",`ExpoCast`="white"))+
  annotation_logticks(sides="b")
print(hdmi.expo.boxplot)

#BEMI blood
bemi_b.bio <- left_join(be.blood.df, bio.data[1:4,2:3], by="CAS")
bemi_b.bio$`5th_quantile` <- apply(bemi_b.bio[,3:10002], 1, quantile, probs=0.05)
bemi_b.bio.df <- bemi_b.bio %>% group_by(CAS) %>% slice(which.min(`5th_quantile`))
bemi_b.bio.df <- bemi_b.bio.df %>% arrange(CAS)
bemi_b.m <- melt(bemi_b.bio.df, id=c(1,2,10003:10005))
bemi_b.m$order <- hdmi.expo_bind$order[1:190000]
bemi_b.m$name <- "Chemical-specific"

be_b.5th <- bemi_b.bio.df[,c(1,10003,10005)]
be_b.5th$name <- "Chemical-specific BEMI"
colnames(be_b.5th)[3] <- "value"
blood.data <- bemi_b.bio.df[,c(1,10003,10004)]
blood.data$name <- "Biomonitoring data in blood"
colnames(blood.data)[3] <- "value"
bemi_b.point <- rbind(be_b.5th, blood.data)

bemi_b.boxplot <- ggplot()+
  stat_summary(fun.data = quantiles_95, data = bemi_b.m, geom="boxplot", aes(x=value, y=reorder(Chemical,order), fill=name), outlier.shape=NA, position = position_dodge(width = 0.75))+
  geom_point(data=blood.data, aes(x=value, y=Chemical, color=name), shape=15, size=4)+
  scale_x_log10(limits=c(1e-07, 1e+03),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme_bw()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), 
        axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 12), 
        legend.title = element_blank(), legend.text = element_text(size = 12), legend.position = "top", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(5.5, 5.5, 5.5, 0),"pt"))+
  xlab(bquote(BE[M]^I~`in`~blood~(mg/L)))+
  scale_fill_manual(name="Biomonitoring data in blood", 
                    breaks="Chemical-specific", 
                    labels=c(bquote(Chemical-specific~BBE[M]^I)),
                    values=c(`Chemical-specific`="red"))+
  scale_color_manual(name="Chemical-specific", 
                     breaks="Biomonitoring data in blood", 
                     labels=c("Biomonitoring data in blood"="Blood data"),
                     values=c(`Biomonitoring data in blood`="#56B4E9"))+
  annotation_logticks(sides="b")
print(bemi_b.boxplot)

#BEMI urine
bemi_u.bio <- left_join(be.urine.df, bio.data[5:6,2:3], by="CAS")
bemi_u.bio$`5th_quantile` <- apply(bemi_u.bio[,3:10002], 1, quantile, probs=0.05)
bemi_u.bio.df <- bemi_u.bio %>% group_by(CAS) %>% slice(which.min(`5th_quantile`))
bemi_u.bio.df <- bemi_u.bio.df %>% arrange(CAS)
bemi_u.m <- melt(bemi_u.bio.df, id=c(1,2,10003:10005))
bemi_u.m$order <- hdmi.expo_bind$order[1:190000]
bemi_u.m$name <- "Chemical-specific"

be_u.5th <- bemi_u.bio.df[,c(1,10003,10005)]
be_u.5th$name <- "Chemical-specific BEMI"
colnames(be_u.5th)[3] <- "value"
urine.data <- bemi_u.bio.df[,c(1,10003,10004)]
urine.data$name <- "Biomonitoring data in urine"
colnames(urine.data)[3] <- "value"
bemi_u.point <- rbind(be_u.5th, urine.data)

bemi_u.boxplot <- ggplot()+
  stat_summary(fun.data = quantiles_95, data = bemi_u.m, geom="boxplot", aes(x=value, y=reorder(Chemical,order), fill=name), outlier.shape=NA, position = position_dodge(width = 0.75))+
  geom_point(data=urine.data, aes(x=value, y=Chemical, color=name), shape=17, size=4)+
  scale_x_log10(limits=c(1e-09, 1),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme_bw()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), 
        axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 12), 
        legend.title = element_blank(), legend.text = element_text(size = 12), legend.position = "top", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(5.5, 5.5, 5.5, 0),"pt"))+
  xlab(bquote(BE[M]^I~`in`~urine~(mg/kg/day)))+
  scale_fill_manual(name="Biomonitoring data in urine",
                    breaks="Chemical-specific", 
                    labels=c(bquote(Chemical-specific~UBE[M]^I)),
                    values=c(`Chemical-specific`="red"))+
  scale_color_manual(name="Chemical-specific",
                     breaks="Biomonitoring data in urine", 
                     labels=c("Biomonitoring data in urine"="Urinary data"),
                     values=c(`Biomonitoring data in urine`="#56B4E9"))+
  annotation_logticks(sides="b")
print(bemi_u.boxplot)

boxplot <- ggarrange(hdmi.expo.boxplot, bemi_b.boxplot, bemi_u.boxplot, ncol=3, widths = c(1.2/3,0.9/3,0.9/3))
ggsave(boxplot, file="Chemical-specific toxicity values with exposure data boxplot.pdf", 
       width=25, height=14, path="Plot output", scale=0.6)
