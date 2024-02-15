library(reshape2)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggridges)
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

hdmi.expo.plot <- ggplot()+
  geom_density_ridges(data = hdmi.expo_bind, aes(x=value, y=reorder(Chemical,order), color=name, fill=name), 
                      scale=0.8, alpha=0.4)+
  geom_point(data=hdmi.expo_point, aes(x=quantile, y=Chemical, color=name, shape=name), size=2)+
  scale_x_log10(limits=c(1e-013, 1e+02),
                breaks = scales::trans_breaks("log10", function(x) 10^x, n=4),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme_bw()+
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 12), 
        axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 12), 
        legend.text = element_text(size = 12), legend.title = element_blank(), legend.position = "top", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(5.5, 5.5, 5.5, 0),"pt"))+
  xlab(bquote(HD[M]^I~or~ExpoCast~(mg/kg/day)))+
  scale_color_manual(name="legend",
                     breaks=c("Chemical-specific","ExpoCast"), 
                     labels=c(bquote(Chemical-specific~HD[M]^I), "ExpoCast"),
                     values=c(`Chemical-specific`="#D55E00",`ExpoCast`="#009E73"))+
  scale_fill_manual(name="legend",
                    breaks=c("Chemical-specific","ExpoCast"), 
                    labels=c(bquote(Chemical-specific~HD[M]^I), "ExpoCast"),
                    values=c(`Chemical-specific`="#D55E00",`ExpoCast`="#009E73"))+
  scale_shape_manual(name="legend",
                     breaks=c("Chemical-specific","ExpoCast"), 
                     labels=c(bquote(Chemical-specific~HD[M]^I), "ExpoCast"),
                     values=c(`Chemical-specific`=16,`ExpoCast`=18))+
  annotation_logticks(sides="b")+
  guides(color = guide_legend(nrow = 2), fill = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))
print(hdmi.expo.plot)

#BEMI blood
bemi_b.bio <- left_join(be.blood.df, bio.data[1:4,2:3], by="CAS")
bemi_b.bio$`5th_quantile` <- apply(bemi_b.bio[,3:10002], 1, quantile, probs=0.05)
bemi_b.bio.df <- bemi_b.bio %>% group_by(CAS) %>% slice(which.min(`5th_quantile`))
bemi_b.bio.df <- bemi_b.bio.df %>% arrange(CAS)
bemi_b.m <- melt(bemi_b.bio.df, id=c(1,2,10003:10005))
bemi_b.m$order <- hdmi.expo_bind$order[1:190000]

be_b.5th <- bemi_b.bio.df[,c(1,10003,10005)]
be_b.5th$name <- "Chemical-specific BEMI"
colnames(be_b.5th)[3] <- "value"
blood.data <- bemi_b.bio.df[,c(1,10003,10004)]
blood.data$name <- "Biomonitoring data in blood"
colnames(blood.data)[3] <- "value"
bemi_b.point <- rbind(be_b.5th, blood.data)

bemi_b.plot <- ggplot()+
  geom_density_ridges(data = bemi_b.m, aes(x=value, y=reorder(Chemical, order)), color="#D55E00", fill="#D55E00", 
                      scale=0.7, alpha=0.4)+
  geom_point(data=bemi_b.point, aes(x=value, y=Chemical, color=name), size=2)+
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
  scale_color_manual(breaks=c("Chemical-specific BEMI","Biomonitoring data in blood"), 
                     labels=c(bquote(Chemical-specific~BE[M]^I~"in"~blood),"Biomonitoring data in blood"),
                     values=c(`Chemical-specific BEMI`="#D55E00",`Biomonitoring data in blood`="#56B4E9"))+
  annotation_logticks(sides="b")+
  guides(color = guide_legend(nrow = 2), fill = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))
print(bemi_b.plot)

#BEMI urine
bemi_u.bio <- left_join(be.urine.df, bio.data[5:6,2:3], by="CAS")
bemi_u.bio$`5th_quantile` <- apply(bemi_u.bio[,3:10002], 1, quantile, probs=0.05)
bemi_u.bio.df <- bemi_u.bio %>% group_by(CAS) %>% slice(which.min(`5th_quantile`))
bemi_u.bio.df <- bemi_u.bio.df %>% arrange(CAS)
bemi_u.m <- melt(bemi_u.bio.df, id=c(1,2,10003:10005))
bemi_u.m$order <- hdmi.expo_bind$order[1:190000]

be_u.5th <- bemi_u.bio.df[,c(1,10003,10005)]
be_u.5th$name <- "Chemical-specific BEMI"
colnames(be_u.5th)[3] <- "value"
urine.data <- bemi_u.bio.df[,c(1,10003,10004)]
urine.data$name <- "Biomonitoring data in urine"
colnames(urine.data)[3] <- "value"
bemi_u.point <- rbind(be_u.5th, urine.data)

bemi_u.plot <- ggplot()+
  geom_density_ridges(data = bemi_u.m, aes(x=value, y=reorder(Chemical, order)), color="#D55E00", fill="#D55E00", 
                      scale=0.8, alpha=0.4)+
  geom_point(data=bemi_u.point, aes(x=value, y=Chemical, color=name), size=2)+
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
  scale_color_manual(breaks=c("Chemical-specific BEMI","Biomonitoring data in urine"), 
                     labels=c(bquote(Chemical-specific~BE[M]^I~"in"~urine),"Biomonitoring data in urine"),
                     values=c(`Chemical-specific BEMI`="#D55E00",`Biomonitoring data in urine`="#56B4E9"))+
  annotation_logticks(sides="b")+
  guides(color = guide_legend(nrow = 2), fill = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))

combine <- ggarrange(hdmi.expo.plot, bemi_b.plot, bemi_u.plot, labels = c("A","B","C"), ncol=3, 
                     widths = c(1.2/3,0.9/3,0.9/3))
ggsave(combine, file="the most sensitive HDMI and BEMI in blood and urine ggridges plot.pdf", 
       width=25, height=14, path="BEMI plots", scale=0.6)

####-------------------------------------------------
##all endpoints
#HDMI and Expocast
cs.hdmi$quantile <- apply(cs.hdmi[,3:10002], 1, quantile, probs=0.05) #5th quantile
cs.hdmi.chem <- merge(chem, cs.hdmi, by="CAS") %>% arrange(endpoint)
cs.hdmi.all.m <- melt(cs.hdmi.chem, id=c(1:3,10004))
cs.hdmi.all.m$name <- "Chemical-specific"

expo.df <- expo.dis %>% arrange(CAS)
expo.df$quantile <- apply(expo.df[,3:10002], 1, quantile, probs=0.95) #95th quantile
expo.endpoint <- merge(cs.hdmi[,c(1,2)],expo.df, by="CAS") %>% arrange(endpoint)
expo.endpoint.m <- melt(expo.endpoint, id=c(1:3,10004))
expo.endpoint.m$name <- "ExpoCast"

hdmi.expo.all_bind <- rbind(cs.hdmi.all.m, expo.endpoint.m)
hdmi.expo.all_bind$order <- hdmi.expo.all_bind$quantile
hdmi.expo.all_bind[350001:700000,8] <- hdmi.expo.all_bind[1:350000,8]

hdmi.expo.all_point <- rbind(cs.hdmi[,c(1,2,10003)], expo.endpoint[,c(1,2,10004)])
hdmi.expo.all_point$name <- "Chemical-specific"
hdmi.expo.all_point$name[36:70] <- "ExpoCast"

hdmi.expo.all.plot <- ggplot()+
  geom_density_ridges(data = hdmi.expo.all_bind, aes(x=value, y=reorder(endpoint,order), color=name, fill=name), 
                      scale=0.8, alpha=0.4)+
  geom_point(data=hdmi.expo.all_point, aes(x=quantile, y=endpoint, color=name, shape=name), size=2)+
  scale_x_log10(limits=c(1e-014, 1e+02),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme_bw()+
  theme(axis.title.y = element_blank(), 
        axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 12), 
        legend.text = element_text(size = 12), legend.title = element_blank(), legend.position = "top", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(5.5, 5.5, 5.5, 0),"pt"))+
  xlab(bquote(HD[M]^I~or~ExpoCast~(mg/kg/day)))+
  scale_color_manual(name="legend",
                     breaks=c("Chemical-specific","ExpoCast"), 
                     labels=c(bquote(Chemical-specific~HD[M]^I), "ExpoCast"),
                     values=c(`Chemical-specific`="#D55E00",`ExpoCast`="#009E73"))+
  scale_fill_manual(name="legend",
                    breaks=c("Chemical-specific","ExpoCast"), 
                    labels=c(bquote(Chemical-specific~HD[M]^I), "ExpoCast"),
                    values=c(`Chemical-specific`="#D55E00",`ExpoCast`="#009E73"))+
  scale_shape_manual(name="legend",
                     breaks=c("Chemical-specific","ExpoCast"), 
                     labels=c(bquote(Chemical-specific~HD[M]^I), "ExpoCast"),
                     values=c(`Chemical-specific`=16,`ExpoCast`=18))+
  annotation_logticks(sides="b")+
  guides(color = guide_legend(nrow = 2), fill = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))

#BEMI blood
bemi_b.bio.all <- left_join(be.blood, bio.data[1:4,2:3], by="CAS")
bemi_b.bio.all$`5th_quantile` <- apply(bemi_b.bio.all[,3:10002], 1, quantile, probs=0.05)
bemi_b.bio.all <- bemi_b.bio.all %>% arrange(endpoint)
bemi_b.all.m <- melt(bemi_b.bio.all, id=c(1,2,10003,10004))
bemi_b.all.m$order <- cs.hdmi.all.m$quantile

bemi_b.all.point <- bemi_b.bio.all[,c(1,2,10003,10004)]

be_b.all.5th <- bemi_b.bio.all[,c(1,2,10004)]
be_b.all.5th$name <- "Chemical-specific BEMI"
colnames(be_b.all.5th)[3] <- "value"
blood.all.data <- bemi_b.bio.all[,c(1,2,10003)]
blood.all.data$name <- "Biomonitoring data in blood"
colnames(blood.all.data)[3] <- "value"
bemi_b.all.point <- rbind(be_b.all.5th, blood.all.data)

bemi_b.all.plot <- ggplot()+
  geom_density_ridges(data = bemi_b.all.m, aes(x=value, y=reorder(endpoint, order)), color="#D55E00", fill="#D55E00", 
                      scale=0.7, alpha=0.4)+
  geom_point(data=bemi_b.all.point, aes(x=value, y=endpoint, color=name), size=2)+
  scale_x_log10(limits=c(1e-06, 1e+03),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme_bw()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), 
        axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 12), 
        legend.title = element_blank(), legend.text = element_text(size = 12), legend.position = "top", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(5.5, 5.5, 5.5, 0),"pt"))+
  xlab(bquote(BE[M]^I~`in`~blood~(mg/L)))+
  scale_color_manual(breaks=c("Chemical-specific BEMI","Biomonitoring data in blood"), 
                     labels=c(bquote(Chemical-specific~BE[M]^I~"in"~blood),"Biomonitoring data in blood"),
                     values=c(`Chemical-specific BEMI`="#D55E00",`Biomonitoring data in blood`="#56B4E9"))+
  annotation_logticks(sides="b")+
  guides(color = guide_legend(nrow = 2), fill = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))

#BEMI urine
bemi_u.bio.all <- left_join(be.urine, bio.data[5:6,2:3], by="CAS")
bemi_u.bio.all$`5th_quantile` <- apply(bemi_u.bio.all[,3:10002], 1, quantile, probs=0.05)
bemi_u.bio.all <- bemi_u.bio.all %>% arrange(endpoint)
bemi_u.all.m <- melt(bemi_u.bio.all, id=c(1,2,10003,10004))
bemi_u.all.m$order <- cs.hdmi.all.m$quantile

bemi_u.all.point <- bemi_u.bio.all[,c(1,2,10003:10004)]

be_u.all.5th <- bemi_u.bio.all[,c(1,2,10004)]
be_u.all.5th$name <- "Chemical-specific BEMI"
colnames(be_u.all.5th)[3] <- "value"
urine.all.data <- bemi_b.bio.all[,c(1,2,10003)]
urine.all.data$name <- "Biomonitoring data in urine"
colnames(urine.all.data)[3] <- "value"
bemi_u.all.point <- rbind(be_u.all.5th, urine.all.data)

bemi_u.all.plot <- ggplot()+
  geom_density_ridges(data = bemi_u.all.m, aes(x=value, y=reorder(endpoint, order)), color="#D55E00", fill="#D55E00", 
                      scale=0.7, alpha=0.4)+
  geom_point(data=bemi_u.all.point, aes(x=value, y=endpoint, color=name), size=2)+
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
  scale_color_manual(breaks=c("Chemical-specific BEMI","Biomonitoring data in urine"), 
                     labels=c(bquote(Chemical-specific~BE[M]^I~"in"~urine),"Biomonitoring data in urine"),
                     values=c(`Chemical-specific BEMI`="#D55E00",`Biomonitoring data in urine`="#56B4E9"))+
  annotation_logticks(sides="b")+
  guides(color = guide_legend(nrow = 2), fill = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))

combine.all <- ggarrange(hdmi.expo.all.plot, bemi_b.all.plot, bemi_u.all.plot, labels = c("A","B","C"), 
                         ncol=3, widths = c(1.4/3,0.8/3,0.8/3))
ggsave(combine.all, file="all endpoints HDMI and BEMI in blood and urine ggridges plot.pdf", 
       width=25, height=14, path="BEMI plots", scale=0.6)
