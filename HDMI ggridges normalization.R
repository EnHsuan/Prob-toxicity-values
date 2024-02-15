library(reshape2)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggridges)
library(ggpubr)

rfd <- read.csv("regulatory reference dose.csv")[,c(1,2,5,7,8)]

bmd <- read.csv("bmd samples CAS by row.csv")
who.bmd.df <- read.csv("WHO bmd samples CAS by row.csv")

hdmi.df <- read.csv("HDMI.csv")
who.hdmi.df <- read.csv("WHO approximate HDMI.csv")

###the most sensitive endpoint
###BMD
#normalized by regulatory POD
bmd.pod <- left_join(bmd, rfd[,c(1:3)], by="CAS")
bmd.nor <- with(bmd.pod, bmd.pod[,1:10000]/bmd.pod[,10004])
bmd.nor.df <- cbind(bmd.pod[,10001:10003], bmd.nor)
bmd.nor.end <- left_join(bmd.nor.df, who.bmd.df[,10001:10002], by="CAS")
bmd.nor.end$end <- paste("CS =", bmd.nor.end$endpoint.x, ", WHO =", bmd.nor.end$endpoint.y)
bmd.nor.end.df <- bmd.nor.end[-c(3,4,21,22,24,27,29),]
bmd.nor.end.df$`5th_quantile` <- apply(bmd.nor.end.df[,4:10003], 1, quantile, probs=0.05)
bmd.nor.sen.end.df <- bmd.nor.end.df %>% group_by(CAS) %>% slice(which.min(`5th_quantile`))
bmd.nor.sen.m <- melt(bmd.nor.sen.end.df, id=c(1:3,10004:10006))
bmd.nor.sen.m$name <- "Bayesian BMD or WHO/IPCS (2018) BMD"

who.bmd.pod <- left_join(who.bmd.df, rfd[,c(1:3)], by="CAS")
who.bmd.nor <- with(who.bmd.pod, who.bmd.pod[,1:10000]/who.bmd.pod[,10004])
who.bmd.nor.df <- cbind(who.bmd.pod[,10001:10003], who.bmd.nor)
who.bmd.nor.end <- left_join(bmd[,10001:10002], who.bmd.nor.df, by="CAS")
who.bmd.nor.end$end <- paste("CS =", who.bmd.nor.end$endpoint.x, ", WHO =", who.bmd.nor.end$endpoint.y)
who.bmd.nor.end.df <- who.bmd.nor.end[-c(3,4,21,22,24,27,29),]
who.bmd.nor.end.df$`5th_quantile` <- apply(who.bmd.nor.end.df[,5:10004], 1, quantile, probs=0.05)
who.bmd.nor.sen.end.df <- who.bmd.nor.end.df %>% group_by(CAS) %>% slice(which.min(`5th_quantile`))
who.bmd.nor.sen.m <- melt(who.bmd.nor.sen.end.df, id=c(1:4,10005:10006))
who.bmd.nor.sen.m$name <- "WHO/IPCS (2018) BMD"

bmd.sen.bind <- rbind(bmd.nor.sen.m, who.bmd.nor.sen.m)
bmd.sen.bind$order <- bmd.sen.bind$`5th_quantile`
bmd.sen.bind[1:190000,10] <- bmd.sen.bind[190001:380000,10]
bmd.sen.bind$name <- factor(bmd.sen.bind$name, levels = c("WHO/IPCS (2018) BMD","Bayesian BMD or WHO/IPCS (2018) BMD"))

bmd.point <- rbind(bmd.nor.sen.end.df[,c(1,3,10006)],who.bmd.nor.sen.end.df[,c(1,4,10006)])
bmd.point$name <- "Bayesian BMD or WHO/IPCS (2018) BMD"
bmd.point$name[20:38] <- "WHO/IPCS (2018) BMD"
bmd.point$name <- factor(bmd.point$name, levels = c("WHO/IPCS (2018) BMD","Bayesian BMD or WHO/IPCS (2018) BMD"))

bmd.sen.plot <- ggplot()+
  geom_density_ridges(data = bmd.sen.bind, aes(x=value, y=reorder(Chemical, order), color=name, fill=name), 
                      scale=0.8, alpha=0.6)+
  geom_point(data=bmd.point, aes(x=`5th_quantile`, y=Chemical, color=name))+
  geom_vline(xintercept = 1, color="black", linetype="dashed")+
  scale_x_log10(limits=c(1e-02, 1e+03),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme_bw()+
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 10), 
        legend.text = element_text(size = 10), legend.title = element_blank(),
        legend.position = "top", legend.key.size = unit(3, 'mm'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab(bquote(BMD[M]^I/Regulatory~POD))+
  scale_color_manual(breaks=c("WHO/IPCS (2018) BMD","Bayesian BMD or WHO/IPCS (2018) BMD"),
                     values=c(`WHO/IPCS (2018) BMD`="#56B4E9", `Bayesian BMD or WHO/IPCS (2018) BMD`="#D55E00"),
                     labels=c("WHO/IPCS (2018) BMD"="WHO/IPCS BMD","Bayesian BMD or WHO/IPCS (2018) BMD"="BBMD or WHO/IPCS BMD"))+
  scale_fill_manual(breaks=c("WHO/IPCS (2018) BMD","Bayesian BMD or WHO/IPCS (2018) BMD"),
                    values=c(`WHO/IPCS (2018) BMD`="#56B4E9", `Bayesian BMD or WHO/IPCS (2018) BMD`="#D55E00"),
                    labels=c("WHO/IPCS (2018) BMD"="WHO/IPCS BMD","Bayesian BMD or WHO/IPCS (2018) BMD"="BBMD or WHO/IPCS BMD"))+
  annotation_logticks(sides="b")+
  guides(color = guide_legend(nrow = 2), fill = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))

###HDMI
#normalized by regulatory RfD
hdmi.rfd <- left_join(hdmi.df, rfd[,c(1,2,4)], by="CAS")
hdmi.nor <- with(hdmi.rfd, hdmi.rfd[,3:10002]/hdmi.rfd[,10004])
hdmi.nor.df <- cbind(hdmi.rfd[,c(1:2,10003)], hdmi.nor)
hdmi.nor.end <- left_join(hdmi.nor.df, who.hdmi.df[,1:2], by="CAS")
hdmi.nor.end$end <- paste("CS =", hdmi.nor.end$endpoint.x, ", WHO =", hdmi.nor.end$endpoint.y)
hdmi.nor.end.df <- hdmi.nor.end[-c(2,5,7,17,18,24,25),]
hdmi.nor.end.df$`5th_quantile` <- apply(hdmi.nor.end.df[,4:10003], 1, quantile, probs=0.05)
hdmi.nor.sen.end.df <- hdmi.nor.end.df %>% group_by(CAS) %>% slice(which.min(`5th_quantile`))
hdmi.nor.sen.m <- melt(hdmi.nor.sen.end.df, id=c(1:3,10004:10006))
hdmi.nor.sen.m$name <- "chemical-specific"

who.hdmi.rfd <- left_join(who.hdmi.df, rfd[,c(1,2,4)], by="CAS")
who.hdmi.nor <- with(who.hdmi.rfd, who.hdmi.rfd[,3:10002]/who.hdmi.rfd[,10004])
who.hdmi.nor.df <- cbind(who.hdmi.rfd[,c(1:2,10003)], who.hdmi.nor)
who.hdmi.nor.end <- left_join(hdmi.df[,1:2], who.hdmi.nor.df, by="CAS")
who.hdmi.nor.end$end <- paste("CS =", hdmi.nor.end$endpoint.x, ", WHO =", hdmi.nor.end$endpoint.y)
who.hdmi.nor.end.df <- who.hdmi.nor.end[-c(2,5,7,17,18,24,25),]
who.hdmi.nor.end.df$`5th_quantile` <- apply(who.hdmi.nor.end.df[,5:10004], 1, quantile, probs=0.05)
who.hdmi.nor.sen.end.df <- who.hdmi.nor.end.df %>% group_by(CAS) %>% slice(which.min(`5th_quantile`))
who.hdmi.nor.sen.m <- melt(who.hdmi.nor.sen.end.df, id=c(1:4,10005,10006))
who.hdmi.nor.sen.m$name <- "WHO"

hdmi.sen.bind <- rbind(hdmi.nor.sen.m, who.hdmi.nor.sen.m)
hdmi.sen.bind$order <- bmd.sen.bind$`5th_quantile`
hdmi.sen.bind[1:190000,10] <- hdmi.sen.bind[190001:380000,10]
hdmi.sen.bind$order2 <- hdmi.sen.bind$`5th_quantile`
hdmi.sen.bind[1:190000,11] <- hdmi.sen.bind[190001:380000,11]
hdmi.sen.bind$name <- factor(hdmi.sen.bind$name, levels = c("WHO","chemical-specific"))

hdmi.point <- rbind(hdmi.nor.sen.end.df[,c(1,3,10006)],who.hdmi.nor.sen.end.df[,c(1,4,10006)])
hdmi.point$name <- "chemical-specific"
hdmi.point$name[20:38] <- "WHO"
hdmi.point$name <- factor(hdmi.point$name, levels = c("WHO","chemical-specific"))

hdmi.sen.plot <- ggplot()+
  geom_density_ridges(data = hdmi.sen.bind, aes(x=value, y=reorder(Chemical,order), color=name, fill=name), 
                      scale=0.8, alpha=0.4)+
  geom_point(data=hdmi.point, aes(x=`5th_quantile`, y=Chemical, color=name))+
  geom_vline(xintercept = 1, color="black", linetype="dashed")+
  scale_x_log10(limits=c(1e-03, 1e+03),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme_bw()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), 
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 10), 
        legend.text = element_text(size = 10), legend.title = element_blank(),
        legend.position = "top", legend.key.size = unit(3, 'mm'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.margin = unit(c(5.5, 5.5, 5.5, 0),"pt"))+
  xlab(bquote(HD[M]^I/Regulatory~RfD))+
  scale_color_manual(breaks=c("WHO","chemical-specific"), 
                     labels=c("WHO/IPCS approximate approach","Chemical-specific approach"),
                     values=c(`chemical-specific`="#D55E00",`WHO`="#56B4E9"))+
  scale_fill_manual(breaks=c("WHO","chemical-specific"), 
                    labels=c("WHO/IPCS approximate approach","Chemical-specific approach"),
                    values=c(`chemical-specific`="#D55E00",`WHO`="#56B4E9"))+
  annotation_logticks(sides="b")+
  guides(color = guide_legend(nrow = 2), fill = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))

sen.plot <- ggarrange(bmd.sen.plot, hdmi.sen.plot, ncol=2, labels=c("A","B"), widths = c(1.1/2,0.9/2), align="h")
ggsave(sen.plot, file="the most sensitive endpoint BBMD HDMI plot.pdf", width=15, height=14, path="Plot output", scale=0.7)

########------------------------------------------------------------
###all endpoints
###BMD
#normalized by regulatory POD
bmd.pod <- left_join(bmd, rfd[,c(1:3)], by="CAS")
bmd.nor <- with(bmd.pod, bmd.pod[,1:10000]/bmd.pod[,10004])
bmd.nor.df <- cbind(bmd.pod[,10001:10003], bmd.nor)
bmd.nor.end <- left_join(bmd.nor.df, who.bmd.df[,10001:10002], by="CAS")
bmd.nor.end$end <- paste("CS =", bmd.nor.end$endpoint.x, ", WHO =", bmd.nor.end$endpoint.y)
bmd.nor.end.df <- bmd.nor.end[-c(3,4,21,22,24,27,29),]
bmd.nor.end.df$`5th_quantile` <- apply(bmd.nor.end.df[,4:10003], 1, quantile, probs=0.05)
bmd.nor.end.df <- bmd.nor.end.df %>% arrange(CAS)
bmd.nor.m <- melt(bmd.nor.end.df, id=c(1:3,10004:10006))
bmd.nor.m$name <- "Bayesian BMD or WHO/IPCS (2018) BMD"

who.bmd.pod <- left_join(who.bmd.df, rfd[,c(1:3)], by="CAS")
who.bmd.nor <- with(who.bmd.pod, who.bmd.pod[,1:10000]/who.bmd.pod[,10004])
who.bmd.nor.df <- cbind(who.bmd.pod[,10001:10003], who.bmd.nor)
who.bmd.nor.end <- left_join(bmd[,10001:10002], who.bmd.nor.df, by="CAS")
who.bmd.nor.end$end <- paste("CS =", who.bmd.nor.end$endpoint.x, ", WHO =", who.bmd.nor.end$endpoint.y)
who.bmd.nor.end.df <- who.bmd.nor.end[-c(3,4,21,22,24,27,29),]
who.bmd.nor.end.df$`5th_quantile` <- apply(who.bmd.nor.end.df[,5:10004], 1, quantile, probs=0.05)
who.bmd.nor.end.df <- who.bmd.nor.end.df %>% arrange(CAS)
who.bmd.nor.m <- melt(who.bmd.nor.end.df, id=c(1:4,10005:10006))
who.bmd.nor.m$name <- "WHO/IPCS (2018) BMD"

bmd.bind <- rbind(bmd.nor.m, who.bmd.nor.m)
bmd.bind$order <- bmd.bind$`5th_quantile`
bmd.bind[1:350000,10] <- bmd.bind[350001:700000,10]
bmd.bind$name <- factor(bmd.bind$name, levels = c("WHO/IPCS (2018) BMD","Bayesian BMD or WHO/IPCS (2018) BMD"))

bmd.all.point <- rbind(bmd.nor.end.df[,c(1,2,10006)],who.bmd.nor.end.df[,c(1,2,10006)])
bmd.all.point$name <- "Bayesian BMD or WHO/IPCS (2018) BMD"
bmd.all.point$name[36:70] <- "WHO/IPCS (2018) BMD"
bmd.all.point$name <- factor(bmd.all.point$name, levels = c("WHO/IPCS (2018) BMD","Bayesian BMD or WHO/IPCS (2018) BMD"))

bmd.plot <- ggplot()+
  geom_density_ridges(data = bmd.bind, 
                      aes(x=value, y=reorder(endpoint.x, order), color=name, fill=name), 
                      scale=0.8, alpha=0.6)+
  geom_point(data=bmd.all.point, aes(x=`5th_quantile`, y=endpoint.x, color=name))+
  geom_vline(xintercept = 1, color="black", linetype="dashed")+
  scale_x_log10(limits=c(1e-02, 1e+03),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme_bw()+
  theme(axis.title.y = element_blank(), legend.title = element_blank(),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12), legend.text = element_text(size = 12),
        legend.position = "top",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab(bquote(BMD[M]^I/Regulatory~POD))+
  scale_color_manual(breaks=c("WHO/IPCS (2018) BMD","Bayesian BMD or WHO/IPCS (2018) BMD"),
                     values=c(`WHO/IPCS (2018) BMD`="#56B4E9", `Bayesian BMD or WHO/IPCS (2018) BMD`="#D55E00"),
                     labels=c("WHO/IPCS (2018) BMD"="WHO/IPCS BMD","Bayesian BMD or WHO/IPCS (2018) BMD"="BBMD or WHO/IPCS BMD"))+
  scale_fill_manual(breaks=c("WHO/IPCS (2018) BMD","Bayesian BMD or WHO/IPCS (2018) BMD"),
                    values=c(`WHO/IPCS (2018) BMD`="#56B4E9", `Bayesian BMD or WHO/IPCS (2018) BMD`="#D55E00"),
                    labels=c("WHO/IPCS (2018) BMD"="WHO/IPCS BMD","Bayesian BMD or WHO/IPCS (2018) BMD"="BBMD or WHO/IPCS BMD"))+
  annotation_logticks(sides="b")+
  guides(color = guide_legend(nrow = 2), fill = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))

###HDMI
#normalized by regulatory RfD
hdmi.rfd <- left_join(hdmi.df, rfd[,c(1,2,4)], by="CAS")
hdmi.nor <- with(hdmi.rfd, hdmi.rfd[,3:10002]/hdmi.rfd[,10004])
hdmi.nor.df <- cbind(hdmi.rfd[,c(1:2,10003)], hdmi.nor)
hdmi.nor.end <- left_join(hdmi.nor.df, who.hdmi.df[,1:2], by="CAS")
hdmi.nor.end$end <- paste("CS =", hdmi.nor.end$endpoint.x, ", WHO =", hdmi.nor.end$endpoint.y)
hdmi.nor.end.df <- hdmi.nor.end[-c(2,5,7,17,18,24,25),]
hdmi.nor.end.df$`5th_quantile` <- apply(hdmi.nor.end.df[,4:10003], 1, quantile, probs=0.05)
hdmi.nor.end.df <- hdmi.nor.end.df %>% arrange(CAS)
hdmi.nor.m <- melt(hdmi.nor.end.df, id=c(1:3,10004:10006))
hdmi.nor.m$name <- "chemical-specific"

who.hdmi.rfd <- left_join(who.hdmi.df, rfd[,c(1,2,4)], by="CAS")
who.hdmi.nor <- with(who.hdmi.rfd, who.hdmi.rfd[,3:10002]/who.hdmi.rfd[,10004])
who.hdmi.nor.df <- cbind(who.hdmi.rfd[,c(1:2,10003)], who.hdmi.nor)
who.hdmi.nor.end <- left_join(hdmi.df[,1:2], who.hdmi.nor.df, by="CAS")
who.hdmi.nor.end$end <- paste("CS =", hdmi.nor.end$endpoint.x, ", WHO =", hdmi.nor.end$endpoint.y)
who.hdmi.nor.end.df <- who.hdmi.nor.end[-c(2,5,7,17,18,24,25),]
who.hdmi.nor.end.df$`5th_quantile` <- apply(who.hdmi.nor.end.df[,5:10004], 1, quantile, probs=0.05)
who.hdmi.nor.end.df <- who.hdmi.nor.end.df %>% arrange(CAS)
who.hdmi.nor.m <- melt(who.hdmi.nor.end.df, id=c(1:4,10005,10006))
who.hdmi.nor.m$name <- "WHO"

hdmi.bind <- rbind(hdmi.nor.m, who.hdmi.nor.m)
hdmi.bind$order <- bmd.bind$`5th_quantile`
hdmi.bind[1:350000,10] <- hdmi.bind[350001:700000,10]
hdmi.bind$name <- factor(hdmi.bind$name, levels = c("WHO","chemical-specific"))

hdmi.allpoint <- rbind(hdmi.nor.end.df[,c(1,2,10006)],who.hdmi.nor.end.df[,c(1,2,10006)])
hdmi.allpoint$name <- "chemical-specific"
hdmi.allpoint$name[36:70] <- "WHO"
hdmi.allpoint$name <- factor(hdmi.allpoint$name, levels = c("WHO","chemical-specific"))

hdmi.plot <- ggplot()+
  geom_density_ridges(data = hdmi.bind,
                      aes(x=value, y=reorder(endpoint.x, order), color=name, fill=name), 
                      scale=0.8, alpha=0.4)+
  geom_point(data=hdmi.allpoint, aes(x=`5th_quantile`, y=endpoint.x, color=name))+
  geom_vline(xintercept = 1, color="black", linetype="dashed")+
  scale_x_log10(limits=c(1e-03, 1e+03),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme_bw()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), 
        axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 12), 
        legend.text = element_text(size = 12), legend.title = element_blank(),
        legend.position = "top", 
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.margin = unit(c(5.5, 5.5, 5.5, 0),"pt"))+
  xlab(bquote(HD[M]^I/Regulatory~RfD))+
  scale_color_manual(breaks=c("WHO","chemical-specific"), 
                     labels=c("WHO/IPCS approximate approach","Chemical-specific approach"),
                     values=c(`chemical-specific`="#D55E00",`WHO`="#56B4E9"))+
  scale_fill_manual(breaks=c("WHO","chemical-specific"), 
                    labels=c("WHO/IPCS approximate approach","Chemical-specific approach"),
                    values=c(`chemical-specific`="#D55E00",`WHO`="#56B4E9"))+
  annotation_logticks(sides="b")+
  guides(color = guide_legend(nrow = 2), fill = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))

allend.plot <- ggarrange(bmd.plot, hdmi.plot, ncol=2, labels=c("A","B"), widths = c(1.25/2,0.75/2))
ggsave(allend.plot, file="all endpoints normalized BBMD HDMI ggridges plot.pdf", width=15, height=14, path="Plot output")
