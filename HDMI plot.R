library(reshape2)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggridges)
library(ggpubr)

yaxis_order <- c("Pentachlorophenol","Chlorpyrifos","Dieldrin","Endosulfan","Azinphos-methyl","Disulfoton","2,4,5-Trichlorophenol",
                 "Heptachlor epoxide","Ethion","Dicofol","p,p'-DDT","Parathion","4,6-Dinitro-o-cresol","Methoxychlor","Dibutyl phthalate","Diazinon",
                 "Endrin","Heptachlor","Aldrin")

rfd <- read.csv("regulatory reference dose.csv")[,c(1,2,5,7,8)]

hdmi.df <- read.csv("HDMI.csv")
who.hdmi.df <- read.csv("WHO approximate HDMI.csv")

###the most sensitive endpoint
###HDMI
#normalized by regulatory rfd
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
hdmi.sen.bind$Chemical <- factor(hdmi.sen.bind$Chemical, levels = yaxis_order)
hdmi.sen.bind$name <- factor(hdmi.sen.bind$name, levels = c("chemical-specific","WHO"))

hdmi.point <- rbind(hdmi.nor.sen.end.df[,c(1,3,10006)],who.hdmi.nor.sen.end.df[,c(1,4,10006)])
hdmi.point$name <- "chemical-specific"
hdmi.point$name[20:38] <- "WHO"
hdmi.point$Chemical <- factor(hdmi.point$Chemical, levels = yaxis_order)
hdmi.point$name <- factor(hdmi.point$name, levels = c("chemical-specific","WHO"))

#set whisker to [5%,95%]
quantiles_95 <- function(x) {
  r <- quantile(x, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

hdmi.sen.plot <- ggplot()+
  stat_summary(fun.data = quantiles_95, data = hdmi.sen.bind, geom="boxplot", aes(x=value, y=Chemical, color=name), outlier.shape=NA, position = position_dodge(width = 0.75))+
  geom_point(data=hdmi.point, aes(x=`5th_quantile`, y=Chemical, color=name), position = position_dodge(width = 0.75))+
  geom_vline(xintercept = 1, color="black", linetype="dashed")+
  scale_x_log10(limits=c(1e-03, 1e+03),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme_bw()+
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 10), 
        legend.title = element_blank(),
        legend.position = "top", legend.key.size = unit(3, 'mm'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab(bquote(HD[M]^I/Regulatory~RfD))+
  scale_color_manual(breaks=c("WHO","chemical-specific"),
                     values=c(`WHO`="#00BFC4", `chemical-specific`="#F8766D"),
                     labels=c("WHO"="WHO/IPCS approximate approach","chemical-specific"="Chemical-specific approach"))+
  guides(color = guide_legend(nrow = 2))+
  annotation_logticks(sides="b")
print(hdmi.sen.plot)

hdmi.sen.total <- ggplot(data=hdmi.point, aes(x=`5th_quantile`, y=name, color=name))+
  geom_boxplot()+
  geom_point(position = position_jitter(width = 0.2))+
  geom_vline(xintercept = 1, color="black", linetype="dashed")+
  scale_x_log10(limits=c(1e-03, 1e+03),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_discrete(limits=c("chemical-specific","WHO"),
                   labels=c("Chemical-specific approach","WHO/IPCS approximate approach"))+
  theme_bw()+
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 10), 
        legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab(bquote("5th"~quantile~of~HD[M]^I/Regulatory~RfD))+
  annotation_logticks(sides="b")
print(hdmi.sen.total)

#HDMI differences
#use the most sensitive endpoints
chem <- read.csv("regulatory reference dose.csv")[,1:2]
bbmd <- read.csv("BBMD quantile.csv")[,4:5]
bbmd.chem <- data.frame(CAS=unique(bbmd$CAS))
bbmd.chem$BBMD <- "BBMD"

who.hdmi.quan <- read.csv("WHO approximate HDMI quantile.csv")
cs.hdmi.quan <- read.csv("HDMI quantile.csv")
hdmi.join <- left_join(cs.hdmi.quan, who.hdmi.quan, by="CAS")
hdmi <- merge(hdmi.join, chem, by="CAS")
hdmi.35 <- hdmi[-c(2,5,7,17,18,24,25),]
hdmi.min <- hdmi.35 %>% group_by(CAS) %>% slice(which.min(`X5..x`))
hdmi.df <- left_join(hdmi.min, bbmd.chem, by="CAS")
hdmi.df <- hdmi.df %>% arrange(CAS)
hdmi.df$log10.cshdmi50 <- log10(hdmi.df$X50..x)
hdmi.df$log10.whohdmi50 <- log10(hdmi.df$X50..y)
hdmi.df$log10.cs_who <- hdmi.df$log10.cshdmi50-hdmi.df$log10.whohdmi50
hdmi.df$Chemical <- factor(hdmi.df$Chemical, levels = yaxis_order)

hdmi.diff <- ggplot(hdmi.df, aes(x=log10.cs_who, y=Chemical))+
  geom_bar(aes(fill = BBMD == BBMD), stat="identity")+
  theme_bw()+
  theme(#axis.title.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.y = element_blank(), 
        legend.text = element_text(size = 10), legend.title = element_blank(),
        legend.position = "top", legend.key.size = unit(3, 'mm'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(5.5, 5.5, 0, 5.5),"pt"))+
  geom_vline(xintercept = 0)+
  xlim(-2,2)+
  scale_fill_manual(breaks = c(F, T), values=c("gray", "gray3"), labels=c("Use WHO/IPCS BMD","Use BBMD"))+
  xlab(bquote(Delta~log[10]~HD[M]^I))
print(hdmi.diff)

hdmi.m <- melt(hdmi.df[,c(1,10,14)])
hdmi.diff.box <- ggplot(hdmi.m, aes(x=value, y=variable))+
  geom_boxplot()+
  geom_jitter(shape=20, position=position_jitter(0.2))+
  theme_bw()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.margin = unit(c(0, 5.5, 5.5, 5.5),"pt"))+
  geom_vline(xintercept = 0, linetype="dashed")+
  xlim(-2,2)+
  xlab(bquote(Delta~log[10]~HD[M]^I))
print(hdmi.diff.box)

#check which one contributes more to the shift in the overall distribution
#use absolute value to compare the difference
abs.df <- merge(bmd.df[,c(1,10,14)], ufh.df[,c(6,13)], by="Chemical")
abs.df$abs.bmd <- abs(abs.df$log10.bmd_whobmd)
abs.df$abs.ufh <- abs(abs.df$log10.cs_who)
abs.df$diff <- abs.df$abs.bmd-abs.df$abs.ufh
#BMD is smaller:  15 chemicals (total 19 chemicals)

#degree of uncertainty
chem <- read.csv("regulatory reference dose.csv")[,1:2]
who.hdmi.quan <- read.csv("WHO approximate HDMI quantile.csv")
cs.hdmi.quan <- read.csv("HDMI quantile.csv")

cs.hdmi.quan.chem <- merge(cs.hdmi.quan, chem, by="CAS")
who.hdmi.quan.chem <- merge(who.hdmi.quan, chem, by="CAS")

cs.hdmi.quan.min <- cs.hdmi.quan.chem %>% group_by(CAS) %>% slice(which.min(`X5.`))
who.hdmi.quan.min <- who.hdmi.quan.chem %>% group_by(CAS) %>% slice(which.min(`X5.`))

who.hdmi.quan.min$degree <- with(who.hdmi.quan.min, X95./X5.)
who.hdmi.quan.min$name <- "WHO/IPCS approximate approach"
cs.hdmi.quan.min$degree <- with(cs.hdmi.quan.min, X95./X5.)
cs.hdmi.quan.min$name <- "Chemical-specific approach"

hdmi.degree <- rbind(who.hdmi.quan.min, cs.hdmi.quan.min)
hdmi.degree$Chemical <- factor(hdmi.degree$Chemical, levels = yaxis_order)
hdmi.degree$name <- factor(hdmi.degree$name, levels = c("WHO/IPCS approximate approach","Chemical-specific approach"))

hdmi.degree.dumb <- ggplot()+
  geom_point(data=hdmi.degree, aes(x=degree, y=Chemical, color=name), size=3)+
  geom_segment(aes(x =who.hdmi.quan.min$degree , y = who.hdmi.quan.min$Chemical, xend = cs.hdmi.quan.min$degree, yend = cs.hdmi.quan.min$Chemical),
               arrow = arrow(length = unit(0.2, "cm")))+
  scale_x_log10(limits=c(1, 1e+03),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme_bw()+
  theme(axis.title.y = element_blank(), 
        legend.title = element_blank(), legend.position = "top", legend.text = element_text(size = 10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("95th / 5th ratio")+
  scale_color_manual(breaks=c("WHO/IPCS approximate approach","Chemical-specific approach"),
                     values=c(`WHO/IPCS approximate approach`="#00BFC4", `Chemical-specific approach`="#F8766D"))+
  guides(color = guide_legend(nrow = 2))+
  annotation_logticks(sides="b")
print(hdmi.degree.dumb)

hdmi.degree.box <- ggplot(data=hdmi.degree, aes(x=degree, y=name, color=name))+
  geom_boxplot()+
  geom_jitter(shape=20, position=position_jitter(0.2))+
  scale_x_log10(limits=c(1, 1e+03),
                breaks = scales::trans_breaks("log10", function(x) 10^x, n=4),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_discrete(limits=c("Chemical-specific approach","WHO/IPCS approximate approach"))+
  theme_bw()+
  theme(axis.title.y = element_blank(), 
        axis.ticks.y=element_blank(),
        legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("95th / 5th ratio")+
  scale_color_manual(breaks=c("WHO/IPCS approximate approach","Chemical-specific approach"),
                     values=c(`WHO/IPCS approximate approach`="#00BFC4", `Chemical-specific approach`="#F8766D"))+
  guides(color = guide_legend(nrow = 2))+
  annotation_logticks(sides="b")
print(hdmi.degree.box)

##check log variance of HDMI
#log variance is proportional to (log10(95th/5th))^2
who.hdmi.log10ratio <- data.frame(chemical=who.hdmi.quan.min$Chemical, CAS=who.hdmi.quan.min$CAS,
                                  uncertainty=who.hdmi.quan.min$degree, log10ratio=log10(who.hdmi.quan.min$degree))
cs.hdmi.log10ratio <- data.frame(chemical=cs.hdmi.quan.min$Chemical, CAS=cs.hdmi.quan.min$CAS,
                                 uncertainty=cs.hdmi.quan.min$degree, log10ratio=log10(cs.hdmi.quan.min$degree))
#fold reduction in log variance
who.median <- median(who.hdmi.log10ratio$log10ratio)
cs.median <- median(cs.hdmi.log10ratio$log10ratio)
hdmifold <- who.median^2/cs.median^2
#1.38 fold reduction in log variance

#combine plots
#HDMI/Reg RfD
norm.hdmi.plot <- ggarrange(hdmi.sen.plot, hdmi.sen.total, nrow=2, heights = c(0.75,0.25), align="v")
#HDMI difference
hdmi.diff.plot <- ggarrange(hdmi.diff, hdmi.diff.box, nrow=2, heights = c(0.75,0.25), align="v")
#HDMI degree of uncertainty
hdmi.degree.plot <- ggarrange(hdmi.degree.dumb, hdmi.degree.box, nrow=2, heights = c(0.75,0.25), align="v")

hdmi.bind <- ggarrange(norm.hdmi.plot, hdmi.diff.plot, hdmi.degree.plot, ncol=3, align="hv")
ggsave(hdmi.bind, file="HDMI plot.pdf", width = 24, height = 12, path = "HDMI plots", scale=0.7)

###-------------------------------------------------------------------------------------------------
#generate csv file for Prism
#A panel normalized distributions
hdmi.nor.sen.end.df$Chemical <- factor(hdmi.nor.sen.end.df$Chemical, levels = yaxis_order)
hdmi.nor.sen.end.df_sort <- hdmi.nor.sen.end.df[order(hdmi.nor.sen.end.df$Chemical),]
hdmi.nor.sen.end.df.t <- data.frame(t(hdmi.nor.sen.end.df_sort[,3:10003]))
colnames(hdmi.nor.sen.end.df.t) <- hdmi.nor.sen.end.df.t[1,]
hdmi.data <- hdmi.nor.sen.end.df.t[-1,]
write.csv(hdmi.data, file="prism file/Fig 4 (A) chemical-specific HDMI normalized data.csv", row.names = FALSE)

who.hdmi.nor.sen.end.df$Chemical <- factor(who.hdmi.nor.sen.end.df$Chemical, levels = yaxis_order)
who.hdmi.nor.sen.end.df_sort <- who.hdmi.nor.sen.end.df[order(who.hdmi.nor.sen.end.df$Chemical),]
who.hdmi.nor.sen.end.df.t <- data.frame(t(who.hdmi.nor.sen.end.df_sort[,4:10004]))
colnames(who.hdmi.nor.sen.end.df.t) <- who.hdmi.nor.sen.end.df.t[1,]
who.hdmi.data <- who.hdmi.nor.sen.end.df.t[-1,]
write.csv(who.hdmi.data, file="prism file/Fig 4 (A) WHO HDMI normalized data.csv", row.names = FALSE)

hdmi.nor.5th <- hdmi.nor.sen.end.df[,c(3,10006)]
hdmi.nor.5th$name <- "Chemical-specific intra"
who.hdmi.nor.5th <- who.hdmi.nor.sen.end.df[,c(4, 10006)]
who.hdmi.nor.5th$name <- "WHO/IPCS intra"
hdmi.5th.bind <- merge(hdmi.nor.5th, who.hdmi.nor.5th, by="Chemical")
hdmi.5th.bind$Chemical <- factor(hdmi.5th.bind$Chemical, levels = yaxis_order)
hdmi.5th.bind_sort <- hdmi.5th.bind[order(hdmi.5th.bind$Chemical),]
write.csv(hdmi.5th.bind_sort, file="prism file/Fig 4 (D) HDMI normalized 5th.csv", row.names = FALSE)

#B panel for differences
hdmi.diff.df <- hdmi.df[,c(10,11,14)]
hdmi.diff.df$Chemical <- factor(hdmi.diff.df$Chemical, levels = yaxis_order)
hdmi.diff.df_sort <- hdmi.diff.df[order(hdmi.diff.df$Chemical),]
write.csv(hdmi.diff.df_sort, file="prism file/Fig 4 (B)(E) HDMI differences.csv", row.names = FALSE)

#C panel for uncertainty degree
hdmi.degree.df <- merge(hdmi.degree[1:19,6:8], hdmi.degree[20:38,6:8], by="Chemical")
hdmi.degree.df$Chemical <- factor(hdmi.degree.df$Chemical, levels = yaxis_order)
hdmi.degree.df_sort <- hdmi.degree.df[order(hdmi.degree.df$Chemical),]
write.csv(hdmi.degree.df_sort, file="prism file/Fig 4 (C)(F) HDMI uncertainty degree.csv", row.names = FALSE)
