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

bmd <- read.csv("bmd samples CAS by row.csv")
who.bmd.df <- read.csv("WHO bmd samples CAS by row.csv")

###normalized BMD between two approaches
###the most sensitive endpoint
#normalized by regulatory POD
bmd.pod <- left_join(bmd, rfd[,c(1:3)], by="CAS")
bmd.nor <- with(bmd.pod, bmd.pod[,1:10000]/bmd.pod[,10004])
bmd.nor.df <- cbind(bmd.pod[,10001:10003], bmd.nor)
bmd.nor.end <- left_join(bmd.nor.df, who.bmd.df[,10001:10002], by="CAS")
bmd.nor.end$end <- paste("CS =", bmd.nor.end$endpoint.x, ", WHO =", bmd.nor.end$endpoint.y)
bmd.nor.end.df <- bmd.nor.end[-c(3,4,21,22,24,27,29),]
bmd.nor.end.df$`5th_quantile` <- apply(bmd.nor.end.df[,4:10003], 1, quantile, probs=0.05)
bmd.nor.sen.end.df <- bmd.nor.end.df %>% group_by(CAS) %>% slice(which.min(`5th_quantile`))
bbmd.nor.sen.end.df <- bmd.nor.sen.end.df[bmd.nor.sen.end.df$Chemical %in% c("Aldrin", "Heptachlor","Methoxychlor",
                                                                             "Dicofol","Azinphos-methyl","Endosulfan",
                                                                             "Dieldrin","Chlorpyrifos","Pentachlorophenol"),]
bbmd.nor.sen.m <- melt(bbmd.nor.sen.end.df, id=c(1:3,10004:10006))
bbmd.nor.sen.m$name <- "BBMD"

#bmd.nor.sen.m <- melt(bmd.nor.sen.end.df, id=c(1:3,10004:10006))
#bmd.nor.sen.m$name <- "BBMD or WHO/IPCS BMD"

who.bmd.pod <- left_join(who.bmd.df, rfd[,c(1:3)], by="CAS")
who.bmd.nor <- with(who.bmd.pod, who.bmd.pod[,1:10000]/who.bmd.pod[,10004])
who.bmd.nor.df <- cbind(who.bmd.pod[,10001:10003], who.bmd.nor)
who.bmd.nor.end <- left_join(bmd[,10001:10002], who.bmd.nor.df, by="CAS")
who.bmd.nor.end$end <- paste("CS =", who.bmd.nor.end$endpoint.x, ", WHO =", who.bmd.nor.end$endpoint.y)
who.bmd.nor.end.df <- who.bmd.nor.end[-c(3,4,21,22,24,27,29),]
who.bmd.nor.end.df$`5th_quantile` <- apply(who.bmd.nor.end.df[,5:10004], 1, quantile, probs=0.05)
who.bmd.nor.sen.end.df <- who.bmd.nor.end.df %>% group_by(CAS) %>% slice(which.min(`5th_quantile`))
who.bmd.nor.sen.m <- melt(who.bmd.nor.sen.end.df, id=c(1:4,10005:10006))
who.bmd.nor.sen.m$name <- "WHO/IPCS BMD"

bmd.sen.bind <- rbind(bbmd.nor.sen.m, who.bmd.nor.sen.m)
bmd.sen.bind$Chemical <- factor(bmd.sen.bind$Chemical, levels = yaxis_order)
bmd.sen.bind$name <- factor(bmd.sen.bind$name, levels = c("BBMD","WHO/IPCS BMD"))

bmd.point <- rbind(bbmd.nor.sen.end.df[,c(1,3,10006)],who.bmd.nor.sen.end.df[,c(1,4,10006)])
bmd.point$name <- "BBMD"
bmd.point$name[10:28] <- "WHO/IPCS BMD"
bmd.point$Chemical <- factor(bmd.point$Chemical, levels = yaxis_order)
bmd.point$name <- factor(bmd.point$name, levels = c("BBMD","WHO/IPCS BMD"))

#set whisker to [5%,95%]
quantiles_95 <- function(x) {
  r <- quantile(x, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

bmd.sen.plot <- ggplot()+
  stat_summary(fun.data = quantiles_95, data = bmd.sen.bind, geom="boxplot", aes(x=value, y=Chemical, fill=name), outlier.shape=NA, position = position_dodge(width = 0.75))+
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
  xlab(bquote(BMD[M]/Regulatory~POD))+
  scale_fill_manual(breaks=c("WHO/IPCS BMD","BBMD"),
                     values=c(`WHO/IPCS BMD`="white", `BBMD`="red"))+
  annotation_logticks(sides="b")
print(bmd.sen.plot)

bmd.sen.total <- ggplot(data=bmd.point, aes(x=`5th_quantile`, y=name, fill=name))+
  geom_boxplot()+
  geom_point(position = position_jitterdodge(jitter.width = 0, dodge.width = 0), size=5, shape=1)+
  geom_vline(xintercept = 1, color="black", linetype="dashed")+
  scale_x_log10(limits=c(1e-02, 1e+03),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_discrete(limits=c("BBMD","WHO/IPCS BMD"))+
  theme_bw()+
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 10), 
        legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab(bquote("5"^th~'%ile'~of~BMD[M]/Regulatory~POD))+
  scale_fill_manual(breaks=c("WHO/IPCS BMD","BBMD"),
                     values=c(`WHO/IPCS BMD`="white", `BBMD`="red"))+
  annotation_logticks(sides="b")
print(bmd.sen.total)

#BMD differences
chem <- read.csv("regulatory reference dose.csv")[,1:2]
bbmd <- read.csv("BBMD quantile.csv")[,4:5]
bbmd.chem <- data.frame(CAS=unique(bbmd$CAS))
bbmd.chem$BBMD <- "BBMD"

bmd.quan <- read.csv("BMD quantile.csv")
who.bmd.quan <- read.csv("WHO BMD quantile.csv")
bmd.join <- left_join(bmd.quan, who.bmd.quan, by="CAS")
bmd <- merge(bmd.join, chem, by="CAS")
bmd.35 <- bmd[-c(2,5,7,17,18,24,25),]
bmd.min <- bmd.35 %>% group_by(CAS) %>% slice(which.min(`X5..x`))
bmd.df <- left_join(bmd.min, bbmd.chem, by="CAS")
bmd.df <- bmd.df %>% arrange(CAS)
bmd.df$log10.bmd50 <- log10(bmd.df$X50..x)
bmd.df$log10.who.bmd50 <- log10(bmd.df$X50..y)
bmd.df$log10.bmd_whobmd <- bmd.df$log10.bmd50-bmd.df$log10.who.bmd50
bmd.df$Chemical <- factor(bmd.df$Chemical, levels = yaxis_order)
bmd.df$log10.bmd_whobmd[is.na(bmd.df$BBMD)] <- 0 #no difference if the chemicals has no available dose-response data for BBMD

bmd.diff <- ggplot(bmd.df, aes(x=log10.bmd_whobmd, y=Chemical))+
  geom_bar(aes(fill = BBMD == BBMD), stat="identity")+
  theme_bw()+
  theme(#axis.title.x = element_blank(), #axis.text.x = element_blank(), 
        axis.title.y = element_blank(), #axis.text.y = element_blank(), 
        legend.text = element_text(size = 10), legend.title = element_blank(),
        legend.position = "top", legend.key.size = unit(3, 'mm'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  geom_vline(xintercept = 0)+
  xlim(-2,2)+
  scale_fill_manual(breaks = c(F, T), values=c(NA, "gray"), labels=c("Use WHO/IPCS BMD","Use BBMD"))+
  xlab(bquote(Delta~Log[10]~BMD[M]))
print(bmd.diff)

bmd.df.m <- melt(bmd.df[,c(1,10,14)])
bmd.diff.box <- ggplot(bmd.df.m, aes(x=value, y=variable))+
  geom_boxplot(fill="gray")+
  geom_point(size=5, shape=1)+
  theme_bw()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  geom_vline(xintercept = 0, linetype="dashed")+
  xlim(-2,2)+
  xlab(bquote(Delta~Log[10]~BMD[M]))
print(bmd.diff.box)

#degree of uncertainty
chem <- read.csv("regulatory reference dose.csv")[,1:2]
bbmd.quan <- read.csv("BBMD quantile.csv")
who.bmd.quan <- read.csv("WHO BMD quantile.csv")

bbmd.quan.chem <- merge(bbmd.quan, chem, by="CAS")
who.bmd.quan.chem <- merge(who.bmd.quan, chem, by="CAS")

bbmd.quan.min <- bbmd.quan.chem %>% group_by(CAS) %>% slice(which.min(`X5.`))
who.bmd.quan.min <- who.bmd.quan.chem %>% group_by(CAS) %>% slice(which.min(`X5.`))

bbmd.quan.min$degree <- with(bbmd.quan.min, X95./X5.)
bbmd.quan.min$name <- "BBMD"
who.bmd.quan.min$degree <- with(who.bmd.quan.min, X95./X5.)
who.bmd.quan.min$name <- "WHO/IPCS BMD"

who.bmd.quan.min.arrow <- who.bmd.quan.min[who.bmd.quan.min$Chemical %in% c("Aldrin", "Heptachlor","Methoxychlor",
                                                                            "Dicofol","Azinphos-methyl","Endosulfan",
                                                                            "Dieldrin","Chlorpyrifos","Pentachlorophenol"),]

bmd.degree <- rbind(who.bmd.quan.min, bbmd.quan.min)
bmd.degree$Chemical <- factor(bmd.degree$Chemical, levels = yaxis_order)
bmd.degree$name <- factor(bmd.degree$name, levels = c("WHO/IPCS BMD","BBMD"))

bmd.degree.dumb <- ggplot()+
  geom_point(data=bmd.degree, aes(x=degree, y=Chemical, color=name, shape=name), size=3)+
  geom_segment(aes(x =who.bmd.quan.min.arrow$degree , y = who.bmd.quan.min.arrow$Chemical, xend = bbmd.quan.min$degree, yend = bbmd.quan.min$Chemical),
               arrow = arrow(length = unit(0.2, "cm")))+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme_bw()+
  theme(axis.title.y = element_blank(), 
        legend.title = element_blank(), legend.position = "top", legend.text = element_text(size = 10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab(bquote("95"^"th"/"5"^"th"~"%ile Ratio"))+
  scale_color_manual(name="label",
                     breaks=c("WHO/IPCS BMD","BBMD"),
                     values=c(`WHO/IPCS BMD`="black", `BBMD`="red"))+
  scale_shape_manual(name="label",
                     breaks=c("WHO/IPCS BMD","BBMD"),
                     values=c(`WHO/IPCS BMD`=1, `BBMD`=16))+
  annotation_logticks(sides="b")
print(bmd.degree.dumb)

bmd.degree.box <- ggplot(data=bmd.degree, aes(x=degree, y=name, fill=name))+
  geom_boxplot()+
  geom_point(position = position_jitterdodge(jitter.width = 0, dodge.width = 0), size=5, shape=1)+
  scale_x_log10(limits=c(1e+0, 1e+3),
                breaks = scales::trans_breaks("log10", function(x) 10^x, n=4),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_discrete(limits=c("BBMD","WHO/IPCS BMD"))+
  theme_bw()+
  theme(axis.title.y = element_blank(), 
        axis.ticks.y=element_blank(),
        legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab(bquote("95"^"th"/"5"^"th"~"%ile Ratio"))+
  scale_fill_manual(breaks=c("WHO/IPCS BMD","BBMD"),
                     values=c(`WHO/IPCS BMD`="white", `BBMD`="red"))+
  annotation_logticks(sides="b")
print(bmd.degree.box)

##check log variance of BMD
#log variance is proportional to (log10(95th/5th))^2
who.bmd.log10ratio <- data.frame(chemical=who.bmd.quan.min.arrow$Chemical, CAS=who.bmd.quan.min.arrow$CAS,
                                 uncertainty=who.bmd.quan.min.arrow$degree, log10ratio=log10(who.bmd.quan.min.arrow$degree))
bbmd.log10ratio <- data.frame(chemical=bbmd.quan.min$Chemical, CAS=bbmd.quan.min$CAS,
                             uncertainty=bbmd.quan.min$degree, log10ratio=log10(bbmd.quan.min$degree))
#fold reduction in log variance
who.bmd.median <- median(who.bmd.log10ratio$log10ratio)
bbmd.median <- median(bbmd.log10ratio$log10ratio)
bmdfold <- who.bmd.median^2/bbmd.median^2
#0.84 fold reduction in log variance

#combine plots
#BMD/Reg POD
norm.bmd.plot <- ggarrange(bmd.sen.plot, bmd.sen.total, nrow=2, heights = c(0.75,0.25), align="v")
#BMD difference
bmd.diff.plot <- ggarrange(bmd.diff, bmd.diff.box, nrow=2, heights = c(0.75,0.25), align="v")
#BMD degree of uncertainty
bmd.degree.plot <- ggarrange(bmd.degree.dumb, bmd.degree.box, nrow=2, heights = c(0.75,0.25), align="v")

bmd.bind <- ggarrange(norm.bmd.plot, bmd.diff.plot, bmd.degree.plot, ncol=3, align="hv")
ggsave(bmd.bind, file="BMD plot.pdf", width = 24, height = 12, path = "Plot output", scale=0.7)
