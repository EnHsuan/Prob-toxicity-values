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
bmd.sen.bind$Chemical <- factor(bmd.sen.bind$Chemical, levels = yaxis_order)
bmd.sen.bind$name <- factor(bmd.sen.bind$name, levels = c("Bayesian BMD or WHO/IPCS (2018) BMD","WHO/IPCS (2018) BMD"))

bmd.point <- rbind(bmd.nor.sen.end.df[,c(1,3,10006)],who.bmd.nor.sen.end.df[,c(1,4,10006)])
bmd.point$name <- "Bayesian BMD or WHO/IPCS (2018) BMD"
bmd.point$name[20:38] <- "WHO/IPCS (2018) BMD"
bmd.point$Chemical <- factor(bmd.point$Chemical, levels = yaxis_order)
bmd.point$name <- factor(bmd.point$name, levels = c("Bayesian BMD or WHO/IPCS (2018) BMD","WHO/IPCS (2018) BMD"))

#set whisker to [5%,95%]
quantiles_95 <- function(x) {
  r <- quantile(x, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

bmd.sen.plot <- ggplot()+
  stat_summary(fun.data = quantiles_95, data = bmd.sen.bind, geom="boxplot", aes(x=value, y=Chemical, color=name), outlier.shape=NA, position = position_dodge(width = 0.75))+
  geom_point(data=bmd.point, aes(x=`5th_quantile`, y=Chemical, color=name), position = position_dodge(width = 0.75))+
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
  scale_color_manual(breaks=c("WHO/IPCS (2018) BMD","Bayesian BMD or WHO/IPCS (2018) BMD"),
                     values=c(`WHO/IPCS (2018) BMD`="#00BFC4", `Bayesian BMD or WHO/IPCS (2018) BMD`="#F8766D"),
                     labels=c("WHO/IPCS (2018) BMD"="WHO/IPCS BMD","Bayesian BMD or WHO/IPCS (2018) BMD"="BBMD or WHO/IPCS BMD"))+
  annotation_logticks(sides="b")
print(bmd.sen.plot)

bmd.sen.total <- ggplot(data=bmd.point, aes(x=`5th_quantile`, y=name, color=name))+
  geom_boxplot()+
  geom_point(position = position_jitter(width = 0.2))+
  geom_vline(xintercept = 1, color="black", linetype="dashed")+
  scale_x_log10(limits=c(1e-02, 1e+03),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_discrete(limits=c("Bayesian BMD or WHO/IPCS (2018) BMD","WHO/IPCS (2018) BMD"),
                   labels=c("BBMD or WHO/IPCS BMD","WHO/IPCS BMD"))+
  theme_bw()+
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 10), 
        legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab(bquote("5th"~quantile~of~BMD[M]/Regulatory~POD))+
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
  scale_fill_manual(breaks = c(F, T), values=c("gray", "gray3"), labels=c("Use WHO/IPCS BMD","Use BBMD"))+
  xlab(bquote(Delta~log[10]~BMD[M]))
print(bmd.diff)

bmd.df.m <- melt(bmd.df[,c(1,10,14)])
bmd.diff.box <- ggplot(bmd.df.m, aes(x=value, y=variable))+
  geom_boxplot()+
  geom_jitter(shape=20, position=position_jitter(0.2))+
  theme_bw()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  geom_vline(xintercept = 0, linetype="dashed")+
  xlim(-2,2)+
  xlab(bquote(Delta~log[10]~BMD[M]))
print(bmd.diff.box)

#degree of uncertainty
chem <- read.csv("regulatory reference dose.csv")[,1:2]
bmd.quan <- read.csv("BMD quantile.csv")
who.bmd.quan <- read.csv("WHO BMD quantile.csv")

bmd.quan.chem <- merge(bmd.quan, chem, by="CAS")
who.bmd.quan.chem <- merge(who.bmd.quan, chem, by="CAS")

bmd.quan.min <- bmd.quan.chem %>% group_by(CAS) %>% slice(which.min(`X5.`))
who.bmd.quan.min <- who.bmd.quan.chem %>% group_by(CAS) %>% slice(which.min(`X5.`))

bmd.quan.min$degree <- with(bmd.quan.min, X95./X5.)
bmd.quan.min$name <- "BBMD or WHO/IPCS BMD"
who.bmd.quan.min$degree <- with(who.bmd.quan.min, X95./X5.)
who.bmd.quan.min$name <- "WHO/IPCS BMD"

bmd.degree <- rbind(who.bmd.quan.min, bmd.quan.min)
bmd.degree$Chemical <- factor(bmd.degree$Chemical, levels = yaxis_order)
bmd.degree$name <- factor(bmd.degree$name, levels = c("WHO/IPCS BMD","BBMD or WHO/IPCS BMD"))

bmd.degree.dumb <- ggplot()+
  geom_point(data=bmd.degree, aes(x=degree, y=Chemical, color=name), size=3)+
  geom_segment(aes(x =who.bmd.quan.min$degree , y = who.bmd.quan.min$Chemical, xend = bmd.quan.min$degree, yend = bmd.quan.min$Chemical),
               arrow = arrow(length = unit(0.2, "cm")))+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme_bw()+
  theme(axis.title.y = element_blank(), 
        legend.title = element_blank(), legend.position = "top", legend.text = element_text(size = 10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("95th / 5th ratio")+
  scale_color_manual(breaks=c("WHO/IPCS BMD","BBMD or WHO/IPCS BMD"),
                     values=c(`WHO/IPCS BMD`="#00BFC4", `BBMD or WHO/IPCS BMD`="#F8766D"))+
  annotation_logticks(sides="b")
print(bmd.degree.dumb)

bmd.degree.box <- ggplot(data=bmd.degree, aes(x=degree, y=name, color=name))+
  geom_boxplot()+
  geom_jitter(shape=20, position=position_jitter(0.2))+
  scale_x_log10(limits=c(1e+0, 1e+3),
                breaks = scales::trans_breaks("log10", function(x) 10^x, n=4),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_discrete(limits=c("BBMD or WHO/IPCS BMD","WHO/IPCS BMD"))+
  theme_bw()+
  theme(axis.title.y = element_blank(), 
        axis.ticks.y=element_blank(),
        legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("95th / 5th ratio")+
  scale_color_manual(breaks=c("WHO/IPCS BMD","BBMD or WHO/IPCS BMD"),
                     values=c(`WHO/IPCS BMD`="#00BFC4", `BBMD or WHO/IPCS BMD`="#F8766D"))+
  annotation_logticks(sides="b")
print(bmd.degree.box)

#combine plots
#BMD/Reg POD
norm.bmd.plot <- ggarrange(bmd.sen.plot, bmd.sen.total, nrow=2, heights = c(0.75,0.25), align="v")
#BMD difference
bmd.diff.plot <- ggarrange(bmd.diff, bmd.diff.box, nrow=2, heights = c(0.75,0.25), align="v")
#BMD degree of uncertainty
bmd.degree.plot <- ggarrange(bmd.degree.dumb, bmd.degree.box, nrow=2, heights = c(0.75,0.25), align="v")

bmd.bind <- ggarrange(norm.bmd.plot, bmd.diff.plot, bmd.degree.plot, ncol=3, align="hv")
ggsave(bmd.bind, file="BMD plot.pdf", width = 24, height = 12, path = "HDMI plots", scale=0.7)

###-------------------------------------------------------------------------------------------------
#generate csv file for Prism
#A panel normalized distributions
bmd.nor.sen.end.df$Chemical <- factor(bmd.nor.sen.end.df$Chemical, levels = rev(yaxis_order))
bmd.nor.sen.end.df_sort <- bmd.nor.sen.end.df[order(bmd.nor.sen.end.df$Chemical),]
bmd.nor.sen.end.df.t <- data.frame(t(bmd.nor.sen.end.df_sort[,3:10003]))
colnames(bmd.nor.sen.end.df.t) <- bmd.nor.sen.end.df.t[1,]
bmd.data <- bmd.nor.sen.end.df.t[-1,]
write.csv(bmd.data, file="prism file/Fig 2 (A) BMD normalized data.csv", row.names = FALSE)

who.bmd.nor.sen.end.df$Chemical <- factor(who.bmd.nor.sen.end.df$Chemical, levels = rev(yaxis_order))
who.bmd.nor.sen.end.df_sort <- who.bmd.nor.sen.end.df[order(who.bmd.nor.sen.end.df$Chemical),]
who.bmd.nor.sen.end.df.t <- data.frame(t(who.bmd.nor.sen.end.df_sort[,4:10004]))
colnames(who.bmd.nor.sen.end.df.t) <- who.bmd.nor.sen.end.df.t[1,]
who.bmd.data <- who.bmd.nor.sen.end.df.t[-1,]
write.csv(who.bmd.data, file="prism file/Fig 2 (A) WHO BMD normalized data.csv", row.names = FALSE)

bmd.nor.5th <- bmd.nor.sen.end.df[,c(3,10006)]
bmd.nor.5th$name <- "BBMD or WHO/IPCS BMD"
who.bmd.nor.5th <- who.bmd.nor.sen.end.df[,c(4, 10006)]
who.bmd.nor.5th$name <- "WHO/IPCS BMD"
bmd.5th.bind <- merge(bmd.nor.5th, who.bmd.nor.5th, by="Chemical")
bmd.5th.bind$Chemical <- factor(bmd.5th.bind$Chemical, levels = rev(yaxis_order))
bmd.5th.bind_sort <- bmd.5th.bind[order(bmd.5th.bind$Chemical),]
write.csv(bmd.5th.bind_sort, file="prism file/Fig 2 (A) BMD normalized 5th.csv", row.names = FALSE)

#B panel for differences
bmd.diff.df <- bmd.df[,c(10,11,14)]
bmd.diff.df$Chemical <- factor(bmd.diff.df$Chemical, levels = rev(yaxis_order))
bmd.diff.df_sort <- bmd.diff.df[order(bmd.diff.df$Chemical),]
write.csv(bmd.diff.df_sort, file="prism file/Fig 2 (B) BMD differences.csv", row.names = FALSE)

#C panel for uncertainty degree
bmd.degree.df <- merge(bmd.degree[1:19,6:8], bmd.degree[20:38,6:8], by="Chemical")
bmd.degree.df$Chemical <- factor(bmd.degree.df$Chemical, levels = rev(yaxis_order))
bmd.degree.df_sort <- bmd.degree.df[order(bmd.degree.df$Chemical),]
write.csv(bmd.degree.df_sort, file="prism file/Fig 2 (C) BMD uncertainty degree.csv", row.names = FALSE)
