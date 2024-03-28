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

cs.ufh <- read.csv("TKTDVF01.csv")
who.ufh <- read.csv("WHO approximate intraspecies factor.csv")

###normalized UFh between two approaches
###the most sensitive endpoint
#normalized by an uncertainty factor of 10

ufh.chem <- left_join(cs.ufh, rfd[,c(1:2)], by="CAS")
ufh.nor <- with(ufh.chem, ufh.chem[,3:10002]/10)
ufh.nor.df <- cbind(ufh.chem[,c(1,2,10003)], ufh.nor)
ufh.nor.df$`95th_quantile` <- apply(ufh.nor.df[,4:10003], 1, quantile, probs=0.95)
ufh.nor.sen.df <- ufh.nor.df %>% group_by(CAS) %>% slice(which.max(`95th_quantile`))
ufh.nor.sen.m <- melt(ufh.nor.sen.df, id=c(1:3,10004))
ufh.nor.sen.m$name <- "Chemical-specific UFh"

who.ufh.rep <- as.data.frame(apply(who.ufh, 2, rep, times=19))
who.ufh.chem <- cbind(rfd[,c(1:2)], who.ufh.rep)
who.ufh.nor <- with(who.ufh.chem, who.ufh.chem[,3:10002]/10)
who.ufh.nor.df <- cbind(rfd[,c(1:2)], who.ufh.nor)
who.ufh.nor.df$`95th_quantile` <- apply(who.ufh.nor.df[,3:10002], 1, quantile, probs=0.95)
who.ufh.nor.m <- melt(who.ufh.nor.df, id=c(1,2,10003))
who.ufh.nor.m$name <- "WHO/IPCS UFh"

ufh.sen.bind <- rbind(ufh.nor.sen.m[,-2], who.ufh.nor.m)
ufh.sen.bind$Chemical <- factor(ufh.sen.bind$Chemical, levels = yaxis_order)
ufh.sen.bind$name <- factor(ufh.sen.bind$name, levels = c("Chemical-specific UFh","WHO/IPCS UFh"))

ufh.point <- rbind(ufh.nor.sen.df[,c(1,3,10004)],who.ufh.nor.df[,c(1,2,10003)])
ufh.point$name <- "Chemical-specific UFh"
ufh.point$name[20:38] <- "WHO/IPCS UFh"
ufh.point$Chemical <- factor(ufh.point$Chemical, levels = yaxis_order)
ufh.point$name <- factor(ufh.point$name, levels = c("Chemical-specific UFh","WHO/IPCS UFh"))

#set whisker to [5%,95%]
quantiles_95 <- function(x) {
  r <- quantile(x, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

ufh.sen.plot <- ggplot()+
  stat_summary(fun.data = quantiles_95, data = ufh.sen.bind, geom="boxplot", aes(x=value, y=Chemical, color=name), outlier.shape=NA, position = position_dodge(width = 0.75))+
  geom_point(data=ufh.point, aes(x=`95th_quantile`, y=Chemical, color=name), position = position_dodge(width = 0.75))+
  geom_vline(xintercept = 1, color="black", linetype="dashed")+
  scale_x_log10(limits=c(1e-1, 1e+02),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme_bw()+
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 10), 
        legend.text = element_text(size = 10), legend.title = element_blank(),
        legend.position = "top", legend.key.size = unit(3, 'mm'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab(bquote(AF[intra]/10))+
  scale_color_manual(breaks=c("WHO/IPCS UFh","Chemical-specific UFh"),
                     values=c(`WHO/IPCS UFh`="#00BFC4", `Chemical-specific UFh`="#F8766D"),
                     labels=c("WHO/IPCS UFh"=bquote(WHO/IPCS~AF[intra]),"Chemical-specific UFh"=bquote(Chemical-specfiic~AF[intra])))
print(ufh.sen.plot)

ufh.sen.total <- ggplot(data=ufh.point, aes(x=`95th_quantile`, y=name, color=name))+
  geom_boxplot()+
  geom_point(position = position_jitter(width = 0.2))+
  geom_vline(xintercept = 1, color="black", linetype="dashed")+
  scale_x_log10(limits=c(1e-1, 1e+02),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_discrete(limits=c("Chemical-specific UFh","WHO/IPCS UFh"),
                   labels=c(bquote(Chemical-specfiic~AF[intra]), bquote(WHO/IPCS~AF[intra])))+
  theme_bw()+
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 10), 
        legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab(bquote(AF[intra]/10))
print(ufh.sen.total)

#UFh differences
chem <- read.csv("regulatory reference dose.csv")[,1:2]
bbmd <- read.csv("BBMD quantile.csv")[,4:5]
bbmd.chem <- data.frame(CAS=unique(bbmd$CAS))
bbmd.chem$BBMD <- "BBMD"

ufh.quan <- read.csv("TKTDVF01 quantile.csv")
ufh.quan.chem <- merge(ufh.quan, chem, by="CAS")
ufh.quan.19 <- ufh.quan.chem %>% group_by(CAS) %>% slice(which.max(`X95.`))
who.ufh.quan <- read.csv("WHO approximate intraspecies factor quantile.csv")
who.ufh.quan.rep <- as.data.frame(apply(who.ufh.quan, 2, rep, times=19))

ufh.bind <- cbind(ufh.quan.19, who.ufh.quan.rep)
ufh.df <- left_join(ufh.bind, bbmd.chem, by="CAS")
ufh.df <- ufh.df %>% arrange(CAS)
ufh.df$log10.csufh50 <- log10(ufh.df$X50....3)
ufh.df$log10.whoufh50 <- log10(ufh.df$X50....7)
ufh.df$log10.cs_who <- ufh.df$log10.csufh50-ufh.df$log10.whoufh50
ufh.df$Chemical <- factor(ufh.df$Chemical, levels = yaxis_order)

ufh.diff <- ggplot(ufh.df, aes(x=log10.cs_who, y=Chemical))+
  geom_bar(stat="identity")+
  theme_bw()+
  theme(#axis.title.x = element_blank(), #axis.text.x = element_blank(), 
    axis.title.y = element_blank(), #axis.text.y = element_blank(), 
    legend.text = element_text(size = 10), legend.title = element_blank(),
    legend.position = "top", legend.key.size = unit(3, 'mm'),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  geom_vline(xintercept = 0)+
  xlim(-2,2)+
  #scale_fill_manual(breaks = c(F, T), values=c("gray", "gray3"), labels=c("Use WHO/IPCS BMD","Use BBMD"))+
  xlab(bquote(Delta~log[10]~AF[intra]))
print(ufh.diff)

ufh.df.m <- melt(ufh.df[,c(1,6,13)])
ufh.diff.box <- ggplot(ufh.df.m, aes(x=value, y=variable))+
  geom_boxplot()+
  geom_jitter(shape=20, position=position_jitter(0.2))+
  theme_bw()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  geom_vline(xintercept = 0, linetype="dashed")+
  xlim(-2,2)+
  xlab(bquote(Delta~log[10]~AF[intra]))
print(ufh.diff.box)

#degree of uncertainty
chem <- read.csv("regulatory reference dose.csv")[,1:2]
ufh.quan <- read.csv("TKTDVF01 quantile.csv")
ufh.quan.19 <- ufh.quan %>% group_by(CAS) %>% slice(which.max(`X95.`))
who.ufh.quan <- read.csv("WHO approximate intraspecies factor quantile.csv")
who.ufh.quan.rep <- as.data.frame(apply(who.ufh.quan, 2, rep, times=19))

ufh.quan.chem <- merge(ufh.quan.19, chem, by="CAS")
ufh.quan.chem <- ufh.quan.chem[order(ufh.quan.chem$CAS),]
who.ufh.quan.chem <- cbind(who.ufh.quan.rep, chem)
who.ufh.quan.chem <- who.ufh.quan.chem[order(who.ufh.quan.chem$CAS),]

who.ufh.quan.chem$degree <- with(who.ufh.quan.chem, X95./X5.)
who.ufh.quan.chem$name <- "WHO/IPCS UFh"
ufh.quan.chem$degree <- with(ufh.quan.chem, X95./X5.)
ufh.quan.chem$name <- "Chemical-specific UFh"

ufh.degree <- rbind(who.ufh.quan.chem, ufh.quan.chem[,-2])
ufh.degree$Chemical <- factor(ufh.degree$Chemical, levels = yaxis_order)
ufh.degree$name <- factor(ufh.degree$name, levels = c("WHO/IPCS UFh","Chemical-specific UFh"))

ufh.degree.dumb <- ggplot()+
  geom_point(data=ufh.degree, aes(x=degree, y=Chemical, color=name), size=3)+
  geom_segment(aes(x =who.ufh.quan.chem$degree , y = who.ufh.quan.chem$Chemical, xend = ufh.quan.chem$degree, yend = ufh.quan.chem$Chemical),
               arrow = arrow(length = unit(0.2, "cm")))+
  scale_x_log10(limits=c(1, 1e+02), 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme_bw()+
  theme(axis.title.y = element_blank(), 
        legend.title = element_blank(), legend.position = "top", legend.text = element_text(size = 10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("95th / 5th ratio")+
  scale_color_manual(breaks=c("WHO/IPCS UFh","Chemical-specific UFh"),
                     values=c(`WHO/IPCS UFh`="#00BFC4", `Chemical-specific UFh`="#F8766D"),
                     labels=c(bquote(WHO/IPCS~AF[intra]), bquote(Chemical-specfiic~AF[intra])))+
  #guides(color = guide_legend(nrow = 2))+
  annotation_logticks(sides="b")
print(ufh.degree.dumb)

ufh.degree.box <- ggplot(data=ufh.degree, aes(x=degree, y=name, color=name))+
  geom_boxplot()+
  geom_jitter(shape=20, position=position_jitter(0.2))+
  scale_x_log10(limits=c(1, 1e+02), 
                breaks = scales::trans_breaks("log10", function(x) 10^x, n=4),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_discrete(limits=c("Chemical-specific UFh","WHO/IPCS UFh"),
                   labels=c(bquote(Chemical-specfiic~AF[intra]), bquote(WHO/IPCS~AF[intra])))+
  theme_bw()+
  theme(axis.title.y = element_blank(), 
        axis.ticks.y=element_blank(),
        legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("95th / 5th ratio")+
  scale_color_manual(breaks=c("WHO/IPCS UFh","Chemical-specific UFh"),
                     values=c(`WHO/IPCS UFh`="#00BFC4", `Chemical-specific UFh`="#F8766D"),
                     labels=c(bquote(WHO/IPCS~AF[intra]), bquote(Chemical-specfiic~AF[intra])))+
  #guides(color = guide_legend(nrow = 2))+
  annotation_logticks(sides="b")
print(ufh.degree.box)

##check log variance of AFintra
#log variance is proportional to (log10(95th/5th))^2
who.ufh.log10ratio <- data.frame(chemical=who.ufh.quan.chem$Chemical, CAS=who.ufh.quan.chem$CAS,
                                 uncertainty=who.ufh.quan.chem$degree, log10ratio=log10(who.ufh.quan.chem$degree))
ufh.log10ratio <- data.frame(chemical=ufh.quan.chem$Chemical, CAS=ufh.quan.chem$CAS,
                             uncertainty=ufh.quan.chem$degree, log10ratio=log10(ufh.quan.chem$degree))
#fold reduction in log variance
who.ufh.median <- median(who.ufh.log10ratio$log10ratio)
ufh.median <- median(ufh.log10ratio$log10ratio)
ufhfold <- who.ufh.median^2/ufh.median^2
#13.54 fold reduction in log variance

#combine plots
#UFh/10
norm.ufh.plot <- ggarrange(ufh.sen.plot, ufh.sen.total, nrow=2, heights = c(0.75,0.25), align="v")
#UFh difference
ufh.diff.plot <- ggarrange(ufh.diff, ufh.diff.box, nrow=2, heights = c(0.75,0.25), align="v")
#UFh degree of uncertainty
ufh.degree.plot <- ggarrange(ufh.degree.dumb, ufh.degree.box, nrow=2, heights = c(0.75,0.25), align="v")

ufh.bind <- ggarrange(norm.ufh.plot, ufh.diff.plot, ufh.degree.plot, ncol=3, align="hv")
ggsave(ufh.bind, file="Intraspecies uncertainty plot.pdf", width = 24, height = 12, path = "HDMI plots", scale=0.7)

###-------------------------------------------------------------------------------------------------
#generate csv file for Prism
#A panel normalized distributions
ufh.nor.sen.df$Chemical <- factor(ufh.nor.sen.df$Chemical, levels = yaxis_order)
ufh.nor.sen.df_sort <- ufh.nor.sen.df[order(ufh.nor.sen.df$Chemical),]
ufh.nor.sen.df.t <- data.frame(t(ufh.nor.sen.df_sort[,3:10003]))
colnames(ufh.nor.sen.df.t) <- ufh.nor.sen.df.t[1,]
ufh.data <- ufh.nor.sen.df.t[-1,]
write.csv(ufh.data, file="prism file/Fig 3 (A) chemical-specific AFintra normalized data.csv", row.names = FALSE)

who.ufh.nor.df$Chemical <- factor(who.ufh.nor.df$Chemical, levels = yaxis_order)
who.ufh.nor.df_sort <- who.ufh.nor.df[order(who.ufh.nor.df$Chemical),]
who.ufh.nor.df.t <- data.frame(t(who.ufh.nor.df_sort[,c(1,3:10002)]))
colnames(who.ufh.nor.df.t) <- who.ufh.nor.df.t[1,]
who.ufh.data <- who.ufh.nor.df.t[-1,]
write.csv(who.ufh.data, file="prism file/Fig 3 (A) WHO AFintra normalized data.csv", row.names = FALSE)

ufh.nor.95th <- ufh.nor.sen.df[,c(3,10004)]
ufh.nor.95th$name <- "Chemical-specific intra"
who.ufh.nor.95th <- who.ufh.nor.df[,c(1, 10003)]
who.ufh.nor.95th$name <- "WHO/IPCS intra"
ufh.95th.bind <- merge(ufh.nor.95th, who.ufh.nor.95th, by="Chemical")
ufh.95th.bind$Chemical <- factor(ufh.95th.bind$Chemical, levels = yaxis_order)
ufh.95th.bind_sort <- ufh.95th.bind[order(ufh.95th.bind$Chemical),]
write.csv(ufh.95th.bind_sort, file="prism file/Fig 3 (D) AFintra normalized 95th.csv", row.names = FALSE)

#B panel for differences
ufh.diff.df <- ufh.df[,c(6,10,13)]
ufh.diff.df$Chemical <- factor(ufh.diff.df$Chemical, levels = yaxis_order)
ufh.diff.df_sort <- ufh.diff.df[order(ufh.diff.df$Chemical),]
write.csv(ufh.diff.df_sort, file="prism file/Fig 3 (B)(E) AFintra differences.csv", row.names = FALSE)

#C panel for uncertainty degree
ufh.degree.df <- merge(ufh.degree[1:19,c(4,6,7)], ufh.degree[20:38,c(4,6,7)], by="Chemical")
ufh.degree.df$Chemical <- factor(ufh.degree.df$Chemical, levels = yaxis_order)
ufh.degree.df_sort <- ufh.degree.df[order(ufh.degree.df$Chemical),]
write.csv(ufh.degree.df_sort, file="prism file/Fig 3 (C)(F) AFintra uncertainty degree.csv", row.names = FALSE)
