##BMD and HDMI differences
library(dplyr)
library(ggplot2)
library(ggpubr)
library(reshape2)

chem <- read.csv("regulatory reference dose.csv")[,1:2]
bbmd <- read.csv("BBMD quantile.csv")[,4:5]
bbmd.chem <- data.frame(CAS=unique(bbmd$CAS))
bbmd.chem$BBMD <- "BBMD"

#use the most sensitive endpoints
#HDMI
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

hdmi.diff <- ggplot(hdmi.df, aes(x=log10.cs_who, y=reorder(Chemical, log10.cs_who)))+
  geom_bar(aes(fill = BBMD == BBMD), stat="identity")+
  theme_bw()+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.y = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(5.5, 5.5, 0, 5.5),"pt"))+
  geom_vline(xintercept = 0)+
  xlim(-2,2)+
  scale_fill_manual(guide = "none", breaks = c(FALSE, TRUE), values=c("gray", "gray3"))
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
  xlab(bquote(Delta~HD[M]^I))
print(hdmi.diff.box)

#BMD
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
bmd.df$order <- hdmi.df$log10.cs_who

bmd.diff <- ggplot(bmd.df, aes(x=log10.bmd_whobmd, y=reorder(Chemical, order)))+
  geom_bar(aes(fill = BBMD == BBMD), stat="identity")+
  theme_bw()+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.y = element_blank(), axis.text.y = element_blank(), 
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.margin = unit(c(5.5, 5.5, 0, 0),"pt"))+
  geom_vline(xintercept = 0)+
  xlim(-2,2)+
  scale_fill_manual(guide = "none", breaks = c(FALSE, TRUE), values=c("gray", "gray3"))
print(bmd.diff)

bmd.df.m <- melt(bmd.df[,c(1,10,14)])
bmd.diff.box <- ggplot(bmd.df.m, aes(x=value, y=variable))+
  geom_boxplot()+
  geom_jitter(shape=20, position=position_jitter(0.2))+
  theme_bw()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.margin = unit(c(0, 5.5, 5.5, 0),"pt"))+
  geom_vline(xintercept = 0, linetype="dashed")+
  xlim(-2,2)+
  xlab(bquote(Delta~BMD["M"]^"I"))
print(bmd.diff.box)

#HDMI-BMD
hdmi.bmd <- merge(bmd.df[,c(1,10,11,14)], hdmi.df[,c(1,14)], by="CAS")
hdmi.bmd$hdmi_bmd <- hdmi.bmd$log10.cs_who-hdmi.bmd$log10.bmd_whobmd
hdmi.bmd$order <- hdmi.df$log10.cs_who

hdmi.bmd.diff <- ggplot(hdmi.bmd, aes(x=hdmi_bmd, y=reorder(Chemical, order)))+
  geom_bar(aes(fill = BBMD == BBMD),  stat="identity")+
  theme_bw()+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.y = element_blank(), axis.text.y = element_blank(), 
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.margin = unit(c(5.5, 5.5, 0, 0),"pt"))+
  geom_vline(xintercept = 0)+
  xlim(-2,2)+
  scale_fill_manual(guide = "none", breaks = c(FALSE, TRUE), values=c("gray", "gray3"))
print(hdmi.bmd.diff)

hdmi.bmd.m <- melt(hdmi.bmd[,c(1,2,6)])
hdmi.bmd.diff.box <- ggplot(hdmi.bmd.m, aes(x=value, y=variable))+
  geom_boxplot()+
  geom_jitter(shape=20, position=position_jitter(0.2))+
  theme_bw()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.margin = unit(c(0, 5.5, 5.5, 0),"pt"))+
  geom_vline(xintercept = 0, linetype="dashed")+
  xlim(-2,2)+
  xlab(bquote(Delta~AF[intra]))
print(hdmi.bmd.diff.box)

hdmi.bind <- ggarrange(hdmi.diff, hdmi.diff.box, nrow=2, heights = c(0.75,0.25), align="v")
bmd.bind <- ggarrange(bmd.diff, bmd.diff.box, nrow=2, heights = c(0.75,0.25), align="v")
hdmi.bmd.bind <- ggarrange(hdmi.bmd.diff, hdmi.bmd.diff.box, nrow=2, heights = c(0.75,0.25), align="v")
diff.bind <- ggarrange(hdmi.bind, bmd.bind, hdmi.bmd.bind, ncol=3, align="hv", widths = c(1.2/3,0.9/3,0.9/3))
ggsave(diff.bind, file="Differences in HDMI BMD UFh bar plot.pdf", width=12, height=6, path="Plot output")

#check which one contributes more to the shift in the overall distribution
#use absolute value to compare the difference
abs.df <- merge(bmd.df[,c(1,10,14)], hdmi.bmd[,c(2,5)], by="Chemical")
abs.df$abs.bmd <- abs(abs.df$log10.bmd_whobmd)
abs.df$abs.ufh <- abs(abs.df$hdmi_bmd)
abs.df$diff <- abs.df$abs.bmd-abs.df$abs.ufh
#BMD is smaller:  15 chemicals (total 19 chemicals)