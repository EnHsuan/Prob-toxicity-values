library(dplyr)
library(ggplot2)
library(ggpubr)

expo <- read.csv("expocast upper bound.csv")

##HDMI and expocast
who.hdmi.quan <- read.csv("WHO approximate HDMI quantile.csv")
cs.hdmi.quan <- read.csv("HDMI quantile.csv")

who.hdmi.quan$name <- "WHO/IPCS"
cs.hdmi.quan$name <- "chemical-specific"

hdmi.quan.bind <- merge(cs.hdmi.quan, who.hdmi.quan, by="CAS")
hdmi.quan.bind.df <- hdmi.quan.bind[-c(5,7,17,18,23,26),]
cs.hdmi.quan.df <- data.frame(CAS=hdmi.quan.bind.df$CAS, endpoint=hdmi.quan.bind.df$endpoint.x,
                              probRfD=hdmi.quan.bind.df$X5..x, name=hdmi.quan.bind.df$name.x)
cs.hdmi.expo <- merge(cs.hdmi.quan.df, expo, by="CAS")
cs.hdmi.expo$MOE <- with(cs.hdmi.expo, probRfD/expocast)
cs.hdmi.expo.df <- cs.hdmi.expo %>% group_by(CAS) %>% slice(which.min(probRfD))

who.hdmi.quan.df <- data.frame(CAS=hdmi.quan.bind.df$CAS, endpoint=hdmi.quan.bind.df$endpoint.x,
                               probRfD=hdmi.quan.bind.df$X5..y, name=hdmi.quan.bind.df$name.y)
who.hdmi.expo <- merge(who.hdmi.quan.df, expo, by="CAS")
who.hdmi.expo$MOE <- with(who.hdmi.expo, probRfD/expocast)
who.hdmi.expo.df <- who.hdmi.expo %>% group_by(CAS) %>% slice(which.min(probRfD))

hdmi.expo <- rbind(cs.hdmi.expo.df, who.hdmi.expo.df)
hdmi.expo$reorder <- hdmi.expo$MOE
hdmi.expo[1:19,8] <- hdmi.expo[20:38, 8]

expo.moe.plot <- ggplot()+
  geom_point(data=hdmi.expo, aes(x=MOE, y=reorder(Chemical, reorder), color=name, shape=name), size=3)+
  geom_segment(aes(x =who.hdmi.expo.df$MOE , y = who.hdmi.expo.df$Chemical, xend = cs.hdmi.expo.df$MOE, yend = cs.hdmi.expo.df$Chemical),
               arrow = arrow(length = unit(0.2, "cm")))+
  geom_rect(aes(xmin = 1e-02, xmax = 1, ymin = -Inf, ymax = Inf), fill="lightpink", alpha=0.3)+
  geom_rect(aes(xmin = 1, xmax = 100, ymin = -Inf, ymax = Inf), fill="orange", alpha=0.3)+
  geom_vline(xintercept=1, color="black", linetype="dashed")+
  geom_vline(xintercept=100, color="black", linetype="dashed")+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme_bw()+
  theme(axis.title.y = element_blank(), 
        legend.title = element_blank(), legend.position = "top", legend.text = element_text(size = 10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Margin of Exposure")+
  scale_color_manual(name="label",
                       breaks=c("WHO/IPCS","chemical-specific"),
                       labels=c("WHO/IPCS Prob. RfD/ExpoCast","Chem.-Specif. Prob. RfD/ExpoCast"),
                       values=c(`WHO/IPCS`=NA,`chemical-specific`="red"))+
  scale_shape_manual(name="label",
                       breaks=c("WHO/IPCS","chemical-specific"),
                       labels=c("WHO/IPCS Prob. RfD/ExpoCast","Chem.-Specif. Prob. RfD/ExpoCast"),
                       values=c(`WHO/IPCS`=1,`chemical-specific`=16))+
  guides(color = guide_legend(nrow = 2))+
  annotation_logticks(sides="b")
print(expo.moe.plot)

#combine HDMI, BE and exposure data
##biomonitoring data MOE
bio.data <- read.csv("bimonitoring data.csv")
be_b.quan <- read.csv("BEMI blood quantile.csv")
be_b.quan.df <- be_b.quan %>% group_by(CAS) %>% slice(which.min(X5.))

be_u.quan <- read.csv("BEMI urine quantile.csv")
be_u.quan.df <- be_u.quan %>% group_by(CAS) %>% slice(which.min(X5.))

bio.be_b <- merge(bio.data[1:4,], be_b.quan.df, by="CAS")
bio.be_u <- merge(bio.data[5:6,], be_u.quan.df, by="CAS")

bio.be_b$MOE <- with(bio.be_b, X5./biomonitoring.value)
bio.be_u$MOE <- with(bio.be_u, X5./biomonitoring.value)

bio.be_b$name <- "BE blood"
bio.be_u$name <- "BE urine"
bio.bind <- rbind(bio.be_b, bio.be_u)

cs.hdmi.moe.subset <- subset(cs.hdmi.expo.df, Chemical %in% c("Heptachlor epoxide","Endosulfan","p,p'-DDT","Dieldrin","Pentachlorophenol","2,4,5-Trichlorophenol"))

bioexpo.moe <- rbind(bio.bind[,c(1,2,9,10)], cs.hdmi.moe.subset[,c(1,5,7,4)])

bioexpo.plot <- ggplot()+
  geom_point(data=bioexpo.moe, aes(x=MOE, y=Chemical, shape=name), color="red", size=3)+
  geom_segment(aes(x =cs.hdmi.moe.subset$MOE , y = cs.hdmi.moe.subset$Chemical, xend = bio.bind$MOE, yend = bio.bind$Chemical),
               arrow = arrow(length = unit(0.2, "cm")))+
  geom_rect(aes(xmin = 1e-02, xmax = 1, ymin = -Inf, ymax = Inf), fill="lightpink", alpha=0.3)+
  geom_rect(aes(xmin = 1, xmax = 100, ymin = -Inf, ymax = Inf), fill="orange", alpha=0.3)+
  geom_vline(xintercept=1, color="black", linetype="dashed")+
  geom_vline(xintercept=100, color="black", linetype="dashed")+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme_bw()+
  theme(axis.title.y = element_blank(), 
        legend.title = element_blank(), legend.position = "top", legend.text = element_text(size = 10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Margin of Exposure")+
  scale_shape_manual(breaks=c("chemical-specific","BE blood","BE urine"),
                     labels=c("chemical-specific"="Chemical-specific Prob RfD/ExpoCast", 
                                "BE blood"="Prob BE in blood/biomonitoring data", 
                                "BE urine"="Prob BE in urine/biomonitoring data"),
                     values=c("chemical-specific"=16, "BE blood"=15, "BE urine"=17))+
  guides(shape = guide_legend(nrow = 3))+
  annotation_logticks(sides="b")
print(bioexpo.plot)

moe.plot <- ggarrange(expo.moe.plot, bioexpo.plot, align = "hv", nrow=2, ncol=1, heights = c(2,1.5))
ggsave(moe.plot, file="MOE plot (Expocast and biomonitoring).pdf", width = 6, height=12, scale=0.8, path="Plot output")
