library(dplyr)
library(ggplot2)
library(ggpubr)

expo <- read.csv("expocast upper bound.csv")

##HDMI uncertainty
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

moe.plot <- ggplot()+
  geom_point(data=hdmi.expo, aes(x=MOE, y=reorder(Chemical, reorder), color=name, shape=name), size=3)+
  geom_point(data=bio.bind, aes(x=MOE, y=Chemical, color=name, shape=name), size=3)+
  geom_segment(aes(x =who.hdmi.expo.df$MOE , y = who.hdmi.expo.df$Chemical, xend = cs.hdmi.expo.df$MOE, yend = cs.hdmi.expo.df$Chemical),
               arrow = arrow(length = unit(0.2, "cm")))+
  geom_rect(aes(xmin = 1e-02, xmax = 1, ymin = -Inf, ymax = Inf), color="gray", alpha=0.3)+
  geom_vline(xintercept=1, color="black", linetype="dashed")+
  geom_vline(xintercept=100, color="black", linetype="dashed")+
  scale_x_log10(limits=c(1e-02,1e+06),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme_bw()+
  theme(axis.title.y = element_blank(), 
        legend.title = element_blank(), legend.position = "top", legend.text = element_text(size = 10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Margin of Exposure (MOE)")+
  scale_color_manual(name="MOE",
                     breaks=c("WHO/IPCS","chemical-specific","BE blood","BE urine"),
                     labels=c("WHO/IPCS Prob RfD/ExpoCast","Chemical-specific Prob RfD/ExpoCast", 
                              "Chemical-specific Prob BE in blood/biomonitoring data", 
                              "Chemical-specific Prob BE in urine/biomonitoring data"),
                     values=c(`WHO/IPCS`="#56B4E9",`chemical-specific`="#D55E00",`BE blood`="#D55E00",`BE urine`="#D55E00"))+
  scale_shape_manual(name="MOE",
                     breaks=c("WHO/IPCS","chemical-specific","BE blood","BE urine"),
                     labels=c("WHO/IPCS Prob RfD/ExpoCast","Chemical-specific Prob RfD/ExpoCast", 
                              "Chemical-specific Prob BE in blood/biomonitoring data", 
                              "Chemical-specific Prob BE in urine/biomonitoring data"),
                     values=c(`WHO/IPCS`=16,`chemical-specific`=16,`BE blood`=15,`BE urine`=17))+
  guides(color = guide_legend(nrow = 2))+# Split the legend into 2 rows
  annotation_logticks(sides="b")
print(moe.plot)
ggsave(moe.plot, file="MOE plot.pdf", width=10, height=8, path="Plot output", scale=0.8)