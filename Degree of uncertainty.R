library(dplyr)
library(ggplot2)
library(ggpubr)
library(reshape2)

##BMD uncertainty
#only use the chemicals and endpoints with available dose-response data
bbmd.quan <- read.csv("BBMD quantile.csv")
who.bmd.quan <- read.csv("WHO BMD quantile.csv")[-c(1:3,7,8,11,12,15,16,18,20,21),]

who.bmd.quan$`95th_5th` <- with(who.bmd.quan, X95./X5.)
who.bmd.quan$name <- "WHO/IPCS BMD"
bbmd.quan$`95th_5th` <- with(bbmd.quan, X95./X5.)
bbmd.quan$name <- "BBMD"

bmd.quan <- left_join(bbmd.quan, who.bmd.quan, by="CAS")
bmd.quan.df <- bmd.quan[-c(13,16,18),] #23 endpoints
who.bmd.bind <- data.frame(endpoint=bmd.quan.df$endpoint.x, CAS=bmd.quan.df$CAS,
                            uncertainty=bmd.quan.df$`95th_5th.y`, name=bmd.quan.df$name.y)
bbmd.bind <- data.frame(endpoint=bmd.quan.df$endpoint.x, CAS=bmd.quan.df$CAS,
                           uncertainty=bmd.quan.df$`95th_5th.x`, name=bmd.quan.df$name.x)
bmd.bind <- rbind(who.bmd.bind, bbmd.bind)
bmd.bind$order <- bmd.bind$uncertainty
bmd.bind[24:46,5] <- bmd.bind[1:23,5]

##check log variance of BMD
#log variance is proportional to (log10(95th/5th))^2
who.bmd.log10ratio <- data.frame(endpoint=bmd.quan.df$endpoint.x, CAS=bmd.quan.df$CAS,
                                  uncertainty=bmd.quan.df$`95th_5th.y`, log10ratio=log10(bmd.quan.df$`95th_5th.y`))
bbmd.log10ratio <- data.frame(endpoint=bmd.quan.df$endpoint.x, CAS=bmd.quan.df$CAS,
                                 uncertainty=bmd.quan.df$`95th_5th.x`, log10ratio=log10(bmd.quan.df$`95th_5th.x`))
#fold reduction in log variance
who.bmd.median <- median(who.bmd.log10ratio$log10ratio)
bbmd.median <- median(bbmd.log10ratio$log10ratio)
bmdfold <- who.bmd.median^2/bbmd.median^2
#0.99 fold reduction in log variance

##HDMI uncertainty
who.hdmi.quan <- read.csv("WHO approximate HDMI quantile.csv")
cs.hdmi.quan <- read.csv("HDMI quantile.csv")

who.hdmi.quan$`95th_5th` <- with(who.hdmi.quan, X95./X5.)
who.hdmi.quan$name <- "WHO"
cs.hdmi.quan$`95th_5th` <- with(cs.hdmi.quan, X95./X5.)
cs.hdmi.quan$name <- "chemical-specific"

hdmi.quan <- left_join(cs.hdmi.quan, who.hdmi.quan, by="CAS")
hdmi.quan.df <- hdmi.quan[-c(2,5,7,17,18,24,25),] #35 endpoints
who.hdmi.bind <- data.frame(endpoint=hdmi.quan.df$endpoint.x, CAS=hdmi.quan.df$CAS,
                            uncertainty=hdmi.quan.df$`95th_5th.y`, name=hdmi.quan.df$name.y)
cs.hdmi.bind <- data.frame(endpoint=hdmi.quan.df$endpoint.x, CAS=hdmi.quan.df$CAS,
                            uncertainty=hdmi.quan.df$`95th_5th.x`, name=hdmi.quan.df$name.x)
hdmi.bind <- rbind(who.hdmi.bind, cs.hdmi.bind)
hdmi.bind$order <- hdmi.bind$uncertainty
hdmi.bind[36:70,5] <- hdmi.bind[1:35,5]

##check log variance
#log variance is proportional to (log10(95th/5th))^2
who.hdmi.log10ratio <- data.frame(endpoint=hdmi.quan.df$endpoint.x, CAS=hdmi.quan.df$CAS,
                                  uncertainty=hdmi.quan.df$`95th_5th.y`, log10ratio=log10(hdmi.quan.df$`95th_5th.y`))
cs.hdmi.log10ratio <- data.frame(endpoint=hdmi.quan.df$endpoint.x, CAS=hdmi.quan.df$CAS,
                                 uncertainty=hdmi.quan.df$`95th_5th.x`, log10ratio=log10(hdmi.quan.df$`95th_5th.x`))
#fold reduction in log variance
who.median <- median(who.hdmi.log10ratio$log10ratio)
cs.median <- median(cs.hdmi.log10ratio$log10ratio)
fold <- who.median^2/cs.median^2
#1.26 fold reduction in log variance

#violin plot
bmd.violin.plot <- ggplot(data=bmd.bind)+
  geom_violin(aes(x=name, y=uncertainty, color=name), draw_quantiles = c(0.05, 0.5, 0.95))+
  geom_point(aes(x=name, y=uncertainty, group=endpoint, color=name))+
  geom_line(aes(x=name, y=uncertainty, group=endpoint), size=0.05)+
  scale_y_log10(limits=c(1e+0, 1e+3),
                breaks = scales::trans_breaks("log10", function(x) 10^x, n=4),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme_bw()+
  theme(axis.title.x = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="none")+
  scale_x_discrete(limits=c("WHO/IPCS BMD","BBMD"))+
  ylab(bquote(95^th/5^th~ratio))+
  annotation_logticks(sides="l")
print(bmd.violin.plot)

hdmi.violin.plot <- ggplot(data=hdmi.bind)+
  geom_violin(aes(x=name, y=uncertainty, color=name), draw_quantiles = c(0.05, 0.5, 0.95))+
  geom_point(aes(x=name, y=uncertainty, group=endpoint, color=name))+
  geom_line(aes(x=name, y=uncertainty, group=endpoint), size=0.05)+
  scale_y_log10(limits=c(1e+0, 1e+3),
                #breaks = scales::trans_breaks("log10", function(x) 10^x),
                #labels = scales::trans_format("log10", scales::math_format(10^.x))
                )+
  theme_bw()+
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), axis.text.y = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="none",
        plot.margin = unit(c(5.5, 5.5, 5.5, 0),"pt"))+
  scale_x_discrete(limits=c("WHO","chemical-specific"), labels=c(bquote(WHO/IPCS~HD[M]^I), bquote(Chemical-specific~HD[M]^I)))+
  ylab(bquote(95^th/5^th~ratio))+
  annotation_logticks(sides="l")
print(hdmi.violin.plot)

violin.plot <- ggarrange(bmd.violin.plot, hdmi.violin.plot, ncol=2, align="h", widths = c(0.54,0.46))
ggsave(violin.plot, file="Degree of uncertainty uncertainty violin plot.pdf", 
       width=10, height=6, path="Plot output", scale=0.7)
