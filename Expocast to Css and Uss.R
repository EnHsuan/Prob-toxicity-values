library(dplyr)
expo.95 <- read.csv("expocast median and upper bound.csv") %>% arrange(CAS)
css.df <- read.csv("Css samples CAS by row.csv") %>% arrange(CAS) #unit: mg/L
gfr_fub.df <- read.csv("GFRxFub samples CAS by row.csv") %>% arrange(CAS) #unit: L/kg/day

set.seed(1234)

css <- data.frame(matrix(NA, nrow=2, ncol=19))
uss <- data.frame(matrix(NA, nrow=2, ncol=19))
for (i in 1:19){
  #Css
  css.temp <- as.numeric(css.df[i, 1:10000])
  css.50 <- quantile(css.temp, probs = 0.5, na.rm=TRUE) #Css median
  css[1,i] <- expo.95$median[i]*css.50
  css[2,i] <- expo.95$upper[i]*css.50
  
  #Uss
  gfr_fub.temp <- as.numeric(gfr_fub.df[i, 1:10000])
  gfr_fub.50 <- quantile(gfr_fub.temp, probs = 0.5, na.rm=TRUE) #gfr_fub median
  uss[1,i] <- css[1,i]*gfr_fub.50
  uss[2,i] <- css[2,i]*gfr_fub.50
}

colnames(css) <- css.df$CAS
colnames(uss) <- gfr_fub.df$CAS
rownames(css) <- c("50th","95th")
rownames(uss) <- c("50th","95th")

css.t <- data.frame(css=t(css))
css.t$CAS <- rownames(css.t)
uss.t <- data.frame(uss=t(uss))
uss.t$CAS <- rownames(uss.t)

bio.data <- read.csv("bimonitoring data.csv")

css.bio.m <- merge(bio.data[1:4,], css.t, by="CAS")
uss.bio.m <- merge(bio.data[5:6,], uss.t, by="CAS")

css.bio_b <- data.frame(CAS=css.bio.m$CAS, Chemical=css.bio.m$Chemical,
                       value=css.bio.m$biomonitoring.value, name="Biomonitoring data")
css.50 <- data.frame(CAS=css.bio.m$CAS, Chemical=css.bio.m$Chemical,
                     value=css.bio.m$css.50th, name="Css median")
css.95 <- data.frame(CAS=css.bio.m$CAS, Chemical=css.bio.m$Chemical,
                     value=css.bio.m$css.95th, name="Css 95th")
uss.bio_u <- data.frame(CAS=uss.bio.m$CAS, Chemical=uss.bio.m$Chemical,
                        value=uss.bio.m$biomonitoring.value, name="Biomonitoring data")
uss.50 <- data.frame(CAS=uss.bio.m$CAS, Chemical=uss.bio.m$Chemical,
                     value=uss.bio.m$uss.50th, name="Uss median")
uss.95 <- data.frame(CAS=uss.bio.m$CAS, Chemical=uss.bio.m$Chemical,
                     value=uss.bio.m$uss.95th, name="Uss 95th")

css.bind <- rbind(css.bio_b, css.50, css.95)
uss.bind <- rbind(uss.bio_u, uss.50, uss.95)

library(ggplot2)
library(ggpubr)

css.bio_b.plot <- ggplot()+
  geom_point(data=css.bind[5:12,], aes(x=value, y=Chemical, shape=name), size=2)+
  geom_segment(aes(x=css.50$value, y=css.50$Chemical, xend=css.95$value, yend=css.95$Chemical), size=0.2)+
  geom_point(data=css.bind[1:4,], aes(x=value, y=Chemical, shape=name), size=2)+
  scale_x_log10(limits=c(1e-11, 1e-1),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  xlab("Css (mg/L)")+
  theme_bw()+
  theme(axis.title.y = element_blank(), legend.title = element_blank(), 
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  annotation_logticks(sides="b")+
  scale_shape_manual(name="legend",
                     breaks=c("Css median","Css 95th","Biomonitoring data"),
                     labels=c(bquote(C[SS]~50^th), bquote(C[SS]~95^th), "Blood data"),
                     values=c(`Css median`=16,`Css 95th`=1,`Biomonitoring data`=0))
print(css.bio_b.plot)

uss.bio_u.plot <- ggplot()+
  geom_point(data=uss.bind[3:6,], aes(x=value, y=Chemical, shape=name), size=2)+
  geom_segment(aes(x=uss.50$value, y=uss.50$Chemical, xend=uss.95$value, yend=uss.95$Chemical), size=0.2)+
  geom_point(data=uss.bind[1:2,], aes(x=value, y=Chemical, shape=name), size=2)+
  scale_x_log10(limits=c(1e-12, 1e-4),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  xlab("Uss (mg/kg/day)")+
  theme_bw()+
  theme(axis.title.y = element_blank(), legend.title = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  annotation_logticks(sides="b")+
  scale_shape_manual(name="legend",
                     breaks=c("Uss median","Uss 95th","Biomonitoring data"),
                     labels=c(bquote(U[SS]~50^th), bquote(U[SS]~95^th), "Urinary data"),
                     values=c(`Uss median`=16,`Uss 95th`=1,`Biomonitoring data`=2))
print(uss.bio_u.plot)

expocast.plot <- ggarrange(css.bio_b.plot, uss.bio_u.plot, align="hv", nrow=2, labels = c("A","B"))
ggsave(expocast.plot, file="Expocast Css and Uss vs biomonitoring data plot.pdf", 
       width=8, height=6, path="BEMI plots", scale=0.7)
