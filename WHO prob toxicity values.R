library(dplyr)

who_hed.df <- read.csv("WHO hed samples CAS by row.csv")

chem <- read.csv("regulatory reference dose.csv")
who_hed.spe.df <- merge(who_hed.df, chem[,c(2,9)], by="CAS")

set.seed(1234)
### Inter-species TK/TD
#P50=1, P95/P50=3
AF_interTKTD <- rlnorm(10000,meanlog = log(1), sdlog=log(3^(1/1.645)))

who_hed.a.df <- who_hed.spe.df[-which(who_hed.spe.df$Species == "Human"),]
who_hd50.a <- who_hed.a.df[,2:10001]/AF_interTKTD
who_hd50.a <- cbind(who_hd50.a, who_hed.a.df[,c(1,10002)])
who_hd50 <- rbind(who_hd50.a, who_hed.spe.df[which(who_hed.spe.df$Species == "Human"),1:10002])
who_hd50.df <- who_hd50 %>% arrange(endpoint)

### Intra-species
# Chiu et al. (2018) Figure 4
# Incidence = 1%
#P50 = 9.7
#P95/P50 = 4.3
AF_intra <- rlnorm(10000, meanlog = log(9.7), sdlog = log(4.3^(1/1.645)))

#22 endpoints
#HD50.data.all <- c()
HDMI.data.all <- c()
for (i in 1:22){ 
  
  cas <- who_hd50.df$CAS[i]
  
  #BMD/AF_inter_BS
  #hed.data <- who_hed.df[i, 1:10000]
  hd50.d <- who_hd50.df[i, 1:10000]
  
  end <- as.data.frame(who_hed.df[i, 10002]) #endpoints
  colnames(end) <- "endpoint"
  num.end <- nrow(hd50.d) #endpoint
  
  if(num.end == 1){
    #HD50.data <- data.frame(CAS = cas, endpoint = end)
    #HD50.data[1, 3:10002] <- as.numeric(hed.data[1, 1:10000])/AF_interTD
    HDMI.data <- data.frame(CAS = cas, endpoint = end)
    HDMI.data[1, 3:10002] <- as.numeric(hd50.d[1, 1:10000])/AF_intra
    
  }else{
    
    #HD50.data <- data.frame(CAS = cas, endpoint = end)
    HDMI.data <- data.frame(CAS = cas, endpoint = end)
    
    for(j in 1:num.end){
      #HD50.data[j, 3:10002] <- as.numeric(hed.data[j, 1:10000])/AF_interTD
      HDMI.data[j, 3:10002] <- as.numeric(hd50.d[j, 1:10000])/AF_intra
    }
    
  }
  
  #HD50.data.all <- rbind(HD50.data.all, HD50.data)
  HDMI.data.all <- rbind(HDMI.data.all, HDMI.data)
  
}

write.csv(HDMI.data.all, file="WHO approximate HDMI.csv", row.names = FALSE)

HD50.data.all <- data.frame(who_hd50.df[,10001:10002],who_hd50.df[,1:10000])
write.csv(HD50.data.all, file="WHO approximate HD50.csv", row.names = FALSE)

HDMI.quan <- as.data.frame(t(apply(HDMI.data.all[3:10002], 1, quantile, probs=c(0.5, 0.05, 0.95))))
HDMI.quan.df <- cbind(HDMI.data.all[,1:2], HDMI.quan)
write.csv(HDMI.quan.df, file="WHO approximate HDMI quantile.csv", row.names = FALSE)

AF_intra.df <- data.frame(t(AF_intra))
write.csv(AF_intra.df, file="WHO approximate intraspecies factor.csv", row.names = FALSE)
AF_intra.quan.df <- as.data.frame(t(apply(AF_intra.df, 1, quantile, probs=c(0.5, 0.05, 0.95))))
write.csv(AF_intra.quan.df, file="WHO approximate intraspecies factor quantile.csv", row.names = FALSE)
