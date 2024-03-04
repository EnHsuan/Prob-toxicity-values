library(dplyr)

#19 chemicals from Ford et al. (2022)
hed.df <- read.csv("hed samples CAS by row.csv")
css.df <- read.csv("Css samples CAS by row.csv") #unit: mg/L
td_h.df <- read.csv("td gsd distribution CAS by row.csv")
gfr_fub.df <- read.csv("GFRxFub samples CAS by row.csv") #unit: L/kg/day

chem <- read.csv("regulatory reference dose.csv")
hed.spe.df <- merge(hed.df, chem[,c(2,9)], by="CAS")

set.seed(1234)
### Inter-species TK/TD
#P50=1, P95/P50=3
AF_interTKTD <- rlnorm(10000,meanlog = log(1), sdlog=log(3^(1/1.645)))

hed.a.df <- hed.spe.df[-which(hed.spe.df$Species == "Human"),]
hd50.a <- hed.a.df[,2:10001]/AF_interTKTD
hd50.a <- cbind(hd50.a, hed.a.df[,c(1,10002)])
hd50 <- rbind(hd50.a, hed.spe.df[which(hed.spe.df$Species == "Human"),1:10002])
hd50.df <- hd50 %>% arrange(CAS)

#35 endpoints
#HD50.data.all <- c()
BEMI_b.data.all <- c()
BEMI_u.data.all <- c()
intra.data.all <- c()
HDMI.data.all <- c()
zrand <- rnorm(n=10000) #variability
for (i in 1:19){
  
  cas <- css.df$CAS[i]
  
  css <- as.numeric(css.df[i, 1:10000]) #uncertainty
  css.50 <- quantile(css, probs = 0.5, na.rm=TRUE) #Css median
  
  #parameters for Uss
  gfr_fub <- as.numeric(gfr_fub.df[i, 1:10000]) #uncertainty
  
  #intraspecies TD
  td_h <- as.numeric(td_h.df[which(td_h.df$CAS == cas), 1:10000]) #uncertainty
  
  #HED=BMD/AF_inter_BS
  #hed.data <- as.data.frame(hed.df[which(hed.df$CAS == cas), 1:10000]) #uncertainty
  
  #HD50=BMD/(AF_inter_BS*AF_inter_TD)
  hd50.d <- as.data.frame(hd50.df[which(hd50.df$CAS == cas), 1:10000]) #uncertainty
  
  end <- as.data.frame(hed.df[which(hed.df$CAS == cas), 10002]) #endpoints
  colnames(end) <- "endpoint"
  num.end <- nrow(hd50.d) #endpoint
  
  if(num.end == 1){
    #HD50.data <- data.frame(CAS = cas, endpoint = end)
    BEMI.blood.data <- data.frame(CAS = cas, endpoint = end)
    BEMI.urine.data <- data.frame(CAS = cas, endpoint = end)
    intra.data <- data.frame(CAS = cas, endpoint = end)
    HDMI.data <- data.frame(CAS = cas, endpoint = end)
    
    #HD50=HED/AF_inter_TD
    #HD50.data[1, 3:10002] <- as.numeric(hed.data[1, 1:10000])/AF_interTD
    
    #human median concentration
    #HCss50=(HED*Css.50)/AF_inter_TD  
    #HCss.50 <- (as.numeric(hd50.data[1, 1:10000])*css.50)/AF_interTD
    HCss.50 <- as.numeric(hd50.d[1, 1:10000])*css.50
    
    #BEMI in blood (unit: mg/L)
    BEMI.blood.data[1, 3:10002] <- HCss.50/(td_h^qnorm(0.99)) #I=1%
    
    #BEMI in urine (unit: mg/kg/day)
    for (k in 1:10000){
      HCss.50.temp <- HCss.50[k] #take one sample from HCss.50
      td.rand <- td_h[k]^zrand #take one sample from td_h and run z score to generate variability
      bemi_u.rand <- (HCss.50.temp*gfr_fub)/td.rand #variability combining human TK and TD
      bemi_u.temp <- quantile(bemi_u.rand, prob=0.01, na.rm=TRUE) #I=1%
      BEMI.urine.data[1, 2+k] <- bemi_u.temp
    }
    
    #HDMI (unit: mg/kg/day)
    for (l in 1:10000){
      HCss.50.temp <- HCss.50[l] #take one sample from HCss.50
      td.rand <- td_h[l]^zrand #take one sample from td_h and run z score to generate variability
      hdmi.rand <- HCss.50.temp/(css*td.rand) #variability combining human TK and TD
      intra.temp <- quantile(css*td.rand, prob=0.99, na.rm=TRUE)/quantile(css*td.rand, prob=0.5, na.rm=TRUE) #TKTDVF01
      intra.data[1, 2+l] <- intra.temp
      hdmi.temp <- quantile(hdmi.rand, prob=0.01, na.rm=TRUE) #I=1%
      HDMI.data[1, 2+l] <- hdmi.temp
    } 
    
    
    }else{
      
      #create dataframe/matrix for loop
      #HD50.data <- data.frame(CAS = cas, endpoint = end)
      BEMI.blood.data <- data.frame(CAS = cas, endpoint = end)
      BEMI.urine.data <- data.frame(CAS = cas, endpoint = end)
      intra.data <- data.frame(CAS = cas, endpoint = end)
      HDMI.data <- data.frame(CAS = cas, endpoint = end)
      HCss.50 <- matrix(nrow = num.end, ncol = 10000)
      
      for(j in 1:num.end){
        
        #HD50=HED/AF_inter_TD
        #HD50.data[j, 3:10002] <- as.numeric(hed.data[j, 1:10000])/AF_interTD
        
        #human median concentration
        #HCss50=(HED*Css.50)/AF_inter_TD
        #HCss.50[j,] <- (as.numeric(hed.data[j, 1:10000])*css.50)/AF_interTD
        HCss.50[j,] <- as.numeric(hd50.d[j, 1:10000])*css.50
        
        #BEMI in blood (unit: mg/L)
        BEMI.blood.data[j, 3:10002] <- HCss.50[j,1:10000]/(td_h^qnorm(0.99)) #I=1% 
        
        #BEMI in urine (unit: mg/kg/day)
        HCss.50.temp <- matrix(nrow = num.end, ncol = 1)
        bemi_u.rand <- matrix(nrow = num.end, ncol = 10000)
        bemi_u.temp <- matrix(nrow = num.end, ncol = 1)
        for (k in 1:10000){
          HCss.50.temp  <- HCss.50[j,k] #take one sample from HCss.50
          td.rand <- td_h[k]^zrand #take one sample from td_h and run z score to generate variability
          bemi_u.rand[j,]  <- (HCss.50.temp*gfr_fub)/td.rand #variability combining human TK and TD
          bemi_u.temp[j,]  <- quantile(bemi_u.rand[j,] , prob=0.01, na.rm=TRUE) #I=1%
          BEMI.urine.data[j, 2+k] <- bemi_u.temp[j,] 
        }
        
        #HDMI (unit: mg/kg/day)
        HCss.50.temp <- matrix(nrow = num.end, ncol = 1)
        hdmi.rand <- matrix(nrow = num.end, ncol = 10000)
        intra.temp <- matrix(nrow = num.end, ncol = 1)
        hdmi.temp <- matrix(nrow = num.end, ncol = 1)
        for (l in 1:10000){
          HCss.50.temp <- HCss.50[j,l] #take one sample from HCss.50
          td.rand <- td_h[l]^zrand #take one sample from td_h and run z score to generate variability
          hdmi.rand[j,] <- HCss.50.temp/(css*td.rand) #variability combining human TK and TD
          
          intra.temp[j,] <- quantile(css*td.rand, prob=0.99, na.rm=TRUE)/quantile(css*td.rand, prob=0.5, na.rm=TRUE)
          intra.data[j, 2+l] <- intra.temp[j,]
          
          hdmi.temp[j,] <- quantile(hdmi.rand[j,], prob=0.01, na.rm=TRUE) #I=1%
          HDMI.data[j, 2+l] <- hdmi.temp[j,]
        }
      
    }
    
  }
  
  #HD50.data.all <- rbind(HD50.data.all, HD50.data)
  BEMI_b.data.all <- rbind(BEMI_b.data.all, BEMI.blood.data)
  BEMI_u.data.all <- rbind(BEMI_u.data.all, BEMI.urine.data)
  HDMI.data.all <- rbind(HDMI.data.all, HDMI.data)
  
  intra.data.all <- rbind(intra.data.all, intra.data)
  
}

write.csv(BEMI_b.data.all, file="BEMI blood.csv", row.names = FALSE)
write.csv(BEMI_u.data.all, file="BEMI urine.csv", row.names = FALSE)
write.csv(HDMI.data.all, file="HDMI.csv", row.names = FALSE)
write.csv(intra.data.all, file="TKTDVF01.csv", row.names = FALSE)

HD50.data.all <- data.frame(hd50.df[,10001:10002],hd50.df[,1:10000])
write.csv(HD50.data.all, file="HD50.csv", row.names = FALSE)

HD50.quan <- as.data.frame(t(apply(HD50.data.all[3:10002], 1, quantile, probs=c(0.5, 0.05, 0.95))))
HD50.quan.df <- cbind(HD50.data.all[,1:2], HD50.quan)
write.csv(HD50.quan.df, file="HD50 quantile.csv", row.names = FALSE)

BEMI_b.quan <- as.data.frame(t(apply(BEMI_b.data.all[3:10002], 1, quantile, probs=c(0.5, 0.05, 0.95))))
BEMI_b.quan.df <- cbind(BEMI_b.data.all[,1:2], BEMI_b.quan)
write.csv(BEMI_b.quan.df, file="BEMI blood quantile.csv", row.names = FALSE)

BEMI_u.quan <- as.data.frame(t(apply(BEMI_u.data.all[3:10002], 1, quantile, probs=c(0.5, 0.05, 0.95))))
BEMI_u.quan.df <- cbind(BEMI_u.data.all[,1:2], BEMI_u.quan)
write.csv(BEMI_u.quan.df, file="BEMI urine quantile.csv", row.names = FALSE)

HDMI.quan <- as.data.frame(t(apply(HDMI.data.all[3:10002], 1, quantile, probs=c(0.5, 0.05, 0.95))))
HDMI.quan.df <- cbind(HDMI.data.all[,1:2], HDMI.quan)
write.csv(HDMI.quan.df, file="HDMI quantile.csv", row.names = FALSE)

intra.quan <- as.data.frame(t(apply(intra.data.all[3:10002], 1, quantile, probs=c(0.5, 0.05, 0.95))))
intra.quan.df <- cbind(intra.data.all[,1:2], intra.quan)
write.csv(intra.quan.df, file="TKTDVF01 quantile.csv", row.names = FALSE)