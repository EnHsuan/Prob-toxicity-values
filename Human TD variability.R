#read LCL data
chemlist <- read.csv("19 overlapping chemical list.csv")
tdvf_ford <- read.csv("tdvf_ford et al 2022.csv")
chemlist.df <- data.frame(Chemical=chemlist$Chemical, CAS=chemlist$CAS.y)

chemtdvf.ford <- merge(chemlist.df, tdvf_ford, by="CAS")
chemtdvf.ford.df <- chemtdvf.ford[c(1:18,20),c(1,3:6)]
colnames(chemtdvf.ford.df)[2] <- "Chemical"
chemtdvf <- chemtdvf.ford.df

rownames(chemtdvf) <- chemtdvf[,1]
chemtdvf.df <- chemtdvf[,-c(1:2)] #generate a dataframe of TDVF quantiles

#TDVF05 = exp(z_0.95*sigma_H)
#qnorm(0.95) = 1.645
#ln(TDVF05) = 1.645 * sigma_H
#sigma_H = ln(TDVF05)/1.645
lntdvf.df <- log(chemtdvf.df)
sigma.df <- lntdvf.df/1.645 #assume lognormal

#ln(sigma_H) ~ Normal distribution
lnsigma.df <- log(sigma.df)
lnsigma.df$mean <- lnsigma.df$Unc50.
lnsigma.df$sd <- (lnsigma.df$mean-lnsigma.df$Unc5.)/1.645

lnsigma.df$CAS <- rownames(lnsigma.df)

lnsigma_d <- data.frame(matrix(NA, nrow=10000, ncol=19))
sigma_d <- data.frame(matrix(NA, nrow=10000, ncol=19))
expsigma_d <- data.frame(matrix(NA, nrow=10000, ncol=19))
for (i in 1:19){
  lnsigma_d[,i] <- rnorm(10000, mean = lnsigma.df$mean[i], sd = lnsigma.df$sd[i]) #ln(sigma) distribution
  sigma_d[,i] <- exp(lnsigma_d[[i]]) #exp(lnsigma)=sigma
  expsigma_d[,i] <- exp(sigma_d[,i]) #exp(sigma)=GSD
}

colnames(lnsigma_d) <- lnsigma.df$CAS
colnames(sigma_d) <- lnsigma.df$CAS
colnames(expsigma_d) <- lnsigma.df$CAS
expsigma_d.t <- as.data.frame(t(expsigma_d))
expsigma_d.t$CAS <- row.names(expsigma_d.t)
write.csv(expsigma_d.t, file="td gsd distribution CAS by row.csv", row.names = FALSE)