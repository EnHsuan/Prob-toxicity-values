library(httk)

chem <- read.csv("19 overlapping chemical list.csv")

set.seed(1234)

DT_pop <- httkpop_generate(method = "d", #direct resampling
                           nsamp = 10000, #number of individuals
                           gfr_resid_var = TRUE,
                           ckd_epi_race_coeff = FALSE)
DT_mc <- list()
css <- data.frame(matrix(NA, nrow=10000, ncol=19))
colnames(css) <- chem$CAS.No.
par.urine <- list()
fub <- data.frame(matrix(NA, nrow=10000, ncol=19))
gfr <- data.frame(matrix(NA, nrow=10000, ncol=19))
bsa <- data.frame(matrix(NA, nrow=10000, ncol=19))
bw <- data.frame(matrix(NA, nrow=10000, ncol=19))
for (i in 1:19){
  cas <- chem[i,58]
  DT_mc[[i]] <- create_mc_samples(chem.cas=cas, species="Human",
                                  model = "3compartmentss", #insert your chosen model here
                                  httkpop.dt = DT_pop,
                                  httkpop = TRUE, #in httk v2.2.2, this will no longer overwrite your user-supplied 
                                                  #httkpop.dt value – and now it's needed to tell create_mc_samples() 
                                                  #that you are using population physiology
                                  samples = nrow(DT_pop) #need to tell create_mc_samples() how many samples you are 
                                                         #giving it, so it can generate the correct number of samples 
                                                         #for other parameters
  )
  css[, i] <- DT_mc[[i]][,
                         calc_analytic_css(chem.cas=cas, output.units='mg/L', species="Human",
                                           parameters = .SD , #.SD = the subset of columns of DT_mc specified in .SDcols argument
                                           model = "3compartmentss",
                                           suppress.messages = TRUE),
                         by = 1:nrow(DT_mc[[i]]), #go row by row
                         .SDcols = names(DT_mc[[i]]) #pass all columns when .SD is used, not just a subset of them
  ]$V1
  
  
  #extract GFR and fraction unbound, also go row by row
  colnames(DT_mc[[i]])[1] <- "weight_adj"
  par.urine[[i]] <- merge(DT_mc[[i]], DT_pop, by="weight_adj")
  fub[, i] <- par.urine[[i]]$Funbound.plasma #Fraction unbound (unitless)
  gfr[, i] <- par.urine[[i]]$gfr_est #GFR (original unit: mL/min/1.73 m^2 body surface area)
  bsa[, i] <- par.urine[[i]]$BSA_adj #BSA adjusted (unit: cm^2)
  bw[, i] <- par.urine[[i]]$weight_adj #Body weight adjusted (unit: kg)
}

#Css
colnames(css) <- chem$CAS.y
css.t <- as.data.frame(t(css))
css.t$CAS <- row.names(css.t)
write.csv(css.t, file="Css samples CAS by row.csv", row.names = FALSE)

#Uss_tk (unit: L/kg/day)
uss_tk <- gfr*(1/1000)*1440*(bsa/10000)*fub*(1/bw)
colnames(uss_tk) <- chem$CAS.y
uss_tk.t <- as.data.frame(t(uss_tk))
uss_tk.t$CAS <- row.names(uss_tk.t)
write.csv(uss_tk.t, file="GFRxFub samples CAS by row.csv", row.names = FALSE)
