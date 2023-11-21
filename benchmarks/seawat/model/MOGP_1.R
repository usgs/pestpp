library(laGP)

#get training data points
training_data <- read.csv("./model/input/trainingdata.csv")

#set training data for objectives
X <- training_data[,2:11]

Y1 <- training_data$mean_conc
Y2 <- training_data$wetland_dd

#get untried location
dv <- read.csv("./model/input/flow.wel_stress_period_data_scenario_base.txt", header=F, sep = "")[86:93,4]
dv <- append(dv,read.csv("./model/input/artrch.dat", header=F, sep = "")[1:2,2])

Xref_obj <- matrix(unlist(dv),nrow=1)

#laGP appoximation
mn_conc.alc <- laGP(Xref_obj, 10, min(round(nrow(X)/2), 150), X, Y1, d=list(mle = T, start = 10, min = 0.1, max = 1500), method = "alc")
wetland_dd.alc <- laGP(Xref_obj, 10, min(round(nrow(X)/2), 150), X, Y2, d=list(mle = T, start = 10, min = 0.1, max = 1500), method = "alc")

total_pump <- -sum(dv[1:8])
ar_rate_total <- sum(dv[9:10])
ar_rate_total_sd <- 0.00001

#save results to pass to MOU
row_names <- c("mean_conc", 
               "salinity", 
               "total_pump_rate", 
               "ar_rate_total", 
               "wetland_dd", 
               "ar_rate_t",
               "mean_conc_sd", 
               "salinity_sd", 
               "ar_rate_total_sd",
               "wetland_dd_sd")

objective <- c(mn_conc.alc$mean, 
               mn_conc.alc$mean, 
               total_pump,
               ar_rate_total,
               wetland_dd.alc$mean,
               ar_rate_total,
               sqrt(mn_conc.alc$s2),
               sqrt(mn_conc.alc$s2),
               ar_rate_total_sd, 
               sqrt(wetland_dd.alc$s2))

write.csv(objective, "./model/output/gp_output.dat", quote = F, row.names = row_names)
