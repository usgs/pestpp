library(laGP)

#get training data points
training_data <- read.csv("./model/input/trainingdata.csv")

#set training data for objectives
X <- cbind(X1 = training_data$x1, X2 = training_data$x2, X3 = training_data$x3)
          
Y1 <- training_data$obj_1
Y2 <- training_data$obj_2

#get untried location
dv <- read.csv("./model/input/dv.dat", header=T)
Xref_obj <- matrix(unlist(dv),nrow=1)

#laGP appoximation
OBJ1.alc <- laGP(Xref_obj, 10, 60, X, Y1, d=NULL, method = "alc")
OBJ2.alc <- laGP(Xref_obj, 10, 60, X, Y2, d=NULL, method = "alc")

Pf_obj_2 <- 1-pnorm(-20-OBJ2.alc$mean/sqrt(OBJ2.alc$s2))

#save results to pass to MOU
objective <- c(OBJ1.alc$mean, OBJ2.alc$mean, sqrt(OBJ1.alc$s2), sqrt(OBJ2.alc$s2), Pf_obj_2)

write.csv(objective, "./model/output/gp_output.dat", quote = F, row.names = F)