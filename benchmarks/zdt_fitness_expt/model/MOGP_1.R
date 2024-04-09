library(laGP)

#get training data points
training_data <- read.csv("./model/input/trainingdata.csv")

#set training data for objectives
X <- cbind(X1 = training_data$x1, X2 = training_data$x2, X3 = training_data$x3, X4 = training_data$x4, X5 = training_data$x5,
           X6 = training_data$x6, X7 = training_data$x7, X8 = training_data$x8, X9 = training_data$x9, X10 = training_data$x10,
           X11 = training_data$x11, X12 = training_data$x12, X13 = training_data$x13, X14 = training_data$x14, X15 = training_data$x15,
           X16 = training_data$x16, X17 = training_data$x17, X18 = training_data$x18, X19 = training_data$x19, X20 = training_data$x20,
           X21 = training_data$x21, X22 = training_data$x22, X23 = training_data$x23, X24 = training_data$x24, X25 = training_data$x25,
           X26 = training_data$x26, X27 = training_data$x27, X28 = training_data$x28, X29 = training_data$x29, X30 = training_data$x30)

Y1 <- training_data$obj_1
Y2 <- training_data$obj_2

#get untried location
dv <- read.csv("./model/input/dv.dat", header=T)
Xref_obj <- matrix(unlist(dv),nrow=1)

#laGP appoximation
OBJ1.alc <- laGP(Xref_obj, 10, 60, X, Y1, d=NULL, method = "alc")
OBJ2.alc <- laGP(Xref_obj, 10, 60, X, Y2, d=NULL, method = "alc")

OBJ1_true <- Xref_obj[1]
OBJ1_SD <- 0

#save results to pass to MOU
objective <- c(OBJ1_true, OBJ2.alc$mean, OBJ1_SD, sqrt(OBJ2.alc$s2), OBJ1_true, OBJ1_true)

write.csv(objective, "./model/output/gp_output.dat", quote = F, row.names = F)