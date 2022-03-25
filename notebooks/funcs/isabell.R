#import library
library("readxl")

#read dataset
influVE = read_excel("influVE.xlsx")
influVE_num = influVE 


#change sex from 1,2 numeric to catercorial to male, female
influVE$sex = as.factor(influVE$sex) 
levels(influVE$sex) <- list(male = "1", female = "2")

#change vacc from 0,1 numeric to catercorial to vaccinated and notvaccinated
influVE$vacc = as.factor(influVE$vacc) 
levels(influVE$vacc) <- list(notvaccinated = "0", vaccinated = "1")

#change sex from 1,2 numeric to catercorial to male, female
influVE$flu = as.factor(influVE$flu) 
levels(influVE$flu) <- list(PCRnegative = "0", PCRpositive = "1")

#change age from numeric age to categeorial overequal60 to under60
test = influVE$age >= 60 #test array, where is influVE$age column bigger or euqal to 60
test = as.factor(test) #factorize to TRUE and FALSE
levels(test) <- list(Overequal60 = "TRUE", Under60 = "FALSE") #change categories to Over60 Under60
influVE$age = test #overwrite in influVE

plot(influVE$vacc, influVE$sex, xlab = "Vaccination Status", ylab = "Gender")

plot(influVE$vacc, influVE$age, xlab = "Vaccination Status", ylab = "Age")

plot(influVE$vacc, influVE$flu, xlab = "Vaccination Status", ylab = "PCR Test Status")

plot(as.factor(influVE$dsample))

library(epitools)
out1 <- oddsratio(with(influVE, table(influVE$vacc, influVE$flu)))
out1

out2 <- riskratio(with(influVE, table(influVE$vacc, influVE$flu)))
out2$measure

