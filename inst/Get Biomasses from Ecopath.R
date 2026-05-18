#Redistribute and Calculate Biomasses from Ecopath

#Load in ecopath data
ecopath_estimates <-
    read.csv(here("inst","extdata","Ecopath-Basic estimates.csv"))

#Use Weightings from Lauria 2012 to split out Red Gurnard and Boarfish from
#grouped columns
Gurnards <- ecopath_estimates[46, ]
Boarfish <- ecopath_estimates[49, ]
num_cols <- sapply(Gurnards, is.numeric)
Red_gurnard<-Gurnards
Other_Gurnards<-Gurnards
Boarfish<-Boarfish
Small_pelag<-Boarfish
#37% of the biomass is from red gurnard
Red_gurnard[num_cols] <- Gurnards[num_cols] * 0.37
Red_gurnard$Group.name<-"Red gurnard"
#63% is for the other gurnards
Other_Gurnards[num_cols] <- Gurnards[num_cols] * 0.63
#97% of the biomass is from boarfish
Boarfish[num_cols]<- Boarfish[num_cols] * 0.97
Boarfish$Group.name<-"Boarfish"
#3% is from the other small pelagic species
Small_pelag[num_cols]<- Boarfish[num_cols] * 0.03

#Doesn't matter that tropic level and hab area are being affected as they are
#not used in other calculations
ecopath_params<-ecopath_estimates[-c(46,49),]
ecopath_params<-rbind(ecopath_params,Other_Gurnards, Red_gurnard, Boarfish,
                      Small_pelag)
