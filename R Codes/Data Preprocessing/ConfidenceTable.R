# this source code is written to generate 
# the table where membership values for each district is collected

GetCorrection <- function(foursquareRatio) {

####################################################
#	                                             #
#	This function is developed to generate       #
#     education table for each district            #
#	                                             #
####################################################	 
	
####################################################
#	                                             #
#	    Education[,1]<-census$egitim_io          #
#Input:   Education[,2]<-census$egitim_lis         #
#	    Education[,3]<-census$egitim_uni         #
####################################################		
	
####################################################
#	                                             #
#	first it generates each districts membership #
#     values and list them in a table              #
#	                                             #
####################################################

# let's find the total number of districts
nDistrict <- length(foursquareRatio)

correctionTable<-matrix(nrow=nDistrict, ncol=6)
for (i in 1:nDistrict) {
 correctionTable[i,]  <- CorrMember(foursquareRatio[i])
 #educationTable[i,1] <- eduMem[[1]]
 #educationTable[i,2] <- eduMem[[2]]
}
 return(correctionTable)
}

GetEducation <- function(education) {

####################################################
#	                                             #
#	This function is developed to generate       #
#     education table for each district            #
#	                                             #
####################################################	 
	
####################################################
#	                                             #
#	    Education[,1]<-census$egitim_io          #
#Input:   Education[,2]<-census$egitim_lis         #
#	    Education[,3]<-census$egitim_uni         #
####################################################		
	
####################################################
#	                                             #
#	first it generates each districts membership #
#     values and list them in a table              #
#	                                             #
####################################################

# let's find the total number of districts
nDistrict <- length(education[,1])
#print(nDistrict)
educationTable<-matrix(nrow=nDistrict, ncol=2)
for (i in 1:nDistrict) {
# educationTable[i,] 
 eduMem  <- MemberFunEdu(education[i,])
 educationTable[i,1] <- eduMem[[1]]
 educationTable[i,2] <- eduMem[[2]]
}
 return(educationTable)
}

GetDevelopment <- function(developmentIndex) {

####################################################
#	                                             #
#	This function is developed to generate       #
#     education table for each district            #
#	                                             #
####################################################	 
	
####################################################
#	                                             #
#	    Education[,1]<-census$egitim_io          #
#Input:   Education[,2]<-census$egitim_lis         #
#	    Education[,3]<-census$egitim_uni         #
####################################################		
	
####################################################
#	                                             #
#	first it generates each districts membership #
#     values and list them in a table              #
#	                                             #
####################################################

# let's find the total number of districts
nDistrict <- length(developmentIndex)
#print(nDistrict)
developmentTable<-matrix(nrow=nDistrict, ncol=3)
for (i in 1:nDistrict) {
# educationTable[i,] 
 developmentTable[i,]<- MemberFunDev(developmentIndex[i])
 # print(devTab)
 #educationTable[i,1] <- eduMem[[1]]
 #educationTable[i,2] <- eduMem[[2]]
}
 return(developmentTable)
}

GetAge <- function(ageDistribution) {

####################################################
#	                                             #
#	This function is developed to generate       #
#     education table for each district            #
#	                                             #
####################################################	 
	
####################################################
#	                                             #
#	    Education[,1]<-census$egitim_io          #
#Input:   Education[,2]<-census$egitim_lis         #
#	    Education[,3]<-census$egitim_uni         #
####################################################		
	
####################################################
#	                                             #
#	first it generates each districts membership #
#     values and list them in a table              #
#	                                             #
####################################################

# let's find the total number of districts
nDistrict <- length(ageDistribution[,1])
#print(nDistrict)
ageTable<-matrix(nrow=nDistrict, ncol=3)
for (i in 1:nDistrict) {
# educationTable[i,] 
 #ageTable[i,]  
 ageTab <- MemberFunAge(ageDistribution[i,])
 ageTable[i,1] <- ageTab[[1]]
 ageTable[i,2] <- ageTab[[2]]
 ageTable[i,3] <- ageTab[[3]]
 }
 return(ageTable)
}

CheckSupport <- function(Table,Support) {

####################################################
#	                                             #
#	This function is developed to generate       #
#     education table for each district            #
#	                                             #
####################################################	 
	
####################################################
#	                                               #
#           Input: Table[nDistrict,2 or 3]         #
#												   #
####################################################		
	
####################################################
#	                                             #
#	first it generates each districts membership #
#     values and list them in a table              #
#	                                             #
####################################################

# let's find the total number of districts
nDistrict <- length(Table[,1])
nClass <- length(Table[1,])

satisfySupp <- numeric(nClass)
for (i in 1:nClass) {
 if(sum(Table[,i])/nDistrict > Support) {
   	satisfySupp[i] <- 1
 } else
    satisfySupp[i] <- 0
}
 return(satisfySupp)
}