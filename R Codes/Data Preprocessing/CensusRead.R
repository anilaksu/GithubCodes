####################################################
#	                                             #
#	       Program main starts here              #
#     						         #
####################################################	

## functions
source('FuzzyFunc.R')
source('ConfidenceTable.R')
source('AssociationRules.R')

### reading data from .dat file 
Census <- read.table("Census.dat", header = FALSE, sep="\t",row.names=NULL) 

### let's see the dimensions of data
 print("The dimensions of data")
 print(dim(Census))



### Cities
City<-as.character(Census[,1])

print("Cities")
print(City)

### Districts
Mahalle<-as.character(Census[,2])

print("Districts")
print(Mahalle)

### Development index 
DevInd<-as.character(Census[,3])

print("Development Index")
print(DevInd)


## age array for all ages
Age<-matrix(nrow=40, ncol=10)

## all ages for given districts
Age<-Census[,4:13]

print("Age Distribution")
print(Age)

#Age[,1]<-census$yas_0_4 

#Age[,2]<-census$yas_5_9

#Age[,3]<-census$yas_10_14
 
#Age[,4]<-census$yas_15_19 

#Age[,5]<-census$yas_20_24

#Age[,6]<-census$yas_25_29

#Age[,7]<-census$yas_30_34

#Age[,8]<-census$yas_35_39

#Age[,9]<-census$yas_40_44

#Age[,10]<-census$yas_45_49 

## all ages for given districts
Edu<-Census[,14:16]

######################################
#                                    #
# By using data given above we wanna # 
# associate development index and age# 
# with foursquare using fraction     #
#                                    #
######################################

## 

#plot(Age[1,])

SampleAge<-numeric(10)

for (i in 1:10){
	
	SampleAge[i]<-Age[1,i]
	
}

  
AgeMem<-MemberFunAge(SampleAge)

print("MemberShip Values for Age Distribution")
print(AgeMem)

DevMem<-MemberFunDev(DevInd[19])

print("MemberShip Values for Development Index")
print(DevMem)

#print(Edu)
educationTable <- GetEducation(Edu)
developmentTable <- GetDevelopment(DevInd)
ageTable <- GetAge(Age)

print(educationTable)
print(developmentTable)
print(ageTable)

ageClass <- CheckSupport(ageTable,0.2)
developmentClass <- CheckSupport(developmentTable,0.2)
educationClass <- CheckSupport(educationTable,0.2)

print(ageClass)
print(educationClass)
print(developmentClass)

## let's generate Correction factors for each district
numDistrict <- length(ageTable[,1])
foursquareRatio <- numeric(numDistrict)
foursquareRatio <- runif(numDistrict, 1, 32)

print(foursquareRatio)

correctionTable <- GetCorrection(foursquareRatio)
correctionClass <- CheckSupport(correctionTable,0.2)

print(correctionClass)

## let's collect all class and tables in list

educationList <- list(table= educationTable, class= educationClass)
ageList <- list(table= ageTable, class= ageClass)
developmentList <- list(table= developmentTable, class= developmentClass)
correctionList <- list(table= correctionTable, class= correctionClass)

sampleList <- GetRules(educationList,ageList,developmentList,correctionList) 



