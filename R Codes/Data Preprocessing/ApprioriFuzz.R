
## library to read census data
library("maptools")
library("pracma")
library("rgeos")
library("sp")
library("rgdal")
library(rmongodb)


## functions
source('FuzzyFunc.R')

## let's read the data
census<-readShapeSpatial("ornek_ses.shp")

## objectid's
ID<-census$OBJECTID

## development index
DevInd<- census$tr_index

## city name
ILADI<-census$ILADI

## district name
Mahalle<-census$Mahalle

Num<-length(ILADI)

## age array for all ages

Age<-matrix(nrow=Num, ncol=10)

Age[,1]<-census$yas_0_4 

Age[,2]<-census$yas_5_9

Age[,3]<-census$yas_10_14
 
Age[,4]<-census$yas_15_19 

Age[,5]<-census$yas_20_24

Age[,6]<-census$yas_25_29

Age[,7]<-census$yas_30_34

Age[,8]<-census$yas_35_39

Age[,9]<-census$yas_40_44

Age[,10]<-census$yas_45_49     

## education data

Education<-matrix(nrow=Num, ncol=3)

Education[,1]<-census$egitim_io 

Education[,2]<-census$egitim_lis

Education[,3]<-census$egitim_uni

## this transforms into Turkish character
#Encoding(DevInd) <- "UTF-8"
print(typeof(DevInd))

tic()
print(DevInd)
toc()

print(ILADI)

print(ID)

print(Mahalle)

print(Age[,1])

print(Education)

Tab<- matrix(nrow=Num, ncol=16)

Tab[,1]<-as.character(ILADI)

Tab[,2]<-as.character(Mahalle)

Tab[,3]<-as.character(DevInd)

Tab[,4:13]<-Age

Tab[,14:16]<-Education

 write.table(Tab, "Census.dat", sep="\t", row.names=FALSE, col.names=FALSE, eol="\n")



mongo <- mongo.create()

cursor = mongo.find(mongo, ns="yerbil.categories");

pipe_1 <- mongo.bson.from.JSON('{"$group":{"_id":"$VenueID","TotUser":{"$max":"$usersCount"},"lat":{"$max":"$lat"},"lon":{"$max":"$lng"}}}');
 
cmd_list<-list(pipe_1)

res <- mongo.aggregation(mongo,"yerbil.venues", cmd_list)

Venues <- mongo.bson.to.list(res)

##  total number of venues used in Correction Analysis
TotNumVenues<-length(Venues$result)

UniCheckCoord<- matrix(nrow=TotNumVenues, ncol=4)

VenueIDs<-character(TotNumVenues)

for(i in 1:TotNumVenues){
# VenueID
VenueIDs[i]<-Venues$result[[i]]$'_id'
# lattitudes
UniCheckCoord[i,4]<-Venues$result[[i]]$TotUser
# lattitudes
UniCheckCoord[i,2]<-Venues$result[[i]]$lon
# longitudes
UniCheckCoord[i,3]<-Venues$result[[i]]$lat
}


 print("Aggregate Coordinates") 
 print(Venues$result[[1]]$lat)

## here we check location of each coordinate
Index<-integer(TotNumVenues) 
for(i in 1:Num){

 poly<-census@polygons[[i]]@Polygons[[1]]@coords

 res<-point.in.polygon(UniCheckCoord[,2], UniCheckCoord[,3], poly[,1], poly[,2], mode.checked=FALSE)

 # to give the index let's multiply it with i
 res<-res*i;

 # let's update index array
 Index<-Index+res;

}

TotCheckDis<-numeric(Num)

# let's calculate number of people in each location
for(i in 1:TotNumVenues){

 if(Index[i]>0) {
   TotCheckDis[Index[i]]<-TotCheckDis[Index[i]]+UniCheckCoord[i,4]
  }

}


## population of each district

Population<-numeric(Num)
 
 for(i in 1:Num){
 Population[i] <-sum(Age[i,])
 }
print(TotCheckDis)

 print(Population)


### let's calculate membership values

#EduMemVal<-MemberFunEdu(Education[1,])
 
#print(EduMemVal)

## let's find membership values for age
 AgeMem <-AgeMember(Age)

## let's find membership values for development index
 DevMem<-DevMember(DevInd)

## let's find membership values for education level
 EduMem<-EduMember(Education)

## let's find membership values for Correction Factor
 CorrMem<-CorrMember(TotCheckDis,Population)

## Here 
print(CorrMem)

  

