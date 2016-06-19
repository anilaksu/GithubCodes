
GetRules <- function(educationList,ageList,developmentList,correctionList) {

####################################################
#	                                             #
#	This function is developed to generate       #
#     education table for each district            #
#	                                             #
####################################################	 
	
####################################################
#	                                             #
#educationList <- list(table= educationTable, class= educationClass)
#ageList <- list(table= ageTable, class= ageClass)
#developmentList <- list(table= developmentTable, class= developmentClass)
#correctionList <- list(table= correctionTable, class= correctionClass)
#
####################################################		
	
####################################################
#	                                             #
#	first it generates each districts membership #
#     values and list them in a table              #
#	                                             #
####################################################

educationTable <- educationList$table
educationClass <- educationList$class

developmentTable <- developmentList$table
developmentClass <- developmentList$class

correctionTable <- correctionList$table
correctionClass <- correctionList$class

ageTable <- ageList$table
ageClass <- ageList$class

nDistrict <- length(ageTable[,1])

## let's total number of combinations
nTotal <- length(educationClass)*length(ageClass)*length(developmentClass)*length(correctionClass)
print(nTotal)
ruleList <- array(list( index = list(education = integer(1), age = integer(1),development = integer(1), correction = integer(1), min = numeric(1))) 
, nTotal )

for (i in 1:nDistrict ){
## here we get values for all possible rules
NewList <- GetRuleValues(educationTable[i,],ageTable[i,],developmentTable[i,],correctionTable[i,])
ruleList <- AddValues(ruleList,NewList) 
#print(ruleList)

}

return(ruleList)

}

GetRuleValues <- function(education,age,development,correction) {

####################################################
#	                                             #
#	This function is developed to find the       #
#     minimum value
#	                                             #
####################################################	 

## let's total number of combinations
nTotal <- length(education)*length(age)*length(development)*length(correction)

## let's store indexes 
minValues <- array(list( index = list(education = integer(1), age = integer(1),development = integer(1), correction = integer(1), min = numeric(1))) 
, nTotal )

nCount <- 1
## first for loop is generated for education combinations
for(i in 1:2){
 ## second for loop is generated for age combinations
 for(j in 1:3) {
  ## third loop is generated for development combinations
   for (k in 1:3) {
    ## forth loop is generated for correction factor combinations
    for (l in 1:6) {
    minValues[[nCount]]$education <- i
    minValues[[nCount]]$age <- j
    minValues[[nCount]]$development <- k
    minValues[[nCount]]$correction <- l
    minValues[[nCount]]$min <- GetMinimum(education[i],age[j],development[k],correction[l])
    nCount <- nCount + 1 
    }
   }
 } 
}
return(minValues)
}

GetMinimum <- function(education,age,development,correction) {

####################################################
#	                                             #
#	This function is developed to find the       #
#     minimum value
#	                                             #
####################################################	 
	
## first let's store them in array

parameterArray <- numeric(4)
parameterArray[1] <- education
parameterArray[2] <- age
parameterArray[3] <- development
parameterArray[4] <- correction

minimum <- min(parameterArray)

return(minimum)
}

AddValues <- function(OldList,NewList) {

####################################################
#	                                             #
#	This function is developed to update the list#
#     list we store values
#	                                             #
####################################################	 
	
nCount <- 1
## first for loop is generated for education combinations
for(i in 1:2){
 ## second for loop is generated for age combinations
 for(j in 1:3) {
  ## third loop is generated for development combinations
   for (k in 1:3) {
    ## forth loop is generated for correction factor combinations
    for (l in 1:6) {
    OldList[[nCount]]$education <- i
    OldList[[nCount]]$age <- j
    OldList[[nCount]]$development <- k
    OldList[[nCount]]$correction <- l
    OldList[[nCount]]$min <- OldList[[nCount]]$min +NewList[[nCount]]$min
    nCount <- nCount + 1 
    }
   }
 } 
}
#print(nCount)
#print(OldList)
return(OldList)
}






