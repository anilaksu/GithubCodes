# In this subroutine we define our membership functions for age
# distribution, development index, education and Correction factor 
# which is the ratio of the total population in a given district to
# total number of people using Foursquare-Swarm

CorrMember<- function(foursquareRatio){
	
####################################################
#	                                             #
#    This function is developed to find membership #
#    values for Correction Factor                  #
#	                                             #
####################################################	 
	
####################################################
#	                                             #
#	                                             #
#Input:   TotFour: Total People using Foursquare   #
#         in given region                          #
#         TotPop: Total Population in given region #
#	                                             #
####################################################		
	
	
####################################################
#	                                             #
#	Output: f_membership                         #
#   Form of f_membership: Log(1): C_r=1            #
#  				  Log(2): C_r=2            #
#       			  Log(4): C_r=4            #
#                         Log(8): C_r=8   	   #
#                         Log(16): C_r=16		   #
#                         Log(32): C_r=32		   #
#	                                             #
####################################################

#### Number of districts are calculated 
numDistrict <- length(foursquareRatio)

## let's write ratio's 
invRatio <- numeric(6)
for (i in 1:6) {
invRatio[i] <- 1/(2^(i-1))
} 
invFourRatio <- 1/foursquareRatio
#print(invRatio)
#print(invFourRatio)
## let's generate membership values
memberCorr <- numeric(6)
for(i in 1:5){
 if((invRatio[i] >= invFourRatio) & (invRatio[i+1] < invFourRatio) ) {
  subInterval <- invRatio[i] - invRatio[i+1];
  distToIntI <- (invRatio[i] - 1/foursquareRatio)/subInterval
  memberCorr[i] <- 1 - distToIntI;
  memberCorr[i+1] <- distToIntI;
  }
}
 return(memberCorr)	

}



EduMember<- function(Education){
	
####################################################
#	                                             #
#  This function is developed to find membership   #
#    values of Age Distribution Table              #
#	                                             #
####################################################	 
	
####################################################
#	                                             #
#	    Education[,1]<-census$egitim_io		   #
#Input:   Education[,2]<-census$egitim_lis	   #	
#	    Education[,3]<-census$egitim_uni	   #	
####################################################		
	
	
####################################################
#	                                             #
#	Output: f_membership                         #
#   Form of f_membership: Low: [0-1]               #
#  				  Medium: [0-1]            #
#       			  High:  [0-1]             #
#	                                             #
####################################################

#### it organizes all Age distribution and returns a matrix

#### Number of districts are calculated 

Ndist<- length(Education[,1])	

### Now let's generate matrix we store member ship values

EduMem<- matrix(nrow=Ndist, ncol=2)

 for(i in 1:Ndist){

  EduMem[i,]<-MemberFunEdu(Education[i,])

 }
  
   return(EduMem)		
	
}



DevMember<- function(DevInd){
	
#####################################################
#	                                              #
#    This function is developed to find membership  #
#    values of Development Table                    #
#	                                              #
#####################################################	 
	
####################################################
#	                                             #
#Input:   DevInd: in form of Alt-Orta-Ust          #
#                                                  #
####################################################		
	
	
####################################################
#	                                             #
#	Output: f_membership                         #
#   Form of f_membership: Low: [0-1]               #
#  				  Medium: [0-1]            #
#       			  High:  [0-1]             #
#	                                             #
####################################################

#### it organizes all Age distribution and returns a matrix

#### Number of districts are calculated 

Ndist<- length(DevInd)	

### Now let's generate matrix we store member ship values

DevMem<- matrix(nrow=Ndist, ncol=3)

 for(i in 1:Ndist){

  DevMem[i,]<-MemberFunDev(DevInd[i])

 }
  
   return(DevMem)		
	
}

AgeMember<- function(Age){
	
####################################################
#	                                             #
#	This function is developed to find membership#
#    values of Age Distribution Table              #
#	                                             #
####################################################	 
	
####################################################
#	                                             #
#         Age[,1]<-census$yas_0_4 
#	    Age[,2]<-census$yas_5_9
#	    Age[,3]<-census$yas_10_14
#	    Age[,4]<-census$yas_15_19 
#Input:   Age[,5]<-census$yas_20_24
#	    Age[,6]<-census$yas_25_29
#         Age[,7]<-census$yas_30_34
#         Age[,8]<-census$yas_35_39
#         Age[,9]<-census$yas_40_44
#         Age[,10]<-census$yas_45_49  	         #
####################################################		
	
	
####################################################
#	                                             #
#	Output: f_membership                         #
#   Form of f_membership: Low: [0-1]               #
#  				  Medium: [0-1]            #
#       			  High:  [0-1]             #
#	                                             #
####################################################

#### it organizes all Age distribution and returns a matrix

#### Number of districts are calculated 

Ndist<- length(Age[,1])	

### Now let's generate matrix we store member ship values

AgeMem<- matrix(nrow=Ndist, ncol=3)

 for(i in 1:Ndist){

  AgeMem[i,]<-MemberFunAge(Age[i,])

 }
  
   return(AgeMem)		
	
}


MemberFunAge<-function(AgeDist){
	
####################################################
#	                                             #
#	This function is developed to find membership#
#    values of given Age Distribution              #
#	                                             #
####################################################


###################################################
#	                                             #
#	Input: AgeDist same age distribution for one #
#   analysis region                                #
#	                                             #
####################################################

# We have to normalize the given distibution with respect to
# total people in a given district
TotPeop<- sum(AgeDist)	;

# Normalized distribution
NormDist<-AgeDist/TotPeop;

# let's define membership
f<- numeric(3)
# f[1]: Low
# f[2]: Medium
# f[3]: High

# Let's calculate Low value for given distribution

for(i in 1:6){

if(i==6){
	f[1]<-f[1]+	0.5*NormDist[i];
} 	else
   f[1]<-f[1]+NormDist[i];
	
}


# Let's calculate Medium value for given distribution

for(i in 6:9){

if(i==6){
	f[2]<-f[2]+	0.5*NormDist[i];
} 	else if(i==9){
	f[2]<-f[2]+	0.5*NormDist[i];
}   else 
   f[2]<-f[2]+NormDist[i];
	
}

	
# Let's calculate High value for given distribution

for(i in 9:10){

if(i==9){
	f[3]<-f[3]+	0.5*NormDist[i];
} 	else 
	f[3]<-f[3]+	NormDist[i];
	
	}
	
	## it only returns f vector
	
	return(f)
	
}

MemberFunEdu<-function(EduDist){
	
###################################################
#	                                               #
#	This function is developed to find membership  #
#    values of given Education Distribution           #
#	                                               #
######################################################


###################################################
#	                                               #
#	Input: EduDist same age distribution for one   #
#   analysis region                                #
#	                                               #
####################################################

# We have to normalize the given distibution with respect to
# total people in a given district
TotPeop<- sum(EduDist[2:3])	;

# Normalized distribution
NormDist<-EduDist[2:3]/TotPeop;

# let's define membership
f<- numeric(2)

# f[2]: Medium
# f[3]: High

# Let's calculate membership for given distribution

 f[1]<-NormDist[1]
 f[2]<-NormDist[2]
 #f[3]<-NormDist[3]
	
	return(f)
	
}


MemberFunDev<-function(DevInd){
	
###################################################
#	                                               #
#	This function is developed to find membership  #
#    values of given developmet index              #
#	                                               #
####################################################


###################################################
#	                                               #
#	Input: Development Index                       #
#	                                               #
####################################################

# let's define membership
f<- numeric(3)
# f[1]: Alt
# f[2]: Orta
# f[3]: Ãœst

if(DevInd=="Alt"){
 f[1]<-1
 f[2]<-0
 f[3]<-0	
} else if(DevInd=="Alt-Orta"){
 f[1]<-0.5
 f[2]<-0.5
 f[3]<-0	
} else if(DevInd=="Orta"){
 f[1]<-0
 f[2]<-1
 f[3]<-0
} else if(DevInd=="Orta-Ãœst"){
 f[1]<-0
 f[2]<-0.5
 f[3]<-0.5
}else if(DevInd=="Ãœst"){
 f[1]<-0
 f[2]<-0
 f[3]<-1
}

	return(f)
	
}

MemberFunCorr<-function(AgeDist){
	
###################################################
#	                                               #
#	This function is developed to find membership  #
#    values of given Age Distribution              #
#	                                               #
####################################################


###################################################
#	                                               #
#	Input: AgeDist same age distribution for one   #
#   analysis region                                #
#	                                               #
####################################################

# We have to normalize the given distibution with respect to
# total people in a given district
TotPeop<- sum(AgeDist)	;

# Normalized distribution
NormDist<-AgeDist/TotPeop;

# let's define membership
f<- numeric(3)
# f[1]: Low
# f[2]: Medium
# f[3]: High

# Let's calculate Low value for given distribution

for(i in 1:6){

if(i==6){
	f[1]<-f[1]+	0.5*NormDist[i];
} 	else
   f[1]<-f[1]+NormDist[i];
	
}


# Let's calculate Medium value for given distribution

for(i in 6:9){

if(i==6){
	f[2]<-f[2]+	0.5*NormDist[i];
} 	else if(i==9){
	f[2]<-f[2]+	0.5*NormDist[i];
}   else 
   f[2]<-f[2]+NormDist[i];
	
}

	
# Let's calculate High value for given distribution

for(i in 9:10){

if(i==9){
	f[3]<-f[3]+	0.5*NormDist[i];
} 	else 
	f[3]<-f[3]+	NormDist[i];
	
	}
	
	## it only returns f vector
	
	return(f)
	
}


