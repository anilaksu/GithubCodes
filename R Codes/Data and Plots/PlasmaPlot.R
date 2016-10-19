# this function created to plot .csv files obtianed by simulation for 
# electric circuits simulations for plasma torch 

resistorData <-  read.csv("EnergyTransferred.csv");
#first column of the data is time in second
time <- resistorData[,1]
#second column of the data is current in ampere
current <- resistorData[,2]
#third column of data is potential in volt
volt <- resistorData[,3]
#forth column of data is power in watt
power <- resistorData[,4]

plot( time[1:45], current[1:45], type="l", col="red",xlab="seconds",ylab="Ampere")
par(new=T)
plot( time[1:45], power[1:45], type="l", col="green" ,axes=FALSE
,ylab="")
legend('topright', "Current in Amperes","Power in Watt" , 
   lty=1, col=c('red',  'green'), bty='n', cex=.75)


#print(length(time));