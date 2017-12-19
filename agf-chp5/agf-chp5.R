#######################################################################
#
#  SIMULATE SIMPLE RANDOM SAMPLE AND AGE LENGTH KEY (ALK) FOR A 
#  FISH POPULATION GIVEN LENGTH AT INFINITY (LINF) AS THE INPUT
#  THE SCRIPT BELOW GENERATES THE REFERENCE POPUALTION AND SAMPLES
#  THE POPULATION DESCRIBED IN BOX 5.2 IN AGE AND GROWTH OF FISHES
#  MICHAEL.COLVIN@MSSTATE.EDU TO REPORT ERRORS OR QUESTIONS 
#
#######################################################################


source("agf-chp5-functions.R")# READ IN THE FUNCTION 
library(FSA) ## NEEDED IN FUNCTION
library(pbapply) ## for progress bar




nsims=100 ## HOW MANY STOCHASTIC REPLICATES
## THIS FUNCTION RUNS THE SAMPLE SIMULATION FOR 
## THE SPECIFICIED NUMBER OF REPLICATES
## IT IS NOT A PARTICULARLY FAST FUNCTION SO START WITH A 
## FEW SIMS BEFORE SET UP TO RUN 1000S
sims<-pblapply(c(1:nsims) ,function(x)
    {
    tmp<-sim_sample(linf=40,    # length at infinity in cm
        populationSize=100000,  # total population size
        sampleSize=100,         # how many fish sampled
        bin=2.5,                # bin size in cm
        bin_n=10)               # how many fish per bin
    tmp[[1]]$rep<-x
    tmp[[2]]$rep<-x
    tmp[[3]]$rep<-x
    return(tmp)
    })

# CALCULATE PERFORMANCE
## MEAN ABSOLUTE DEVIATION (MAD)
## ABSOLUTE PERCENT DEVIATION (APD)
## MEAN ABSOLUTE PERCENT DEVIATION (MAPD)    

    
## MEAN AGE
ma<-lapply(sims,function(x) x[[1]] )  
ma<-do.call("rbind",ma)    
ma$MAD_srs<- abs(ma$true-ma$srs)
ma$MAD_alk<- abs(ma$true-ma$alk)
ma$APD_srs<- abs((ma$true-ma$srs)/ma$true)*100
ma$APD_alk<- abs((ma$true-ma$alk)/ma$true)*100


## AGE FREQUENCY
af<-lapply(sims,function(x) x[[2]] )  
af<-do.call("rbind",af)
af$MAD_srs<- abs(af$true-af$srs)
af$MAD_alk<- abs(af$true-af$alk)
af$MAPD_srs<- abs((af$true-af$srs)/af$true)*100
af$MAPD_alk<- abs((af$true-af$alk)/af$true)*100

### SUMMARY
af<-aggregate(cbind(MAD_srs,MAD_alk,MAPD_srs,MAPD_alk)~rep,
    data=af,
    FUN=mean)

    
## MEAN LENGTH AT AGE
mla<-lapply(sims,function(x) x[[3]])   
mla<-do.call("rbind",mla) 
mla$MAD_srs<-abs(mla$true-mla$srs)
mla$MAD_alk<-abs(mla$true-mla$alk)
mla$MAPD_srs<-abs((mla$true-mla$srs)/mla$true)*100
mla$MAPD_alk<-abs((mla$true-mla$alk)/mla$true)*100

### SUMMARY
mla<-aggregate(cbind(MAD_srs,MAD_alk,MAPD_srs,MAPD_alk)~rep,
    data=mla,
    FUN=mean)