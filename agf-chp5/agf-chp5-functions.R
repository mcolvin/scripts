

sim_sample<-function(linf=40,
    populationSize=100000,
    sampleSize=100,
    bin=2.5,
    bin_n=10)# bin size in cm
    {
    # STAGE 1

    ## SET TOTAL POPULATION SIZE
    N<- populationSize
    
    ## VBGF 
    Linf<- linf
    k<- exp(1.2434+-0.66*log(Linf))# DEVELOPED FROM FISHBASE
    t0<- 0	
    
    ## MAXIMUM AGE AS A FUNCTION OF VBGF PARAMETERS
    maxAge<-round(log(1-(((Linf*0.95)/Linf)))/-k + t0)
    
    ## RELATIONSHIP OF ANNUAL MORTALITY AND MAXIMUM AGE
    b = -1.01
    a = exp(1.46)
    A<- 1-exp(-a*maxAge^b )
    
    ## SET CUMULATIVE SURVIVAL ASSUMING EQUILIBRIUM
    ## AND DRAW FROM A MULTINOMIAL DISTRIBUTION
    survival<- cumprod(rep(1-A,maxAge))
    dat<-data.frame(age=c(1:maxAge),
        N=rmultinom(1,N,survival))
    
    ## ASSIGN EXPECTED LENGHT AT AGE GIVEN VBGF PARAMETERS
    dat$exp_length<- Linf  * (1 - exp(-k * (dat$age-t0))) 
    
    ## SET UP A POPULATION TO SAMPLE
    pop<- as.data.frame(lapply(dat,function(x) rep(x,dat$N)))
    pop$age_f<- as.factor(pop$age)# MAKE A FACTOR
    
    ## ASSIGN LENGTH TO EXPANDED DATASET AND 
    ## ADD NOISE
    pop$len_true<- rnorm(length(pop$exp_length),
        pop$exp_length,pop$exp_length*0.2)
        
    ## ROUND TO NEAREST MM FOR OBSERVED LENGTH
    pop$len_obs<- round(pop$len_true,0)
    
    ## ASSIGN LENGTH BINS
    bins<-seq(from=0,to=max(pop$len_obs)+bin,
        by=bin)
    pop$bin<-cut(pop$len_obs,
        breaks=bins,
        labels=(bins+(bin*0.5))[-length(bins)],
        right=TRUE)
    pop$bin<- factor(pop$bin, levels=c((bins+(bin*0.5))[-length(bins)]))
    
    
    # TRUE VALUES: POPULATION LEVEL
    
    ## CALCULATE TRUE AGE FREQUENCY
    ## OF THE POPULATION
    af<-as.data.frame(prop.table(table(pop$age)))
    
    ## CALCULATE TRUE MEAN AGE
    ## OF THE POPULATION 
    ma<- mean(pop$age)
    
    ## CALCULATE TRUE MEAN LENGTH AT AGE
    ## OF THE POPULATION      
    mla<-aggregate(len_obs~age_f,pop,
        FUN=mean)

        
    # SAMPLE THE POPULATION
    ## SIMPLE RANDOM SAMPLE
    srs<- pop[sample(c(1:N),sampleSize,replace=FALSE),]
    
    ## ESTIMATES
    ### AGE FREQUENCY
    srs_af<-as.data.frame(prop.table(table(srs$age_f)))
    
    ### MEAN AGE
    srs_ma<- mean(srs$age)
    
    ### MEAN LENGTH AT AGE
    srs_mla<- tapply(srs$len_obs,srs$age_f,mean)
    
    ## AGE LENGTH KEY
    alk<- srs[order(srs$bin),]
    alk$num <- ave(alk$len_obs, alk$bin, FUN = seq_along)
    alk<-subset(alk,num<=bin_n)
    len.n <- xtabs(~bin,data=srs)
	raw <- xtabs(~bin+age_f,data=alk)
	key <- prop.table(raw, margin=1)
    key[is.na(key)]<-0
    
    ## ESTIMATES
    ### AGE FREQUENCY
    alk_af<-c(prop.table(t(key) %*% table(srs$bin)))
    
    ### MEAN AGE
    alk_ma<-c(prop.table(t(key) %*% table(srs$bin)))
    alk_ma<-sum(as.numeric(colnames(key))*alk_ma) 

    ### CALCULATE MEAN LENGTH AT AGE FROM BETTOLI AND MIRANDA
	alk_mla<- suppressWarnings(FSA::alkMeanVar(key,
        len_obs~bin+age_f,
        alk,len.n)[,c(1:2)])
        
    ## CALCULATE UP MAD, APD, MAPD
    meanAge<-data.frame(true=ma,
        srs=srs_ma,
        alk=alk_ma)
    ageFrequency<-data.frame(age=af[,1],
        true=af[,2],
        srs=srs_af[,2],
        alk=alk_af)
    meanLengthAtAge<-data.frame(age=mla[,1],
        true=mla[,2],
        srs=srs_mla,
        alk=alk_mla[,2])
    return(list(meanAge=meanAge,
        ageFrequency=ageFrequency,
        meanLengthAtAge=meanLengthAtAge))
    }