#update date: 11/17/2019

##############DEPENDENCIES
usePackage <- function(p) 
{
  p = deparse(substitute(p)) 
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}


deparse(substitute(lmerTest)) 

usePackage(lmerTest)
usePackage(boot)
usePackage(LMERConvenienceFunctions)
usePackage(simpleboot)
usePackage(broom) #for visualizations
usePackage(merTools)
usePackage(foreach)
usePackage(doParallel)
usePackage(gridExtra)
usePackage(tcltk)
usePackage(data.table)
usePackage(broom) #for residual plots
usePackage(ggplot2) #for plot visualization
usePackage(cowplot) #for visualization
usePackage(plyr)
usePackage(zoo) #for rollmean function in integration



#################################################################################
####################FUNCTIONS ##################################################
################################################################################

###Preprocesses raw data by parsing and filtering inputs


preprocess = function(rawdat2){
  #rawdat2 = unsorted
  rawdat = lapply(c(1:length(rawdat2)), function(y){
    dt = rawdat2[[y]]
    if(ncol(dt)>14){
      dt$R2_S1[is.na(dt$R2_S1)] = .001
      dt$R2_S2[is.na(dt$R2_S2)] = .001
    }
    t = dt[complete.cases(dt),] #filter out any rows containing NAs
    trash = dt[!complete.cases(dt),] #filter out any rows containing NAs
    write.table(trash, file = paste(dirs,"trash.csv",sep="/"),row.names=FALSE,col.names=TRUE,sep=",")
    if (ncol(dt) < 15){
    t$Validation_Training[t$Validation_Training=="V"] = "T" #makes everything change to T
    t$Validation_Training[t$Validation_Training==TRUE] = "T"} #makes everything change to T
    return(t)})
        
  #test_dat = rawdat[which(sapply(rawdat, function(y) "V" %in% y$Validation_Training))] #extract testing/validation data
  trn_dat = rawdat[which(sapply(rawdat, function(y) "T" %in% y$Validation_Training))] #extract calibration data; r converts "T" to "TRUE 
  meas_dat = rawdat[which(sapply(rawdat, function(y) "TRUE" %in% (y$WL > 0)))] #extract measurement data
  ###check data###
  ##add whale field and check for missing data
  colnames=c("Name","Sighting","Flight","Whale","Flight_Date","Date","Original_Filename","Output_Image_Name","Baro_Alt","Adj_Baro_Alt","GSD","Focal_Length","Cam_Pix_mm",
              "R2_S1","R2_S2","SA","WL","FW","MW","OW","Pcnt20Wd","Pcnt30Wd","Pcnt40Wd","Pcnt50Wd","Pcnt60Wd")
  meas_dat = lapply(c(1:length(meas_dat)), function(yy) {
    y = meas_dat[[yy]]
        if (("Whale" %in% names(y))=="FALSE"){
          y$Whale = "A"
        }
        y$Whale = gsub('([[:punct:]])|\\s+','_',y$Whale)
        y$Name = factor(paste(y$Flight_Date,y$Flight,y$Whale,sep="_"))
        y = y[y$GSD>0,] #filters rows with bad GSD
        y = y[y$Baro_Alt>0,] #filters rows with bad altitude
        y = y[y$WL>10,] #filters rows with bad length measurement
        y = y[y$OW>10,] #filters rows with bad width measurement
        y$Date = NA
        ##iterate through dates
        for (iii in 1:nrow(y)){
          spltdt = strsplit(y$Flight_Date[iii], "/")
          month_val = as.numeric(spltdt[[1]][1])
          day_val = as.numeric(spltdt[[1]][2])
          year_val = as.numeric(spltdt[[1]][3])
          try(if (year_val < 2000){
            year_val = year_val + 2000
          },TRUE)
          newdt = paste(month_val,day_val,year_val,sep="/")
          newdt2 = as.character(as.Date(newdt,"%m/%d/%Y"))
          y$Date[iii] = newdt2
        }
        y$Flight_Date = paste(y$Flight,y$Date,sep="_")
        y = na.omit(y) #removes any row with NA values
        y = y[,colnames,with=FALSE]
        y$GSD = y$Adj_Baro_Alt * (y$Cam_Pix_mm/y$Focal_Length)
        y$BL2 = (y$WL*0.4)^2
        y$BAI = y$SA *100 / y$BL2
        data.frame(y)
  })

  meas_dat2 = do.call("rbind",meas_dat)
  nameslist = unique(meas_dat2$Name)
  meas_dat3 = lapply(nameslist,function(x){
    meas_dat2[meas_dat2$Name == x,]
  })
  meas_dat = meas_dat3
  trn_dat = lapply(trn_dat, function(y) {
    y = y[y$GSD>0,] #filters rows with bad GSD
    y = y[y$Baro_Alt>0,] #filters rows with bad altitude
    y = y[y$Object_Length_Pixels>10,] #filters rows with bad length measurement
    y = na.omit(y)
    y$Flight = factor(y$Flight)
    y$Sighting = factor(y$Sighting)
    y$Date = NA
    ##iterate through dates
    for (iii in 1:nrow(y)){
      spltdt = strsplit(y$Flight_Date[iii], "/")
      month_val = as.numeric(spltdt[[1]][1])
      day_val = as.numeric(spltdt[[1]][2])
      year_val = as.numeric(spltdt[[1]][3])
      if (year_val < 2000){
        year_val = year_val + 2000
      }
      newdt = paste(month_val,day_val,year_val,sep="/")
      newdt2 = as.character(as.Date(newdt,"%m/%d/%Y"))
      y$Date[iii] = newdt2
    }
    y$Flight_Date = paste(y$Flight,y$Date,sep="_")
    y$GSD = y$Adj_Baro_Alt * (y$Cam_Pix_mm/y$Focal_Length) #calc GSD from adj baro alt
    data.frame(y) #convert to data frame so program behaves normally
  })
  ##name the lists according to sighting and flight-date
  names1 = unique(unlist(lapply(trn_dat, function(y) paste(y$Date,y$Flight_Date,sep="_"))))
  names(trn_dat) = names1

  return(list(meas_dat,trn_dat))
}


###Different modeling Methods for training data (e.g. GSD calculation)
modelbattery = function(x,method){
  #note method 1 is NO correction
  mod2 = NA
  mod3 = NA
  mod4 = NA
  mod5 = NA
  mod6 = NA
  mod7 = NA
  
  if (4 %in% method){
    mod4 <- lm(log10(Object_Length_Pixels) ~ log10(Adj_Baro_Alt), data = x) #Method 4 smoothing
    mod6 <- lm(Object_Length_mm/1000/10^predict(mod4) ~ Adj_Baro_Alt, data = x) #Method 4 correction
  }
  
  if (5 %in% method){
    mod5 <- lmer(log10(Object_Length_Pixels) ~ log10(Adj_Baro_Alt) + (1 | Flight_Date), data = x,REML = F) #Method 5 smoothing
    mod7 <- lmer(Object_Length_mm/1000/10^predict(mod5) ~ Adj_Baro_Alt + (1| Flight_Date), data = x,REML = F) #Method 5 correction
  }
  
  if (3 %in% method){
    mod3 <- lmer(Object_Length_mm/1000/Object_Length_Pixels ~ Adj_Baro_Alt + (1 | Flight_Date),data = x,REML = F) #Method 3
  }
  
  if (2 %in% method){
    mod2 <- lm(Object_Length_mm/1000/Object_Length_Pixels ~ Adj_Baro_Alt, data=x) #Method 2
  }
  return(list(mod2,mod3,mod4,mod5,mod6,mod7))
}

###Bootstrapped prediction with lmer
lmerboots = function(model,new.data){
  #cl<-makeCluster(6) #  #set parallel backend to use 8 processors
  #registerDoParallel(cl)
  pred = predictInterval(model,new.data,n.sims=2000,stat=c("mean"),level=0.68,which=c("full"),.parallel=FALSE)
}


## 
CIboot = setClass("CIboot", slots = c(lwr="numeric", upr="numeric",mean="numeric"))

###bootstrapped CIs####
# ####boot strapping confidence intervals and error propogation
# a function to perform bootstrapping

bootSE = function(raw.data, B=1000){
  # this function will take 1,000 (by default) bootsamples calculate the mean of
  # each one, store it, & return the bootstrapped sampling distribution of the mean

  boot.dist = vector(length=B)     # this will store the means
  N         = length(raw.data)     # this is the N from your data
  for(i in 1:B){
    boot.sample  = sample(x=raw.data, size=N, replace=TRUE)
    boot.dist[i] = mean(boot.sample)
  }
  boot.dist = sort(boot.dist)
  min = boot.dist[25]
  max = boot.dist[976]
  mean = mean(raw.data)
  val = CIboot(lwr=min,upr=max,mean=mean)
  return(val)
}


###Bootstrapped prediction with lm
lmboots = function(model,new.data,old.data){
  #references: http://www.statoo.com/en/publications/bootstrap_scgn_v131.pdf
  #http://stats.stackexchange.com/questions/64813/two-ways-of-using-bootstrap-to-estimate-the-confidence-interval-of-coefficients
  #http://blog.minitab.com/blog/adventures-in-statistics-2/when-should-i-use-confidence-intervals-prediction-intervals-and-tolerance-intervals
  #http://www.ats.ucla.edu/stat/r/faq/boot.htm
  #http://stats.stackexchange.com/questions/17065/estimate-confidence-interval-of-mean-by-bootstrap-t-method-or-simply-by-bootstra
  #formula is directly from ?boot
  
  n.sim = 1000
  lmdat.data <- data.frame(old.data, resid = residuals(model), fit = fitted(model)) #collect pertinent metrics
  new.fit <- predict(model, new.data) #make prediction
  lmdat.fun <- function(dat, inds, i.pred, fit.pred, x.pred)
  {
    y = dat$fit+dat$resid[inds]
    lm.b <- lm(as.formula(paste("y~", paste(deparse(model$call[[2]][[3]])))),data=dat)
    pred.b <- predict(lm.b, x.pred) #make 
    c(pred.b - (fit.pred + dat$resid[i.pred]),coef(lm.b))
  }
  lmdat.boot <- boot(lmdat.data, lmdat.fun, R = n.sim, m = 1, 
                     fit.pred = new.fit, x.pred = new.data)
  pred = data.frame(new.fit,t(sapply(c(1:length(new.fit)),function(x) new.fit[x] - sort(lmdat.boot$t[,x])[c(.16*n.sim, .84*n.sim)])))
  names(pred) = c("fit","upr","lwr")
  return(pred)
}

### Function to handle prediction of different model types (LM vs LMER)
predictions = function(model,new.data,old.data){
  if(is.na(model)){
    preds = NA
  }else if (attributes(model)$class == "lm"){
      preds = lmboots(model,new.data,old.data)
  }else{preds = lmerboots(model,new.data)}
  return(preds)
}

###Length Estimation for comparing training model Methods
calc_length = function(dat,predictions,method){
  ids=c("Flight_Date","Date","Method","ID","Name")
  uid = paste(dat$Flight_,dat$Date,sep="_")  
  #set everything as NA value, then replace by method
  s1 = NA
  s2 = NA
  s3 = NA
  s4 = NA
  s5 = NA
  
  if (4 %in% method){
    s4 = data.frame(dat$Flight_Date,dat$Date,"4",uid,"Cal_Obj",dat$Object_Length_Pixels* predictions[[5]])
    names(s4)[1:5]=ids; s4[1:5]=lapply(s4[1:5],function(x) factor(x))
    names(s4)[6:8] = c("Cal.fit","Cal.upr","Cal.lwr")
  }
  
  if (5 %in% method){
    s5 = data.frame(dat$Flight_Date,dat$Date,"5",uid,"Cal_Obj",dat$Object_Length_Pixels * predictions[[6]])
    names(s5)[1:5]=ids; s5[1:5]=lapply(s5[1:5],function(x) factor(x))
    names(s5)[6:8] = c("Cal.fit","Cal.upr","Cal.lwr")
  }
  
  if (2 %in% method){
    #calibrated
    s2 = data.frame(dat$Flight_Date,dat$Date,"2",uid,"Cal_Obj",predictions[[1]] * dat$Object_Length_Pixels)
    names(s2)[1:5]=ids; s2[1:5]=lapply(s2[1:5],function(x) factor(x))
    names(s2)[6:8] = c("Cal.fit","Cal.upr","Cal.lwr")
  }
  
  if (3 %in% method){
    s3 = data.frame(dat$Flight_Date,dat$Date,"3",uid,"Cal_Obj",predictions[[2]] * dat$Object_Length_Pixels)
    names(s3)[1:5]=ids; s3[1:5]=lapply(s3[1:5],function(x) factor(x))
    names(s3)[6:8] = c("Cal.fit","Cal.upr","Cal.lwr")
  }
  
  if (1 %in% method){
    #no adjustment
    s1 = data.frame(dat$Flight_Date,dat$Date,"1",uid,"Cal_Obj",dat$Object_Length_Pixels * dat$GSD) #replicates 3x so it has same # of columns as the others
    names(s1)[1:5]=ids; s1[1:5]=lapply(s1[1:5],function(x) factor(x))
    names(s1)[6] = "Cal.fit"
    s1$Cal.upr = NaN
    s1$Cal.lwr = NaN
  }
    
  return(list(s1,s2,s3,s4,s5))
}

#determines if individual observations are more extreme than the others
#only remove very extreme values because SD only changes as a function of the observation being examined rather than changing SD as a function of previously removed values
filt_dat = function(dat,idx,filt_intens){
  dat[idx] = sapply(dat[idx],function(x){
    tmpdf = data.frame(x * dat$GSD)
    tmpdf2 = tmpdf
    for (ii in 1:dim(tmpdf)[1]){
      testdat = tmpdf[,1][! tmpdf[,1] %in% tmpdf[ii,]]
      sds = sd(testdat,na.rm=TRUE)
      means = mean(testdat,na.rm=TRUE)
      num_obs2 = dim(tmpdf)[1]
      min = means - abs(qt(filt_intens/2,num_obs2)) * sds
      max = means + abs(qt(filt_intens/2,num_obs2)) * sds
        #min = means - abs(qt(0.025,5)) * sds
      #max = means + abs(qt(0.025,5)) * sds
      tmpdf2[ii,][tmpdf2[ii,]>max]=NA
      tmpdf2[ii,][tmpdf2[ii,]<min]=NA
    }
    tmpdf = tmpdf2
    x[is.na(tmpdf)]=NA 
    x
  })
  dat = dat[complete.cases(dat),] #removes any row containing NA
  return(dat)
}


#determines if individual observations are more extreme than the others
#uses more intense filter by removing extreme values and reducing SD and repeating 
filt_valdat = function(dat,idx){
  dat[idx] = sapply(dat[idx],function(x){
    tmpdf = data.frame(x * dat$GSD)
    for (ii in 1:dim(tmpdf)[1]){
      testdat = tmpdf[,1][! tmpdf[,1] %in% tmpdf[ii,]]
      sds = sd(testdat,na.rm=TRUE)
      means = mean(testdat,na.rm=TRUE)
      min = means - abs(qt(filt_intens/2,dim(tmpdf)[1])) * sds
      max = means + abs(qt(filt_intens/2,dim(tmpdf)[1])) * sds
      tmpdf[ii,][tmpdf[ii,]>max]=NA
      tmpdf[ii,][tmpdf[ii,]<min]=NA
    }
    x[is.na(tmpdf)]=NA 
    x
  })
  dat = dat[complete.cases(dat),] #removes any row containing NA
  return(dat)
}
###smoothing model specific to each morphometric parameter
whalesmooth = function(dat,idx){
  #independent smoothing --this model is specific to each 'object'
  lma <- lapply(dat[idx],function(x) lm(log10(x) ~ log10(Adj_Baro_Alt), data = dat)) #FE model
}

###Calculations on actual whale morphometric parameters
morphometrics = function(predictions,dat,idx,cntr,method){ #calculates lengths from estimates

  # predictions = whalepreds[[ii]]
  # dat = meas_dat[[ii]]
  # idx = idx
  # cntr = ii
  # method = method
  
  ids=c("Flight","Date","Whale","Method","ID")
  #uid = paste(dat$Flight[1],dat$Date[1],dat$Whale[1],sep="_")
  uid = paste("Whale",cntr,sep="_")
  p2 = NA
  p3 = NA
  p1 = NA
  p6 = NA
  p7 = NA
  
  if (2 %in% method){
  #calibrated only;
    p2 = data.frame(dat$Flight_Date,dat$Date,dat$Whale,"2",uid,lapply(dat[idx][2:10],function(x) x * predictions[[1]]))
    ndf = data.frame(dat[idx][[1]] * predictions[[1]]^2,dat[idx][[11]],NaN,NaN,dat[idx][[12]],NaN,NaN)
    names(ndf) = c("SA.fit","SA.upr","SA.lwr","BL2.fit","BL2.upr","BL2.lwr","BAI.fit","BAI.upr","BAI.lwr")
    p2 = cbind(p2,ndf)
    names(p2)[1:5]=ids; p2[1:5]=lapply(p2[1:5],function(x) factor(x))
  }
  
  if (3 %in% method){
    p3 = data.frame(dat$Flight_Date,dat$Date,dat$Whale,"3",uid,lapply(dat[idx][2:10],function(x) x * predictions[[2]]))
    ndf = data.frame(dat[idx][[1]] * predictions[[2]]^2,dat[idx][[11]],NaN,NaN,dat[idx][[12]],NaN,NaN)
    names(ndf) = c("SA.fit","SA.upr","SA.lwr","BL2.fit","BL2.upr","BL2.lwr","BAI.fit","BAI.upr","BAI.lwr")
    p3 = cbind(p3,ndf)
    names(p3)[1:5]=ids; p3[1:5]=lapply(p3[1:5],function(x) factor(x))
  }
  
  if (4 %in% method){
  #calibrated from smoothed+calibrated model
    p6 = data.frame(dat$Flight_Date,dat$Date,dat$Whale,"4",uid,lapply(dat[idx][2:10],function(x) (x * predictions[[3]]))) #based on lmer
    ndf = data.frame(dat[idx][[1]] * predictions[[3]]^2,dat[idx][[11]],NaN,NaN,dat[idx][[12]],NaN,NaN)
    names(ndf) = c("SA.fit","SA.upr","SA.lwr","BL2.fit","BL2.upr","BL2.lwr","BAI.fit","BAI.upr","BAI.lwr")
    p6 = cbind(p6,ndf)
    names(p6)[1:5]=ids; p6[1:5]=lapply(p6[1:5],function(x) factor(x))
  }
  
  if (5 %in% method){
    p7 = data.frame(dat$Flight_Date,dat$Date,dat$Whale,"5",uid,lapply(dat[idx][2:10],function(x) (x * predictions[[4]]))) #based on lm
    names(p7)[6:8] = c("WL.fit","WL.upr","WL.lwr")
    ndf = data.frame(dat[idx][[1]] * predictions[[4]]^2,dat[idx][[11]],NaN,NaN,dat[idx][[12]],NaN,NaN)
    names(ndf) = c("SA.fit","SA.upr","SA.lwr","BL2.fit","BL2.upr","BL2.lwr","BAI.fit","BAI.upr","BAI.lwr")
    p7 = cbind(p7,ndf)
    names(p7)[1:5]=ids; p7[1:5]=lapply(p7[1:5],function(x) factor(x))
  }
  
  if (1 %in% method){
  #no adjustment
  p1 = data.frame(dat$Flight_Date,dat$Date,dat$Whale,"1",uid,lapply(dat[idx][2:10],function(x) x * dat$GSD));
  ndf = data.frame(dat[idx][[1]] * dat$GSD^2,dat[idx][[11]],dat[idx][[12]])
  names(ndf) = c("SA","BL2","BAI")
  p1 = cbind(p1,ndf)
  names(p1)[1:5]=ids; p1[1:5]=lapply(p1[1:5],function(x) factor(x))
  idx4 = 6:length(names(p1))
  p1[,paste(names(p1[idx4]),"upr",sep=".")] =  NA
  p1[,paste(names(p1[idx4]),"lwr",sep=".")] =  NA
  names(p1)[idx4] = paste(names(p1[idx4]),"fit",sep=".")
  }

  return(list(p1,p2,p3,p6,p7))
}

#model LENGTH 
lengthmdl = function(pixpreds,gsdpreds,old.dat){
  newdat = do.call("rbind",lapply(seq(1:100),function(y){
    OL = old.dat$Object_Length_mm[1]/1000 #object length in mm
    newOL = y/2
    pix = 10^pixpreds$fit/OL * newOL #convert to pix/meter then multiply
    GSD = gsdpreds$fit
    df=data.frame(pix,newOL,GSD)
    names(df) = c("pix","newOL","GSD")
    df}))
  nlm = lm(newOL ~ pix * GSD, data=newdat)
}

uncormeans = function(x){
  ids = unique(x$ID)
  dat = do.call("rbind",lapply(ids,function(y){ #loop by id
    df2 = subset(x,ID == y)
    lwrpos = which(grepl('.lwr', names(df2))) 
    uprpos = which(grepl('.upr', names(df2))) 
    datpos = which(grepl('.fit', names(df2)))
    df3 = data.frame(mapply(function(a,b,c){
      if (anyNA(df2[[a]],recursive=TRUE)==TRUE){
        df = data.frame(matrix(data=NaN,nrow=dim(df2)[1],ncol=3)) #keeps NaNs from going into bootSE
        names(df) = c("dat.upr","dat.lwr","dat.mean")
        df
      }else{
        if (dim(df2)[1]==1){
          df = data.frame(matrix(data=NaN,nrow=dim(df2)[1],ncol=3)) #keeps NaNs from going into bootSE
          names(df) = c("dat.upr","dat.lwr","dat.mean")
          df$dat.mean = df2[[a]]
          df
      }else{
        dat = bootSE(df2[[a]])
        df = data.frame(dat@upr,dat@lwr,dat@mean)
        df
      }}},datpos,uprpos,lwrpos))
    names(df3) = names(df2[datpos])
    df3 = unlist(df3)
    newnames = unlist(lapply(names(df3),function(x){
      paste(strsplit(x,split="\\.")[[1]][1],strsplit(x,split="\\.")[[1]][4],sep=".")
    }))
    names(df3) = newnames
    df3 = data.frame(as.list(df3))
    df3 = cbind(df2[1,1],unique(df2[2]),unique(df2[3]),unique(df2[4]),unique(df2[5]),df3)
    names(df3)[1] = "Flight_Date"
    df3
  }))
  return(dat)
}

cormeans = function(x){
  ids = unique(x$ID)
  dat = do.call("rbind",lapply(ids,function(y){ #loop by id
    df2 = subset(x,ID == y)
    lwrpos = which(grepl('.lwr', names(df2))) 
    uprpos = which(grepl('.upr', names(df2))) 
    datpos = which(grepl('.fit', names(df2)))
    df3 = data.frame(mapply(function(a,b,c){
      if (anyNA(df2[[a]],recursive=TRUE)==TRUE){
        df = data.frame(matrix(data=NaN,nrow=1,ncol=3)) #keeps NaNs from going into bootSE
        names(df) = c("dat.upr","dat.lwr","dat.mean")
        df
      }else{
        if (a>33 & dim(df2)[1]>1){
          dat = bootSE(df2[[a]])
          df = data.frame(dat@mean,dat@upr,dat@lwr)
          df
        }else{
          if (dim(df2)[1]==1){
            df = data.frame(matrix(data=NaN,nrow=dim(df2)[1],ncol=3)) #keeps NaNs from going into bootSE
            names(df) = c("dat.upr","dat.lwr","dat.mean")
            df$dat.mean = df2[[a]]
            df
        }else{
        dat = bootSE(df2[[a]])
        dat2 = bootSE(df2[[b]])
        dat3 = bootSE(df2[[c]])
        df = data.frame(dat@mean,dat2@upr,dat3@lwr)
        df
        }}}},datpos,uprpos,lwrpos))
    names(df3) = names(df2[datpos])
    df3 = unlist(df3)
    newnames = unlist(lapply(names(df3),function(x){
      paste(strsplit(x,split="\\.")[[1]][1],strsplit(x,split="\\.")[[1]][4],sep=".")
    }))
    names(df3) = newnames
    df3 = data.frame(as.list(df3))
    df3 = cbind(df2[1,1],unique(df2[2]),unique(df2[3]),unique(df2[4]),unique(df2[5]),df3)
    names(df3)[1] = "Flight_Date"
    df3
  }))
  return(dat)
}

my.file.rename <- function(from, to) {
  #from https://stackoverflow.com/questions/10266963/moving-files-between-folders
  todir <- dirname(to)
  if (!isTRUE(file.info(todir)$isdir)) dir.create(todir, recursive=TRUE)
  file.rename(from = from,  to = to)
}

#input and preprocessing 
inproc = function(dirs,filt_cal,filt,filt_intens){
  print("Reading and Preprocessing raw data")
  
  meas_dat.drop = NA #placeholder in case meas_dat.drop doesn't exist
  
  if(!dir.exists(dirs)){
    stop(sprintf("The data directory %s does not exist", dirs))
  }
  
  #see if old procesing files exist and move them if they do
  if(file.exists("Validation.csv")){
    my.file.rename(from = "Validation.csv",
                   to = paste(dirs,"processing_archive","Validation.csv",sep="/"))
  }
  if(file.exists("Morph_Measures.csv")){
    my.file.rename(from = "Morph_Measures.csv",
                   to = paste(dirs,"processing_archive","Morph_Measures.csv",sep="/"))
  }
  datfiles = {}
  datfiles = list.files(path = dirs ,pattern=".csv",full.name=TRUE,recursive=TRUE) #lists all .csv files in the selected folder
  if (length(datfiles)<1){sprintf("There are no datfiles in directory: %s ", dirs)
    stop(sprintf("There are no datfiles in directory: %s ", dirs))
  }
  
  z2 = c()
  for(z in c(1:length(datfiles))){
    dt = datfiles[[z]]
    if (grepl("Morph_Measures.csv",dt[1])==TRUE){
      #datfiles2 = datfiles[-z]
      z2 = append(z2,z)
    }else if (grepl("Validation.csv",dt[1])==TRUE){
      #datfiles2 = datfiles[-z]
      z2 = append(z2,z)
    }}
  
  #only remove files if they exist
  if (length(z2)>1){
    datfiles = datfiles[-z2]
  }
  
  if (length(datfiles) < 1){
    stop("No measurement files were detected")
  }else{sprintf("There were %s measurement files detected for processing", length(datfiles))}
  #begin sorting data
  unsorted = lapply(datfiles,function(i){fread(i)}) #read in the data
  datlists = preprocess(unsorted) #process raw data
  
  ##training and validation data
  valdat=do.call("rbind",datlists[[2]])
  if (length(datlists[[2]]) <1){
    print("WARNING!! Zero calibation measurements were detected; all data processed under method 1. Make sure that .csv files with following prefix are present: whale_measurements_")
  }
  idx = which(names(valdat)=="Object_Length_Pixels")
  if (filt_cal == "yes"){
    valdat = filt_valdat(valdat,idx) #more extreme filtering than on measurement data
    
  }
  
  
  ##filtering measurement data
  meas_dat = datlists[[1]]
  if (length(datlists[[1]]) <1){
    stop("Zero measured whales were detected; make sure that .csv files with following prefix are present: whale_measurements_")
  }
  idx3 = which(names(meas_dat[[1]])=="WL"):which(names(meas_dat[[1]])=="WL")
  if (filt == "yes"){
    meas_dat = lapply(meas_dat,function(x) filt_dat(x,idx3,filt_intens))}
  newdat = do.call("rbind",meas_dat) #collapse the data
  
  #remove measurement not having matching calibration data"
  newdat$fltdt = paste(newdat$Flight,newdat$Date,sep="_")
  valdat$fltdt = paste(valdat$Flight,valdat$Date,sep="_")
  test2 = match(newdat$fltdt,valdat$fltdt) 
  newdat.drop = newdat[is.na(test2),]
  newdat.keep = newdat[!is.na(test2),]
  
  if (nrow(newdat.drop )>0){
    print("Unique number of dates and flights in calibration data are not the same as the unique number of dates and flights in the measurement data")
  }
  
  meas_dat.keep = vector(0,mode="list")
  meas_dat.keep = lapply(unique(newdat.keep$Name),function(pp){
    subset(newdat.keep,Name == pp)
  })
  
  meas_dat.drop = vector(0,mode="list")
  if (nrow(newdat.drop)>0){
    meas_dat.drop = lapply(unique(newdat.drop$Name),function(pp){
      subset(newdat.drop,Name == pp)
    })
  }
  
  meas_dat = meas_dat.keep
  
  if (length(datlists[[2]]) <1){
      valdat = NA
  }
  
  
  return(list(meas_dat,meas_dat.drop,valdat))
}

#function that controls calibration data processing
calproc = function(valdat,method){
  
  #valdat = inprocdat[[3]] #debugging line; comment when not using
  
  print("Processing Calibration Data and Creating Models")
  modlist = modelbattery(valdat,method) #creates model Methods based on validation data
  preds = lapply(modlist,function(x) predictions(x,valdat,valdat)) #makes bootstrapped predictions of means and prediction intervals
  ##NOTE## Warning "is.na(model)..." is ok because if there was an NA, there wouldn't be a warning and behavior would be as expected
  lengths = calc_length(valdat,preds,method) #calculates lengths from estimates
  lengths = lengths[method] #filter out lengths that aren't part of method
  ##calculate metrics for comparison
  #mean by Method and ID (ID is unique flight and Flight_Date combination)
  
  #calculates SE from length estimates
  if ((1 %in% method) & (length(method) == 1)){
    groupmeanuncor = vector(mode="list",1)
    groupmeanuncor[[1]] = uncormeans(lengths[[1]])
    groupmean=groupmeanuncor
  }else if ((1 %in% method)& (length(method) >1)){#uncorrected method and more than 1
    groupmeanuncor = vector(mode="list",1)
    groupmeanuncor[[1]] = uncormeans(lengths[[1]])
    groupmeancor= lapply(lengths[2:length(lengths)],function(x){
      cormeans(x)})
    groupmean=append(groupmeanuncor,groupmeancor)
  }else{
    groupmeancor= lapply(lengths,function(x){
      cormeans(x)})
    groupmean=groupmeancor
  }
  #RMSE by Method and ID
  grouprmse = lapply(c(1:length(lengths)),function(x){
    lengthdat = lengths[[x]]
    rmse = function(y){mean((y-valdat$Object_Length_mm[1]/1000)^2)^(1/2)};
    dat = aggregate(Cal.fit ~ ID,data=lengthdat,rmse)
    data.frame(unique(lengthdat[1:5]),dat)
  })
  
  #Sq Error by Method
  populationsqe = lapply(lengths,function(x) {
    dt = data.frame(x$Flight_Date,valdat$Adj_Baro_Alt,(x$Cal.fit-valdat$Object_Length_mm[1]/1000)^2);
    names(dt) = c("Flight_Date","Adj_Baro_Alt","Sq_Error");
    dt
  })
  
  #calculates CV from length estimates
  groupCV = lapply(lengths,function(x){
    cv = function(y){sd(y)/mean(y)} * 100;
    dt = aggregate(Cal.fit ~ ID,data=x,cv);
    dt = data.frame(unique(x[1:2]),dt)
    dt
  })
  
  popcv = data.frame(unlist(lapply(groupCV,function(x) mean(x$Cal.fit,na.rm=TRUE))))
  popmean = do.call("rbind.fill",lapply(groupmean,function(x) data.frame(mean(x$Cal.mean,na.rm=TRUE),mean(x$Cal.lwr,na.rm=TRUE),mean(x$Cal.upr,na.rm=TRUE))))
  poprmse = data.frame(unlist(lapply(grouprmse,function(x) mean(x$Cal.fit))))
  val_summ = data.frame(popmean,poprmse,popcv,popmean[1]-valdat$Object_Length_mm[1]/1000)
  names(val_summ) = c("Mean","Mean.lwr","Mean.upr","RMSE","CV","Bias")
  val_summ$Method = method
  val_summ = val_summ[,c(ncol(val_summ),1:(ncol(val_summ)-1))]
  
  return(list(val_summ,modlist))
}

#processing chain for morphometric measurement analysis
main_morph_proc = function(meas_dat,meas_dat.drop,valdat,modlist,method){
  # meas_dat = inprocdat[[1]]
  # meas_dat.drop = inprocdat[[2]]
  # valdat = inprocdat[[3]]
  # modlist = calprocdat[[2]]
  
  if (is.data.frame(valdat)){
    print("Analyzing Morphometric Measurements")
    idx = which(names(meas_dat[[1]])=="SA"):length(meas_dat[[1]])
    ##prediction to correct values
    cl<-makeCluster((detectCores()-2),timeout=600)
    registerDoParallel(cl)
    modlist2 = modlist[c(1,2,5,6)] #remove 3 and  because they don't apply in this context
    whalepreds = foreach(ii = 1:length(meas_dat)) %dopar% {
      source("Whale_Quant_Funcs.R")
      lapply(modlist2,function(x) predictions(x,meas_dat[[ii]],valdat))}
    ## calculation of measurements by observation
    morphs = foreach(ii=1:length(meas_dat)) %do% {
      source("Whale_Quant_Funcs.R")
      morphometrics(whalepreds[[ii]],data.frame(meas_dat[[ii]]),idx,ii,method)} #must convert meas_dat to dataframe for proper operation

    #sort by Method and flight as opposed to flight then Method
    morphs1 = as.list(data.frame(do.call(rbind, morphs)))# reorder the list
    morphs = lapply(c(1:length(morphs1)),function(x){
      dt = morphs1[[x]]
      dt2 = do.call("rbind",dt)
    })
    morphs = morphs[method] #only keep morphs matching desired method; rest are NA nayway
      
    #incorporate dropped data for method 1
    if ((1 %in% method) & (length(meas_dat.drop) > 0)){
      morphs = vector(0,mode="list")
      
      morphs.drop = foreach(ii=1:length(meas_dat.drop)) %do% {
        #source('Whale_Quant_Funcs.R')
        morphometrics(whalepreds[[ii]],meas_dat.drop[[ii]],idx,ii+length(unique(morphs$X1$ID)),c(1))}
      morphs1.drop = as.list(data.frame(do.call(rbind, morphs.drop)))# reorder the list
      morphs.drop = lapply(morphs1.drop,function(x) do.call("rbind",x))
      morphs.drop = morphs.drop[[1]]
      #combind morphs and morphs.dropped
      morphs[[1]] = rbind(morphs,morphs.drop)
    }
  #when present with all or some calibration data
  }else{
    print("Analyzing Morphometric Measurements Without CAlibration data")
    idx = which(names(meas_dat.drop[[1]])=="SA"):length(meas_dat.drop[[1]])
    ##prediction to correct values
    cl<-makeCluster((detectCores()-2),timeout=600)
    registerDoParallel(cl)
    morphs = vector(0,mode="list")
    morphs.drop = foreach(ii=1:length(meas_dat.drop)) %do% {
      #source('Whale_Quant_Funcs.R')
      morphometrics(whalepreds[[ii]],meas_dat.drop[[ii]],idx,ii,c(1))
      }
    
    morphs1.drop = as.list(data.frame(do.call(rbind, morphs.drop)))# reorder the list
    morphs.drop = lapply(morphs1.drop,function(x) do.call("rbind",x))
    morphs.drop = morphs.drop[[1]]
    morphs[[1]]=morphs.drop
  } #end exception loop for no calibration data
  
  #calc coefficient of variation
  morph_cv = do.call("rbind",lapply(morphs,function(x){
    cv = function(y){sd(y)/mean(y)} * 100;
    dt = aggregate(. ~ ID,data=x,cv,na.action=na.pass)[6:dim(x)[2]];
    idx2 = which(grepl('.fit',names(dt)))
    dt = dt[idx2] #extract .fit columns
    names(dt) = gsub(".fit", ".cv", names(dt))
    data.frame(unique(x[3:5]),dt)
  }))
  
  ## fix morph mean data for the null Method so that it has SE for CI
  if (1 %in% method & length(method) == 1){ #handles uncorrected method selection only
    uncormorphs = uncormeans(morphs[[1]]) #CI and mean for uncorrected Method
    morph_mean = uncormorphs
  }else if (1 %in% method & length(method) > 1){#handles uncorrected method combined with other methods
    uncormorphs = uncormeans(morphs[[1]]) #CI and mean for uncorrected Method
    morphs12 = foreach(ii=2:length(morphs)) %dopar% {
      source("Whale_Quant_Funcs.R")
      dt = morphs[[ii]]
      cormeans(dt)}
    cormorphs = do.call("rbind",morphs12)
    morph_mean = rbind(uncormorphs,cormorphs)
  }else{#other correction methods only
    morphs12 = foreach(ii=1:length(morphs)) %dopar% {
      source("Whale_Quant_Funcs.R")
      dt = morphs[[ii]]
      cormeans(dt)}
    cormorphs = do.call("rbind",morphs12)
    morph_mean = cormorphs
  }
  morph_summs = merge(data.table(morph_mean),data.table(morph_cv),by=c("ID","Method"),sort="TRUE")
  morph_summ = morph_summs[order(Method,ID)]
  stopCluster(cl)
  return(morph_summ)
}
