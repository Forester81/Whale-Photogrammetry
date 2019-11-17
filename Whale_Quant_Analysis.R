###########################################################
# Title: Whale Quantitative Analysis
# Objective: Creates summary tables of morphometric data measured in the complementary MATLAB programs
# Author: Jonathan Burnett
# Date Modified: November 19 2019
# Version 1.0
# Tested in R version 3.4.3 (2017-11-30) on Windows 10 x64 Platform
###########################################################

##############INSTRUCTIONS#################################
# SEE ATTACHED 'read.me' for detailed instructions
# make sure 'Whale_Quant_Funcs.R' is in the same folder as the 'Whale_Quant_Analysis
# Set configuration options under the 'Set Options' heading
# From R Studio window click 'source' to run the program
# Output is written to the parent folder containing the data
###########################################################
#Environmental Settings #
rm(list = lsf.str())
rm(list = ls.str(mode = 'numeric'))
library(rstudioapi)    
wd = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)
source('Whale_Quant_Funcs.R')

################################################################################
########################### Set Options########################################
################################################################################


############## Calibration Methods ##################
#refer to Burnett et al. 2018 for additional details
#method 1 = uncorrected
#method 5 = mixed effect correction with smoothing
#method 2 = basic linear correction
#method 3 = basic linear correction with smoothing
#method 4 = mixed effect correction
method = c(1) ###SET METHODS HERE; c(1:5) will run them all like the manuscript but takes a long time, c(1,5) runs method 1 and 5
############# end calibration methods input#########

############## FILTERING ###################
###turn on measurement data filter to remove extreme observations
filt = "no" #change to "yes" to activate filtering
filt_cal = "no" #change to "yes" to activate filtering of calibration data; filtering may remove observations that prevents scaling of data from entire flights
filt_intens = .01 #the outlier threshold e.g. .05 corresponds with a value outside the 95% confidence 
#############END FILTERING######################

############## Input Data Directory Selection ###################
#identify the directory containing all of the measurement data.  It's looking a folder containing .csv outputs from the two MATLAB programs.  

ingui = "yes" #change to 'no' to manually set input directory; input 'yes' to use interactive directory selection
if (ingui == "yes"){ #use no if you are a mac user because this doesn't work in mac
  print("Use GUI to select directory")
  dirs = tk_choose.dir(default = "", "Select directory containing whale data") #user interactive selection of data folder
}else{
  #for mac users be sure to use \\Users if data resides on internal HD or use \\Volumes if data resides on an external drive
  dirs = "G:\\hale_taylor\\pcformat\\" #this line is for debugging purposes
}
#############END DIRECTORY SELECTION ############################

#############END OPTIONS SECTION ################################


############# Main Processing Chain ###############

#run this section using either 'run' or 'source'

#version check
if ((as.numeric(R.version$major) < 3) | (as.numeric(R.version$minor) < 4.3)) {
  stop(sprintf("R version is too low: %s , requires version 3.4.3 or higher", paste(R.version$major,R.version$minor,sep=".")))
}

#data input and preprocessing 
inprocdat = inproc(dirs,filt_cal,filt,filt_intens)
write.table(rbindlist(inprocdat[[1]]), file = paste(dirs,"formatted_measurement_data.csv",sep="/"),row.names=FALSE,col.names=TRUE,sep=",")
#diagnostics
scaletable = inprocdat[[3]]
if (is.list(scaletable)){
  scaletable$Est_Len_mm = scaletable$GSD * scaletable$Object_Length_Pixels*1000
  scaletable$Perc_Diff_Len = abs((scaletable$Est_Len_mm - scaletable$Object_Length_mm)/scaletable$Object_Length_mm * 100)
  write.table(scaletable, file = paste(dirs,"formatted_scaling_data.csv",sep="/"),row.names=FALSE,col.names=TRUE,sep=",")
  #calibration data processing
  calprocdat = calproc(inprocdat[[3]],method)
}else{
  calprocdat=as.list(c(NA,NA))
  }#only run if there is calibrationdata

#measurement data processing
morph_summ = main_morph_proc(inprocdat[[1]],inprocdat[[2]],inprocdat[[3]],calprocdat[[2]],method)

#data writing
write.table(calprocdat[[1]], file = paste(dirs,"Validation.csv",sep="/"),row.names=FALSE,col.names=TRUE,sep=",")
cat(sprintf("calibration summary results written to %s\n",(paste(dirs,"Validation.csv",sep="/"))))

write.table(morph_summ, file = paste(dirs,"Morph_Measures.csv",sep="/"),row.names=FALSE,col.names=TRUE,sep=",")
cat(sprintf("Morphometric results written to %s\n",paste(dirs,"Morph_Measures.csv",sep="/")))

if (length(inprocdat[[2]])>0){#if there are dropped data
  write.table(do.call("rbind",inprocdat[[2]]), file = paste(dirs,"Non_Calibrated_Whales.csv",sep="/"),row.names=FALSE,col.names=TRUE,sep=",")
  cat(sprintf("%s whales were missing associated calibration data. Raw data for these whales were output to Non_Calibrated_Whales.csv, however, these whales were proccessed under calibration method 1 (no correction) and appear in the Morph_Measures.csv file written to the data folder prvided.\n", length(inprocdat[[2]])))
}

############ END MAIN PROGRAM #####################


