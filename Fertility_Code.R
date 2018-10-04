# Calculate asfr and tfr for MICS surveys
##Using both birth histories (modified from spss codes from Patrick) and last child birth questions

rm(list = ls())

library(tidyr)
library(dplyr) 
library(tools)

#List of MICS Women(WM), Birth History(BH) and Houselisting(HL data)
Path_WM <- "V:/FertilitySection/Alkema_Joint project on contraceptive use trends/1_All-women-estimates&projections/Data-tabulations/MICS/Translated_RDataFiles"
WM.data <- dir(Path_WM)

Path_BH <- "V:/MICS/Documentation/scripts_Natalie/BirthHistory_microdata/Tramslated_BHData"
BH.data <- dir(Path_BH)

Path_HL <- "V:/FertilitySection/Data-sets_Marriage/World Marriage Data 2017/MICS microdata/MICS_HH_Translate/Translated_RData"
HL.data <- dir(Path_HL)

MICS.Master<-read.csv("V:/MICS/MICSMaster.csv")

#Period - Define the numbers of years in a period
# 5 years = 60
# 3 years = 36
period <- 60

#maxperiod - Define upper limit for periods
# 5 years and 3 years = 5
maxperiod <- 5
  
#NoBH_method - Choose method to calculate no birth history 
#Variable = using last birth variable - set as default (decided this ensure population of numerator and denominator match)
#Calculate = using calculation to extract # of children
NoBH_method <- "Variable"

#Output File
filename <- paste0("V:/MICS/Documentation/scripts_Natalie/Fertility_MICS_",(period/12),"years_",NoBH_method,".csv")
file.create(filename) #Overwrite existing output

#Only for bh file
FertB <- function(df){
  ##Get Limit based on desired observation years
  df$UpperLim <- df$INTDATECMC-1
  df$LowerLim <- df$INTDATECMC-(maxperiod*period) 
  
  df$AgeAtBirth <- df$KIDDOBCMC - df$DOBCMC #Children's birth - women's birth = age of woman when gave birth
  df$age5 <- trunc(df$AgeAtBirth/60) - 2 #Fit in 5 years age group
  df$age5_LAB <- factor(df$age5,
                        levels = c(0,1,2,3,4,5,6,7),
                        labels = c("<15","15-19","20-24","25-29","30-34","35-39","40-44","45-49"))
  
  df$colper <- trunc((df$UpperLim - df$KIDDOBCMC)/period) #Gets number of group depending on the period desired for analysis 
  
  df.3years <- df[which(df$KIDDOBCMC >= df$LowerLim & df$UpperLim>= df$KIDDOBCMC & df$age5 >0),] #Restrict dataset to be within analysis period
  
  #Aggregate for the number of births by woman, age and period
  BirthCount <- group_by(df.3years, ID,age5,colper)%>%
    summarise(BirthCount = n())
  
  BirthCount <- merge(df.3years[,c("ID","WMWEIGHT")],BirthCount,by="ID")
  BirthCount <- BirthCount[!duplicated(BirthCount),]
  
  #Summarize by age, number of birth, and period
  BirthNum <- as.data.frame(xtabs(WMWEIGHT~age5+BirthCount+colper,BirthCount))
  
  #Weighted number of birth = number of birth (by aggregation) * Frequency of such number of birth from xtabs
  BirthNum$WeightedNum <- as.numeric(BirthNum$BirthCount) * as.numeric(BirthNum$Freq)
  
  #Aggregate number of birth by age group and period
  BIRTHS_TABLE <- aggregate(BirthNum$WeightedNum, by = list(age5=BirthNum$age5,colper = BirthNum$colper), FUN = sum)
  BIRTHS_TABLE$Number_Births <- BIRTHS_TABLE$x
  BIRTHS_TABLE<-BIRTHS_TABLE[,c("age5","colper","Number_Births")]
  
  #Set indicator variable if 45-49 is in age group
  if(7%in%df$age5){
    IndicatorAllo <- T  
  }
  
  assign("IndicatorAllo",IndicatorAllo,.GlobalEnv)
  
  return (BIRTHS_TABLE)
}

FertE <- function(df){
  #c - loop through 4 groups of years preceding the survey
  for(c in 0:4){
    
    #Get Time Limit based on desired observation years
    df$UpperLim <- df$INTDATECMC-(period*c)-1
    #Get age group 
    df$age <- df$UpperLim - df$DOBCMC
    df$age5 <- trunc(df$age/60)
    
    #Get Exposure in months for current age group and previous age group within analysis period
    #Exp_CurAge = Exposure in current age group
    #Exp_PreAge = Exposure in previous age group
    df$Exp_CurAge <- df$age - (df$age5*60)+1
    df$Exp_CurAge<- ifelse(df$Exp_CurAge>period,period,df$Exp_CurAge) #Bound Exposure of current age within the months of analysis period
    
    df$Exp_PreAge <- ifelse(df$Exp_CurAge<period,(period-df$Exp_CurAge),0)
    
    
    #Recast age to group from 1-7
    df$age5 <- df$age5-2
    df.CurAge <- subset(df[,c("age","age5","Exp_CurAge","wmweight")], subset = df$age5>0)
    
    #Getting sum of exposure (by high and low refers to age group)
    HIGH.EXP<-as.data.frame(xtabs((wmweight*Exp_CurAge)~age5,df.CurAge))
    HIGH.EXP$colper <- c
    
    #Reduce age group by 1 (for possible of respondent's age 3 years ago)
    df$age5 <- df$age5-1
    df.PreAge <- subset(df[,c("age","age5","Exp_PreAge","wmweight")],subset=df$age5>0&df$Exp_PreAge>0)
    
    LOW.EXP <- as.data.frame(xtabs((wmweight*Exp_PreAge)~age5, df.PreAge))
    LOW.EXP$colper <- c
    
    #Merge all the time period to 1 table
    if(c==0){
      tab.exp <- rbind(HIGH.EXP,LOW.EXP)
    }else{
      exp <- rbind(HIGH.EXP, LOW.EXP)
      tab.exp <- rbind(tab.exp,exp)
    }
  }
  
  tab.exp$Freq <- tab.exp$Freq / 12
  EXP <- aggregate(tab.exp$Freq, by = list(age5=tab.exp$age,colper=tab.exp$colper),FUN=sum)
  EXP$Exposure <- EXP$x
  EXP<-EXP[,c("age5","colper","Exposure")]
  
  return(EXP)
}

FertP <- function(df){
  #Get total number of woman per age group from household file 
  df$age5 <- ifelse((df$AGE5YEAR>=15 & df$AGE5YEAR<=19),1,
                    ifelse((df$AGE5YEAR>=20 & df$AGE5YEAR<=24),2,
                           ifelse((df$AGE5YEAR>=25 & df$AGE5YEAR<=29),3,
                                  ifelse((df$AGE5YEAR>=30 & df$AGE5YEAR<=34),4,
                                         ifelse((df$AGE5YEAR>=35 & df$AGE5YEAR<=39),5,
                                                ifelse((df$AGE5YEAR>=40 & df$AGE5YEAR <=44),6,
                                                       ifelse((df$AGE5YEAR>=45 & df$AGE5YEAR <=49),7,100)))))))
  
  tab.numW <- as.data.frame(xtabs(hhweight~age5+SEX,df))
  tab.numW$Population <- sum(tab.numW$Freq)
  
  tab.numW$nwm <- tab.numW$Freq
  #Subset to woman and age bewtween 15-49
  NWM.Table <- tab.numW[which(tab.numW$SEX==2&tab.numW$age5!=100),c("age5","nwm","Population")]
  NWM.Table$colper <- 0 #Force this to only display once in output
  
  return(NWM.Table)
}

Fert_WithLastBirth <- function (df){
  var_indicator <- NA
  #Addressing difference in Last Birth variable
  if("BIRTHLAST2YRS"%in% names(df)){
    var_indicator <- "2years"
    #Number of woman that has birth in the last 2years (By age)
    tab_Births<-as.data.frame(xtabs(wmweight~AGE5YEAR+BIRTHLAST2YRS,df))
    tab_Births <- tab_Births[which(tab_Births$BIRTHLAST2YRS==1),] 
  }else if("BIRTHLASTYR" %in% names(df)){
    var_indicator <- "1year"
    #Number of woman that has birth in the last year (By age)
    tab_Births<-as.data.frame(xtabs(wmweight~AGE5YEAR+BIRTHLASTYR,df))
    tab_Births <- tab_Births[which(tab_Births$BIRTHLASTYR==1),]
  }else{
    return (NULL)
  }
  
  tab_Births$Number_Births <- tab_Births$Freq
  
  #Exposed woman = number of woman in wm survey per age group
  tab_NumW<-as.data.frame(xtabs(wmweight~AGE5YEAR,df))
  tab_NumW$Exposure <- tab_NumW$Freq
  
  Fert.Table <- merge(tab_Births[,c("AGE5YEAR","Number_Births")], tab_NumW[,c("AGE5YEAR","Exposure")], by = "AGE5YEAR",all=T)
  Fert.Table$asfr <- (Fert.Table$Number_Births/Fert.Table$Exposure) * 1000
  #Regroup age to numeric (for aggregation later)
  Fert.Table$age5 <- ifelse(Fert.Table$AGE5YEAR == 113, 0,
                            ifelse(Fert.Table$AGE5YEAR == 120, 1,
                                   ifelse(Fert.Table$AGE5YEAR == 130,2,
                                          ifelse(Fert.Table$AGE5YEAR == 140,3,
                                                 ifelse(Fert.Table$AGE5YEAR == 150,4,
                                                        ifelse(Fert.Table$AGE5YEAR == 160, 5,
                                                               ifelse(Fert.Table$AGE5YEAR == 170, 6,
                                                                      ifelse(Fert.Table$AGE5YEAR == 180,7,NA))))))))
  Fert.Table <- Fert.Table[,-1]
  assign("var_indicator",var_indicator,.GlobalEnv)
  return(Fert.Table)
}


Output <- function(Tab){
  if (file.size(filename) == 0){
    # if the csv output file is empty append the computed values to the output file but output the column names first, that is,in the first row
    write.csv(Tab, file = filename, quote = TRUE, append = TRUE, sep = ",", row.names = FALSE, col.names = TRUE)
  } else {
    # if the csv output file already has observations in it append the results to the output file without displaying column names each time data is outputted
    write.table(Tab, file = filename, quote = TRUE, append = TRUE, sep = ",", row.names = FALSE, col.names = FALSE)
  }
}

for (i in 1:length(WM.data)){
  #Get surveycode for searching through corresponding bh and hl data
  hl.name <- paste0(substr(WM.data[i],1,2),"hl",substr(WM.data[i],5,5),".RData")
  bh.name <- paste0(substr(WM.data[i],1,2),"bh",substr(WM.data[i],5,5),".RData")
  
  #Survey Information from Master
  SurveyInfo <- MICS.Master[which(MICS.Master$Individual.Recode == file_path_sans_ext(WM.data[i])),]
  
  if(bh.name %in% BH.data){
    #When Birth History data exists
    ####BH FILE TO CALCULATE NUMBER OF BIRTHS####
    load(paste(Path_BH,bh.name,sep="/"))
    
    #Missing cluster number - cannot create unique ID for merge
    if(!"CLUSTERNO"%in%names(df)){
        print(paste0("NO ID: ",bh.name))
      next()
    }
    print(bh.name)

    df$CLUSTERID <- sprintf(paste0("%0",max(nchar(df$CLUSTERNO)),"d"),as.numeric(df$CLUSTERNO))
    df$HHID <- sprintf(paste0("%0",max(nchar(df$HHNUM)),"d"),as.numeric(df$HHNUM))
    df$WMID <- sprintf(paste0("%0",max(nchar(df$LINENO)),"d"),as.numeric(df$LINENO))

    #Unique ID for woman
    df$ID <- paste0(df$CLUSTERID,df$HHID,df$WMID)
    
    #Set CMC Date to Numeric
    df$KIDDOBCMC <- as.numeric(df$KIDDOBCMC)
    df$INTDATECMC <- as.numeric(df$INTDATECMC)
    df$DOBCMC <- as.numeric(df$DOBCMC)
    df$WMWEIGHT <- as.numeric(df$WMWEIGHT)
    
    BIRTHS_TABLE<-FertB(df) #Return number of births
    ##############################################
    rm(df)
    ####WM FILE TO CALCULATE EXPOSURES####
    load(paste(Path_WM,WM.data[i],sep="/"))
    
    #Set CMC Date to Numeric
    df$INTDATECMC <- as.numeric(df$INTDATECMC)
    df$DOBCMC <- as.numeric(df$DOBCMC)
    df$wmweight <- as.numeric(df$wmweight)
    
    EXP <- FertE(df) #Return number of exposed women
    ##############################################
    rm(df)
    ####HL FILE TO CALCULATE NUMBER OF WOMEN and POPULATION####
    load(paste(Path_HL,hl.name,sep="/"))
    
    #Set numeric variables to correct datatype
    df$hhweight <- as.numeric(df$hhweight)
    if(hl.name == "sahl4.RData"){
      df$AGE5YEAR <- df$WA
    }else{
      df$AGE5YEAR <- as.numeric(df$AGE5YEAR)
    }
    
    NWM.Table <- FertP(df) #Return number of women in population
    #############################################
    Fert.Table <- merge(BIRTHS_TABLE, EXP, by = c("age5","colper"),all=T)
    Fert.Table <- merge(Fert.Table, NWM.Table, by=c("age5","colper"),all=T)
    Fert.Table$asfr <- Fert.Table$Number_Births / Fert.Table$Exposure *1000
    Fert.Table$methods <- "Birth Histories"
    
    Fert.Table <- Fert.Table[!is.na(Fert.Table$colper),]
   
  }else{
    #Cases when birth history file does not exists
    #Load woman file
    load(paste(Path_WM,WM.data[i],sep="/"))
    
    print(WM.data[i])
    
    var.List <- c("BIRTHLAST2YRS" , "BIRTHLASTYR")
    
    #Must have recent birth questions to be able to calculate
    if(all(!var.List %in% names(df))){
       next()
    }
    
    Fert.Table <- Fert_WithLastBirth (df) #Return df with number of births, number of woman (both in population and wm survey) and asfr
    
    IndicatorAllo <- F
    
    #If number of births = NA, set to 0 - suggesting no women in such age group had a birth recently
    Fert.Table$Number_Births[which(!is.na(Fert.Table$Exposure) & Fert.Table$Exposure>0 & is.na(Fert.Table$Number_Births))] <- 0

    #Calculate for asfr
    if("BIRTHLAST2YRS"%in% names(df)){
      Fert.Table$asfr <- 0.5*(Fert.Table$Number_Births/Fert.Table$Exposure) * 1000
    }else{
      Fert.Table$asfr <- (Fert.Table$Number_Births/Fert.Table$Exposure) * 1000
    }
    ####HL FILE TO CALCULATE NUMBER OF WOMEN and POPULATION####
    load(paste(Path_HL,hl.name,sep="/"))
    
    #Set numeric variables to correct datatype
    df$hhweight <- as.numeric(df$hhweight)
    NWM.Table <- FertP(df)
    #############################################
    Fert.Table <- merge(Fert.Table, NWM.Table, by=c("age5"),all=T)
    
    Fert.Table$methods <- "Recent Births"
    Fert.Table <- Fert.Table[!is.na(Fert.Table$asfr),]
  }
  
  #For both recent birth and birth history
  ##Aggregate 
  cumSum <- Fert.Table %>%
    group_by(colper,age5)%>%
    summarise(value = sum(asfr),
              births = sum(Number_Births),
              exp = sum(Exposure))%>%
    mutate(sumasfr=cumsum(value),
           cum_births = cumsum(births),
           cum_exp = cumsum(exp))

  Fert.Table <- merge(Fert.Table, cumSum[,c("colper","age5","sumasfr","cum_births","cum_exp")], by = c("colper","age5"),all.x=T)

  Fert.Table <- Fert.Table[order(Fert.Table$colper,Fert.Table$age5),]
  
  #Calculate asfr and tfr when specific age groups are missing - Imputing based on previous period and age groups
  if(IndicatorAllo == T){
    Allocate.list <- which(is.na(Fert.Table$sumasfr))
    for(r in Allocate.list){
      #Assumed the missing age group has an asfr of the difference between previous year breakdown and the 5years younger of the 
      #current year break down divided by the 5years younger group of the previous year breakdown. 
      Fert.Table[r,"sumasfr"] <- (Fert.Table[r-7,"sumasfr"]*Fert.Table[r-1,"sumasfr"])/Fert.Table[r-8,"sumasfr"]
      Fert.Table[r,"asfr"] <- Fert.Table[r,"sumasfr"] - Fert.Table[r-1,"sumasfr"]
    }
  }else{
    Fert.Table <- Fert.Table[!is.na(Fert.Table$asfr),]
  }
  
  Fert.Table$tfr <- NA
  Fert.Table$tfr[which(Fert.Table$age5 == 7)] <- (Fert.Table$sumasfr[which(Fert.Table$age5==7)]*5)/1000
  
  Fert.Table$age5 <- factor(Fert.Table$age5,
                            levels = c(0,1,2,3,4,5,6,7),
                            labels = c("[<15]","[15-19]","[20-24]","[25-29]","[30-34]","[35-39]","[40-44]","[45-49]"))
  
  Fert.Table$cbr <- NA
  
  Fert.Table <- Fert.Table%>%
    mutate(Country.Name = SurveyInfo$CountryName.UN,
           ISO.Code = SurveyInfo$LocID,
           Catalog.ID = SurveyInfo$CatalogID,
           Start.Year = SurveyInfo$StartYear,
           End.Year = SurveyInfo$EndYear,
           Survey.ShortName = SurveyInfo$ShortName,
           Survey.LongName = SurveyInfo$SurveyName,
           FieldTime.Mid = as.numeric(as.character(SurveyInfo$FieldWorkMiddle))
    )
  
  #Calculate midpoint of time
  Fert.Table$TimeMid <- NA
  #For recent birth
  if(NoBH_method == "Variable"){
    if(var_indicator == "1year"){
      Fert.Table$TimeMid <- Fert.Table$FieldTime.Mid - 0.5
    }else{
      Fert.Table$TimeMid <- Fert.Table$FieldTime.Mid - 1
    }
  }else{
    Fert.Table$TimeMid[which(Fert.Table$colper==0)] <- Fert.Table$FieldTime.Mid - 0.5*(period/12)
    Fert.Table$TimeMid[which(Fert.Table$colper!=0)] <- Fert.Table$FieldTime.Mid - (0.5*(period/12)+(as.numeric(Fert.Table$colper)*(period/12)))
  }
  
  Fert.Table <- select(Fert.Table, Country.Name:End.Year,TimeMid, colper:age5, Number_Births:Exposure,nwm:Population,asfr,sumasfr:tfr, Survey.ShortName:Survey.LongName, methods)
  
  Output(Fert.Table) 
  rm(Fert.Table, cumSum, df)
}
