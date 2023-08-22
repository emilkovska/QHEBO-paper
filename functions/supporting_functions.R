#################################### MAIN FUNCTIONS   ####################################################

getModelTransitions <- function(sens_analysis = FALSE, what = NULL, change = NULL, inc_adjust = NULL) { #what: 1- mortality, 2-diab prevalence 3-diab incidence 4-RR 5-SMR
  
  bmi <- read.csv(paste(wd,"data/NCD_RisC_Lancet_smooth_2030_Bulgaria.csv",sep=""))
  names(bmi) <- c("X","year","sex","age","p_ow","p_ob","log.se.ob","log.se.ow","p_no")
  H.dimnames <- c("healthy","diabetic","dead")
  v_years <- unique(bmi$year)
  v_bmis <- c("no","ow","ob")
  v_sex <- c("Female","Male")
  bmi$X <- NULL
  
  df <- get_dfall_template(sens_analysis = sens_analysis, what = what, change = change, inc_adjust = inc_adjust)
  if (any(is.na(df)) != FALSE)
    stop("NaNs produced in df")
  
  df <- merge(df,bmi,by=c("year","age","sex"))
  
  df$sex <- as.factor(df$sex)
  df <- df[order(df$year,df$age),]
  df <- df[order(df$sex),]
  rownames(df) <- NULL 
  
  ############################ Populate the main data frame ##################################################
  P_array <- array(dim = c(3,3,81,15,3,2), dimnames = list(H.dimnames, 
                                                        H.dimnames,
                                                        20:100,
                                                        v_years, 
                                                        v_bmis,
                                                        v_sex))
  
  
  #Disaggregate mortality and diabetes prevalence into risk factor-specific ones
  df$mort_no <- df$mort/(df$p_no + 
                           df$p_ow*df$RR_m_ow +
                           df$p_ob*df$RR_m_ob)
  df$mort_ow <- df$mort_no*df$RR_m_ow
  df$mort_ob <- df$mort_no*df$RR_m_ob
  
  df$diab_no <- df$diab/(df$p_no + 
                           df$p_ow*df$RR_d_ow +
                           df$p_ob*df$RR_d_ob)
  df$diab_ow <- df$diab_no*df$RR_d_ow
  df$diab_ob <- df$diab_no*df$RR_d_ob
  
  #Disaggregate mortality into mort due to disease and not due to disease
  df$m_d_0 <- df$mort/(df$diab*df$smr + (1-df$diab))
  df$m_d_1 <- df$m_d_0*df$smr
  df$em <- df$m_d_1 - df$m_d_0        #m_d_1 - m_d_0 = em(YODO) 
  
  
  
  if (any(df$em < 0)) { #if em turns out negative, get the last positive one within the same year and use until age 100
    df$em <- repeat.neg(df$em)
  } 
  
  p_no_d_1 <- df$diab_no*df$p_no/df$diab
  p_ow_d_1 <- df$diab_ow*df$p_ow/df$diab
  p_ob_d_1 <- df$diab_ob*df$p_ob/df$diab
  numerator <- (1-df$diab)*df$em - ((p_no_d_1*(df$mort_no)+
                                       p_ow_d_1*(df$mort_ow)+
                                       p_ob_d_1*(df$mort_ob)) - df$mort)
  denominator <- p_no_d_1*(1-df$diab_no) + p_ow_d_1*(1-df$diab_ow) + p_ob_d_1*(1-df$diab_ob)
  df$am <- numerator/denominator
  
  if (any(df$am < 0)) { #if am turns out negative, get the last positive one within the same year and use until age 100
    df$am <- repeat.neg(df$am)
  } 
  
  
  df$d_inc_no <- df$d_inc_all/(df$p_no + df$p_ow*df$RR_d_ow + df$p_ob*df$RR_d_ob)
  df$d_inc_ow <- df$d_inc_no*df$RR_d_ow
  df$d_inc_ob <- df$d_inc_no*df$RR_d_ob

  
  #Used as transition probability to death state due to other causes 
  #for healthy individuals
  #Following equation (3) from YODO, we can calculate RF-specific other-cause mortality moc(R)
  # m(R) = moc(R)+p(d=1 | R) * am =>
  # moc(R) = m(R) - p(d=1 | R) * am 
  # moc(R) == m(d=0,R) - used as transition probability for both healthy and diseased to the death state (1,4) & (2,4) in the matrix
  df$m_d_0_no <- df$mort_no - df$diab_no*df$am
  df$m_d_0_ow <- df$mort_ow - df$diab_ow*df$am
  df$m_d_0_ob <- df$mort_ob - df$diab_ob*df$am
  
  #Then, following equation (1) from YODO where:
  #for diabetic individuals
  #m(d=0,R) = moc(R) & 
  #m(d=1,R) = m(d=0,R) + am
  #we can find RF-specific other-cause mortality for diseased & non-diseased. 
  #m(d=0,R) Used as transition probability from healthy --> death (1,4)
  #m(d=1,R) Used as transition probability from diabetic --> other-cause death (1,4) + death from diabetes (1,3) (AM)
  df$m_d_1_no <- df$m_d_0_no + df$am
  df$m_d_1_ow <- df$m_d_0_ow + df$am
  df$m_d_1_ob <- df$m_d_0_ob + df$am
  
  
  df <- df[order(df$year,df$age,df$sex),]
  df <- df[order(df$sex),]
  rownames(df) <- NULL
  
  
  #DF  for normal
  keepcol <- c("year","age","sex","p_no",
               "m_d_0_no", #rate to die from all-cause | BMI, healthy
               "m_d_1_no", #rate to die from all-cause | BMI, diabetic
               "am",       #good to have
               "d_inc_no") #rate to become diabetic
  df_no <- df[,keepcol]
  
  
  #DF  for overweight
  keepcol <- c("year","age","sex","p_ow",
               "m_d_0_ow", #rate to die from all-cause | BMI, healthy
               "m_d_1_ow", #rate to die from all-cause | BMI, diabetic
               "am",       #good to have
               "d_inc_ow") #rate to become diabetic
  df_ow <- df[,keepcol]
  
  #DF  for obese
  keepcol <- c("year","age","sex","p_ob",
               "m_d_0_ob", #rate to die from all-cause | BMI, healthy
               "m_d_1_ob", #rate to die from all-cause | BMI, diabetic
               "am",       #good to have
               "d_inc_ob") #rate to become diabetic
  df_ob <- df[,keepcol]
  
  names(df_no) <- names(df_ow) <- names(df_ob) <- c("year","age","sex","p_BMI","m_d_0","m_d_1","am","d_inc")
  
  df <- list(df_no,df_ow,df_ob)
  names(df) <- c("no","ow","ob")
  
  if (sens_analysis != TRUE) {
    saveRDS(df, file= paste(wd,"BMIdfs.RData",sep=""))
  } else {
    saveRDS(df, file= paste(wd,"BMIdfs_sens.RData",sep=""))
  }
  
  ###########################       Populate the transition matrices     ########################### 
 
  
  for (i in 1:2) {                                    #For each sex; 1- Women, 2- Men
    for (ii in 1:3) {                                 #For each BMI state
      for (iii in 1:length(v_years)) {                #For each year
        for (iiii in 20:100) {                        #For each age
          
          sub <- df[[ii]] %>%
            filter(year == v_years[iii],
                   age == iiii,
                   sex == v_sex[i])
          
          p_all_d_0 <- sub$m_d_0      #rate to die from all-cause | BMI
          p_all_d_1 <- sub$m_d_1
          p_inc <- sub$d_inc          #probability to develop diabetes | BMI 
          
          ###########################       Rates to Prob Normal     ########################### 
          
          probs <- as.numeric(unlist(c(p_all_d_0,p_all_d_1)))
          if(any(probs<0)) 
            stop(paste("NaNs produced in prob_t=",t,sep=""))
          probs <- rate_to_prob(probs,t=1)
          p_all_d_0 <- probs[1]
          p_all_d_1 <- probs[2]
          
          #compile matrix
          P <- matrix(0,3,3)
          rownames(P) <- colnames(P) <- H.dimnames
          
          #From Healthy
          P["healthy","healthy"] <- (1-p_all_d_0)*(1-p_inc) #Conditionally on not dying, then remain healthy 
          P["healthy","diabetic"] <- (1-p_all_d_0)*p_inc    #Conditionally on not dying, then p(inc)
          P["healthy","dead"] <- p_all_d_0
          
          #From Sick
          P["diabetic","diabetic"] <- (1-p_all_d_1) # AM and moc(R) are mutually exclusive because we care about the probability of
          P["diabetic","dead"] <- p_all_d_1             # death either from AM or other causes. P(A or B) = P(A) + P(B). Either can occur
                                                    #but the prob to remain alive and diabetic is contingent on neither of the deaths occurring at once
                                                         
          #From Dead
          P["dead","dead"] <- 1
         
          
          #Record transition probs
          P_array[,,paste(iiii),
                  paste(v_years[iii]),
                  v_bmis[ii],
                  v_sex[i]] <- P
        }
      }
    }
  }
  
  if (sens_analysis != TRUE) {
    saveRDS(P_array, file= paste(wd,"P_array.RData",sep="")) 
  } else {
    saveRDS(P_array, file= paste(wd,"P_array_sens.RData",sep=""))
  }
  
 
  
  return(P_array)
}

runCohortBMI <- function(sens_analysis,start_age,input_year,prev_change = NULL,risk_change = NULL) {
  
  
  if (!isTRUE(sens_analysis)) {
    BMIdf <- readRDS(file=paste(wd,"BMIdfs.RData",sep=""))
    P_array <- readRDS(file= paste(wd,"P_array.RData",sep="")) #dim=c(4,4,age,year,BMI,sex)
  } else {
    BMIdf <- readRDS(file=paste(wd,"BMIdfs_sens.RData",sep=""))
    P_array <- readRDS(file= paste(wd,"P_array_sens.RData",sep="")) #dim=c(4,4,age,year,BMI,sex)
  }
  
  v_age <- seq(from = start_age, to = 100)
  v_sex <- c("Female","Male")
  M_f <- readRDS(file=paste(wd,"P1matrix_Female_",input_year,".RData",sep=""))
  M_m <- readRDS(file=paste(wd,"P1matrix_Male_",input_year,".RData",sep=""))
  n_t <- length(v_age)
  v_bmi <- c("no","ow","ob")
  H.dimnames <- c("healthy","diabetic","dead")
  cohort <- array(dim = c(n_t,4,3,2),dimnames = list(v_age,
                                                     c("healthy","diabetic","dead","incidence"),
                                                     v_bmi,
                                                     v_sex))
  
  
  #Get the baseline RF vector - first column is Women, second - Men
  m_RF <- rbind(BMIdf[["no"]]$p_BMI[BMIdf[["no"]]$year == input_year & BMIdf[["no"]]$age == start_age],
                BMIdf[["ow"]]$p_BMI[BMIdf[["ow"]]$year == input_year & BMIdf[["ow"]]$age == start_age],
                BMIdf[["ob"]]$p_BMI[BMIdf[["ob"]]$year == input_year & BMIdf[["ob"]]$age == start_age])
  
  
  
  #Apply childhood intervention if any
  if (!is.null(prev_change) & prev_change != 1) {
    m_RF[2:3,] <- m_RF[2:3,] * prev_change
    m_RF[1,] <- 1 - colSums(m_RF[2:3,])
  }
  
  
  #Assign cohort start vector 
  cohort[paste(start_age),1:3,,] <- matrix(0,6,3)
  cohort[paste(start_age),"healthy",,] <- m_RF
  
  #P_array: dim=c(3,3,age,year,BMI,sex)
  
  for (s in 1:2) {
    for (t in (start_age+1):100) {
      for (b in 1:3) {
        
        cohort[paste(t),1:3,v_bmi[b],v_sex[s]] <- cohort[paste(t-1),1:3,
                                                         v_bmi[b],v_sex[s]]%*%P_array[,,paste(t-1),
                                                                                             paste(input_year),v_bmi[b],v_sex[s]]
        
        cohort[paste(t),"incidence",v_bmi[b],v_sex[s]] <- cohort[paste(t-1),"healthy",
                                                                 v_bmi[b],v_sex[s]] * P_array["healthy","diabetic",
                                                                                                           paste(t-1),
                                                                                                           paste(input_year),v_bmi[b],v_sex[s]]
      }
      
      #Get BMI Matrix
      if (s == 1) {
        M <- M_f
      } else {
        M <- M_m
      }
      
      #(t-1)-19
      
      M1 <- M[["E.trans.prob"]][,,paste("age",t-1,sep="_")]
      M1 <- t(M1) 
      
      #Apply adulthood intervention if any
      if (!is.null(risk_change) & risk_change != 1 & (20 <= (t-1) & (t-1) < 65)) {
       M1["no","ob"] <- M1["no","ob"] * risk_change
       M1["no","ow"] <- M1["no","ow"] * risk_change
       M1["ow","ob"] <- M1["ow","ob"] * risk_change
       M1[1,1]       <- 1 - (M1[1,2] + M1[1,3])
       M1[2,2]       <- 1 - (M1[2,1] + M1[2,3])
      }
      
      #Adjust current vector - healthy
      cohort[paste(t),"healthy",,v_sex[s]]  <- cohort[paste(t),"healthy",,v_sex[s]]%*%M1

      #Adjust current vector - diabetic
      cohort[paste(t),"diabetic",,v_sex[s]]  <- cohort[paste(t),"diabetic",,v_sex[s]]%*%M1
      
    }
  }

  return(cohort)
} 

#################################### SUPPORTING FUNCTIONS   ####################################################

#01 # Create function to clear data ############################################################################
#Used in function break_down

repeat.before <- function(x) {
  ind <- which(!is.na(x))
  ind_rep <- ind
  if (is.na(x[1])) {
    ind_rep <- c(min(ind), ind)
    ind <- c(1, ind)
  }
  rep(x[ind_rep], times = diff(c(ind, length(x) + 1)))
}

#02 # Create function to substitute negative AMs ###############################################################
#Used in function calc_em_am

repeat.neg <- function(x) {
  neg <- which(x < 0)
  l_neg <- length(neg)
  for (i in 1:l_neg) {
    x[neg[i]] <- ifelse(neg[i]==1,
                        x[min(which(x>0))],
                        x[neg[i]-1])
  }
  return(x)
}



#03 # Create function to covert rates to probabilities  ################################################
#Previously this function was built in package "darthtools" but it seems the package does not work in the newest version of R

rate_to_prob <- function(r, t = 1){
  if ((sum(r < 0) > 0)) {
    stop("rate not greater than or equal to 0")
  }
  p <- 1 - exp(-r * t)
  return(p)
}

#04 # Create function equivalent to RIGHT and LEFT in Excel  ################################################
right <- function (string, char) {
  substr(string,nchar(string)-(char-1),nchar(string))
}

left <- function (string,char) {
  substr(string,1,char)
}

#05 # Create function create the template df_all  ################################################
#Used in get_dfall()


get_dfall_template <- function(sens_analysis,what, change,inc_adjust) { 
  
  df <- data.frame(year = rep(rep(2000:2050, length(20:100)),2),
                   age = rep(rep(20:100, each = length(2000:2050)),2),
                   sex = rep(c("Female","Male"),each = 4131))
  
  which_col <- list("mort","diab","d_inc_all",
                    c("RR_m_ow","RR_m_ob","RR_d_ow","RR_d_ob","RR_m_ow_diab","RR_m_ob_diab",
                      "RR_m_ow_nondiab","RR_m_ob_nondiab"),
                    "smr")
  
  diab <- read.csv(paste(wd,"data/pred_prev_diab.csv",sep=""))            #Import diabetes prevalence projections to 2050
  smr <- read_excel(paste(wd,"data/source/smr.xlsx",sep=""))              #Import SMR
  mort <- read.csv(paste(wd,"data/pred_allcause_mort_LC_GAM.csv",sep="")) #Import mortality estimates
 
  
  
  #Generate RR df
  RR <- all_RRs(20)
  RR <- RR %>% pivot_longer(cols = m_ow_women:d_ob_men,
                            names_to = "TypeRR",
                            values_to = "RR")
  RR$sex <- RR$type <- 0
  RR$sex[grep("_men",RR$TypeRR)] <- "Male"
  RR$sex[grep("_women",RR$TypeRR)] <- "Female"
  RR$type[grep("m_ow_",RR$TypeRR)] <- "RR_m_ow"
  RR$type[grep("m_ob_",RR$TypeRR)] <- "RR_m_ob"
  RR$type[grep("d_ow_",RR$TypeRR)] <- "RR_d_ow"
  RR$type[grep("d_ob_",RR$TypeRR)] <- "RR_d_ob"
  RR <- RR[,-2]
  RR <- pivot_wider(RR,names_from=type,values_from = RR)
  
  
  
  #Generate rr df
  rr <- read_excel(paste(wd,"data/source/rr.xlsx",sep=""),sheet = "RR")
  rr$sex[rr$sex==2] <- "Female"
  rr$sex[rr$sex==1] <- "Male"
  rr <- pivot_wider(rr,names_from=c(bmi_gr,diab),values_from = RR_det)
  names(rr) <- c("sex","RR_m_ow_diab","RR_m_ob_diab","RR_m_ow_nondiab","RR_m_ob_nondiab")

  
  
  #Generate smr df
  names(smr) <- c("age","Female","Male")
  smr <- smr %>% pivot_longer(cols = Female:Male,
                              names_to = "sex",
                              values_to = "smr")
  
  
  #Generate mort_diab df
  mort_diab <- merge(mort,diab, by=c("year","age","sex"))
  drops <- c("X.x","X.y")
  mort_diab <- mort_diab[ ,!(names(mort_diab) %in% drops)]
  names(mort_diab) <- c("year","age","sex","mort","se.mort","diab",
                        "se.diab","lower.d","upper.d")
  
  
  #Generate d_inc df
  d_inc <- read.csv(paste(wd,"data/source/diab_inc_smooth.csv",sep=""))
  d_inc$X <- NULL
  names(d_inc) <- c("age","Male","Female")
  
  
  if (!is.null(inc_adjust)) {
    d_inc$Male   <- d_inc$Male * inc_adjust[1]
    d_inc$Female <- d_inc$Female * inc_adjust[2]
  }
  
  d_inc <- d_inc %>% pivot_longer(cols = Male:Female,
                                  names_to = "sex",
                                  values_to = "d_inc_all")
  
  
  df <- merge(df,mort_diab,by= c("year","age","sex"))
  df <- merge(df,RR, by=c("age","sex"))
  df <- merge(df,rr,by="sex")
  df <- merge(df,smr,by=c("age","sex"))
  df <- merge(df,d_inc,by.x=c("age","sex"))
  
  
  df$sex <- as.factor(df$sex)
  df <- df[order(df$year,df$age),]
  df <- df[order(df$sex),]
  rownames(df) <- NULL 
  
  df$RR_m_ow_nondiab <- df$RR_m_ow
  df$RR_m_ob_nondiab <- df$RR_m_ob
  
  
  
  if (sens_analysis == TRUE) {
        df[,which_col[[what]]] <- df[,which_col[[what]]] * change
        write.table(df,file = paste(wd,"data/source/df_all_sens_template.txt",sep=""))
  } else {
        write.table(df,file = paste(wd,"data/source/df_all_template.txt",sep="")) 
  }

  
  return(df)
}

#06 # Create function to get all the RRs #######################################################################
#Used in function break_down     

all_RRs <- function(start_age) {
  rr_default <- c(1.15,1.5,1.2,1.55) #Default RR for ow women, ob women, ow men, ob men. Same as in file rr_bmi_m_allcause_DYNAMO-HIA_BG.xlsx
  v_adjust <- c(0.98,0.95,0.90)      #Adjustment coefficient per DYNAMO-HIA
  v_asofage <- c(50,60,70)           #Apply adjustment coefficient as of age vector
  RR <- data.frame(age <- seq(start_age,to=100),
                   matrix(NA,101-start_age,
                          length(rr_default)*2))    #create blank data frame to store final RR
  names(RR) <- c("age","m_ow_women","m_ob_women","m_ow_men","m_ob_men", #all-cause mortality RRs
                 "d_ow_women","d_ob_women","d_ow_men","d_ob_men") #diabetes prevalence RRs
  
  
  #get All-cause mortality RRs for population aged start_age +
  for (i in 1:length(rr_default)) {
    RR[,i+1] <- get_RR(start_age,rr_default = rr_default[i] ,v_adjust,v_asofage)[,2]
  }
  
  #get Diabetes RRs for population aged start_age +
  rr_default <- c(2.3,7,2.25,5.5) #ow women, ob women, ow men, ob men
  v_adjust <- c(0.92,0.90)
  v_asofage <- c(60,75)
  
  for (i in 1:length(rr_default)) {
    RR[,i+5] <- get_RR(start_age,rr_default = rr_default[i],v_adjust,v_asofage)[,2]
  }
  
  return(RR)
}

#07 #   Create function to generate RRs ########################################################################
#Used in function all_RRs

get_RR <- function(start_age,rr_default,v_adjust,v_asofage) { #v_adjust and v_asofage need to be the same length
  n <- 101 - start_age
  data <- data.frame(Age = seq(start_age,length.out = n, by=1))
  x <- length(v_asofage)
  data$RR <- rep(rr_default,n)
  if (x > 0) {
    for (i in 1:x) {
      start_position <- v_asofage[i] - start_age + 1
      data$RR[start_position:n] <- 1 + v_adjust[i]*(rr_default-1)
    }
  }
  return(data)
}

