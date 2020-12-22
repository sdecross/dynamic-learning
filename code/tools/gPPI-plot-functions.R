#######################################################
# FUNCTIONS FOR DATA WRANGLING FOR MT FEAR PPI PLOTS
#######################################################


# MANIPULATE ROI FILES ----------------------------------------------------------------------------
# Usage: roi1_ctl <- lapply(roi1_ctl, SetUpROI)

SetUpROI <- function(x){
  output <- x %>%
    dplyr::select(timecourse = V1) %>%
    mutate(TR = seq(2, 320, by = 2)) # TR counter
  mean <- mean(output$timecourse)
  output <- output %>% # demeaning the timecourse to enable cross-subject comparison
    mutate(timecourse_dm = timecourse - mean)
  return(output)
}




# MANIPULATE TASK ---------------------------------------------------------------------------------
# Usage: task_ctl <- lapply(task_ctl, SetUpTask)           # lapply(data, function)

SetUpTask <- function(x){
  
  # Onset, duration, and CS+/- indicator value
  output <- x %>%
    dplyr::select(onset = V1, duration = V2, value = V3) %>%
    mutate(stim = ifelse(value == -1, "CS-", "CS+"))
  
  # Column indicating whether it's the 1st, 2nd, 3rd, or 4th block of each type
  output[ , "num"] <- 0  # make additional last row of 0s, to fill in
  output[which(output[ ,3] == -1), 5] <- seq(1,4)
  output[which(output[ ,3] == 1), 5] <- seq(1,4)
  output$num <- as.character(output$num)
  
  # Column with CS+/- and block number in same cell
  output$block <- with(output, paste0(stim, num))
  
  # Sort to have CS-s in order then CS+s in order
  output <- output %>% arrange(block)
  
  return(output)
  
}




# FIND TR START OF BLOCKS IN TASK -------------------------------------------------------------------
# Usage: task_ctl <- lapply(task_ctl, StartBlocks)

StartBlocks <- function(x){
  
  output <- x
  
  # Vector of the TR numbers that begin each new "block" throughout the whole task
  #startbreaks <- c(2, 22, 42, 62, 82, 102, 122, 142, 162, 182, 202, 222, 242, 262, 282, 302)
  
  startbreaks <- seq(2, 320, by = 2) # NOW WILL LINE UP "ONSET" TIME BY CLOSEST TR, NOT "PARTICULAR START BLOCK TIMES"
  
  # Finding which TR numbers begin each CS+/- block for this individual
  output$startblock <- 0 # add column with this information
  i <- 1 # counter
  for (i in 1:8){ # 4 blocks of each stim type, so will do this 8 times
    y <- output[i, 1] # save the onset time
    absvals <- (abs(y - startbreaks)) # find the absolute value of the distance of the onset time to each potential startbreak time
    ind <- which.min(absvals) # find the minimum distance and save its index 
    output$startblock[i] <- startbreaks[ind] # in startblock col, save the TR at which that relevant block begins
    i <- i + 1 # add one to the counter to do for the next subsequent row
  }
  
  return(output)
  
}



# SUBSET PARTS OF TIMECOURSES TO FIND ACTIVITY PER TR IN EACH BLOCK ------------------------------------------
# Usage: 
#          - arg x is task list; arg y is roi1 list; arg z is roi2 list
#          - output is data list, which contains the list of activity per TR in each block
#          - need to apply the function over 3 lists and iterate simultaneously >> map2 is for 2; pmap is for >2
#          - data_ctl <- pmap(task_ctl, roi1_ctl, roi2_ctl, ActivityDuringBlocks)


ActivityDuringBlocks <- function(x, y, z){
  
  # Create new dataframe to put the values for each TR from each ROI
  stim_m <- rep("CS-", 28)
  stim_p <- rep("CS+", 28)
  stim <- c(stim_m, stim_p) # make stim col
  
  num_m <- rep(1:4, each = 7)
  num_p <- rep(1:4, each = 7)
  num <- c(num_m, num_p) # make num of block col (don't need, but might as well keep)
  
  roi1 <- rep(NA, 56) # make empty vals col
  roi2 <- rep(NA, 56) # make empty vals col
  
  data <- data.frame(stim, num, roi1, roi2) # make empty df to put vals of each relevant TR into
  
  # DEAL WITH ***ROI1*** VALUES (arg y) - TAKE DEMEANED TIMESERIES AS SUBSET RELEVANT TRS INTO THE "DATA" DF
 
  i <- 1 # counter for CS+/CS- blocks (8 total)
  f <- 1 # f and g are counters for the rows in "data" df you're going to stick the relevant data
  g <- 7
  
  for (i in 1:8){ 
    
    a <- which(y[ , 2] == x$startblock[i]) # which TR in roi1 df # equals "starblock" TR for that block from task df - "a" is the ROW NUMBER in roi1 df
    
    data[f:g, "roi1"] <- y[(a+3):(a+9), 3]  # take roi1 df "y" and take TRs 4-10 of that block for col 3 (demeaned timeseries) and stick in "data" df
    f <- f + 7
    g <- g + 7 
    
    i <- i + 1 # CS+/CS- block counter
  }
  
  
  
  # DEAL WITH ***ROI2*** VALUES (arg z) - TAKE DEMEANED TIMESERIES AS SUBSET RELEVANT TRS INTO THE "DATA" DF
  
  i <- 1 # counter for CS+/CS- blocks (8 total)
  f <- 1 # f and g are counters for the rows in "data" df you're going to stick the relevant data
  g <- 7
  
  for (i in 1:8){ 
    
    a <- which(z[ , 2] == x$startblock[i]) # which TR in roi2 df # equals "starblock" TR for that block from task df - "a" is the ROW NUMBER in roi2 df
    
    data[f:g, "roi2"] <- z[(a+3):(a+9), 3]  # take roi1 df "z" and take TRs 4-10 of that block for col 3 (demeaned timeseries) and stick in "data" df
    f <- f + 7
    g <- g + 7
    
    i <- i + 1 # CS+/CS- block counter
  }
  
  return(data)
}







