###############################################################################
# DIFFERENTIAL GRAPH FUNCTIONS FOR FEAR COND
###############################################################################

# 1) get just the differentials and separate by group... reshaping data 

diff_ctl <- function(x){
  output <- x %>%
    dplyr::select(subject, trauma, CSPRvPN_1, CSPRvPN_2, CSPRvPN_3, CSPRvPN_4) %>%
    gather(key = "block", value = "response", 3:6 ) %>%
    filter(trauma == 0)
  return(output)
}

diff_tr <- function(x){
  output <- x %>%
    dplyr::select(subject, trauma, CSPRvPN_1, CSPRvPN_2, CSPRvPN_3, CSPRvPN_4) %>%
    gather(key = "block", value = "response", 3:6 ) %>%
    filter(trauma == 1)
  return(output)
}

# 2) summaries of differentials by group

diff_ctl_sum <- function(x){
  output <- data.frame(mean <- tapply(x$response,x$block,mean), # calculating the stats and putting it in a dataframe
                       n <- tapply(x$response,x$block,length),
                       sd <- tapply(x$response,x$block,sd))
  output$se <- output$sd/sqrt(output$n) # calculating standard error
  output["trauma"] <- 0 #adding trauma/grouping col into summary file
  output2 <- setDT(output, keep.rownames = TRUE)[] #moves the rownames to the first col
  colnames(output2)=c("block", "mean", "n", "sd", "se", "trauma") # naming
  return(output2)
}

diff_tr_sum <- function(x){
  output <- data.frame(mean <- tapply(x$response,x$block,mean),
                       n <- tapply(x$response,x$block,length),
                       sd <- tapply(x$response,x$block,sd))
  output$se <- output$sd/sqrt(output$n)
  output["trauma"] <- 1 #adding trauma/grouping col into summary file
  output2 <- setDT(output, keep.rownames = TRUE)[] #moves the rownames to the first col
  colnames(output2)=c("block", "mean", "n", "sd", "se", "trauma") 
  return(output2)
}

# 3) combine summaries

combine_summary_diff <- function(x,y){
  output <- rbind(x, y) # combines
  return(output)
}

# 4) plot function

diff_plot <- function(x){
  ggplot(data = x[[1]], aes(x=block, y=mean)) +
    theme_classic(base_size = 14) + # was 12; made bigger
    geom_hline(yintercept = 0, color = "gray", linetype = "dashed") +
    theme(panel.border = element_rect(fill=NA)) +
    geom_point(size=4.5) +
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width = 0.3, size = 1.5) +
    xlab("Block") +
    theme(axis.text = element_text(size = 12), # size of axis labels and titles
          axis.title = element_text(size = 14)) +
    scale_x_discrete(labels=block_names) +
    ylab("BOLD response") +
    ggtitle(paste0(x[[2]],": CS+R minus CS+")) +
    theme(plot.title = element_text(hjust = 0.5)) +                # centers title, to make bold,
    facet_grid(~trauma, labeller = as_labeller(trauma_names)) +    # face = "bold" (goes after hjust)
    theme(strip.text = element_text(size = 14)) # facet label sizes
}




###############################################################################
# STIMULI GRAPH FUNCTIONS FOR FEAR COND
###############################################################################
# 1) get just the stim and separate by group

stim_ctl <- function(x){
  output <- x %>%
    dplyr::select(subject, trauma, CSPR_1, CSPR_2, CSPR_3, CSPR_4, CSPN_1, CSPN_2, CSPN_3, CSPN_4) %>%
    gather(key = "stim", value = "response", 3:10 ) %>%
    filter(trauma == 0)
  return(output)
}

stim_tr <- function(x){
  output <- x %>%
    dplyr::select(subject, trauma, CSPR_1, CSPR_2, CSPR_3, CSPR_4, CSPN_1, CSPN_2, CSPN_3, CSPN_4) %>%
    gather(key = "stim", value = "response", 3:10 ) %>%
    filter(trauma == 1)
  return(output)
}

# 2) summaries of differentials by group

stim_ctl_sum <- function(x){
  output <- data.frame(mean <- tapply(x$response,x$stim,mean),
                       n <- tapply(x$response,x$stim,length),
                       sd <- tapply(x$response,x$stim,sd))
  output$se <- output$sd/sqrt(output$n)
  output["trauma"] <- 0 #adding trauma/grouping col into summary file
  output2 <- setDT(output, keep.rownames = TRUE)[] #moves the rownames to the first col
  colnames(output2)=c("stim", "mean", "n", "sd", "se", "trauma")
  return(output2)
}

stim_tr_sum <- function(x){
  output <- data.frame(mean <- tapply(x$response,x$stim,mean),
                       n <- tapply(x$response,x$stim,length),
                       sd <- tapply(x$response,x$stim,sd))
  output$se <- output$sd/sqrt(output$n)
  output["trauma"] <- 1 #adding trauma/grouping col into summary file
  output2 <- setDT(output, keep.rownames = TRUE)[] #moves the rownames to the first col
  colnames(output2)=c("stim", "mean", "n", "sd", "se", "trauma") 
  return(output2)
}

# 3) combine summaries

combine_summary_stim <- function(x,y){
  output <- rbind(x, y) 
  return(output)
}

#add stimsign col
addstimsign <- function(x){
  output <- x %>%
    mutate(stimsign = 
             if_else(substr(stim, 1, 4) == "CSPR", "CS+R", "CS+"))
  return(output)
}

#add block col
addblockcol <- function(x){
  output <- x %>%
    mutate(block = 
             if_else(substr(stim, 5, 6) == "_1", 1,
                     if_else(substr(stim, 5, 6) == "_2", 2,
                             if_else(substr(stim, 5, 6) == "_3", 3, 4))))
  return(output)
}

# 4) plot function

stim_plot <- function(x){
  ggplot(data = x[[1]], aes(x=block, y=mean, color=stimsign)) +
    theme_classic(base_size = 14) +
    geom_hline(yintercept = 0, color = "gray", linetype = "dashed") +
    theme(panel.border = element_rect(fill=NA)) +
    geom_point(size=4.5) +
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width = 0.3, size = 1.5, show.legend=FALSE) +
    xlab("Block") +
    theme(axis.text = element_text(size = 12), # size of axis labels and titles
          axis.title = element_text(size = 14)) +
    ylab("BOLD response") +
    ggtitle(paste0(x[[2]],": CS+R and CS+")) +
    theme(plot.title = element_text(hjust = 0.5)) + #centers title
    scale_color_manual(name = "Stimuli", values = c("CS+" = "red", "CS+R" = "red4")) +
    guides(color = guide_legend(override.aes = list(size=2, shape=15), reverse=TRUE)) + #legend squares
    #guides(color = guide_legend(override.aes = list(size=3), reverse=TRUE)) + #legend circles
    theme(legend.title=element_text(size=12)) + # was 10
    theme(legend.text=element_text(size=12)) + # was 8
    facet_grid(~trauma, labeller = as_labeller(trauma_names)) +
    theme(strip.text =  element_text(size = 14)) # facet label sizes
}




