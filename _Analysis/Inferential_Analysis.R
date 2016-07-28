library(ggplot2)
library(plyr)
library(reshape)

setwd("~/Documents/Multimodal_IOR/P_averages")


##################################################################
####                  Grand Averages                          ####
##################################################################

b_list = list.files(
  path = "."
  , pattern = "grand_averages"
)
b_data = b_list[grep(".csv", b_list)]
b = ldply(
  .data = b_data
  , .fun = function(file) {
    df = read.csv(
      file
      , header = TRUE
    )
    P_col = data.frame(id = substr(file,16,18))
    to_return = cbind(P_col, df)
    return(to_return)
  }
)


#---------------------------- Function --------------------------#
grand_batches = function(df = b, batch = 0) {
  if (batch == 0) {
    batch_list = unique(b$id)
  } else if (batch == 1) {
    batch_list = c("e01", "e02", "e03", "e04", "e05", "e06", "p06", "e27")
  } else if (batch == 2) {
    batch_list =  c("e12", "e16", "e17", "e20", "e21", "e22")
  } else {
    stop("Unknow 'batch' argument. Must be 0, 1, or 2")
  }
  
  # batch 1 df
  b2 = b[b$id %in% batch_list,]
  
  # melt 
  b2 = melt(b2[,-c(2:3)])
  
  # add time variable
  b2$X = -200:500
  
  # average across participants
  b2_agg = aggregate(value ~ variable + X, data = b2, FUN = mean)
  
  # add cueing/cue/target/laterality columns
  b2_agg$target_modality = ifelse(
    substr(b2_agg$variable,4,4) == "T"
    , "tactile"
    , "visual"  
  )
  
  # convert to factor
  b2_agg$target_modality = factor(b2_agg$target_modality)
  
  # change names
  levels(b2_agg$target_modality) = c('Tactile Target\n(C3/4)', 'Visual Target\n(PO7/8)')
  
  gg = ggplot(
    b2_agg
    , aes(x = X, y = value)
  ) +
    geom_line() +
    facet_grid(. ~ target_modality)+
    geom_vline(
      xintercept =  0
      , linetype = 2
    )+ 
    geom_hline(
      yintercept = 0
    )+
    scale_x_continuous("Time (ms)")+
    scale_y_reverse("Voltage (V)")
  
  return(gg)
}
#---------------------------- Function --------------------------#


#----------------------------- Batch 1 --------------------------#
grand_batches(b, 1)
#----------------------------- Batch 1 --------------------------#


#----------------------------- Batch 2 --------------------------#
grand_batches(b, 2)
#----------------------------- Batch 2 --------------------------#


#----------------------------- Together -------------------------#
grand_batches(b, 0)
#----------------------------- Together -------------------------#



##################################################################
####                 Condition Averages                       ####
##################################################################

a_list = list.files(
  path = "."
  , pattern = "condition_averages"
)
a_data = a_list[grep(".csv", a_list)]
a = ldply(
  .data = a_data
  , .fun = function(file) {
      df = read.csv(
      file
      , header = TRUE
      )
      P_col = data.frame(id = substr(file,20,22))
      to_return = cbind(P_col, df)
    return(to_return)
  }
)


#---------------------------- Function --------------------------#
condition_batches = function(df = a, batch = 0) {
  if (batch == 0) {
    batch_list = unique(b$id)
  } else if (batch == 1) {
    batch_list = c("e01", "e02", "e03", "e04", "e05", "e06", "p06", "e27")
  } else if (batch == 2) {
    batch_list =  c("e12", "e16", "e17", "e20", "e21", "e22")
  } else {
    stop("Unknow 'batch' argument. Must be 0, 1, or 2")
  }
  
  # batch 1 df
  a2 = a[a$id %in% batch_list,]
  
  # melt 
  a2 = melt(a2[,-c(2:3)])
  
  # add time variable
  a2$X = -200:500
  
  # average across participants
  a2_agg = aggregate(value ~ variable + X, data = a2, FUN = mean)
  
  # add cueing/cue/target/laterality columns
  a2_agg$cue_modality = ifelse(
    substr(a2_agg$variable,3,3) == "T"
    , "tactile"
    , "visual"
  )
  a2_agg$target_modality = ifelse(
    substr(a2_agg$variable,4,4) == "T"
    , "tactile"
    , "visual"  
  )
  a2_agg$laterality = ifelse(
    substr(a2_agg$variable,2,2) == "C"
    , "contra"
    , "ipsi"  
  )
  a2_agg$cueing = ifelse(
    substr(a2_agg$variable,1,1) == "C"
    , "cued"
    , "uncued"  
  )
  
  # change factor order (note: this does not change DF)
  a2_agg$cueing = factor(a2_agg$cueing, c("uncued", "cued"))
  
  # convert to factor
  a2_agg$cue_modality = factor(a2_agg$cue_modality)
  a2_agg$target_modality = factor(a2_agg$target_modality)
  a2_agg$laterality = factor(a2_agg$laterality)
  
  # change names
  levels(a2_agg$cue_modality) = c('Tactile Cue', 'Visual Cue')
  levels(a2_agg$target_modality) = c('Tactile Target\n(C3/4)', 'Visual Target\n(PO7/8)')
  levels(a2_agg$laterality) = c("Contralateral", "Ipsilateral")
  levels(a2_agg$cueing) = c("Uncued", "Cued")
  
  gg = ggplot(
    a2_agg
    , aes(x = X, y = value, group = cueing, color = cueing)
  ) +
    geom_line() +
    facet_grid(cue_modality ~ target_modality + laterality)+
    geom_vline(
      xintercept =  0
      , linetype = 2
    )+ 
    geom_hline(
      yintercept = 0
    )+
    scale_x_continuous("Time (ms)")+
    scale_y_reverse("Voltage (V)")+
    labs(color = "Cueing")
  
  return(gg)
}
#---------------------------- Function --------------------------#


#----------------------------- Batch 1 --------------------------#
condition_batches(b, 1)
#----------------------------- Batch 1 --------------------------#


#----------------------------- Batch 2 --------------------------#
condition_batches(b, 2)
#----------------------------- Batch 2 --------------------------#


#----------------------------- Together -------------------------#
condition_batches(b, 0)
#----------------------------- Together -------------------------#

