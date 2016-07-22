library(ggplot2)
library(plyr)
library(reshape)

setwd("~/Documents/Multimodal_IOR/condition_averages/")

a = ldply(
  .data = list.files(
    path = "."
  )
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

batch1_list = c("e01", "e02", "e03", "e04", "e05", "e06", "p06", "e27")

# batch 1 df
a1 = a[a$id %in% batch1_list,]

# melt 
a1 = melt(a1[,-c(2:3)])

# add time variable
a1$X = -200:500

# average across participants
a1_agg = aggregate(value ~ variable + X, data = a1, FUN = mean)

# add cueing/cue/target/laterality columns
a1_agg$cue_modality = ifelse(
  substr(a1_agg$variable,3,3) == "T"
  , "tactile"
  , "visual"
)
a1_agg$target_modality = ifelse(
  substr(a1_agg$variable,4,4) == "T"
  , "tactile"
  , "visual"  
)
a1_agg$laterality = ifelse(
  substr(a1_agg$variable,2,2) == "C"
  , "contra"
  , "ipsi"  
)
a1_agg$cueing = ifelse(
  substr(a1_agg$variable,1,1) == "C"
  , "cued"
  , "uncued"  
)



##################################################################
####                    By Condition                          ####
##################################################################

ggplot(
  a1_agg
  , aes(x = X, y = value, group = cueing, color = cueing)
) +
  geom_line() +
  facet_grid(cue_modality ~ target_modality + laterality)+
  scale_y_reverse()
  
