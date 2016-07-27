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



##################################################################
####                    Batch 1                               ####
##################################################################

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

ggplot(
  a1_agg
  , aes(x = X, y = value, group = cueing, color = cueing)
) +
  geom_line() +
  facet_grid(cue_modality ~ target_modality + laterality)+
  scale_y_reverse()
  



##################################################################
####                    Batch 2                               ####
##################################################################

batch2_list = c("e12", "e16", "e17", "e20", "e21", "e22")

# batch 1 df
a2 = a[a$id %in% batch2_list,]

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

ggplot(
  a2_agg
  , aes(x = X, y = value, group = cueing, color = cueing)
) +
  geom_line() +
  facet_grid(cue_modality ~ target_modality + laterality)+
  scale_y_reverse()



##################################################################
####                    Together                              ####
##################################################################

aa = rbind(a1, a2)

aa_agg = aggregate(value ~ variable + X, data = aa, FUN = mean)

# add cueing/cue/target/laterality columns
aa_agg$cue_modality = ifelse(
  substr(aa_agg$variable,3,3) == "T"
  , "tactile"
  , "visual"
)
aa_agg$target_modality = ifelse(
  substr(aa_agg$variable,4,4) == "T"
  , "tactile"
  , "visual"  
)
aa_agg$laterality = ifelse(
  substr(aa_agg$variable,2,2) == "C"
  , "contra"
  , "ipsi"  
)
aa_agg$cueing = ifelse(
  substr(aa_agg$variable,1,1) == "C"
  , "cued"
  , "uncued"  
)

# change factor order (note: this does not change DF)
aa_agg$cueing = factor(aa_agg$cueing, c("uncued", "cued"))

# convert to factor
aa_agg$cue_modality = factor(aa_agg$cue_modality)
aa_agg$target_modality = factor(aa_agg$target_modality)
aa_agg$laterality = factor(aa_agg$laterality)

# change names
levels(aa_agg$cue_modality) = c('Tactile Cue', 'Visual Cue')
levels(aa_agg$target_modality) = c('Tactile Target\n(C3/4)', 'Visual Target\n(PO7/8)')
levels(aa_agg$laterality) = c("Contralateral", "Ipsilateral")
levels(aa_agg$cueing) = c("Uncued", "Cued")

ggplot(
  aa_agg
  , aes(x = X, y = value, group = cueing, color = cueing)
) +
  geom_line() +
  facet_grid(cue_modality ~ target_modality + laterality)+
  scale_y_reverse()+
  geom_vline(
    xintercept =  0
    , linetype = 2
  )+ 
  geom_hline(
    yintercept = 0
  )+
  labs(color = "Cueing")