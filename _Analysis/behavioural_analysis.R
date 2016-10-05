# libraries
library(grid)
library(gridExtra)
library(plyr)
library(ggplot2) 
library(stringr)

# directory
setwd("/Volumes/Seagate Backup Plus Drive/Experiments/multimodal_ior/_Data/forR/BeforeSummer_ForAnalysis")

# old data 
a = ldply(
  .data = list.files(
    pattern = "_data.txt"
  )
  , .fun = function(x) {
    read.table(
      x
      , header = T
    )
  }
)

# make copy
b=a

# switch left from right 
b$is_right = b$cue_location == 'right'

b$cue_location = factor(
  ifelse( b$is_right
          , 'left'
          , 'right'
  )
)

# now do it for target 
b$is_right_target = b$target_location == 'right' 

b$target_location = factor(
  ifelse( b$is_right_target
          , 'left'
          , 'right'
  )
)  

# make another copy
d = b 

# factor this 
d$id = as.factor(d$id)

# define recovered data files
d$recovered = FALSE
d[d$id == 'e03' | d$id == 'e04' | d$id == 'e05' |d$id == 'e06' | d$id == 'e07',]$recovered = TRUE

# add some columns 
d[d$recovered,]$target_response_rt = d[d$recovered,]$response_time - d[d$recovered,]$target_on_time 

# create cued factor 
d$cued = FALSE
d$cued[(d$target_location == "right" & d$cue_location == "right") | (d$target_location == "left" & d$cue_location == "left")] = TRUE

# make NA factor 
d$block = as.factor(ifelse(is.na(d$block), "NA", d$block))

print("getting rid of practice trials...")
d = d[d$block != "practice",]

summarize_d = summary(d) 
print(summarize_d)

# overall RT distribution
hist(d$target_response_rt, breaks = 50)
abline(v = 100)


# proportions of exclusions
length_id = aggregate(critical_blink ~ id, data = d, FUN = length)

blink_excl = aggregate(critical_blink ~ id,data = d, FUN = sum)
blink_props_id = blink_excl$critical_blink / length_id$critical_blink
print("Proportion of blinks: ")
mean(blink_props_id)

saccade_excl = aggregate(critical_saccade ~ id,data = d, FUN = sum)
saccade_props_id = saccade_excl$critical_saccade / length_id$critical_blink
print("Proportion of saccades: ")
mean(saccade_props_id)

# no recovered data sets here
pre_target_response_excl = aggregate(pre_target_response ~ id,data = d, FUN = sum)
pre_target_response_props_id = pre_target_response_excl$pre_target_response / length_id$critical_blink
print("Proportion of pre target responses: ")
mean(pre_target_response_props_id)


# exclusions
e = d
e[is.na(e$pre_target_response),]$pre_target_response = FALSE
e = e[!e$pre_target_response,] 
e = e[!e$critical_saccade,]
e = e[!e$critical_blink,]
e[is.na(e$target_type),]$target_type = "target"
e = e[e$target_type == "target",]
e = e[!is.na(e$target_response_rt),]
e = e[e$target_response_rt >= 100,]
e$count = TRUE

# count trials per condition 
trial_per_condition_post = aggregate(count ~ cue_modality + cue_location + target_modality + target_location + id, data = e, FUN = sum)
# who has fewer than 10 trials in one condition or more?
too_few = unique(trial_per_condition_post[trial_per_condition_post$count < 10,]$id)
print(too_few)
# get rid of them
f = e
f = f[!(f$id %in% too_few), ]
# who's left?
unique(f$id)
# how many?
length(unique(f$id))

# look at IOR by P
IOR = aggregate(target_response_rt ~ cued + cue_modality + target_modality + id, data = f, FUN = mean)

# IOR 
# cued - uncued
IOR_effects = aggregate(target_response_rt ~ cue_modality + target_modality + id, data = IOR, FUN = diff)
names(IOR_effects)[4] = "IOR"

# group IOR
IOR_effecs_group = aggregate(IOR ~ cue_modality + target_modality, data = IOR_effects, FUN = mean)

