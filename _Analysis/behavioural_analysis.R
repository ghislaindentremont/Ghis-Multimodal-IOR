# libraries
library(grid)
library(gridExtra)
library(plyr)
library(ggplot2) 
library(stringr)

# old or new data
NEW = TRUE

if (NEW) {
  # directory
  setwd("/Volumes/Seagate Backup Plus Drive/Experiments/multimodal_ior/_Data/forR/new_data")
  
  # old data 
  a = ldply(
    .data = list.files(
      pattern = "_data.txt"
    )
    , .fun = function(x) {
      read.table(
        x
        , header = T
        , stringsAsFactors = F
      )
    }
  )
  
  # copy
  d = a
  
} else{
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
  
  # make NA factor 
  d$block = as.factor(ifelse(is.na(d$block), "NA", d$block))
  
}



# create cued factor 
d$cued = FALSE
d$cued[(d$target_location == "right" & d$cue_location == "right") | (d$target_location == "left" & d$cue_location == "left")] = TRUE



#### Exclusions ####
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
blink_excl$prop = blink_excl$critical_blink / length_id$critical_blink
print("Proportion of blinks: ")
print(blink_excl)
mean(blink_excl$prop)

saccade_excl = aggregate(critical_saccade ~ id,data = d, FUN = sum)
saccade_excl$prop = saccade_excl$critical_saccade / length_id$critical_blink
print("Proportion of saccades: ")
print(saccade_excl)
mean(saccade_excl$prop)

# no recovered data sets here
pre_target_response_excl = aggregate(pre_target_response ~ id,data = d, FUN = sum)
pre_target_response_excl$prop = pre_target_response_excl$pre_target_response / length_id$critical_blink
print("Proportion of pre target responses: ")
print(pre_target_response_excl)
mean(pre_target_response_excl$prop)

if (!NEW){
  e[is.na(e$pre_target_response),]$pre_target_response = FALSE
  e[is.na(e$target_type),]$target_type = "target" 
}

# exclusions
e = d
e = e[!e$pre_target_response,] 
e = e[!e$critical_saccade,]
e = e[!e$critical_blink,]
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



#### IOR ####
# look at IOR by P
IOR = aggregate(target_response_rt ~ cued + cue_modality + target_modality + id, data = f, FUN = mean)

# IOR 
# cued - uncued
IOR_effects = aggregate(target_response_rt ~ cue_modality + target_modality + id, data = IOR, FUN = diff)
names(IOR_effects)[4] = "IOR"
print(IOR_effects)

# group IOR
IOR_effecs_group = aggregate(IOR ~ cue_modality + target_modality, data = IOR_effects, FUN = mean)
print(IOR_effecs_group)

