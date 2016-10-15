# libraries
library(grid)
library(gridExtra)
library(plyr)
library(ggplot2) 
library(stringr)
library(ez)
library(nlme)

# old or new data
NEW = TRUE
BOTH = TRUE

if (NEW | BOTH) {
  # directory
  setwd("/Users/ghislaindentremont/Documents/Multimodal_IOR/Ghis/TXT/new_data")
  
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
  d2 = a
  
  # get rid of participants with no EEG
  d2 = d2[d2$id != "e37" & d2$id != "e38",]
  
}

if (!NEW | BOTH){
  # directory
  setwd("/Users/ghislaindentremont/Documents/Multimodal_IOR/Ghis/TXT/BeforeSummer_ForAnalysis")
  
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
  d1 = b 
  
  # factor this 
  d1$id = as.factor(d1$id)
  
  # define recovered data files
  d1$recovered = FALSE
  d1[d1$id == 'e03' | d1$id == 'e04' | d1$id == 'e05' |d1$id == 'e06' | d1$id == 'e07',]$recovered = TRUE
  
  # add some columns
  d1[d1$recovered,]$target_response_rt = d1[d1$recovered,]$response_time - d1[d1$recovered,]$target_on_time 
  
  # make NA factor 
  d1$block = as.factor(ifelse(is.na(d1$block), "NA", d1$block))
  
}

# join new and old data if applicable otherwise create copy of single df
if (BOTH) {
  d = rbind.fill(d1, d2)
} else if (NEW) {
  d = d2
} else {
  d = d1
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

if (!NEW | BOTH){
  d[is.na(d$pre_target_response),]$pre_target_response = FALSE
  d[is.na(d$target_type),]$target_type = "target" 
  d[is.na(d$recovered),]$recovered = FALSE
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
P_list = unique(f$id)
print(P_list)
save(P_list, file = "../P_list.Rdata")
# how many?
length(P_list)

# create data frame for people that actually end up being used in analysis
d_left = d[d$id %in% P_list,]

# proportions of exclusions
length_id = aggregate(critical_blink ~ id, data = d_left, FUN = length)
length_id2 = aggregate(critical_blink ~ id, data = d_left[!d_left$recovered,], FUN = length)

blink_excl = aggregate(critical_blink ~ id,data = d_left, FUN = sum)
blink_excl$prop = blink_excl$critical_blink / length_id$critical_blink
print("Proportion of blinks: ")
print(blink_excl)
mean(blink_excl$prop)

saccade_excl = aggregate(critical_saccade ~ id,data = d_left, FUN = sum)
saccade_excl$prop = saccade_excl$critical_saccade / length_id$critical_blink
print("Proportion of saccades: ")
print(saccade_excl)
mean(saccade_excl$prop)

# no recovered data sets here
pre_target_response_excl = aggregate(pre_target_response ~ id,data = d_left[!d_left$recovered,], FUN = sum)
pre_target_response_excl$prop = pre_target_response_excl$pre_target_response / length_id2$critical_blink
print("Proportion of pre target responses: ")
print(pre_target_response_excl)
mean(pre_target_response_excl$prop)



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



#### Analysis ####
# ezA = ezANOVA(
#   IOR
#   , target_response_rt
#   , id
#   , .(cued, cue_modality, target_modality)
#   , return_aov = T
# )
# print(ezA)

# ezP = ezPlot(
#   IOR
#   , target_response_rt
#   , id
#   , .(cued, cue_modality, target_modality)
#   , split = .(cued)
#   , col = .(target_modality)
#   # , diff = factor(cued)
#   , x = .(cue_modality)
# )
# print(ezP)

# ezS = ezStats(
#   IOR
#   , target_response_rt
#   , id
#   , .(cued, cue_modality, target_modality)
#   # , diff = cued
# )
# print(ezS)


m = aov(target_response_rt ~
          cue_modality*target_modality*cued 
        + Error(id/(cue_modality*target_modality*cued))
        , data = IOR)
m_summary = summary(m)

# # Normality Assumption - NOT SURE THIS IS RIGHT
# residz = NULL
# for (i in 1:7) {
#   m_stuff = proj(m)
#   temp = m_stuff[[2+i]][,"Residuals"]
#   temp = temp[seq(1,184,8)]
#   residz = c(residz, temp)
# }
# qqnorm(residz)
# qqline(residz)

# # residual standard error (from ez table) is given by sqrt of MSE of of cueing of full model
# res_st_err = sqrt(318)

# equaltion: LSD = t(alpha/2, N-a) * sqrt(2*MSE/n), when n1=n2
alpha_over_2 = 0.025

# GET MSE from one-way ANOVA
IOR_cueing = aggregate(target_response_rt ~ cued + id, data = IOR, FUN = mean) 
# model
m2 = aov(target_response_rt ~ cued
         + Error(id/cued)
         , data = IOR_cueing)
m2_summary = summary(m2)

MSE = m2_summary$`Error: id:cued`[1][[1]][[3]][2]
df = m2_summary$`Error: id:cued`[1][[1]][[1]][2]
n = df + 1
# double check 
n == length(unique(IOR$id))

# get critical t value and LSD
t_crit = qt(alpha_over_2, df)
LSD = abs(t_crit * sqrt(2*MSE/n))

# now identify your IOR effects 
IOR_effects = aggregate(target_response_rt ~ cue_modality + target_modality + id ,data = IOR, FUN = diff)
IOR_CIs = aggregate(target_response_rt ~ cue_modality + target_modality,data = IOR_effects, FUN = mean)
# change name
names(IOR_CIs)[3] = "M"

# plot 3 way 
# add ci to data frame
IOR_CIs$CI = LSD

# plot 3 way interaction 
IOR_CIs$cue_modality = as.factor(IOR_CIs$cue_modality)
levels(IOR_CIs$cue_modality) = c("Tactile Cue", "Visual Cue")
IOR_CIs$target_modality = as.factor(IOR_CIs$target_modality)
levels(IOR_CIs$target_modality) = c("Tactile Target", "Visual Target")

# generate plot
gg = ggplot(IOR_CIs, aes(x = cue_modality
                                 ,y = M
                                 , group = target_modality
                                 , fill = target_modality
                                 , color = target_modality
)
)+
  geom_line(size = 1)+
  geom_errorbar(aes(ymin = M - CI, ymax = M + CI)
                , width = 0.1
                , size = 1
  )+
  labs(x = "Cue Modality",y = "IOR Effect: Cued - Uncued (ms)", color = "Target Modality")+
  geom_hline(yintercept = 0, size = 1, linetype = "dashed")+
  theme_gray(base_size = 30)+
  theme(panel.grid.major = element_line(size = 1.5)
   ,panel.grid.minor = element_line(size = 1)) 

print(gg)



