# libraries
library(grid)
library(gridExtra)
library(plyr)
library(ggplot2) 
library(stringr)
library(ez)


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
  
  # get rid of participants with no EEG or other exclusions
  d2 = d2[d2$id != "e37" & d2$id != "e38" & d2$id != "e43",]
  
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

# demographics
sex = aggregate(sex ~ id, d_left, unique)
mean(sex$sex == "f")  # approximate proportion of women
age = aggregate(age ~ id, d_left, unique)
median(age$age)  # approximate proportion 

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

print("Proportion 0-100 ms responses: ")
mean(d_left[!is.na(d_left$target_response_rt),]$target_response_rt < 100)

# no recovered data sets here
pre_target_response_excl = aggregate(pre_target_response ~ id,data = d_left[!d_left$recovered,], FUN = sum)
pre_target_response_excl$prop = pre_target_response_excl$pre_target_response / length_id2$critical_blink
print("Proportion of pre target responses: ")
print(pre_target_response_excl)
mean(pre_target_response_excl$prop)



#### IOR ####
# look at IOR by P
IOR = aggregate(target_response_rt ~ cued + cue_modality + target_modality + id, data = f, FUN = mean)
IOR_summ = aggregate(target_response_rt ~ cued + cue_modality + target_modality, data = IOR, FUN = mean)
IOR_summ$SD = aggregate(target_response_rt ~ cued + cue_modality + target_modality, data = IOR, FUN = sd)$target_response_rt

# IOR 
# cued - uncued
IOR_effects = aggregate(target_response_rt ~ cue_modality + target_modality + id, data = IOR, FUN = diff)
names(IOR_effects)[4] = "IOR"
print(IOR_effects)

# group IOR
IOR_effecs_group = aggregate(IOR ~ cue_modality + target_modality, data = IOR_effects, FUN = mean)
print(IOR_effecs_group)



#### Analysis ####
ezA = ezANOVA(
  IOR
  , target_response_rt
  , id
  , .(cued, cue_modality, target_modality)
  , return_aov = T
)
print(ezA)

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

ezS = ezStats(
  IOR
  , target_response_rt
  , id
  , within = .(cued, cue_modality, target_modality)
  # , diff = cued
)
print(ezS)

m = aov(target_response_rt ~
          cue_modality*target_modality*cued 
        + Error(id/(cue_modality*target_modality*cued))
        , data = IOR)
m_summary = summary(m)

# Normality Assumption - NOT SURE THIS IS RIGHT
residz = NULL
for (i in 1:7) {
  m_stuff = proj(m)
  temp = m_stuff[[2+i]][,"Residuals"]
  temp = temp[seq(1,184,8)]
  residz = c(residz, temp)
}
qqnorm(residz)
qqline(residz)


#--------------------------------- SEM (L & M) ----------------------------------------#
alpha_over_2 = 0.025

IOR$mix_factor = paste(IOR$cued, IOR$cue_modality, IOR$target_modality)
# model
m2 = aov(target_response_rt ~ mix_factor
         + Error(id/mix_factor)
         , data = IOR)
m2_summary = summary(m2)

MSE = m2_summary$`Error: id:mix_factor`[1][[1]][[3]][2]
df = m2_summary$`Error: id:mix_factor`[1][[1]][[1]][2]
n = length(unique(IOR$id))

SEM_LM = sqrt(MSE/n)

# get critical t value and LSD
t_crit = qt(1-alpha_over_2, df)
ME_LM = abs(t_crit * SEM_LM)
LSD = sqrt(2) * ME_LM
#--------------------------------- SEM (L & M) ----------------------------------------#


#--------------------------------- Circularity ----------------------------------------#
# NOTE: ciruclarity = sphericity
IOR$mix_factor_num = as.numeric(factor(IOR$mix_factor))
pairings = combn(unique(IOR$mix_factor_num),2)

alpha = 0.05
SEs = NULL
MEs = NULL
ids = NULL
MEANs = NULL
for (i in 1:ncol(pairings)) {
  pair = pairings[,i]
  
  level1 = pair[1]
  level2 = pair[2]
  
  data1 = IOR[IOR$mix_factor_num == level1,]$target_response_rt
  data2 = IOR[IOR$mix_factor_num == level2,]$target_response_rt
  
  dat = data1-data2
  
  MEAN = mean(dat)
  SD = sd(dat)
  N = length(dat)
  DF = N - 1
  T_CRIT = qt(1-alpha/2, DF)
  
  SE = SD/sqrt(N)
  
  ME = T_CRIT * SE
  
  MEANs = c(MEANs, MEAN)
  MEs = c(MEs, ME)
  SEs = c(SEs, SE)
  ids = c(ids, paste(level1, level2))
}

hist(SEs/sqrt(2)) # scale!
abline(v = SEM_LM)

# also look at epsilon
ezEpsi = ezANOVA(data = IOR
        , dv = target_response_rt
        , wid = id
        , within = mix_factor)
ezEpsi$`Sphericity Corrections`$GGe

pairdiffs = data.frame(ids, MEs, SEs, MEANs)
gg2 = ggplot(pairdiffs, aes(x = ids, y = MEANs*1e6))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(ymin = (MEANs - SEs)*1e6, ymax = (MEANs + SEs)*1e6)
                , width = 0.1
                , size = 1
  )+
  theme_gray(base_size = 30)+
  theme(panel.grid.major = element_line(size = 1.5)
        ,panel.grid.minor = element_line(size = 1)) 
print(gg2)

pSEs = pairdiffs$SEs
SEM_LM_estimate = sqrt( mean((pSEs/sqrt(2))^2) )
abs(SEM_LM - SEM_LM_estimate) < 0.0001
#--------------------------------- Circularity ----------------------------------------#


#--------------------------------- condition-wise--------------------------------------#
# get individual confidence intervals
IOR_cueing = aggregate(target_response_rt ~ id + cue_modality + target_modality, data = IOR, FUN = diff) 
IOR_cueing$mix_factor = paste(IOR_cueing$cue_modality, IOR_cueing$target_modality)

# TEST main effect of cueing 
test1 = aggregate(target_response_rt ~ id + cued, data = IOR, FUN = mean) 
test1m = summary(aov(target_response_rt ~ cued
         + Error(id/cued)
         , data = test1))
MSEtest1 = test1m$`Error: id:cued`[1][[1]][[3]][2]
m_summary$`Error: id:cued`[1][[1]][[3]][2]/4  # devide by 4

# TEST: cue modality X cueing
test2 = aggregate(target_response_rt ~ id + cue_modality, data = IOR_cueing, FUN = mean) 
test2m = summary(aov(target_response_rt ~ cue_modality
            + Error(id/cue_modality)
            , data = test2))
MSEtest2 = test2m$`Error: id:cue_modality`[1][[1]][[3]][2]
m_summary$`Error: id:cue_modality:cued`[1][[1]][[3]][2]  # no need to devide 

# TEST: cue modality X cueing X target_modality
test3 = aggregate(target_response_rt ~ id + cue_modality, data = IOR_cueing, FUN = diff) 
test3m = summary(aov(target_response_rt ~ cue_modality
                     + Error(id/cue_modality)
                     , data = test3))
MSEtest3 = test3m$`Error: id:cue_modality`[1][[1]][[3]][2]
m_summary$`Error: id:cue_modality:target_modality:cued`[1][[1]][[3]][2]*4  # need to multiply by 4


CW = ddply(
  .data = IOR_cueing
  , .variables = .(mix_factor, cue_modality, target_modality)
  , .fun = function(x,alpha = 0.05){
    dat = x$target_response
    
    MEAN = mean(dat)
    SD = sd(dat)
    n = length(dat)
    df = n - 1
    t_crit = qt(1-alpha/2, df)
    
    SE = SD/sqrt(n)
    
    ME = t_crit * SE
    
    return(c(ME = ME, M = MEAN, SE = SE))
  }
)

CW$LSD = LSD

# # plot 3 way interaction 
# NOTE: These are not necessarily the error bars I would want for interaction
# CW$cue_modality = as.factor(CW$cue_modality)
# levels(CW$cue_modality) = c("Tactile Cue", "Visual Cue")
# CW$target_modality = as.factor(CW$target_modality)
# levels(CW$target_modality) = c("Tactile Target", "Visual Target")

# focus on IOR effects
CW$mix_factor = as.factor(CW$mix_factor)
levels(CW$mix_factor) = c("Tactile\nTactile", "Tactile\nVisual", "Visual\nTactile", "Visual\nVisual")

# generate plot
gg = ggplot(CW, 
            aes(x = mix_factor
              # x = cue_modality
                         ,y = M
#                          , group = target_modality
#                          , fill = target_modality
#                          , color = target_modality
            )
)+
  # geom_line(size = 1)+
  geom_errorbar(aes(ymin = M - ME, ymax = M + ME)
                , width = 0.1
                , size = 1
  )+
  geom_point(size = 2)+
  # labs(x = "Cue Modality",y = "IOR Effect: Cued - Uncued (ms)", color = "Target Modality")+
  labs(x = "Cue/Target Pairing",y = "IOR Effect: Cued - Uncued (ms)")+
  geom_hline(yintercept = 0, size = 1, linetype = "dashed")+
  theme_gray(base_size = 30)+
  theme(panel.grid.major = element_line(size = 1.5)
        ,panel.grid.minor = element_line(size = 1)) 

print(gg)
#--------------------------------- condition-wise -------------------------------------#


#------------------------------------- pooled -----------------------------------------#
SEM_LM_cueing = sqrt( mean((CW$SE)^2) ) # don't devide by square 2 - because already differences
t_crit_cueing = qt(1-alpha_over_2, n-1)
ME_cueing = t_crit_cueing * SEM_LM_cueing

CW$ME_cueing = ME_cueing

# generate plot
gg = ggplot(CW, 
            aes(x = mix_factor
                ,y = M
            )
)+
  geom_errorbar(aes(ymin = M - ME_cueing, ymax = M + ME_cueing)
                , width = 0.1
                , size = 1
  )+
  geom_point(size = 2)+
  labs(x = "Cue/Target Pairing",y = "IOR Effect: Cued - Uncued (ms)")+
  geom_hline(yintercept = 0, size = 1, linetype = "dashed")+
  theme_gray(base_size = 30)+
  theme(panel.grid.major = element_line(size = 1.5)
        ,panel.grid.minor = element_line(size = 1)) 

print(gg)
#------------------------------------- pooled -----------------------------------------#


# NOTE: why is this not the same as taking the MSE from interaction term of full ANOVA?
# is it meaningful that it is roughly a thrid of it? 
#--------------------------------- interaction ----------------------------------------#
IOR_same = IOR_cueing
IOR_same$cue_modality = ifelse(IOR_cueing$cue_modality == IOR_cueing$target_modality, TRUE, FALSE)

# same modality minus different modality
IOR_int_id = aggregate(target_response_rt ~ target_modality + id, data = IOR_same, FUN =diff)
IOR_int_id$target_response_rt = -IOR_int_id$target_response_rt

# model
m3 = aov(target_response_rt ~ target_modality
         + Error(id/target_modality)
         , data = IOR_int_id)
m3_summary = summary(m3)

MSE3 = m3_summary$`Error: id:target_modality`[1][[1]][[3]][2]
df3 = m3_summary$`Error: id:target_modality`[1][[1]][[1]][2]

SEM_LM3 = sqrt(MSE3/n)
t_crit3 = qt(1-alpha_over_2, n-1)
ME3 = t_crit3 * SEM_LM3

IOR_int = aggregate(target_response_rt ~ target_modality, data = IOR_int_id, FUN =mean)
names(IOR_int)[2] = "M"
IOR_int$ME3 = ME3

IOR_int$target_modality = as.factor(IOR_int$target_modality)
levels(IOR_int$target_modality) = c("Tactile", "Visual")

# generate plot
gg = ggplot(IOR_int, 
            aes(x = target_modality
                ,y = M
            )
)+
  geom_errorbar(aes(ymin = M - ME3, ymax = M + ME3)
                , width = 0.1
                , size = 1
  )+
  geom_point(size = 2)+
  labs(x = "Target Modality",y = "Same Modality IOR - Different Modality IOR (ms)")+
  geom_hline(yintercept = 0, size = 1, linetype = "dashed")+
  theme_gray(base_size = 30)+
  theme(panel.grid.major = element_line(size = 1.5)
        ,panel.grid.minor = element_line(size = 1)) 

print(gg)

# MSE_int = m_summary$`Error: id:cue_modality:target_modality:cued`[1][[1]][[3]][2]
#---------------------------------- interaction ----------------------------------------#


