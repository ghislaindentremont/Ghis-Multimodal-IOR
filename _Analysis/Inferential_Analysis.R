library(ggplot2)
library(plyr)
library(reshape)
library(scales)
library(ggthemes)

setwd("/Users/ghislaindentremont/Documents/Experiments/Multimodal_IOR/Ghis/P_analysis2")


##################################################################
####                  Grand Averages                          ####
##################################################################

b_list2 = list.files(
  path = "."
  , pattern = "grand_averages"
  , recursive = T
)
b_list = b_list2[grep(".csv", b_list2)]

load("/Users/ghislaindentremont/Documents/Experiments/Multimodal_IOR/Ghis/TXT/P_list.Rdata")
b_data = NULL
for (id in P_list) {
  b_data = c(b_data, b_list[grep(id, b_list)])
}

b = ldply(
  .data = b_data
  , .fun = function(file) {
    df = read.csv(
      file
      , header = TRUE
    )
    P_col = data.frame(id = substr(file,1,3))
    to_return = cbind(P_col, df)
    return(to_return)
  }
)

batch = 0
pre_honours = c("e12", "e16", "e17", "e20", "e22", "e02", "e03", "e04", "e05", "e06", "p06", "e27")
select_batch = function(batch) {
  if (batch == 0) {
    batch_list = unique(b$id)
  } else if (batch == 1) {
    batch_list = c("e02", "e03", "e04", "e05", "e06", "p06", "e27")
  } else if (batch == 2) {
    batch_list =  c("e12", "e16", "e17", "e20", "e22")
  } else if (batch == 3) {
    batch_list =  pre_honours
  } else {
    batch_list = P_list[!(P_list %in% pre_honours)]
  }
  
  return(batch_list)
}
batch_list = select_batch(batch)

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
  , aes(x = X, y = value*1e6)
) +
  geom_line(size = 1) +
  facet_grid(. ~ target_modality, scales = "free")+
  geom_vline(
    xintercept =  0
    , linetype = 2
  )+ 
  geom_hline(
    yintercept = 0
  )+
  scale_x_continuous("Time (ms)")+
  scale_y_reverse(paste("Voltage (", expression(u), "V)", sep = ""))+
  theme_gray(base_size = 30)+
  theme(panel.grid.major = element_line(size = 1.5)
        ,panel.grid.minor = element_line(size = 1)) 

print(gg)



##################################################################
####                      For Presentation                    ####
##################################################################

plot_grands = function(dat, vline1, vline2, vline3, vline4) {
  gg = ggplot(
    dat
    , aes(x = X, y = value*1e6)
  ) +
    geom_line(size = 1) +
    geom_vline(
      xintercept =  0
      , linetype = 2
    )+ 
    geom_hline(
      yintercept = 0
    )+
    geom_vline(xintercept = vline1, color = "red", size = 1)+
    geom_vline(xintercept = vline2, color = "red", size = 1)+
    geom_vline(xintercept = vline3, color = "red", size = 1)+
    geom_vline(xintercept = vline4, color = "red", size = 1)+
    scale_x_continuous("Time (ms)", limits = c(-100, 400))+
    scale_y_reverse(paste("Voltage (", expression(u), "V)", sep = ""))+
    theme_gray(base_size = 30)+
    theme(panel.grid.major = element_line(size = 1.5)
          ,panel.grid.minor = element_line(size = 1)) 
  
  return(gg)
}

# tactile
tactile_only = b2_agg[b2_agg$target_modality == 'Tactile Target\n(C3/4)',]
Plo = 25
Phi = 60
Nlo = 60
Nhi = 120
plot_grands(tactile_only, Plo, Phi, Nlo, Nhi)


# visual
visual_only = b2_agg[b2_agg$target_modality == 'Visual Target\n(PO7/8)',]
P1lo = 70
P1hi = 170
N1lo = 170
N1hi = 210
plot_grands(visual_only,P1lo, P1hi, N1lo, N1hi)



##################################################################
####                 Condition Averages                       ####
##################################################################

a_list2 = list.files(
  path = "."
  , pattern = "condition_averages"
  , recursive = T
)

a_list = a_list2[grep(".csv", a_list2)]

a_data = NULL
for (id in P_list) {
  a_data = c(a_data, a_list[grep(id, a_list)])
}

a = ldply(
  .data = a_data
  , .fun = function(file) {
    df = read.csv(
      file
      , header = TRUE
    )
    P_col = data.frame(id = substr(file,1,3), stringsAsFactors=FALSE)
    to_return = cbind(P_col, df)
    return(to_return)
  }
)

select_batch(batch)

# batch 1 df
a2 = a[a$id %in% batch_list,]

# melt 
a2 = melt(a2[,-c(2:3)])

# add time variable
a2$time = -200:500

# add cueing/cue/target/laterality columns
a2$cue_modality = ifelse(
  substr(a2$variable,3,3) == "T"
  , "tactile"
  , "visual"
)
a2$target_modality = ifelse(
  substr(a2$variable,4,4) == "T"
  , "tactile"
  , "visual"  
)
a2$laterality = ifelse(
  substr(a2$variable,2,2) == "C"
  , "contra"
  , "ipsi"  
)
a2$cueing = ifelse(
  substr(a2$variable,1,1) == "C"
  , "cued"
  , "uncued"  
)

# average across participants
a2_agg = aggregate(value ~ variable + time + cue_modality + target_modality + laterality + cueing, data = a2, FUN = mean)

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
  , aes(x = time, y = value*1e6, group = cueing, color = cueing)
) +
  geom_line(size = 1) +
  facet_grid(laterality ~ target_modality + cue_modality)+
  geom_vline(
    xintercept =  0
    , linetype = 2
  )+ 
  geom_hline(
    yintercept = 0
  )+
  scale_x_continuous("Time (ms)", limits = c(-100, 400))+
  scale_y_reverse(paste("Voltage (", expression(u), "V)", sep = ""))+
  labs(color = "Cueing")+
  theme_gray(base_size = 30)+
  theme(panel.grid.major = element_line(size = 1.5)
        ,panel.grid.minor = element_line(size = 1)) 

print(gg)



##################################################################
####           Collapse Across Laterality                     ####
##################################################################

a2_agg_visual = a2_agg[a2_agg$target_modality == 'Visual Target\n(PO7/8)',]
a2_agg2 = aggregate(value ~ time + cue_modality + cueing, data = a2_agg_visual, FUN = mean)

gg = ggplot(
  a2_agg2
  , aes(x = time, y = value*1e6, group = cueing, color = cueing)
) +
  geom_line(size = 1) +
  facet_grid(~cue_modality)+
  geom_vline(
    xintercept =  0
    , linetype = 2
  )+ 
  geom_hline(
    yintercept = 0
  )+
  scale_x_continuous("Time (ms)")+
  scale_y_reverse(paste("Voltage (", expression(u), "V)", sep = ""))+
  labs(color = "Cueing")+
  theme_gray(base_size = 30)+
  theme(panel.grid.major = element_line(size = 1.5)
        ,panel.grid.minor = element_line(size = 1)) 

print(gg)

a2_agg_tactile = a2_agg[a2_agg$target_modality == 'Tactile Target\n(C3/4)',]
a2_agg22 = aggregate(value ~ time + cue_modality + cueing, data = a2_agg_tactile, FUN = mean)

gg = ggplot(
  a2_agg22
  , aes(x = time, y = value*1e6, group = cueing, color = cueing)
) +
  geom_line(size = 1) +
  facet_grid(~cue_modality)+
  geom_vline(
    xintercept =  0
    , linetype = 2
  )+ 
  geom_hline(
    yintercept = 0
  )+
  scale_x_continuous("Time (ms)")+
  scale_y_reverse(paste("Voltage (", expression(u), "V)", sep = ""))+
  labs(color = "Cueing")+
  theme_gray(base_size = 30)+
  theme(panel.grid.major = element_line(size = 1.5)
        ,panel.grid.minor = element_line(size = 1)) 

print(gg)



##################################################################
####         Collapse Across Cue Modality                     ####
##################################################################

a2_agg3 = aggregate(value ~ time + cueing, data = a2_agg_visual, FUN = mean)

gg = ggplot(
  a2_agg3
  , aes(x = time, y = value*1e6, group = cueing, color = cueing)
) +
  geom_line(size = 1) +
  geom_vline(
    xintercept =  0
    , linetype = 2
  )+ 
  geom_hline(
    yintercept = 0
  )+
  scale_x_continuous("Time (ms)")+
  scale_y_reverse(paste("Voltage (", expression(u), "V)", sep = ""))+
  labs(color = "Cueing")+
  theme_gray(base_size = 30)+
  theme(panel.grid.major = element_line(size = 1.5)
        ,panel.grid.minor = element_line(size = 1)) 

print(gg)




##################################################################
####                     Analysis                             ####
##################################################################

load("/Users/ghislaindentremont/Documents/Experiments/Multimodal_IOR/Ghis/TXT/IOR_effects.Rdata")

do_aov = function(component, target_modality, lower_bound, upper_bound) {
  # P45
  eeg_comp = a2[a2$target_modality == target_modality
                & a2$time >= lower_bound
                & a2$time <= upper_bound, ]
  
  # get average voltage in window for each condition (4) x cueing (2) x participant (n)
  eeg_comp_agg = aggregate(value ~ cue_modality + laterality + cueing + id, data = eeg_comp, FUN = mean)
  
  # get means 
  eeg_comp_agg_means = aggregate(value ~ cue_modality + laterality + cueing, data = eeg_comp_agg, FUN = mean)
  print(eeg_comp_agg_means)
  
  # get sds
  eeg_comp_agg_sds = aggregate(value ~ cue_modality + laterality + cueing, data = eeg_comp_agg, FUN = sd)
  print(eeg_comp_agg_sds)
  
  # get cueing effects
  cueing_effects = aggregate(value ~ cue_modality + laterality + id, data = eeg_comp_agg, FUN = diff)
  
  if (component == "N1" | component == "N80.P100.Complex") {
    cueing_effects$value = -cueing_effects$value
  }
  
  cueing_effects_bilateral = aggregate(value ~ cue_modality + id, data = cueing_effects, FUN = mean)
  
  IOR_effects_one_target_type = IOR_effects[IOR_effects$target_modality == target_modality,]
  
  # merge to create correlation data frame 
  corr_df = merge(cueing_effects_bilateral, IOR_effects_one_target_type)
  
  corr_tactile = corr_df[corr_df$cue_modality == "tactile",]
  corr_T = cor(corr_tactile$IOR, corr_tactile$value)
  corr_T_test = cor.test(corr_tactile$IOR, corr_tactile$value)
  print(corr_T_test)
  
  corr_visual = corr_df[corr_df$cue_modality == "visual",]
  corr_V = cor(corr_visual$IOR, corr_visual$value)
  corr_V_test = cor.test(corr_visual$IOR, corr_visual$value)
  print(corr_V_test)
  
  corr_df$cue_modality = as.factor(corr_df$cue_modality)
  levels(corr_df$cue_modality) = c("Tactile", "Visual")
  
  gg3 = ggplot(data = corr_df, aes(x = IOR, y = value*1e6, color = cue_modality))+
    geom_point(size = 2)+
    geom_smooth(method = "lm", se = FALSE)+
    labs(x = "IOR Effect", y = paste("Cueing Reduction (", expression(u),"V)", sep = ""), color = "Cue Modality")+
    geom_hline(yintercept = 0, size = 1, linetype = "dashed")+
    geom_vline(xintercept = 0, size = 1, linetype = "dashed")+
    
    theme_gray(base_size = 30)+
    theme(panel.grid.major = element_line(size = 1.5)
          ,panel.grid.minor = element_line(size = 1)) 
  
  if (component == "P1") {
    gg3 = gg3 + 
      annotate("text", x = 82, y=3, label = sprintf("%f", corr_V), size = 7, color = "#00BFC4")+
      annotate("text", x = 82, y=3.5, label = sprintf("%f", corr_T), size = 7, color = "#F8766D")+
      ggtitle(component)
  } else if (component == "N80.P100.Complex") {
    gg3 = gg3 + 
      annotate("text", x = 82, y=1.5, label = sprintf("%f", corr_V), size = 7, color = "#00BFC4")+
      annotate("text", x = 82, y=2, label = sprintf("%f", corr_T), size = 7, color = "#F8766D")+
      ggtitle(component)
  }
  
  
  
  print(gg3)
  
  # data frame of summary stats 
  eeg_comp_agg_sum = cbind(eeg_comp_agg_means, eeg_comp_agg_sds$value)
  names(eeg_comp_agg_sum)[4:5] = c("M", "SD")
  
  # ANOVA
  m = aov(value ~
            cue_modality*laterality*cueing 
          + Error(id/(cue_modality*laterality*cueing))
          , data = eeg_comp_agg)
  m_summary = summary(m)
  print(m_summary)
  
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
  
  eeg_comp_agg$mix_factor = paste(eeg_comp_agg$cueing, eeg_comp_agg$cue_modality, eeg_comp_agg$laterality)
  # model
  m2 = aov(value ~ mix_factor
           + Error(id/mix_factor)
           , data = eeg_comp_agg)
  m2_summary = summary(m2)
  
  MSE = m2_summary$`Error: id:mix_factor`[1][[1]][[3]][2]
  df = m2_summary$`Error: id:mix_factor`[1][[1]][[1]][2]
  n = length(unique(eeg_comp_agg$id))
  
  SEM_LM = sqrt(MSE/n)
  
  # get critical t value and LSD
  t_crit = qt(1-alpha_over_2, df)
  ME_LM = abs(t_crit * SEM_LM)
  LSD = sqrt(2) * ME_LM
  #--------------------------------- SEM (L & M) ----------------------------------------#
  
  
  #--------------------------------- Circularity ----------------------------------------#
  # NOTE: ciruclarity = sphericity
  eeg_comp_agg$mix_factor_num = as.numeric(factor(eeg_comp_agg$mix_factor))
  pairings = combn(unique(eeg_comp_agg$mix_factor_num),2)
  
  print( cbind(unique(eeg_comp_agg[6]), unique(eeg_comp_agg[7])))
  
  alpha = 0.05
  SEs = NULL
  MEs = NULL
  ids = NULL
  MEANs = NULL
  for (i in 1:ncol(pairings)) {
    pair = pairings[,i]
    
    level1 = pair[1]
    level2 = pair[2]
    
    data1 = eeg_comp_agg[eeg_comp_agg$mix_factor_num == level1,]$value
    data2 = eeg_comp_agg[eeg_comp_agg$mix_factor_num == level2,]$value
    
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
  ezEpsi = ezANOVA(data = eeg_comp_agg
                   , dv = value
                   , wid = id
                   , within = mix_factor)
  print(ezEpsi$`Sphericity Corrections`$GGe)
  
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
  
  
  #----------------------------- Look at interaction ------------------------------------#
  eeg_comp_agg_cueing2 = aggregate(value ~ laterality + cue_modality + id,data = eeg_comp_agg, FUN = diff)
  eeg_comp_agg_cueing2$value = -eeg_comp_agg_cueing2$value  # flip sign
  eeg_comp_agg_cueing = aggregate(value ~ cue_modality + id,data = eeg_comp_agg_cueing2, FUN = mean)
  
  # calculate paired differences (from zero!)
  mods = unique(eeg_comp_agg_cueing$cue_modality)
  MEs = NULL
  for (mod in mods) {
    dat = eeg_comp_agg_cueing[eeg_comp_agg_cueing$cue_modality == mod,]$value
    
    MEAN = mean(dat)
    SD = sd(dat)
    N = length(dat)
    DF = N - 1
    T_CRIT = qt(1-alpha/2, DF)
    
    SE = SD/sqrt(N)
    
    ME = T_CRIT * SE
    MEs = c(MEs, ME)
  }
  
  eeg_comp_agg_cueing$mix_factor = paste(eeg_comp_agg_cueing$cue_modality)
  # model
  m3 = aov(value ~ mix_factor
           + Error(id/mix_factor)
           , data = eeg_comp_agg_cueing)
  m3_summary = summary(m3)
  
  MSE3 = m3_summary$`Error: id:mix_factor`[1][[1]][[3]][2]
  # which is the same as taking the MSE of the interaction from the full ANOVA
  MSE = m_summary$`Error: id:cue_modality:cueing`[1][[1]][[3]][2]
  abs(MSE - MSE3) < 0.01
  df3 = m3_summary$`Error: id:mix_factor`[1][[1]][[1]][2]
  
  SEM_LM3 = sqrt(MSE3/n)
  
  # get critical t value and LSD
  t_crit3 = qt(1-alpha_over_2, df3)
  ME_LM3 = abs(t_crit3 * SEM_LM3)
  LSD3 = sqrt(2) * ME_LM3
  
  EEG_CIs2 = aggregate(value ~ cue_modality + cueing, data = eeg_comp_agg_means, FUN = mean)
  EEG_CIs = aggregate(value ~ cue_modality, data = EEG_CIs2, FUN = diff)
  EEG_CIs$value = -EEG_CIs$value
  EEG_CIs$CI = ME_LM3  # use margin of error !
  
  names(EEG_CIs)[2] = "M"
  
  # plot 3 way interaction 
  EEG_CIs$cue_modality = as.factor(EEG_CIs$cue_modality)
  levels(EEG_CIs$cue_modality) = c("Tactile Cue", "Visual Cue")
  
  # add pair-wise
  EEG_CIs$from_zero = MEs
  
  if(substr(component, 1,1) == "P"){
    EEG_CIs$M = -EEG_CIs$M 
  }
  
  print(EEG_CIs)
  
  # plot!
  gg = ggplot(EEG_CIs, aes(x = cue_modality
                           ,y = M*1e6  # get micro volts
  )
  )+
    geom_point(size = 2)+
    ggtitle(component)+
    geom_errorbar(aes(ymin = (M - CI)*1e6, ymax = (M + CI)*1e6)
                  , width = 0.1
                  , size = 1
    )+
    
    geom_errorbar(aes(ymin = (M - from_zero)*1e6, ymax = (M + from_zero)*1e6)
                  , width = 0.1
                  , size = 1
                  , color = "red"
    )+
    labs(x = "Cue Modality",y = paste("Cueing Reduction (", expression(u),"V)", sep = ""), color = "Laterality")+
    geom_hline(yintercept = 0, size = 1, linetype = "dashed")+
    theme_gray(base_size = 30)+
    theme(panel.grid.major = element_line(size = 1.5)
          ,panel.grid.minor = element_line(size = 1)) 
  
  print(gg)
  
  
  #----------------------------- Look at interaction ------------------------------------#
  
  
  return(eeg_comp_agg_sum)
}



#### P45 ####
P45 = do_aov('P45', 'tactile', Plo, Phi)

#### N80/P100 ####
N80_P100 = do_aov('N80.P100.Complex', 'tactile', Nlo, Nhi)

#### P1 ####
P1 = do_aov("P1", 'visual', P1lo, P1hi)

#### N1 ####
N1 = do_aov("N1", 'visual', N1lo, N1hi)







