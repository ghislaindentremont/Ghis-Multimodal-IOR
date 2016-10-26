library(ggplot2)
library(plyr)
library(reshape)
library(scales)
library(ggthemes)

setwd("/Users/ghislaindentremont/Documents/Multimodal_IOR/Ghis/P_analysis2")


##################################################################
####                  Grand Averages                          ####
##################################################################

b_list2 = list.files(
  path = "."
  , pattern = "grand_averages"
  , recursive = T
)
b_list = b_list2[grep(".csv", b_list2)]

load("/Users/ghislaindentremont/Documents/Multimodal_IOR/Ghis/TXT/P_list.Rdata")
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
    scale_x_continuous("Time (ms)")+
    scale_y_reverse(paste("Voltage (", expression(u), "V)", sep = ""))+
    theme_gray(base_size = 30)+
    theme(panel.grid.major = element_line(size = 1.5)
          ,panel.grid.minor = element_line(size = 1)) 
  
  return(gg)
}

# tactile
tactile_only = b2_agg[b2_agg$target_modality == 'Tactile Target\n(C3/4)',]
plot_grands(tactile_only, 25, 60, 60, 120)


# visual
visual_only = b2_agg[b2_agg$target_modality == 'Visual Target\n(PO7/8)',]
plot_grands(visual_only, 90, 175, 175, 205)



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
  facet_grid(cue_modality ~ target_modality + laterality)+
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




# ##################################################################
# ####                     Analysis                             ####
# ##################################################################
# 
# do_aov = function(component, target_modality, lower_bound, upper_bound) {
#   # P45
#   eeg_comp = a2[a2$target_modality == target_modality
#                           & a2$time >= lower_bound
#                           & a2$time <= upper_bound, ]
#   
#   # get average voltage in window for each condition (4) x cueing (2) x participant (n)
#   eeg_comp_agg = aggregate(value ~ cue_modality + laterality + cueing + id, data = eeg_comp, FUN = mean)
#   
#   # get means 
#   eeg_comp_agg_means = aggregate(value ~ cue_modality + laterality + cueing, data = eeg_comp_agg, FUN = mean)
#   
#   # get sds
#   eeg_comp_agg_sds = aggregate(value ~ cueing + cue_modality + laterality + cueing, data = eeg_comp_agg, FUN = sd)
#   
#   # data frame of summary stats 
#   eeg_comp_agg_sum = cbind(eeg_comp_agg_means, eeg_comp_agg_sds$value)
#   names(eeg_comp_agg_sum)[4:5] = c("M", "SD")
#   
#   
#   # Run ANOVA
#   m <- aov(value ~
#               cue_modality*laterality*cueing
#             + Error(id/(cue_modality*laterality*cueing))
#             , data = eeg_comp_agg)
#   print(summary(m))
#   
# #   # ONLY TEST NORMALITY OF RESIDUALS (CANNOT TEST SPHERICITY WITH TWO LEVELS)
# #   residz = NULL
# #   for (i in 1:7) {
# #     m_stuff = proj(m)
# #     temp = m_stuff[[2+i]][,"Residuals"]
# #     temp = temp[seq(1,184,8)]
# #     residz = c(residz, temp)
# #   }
# #   
# #   qqnorm(residz)
# #   qqline(residz)
#   
#   
#   # let's generate CIs
#   # equaltion: LSD = t(alpha/2, N-a) * sqrt(2*MSE/n), when n1=n2
#   alpha_over_2 = 0.025
#   
#   # GET MSE from one-way ANOVA
#   EEG_cueing = aggregate(value ~ cueing + id, data = eeg_comp_agg, FUN = mean) 
#   # model
#   m2 = aov(value ~ cueing
#            + Error(id/cueing)
#            , data = EEG_cueing)
#   m2_summary = summary(m2)
#   
#   MSE = m2_summary$`Error: id:cueing`[1][[1]][[3]][2]
#   df = m2_summary$`Error: id:cueing`[1][[1]][[1]][2]
#   n = df + 1
#   # double check 
#   n == length(unique(eeg_comp_agg$id))
#   
#   # get critical t value and LSD
#   t_crit = qt(alpha_over_2, df)
#   LSD = abs(t_crit * sqrt(2*MSE/n))
#   
#   # now identify your IOR effects 
#   EEG_effects = aggregate(value ~ cue_modality + laterality + id ,data = eeg_comp_agg, FUN = diff)
#   EEG_CIs = aggregate(value ~ cue_modality + laterality, data = EEG_effects, FUN = mean)
#   # change name
#   names(EEG_CIs)[3] = "M"
#   
#   # plot 3 way 
#   # add ci to data frame
#   EEG_CIs$CI = LSD
#   
#   # plot 3 way interaction 
#   EEG_CIs$cue_modality = as.factor(EEG_CIs$cue_modality)
#   levels(EEG_CIs$cue_modality) = c("Tactile Cue", "Visual Cue")
#   EEG_CIs$laterality = as.factor(EEG_CIs$laterality)
#   levels(EEG_CIs$laterality) = c("Contralateral", "Ipsilateral")
#   
#   # plot!
#   gg = ggplot(EEG_CIs, aes(x = cue_modality
#                            ,y = M*1e6  # get micro volts 
#                            , group = laterality
#                            , fill = laterality
#                            , color = laterality
#   )
#   )+
#     geom_line(size = 1)+
#     geom_errorbar(aes(ymin = (M - CI)*1e6, ymax = (M + CI)*1e6)
#                   , width = 0.1
#                   , size = 1
#     )+
#     labs(x = "Cue Modality",y = paste("Cueing Effect: Cued - Uncued ", expression(u),"V", sep = ""), color = "Laterality")+
#     geom_hline(yintercept = 0, size = 1, linetype = "dashed")+
#     theme_gray(base_size = 30)+
#     theme(panel.grid.major = element_line(size = 1.5)
#           ,panel.grid.minor = element_line(size = 1)) 
#   
#   print(gg)
#   
#   
#   return(eeg_comp_agg_sum)
# }
# 
# 
# 
# #### P45 ####
# P45 = do_aov('P45', 'tactile', 30, 65)
# 
# #### N80/P100 ####
# N80_P100 = do_aov('N80.P100.Complex', 'tactile', 65, 120)
# 
# #### P1 ####
# P1 = do_aov("P1", 'visual', 110, 170)
# 
# #### N1 ####
# N1 = do_aov("N1", 'visual', 170, 205)
# 
# 
# 
