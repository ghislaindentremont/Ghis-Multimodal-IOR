print("
R script that read in and summarizes participant behavioural data.
Get RT distribution, and exclusion proportions, and overall usable trial count.
")

# libraries
library(grid)
library(gridExtra)
library(plyr)
library(ggplot2) 
library(stringr)

# directory
filedir = readline("Where is the txt file found? ")
# filedir = "/Users/ray/Experiments/Ghis-Multimodal-IOR/_Data"

# participant
participant =readline("What is the participant id? ")

# where to save info to later
# savedir = readline("Where should I save results? ")
# savedir = "/Users/ray/Desktop/Multimodal Quick Results"

file_ls = list.files(
    path = filedir
    , pattern = participant
    , full.names = T 
    , recursive = T
  )

filename = file_ls[ grep(".txt", file_ls) ]

a = read.table(filename, header = T)

print("getting rid of practice trials...")
b = a[a$block != "practice",]

summarize_b = summary(b)

print(summarize_b)

hist(b$target_response_rt, breaks = 50)
abline(v = 100)
# dev.copy(png, sprintf("%s/%s.png",savedir,participant))
# dev.off()

print("Proportion of blinks: ")
blink_prop = sum(b$critical_blink)/length(b$critical_blink)
print(blink_prop)

print("Proportion of saccades: ")
sacc_prop = sum(b$critical_saccade)/length(b$critical_saccade)
print(sacc_prop)

print("Proportion of pre target responses: ")
pre_prop = sum(b$pre_target_response)/length(b$pre_target_response)
print(pre_prop)

# how many left in total?
b$keep = TRUE
b[b$pre_target_response,]$keep = FALSE
b[b$critical_saccade,]$keep = FALSE
b[b$critical_blink,]$keep = FALSE
b[b$target_type == "catch",]$keep = FALSE  # can't use these for analysis
b[is.na(b$target_response_rt),]$keep = FALSE
b[b$target_response_rt < 100 & !is.na(b$target_response_rt),]$keep = FALSE


print("Number of trials we began with: ")
init_num = nrow(b)
print(init_num)

print("Number of usable trials left over after all exclusions: ")
usable_num = sum(b$keep)
print(usable_num)

print("Proportion of trials we KEPT: ")
kept_prop_num = sum(b$keep)/nrow(b)
print(kept_prop_num)

df = data.frame(blink_prop
	, sacc_prop
	, pre_prop
	, init_num
	, usable_num
	, kept_prop_num
	)

# write.csv(df, file = sprintf("%s/%s.csv",savedir,participant))

