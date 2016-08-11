# libraries
library(grid)
library(gridExtra)
library(plyr)
library(ggplot2) 
library(stringr)

# directory
filedir = readline("Where is the txt file found? ")

# participant
participant =readline("What is the participant id? ")

# where to save info to later
savefir = readline("Where should I save graphs? ")

file_ls = list.files(
    path = filedir
    , pattern = participant
    , full.names = T 
  )

filename = file_ls[ grep(".txt", dat_ls) ]

a = read.table(filename, header = T)

print("getting rid of practice trials...")
b = a[a$block != "practice",]

summarize_b = summary(b)

print(summarize_b)

print("Proportion of blinks: ")
sum(b$critical_blink)/length(b$critical_blink)
print("Proportion of saccades: ")
sum(b$critical_saccade)/length(b$critical_saccade)
print("Proportion of pre target responses: ")
sum(b$pre_target_response)/length(b$pre_target_response)

hist(b$target_response_rt, breaks = 50)
abline(v = 100)

# how many left in total?
b$keep = TRUE
b[b$pre_target_response,]$keep = FALSE
b[b$critical_saccade,]$keep = FALSE
b[b$critical_blink,]$keep = FALSE
b[b$target_type == "catch",]$keep = FALSE  # can't use these for analysis
b[is.na(b$target_response_rt),]$keep = FALSE
b[!is.na(b$target_response_rt) & b$target_response_rt < 100,]$keep = FALSE


print("Number of trials we began with: ")
print( nrow(b) )
print("Number of usable trials left over after all exclusions: ")
print( sum(b$keep) )
print("Proportion of trials we KEPT: ")
print(sum(b$keep)/nrow(b))

