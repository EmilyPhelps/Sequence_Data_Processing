#This is an R script to calculate the read_depth of resequencing samples accross the genome. 
#Author= Emily Phelps
#Last Update= 29/8/2022

#Packages

library(tidyverse)

#Fixed variables
#Edit DIR, BAM_LIST and OUT

DIR <- " "
OUT <- " "
setwd(DIR)
BAM_LIST <- read.table("depth.list", header=FALSE)
BAM_ARR <- BAM_LIST$V1

for (i in 1:length(BAM_ARR)){
  # Compute depth stats
  bamfile = BAM_ARR[i]
  depth <- read_tsv(paste0(BAM_ARR[i]), col_names= F)$X1
  mean_depth <- mean(depth)
  sd_depth <- sd(depth)
  mean_depth_nonzero <- mean(depth[depth > 0])
  mean_depth_within2sd <- mean(depth[depth < mean_depth + 2 * sd_depth])
  median <- median(depth)
  presence <- as.logical(depth)
  proportion_of_reference_covered <- mean(presence)
  
  
  # Bind stats into dataframe and store sample-specific per base depth and presence data
  if (i==1){
    output <- data.frame(bamfile, mean_depth, sd_depth, mean_depth_nonzero, mean_depth_within2sd, median, proportion_of_reference_covered)
    total_depth <- depth
    total_presence <- presence
	} 

else {
    output <- rbind(output, cbind(bamfile, mean_depth, sd_depth, mean_depth_nonzero, mean_depth_within2sd, median, proportion_of_reference_covered))
    total_depth <- total_depth + depth
    total_presence <- total_presence + presence
  }
}
print(output) 
output %>%
  mutate(across(where(is.numeric), round, 3))

write.table(output, file=paste(OUT, ".tab", sep="\t"), row.names=FALSE, quote=FALSE)
  
# Plot the depth distribution
depth_dist <- 
  tibble(total_depth = total_depth, position = 1:length(total_depth))  %>%
  ggplot(aes(x = position, y = total_depth)) +
  geom_point(size = 0.1)

#Total depth per site accross all individuals
total_depth_summary <- count(tibble(total_depth = total_depth), total_depth)
total_presence_summary <- count(tibble(total_presence = total_presence), total_presence)

total_depth <- 
  total_depth_summary %>%
  ggplot(aes(x = log(total_depth), y = n)) +
  geom_point()
  
total_presence <-
  total_presence_summary %>%
  ggplot(aes(x = total_presence, y = n)) +
  geom_col()

pdf(file=paste(OUT,"_depth_dist.pdf", sep=""))
    depth_dist 
dev.off()
  
pdf(file=paste(OUT,"_total_depth.pdf", sep=""))
    total_depth
dev.off()

pdf(file=paste(OUT,"_depth_presence.pdf", sep=""))  
  total_presence
dev.off()
