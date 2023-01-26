#************* Pairwise coinfection analyses *************#
#### 1. Importing Data and Initial Preparation of Data Frame
#### 2. Antigenic Segment Processing
#### 3. Coinfection Supernatant Titers
#### 4. Population Genetic Statistics
#### 5. Testing Strain Specificity

# Pairwise coinfection analyses Script for Influenza LMGSeq
# Script: pairwise_coinfections.R
# This file is an R script that analyzes strain assignment data from a LMGSeq experiment
# Requires running demultiplexing.sh and amplicon_curation_strain_assignment.sh to generate output files for *pairwise coinfections*:
#   ./demultiplexing.sh runA "shared/cross_list_runA_pairwise.txt" pairwise_infections
#   ./amplicon_curation_strain_assignment.sh runA "shared/cross_list_runA_pairwise.txt" pairwise_infections

#This script is part of the following manuscript:
#Influenza A virus reassortment is strain dependent
#Kishana Y. Taylor | Ilechukwu Agu  | Ivy José | Sari Mäntynen | A.J. Campbell | Courtney Mattson |  Tsui-wen Chou | Bin Zhou | David Gresham | Elodie Ghedin |  Samuel L. Díaz Muñoz


#AFTER SENDING DRAFT - Search for this below

#Load required libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape2)
library(readr)
library(gridExtra)
library(ggthemes)

# 1. Importing Data and Initial Preparation of Data Frame ----

#List output file names
file_names <- list.files("outputs/clinical", "*_strain_database17_98_locus.txt", full.names = T)

#Make a dataframe by recursively reading each output file in the list
two_strain_database17_98_locus <- NULL
for (file in file_names) {
  df <- read.table(file, quote="\"", comment.char="")
  two_strain_database17_98_locus <- rbind(two_strain_database17_98_locus, df)
}

#Assign names to columns
colnames(two_strain_database17_98_locus) <- c("cross_sample_locus", "total_reads", "majority_assigned_reads", "majority_strain_locus")

#Calculate Proportion of Reads that were assigned to one strain in the coinfection or the other
two_strain_database17_98_locus <- mutate(two_strain_database17_98_locus, proportion_assigned = majority_assigned_reads/total_reads)

#Quick quality control plot looking at relationship between total number of reads and proportion assigned to one strain 
qplot(two_strain_database17_98_locus$proportion_assigned, two_strain_database17_98_locus$total_reads)

#Cleaning up text from output files 
cross_sample_locus <- gsub("usearch/all/","", two_strain_database17_98_locus$cross_sample_locus)
two_strain_database17_98_locus$cross_sample_locus <- gsub("_98_merged.b6","", cross_sample_locus)

#Then split cross sample locus into different columns, but keep original
two_strain_database17_98_locus <- separate(two_strain_database17_98_locus, cross_sample_locus, c("cross", "sample", "locus"), sep = "_", remove=FALSE)

#Now create a cross_sample column to identify unique samples more easily
two_strain_database17_98_locus <- unite(two_strain_database17_98_locus, cross_sample, c("cross", "sample"), sep = "_", remove=FALSE)

#Now split majority strain and majority locus
two_strain_database17_98_locus <- separate(two_strain_database17_98_locus, majority_strain_locus, c("majority_strain", "majority_locus"), sep = "_", remove=FALSE)

#Check number of rows in two_strain_database17_98_locus
nrow(two_strain_database17_98_locus)
#[1] 974

#Let's try our first plot
ggplot(two_strain_database17_98_locus, aes(x = locus, y = sample)) + geom_point(aes(col=majority_strain, alpha=proportion_assigned), shape=15, size = 6) + theme(axis.text.x = element_text(angle=90)) + facet_wrap(~ cross)


#Because not all loci were informative for reassortment (identical sequences among pairs), we'll subset the dataframe to reflect only the informative loci for each cross
cross10_subset <- filter(two_strain_database17_98_locus, cross == "cross10") %>% filter(locus != "HA3a" & locus != "NS1d" & locus != "PAc")
cross10_subset <- filter(cross10_subset, sample != "sample87") #This sample was not included in this cross, just easier to remove here than in shell script
cross11_subset <- filter(two_strain_database17_98_locus, cross == "cross11") %>% filter(locus == "Mg" | locus == "NS1d")
cross12_subset <- filter(two_strain_database17_98_locus, cross == "cross12") %>% filter(locus != "HA3a" & locus != "PAc")

two_strain_database17_98_locus <- rbind(cross10_subset, cross11_subset, cross12_subset)

#Now try again
ggplot(two_strain_database17_98_locus, aes(x = locus, y = sample)) + geom_point(aes(col=majority_strain, alpha=proportion_assigned), shape=15, size = 6) + theme(axis.text.x = element_text(angle=90)) + facet_wrap(~ cross)

#Okay, now add the information on what each cross (experimental coinfection) is about. This data sheet is included in the repo. 
cross_data_clinical_samples <- read.csv("~/Dropbox/mixtup/Documentos/ucdavis/papers/influenza_GbBSeq_human/human_influenza_GbBSeq/data/cross_data_clinical_samples.csv")


#Prepare the cross_data_clinical_samples df for a merge with the big data set, so we can plot based on infection conditions
cross_data_clinical_samples <- mutate(cross_data_clinical_samples, cross = paste("cross", cross_id, sep = ""))

#Make a strain combination column for more convenient facetting, note that order will be strainA_strainB
cross_data_clinical_samples <- mutate(cross_data_clinical_samples, strainAB = paste(strainA, strainB, sep = "_"))

#Join sample information
two_strain_database17_98_locus <- left_join(two_strain_database17_98_locus, cross_data_clinical_samples, by = "cross")

#Arrange facets by strain
ggplot(two_strain_database17_98_locus, aes(x = locus, y = sample)) + geom_point(aes(col=majority_strain, alpha=proportion_assigned), shape=15, size = 6) + theme(axis.text.x = element_text(angle=90)) + facet_grid(strainA ~ strainB)

#Summarise number of reads per cross 
two_strain_database17_98_locus %>% group_by(strainAB) %>%
  summarise(
    total_cross_reads = sum(total_reads),
  )
# A tibble: 3 × 2
#strainAB                total_cross_reads
#<chr>                               <int>
#1 H3N2-NYU-19_H3N2-NYU-23            103082
#2 H3N2-NYU-19_H3N2-NYU-26              5643
#3 H3N2-NYU-23_H3N2-NYU-26              2324

#NOTE DIFFERENT SAMPLE SIZES AND LOCI IN EACH CROSS

#Summarize average number of reads per sample in each cross 
reads_sample <- two_strain_database17_98_locus %>% group_by(cross_sample, cross) %>%
  summarise(
    total_cross_reads = sum(total_reads),
  )

summary(reads_sample$total_cross_reads)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.0    99.0   353.0   810.6  1600.0  2702.0  

#Summarise number of reads per locus 
two_strain_database17_98_locus %>% group_by(strainAB, locus) %>%
  summarise(
    total_locus_reads = sum(total_reads),
    average_locus_reads = mean(total_reads),
  )

# # A tibble: 13 × 4
# Groups:   strainAB [3]
#strainAB                locus total_locus_reads average_locus_reads
#<chr>                   <chr>             <int>               <dbl>
#1 H3N2-NYU-19_H3N2-NYU-23 Mg               102464             1314.  
#2 H3N2-NYU-19_H3N2-NYU-23 NS1d                618                7.82
#3 H3N2-NYU-19_H3N2-NYU-26 Mg                  471               19.6 
#4 H3N2-NYU-19_H3N2-NYU-26 NA3a               2335               97.3 
#5 H3N2-NYU-19_H3N2-NYU-26 NPd                2700              112.  
#6 H3N2-NYU-19_H3N2-NYU-26 NS1d                116                5.27
#7 H3N2-NYU-19_H3N2-NYU-26 PB1c                  8                1   
#8 H3N2-NYU-19_H3N2-NYU-26 PB2f                 13                1.62
#9 H3N2-NYU-23_H3N2-NYU-26 Mg                  339               14.7 
#10 H3N2-NYU-23_H3N2-NYU-26 NA3a                647               28.1 
#11 H3N2-NYU-23_H3N2-NYU-26 NPd                1273               55.3 
#12 H3N2-NYU-23_H3N2-NYU-26 PB1c                  3                1   
#13 H3N2-NYU-23_H3N2-NYU-26 PB2f                 28               


#Total Rows in data frame
nrow(two_strain_database17_98_locus)
#[1] 353

#Locus by cross
locus_by_cross <- two_strain_database17_98_locus %>% group_by(strainAB, locus) %>%
  summarise(
    number_samples = length(cross_sample),
  )

#No subtype processing required, all are H3N2

# 3. Quality measures of strain assignments and calculation of reassortment frequencies

#Checking how good the assignments are
mean(two_strain_database17_98_locus$proportion_assigned)
#[1] 0.9028826
sd(two_strain_database17_98_locus$proportion_assigned)
#[1] 0.1394683

#Pretty good overall
#How many sample/locus combinations are over a certain proportion 
nrow(subset(two_strain_database17_98_locus, proportion_assigned > .80)) / nrow(two_strain_database17_98_locus)
#0.7903683
nrow(subset(two_strain_database17_98_locus, proportion_assigned > .75)) / nrow(two_strain_database17_98_locus)
#0.8470255
nrow(subset(two_strain_database17_98_locus, proportion_assigned > .70)) / nrow(two_strain_database17_98_locus)
#0.8895184
nrow(subset(two_strain_database17_98_locus, proportion_assigned > .65)) / nrow(two_strain_database17_98_locus)
#0.9235127

#Histogram of the proportion of reads assigned to one strain (note this is expected to be close to 1, because plaques should be clonal but sequencing errors particularly when read count is small can lead to lower proportions)
ggplot(two_strain_database17_98_locus, aes(x = proportion_assigned)) + geom_histogram()

#Same reassortment plot, but y scales free because different sample sizes
ggplot(two_strain_database17_98_locus, aes(x = locus, y = sample)) + geom_point(aes(col=majority_strain, alpha=proportion_assigned), shape=15, size = 6) + theme(axis.text.x = element_text(angle=90)) + facet_wrap(~ cross, scales = "free_y")

#Rename data frame
controls_df <- two_strain_database17_98_locus

#First need to get the number of typed loci for each sample
controls_df_loci_numbers <- controls_df %>% 
  group_by(cross_sample) %>% 
  summarise(
    total_segments = length(locus),
    max_majority_strain = max(table(majority_strain))
  )

controls_df <- right_join(controls_df, controls_df_loci_numbers)

#Once we have numbers, sort by parentals and plot
ggplot(controls_df, aes(x = locus, y = reorder(sample, max_majority_strain))) + geom_point(aes(col=majority_strain, alpha=proportion_assigned), shape=15, size = 6) + theme(axis.text.x = element_text(angle=90)) + facet_wrap(~ cross, scales = "free_y")

#Once we have numbers, sort by parentals and pairwise grid and plot
ggplot(controls_df, aes(x = locus, y = sample)) + geom_point(aes(col=majority_strain, alpha=proportion_assigned), shape=15, size = 6) + theme(axis.text.x = element_text(angle=90), axis.text.y = element_blank()) + facet_grid(strainA ~ strainB, scales = "free_y")

#Pairwise grid
ggplot(controls_df, aes(x = locus, y = sample)) + geom_point(aes(col=majority_strain, alpha=proportion_assigned), shape=15, size = 6) + theme(axis.text.x = element_text(angle=90), axis.text.y = element_blank()) + facet_grid(strainA ~ strainB, scales = "free_y")

#Not ordered by parentals/reassortants
ggplot(controls_df, aes(x = locus, y = sample)) + geom_point(aes(col=majority_strain, alpha=proportion_assigned), shape=15, size = 6) + theme(axis.text.x = element_text(angle=90)) + facet_wrap(~ cross)

#Not ordered
ggplot(controls_df, aes(x = locus, y = sample)) + geom_point(aes(col=majority_strain, alpha=proportion_assigned), shape=15, size = 6) + theme(axis.text.x = element_text(angle=90)) + theme(axis.text.y = element_blank()) + facet_grid(strainA ~ strainB, scales = "free_y")


#Figure 2A - Not ordered, all samples genotypes, order segments in proper influenza order
ggplot(controls_df, aes(x = factor(locus, level = c('PB2f', 'PB1c', 'PAc', 'HA', 'NPd', 'NA', 'Mg', 'NS1d')), y = sample)) + geom_point(aes(col=majority_strain), shape=15, size = 6) + theme(axis.text.x = element_text(angle=90)) + xlab("Segment") + ylab("Plaque Isolate") + theme(axis.text.y = element_blank()) + facet_grid(strainA ~ strainB, scales = "free_y")
#Let's spruce up re: Elodie comments, this is the final figure 2A
ggplot(controls_df, aes(x = factor(locus, level = c('PB2f', 'PB1c', 'PAc', 'HA', 'NPd', 'NA', 'Mg', 'NS1d')), y = sample)) + geom_point(aes(col=majority_strain), shape=15, size = 6) + xlab("Segment") + ylab("Plaque Isolate") + facet_grid(strainA ~ strainB, scales = "free_y") + theme_tufte() +  theme(text = element_text(size = 20,  family="Helvetica")) + theme(axis.text.x = element_text(angle=90, size = 17)) + theme(axis.text.y = element_blank()) + theme(strip.text = element_text(face = "bold")) + theme(legend.position = "none")      
#ggsave("outputs/figure2A.pdf", width = 8, height = 11)

#Now let's generate some numbers!
parentals <- controls_df %>% 
  group_by(cross_sample, cross, sample) %>% 
  summarise(
    parental = max_majority_strain == total_segments,
    #total_segments == max_majority_strain
  )

parentals <- distinct(parentals)

#This table gives a quick look at the number of samples per cross that are parentals (i.e. not reassortant)
table(parentals$cross, parentals$parental)

#Group by cross_sample and calculate number of parents for each
controls_df_number_parents <- controls_df %>% 
  group_by(cross, cross_sample) %>% 
  summarise(
    number_parents = length(unique(majority_strain)),
    #constellation = paste(majority_strain_locus, sep = ","),
    max_majority_strain = max(table(majority_strain))
  )

nrow(subset(controls_df_number_parents, number_parents == 2))
#92
nrow(controls_df_number_parents)
#[1] 136

#Make a reasortant column 
controls_df_number_parents <- mutate(controls_df_number_parents, reassortant = ifelse(number_parents > 1, yes = 1, no = 0))

#Make a data frame that has per-cross (i.e. experimental coinfection) statistics
cross_stats <- controls_df_number_parents %>% 
  group_by(cross) %>% 
  summarise(
    clones = length(number_parents),
    reassortants = sum(reassortant)
    #constellation = paste(majority_strain_locus, sep = ","),
    #max_majority_strain = max(table(majority_strain))
  )

#Add proportion reassortant to data frame that has per-cross stats
cross_stats <- mutate(cross_stats, prop_reassortant = reassortants / clones)

#Going forward with allowing incomplete genotypes, because this is conservative with respect to reassortemnt, because missing segments decrease the chance of detecting a reassortant
ggplot(controls_df, aes(x = locus, y = reorder(sample, max_majority_strain))) + geom_point(aes(col=majority_strain, alpha=proportion_assigned), shape=15, size = 6) + theme(axis.text.x = element_text(angle=90)) + facet_wrap(~ cross)

#Plot reasortment proportions
ggplot(cross_stats, aes(x = cross, y = prop_reassortant)) + geom_col()

#Same but reorder by height
ggplot(cross_stats, aes(x = reorder(cross, prop_reassortant), y = prop_reassortant)) + geom_col()

#Same but reorder by height with line for 40% reassortment but then also for the theoretical free reassortment
ggplot(cross_stats, aes(x = reorder(cross, prop_reassortant), y = prop_reassortant)) + geom_col() + geom_hline(yintercept = 0.40, linetype = 2, color = "grey", alpha = 0.75) + geom_hline(yintercept = 0.9921875, linetype = 2, color = "red", alpha = 0.75)

#Add strainAB information to cross_stats DF
cross_stats <- left_join(cross_stats, cross_data_clinical_samples, by = "cross")

#Same but reorder by height with line for 40% reassortment but then also for the theoretical free reassortment and add strain info on x axis 
ggplot(cross_stats, aes(x = reorder(strainAB, prop_reassortant), y = prop_reassortant)) + geom_col() + ylab("Proportion of Reassortant Plaque Isolates") + xlab("Strains in Experimental Coinfection") + geom_hline(yintercept = 0.40, linetype = 2, color = "grey", alpha = 0.75) + geom_hline(yintercept = 0.9921875, linetype = 2, color = "red", alpha = 0.75)
#Update for better publication quality. This is the final Figure 2B
ggplot(cross_stats, aes(x = reorder(strainAB, prop_reassortant), y = prop_reassortant)) + geom_col() + ylab("Proportion of Reassortant Plaque Isolates") + xlab("Strains in Experimental Coinfection") + geom_hline(yintercept = 0.40, linetype = 2, color = "grey", alpha = 0.75) + geom_hline(yintercept = 0.9921875, linetype = 2, color = "red", alpha = 0.75) + theme_tufte() +  theme(text = element_text(size = 20,  family="Helvetica")) + theme(axis.text.x = element_text(size = 12)) + theme(legend.position = "none") + theme(panel.grid.major.y = element_line(color = "lightgray",size = 0.5)) 
#ggsave("outputs/figure2B.pdf", width = 12, height = 8.5)

#Now more stats across all crosses 
mean(cross_stats$prop_reassortant)
#[1] 0.8244586

sd(cross_stats$prop_reassortant)
#[1] 0.267279

min(cross_stats$prop_reassortant)
#[1] 0.5168539

max(cross_stats$prop_reassortant)
#[1] 1

#Let's calculate some stats on the expected vs observed number of reasortants

#Reassortment as success. Perfect reassortment expectation is calculated as 254/256, i.e. 2^8 = 256 unique combinations, two of those subtracted becausecause they are parentals
#So first let's do an example "manually" with H3N2-NYU-19_H3N2-NYU-26/cross12
binom.test(c(24, 0), p = 254/256)
#	Exact binomial test
#data:  c(24, 0)
#number of successes = 24, number of trials = 24, p-value = 1
#alternative hypothesis: true probability of success is not equal to 0.9921875
#95 percent confidence interval:
#  0.8575264 1.0000000
#sample estimates:
#  probability of success 
#1

#Then H3N2-NYU-23_H3N2-NYU-26/cross10
binom.test(c(22, 1), p = 254/256)
#Exact binomial test

#data:  c(22, 1)
#number of successes = 22, number of trials = 23, p-value = 0.1651
#alternative hypothesis: true probability of success is not equal to 0.9921875
#95 percent confidence interval:
#  0.7805134 0.9988998
#sample estimates:
#  probability of success 
#0.9565217 

#Then H3N2-NYU-19_H3N2-NYU-23/cross11
binom.test(c(46, 43), p = 254/256)
#Exact binomial test

#data:  c(46, 43)
#number of successes = 46, number of trials = 89, p-value < 2.2e-16
#alternative hypothesis: true probability of success is not equal to 0.9921875
#95 percent confidence interval:
#  0.4084187 0.6241361
#sample estimates:
#  probability of success 
#0.5168539 

#In highest two reassortment experimental coinfections they were not different from random

controls_df <- right_join(controls_df, cross_stats)

#### Quality control effects on strain assignments ####
#We'll make a new DF called controls_df_qc and add modifications sequentially and check on stats

#First let's take a close look at the proportion of reads assigned to one strain or the other.
# To recap overall the assignments are pretty good, overall if a bit lower than the main data set

nrow(subset(two_strain_database17_98_locus, proportion_assigned > .80)) / nrow(two_strain_database17_98_locus)
#0.7903683

ggplot(two_strain_database17_98_locus, aes(x = proportion_assigned)) + geom_histogram()

nrow(two_strain_database17_98_locus)
#[1] 353

#Clearly should get rid of strain assignments of segments based on a single read, even though LMGSeq is designed to work under this scenario, 
# but because PID between the strains in these loci are very high excercising extra caution
controls_df_qc <- filter(two_strain_database17_98_locus, majority_assigned_reads > 1)
nrow(controls_df_qc)
#[1] 316

#Now getting rid of assignments that are based on less than a 66.667% (2/3) of reads when sample size is under 50 (that would give p < 0.05 in an exact binomial)  
remove_low_qual <- filter(controls_df_qc, proportion_assigned < 0.67) %>% filter(total_reads < 50)
controls_df_qc <- anti_join(controls_df_qc, remove_low_qual)
nrow(controls_df_qc)
#[1] 293

#Now let's look at the proportion stats and plots
nrow(subset(controls_df_qc, proportion_assigned > .80)) / nrow(controls_df_qc)
#[1] 0.8464164
ggplot(controls_df_qc, aes(x = proportion_assigned)) + geom_histogram()

#Much better

#Now let's see the effect on reassortment estimates
#First, original plot
ggplot(two_strain_database17_98_locus, aes(x = locus, y = sample)) + geom_point(aes(col=majority_strain, alpha=proportion_assigned), shape=15, size = 6) + theme(axis.text.x = element_text(angle=90)) + facet_wrap(~ cross, scales = "free")
#Second corrected plot
ggplot(controls_df_qc, aes(x = locus, y = sample)) + geom_point(aes(col=majority_strain, alpha=proportion_assigned), shape=15, size = 6) + theme(axis.text.x = element_text(angle=90)) + facet_wrap(~ cross, scales = "free")
#Looks like that basically took out the PB1c locus so now cross 10 and 12 are based on 4 loci (though not complete)

##### Now calculate reasortment stats ####
#Now let's generate some numbers!
#First need to get the number of typed loci for each sample
controls_df_loci_numbers_qc <- controls_df_qc %>% 
  group_by(cross_sample) %>% 
  summarise(
    total_segments = length(locus),
    max_majority_strain = max(table(majority_strain))
  )

controls_df_qc <- right_join(controls_df_qc, controls_df_loci_numbers_qc)

parentals_qc <- controls_df_qc %>% 
  group_by(cross_sample, cross, sample) %>% 
  summarise(
    parental = max_majority_strain == total_segments,
    #total_segments == max_majority_strain
  )

parentals_qc <- distinct(parentals_qc)

#This table gives a quick look at the number of samples per cross that are parentals (i.e. not reassortant)
table(parentals_qc$cross, parentals_qc$parental)
#FALSE TRUE
#cross10    11   12
#cross11    42   41
#cross12    24    0

#So cross 10 took a hit on sample size and reassortants 

#Group by cross_sample and calculate number of parents for each
controls_df_number_parents_qc <- controls_df_qc %>% 
  group_by(cross, cross_sample) %>% 
  summarise(
    number_parents = length(unique(majority_strain)),
    #constellation = paste(majority_strain_locus, sep = ","),
    max_majority_strain = max(table(majority_strain))
  )

nrow(subset(controls_df_number_parents_qc, number_parents == 2))
#77
nrow(controls_df_number_parents_qc)
#[1] 130

#Make a reasortant column 
controls_df_number_parents_qc <- mutate(controls_df_number_parents_qc, reassortant = ifelse(number_parents > 1, yes = 1, no = 0))

#Make a data frame that has per-cross (i.e. experimental coinfection) statistics
cross_stats_qc <- controls_df_number_parents_qc %>% 
  group_by(cross) %>% 
  summarise(
    clones = length(number_parents),
    reassortants = sum(reassortant)
    #constellation = paste(majority_strain_locus, sep = ","),
    #max_majority_strain = max(table(majority_strain))
  )

#Add proportion reassortant to data frame that has per-cross stats
cross_stats_qc <- mutate(cross_stats_qc, prop_reassortant = reassortants / clones)

#Going forward with allowing incomplete genotypes, because this is conservative with respect to reassortemnt, because missing segments decrease the chance of detecting a reassortant
ggplot(controls_df_qc, aes(x = locus, y = reorder(sample, max_majority_strain))) + geom_point(aes(col=majority_strain, alpha=proportion_assigned), shape=15, size = 6) + theme(axis.text.x = element_text(angle=90)) + facet_wrap(~ cross, scales = "free")

#Plot reasortment proportions
ggplot(cross_stats_qc, aes(x = cross, y = prop_reassortant)) + geom_col()

#Same but reorder by height
ggplot(cross_stats_qc, aes(x = reorder(cross, prop_reassortant), y = prop_reassortant)) + geom_col()

#Same but reorder by height with line for 40% reassortment but then also for the theoretical free reassortment
ggplot(cross_stats_qc, aes(x = reorder(cross, prop_reassortant), y = prop_reassortant)) + geom_col() + geom_hline(yintercept = 0.40, linetype = 2, color = "grey", alpha = 0.75) + geom_hline(yintercept = 0.9921875, linetype = 2, color = "red", alpha = 0.75)

#Add strainAB information to cross_stats DF
cross_stats_qc <- left_join(cross_stats_qc, cross_data_clinical_samples, by = "cross")

#Same but reorder by height with line for 40% reassortment but then also for the theoretical free reassortment and add strain info on x axis 
ggplot(cross_stats_qc, aes(x = reorder(strainAB, prop_reassortant), y = prop_reassortant)) + geom_col() + ylab("Proportion of Reassortant Plaque Isolates") + xlab("Strains in Experimental Coinfection") + geom_hline(yintercept = 0.40, linetype = 2, color = "grey", alpha = 0.75) + geom_hline(yintercept = 0.9921875, linetype = 2, color = "red", alpha = 0.75)
#Update for better publication quality. This is the final Figure 2B
ggplot(cross_stats_qc, aes(x = reorder(strainAB, prop_reassortant), y = prop_reassortant)) + geom_col() + ylab("Proportion of Reassortant Plaque Isolates") + xlab("Strains in Experimental Coinfection") + geom_hline(yintercept = 0.40, linetype = 2, color = "grey", alpha = 0.75) + geom_hline(yintercept = 0.9921875, linetype = 2, color = "red", alpha = 0.75) + theme_tufte() +  theme(text = element_text(size = 20,  family="Helvetica")) + theme(axis.text.x = element_text(size = 12)) + theme(legend.position = "none") + theme(panel.grid.major.y = element_line(color = "lightgray",size = 0.5)) 
#ggsave("outputs/figure2B.pdf", width = 12, height = 8.5)

#Now more stats across all crosses 
mean(cross_stats_qc$prop_reassortant)
#[1] 0.6614283

sd(cross_stats_qc$prop_reassortant)
#[1] 0.2935401

min(cross_stats_qc$prop_reassortant)
#[1] 0.4782609

max(cross_stats_qc$prop_reassortant)
#[1] 1

#Let's calculate some stats on the expected vs observed number of reasortants

#Reassortment as success. Perfect reassortment expectation is calculated as 254/256, i.e. 2^8 = 256 unique combinations, two of those subtracted becausecause they are parentals
#So first let's do an example "manually" with H3N2-NYU-19_H3N2-NYU-26/cross12
binom.test(c(24, 0), p = 254/256)
#	Exact binomial test
#data:  c(24, 0)
#number of successes = 24, number of trials = 24, p-value = 1
#alternative hypothesis: true probability of success is not equal to 0.9921875
#95 percent confidence interval:
#  0.8575264 1.0000000
#sample estimates:
#  probability of success 
#1

#Then H3N2-NYU-23_H3N2-NYU-26/cross10
binom.test(c(11, 12), p = 254/256)
#Exact binomial test

#	Exact binomial test

#data:  c(11, 12)
#number of successes = 11, number of trials = 23, p-value < 2.2e-16
#alternative hypothesis: true probability of success is not equal to 0.9921875
#95 percent confidence interval:
#  0.2681962 0.6941220
#sample estimates:
#  probability of success 
#0.4782609 

#Then H3N2-NYU-19_H3N2-NYU-23/cross11
binom.test(c(42, 41), p = 254/256)
#Exact binomial test

#	Exact binomial test

#data:  c(42, 41)
#number of successes = 42, number of trials = 83, p-value < 2.2e-16
#alternative hypothesis: true probability of success is not equal to 0.9921875
#95 percent confidence interval:
#  0.3939794 0.6176265
#sample estimates:
#  probability of success 
#0.5060241 

#After QC only cross12 is not different from random

#This is final figure on supplementary materials for supplement
ggplot(subset(controls_df_qc, cross != "cross11"), aes(x = factor(locus, level = c('PB2f', 'PB1c', 'PAc', 'HA', 'NPd', 'NA', 'Mg', 'NS1d')), y = sample)) + geom_point(aes(col=majority_strain), shape=15, size = 6) + xlab("Segment") + ylab("Plaque Isolate") + facet_wrap( ~ strainAB, scales = "free") + theme(axis.text.x = element_text(angle=90)) + theme(axis.text.y = element_blank()) + theme(strip.text = element_text(face = "bold")) + theme(legend.position = "none")      
ggsave("outputs/figureS5.pdf", width = 5, height = 6)
