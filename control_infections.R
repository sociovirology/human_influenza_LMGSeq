################ Control Experiment Infections  ################
#### 1. Data Sources and Preparing Data
#### 2. Antigenic Segment Processing
#### 3. Basic Reassortment Plots
#### 4. Calculation of Reassortment Rates, using 8 segment
#### 5. Figures Plotting Reasortment Proportions and Replicate Deviation Data 

# Pairwise coinfection analyses Script for Influenza LMGSeq
# Script: pairwise_coinfections.R
# This file is an R script that conducts analyses on experiments to conduct control infections to test the LMGSeq approach
# Requires running demultiplexing.sh and amplicon_curation_strain_assignment.sh to generate output files for *control coinfections*:
#   ./demultiplexing.sh runA "shared/cross_list_runA_control.txt" control_infections
#   ./amplicon_curation_strain_assignment.sh runA "shared/cross_list_runA_control.txt" control_infections

#This script is part of the following manuscript:
#Influenza A virus reassortment is strain dependent
#Kishana Y. Taylor | Ilechukwu Agu  | Ivy José | Sari Mäntynen | A.J. Campbell | Courtney Mattson |  Tsui-wen Chou | Bin Zhou | David Gresham | Elodie Ghedin |  Samuel L. Díaz Muñoz

## Experiments
#We conducted several control experiments:
#1. Biological replicates of HK68/CA09, HK68/PAN99, and CA09/PAN99 
#-HK68/CA09: Run A Cross 2 and 4
#-HK68/PAN99: Run A Cross 3 and 5
#-CA09/PAN99: Run A Cross 11, 13, and 17

#2. Varying MOI of HK/CA09
#-Run A Cross 2 MOI = 10
#-Run A Cross 10 MOI = 0.01
#-Run A Cross 6 MOI = 1


#Load libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape2)
library(readr)

#### 1. Data Sources and Preparing Data ####

#Read files
file_names <- list.files("outputs/control_infections", "*_strain_database17_98_locus.txt", full.names = T)
#file_names <- file_names[c(3,4)]

#Going to add 

two_strain_database17_98_locus <- NULL
for (file in file_names) {
  df <- read.table(file, quote="\"", comment.char="")
  two_strain_database17_98_locus <- rbind(two_strain_database17_98_locus, df)
}

#Assign names to columns
colnames(two_strain_database17_98_locus) <- c("cross_sample_locus", "total_reads", "majority_assigned_reads", "majority_strain_locus")

#Calculate Proportions
two_strain_database17_98_locus <- mutate(two_strain_database17_98_locus, proportion_assigned = majority_assigned_reads/total_reads)

#Quick quality control plot
qplot(two_strain_database17_98_locus$proportion_assigned, two_strain_database17_98_locus$total_reads)

cross_sample_locus <- gsub("usearch/all/","", two_strain_database17_98_locus$cross_sample_locus)

two_strain_database17_98_locus$cross_sample_locus <- gsub("_98_merged.b6","", cross_sample_locus)

#Then split cross sample locus
two_strain_database17_98_locus <- separate(two_strain_database17_98_locus, cross_sample_locus, c("cross", "sample", "locus"), sep = "_", remove=FALSE)

#Now create a cross_sample column to identify unique samples more easily
two_strain_database17_98_locus <- unite(two_strain_database17_98_locus, cross_sample, c("cross", "sample"), sep = "_", remove=FALSE)

#Now split majority strain and majority locus
two_strain_database17_98_locus <- separate(two_strain_database17_98_locus, majority_strain_locus, c("majority_strain", "majority_locus"), sep = "_", remove=FALSE)

#Check number of rows in two_strain_database17_98_locus
nrow(two_strain_database17_98_locus)
#[1] 6503

#Let's try our first plot
ggplot(two_strain_database17_98_locus, aes(x = locus, y = sample)) + geom_point(aes(col=majority_strain, alpha=proportion_assigned), shape=15, size = 6) + theme(axis.text.x = element_text(angle=90)) + facet_wrap(~ cross)

no_subtype_processing_plot <- ggplot(two_strain_database17_98_locus, aes(x = locus, y = sample)) + geom_point(aes(col=majority_strain, alpha=proportion_assigned), shape=15, size = 6) + theme(axis.text.x = element_text(angle=90)) + facet_wrap(~ cross)

#Okay, now add the information on what each cross is about
cross_data_runA <- read.csv("data/cross_data_runA.csv")

#Prepare the cross_data_runA df for a merge with the big data set, so we can plot based on infection conditions
cross_data_runA <- mutate(cross_data_runA, cross = paste("cross", cross_id, sep = ""))

#Make a strain combination column for more convenient facetting
cross_data_runA <- mutate(cross_data_runA, strainAB = paste(strainA, strainB, sep = "_"))

#Join sample information
two_strain_database17_98_locus <- left_join(two_strain_database17_98_locus, cross_data_runA, by = "cross")

#Arrange facets by strain combination
ggplot(two_strain_database17_98_locus, aes(x = locus, y = sample)) + geom_point(aes(col=majority_strain, alpha=proportion_assigned), shape=15, size = 6) + theme(axis.text.x = element_text(angle=90)) + facet_wrap(~ strainAB)

#Summarise number of reads per cross 
two_strain_database17_98_locus %>% group_by(cross) %>%
  summarise(
    total_cross_reads = sum(total_reads),
  )
# A tibble: 9 x 2
#cross   total_cross_reads
#<chr>               <int>
#cross10            865786
#cross11            664352
#cross13            740611
#cross17           1199638
#cross2            1152974
#cross3            1197359
#cross4             859464
#cross5            1176569
#cross6             583357

#Summarise average number of reads per sample in each cross 
reads_sample <- two_strain_database17_98_locus %>% group_by(cross_sample, cross) %>%
  summarise(
    total_cross_reads = sum(total_reads),
  )

summary(reads_sample$total_cross_reads)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#30    5205    9378    9769   13690   24135

#Summarise number of reads per locus 
reads_locus <- two_strain_database17_98_locus %>% group_by(locus) %>%
  summarise(
    total_locus_reads = sum(total_reads),
    average_locus_reads = mean(total_reads),
  )

summary(reads_locus$average_locus_reads)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#14.63  121.79  686.32 1261.37 2575.73 3323.07

mean(reads_locus$average_locus_reads)
#1261.366
sd(reads_locus$average_locus_reads)
#1353.948

#Looking at the plots, need to clean-up HA/NA for SI x HK which is intersubtype
#For now, just going to take any HA1b or NA1c as evidence that HA should be the one

#Need to make this into a function:

## 2. Antigenic Segment Processing ----
# This section uses information from each of two PCR loci (HA3a, HA1b, NA1c, NA3a) for each antigenic segment (HA and NA)
# to decide which strain to assign to that given segment in each progeny isolate clone (sample).
# Aside from using read data, we also used information from control LMGSeq runs (including single strain LMGSeq to establish thresholds for how  primers work behave with just one strain background)

#Total Rows in DF
nrow(two_strain_database17_98_locus)
#[1] 6503

#Locus by cross
locus_by_cross <- two_strain_database17_98_locus %>% group_by(strainAB, locus) %>%
  summarise(
    number_samples = length(cross_sample),
  )

#NA's
View(locus_by_cross[grep("NA", locus_by_cross$locus), ])

#HA's
View(locus_by_cross[grep("HA", locus_by_cross$locus), ])

## 2.1 NA Segment  ----
#Lets check how many samples have any NA segments
nrow(two_strain_database17_98_locus[grep("NA", two_strain_database17_98_locus$locus), ])
#[1] 1193

#Put them in a data frame
na <- two_strain_database17_98_locus[grep("NA", two_strain_database17_98_locus$locus), ]

#Now check how many NA's called per sample/cross
na_by_sample <- group_by(na, cross_sample) %>%
  summarise(
    count = n()
  )

#How many have 2 NA loci called?
nrow(subset(na_by_sample, count == 2))
#[1] 331

#Join that information into the na data frame
na <- left_join(na, na_by_sample)

#How many have 2 loci called
nrow(subset(na, count == 2))
#[1] 662
#Remember this is twice the number above, because here not looking by sample

#Let's subset to these only, the others we can leave untouched
na_two <- subset(na, count == 2)

#First, there should be no HK68_PAN99
filter(na_two, strainAB == "HK68_PAN99")
#<0 rows> (or 0-length row.names)
#Looks good.

#Now down the line
#Using a single rule for all NA1/NA3 decisions based on manual inspection
na_two_wider <- na_two %>% group_by(cross_sample, locus) %>%
  summarise(
    majority_assigned_reads = sum(majority_assigned_reads),
  )

#View
#View(pivot_wider(na_two_wider, names_from = locus, values_from =  majority_assigned_reads))

#Implement
na_two_wider <- pivot_wider(na_two_wider, names_from = locus, values_from =  majority_assigned_reads)

#Number of rows
nrow(na_two_wider)
#331
#So one row per sample

ggplot(na_two_wider, aes(x = NA1c, y = NA3a)) + geom_point()

#So our strategy is to add a "remove" data frame, which will have a 'cross_sample_locus' column,
# that will then be used as a key to remove rows from the larger database

#Our starting rule is going to be NA3a's < 100 we want to give to NA1c, and NA1c's <3 we want to give to NA3a

#But first, let's just agree that any where NA1c is equal or majority, should be NA1c
#View(subset(na_two_wider, NA3a < NA1c))
nrow(subset(na_two_wider, NA3a <= NA1c))
#10

#So go ahead and mark the NA3a's in those samples as remove
#Create a data frame that I can then join to na_two_wider
cross_sample <- subset(na_two_wider, NA3a <= NA1c)$cross_sample
cross_sample_locus <- paste(subset(na_two_wider, NA3a <= NA1c)$cross_sample, "NA3a", sep = "_")

#Make temp remove data frame
remove_na <- data.frame(cross_sample,cross_sample_locus)
nrow(remove_na)
#10

#Now remove those from the na_two_wider data frame so I can whittle things down
#View(anti_join(na_two_wider, remove_na)) #321 looks good!
na_two_wider <- anti_join(na_two_wider, remove_na)

#Now among the ones where NA3a is greater than NA1c, let's assign remaining:
#We will recursively add to the remove DF, as we go through each condition, may circle back to this later

#First, remove those that do not meet criteria to assign to NA1c and add to 
cross_sample <- c(cross_sample, subset(na_two_wider, NA1c < 3)$cross_sample)
cross_sample_locus <- c(cross_sample_locus, paste(subset(na_two_wider, NA1c < 3)$cross_sample, "NA1c", sep = "_"))

remove_na <- data.frame(cross_sample,cross_sample_locus)
nrow(remove_na)
#168

#Now again remove those from the na_two_wider data frame so I can whittle things down
#View(anti_join(na_two_wider, remove_na)) #163, moving down the line
na_two_wider <- anti_join(na_two_wider, remove_na)

#Second, now for those sample that meet criteria to assign to NA3a, I need to add NA1c to the remove for those loci
#View(subset(na_two_wider, NA3a > 100))
nrow(subset(na_two_wider, NA3a > 100))
#159

#Add recursively
cross_sample <- c(cross_sample, subset(na_two_wider, NA3a > 100)$cross_sample)
cross_sample_locus <- c(cross_sample_locus, paste(subset(na_two_wider, NA3a > 100)$cross_sample, "NA1c", sep = "_"))

#Build remove_na df again
remove_na <- data.frame(cross_sample,cross_sample_locus)
nrow(remove_na)
#327

#Now for prob final time remove those from the na_two_wider data frame so I can whittle things down
#View(anti_join(na_two_wider, remove_na)) #4, just these left!
na_two_wider <- anti_join(na_two_wider, remove_na)

#Finally, these remaining 5 ones are clearly NA1c, so need to label them as NA3a for removal
#No subsetting needed
cross_sample <- c(cross_sample, na_two_wider$cross_sample)
cross_sample_locus <- c(cross_sample_locus, paste(na_two_wider$cross_sample, "NA3a", sep = "_"))

#Build remove_na df for a final time
remove_na <- data.frame(cross_sample,cross_sample_locus)
nrow(remove_na)
#331 goal!!!

#Now let's do a quick sanity check, visually with view
View(
  na_two %>% group_by(cross_sample, locus) %>%
    summarise(
      majority_assigned_reads = sum(majority_assigned_reads),
    ) %>% pivot_wider(names_from = locus, values_from =  majority_assigned_reads) %>% left_join(remove_na)
)
#Gosh I love the tidyverse!
#Looks good enough especially at the extremes

#Now let's remove these offending rows from the big data frame!

#First check what are big data frame looks like before modifying
nrow(two_strain_database17_98_locus)
#6503

#And recall how many of the samples had duplicate NA calls
nrow(subset(na_by_sample, count == 2))
#[1] 331

#Recall that is the number of samples, so I had 331*2 = 662 rows in the data frame assoc w/ 2 NA's,
# so want to get that down to this amount:
nrow(two_strain_database17_98_locus) - nrow(subset(na_by_sample, count == 2))
#6172

#Prepare to actually remove rows
#View(two_strain_database17_98_locus[!two_strain_database17_98_locus$cross_sample_locus %in% remove_na$cross_sample_locus, ])
#Beautiful, isn't it?

#Actually change data frame
two_strain_database17_98_locus <- two_strain_database17_98_locus[!two_strain_database17_98_locus$cross_sample_locus %in% remove_na$cross_sample_locus, ]


#Now check 
filter(two_strain_database17_98_locus, locus == "NA3a" | locus == "NA1c") %>% 
  group_by(cross_sample) %>%
  summarise(
    count = n()
  )
#All have only 1 NA locus

#Now rename all the NA loci
two_strain_database17_98_locus$locus[grep("NA", two_strain_database17_98_locus$locus)] <- "NA"
#DONE!!!!! 
## End NA


#### 2.2 HA Segment  ---- 
#Lets check how many samples have any HA segments
nrow(two_strain_database17_98_locus[grep("HA", two_strain_database17_98_locus$locus), ])
#[1] 843

#Put them in a data frame
ha <- two_strain_database17_98_locus[grep("HA", two_strain_database17_98_locus$locus), ]

#Now check how many HA's called per sample/cross
ha_by_sample <- group_by(ha, cross_sample) %>%
  summarise(
    count = n()
  )

#How many have 2 HA loci called?
nrow(subset(ha_by_sample, count == 2))
#[1] 52

#Join that information into the ha data frame
ha <- left_join(ha, ha_by_sample)

#How many have 2 loci called
nrow(subset(ha, count == 2))
#[1] 104
#Remember this is twice the number above, because here not looking by sample

#Let's subset to these only, the others we can leave untouched
ha_two <- subset(ha, count == 2)

#First, there should be no CA09_SI86 samples or HK68_PAN99
filter(ha_two, strainAB == "HK68_PAN99")
#<0 rows> (or 0-length row.names)
#Looks good.

#Now down the line
#Using a single rule for all HA1/HA3 decisions based on manual inspection
ha_two_wider <- ha_two %>% group_by(cross_sample, locus) %>%
  summarise(
    majority_assigned_reads = sum(majority_assigned_reads),
  )

#View
#View(pivot_wider(ha_two_wider, names_from = locus, values_from =  majority_assigned_reads))

#Implement
ha_two_wider <- pivot_wider(ha_two_wider, names_from = locus, values_from =  majority_assigned_reads)

#Number of rows
nrow(ha_two_wider)
#52
#So one row per sample

ggplot(ha_two_wider, aes(x = HA1b, y = HA3a)) + geom_point()

#So our strategy is to add a "remove" data frame, which will have a 'cross_sample_locus' column,
# that will then be used as a key to remove rows from the larger database

#First, let's just agree that any sample where HA1b is equal or majority, the HA3a locus should be marked for removal
#View(subset(ha_two_wider, HA3a <= HA1b))
nrow(subset(ha_two_wider, HA3a <= HA1b))
#2

#So go ahead and mark the HA3a's in those samples as remove
#Create a data frame that I can then join to ha_two_wider
cross_sample <- subset(ha_two_wider, HA3a <= HA1b)$cross_sample
cross_sample_locus <- paste(subset(ha_two_wider, HA3a <= HA1b)$cross_sample, "HA3a", sep = "_")

#Make temp remove data frame
remove_ha <- data.frame(cross_sample,cross_sample_locus)
nrow(remove_ha)
#2

#Now remove those from the na_two_wider data frame so I can whittle things down
#View(anti_join(ha_two_wider, remove_ha)) #50 made a small dent
ha_two_wider <- anti_join(ha_two_wider, remove_ha)

#Now among the ones where HA3a is greater than HA1b, let's assign remaining:
#We will recursively add to the remove DF, as we go through each condition, may circle back to this later

#Controls show that HA1b can't really be amplified by HA3a's. But very conservatively I will exclude 
# singleton HA1b's across the board
#View(subset(ha_two_wider, HA1b == 1))
nrow(subset(ha_two_wider, HA1b == 1))
#31

#First, remove singleton HA1b's
cross_sample <- c(cross_sample, subset(ha_two_wider, HA1b == 1)$cross_sample)
cross_sample_locus <- c(cross_sample_locus, paste(subset(ha_two_wider, HA1b == 1)$cross_sample, "HA1b", sep = "_"))

remove_ha <- data.frame(cross_sample,cross_sample_locus)
nrow(remove_ha)
#33

#Now again remove those from the ha_two_wider data frame so I can whittle things down
#View(anti_join(ha_two_wider, remove_ha)) #33, making progress now
ha_two_wider <- anti_join(ha_two_wider, remove_ha)

#Second, will accept all HA1b's >2, when HA3a's are < 1000. So remove HA3a's under 1000. 
#View(subset(ha_two_wider,  HA3a < 1000))
nrow(subset(ha_two_wider,  HA3a < 1000))
#8

#Add recursively
cross_sample <- c(cross_sample, subset(ha_two_wider,  HA3a < 1000)$cross_sample)
cross_sample_locus <- c(cross_sample_locus, paste(subset(ha_two_wider,  HA3a < 1000)$cross_sample, "HA3a", sep = "_"))

#Build remove_ha df again
remove_ha <- data.frame(cross_sample,cross_sample_locus)
nrow(remove_ha)
#41

#Now for prob final time remove those from the ha_two_wider data frame so I can whittle things down
#View(anti_join(ha_two_wider, remove_ha)) #11, just these left!
ha_two_wider <- anti_join(ha_two_wider, remove_ha)

#These are pretty clear HA3a's maaaaybe one or two potential toss ups close to 1000 HA3a reads but giving to HA3a
#No subsetting needed
cross_sample <- c(cross_sample, ha_two_wider$cross_sample)
cross_sample_locus <- c(cross_sample_locus, paste(ha_two_wider$cross_sample, "HA1b", sep = "_"))

#Build remove_ha df for a final time
remove_ha <- data.frame(cross_sample,cross_sample_locus)
nrow(remove_ha)
#52 goal!!!

#Now let's do a quick sanity check, visually with view
View(
  ha_two %>% group_by(cross_sample, locus) %>%
    summarise(
      majority_assigned_reads = sum(majority_assigned_reads),
    ) %>% pivot_wider(names_from = locus, values_from =  majority_assigned_reads) %>% left_join(remove_ha)
)
#Gosh I love the tidyverse!
#Looks good enough especially at the extremes

#Now let's remove these offending rows from the big data frame!

#First check what are big data frame looks like before modifying
nrow(two_strain_database17_98_locus)
#6172

#And recall how many of the samples had duplicate NA calls
nrow(subset(ha_by_sample, count == 2))
#[1] 52

#Recall that is the number of samples, so I had 52*2 = 104 rows in the data frame assoc w/ 2 HA's,
# so want to get that down to this amount:
nrow(two_strain_database17_98_locus) - nrow(subset(ha_by_sample, count == 2))
#6120

#Prepare to actually remove rows
#View(two_strain_database17_98_locus[!two_strain_database17_98_locus$cross_sample_locus %in% remove_ha$cross_sample_locus, ])
#Beautiful, isn't it?

#Actually change data frame
two_strain_database17_98_locus <- two_strain_database17_98_locus[!two_strain_database17_98_locus$cross_sample_locus %in% remove_ha$cross_sample_locus, ]


#Now check 
filter(two_strain_database17_98_locus, locus == "HA3a" | locus == "HA1b") %>% 
  group_by(cross_sample) %>%
  summarise(
    count = n()
  )
#All have only 1 HA locus

#Now rename all the HA loci
two_strain_database17_98_locus$locus[grep("HA", two_strain_database17_98_locus$locus)] <- "HA"
# End HA ##

#### 3. Basic Reassortment Plots ####
#Now should be able to generate the plots for each of the conditions
#First, replicate of standard (using crosess to subset because it's a bit complicated)
ggplot(subset(two_strain_database17_98_locus, cross == "cross2" | cross == "cross4" | cross == "cross3" | cross == "cross5" | cross == "cross13" | cross == "cross17"), aes(x = locus, y = sample)) + geom_point(aes(col=majority_strain, alpha=proportion_assigned), shape=15, size = 6) + theme(axis.text.x = element_text(angle=90)) + facet_wrap(strainAB ~ cross, nrow = 3)
#Leaving out cross == "cross11" to make it pairs

#Second by MOI
ggplot(subset(two_strain_database17_98_locus, cross == "cross4" | cross == "cross10" | cross == "cross6"), aes(x = locus, y = sample)) + geom_point(aes(col=majority_strain, alpha=proportion_assigned), shape=15, size = 6) + theme(axis.text.x = element_text(angle=90)) + facet_wrap(strainAB ~ as.numeric(moi), ncol =3)
#-Run A Cross 4 MOI = 10
#-Run A Cross 10 MOI = 0.01
#-Run A Cross 6 MOI = 1

#Now run some numbers and stats tests
  #t-test for trypsin no-trypsin (may need to control for strains)
  #regression? ordered ANOVA for timepoints OR binomial regression? Reassortant/not 

#Plot
ggplot(two_strain_database17_98_locus, aes(x = locus, y = sample)) + geom_point(aes(col=majority_strain, alpha=proportion_assigned), shape=15, size = 6) + theme(axis.text.x = element_text(angle=90)) + facet_wrap(~ cross)

controls_df <- two_strain_database17_98_locus

#First need to get the number of typed loci for each sample
controls_df_loci_numbers <- controls_df %>% 
  group_by(cross_sample) %>% 
  summarise(
    total_segments = length(locus),
    max_majority_strain = max(table(majority_strain))
  )

controls_df <- right_join(controls_df, controls_df_loci_numbers)

#Now should be able to generate the plots for each of the conditions
#First, replicate of standard (using crosess to subset because it's a bit complicated)
ggplot(subset(controls_df, cross == "cross2" & total_segments > 6 | cross == "cross4" & total_segments > 6 | cross == "cross3" & total_segments > 6 | cross == "cross5"  & total_segments > 6 | cross == "cross11"  & total_segments > 6 | cross == "cross13"  & total_segments > 6 | cross == "cross15"  & total_segments > 6), aes(x = locus, y = sample)) + geom_point(aes(col=majority_strain, alpha=proportion_assigned), shape=15, size = 6) + theme(axis.text.x = element_text(angle=90)) + facet_wrap(strainAB ~ cross, nrow =3)


#Second by MOI
ggplot(subset(controls_df, cross == "cross4" & total_segments > 6 | cross == "cross10" & total_segments > 6 | cross == "cross6" & total_segments > 6), aes(x = locus, y = sample)) + geom_point(aes(col=majority_strain, alpha=proportion_assigned), shape=15, size = 6) + theme(axis.text.x = element_text(angle=90)) + facet_wrap(strainAB ~ as.numeric(moi), ncol =3)
#-Run A Cross 4 MOI = 10
#-Run A Cross 10 MOI = 0.01
#-Run A Cross 6 MOI = 1

#Lose a ton of samples, should look into removing PB1, if no different keep all data

#Now sort by parentals
#First, replicate of standard (using crosess to subset because it's a bit complicated)
ggplot(subset(controls_df, cross == "cross2" | cross == "cross4" | cross == "cross3" | cross == "cross5" | cross == "cross11" | cross == "cross13" | cross == "cross17"), aes(x = locus, y = reorder(sample, max_majority_strain))) + geom_point(aes(col=majority_strain, alpha=proportion_assigned), shape=15, size = 6) + theme(axis.text.x = element_text(angle=90)) + facet_wrap(strainAB ~ cross, nrow =3)

#Second by MOI
ggplot(subset(controls_df, cross == "cross4" | cross == "cross10" | cross == "cross6"), aes(x = locus, y = reorder(sample, max_majority_strain))) + geom_point(aes(col=majority_strain, alpha=proportion_assigned), shape=15, size = 6) + theme(axis.text.x = element_text(angle=90)) + facet_wrap(strainAB ~ as.numeric(moi), ncol =3)
#-Run A Cross 2 MOI = 10
#-Run A Cross 10 MOI = 0.01
#-Run A Cross 6 MOI = 1

#### 4. Calculation of Reassortment Rates, using 8 segments ####
#Now let's generate some numbers!
parentals <- subset(controls_df, total_segments > 7) %>% 
  group_by(cross_sample, cross, sample) %>% 
  summarise(
  parental = max_majority_strain == total_segments,
  #total_segments == max_majority_strain
)

parentals <- distinct(parentals)

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
#305
nrow(controls_df_number_parents)
#[1] 864

#Make a reasortant column 
controls_df_number_parents <- mutate(controls_df_number_parents, reassortant = ifelse(number_parents > 1, yes = 1, no = 0))

cross_stats <- controls_df_number_parents %>% 
  group_by(cross) %>% 
  summarise(
    clones = length(number_parents),
    reassortants = sum(reassortant)
    #constellation = paste(majority_strain_locus, sep = ","),
    #max_majority_strain = max(table(majority_strain))
  )

#Add proportion reassortant to data frame
cross_stats <- mutate(cross_stats, prop_reassortant = reassortants / clones)

#Now add cross data
cross_stats <- left_join(cross_stats, cross_data_runA, by = "cross")

#### 5. Figures Plotting Reasortment Proportions and Replicate Deviation Data ####
#First, replicate of standard (using crosess to subset because it's a bit complicated)
#Figure 1A
ggplot(subset(cross_stats, cross == "cross2" | cross == "cross4" | cross == "cross3" | cross == "cross5" | cross == "cross11" | cross == "cross13" | cross == "cross17"), aes(x = strainAB, y = prop_reassortant)) + xlab("Experimental Coinfection: Strain Combination") + ylab("Proportion of Reassortant Plaque Isolates") + geom_point(size = 3) + theme_tufte() + theme(text = element_text(size = 20,  family="Helvetica")) + theme(axis.text.x = element_text(size = 12)) + theme(legend.position = "none") + theme(panel.grid.major.y = element_line(color = "lightgray",size = 0.5)) 

#Difference in each reassortment rate estimate within 
reassortment_rate_replicate_estimate <- subset(cross_stats, cross == "cross2" | cross == "cross4" | cross == "cross3" | cross == "cross5" | cross == "cross11" | cross == "cross13" | cross == "cross17")

subset(reassortment_rate_replicate_estimate, strainAB == "CA09_HK68")
#0.156 - 0.0417 = 0.1143

subset(reassortment_rate_replicate_estimate, strainAB == "CA09_PAN99")
#0.5-0.438 = 0.062
#0.5 - 0.406 = 0.094
#0.438 - 0.406 = 0.032

subset(reassortment_rate_replicate_estimate, strainAB == "HK68_PAN99")
#0.781-0.771 = 0.01

#Average deviation 
mean(c(0.1143, 0.062, 0.094, 0.032, 0.01))
#0.06246

#Figure 1B
#Second by MOI
ggplot(subset(cross_stats, cross == "cross10" | cross == "cross6" | cross == "cross4"), aes(x = moi, y = prop_reassortant)) + xlab("Multiplicity of Infection") + ylab("Proportion of Reassortant Plaque Isolates") + scale_x_log10() + geom_point(size = 3) + geom_smooth(method = "lm", formula = y ~ exp(x), se = T) + theme_tufte() + theme(text = element_text(size = 20,  family="Helvetica")) + theme(axis.text.x = element_text(size = 12)) + theme(legend.position = "none") + theme(panel.grid.major.y = element_line(color = "lightgray",size = 0.5)) 