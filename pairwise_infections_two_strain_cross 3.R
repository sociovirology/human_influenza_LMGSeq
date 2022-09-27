################ Culture coinfection analyses ################
#### 1. Data Sources and Preparing Data
#### 2. 
#### 3. Coinfection Supernatant Titers
#### 4. Population Genetic Statistics
#### 5. Testing Strain Specificity

#AFTER SENDING DRAFT - Search for this below

#Load libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape2)
library(readr)
library(gridExtra)
library(ggthemes)

#### 1. Data Sources and Preparing Data ####

#Borrow tick from strain assignment function I developed.

#Read files
file_names <- list.files("outputs/pairwise_infections", "*_strain_database17_98_locus.txt", full.names = T)
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
#[1] 7538

#Let's try our first plot
ggplot(two_strain_database17_98_locus, aes(x = locus, y = sample)) + geom_point(aes(col=majority_strain, alpha=proportion_assigned), shape=15, size = 6) + theme(axis.text.x = element_text(angle=90)) + facet_wrap(~ cross)

no_subtype_processing_plot <- ggplot(two_strain_database17_98_locus, aes(x = locus, y = sample)) + geom_point(aes(col=majority_strain, alpha=proportion_assigned), shape=15, size = 6) + theme(axis.text.x = element_text(angle=90)) + facet_wrap(~ cross)


#Okay, now add the information on what each cross is about
cross_data_runA <- read.csv("~/Dropbox/mixtup/Documentos/ucdavis/papers/influenza_GbBSeq_human/influenza_GbBSeq/data/cross_data_runA.csv")

#Prepare the cross_data_runA df for a merge with the big data set, so we can plot based on infection conditions
cross_data_runA <- mutate(cross_data_runA, cross = paste("cross", cross_id, sep = ""))

#Make a strain combination column for more convenient facetting
cross_data_runA <- mutate(cross_data_runA, strainAB = paste(strainA, strainB, sep = "_"))

#Join sample information
two_strain_database17_98_locus <- left_join(two_strain_database17_98_locus, cross_data_runA, by = "cross")

#Arrange facets by strain
ggplot(two_strain_database17_98_locus, aes(x = locus, y = sample)) + geom_point(aes(col=majority_strain, alpha=proportion_assigned), shape=15, size = 6) + theme(axis.text.x = element_text(angle=90)) + facet_grid(strainA ~ strainB)

#Summarise number of reads per cross 
two_strain_database17_98_locus %>% group_by(strainAB) %>%
  summarise(
    total_cross_reads = sum(total_reads),
  )
#strainAB   total_cross_reads
#<chr>                  <int>
#1 CA09_HK68            1,154,666
#2 CA09_PAN99           664,562
#3 CA09_SI86            907,760
#4 CA09_TX12            1,128,965
#5 HK68_PAN99           1,197,359
#6 HK68_SI86            1,053,568
#7 HK68_TX12            1,383,315
#8 PAN99_SI86           891,395
#9 PAN99_TX12           1,544,278
#10 SI86_TX12           221,086

#Summarise average number of reads per sample in each cross 
reads_sample <- two_strain_database17_98_locus %>% group_by(cross_sample, cross) %>%
  summarise(
    total_cross_reads = sum(total_reads),
  )

summary(reads_sample$total_cross_reads)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#44    5595   10650   10570   15142   33184

#Summarise number of reads per locus 
two_strain_database17_98_locus %>% group_by(strainAB, locus) %>%
  summarise(
    total_locus_reads = sum(total_reads),
    average_locus_reads = mean(total_reads),
  )

#################### New Antigenic Segment Processing - Refer ha_na.Rmd ####################
#Total Rows in DF
nrow(two_strain_database17_98_locus)
#[1] 7538

#Locus by cross
locus_by_cross <- two_strain_database17_98_locus %>% group_by(strainAB, locus) %>%
  summarise(
    number_samples = length(cross_sample),
  )

#NA's
View(locus_by_cross[grep("NA", locus_by_cross$locus), ])

#HA's
View(locus_by_cross[grep("HA", locus_by_cross$locus), ])

##### First NA #####
#Lets check how many samples have any NA segments
nrow(two_strain_database17_98_locus[grep("NA", two_strain_database17_98_locus$locus), ])
#[1] 1300

#Put them in a data frame
na <- two_strain_database17_98_locus[grep("NA", two_strain_database17_98_locus$locus), ]

#Now check how many NA's called per sample/cross
na_by_sample <- group_by(na, cross_sample) %>%
  summarise(
    count = n()
  )

#How many have 2 NA loci called?
nrow(subset(na_by_sample, count == 2))
#[1] 347

#Join that information into the na data frame
na <- left_join(na, na_by_sample)

#How many have 2 loci called
nrow(subset(na, count == 2))
#[1] 694
#Remember this is twice the number above, because here not looking by sample

#Let's subset to these only, the others we can leave untouched
na_two <- subset(na, count == 2)

#First, there should be no CA09_SI86 samples or HK68_PAN99
filter(na_two, strainAB == "CA09_SI86" | strainAB == "HK68_PAN99")
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
#347
#So one row per sample

ggplot(na_two_wider, aes(x = NA1c, y = NA3a)) + geom_point()

#So our strategy is to add a "remove" data frame, which will have a 'cross_sample_locus' column,
# that will then be used as a key to remove rows from the larger database

#Our starting rule is going to be NA3a's < 100 we want to give to NA1c, and NA1c's <3 we want to give to NA3a

#But first, let's just agree that any where NA1c is equal or majority, should be NA1c
#View(subset(na_two_wider, NA3a < NA1c))
nrow(subset(na_two_wider, NA3a <= NA1c))
#141

#So go ahead and mark the NA3a's in those samples as remove
#Create a data frame that I can then join to na_two_wider
cross_sample <- subset(na_two_wider, NA3a <= NA1c)$cross_sample
cross_sample_locus <- paste(subset(na_two_wider, NA3a <= NA1c)$cross_sample, "NA3a", sep = "_")

#Make temp remove data frame
remove_na <- data.frame(cross_sample,cross_sample_locus)
nrow(remove_na)
#141

#Now remove those from the na_two_wider data frame so I can whittle things down
#View(anti_join(na_two_wider, remove_na)) #206 looks good!
na_two_wider <- anti_join(na_two_wider, remove_na)

#Now among the ones where NA3a is greater than NA1c, let's assign remaining:
#We will recursively add to the remove DF, as we go through each condition, may circle back to this later

#First, remove those that do not meet criteria to assign to NA1c and add to 
cross_sample <- c(cross_sample, subset(na_two_wider, NA1c < 3)$cross_sample)
cross_sample_locus <- c(cross_sample_locus, paste(subset(na_two_wider, NA1c < 3)$cross_sample, "NA1c", sep = "_"))

remove_na <- data.frame(cross_sample,cross_sample_locus)
nrow(remove_na)
#201

#Now again remove those from the na_two_wider data frame so I can whittle things down
#View(anti_join(na_two_wider, remove_na)) #146, moving down the line
na_two_wider <- anti_join(na_two_wider, remove_na)

#Second, now for those sample that meet criteria to assign to NA3a, I need to add NA1c to the remove for those loci
#View(subset(na_two_wider, NA3a > 100))
nrow(subset(na_two_wider, NA3a > 100))
#141

#Add recursively
cross_sample <- c(cross_sample, subset(na_two_wider, NA3a > 100)$cross_sample)
cross_sample_locus <- c(cross_sample_locus, paste(subset(na_two_wider, NA3a > 100)$cross_sample, "NA1c", sep = "_"))

#Build remove_na df again
remove_na <- data.frame(cross_sample,cross_sample_locus)
nrow(remove_na)
#342

#Now for prob final time remove those from the na_two_wider data frame so I can whittle things down
#View(anti_join(na_two_wider, remove_na)) #5, just these left!
na_two_wider <- anti_join(na_two_wider, remove_na)

#Finally, these remaining 5 ones are toss ups - I will simply give to NA3a, so need to label them as NA1c
#No subsetting needed
cross_sample <- c(cross_sample, na_two_wider$cross_sample)
cross_sample_locus <- c(cross_sample_locus, paste(na_two_wider$cross_sample, "NA1c", sep = "_"))

#Build remove_na df for a final time
remove_na <- data.frame(cross_sample,cross_sample_locus)
nrow(remove_na)
#347 goal!!!

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
#7538

#And recall how many of the samples had duplicate NA calls
nrow(subset(na_by_sample, count == 2))
#[1] 347

#Recall that is the number of samples, so I had 347*2 = 694 rows in the data frame assoc w/ 2 NA's,
# so want to get that down to this amount:
nrow(two_strain_database17_98_locus) - nrow(subset(na_by_sample, count == 2))
#7191

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

##### End NA #####


##### Second HA #####
#Lets check how many samples have any HA segments
nrow(two_strain_database17_98_locus[grep("HA", two_strain_database17_98_locus$locus), ])
#[1] 1124

#Put them in a data frame
ha <- two_strain_database17_98_locus[grep("HA", two_strain_database17_98_locus$locus), ]

#Now check how many HA's called per sample/cross
ha_by_sample <- group_by(ha, cross_sample) %>%
  summarise(
    count = n()
  )

#How many have 2 HA loci called?
nrow(subset(ha_by_sample, count == 2))
#[1] 235

#Join that information into the ha data frame
ha <- left_join(ha, ha_by_sample)

#How many have 2 loci called
nrow(subset(ha, count == 2))
#[1] 470
#Remember this is twice the number above, because here not looking by sample

#Let's subset to these only, the others we can leave untouched
ha_two <- subset(ha, count == 2)

#First, there should be no CA09_SI86 samples or HK68_PAN99
filter(ha_two, strainAB == "CA09_SI86" | strainAB == "HK68_PAN99")
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
#235
#So one row per sample

ggplot(ha_two_wider, aes(x = HA1b, y = HA3a)) + geom_point()

#So our strategy is to add a "remove" data frame, which will have a 'cross_sample_locus' column,
# that will then be used as a key to remove rows from the larger database

#Our starting rule is going to be NA3a's < 100 we want to give to NA1c, and NA1c's <3 we want to give to NA3a

#But first, let's just agree that any sample where HA1b is equal or majority, the HA3a locus should be marked for removal
#View(subset(ha_two_wider, HA3a <= HA1b))
nrow(subset(ha_two_wider, HA3a <= HA1b))
#36

#So go ahead and mark the HA3a's in those samples as remove
#Create a data frame that I can then join to ha_two_wider
cross_sample <- subset(ha_two_wider, HA3a <= HA1b)$cross_sample
cross_sample_locus <- paste(subset(ha_two_wider, HA3a <= HA1b)$cross_sample, "HA3a", sep = "_")

#Make temp remove data frame
remove_ha <- data.frame(cross_sample,cross_sample_locus)
nrow(remove_ha)
#36

#Now remove those from the na_two_wider data frame so I can whittle things down
#View(anti_join(ha_two_wider, remove_ha)) #199 made a small dent
ha_two_wider <- anti_join(ha_two_wider, remove_ha)

#Now among the ones where HA3a is greater than HA1b, let's assign remaining:
#We will recursively add to the remove DF, as we go through each condition, may circle back to this later

#Controls show that HA1b can't really be amplified by HA3a's. But very conservatively I will exclude 
# singleton HA1b's across the board
#View(subset(ha_two_wider, HA1b == 1))
nrow(subset(ha_two_wider, HA1b == 1))
#47

#First, remove singleton HA1b's
cross_sample <- c(cross_sample, subset(ha_two_wider, HA1b == 1)$cross_sample)
cross_sample_locus <- c(cross_sample_locus, paste(subset(ha_two_wider, HA1b == 1)$cross_sample, "HA1b", sep = "_"))

remove_ha <- data.frame(cross_sample,cross_sample_locus)
nrow(remove_ha)
#83

#Now again remove those from the ha_two_wider data frame so I can whittle things down
#View(anti_join(ha_two_wider, remove_ha)) #152, making progress now
ha_two_wider <- anti_join(ha_two_wider, remove_ha)

#Second, will accept all HA1b's >2, when HA3a's are < 1000. So remove HA3a's under 1000. 
#View(subset(ha_two_wider,  HA3a < 1000))
nrow(subset(ha_two_wider,  HA3a < 1000))
#127

#Add recursively
cross_sample <- c(cross_sample, subset(ha_two_wider,  HA3a < 1000)$cross_sample)
cross_sample_locus <- c(cross_sample_locus, paste(subset(ha_two_wider,  HA3a < 1000)$cross_sample, "HA3a", sep = "_"))

#Build remove_ha df again
remove_ha <- data.frame(cross_sample,cross_sample_locus)
nrow(remove_ha)
#210

#Now for prob final time remove those from the ha_two_wider data frame so I can whittle things down
#View(anti_join(ha_two_wider, remove_ha)) #25, just these left!
ha_two_wider <- anti_join(ha_two_wider, remove_ha)

#These are pretty clear HA3a's just a few potential toss ups close to 1000 HA3a reads but giving to HA3a
#No subsetting needed
cross_sample <- c(cross_sample, ha_two_wider$cross_sample)
cross_sample_locus <- c(cross_sample_locus, paste(ha_two_wider$cross_sample, "HA1b", sep = "_"))

#Build remove_ha df for a final time
remove_ha <- data.frame(cross_sample,cross_sample_locus)
nrow(remove_ha)
#235 goal!!!

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
#7191

#And recall how many of the samples had duplicate NA calls
nrow(subset(ha_by_sample, count == 2))
#[1] 235

#Recall that is the number of samples, so I had 235*2 = 470 rows in the data frame assoc w/ 2 HA's,
# so want to get that down to this amount:
nrow(two_strain_database17_98_locus) - nrow(subset(ha_by_sample, count == 2))
#6956

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


##### End HA #####

#Checking how good the assignments are
mean(two_strain_database17_98_locus$proportion_assigned)
#[1] 0.9416985
sd(two_strain_database17_98_locus$proportion_assigned)
#[1] 0.1118182

#Pretty good overall
#How many sample/locus combinations are over a certain proportion 
nrow(subset(two_strain_database17_98_locus, proportion_assigned > .80)) / nrow(two_strain_database17_98_locus)
#0.8858539
nrow(subset(two_strain_database17_98_locus, proportion_assigned > .75)) / nrow(two_strain_database17_98_locus)
#0.9097182
nrow(subset(two_strain_database17_98_locus, proportion_assigned > .70)) / nrow(two_strain_database17_98_locus)
#0.93387
nrow(subset(two_strain_database17_98_locus, proportion_assigned > .65)) / nrow(two_strain_database17_98_locus)
#0.9544278

#Histogram
ggplot(two_strain_database17_98_locus, aes(x = proportion_assigned)) + geom_histogram()


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

#Once we have numbers, sort by parentals and limit to complete genotypes
ggplot(subset(controls_df, total_segments > 7), aes(x = locus, y = reorder(sample, max_majority_strain))) + geom_point(aes(col=majority_strain, alpha=proportion_assigned), shape=15, size = 6) + theme(axis.text.x = element_text(angle=90)) + facet_wrap(~ cross, scales = "free_y")

#Once we have numbers, sort by parentals and limit to complete genotypes and pairwise grid
ggplot(subset(controls_df, total_segments > 7), aes(x = locus, y = sample)) + geom_point(aes(col=majority_strain, alpha=proportion_assigned), shape=15, size = 6) + theme(axis.text.x = element_text(angle=90), axis.text.y = element_blank()) + facet_grid(strainA ~ strainB, scales = "free_y")

#Pairwise grid
ggplot(controls_df, aes(x = locus, y = sample)) + geom_point(aes(col=majority_strain, alpha=proportion_assigned), shape=15, size = 6) + theme(axis.text.x = element_text(angle=90), axis.text.y = element_blank()) + facet_grid(strainA ~ strainB, scales = "free_y")

#Not ordered
ggplot(controls_df, aes(x = locus, y = sample)) + geom_point(aes(col=majority_strain, alpha=proportion_assigned), shape=15, size = 6) + theme(axis.text.x = element_text(angle=90)) + facet_wrap(~ cross)

#Not ordered
ggplot(controls_df, aes(x = locus, y = sample)) + geom_point(aes(col=majority_strain, alpha=proportion_assigned), shape=15, size = 6) + theme(axis.text.x = element_text(angle=90)) + theme(axis.text.y = element_blank()) + facet_grid(strainA ~ strainB, scales = "free_y")
#Not ordered at least 7 segments
ggplot(subset(controls_df, total_segments > 6), aes(x = locus, y = sample)) + geom_point(aes(col=majority_strain, alpha=proportion_assigned), shape=15, size = 6) + theme(axis.text.x = element_text(angle=90)) + theme(axis.text.y = element_blank()) + facet_grid(strainA ~ strainB, scales = "free_y")

#Figure 2A - Not ordered, all samples genotypes, order segments 
ggplot(controls_df, aes(x = factor(locus, level = c('PB2f', 'PB1c', 'PAc', 'HA', 'NPd', 'NA', 'Mg', 'NS1d')), y = sample)) + geom_point(aes(col=majority_strain), shape=15, size = 6) + theme(axis.text.x = element_text(angle=90)) + xlab("Segment") + ylab("Plaque Isolate") + theme(axis.text.y = element_blank()) + facet_grid(strainA ~ strainB, scales = "free_y")
#Let's spruce up re: Elodie comments 
ggplot(controls_df, aes(x = factor(locus, level = c('PB2f', 'PB1c', 'PAc', 'HA', 'NPd', 'NA', 'Mg', 'NS1d')), y = sample)) + geom_point(aes(col=majority_strain), shape=15, size = 6) + xlab("Segment") + ylab("Plaque Isolate") + facet_grid(strainA ~ strainB, scales = "free_y") + theme_tufte() +  theme(text = element_text(size = 20,  family="Helvetica")) + theme(axis.text.x = element_text(angle=90, size = 17)) + theme(axis.text.y = element_blank()) + theme(strip.text = element_text(face = "bold")) + theme(legend.position = "none")      


#Now let's generate some numbers!
parentals <- controls_df %>% 
  group_by(cross_sample, cross, sample) %>% 
  summarise(
    parental = max_majority_strain == total_segments,
    #total_segments == max_majority_strain
  )

#Complete Genotypes Only
#parentals <- subset(controls_df, total_segments > 7) %>% 
#  group_by(cross_sample, cross, sample) %>% 
#  summarise(
#  parental = max_majority_strain == total_segments,
#  #total_segments == max_majority_strain
#)

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

#Complete Genotypes Only 
#Group by cross_sample and calculate number of parents for each
#controls_df_number_parents <- subset(controls_df, total_segments > 7) %>% 
#  group_by(cross, cross_sample) %>% 
#  summarise(
#    number_parents = length(unique(majority_strain)),
#    #constellation = paste(majority_strain_locus, sep = ","),
#    max_majority_strain = max(table(majority_strain))
#  )

nrow(subset(controls_df_number_parents, number_parents == 2))
#470
nrow(controls_df_number_parents)
#[1] 960

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

#AFTER SENDING DRAFT - APRIL 28 
# Make a supplementary figure that has only 8 segment, 7 segment, or any amount of segments 

#Going forward with allowing incomplete genotypes, because this is conservative with respect to reassortemnt, because missing segments decrease the chance of detecting a reassortant
ggplot(controls_df, aes(x = locus, y = reorder(sample, max_majority_strain))) + geom_point(aes(col=majority_strain, alpha=proportion_assigned), shape=15, size = 6) + theme(axis.text.x = element_text(angle=90)) + facet_wrap(~ cross)

#Plot reasortment proportions
ggplot(cross_stats, aes(x = cross, y = prop_reassortant)) + geom_col()

#Same but reorder by height
ggplot(cross_stats, aes(x = reorder(cross, prop_reassortant), y = prop_reassortant)) + geom_col()

#Same but reorder by height with line for 40% reassortment but then also for the theoretical free reassortment
ggplot(cross_stats, aes(x = reorder(cross, prop_reassortant), y = prop_reassortant)) + geom_col() + geom_hline(yintercept = 0.40, linetype = 2, color = "grey", alpha = 0.75) + geom_hline(yintercept = 0.9921875, linetype = 2, color = "red", alpha = 0.75)

#Add strainAB information to cross_stats DF
cross_stats <- left_join(cross_stats, cross_data_runA, by = "cross")

#Same but reorder by height with line for 40% reassortment but then also for the theoretical free reassortment and add strain info on x axis 
ggplot(cross_stats, aes(x = reorder(strainAB, prop_reassortant), y = prop_reassortant)) + geom_col() + ylab("Proportion of Reassortant Plaque Isolates") + xlab("Strains in Experimental Coinfection") + geom_hline(yintercept = 0.40, linetype = 2, color = "grey", alpha = 0.75) + geom_hline(yintercept = 0.9921875, linetype = 2, color = "red", alpha = 0.75)
#Update for better publication quality
ggplot(cross_stats, aes(x = reorder(strainAB, prop_reassortant), y = prop_reassortant)) + geom_col() + ylab("Proportion of Reassortant Plaque Isolates") + xlab("Strains in Experimental Coinfection") + geom_hline(yintercept = 0.40, linetype = 2, color = "grey", alpha = 0.75) + geom_hline(yintercept = 0.9921875, linetype = 2, color = "red", alpha = 0.75) + theme_tufte() +  theme(text = element_text(size = 20,  family="Helvetica")) + theme(axis.text.x = element_text(size = 12)) + theme(legend.position = "none") + theme(panel.grid.major.y = element_line(color = "lightgray",size = 0.5)) 

#Now more stats
mean(cross_stats$prop_reassortant)
#[1] 0.4895833

sd(cross_stats$prop_reassortant)
#[1] 0.3007834

min(cross_stats$prop_reassortant)
#[1] 0.04166667

max(cross_stats$prop_reassortant)
#[1] 0.9270833

#Let's calculate some stats on the expected vs observed number of reasortants

#Reassortment as success. Perfect reassortment expectation is calculated as 254/256, i.e. 2^8 = 256 unique combinations, two of those subtracted becausecause they are parentals
#So first let's do an example "manually" with CA09_SI86/cross18
binom.test(c(89, 7), p = 254/256)
#Exact binomial test
#data:  c(89, 7)
#number of successes = 89, number of trials = 96, p-value = 1.153e-05
#alternative hypothesis: true probability of success is not equal to 0.9921875
#95 percent confidence interval:
#  0.8555203 0.9701822
#sample estimates:
#  probability of success 
#0.9270833

controls_df <- right_join(controls_df, cross_stats)

#Does strain identity affect the proportion of reassortants?
ggplot(controls_df, aes(x = cross, y = prop_reassortant)) + geom_point(aes(color = strainA))

#Does strain identity affect the proportion of reassortants? Heatmap attempt
ggplot(controls_df, aes(x = strainA, y = strainB, fill= prop_reassortant)) + geom_tile()

#Heatmap! This is Figure 3A
ggplot(controls_df, aes(x = strainA, y = strainB, fill= prop_reassortant)) + geom_tile() + geom_text(aes(label = round(prop_reassortant, 2))) + theme(legend.position = "none") + scale_fill_gradient(low = "white", high = "red") 

#This doesn't work, because of double counting of strains 
#model1 <- lm(prop_reassortant ~ strainA + strainB, controls_df)

#Let's calculate the averages and SD's of each
#SI86
mean(subset(cross_stats, strainA == "SI86" | strainB == "SI86", select = prop_reassortant, drop = T))
#[1] 0.6484375
sd(subset(cross_stats, strainA == "SI86" | strainB == "SI86", select = prop_reassortant, drop = T))
#0.2102989

#TX12
mean(subset(cross_stats, strainA == "TX12" | strainB == "TX12", select = prop_reassortant, drop = T))
#[1] 0.4192708
sd(subset(cross_stats, strainA == "TX12" | strainB == "TX12", select = prop_reassortant, drop = T))
#0.3128904

#PAN99
mean(subset(cross_stats, strainA == "PAN99" | strainB == "PAN99", select = prop_reassortant, drop = T))
#[1] 0.4921875
sd(subset(cross_stats, strainA == "PAN99" | strainB == "PAN99", select = prop_reassortant, drop = T))
#0.2761727

#HK68
mean(subset(cross_stats, strainA == "HK68" | strainB == "HK68", select = prop_reassortant, drop = T))
la#[1] 0.3515625
sd(subset(cross_stats, strainA == "HK68" | strainB == "HK68", select = prop_reassortant, drop = T))
#0.3261348

#CA09
mean(subset(cross_stats, strainA == "CA09" | strainB == "CA09", select = prop_reassortant, drop = T))
#[1] 0.5364583
sd(subset(cross_stats, strainA == "CA09" | strainB == "CA09", select = prop_reassortant, drop = T))
#0.3866347

#Alternative to heatmap, potential new Figure 3A
cross_stats_manual <- read.csv("cross_stats_manual.csv")
View(cross_stats_manual)
ggplot(cross_stats_manual, aes(x = strain, y = prop_reassortant)) + geom_point()
ggplot(cross_stats_manual, aes(x = strain, y = prop_reassortant)) + geom_boxplot() + geom_point(aes(colour = as.factor(prop_reassortant), size = 8)) + scale_color_brewer(palette = 'Set3') + xlab("Strain") + ylab("Proportion of Reassortant Plaque Isolates") + theme(legend.position = "none") + theme_tufte() + theme(text = element_text(size = 20,  family="Helvetica")) + theme(axis.text.x = element_text(size = 12)) + theme(legend.position = "none") + theme(panel.grid.major.y = element_line(color = "lightgray",size = 0.5))


#Let's do a quick proportion test on the reassortment data
prop.test(cross_stats$reassortants, cross_stats$clones)
#	10-sample test for equality of proportions without continuity correction

#data:  cross_stats$reassortants out of cross_stats$clones
#X-squared = 312.8, df = 9, p-value < 2.2e-16
#alternative hypothesis: two.sided
#sample estimates:
#  prop 1     prop 2     prop 3     prop 4     prop 5     prop 6     prop 7     prop 8     prop 9   prop 10
#0.16666667 0.43750000 0.73958333 0.61458333 0.92708333 0.63541667 0.04166667 0.78125000 0.13541667  0.41666667

pairwise.prop.test(cross_stats$reassortants, cross_stats$clones)
#Pairwise comparisons using Pairwise comparison of proportions 

#data:  cross_stats$reassortants out of cross_stats$clones 

#   1       2       3       4       5       6       7       8       9      
#2  0.00153 -       -       -       -       -       -       -       -      
#3  1.7e-13 0.00076 -       -       -       -       -       -       -      
#4  1.5e-08 0.18673 0.53732 -       -       -       -       -       -      
#5  < 2e-16 3.3e-11 0.01591 1.6e-05 -       -       -       -       -      
#6  2.8e-09 0.11943 0.80562 1.00000 5.6e-05 -       -       -       -      
#7  0.11943 1.1e-08 < 2e-16 3.9e-15 < 2e-16 5.2e-16 -       -       -      
#8  2.0e-15 5.3e-05 1.00000 0.18392 0.10997 0.31207 < 2e-16 -       -      
#9  1.00000 0.00017 4.0e-15 6.1e-10 < 2e-16 1.0e-10 0.31207 < 2e-16 -      
#10 0.00443 1.00000 0.00024 0.11943 5.5e-12 0.05766 5.0e-08 1.4e-05 0.00054

#P value adjustment method: holm

#Only 1x7, 1X9, 2x4, 2x6, 2x10, 3x4, 3x6, 3x8, 4x6, 4x8, 4x10, 5x8, 6x8, 6x10, 7x9 = 15 out of 44 
select(cross_stats, c(prop_reassortant, strainAB))
#1           0.167  HK68_TX12 
#2           0.438  CA09_PAN99
#3           0.740  CA09_TX12 
#4           0.615  PAN99_SI86
#5           0.927  CA09_SI86 
#6           0.635  SI86_TX12 
#7           0.0417 CA09_HK68 
#8           0.781  HK68_PAN99
#9           0.135  PAN99_TX12
#10           0.417  HK68_SI86

#Bring up ANOVA analysis, buried below
#Import matrix with crosses and strains for an ANOVA
cross_data_runA_matrix <- read.csv("~/Dropbox/mixtup/Documentos/ucdavis/papers/influenza_GbBSeq_human/influenza_GbBSeq/data/cross_data_runA_matrix.csv")

cross_matrix <- as.matrix(cross_data_runA_matrix[2:6])

#Now add strain information to cross_stats
cross_stats <- right_join(cross_stats, unique(subset(controls_df, select = c(cross, strainA, strainB))))

table(cross_stats$strainA, cross_stats$prop_reassortant)

model5 = lm(prop_reassortant ~ cross_matrix-1, cross_stats)
#Adding -1 to get all coefficients (assuming intercept doesn't make sense because there is no reassortment without strains)
summary(model5)
#Call:
#  lm(formula = prop_reassortant ~ cross_matrix - 1, data = cross_stats)

#Residuals:
#  1        2        3        4        5        6        7        8        9       10 
#-0.18750  0.07986 -0.12500  0.09375  0.23264  0.11111  0.02778 -0.07986 -0.21875  0.06597 

#Coefficients:
#                   Estimate Std. Error t value Pr(>|t|)   
#cross_matrixCA09  0.008681   0.105766   0.082  0.93777   
#cross_matrixHK68  0.348958   0.105766   3.299  0.02149 * 
# cross_matrixPAN99 0.515625   0.105766   4.875  0.00457 **
# cross_matrixSI86  0.345486   0.105766   3.267  0.02228 * 
# cross_matrixTX12  0.005208   0.105766   0.049  0.96263   
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.1958 on 5 degrees of freedom
#Multiple R-squared:  0.9403,	Adjusted R-squared:  0.8806 
#F-statistic: 15.75 on 5 and 5 DF,  p-value: 0.004435


#ANOVA Table
#Bring plots together
#Requires package GridExtra

summary(model5)
as.data.frame(exp(coef(model5))) #exp doesn't make sense here
as.data.frame(coef(model5))

#subset(summary(model5)$coef, "Pr(>|z|)" < 0.05)

model5_coef <- data.frame(summary(model5)$coef)
colnames(model5_coef) <- c("Coef", "Std Error", "t value", "p")
model5_coef

model5_coef_rtable <- tableGrob(round(model5_coef, digits=3))
grid.newpage()
grid.draw(model5_coef_rtable)

grid.arrange(model5_coef)

#Now lets look at H1N1 vs H3N2
#Plot by subtype
ggplot(controls_df, aes(x = subtype, y = prop_reassortant)) + geom_point() + geom_boxplot(alpha = 0.25) 

#Yes or no? This is Figure 4
ggplot(controls_df, aes(x = same_subtype, y = prop_reassortant)) + xlab("Are coinfecting strains the same subtype?") + ylab("Proportion of Reassortant Plaque Isolates") + geom_point(aes(color=cross, size = 7)) + geom_boxplot(alpha = 0.25) + theme(legend.position = "none") + theme_tufte() +  theme(text = element_text(size = 20,  family="Helvetica")) + theme(axis.text.x = element_text(size = 12)) + theme(legend.position = "none") + theme(panel.grid.major.y = element_line(color = "lightgray",size = 0.5))  

summary.lm(aov(prop_reassortant ~ subtype, controls_df))
#Call:
#aov(formula = prop_reassortant ~ subtype, data = controls_df)

#Residuals:
#  Min      1Q  Median      3Q     Max 
#-0.4485 -0.1984  0.0000  0.1452  0.4161 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)       0.927083   0.009296   99.73   <2e-16 ***
#  subtypeH3N2      -0.561963   0.010641  -52.81   <2e-16 ***
#  subtypeH3N2/H1NI -0.436892   0.010020  -43.60   <2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.2404 on 6953 degrees of freedom
#Multiple R-squared:  0.2864,	Adjusted R-squared:  0.2862 
#F-statistic:  1395 on 2 and 6953 DF,  p-value: < 2.2e-16

#This should be a t-test
summary.lm(aov(prop_reassortant ~ same_subtype, controls_df))
#Call:
#aov(formula = prop_reassortant ~ same_subtype, data = controls_df)

#Residuals:
#  Min      1Q  Median      3Q     Max 
#-0.4485 -0.3316  0.1244  0.2494  0.4288 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)   0.490191   0.004427  110.73   <2e-16 ***
#  same_subtypeY 0.008057   0.006948    1.16    0.246    
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.2846 on 6954 degrees of freedom
#Multiple R-squared:  0.0001933,	Adjusted R-squared:  4.956e-05 
#F-statistic: 1.345 on 1 and 6954 DF,  p-value: 0.2462

#Now t-test
intrasubtype <- as.vector(unique(subset(controls_df, same_subtype == "Y", select = "prop_reassortant")))
intersubtype <- as.vector(unique(subset(controls_df, same_subtype == "N", select = "prop_reassortant")))

t.test(intrasubtype, intersubtype)
#Welch Two Sample t-test

#data:  intrasubtype and intersubtype
#t = 0.094822, df = 4.4789, p-value = 0.9285
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.5878083  0.6312111
#sample estimates:
#  mean of x mean of y 
#0.5026042 0.480902

mean(intrasubtype$prop_reassortant)
#[1] 0.5026042
sd(intrasubtype$prop_reassortant)
#[1] 0.4104902

mean(intersubtype$prop_reassortant)
#[1] 0.4809028
sd(intersubtype$prop_reassortant)
#[1] 0.2480319

#Now lets look at genetic similarity, pairwise genetic identity (PID)
pairwise_genetic_identity <- read.csv("~/Dropbox/mixtup/Documentos/ucdavis/papers/influenza_GbBSeq_human/influenza_GbBSeq/data/pairwise_genetic_identity.csv", na.strings="")

#Make a strainAB column
pairwise_genetic_identity <- unite(pairwise_genetic_identity, strainAB, c("strainB", "strainA",), sep = "_", remove = F)
#Note need to flip order of strainB and strainA to match cross_stats

#Let's just calculate some basic stats on pairwise genetic identity

#First average percent identity
average_pid <- group_by(pairwise_genetic_identity, strainAB) %>%
  summarise(count = n(),
            average_pid = mean(pairwise_id)
            )
# A tibble: 10 × 3
#strainAB   count average_pid
#<chr>      <int>       <dbl>
#1 CA09_HK68      8        78.1
#2 CA09_PAN99     8        77.6
#3 CA09_SI86      8        82.8
#4 CA09_TX12      8        77.0
#5 HK68_PAN99     8        93.1
#6 HK68_SI86      8        82.1
#7 HK68_TX12      8        91.3
#8 PAN99_SI86     8        80.4
#9 PAN99_TX12     8        96.4
#10 SI86_TX12      8        81.0

#Second, average percent identity 
average_pid_no_antigenic <- group_by(subset(pairwise_genetic_identity, segment != "HA" & segment != "NA"), strainAB) %>%
  summarise(count = n(),
            average_pid_no_antigenic = mean(pairwise_id)
  )
# A tibble: 10 × 3
#strainAB   count average_pid_no_antigenic
#<chr>      <int>                    <dbl>
#1 CA09_HK68      6                     86.6
#2 CA09_PAN99     6                     85.6
#3 CA09_SI86      6                     84.2
#4 CA09_TX12      6                     85.0
#5 HK68_PAN99     6                     94.2
#6 HK68_SI86      6                     92.1
#7 HK68_TX12      6                     92.8
#8 PAN99_SI86     6                     89.3
#9 PAN99_TX12     6                     97.0
#10 SI86_TX12      6                     88.4

average_pid <- cbind(average_pid, average_pid_no_antigenic$average_pid_no_antigenic)
average_pid <- average_pid[, c(1,3,4)]
colnames(average_pid) <- c("strainAB", "average_pid", "average_pid_no_antigenic")

cross_stats <- right_join(cross_stats, average_pid)

#Now let's test for correlation between PID and reassortment rate
ggplot(cross_stats, aes(x = prop_reassortant, y = average_pid)) + geom_point(aes(colour = strainAB), size = 6) + geom_smooth(method = "lm") + theme(text = element_text(size = 20,  family="Helvetica")) + theme(axis.text.x = element_text(size = 15)) + theme(panel.grid.major.y = element_line(color = "lightgray",size = 0.5))       

reassortment_pid_model <- lm(prop_reassortant ~ average_pid, cross_stats)
summary(reassortment_pid_model)
#Call:
#  lm(formula = prop_reassortant ~ average_pid, data = cross_stats)

#Residuals:
#  Min       1Q   Median       3Q      Max 
#-0.50354 -0.20585  0.00032  0.16816  0.42617 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)  1.278944   1.248935   1.024    0.336
#average_pid -0.009399   0.014826  -0.634    0.544

#Residual standard error: 0.3113 on 8 degrees of freedom
#Multiple R-squared:  0.04784,	Adjusted R-squared:  -0.07118 
#F-statistic: 0.402 on 1 and 8 DF,  p-value: 0.5438

#No significant relationship! 

#Same but excluding antigenic segments
ggplot(cross_stats, aes(x = prop_reassortant, y = average_pid_no_antigenic)) + geom_point(aes(colour = strainAB), size = 6) + geom_smooth(method = "lm") + theme(text = element_text(size = 20,  family="Helvetica")) + theme(axis.text.x = element_text(size = 15)) + theme(panel.grid.major.y = element_line(color = "lightgray",size = 0.5))
reassortment_pid_model_no_antigenic <- lm(prop_reassortant ~ average_pid_no_antigenic, cross_stats)
summary(reassortment_pid_model_no_antigenic)
#Call:
#  lm(formula = prop_reassortant ~ average_pid_no_antigenic, data = cross_stats)
#
#Residuals:
#  Min       1Q   Median       3Q      Max 
#-0.52866 -0.15740  0.05647  0.12578  0.42021 
#
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)               2.94342    2.01936   1.458    0.183
#average_pid_no_antigenic -0.02742    0.02254  -1.216    0.258

#Residual standard error: 0.2931 on 8 degrees of freedom
#Multiple R-squared:  0.1561,	Adjusted R-squared:  0.05061 
#F-statistic:  1.48 on 1 and 8 DF,  p-value: 0.2585

#No significant relationship, either!

#Finally lets look at other measures of reassortment
#Now let's test for correlation between PID and proportion unique genotypes
ggplot(cross_stats, aes(x = proportion_unique_genotypes, y = average_pid)) + geom_point(aes(colour = strainAB)) + geom_smooth(method = "lm") 
unique_genotypes_pid <- lm(proportion_unique_genotypes ~ average_pid, cross_stats)
summary(unique_genotypes_pid)
#Multiple R-squared:  0.1165,	Adjusted R-squared:  0.00603 
#F-statistic: 1.055 on 1 and 8 DF,  p-value: 0.3345

#Not surprising as reassortment rate and unique genotypes are correlated

#Now let's test for correlation between PID and pairwise_heterologous_mean_frequency
ggplot(cross_stats, aes(x = pairwise_heterologous_mean_frequency, y = average_pid)) + geom_point(aes(colour = strainAB)) + geom_smooth(method = "lm") 
pairwise_heterologous_mean_frequency_pid <- lm(pairwise_heterologous_mean_frequency ~ average_pid, cross_stats)
summary(pairwise_heterologous_mean_frequency_pid)
#Multiple R-squared:  0.1561,	Adjusted R-squared:  0.05056 
#F-statistic: 1.479 on 1 and 8 DF,  p-value: 0.2586

#Also uncorrelated!

#Now let's look at the segments involved in reasortment

#What is the bias in segment assignment per locus?
locus_strain <- group_by(controls_df, cross, locus, majority_strain) %>%
  summarise(
    total_samples = length(unique(sample)),
  )

#Same column heights to compare
ggplot(locus_strain, aes(x = locus, y = total_samples, fill = majority_strain)) + geom_col(position = "fill") + facet_wrap(~cross)
#This graph could be really useful for the single cycle conditions experiments

#Now let's rearrange in pairwise combination
#First need to add strainA/strainB information
subset(controls_df, select = c("cross", "strainA", "strainB", "strainAB"))
locus_strain <- right_join(locus_strain, subset(controls_df, select = c("cross", "strainA", "strainB", "strainAB")))
#Plot! This is Figure 6
ggplot(locus_strain, aes(x = factor(locus, level = c('PB2f', 'PB1c', 'PAc', 'HA', 'NPd', 'NA', 'Mg', 'NS1d')), y = total_samples, fill = majority_strain)) + geom_col(position = "fill") + facet_grid(strainA~strainB) + xlab("Segment") + ylab("Proportion of Plaque Isolates Assigned to Each Strain") + theme(axis.text.x = element_text(angle=90)) + geom_hline(yintercept = 0.50, alpha = 0.75, linetype = 2) + theme_tufte() +  theme(text = element_text(size = 20,  family="Helvetica")) + theme(axis.text.x = element_text(angle=90, size = 17)) + theme(strip.text = element_text(face = "bold")) + theme(legend.position = "none")

#To address Anice Lowen comment on parental strain reproduction (solo infections), going to restrict only to reassortant genotypes
controls_df_reassortant <- right_join(controls_df, subset(controls_df_number_parents, select = c("cross_sample", "reassortant")), by="cross_sample")
controls_df_reassortants_only  <- subset(controls_df_reassortant, reassortant == 1)
locus_strain_reassortants_only <- group_by(controls_df_reassortants_only, cross, locus, majority_strain) %>%
  summarise(
    total_samples = length(unique(sample)),
  )
locus_strain_reassortants_only <- right_join(locus_strain_reassortants_only, subset(controls_df_reassortants_only, select = c("cross", "strainA", "strainB", "strainAB")))
#Plot! This is Figure 6, but with reassortants only
ggplot(locus_strain_reassortants_only, aes(x = factor(locus, level = c('PB2f', 'PB1c', 'PAc', 'HA', 'NPd', 'NA', 'Mg', 'NS1d')), y = total_samples, fill = majority_strain)) + geom_col(position = "fill") + facet_grid(strainA~strainB) + xlab("Segment") + ylab("Proportion of Plaque Isolates Assigned to Each Strain") + theme(axis.text.x = element_text(angle=90)) + geom_hline(yintercept = 0.50, alpha = 0.75, linetype = 2)



#Putting in one row
ggplot(locus_strain, aes(x = locus, y = total_samples, fill = majority_strain)) + geom_col(position = "fill") + facet_wrap(~cross, nrow = 1) + geom_hline(yintercept = 0.50)

#
ggplot(locus_strain, aes(x = locus, y = total_samples, fill = majority_strain)) + geom_col(position = "fill") + coord_flip() + facet_wrap(~cross, ncol = 1) + geom_hline(yintercept = 0.50)

#Segments involved in reassortment
#Doing visually and manually here
segments_reassortment <- c(3,8,8,8,8,7,6,8,7,8)
mean(segments_reassortment)
# 7.1
sd(segments_reassortment)
#1.595131

#So first thing I have to do here is calculate the expected 95% CI for binomial. This treats each segment as independent
#I can plot the 95% CI and draw lines on the graph

#So in this case we have 96 samples. What is the confidence interval around a 50-50 match
binom.test(c(48, 48), p = 0.5)
#Exact binomial test
#data:  c(48, 48)
#number of successes = 48, number of trials = 96, p-value = 1
#alternative hypothesis: true probability of success is not equal to 0.5
#95 percent confidence interval:
#  0.3961779 0.6038221
#sample estimates:
#  probability of success 
#   0.5

#So then our 95% confidence interval (for 96 trials) is 0.3961779 0.6038221, 
#Redraw with these limits
ggplot(locus_strain, aes(x = locus, y = total_samples, fill = majority_strain)) + geom_col(position = "fill") + facet_grid(strainA~strainB) + geom_hline(yintercept = 0.50, alpha = 0.75, linetype = 2) + geom_hline(yintercept = c(0.396, 0.604), alpha = 0.75, linetype = 2, color = "grey")

#But I don't think I have equal samples, per cross, argh. Need to adjust, per sample size 



get_CI <- function(sample_size) {
  half_integer <- round(sample_size/2, digits = 0)
  ci_low <- binom.test(c(half_integer, half_integer), p = 0.5)[4]$conf.int[1]
  ci_low <- as.numeric(ci_low)
  ci_high<- binom.test(c(half_integer, half_integer), p = 0.5)[4]$conf.int[2]
  ci_high <- as.numeric(ci_high)
  ci_df <- data.frame(ci_low, ci_high)
  return(ci_df)
}

get_low_CI <- function(sample_size) {
  half_integer <- round(sample_size/2, digits = 0)
  ci_low <- binom.test(c(half_integer, half_integer), p = 0.5)[4]$conf.int[1]
  ci_low <- as.numeric(ci_low)
  return(ci_low)
}

get_high_CI <- function(sample_size) {
  half_integer <- round(sample_size/2, digits = 0)
  ci_high <- binom.test(c(half_integer, half_integer), p = 0.5)[4]$conf.int[2]
  ci_high <- as.numeric(ci_high)
  return(ci_high)
}

#### All Data ####
#Need to add the proportion we see in the plot, but I have questions about how I originally made locus_strain (Line 399), so going with this
binomial_segment <- group_by(controls_df, cross, locus, majority_strain) %>%
  summarise(count = n()) %>% 
    group_by(cross, locus) %>% 
      summarise(proportion = max(count) / sum(count), total = sum(count))
#Here we go!
binomial_segment <- data.frame(binomial_segment, ci_low = sapply(binomial_segment$total, get_low_CI), ci_high = sapply(binomial_segment$total, get_high_CI))

#Now add strains for plotting 
unique(subset(controls_df, select = c("cross", "locus", "strainA", "strainB", "strainAB")))
binomial_segment <- right_join(binomial_segment, unique(subset(controls_df, select = c("cross", "locus", "strainA", "strainB", "strainAB"))))

#Last thing!
#Let's make a column for the ones that fall within the confidence interval
inside_ci <- binomial_segment$proportion >= binomial_segment$ci_low & binomial_segment$proportion <= binomial_segment$ci_high

binomial_segment <- data.frame(binomial_segment, inside_ci)

#Slightly Different plot than above, but very nice as well 
ggplot(binomial_segment, aes(x = locus, y = proportion, colour = locus)) + geom_point(aes(shape = inside_ci), size = 2) + 
  geom_point(colour = "grey90", size = 0.5) + geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.1) + facet_grid(strainA~strainB) + geom_hline(yintercept = 0.50, alpha = 0.75, linetype = 2)
#This one does a better job of highlighting discordant loci. This is Figure 6B
ggplot(binomial_segment, aes(x = factor(locus, level = c('PB2f', 'PB1c', 'PAc', 'HA', 'NPd', 'NA', 'Mg', 'NS1d')), y = proportion, colour = inside_ci)) + geom_point(aes(shape = inside_ci), size = 2) + theme(axis.text.x = element_text(angle=90)) + ylab("Proportion of Plaque Isolates Assigned to Each Strain") + xlab("Segment") +
  geom_point(colour = "grey90", size = 0.5) + geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.1) + facet_grid(strainA~strainB) + geom_hline(yintercept = 0.50, alpha = 0.75, linetype = 2) + theme_tufte() +  theme(text = element_text(size = 20,  family="Helvetica")) + theme(axis.text.x = element_text(angle=90, size = 17)) + theme(strip.text = element_text(face = "bold")) + theme(legend.position = "none")


#Doesn't really look like the ones that fall inside the CI are of a specific locus or strain combination
ggplot(binomial_segment, aes(x = inside_ci, y = proportion)) + geom_boxplot() + geom_point(aes(colour = strainAB)) + facet_wrap(~ locus) 

#Some summary stats for paper and we should be good!

#Take a look by locus
binomial_segment %>% group_by(locus) %>% summarise(
  count = n(),
  mean_proportion_max = mean(proportion),
  sd_proportion_max = sd(proportion),
  number_inside_ci = length(which(inside_ci==TRUE)),
  )
# A tibble: 8 x 5
#locus count mean_proportion_max sd_proportion_max number_inside_ci
#<chr> <int>               <dbl>             <dbl>            <int>
#1 HA       10               0.841             0.156                1
#2 Mg       10               0.864             0.156                1
#3 NA       10               0.829             0.193                3
#4 NPd      10               0.866             0.125                1
#5 NS1d     10               0.803             0.186                3
#6 PAc      10               0.779             0.194                2
#7 PB1c     10               0.827             0.137                2
#8 PB2f     10               0.865             0.111                0

#Not much going on here, tends to be a pretty good bias towards one strain or the other

#Now by cross
binomial_segment %>% group_by(cross, strainAB) %>% summarise(
  count = n(),
  mean_proportion = mean(proportion),
  sd_proportion = sd(proportion),
  number_inside_ci = length(which(inside_ci==TRUE)),
)
# # A tibble: 10 x 6
# Groups:   cross [10]
# cross   strainAB   count mean_proportion sd_proportion number_inside_ci
#1 cross1  HK68_TX12      8           0.967        0.0289                0
#2 cross11 CA09_PAN99     8           0.888        0.147                 1
#3 cross14 CA09_TX12      8           0.708        0.138                 2
#4 cross15 PAN99_SI86     8           0.615        0.0839                6
#5 cross18 CA09_SI86      8           0.727        0.118                 1
#6 cross19 SI86_TX12      8           0.825        0.131                 1
#7 cross2  CA09_HK68      8           0.993        0.0115                0
#8 cross3  HK68_PAN99     8           0.757        0.131                 2
#9 cross7  PAN99_TX12     8           0.968        0.0221                0
#10 cross8  HK68_SI86      8           0.895        0.0778               0

#Okay there is something interesting going on here, which partially supports strain specificity,
# PAN99_SI86 (0.61458333) has 6(!) segments that show free assortment! There's two strains that have 2 free segments.
# These are CA09_TX12 (0.73958333) and HK68_PAN99 (0.7812500). All three are relatively high reassortment, but the highest reassortment coinfection,
# has only 1 conforming

#Any relationship?
segment_assortment_reasortment <- right_join(binomial_segment %>% group_by(cross, strainAB) %>% summarise(
  count = n(),
  mean_proportion = mean(proportion),
  sd_proportion = sd(proportion),
  number_inside_ci = length(which(inside_ci==TRUE)),
), cross_stats)

#No relationship... Would have to be audacious to draw a straight line through that... (but R will do it!)
# This should be a Supplementary figure
ggplot(segment_assortment_reasortment, aes(x = prop_reassortant, y = number_inside_ci)) + geom_point(aes(colour = strainAB), size = 6)

#### END All Data Segment representation ####


#### Reassortants only ####
#Again this is Figure 6, but with reassortants only
ggplot(locus_strain_reassortants_only, aes(x = factor(locus, level = c('PB2f', 'PB1c', 'PAc', 'HA', 'NPd', 'NA', 'Mg', 'NS1d')), y = total_samples, fill = majority_strain)) + geom_col(position = "fill") + facet_grid(strainA~strainB) + xlab("Segment") + ylab("Proportion of Plaque Isolates Assigned to Each Strain") + theme(axis.text.x = element_text(angle=90)) + geom_hline(yintercept = 0.50, alpha = 0.75, linetype = 2)

#Putting in one row
ggplot(locus_strain_reassortants_only, aes(x = locus, y = total_samples, fill = majority_strain)) + geom_col(position = "fill") + facet_wrap(~cross, nrow = 1) + geom_hline(yintercept = 0.50)

#
ggplot(locus_strain_reassortants_only, aes(x = locus, y = total_samples, fill = majority_strain)) + geom_col(position = "fill") + coord_flip() + facet_wrap(~cross, ncol = 1) + geom_hline(yintercept = 0.50)

#Segments involved in reassortment
#Doing visually and manually here
segments_reassortment <- c(3,8,8,8,8,7,6,8,7,8)
mean(segments_reassortment)
# 7.1
sd(segments_reassortment)
#1.595131

#So first thing I have to do here is calculate the expected 95% CI for binomial. This treats each segment as independent
#I can plot the 95% CI and draw lines on the graph

#So in this case we have 96 samples. What is the confidence interval around a 50-50 match
binom.test(c(48, 48), p = 0.5)
#Exact binomial test
#data:  c(48, 48)
#number of successes = 48, number of trials = 96, p-value = 1
#alternative hypothesis: true probability of success is not equal to 0.5
#95 percent confidence interval:
#  0.3961779 0.6038221
#sample estimates:
#  probability of success 
#   0.5

#So then our 95% confidence interval (for 96 trials) is 0.3961779 0.6038221, 
#Redraw with these limits
ggplot(locus_strain_reassortants_only, aes(x = locus, y = total_samples, fill = majority_strain)) + geom_col(position = "fill") + facet_grid(strainA~strainB) + geom_hline(yintercept = 0.50, alpha = 0.75, linetype = 2) + geom_hline(yintercept = c(0.396, 0.604), alpha = 0.75, linetype = 2, color = "grey")
#But I don't think I have equal samples, per cross, argh. Need to adjust, per sample size 



get_CI <- function(sample_size) {
  half_integer <- round(sample_size/2, digits = 0)
  ci_low <- binom.test(c(half_integer, half_integer), p = 0.5)[4]$conf.int[1]
  ci_low <- as.numeric(ci_low)
  ci_high<- binom.test(c(half_integer, half_integer), p = 0.5)[4]$conf.int[2]
  ci_high <- as.numeric(ci_high)
  ci_df <- data.frame(ci_low, ci_high)
  return(ci_df)
}

get_low_CI <- function(sample_size) {
  half_integer <- round(sample_size/2, digits = 0)
  ci_low <- binom.test(c(half_integer, half_integer), p = 0.5)[4]$conf.int[1]
  ci_low <- as.numeric(ci_low)
  return(ci_low)
}

get_high_CI <- function(sample_size) {
  half_integer <- round(sample_size/2, digits = 0)
  ci_high <- binom.test(c(half_integer, half_integer), p = 0.5)[4]$conf.int[2]
  ci_high <- as.numeric(ci_high)
  return(ci_high)
}

#Need to add the proportion we see in the plot, but I have questions about how I originally made locus_strain (Line 399), so going with this
binomial_segment_reassortants_only <- group_by(controls_df_reassortants_only, cross, locus, majority_strain) %>%
  summarise(count = n()) %>% 
  group_by(cross, locus) %>% 
  summarise(proportion = max(count) / sum(count), total = sum(count))
#Here we go!
binomial_segment_reassortants_only <- data.frame(binomial_segment_reassortants_only, ci_low = sapply(binomial_segment$total, get_low_CI), ci_high = sapply(binomial_segment$total, get_high_CI))

#Now add strains for plotting 
unique(subset(controls_df_reassortants_only, select = c("cross", "locus", "strainA", "strainB", "strainAB")))
binomial_segment_reassortants_only <- right_join(binomial_segment_reassortants_only, unique(subset(controls_df_reassortants_only, select = c("cross", "locus", "strainA", "strainB", "strainAB"))))

#Last thing!
#Let's make a column for the ones that fall within the confidence interval
inside_ci_reassortants_only <- binomial_segment_reassortants_only$proportion >= binomial_segment_reassortants_only$ci_low & binomial_segment_reassortants_only$proportion <= binomial_segment_reassortants_only$ci_high

binomial_segment_reassortants_only <- data.frame(binomial_segment_reassortants_only, inside_ci_reassortants_only)

#Slightly Different plot than above, but very nice as well 
ggplot(binomial_segment_reassortants_only, aes(x = locus, y = proportion, colour = locus)) + geom_point(aes(shape = inside_ci_reassortants_only), size = 2) + 
  geom_point(colour = "grey90", size = 0.5) + geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.1) + facet_grid(strainA~strainB) + geom_hline(yintercept = 0.50, alpha = 0.75, linetype = 2)
#This one does a better job of highlighting discordant loci. This is Figure 6B
ggplot(binomial_segment_reassortants_only, aes(x = factor(locus, level = c('PB2f', 'PB1c', 'PAc', 'HA', 'NPd', 'NA', 'Mg', 'NS1d')), y = proportion, colour = inside_ci_reassortants_only)) + geom_point(aes(shape = inside_ci), size = 2) + theme(axis.text.x = element_text(angle=90)) + ylab("Proportion of Plaque Isolates Assigned to Each Strain") + xlab("Segment") +
  geom_point(colour = "grey90", size = 0.5) + geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.1) + facet_grid(strainA~strainB) + geom_hline(yintercept = 0.50, alpha = 0.75, linetype = 2)
#The pattern remains pretty similar using reassortant data only, some blinking in and out

#Doesn't really look like the ones that fall inside the CI are of a specific locus or strain combination
ggplot(binomial_segment_reassortants_only, aes(x = inside_ci_reassortants_only, y = proportion)) + geom_boxplot() + geom_point(aes(colour = strainAB)) + facet_wrap(~ locus) 

#Some summary stats for paper and we should be good!

#Take a look by locus
binomial_segment_reassortants_only %>% group_by(locus) %>% summarise(
  count = n(),
  mean_proportion_max = mean(proportion),
  sd_proportion_max = sd(proportion),
  number_inside_ci = length(which(inside_ci_reassortants_only==TRUE)),
)
# A tibble: 8 × 5
#locus count mean_proportion_max sd_proportion_max number_inside_ci
#<chr> <int>               <dbl>             <dbl>            <int>
#1 HA       10               0.787             0.188                2
#2 Mg       10               0.815             0.152                1
#3 NA       10               0.755             0.195                3
#4 NPd      10               0.767             0.150                2
#5 NS1d     10               0.793             0.139                1
#6 PAc      10               0.744             0.170                2
#7 PB1c     10               0.793             0.180                2
#8 PB2f     10               0.766             0.109                1

#Pretty similar to all data and reassortants only 

#Now by cross
binomial_segment_reassortants_only %>% group_by(cross, strainAB) %>% summarise(
  count = n(),
  mean_proportion = mean(proportion),
  sd_proportion = sd(proportion),
  number_inside_ci = length(which(inside_ci_reassortants_only==TRUE)),
)
# A tibble: 10 × 6
# Groups:   cross [10]
#cross   strainAB   count mean_proportion sd_proportion number_inside_ci
#<chr>   <chr>      <int>           <dbl>         <dbl>            <int>
#1 cross1  HK68_TX12      8           0.875         0.111                0
#2 cross11 CA09_PAN99     8           0.858         0.146                1
#3 cross14 CA09_TX12      8           0.725         0.119                1
#4 cross15 PAN99_SI86     8           0.635         0.115                4
#5 cross18 CA09_SI86      8           0.743         0.129                1
#6 cross19 SI86_TX12      8           0.796         0.188                2
#7 cross2  CA09_HK68      8           0.896         0.146                0
#8 cross3  HK68_PAN99     8           0.712         0.162                3
#9 cross7  PAN99_TX12     8           0.765         0.151                1
#10 cross8  HK68_SI86      8           0.770         0.163                1

#Some changes compared to all data, some loci blinking in or out (due to reduced sample size), but overall no dramatic differences 

#Any relationship?
segment_assortment_reasortment_only <- right_join(binomial_segment_reassortants_only %>% group_by(cross, strainAB) %>% summarise(
  count = n(),
  mean_proportion = mean(proportion),
  sd_proportion = sd(proportion),
  number_inside_ci = length(which(inside_ci_reassortants_only==TRUE)),
), cross_stats)

#No relationship, either...
# This should be a Supplementary figure
ggplot(segment_assortment_reasortment_only, aes(x = prop_reassortant, y = number_inside_ci)) + geom_point(aes(colour = strainAB))

#### END Reassortants only ####


#### 2. Linkage of Segments ####
#Then we go to look at linkage disequilibrium 

controls_df_full <- filter(controls_df, total_segments == 8)

#Calculate the number of samples that have either strain at each locus
locus_strain_full <- group_by(controls_df_full, cross, locus, majority_strain) %>%
  summarise(
    sample_by_strain = length(unique(sample)),
  )

#Calculate total samples, per locus (for each cross)
locus_strain_full_totals <- group_by(locus_strain_full, cross, locus) %>%
  summarise(
    total_samples = sum(sample_by_strain),
  ) 
#COME BACK TO THIS - This is redundant

#Add totals to data frame
locus_strain_full <- right_join(locus_strain_full, locus_strain_full_totals)

#Add proprtions
locus_strain_full <- mutate(locus_strain_full, proportion = sample_by_strain/total_samples)

#First generate list of loci pairings
#New approach with no repeats from, except I used combinations(): 
#https://davetang.org/muse/2013/09/09/combinations-and-permutations-in-r/

library(gtools)
#urn with 3 balls
loci <- unique(locus_strain_full$locus)
combinations(n = 8,r = 2,v = loci)
nrow(combinations(n = 8,r = 2,v = loci))
#[1] 28

#Now implement into data frame
loci_pairings <- as.data.frame(combinations(n = 8,r = 2,v = loci))
#Rename columns
colnames(loci_pairings) <- c("locusA", "locusB")

loci_pairings[1, 2]

#FYI I should be able to automate this...
linkage_disequilibrium <- function(df, target_cross, loci, ref_strain) {
  cross_strain <- subset(df, cross == target_cross & majority_strain == ref_strain)
  #locusA <- loci[1]
  #locusB <- loci[2]
  pi <- cross_strain[cross_strain$locus == loci[1], "proportion"]
  pj <- cross_strain[cross_strain$locus == loci[2], "proportion"]
  pij <- (cross_strain[cross_strain$locus == loci[1], "sample_by_strain"] + cross_strain[cross_strain$locus == loci[2], "sample_by_strain"]) / (cross_strain[cross_strain$locus == loci[1], "total_samples"] + cross_strain[cross_strain$locus == loci[2], "total_samples"])
  d <- pij - (pi*pj) 
  return(d)
}


ld_df <- NULL
d_list <- NULL
#Run for loop
for (i in 1:nrow(loci_pairings)) {
  #print(unlist(loci_pairings[i, ], use.names=FALSE))
  d_list <- linkage_disequilibrium(locus_strain_full, "cross18", unlist(loci_pairings[i, ], use.names=FALSE), "CA09")
  #print(loci_pairings[i, ])
  #print(d_list)
  print(cbind(loci_pairings[i, ], d_list))
  #ld_df <- cbind(loci_pairings[i, ], d_list)
  ld_df <- rbind(ld_df, cbind(loci_pairings[i, ], d_list, cross = "cross18"))
}  

#make a locus pair column
ld_df <- unite(ld_df, locus_pair, c("locusA", "locusB"), sep = "_", remove = FALSE)

#Graph!
ggplot(ld_df, aes(reorder(locus_pair, sample_by_strain), sample_by_strain)) + geom_col(aes(fill=locus_pair), position = "dodge") + ylim(0, 0.3) + geom_hline(yintercept = 0, linetype='dotted') + theme(axis.text.x = element_text(angle = 90), legend.position = "none") + coord_flip() + facet_wrap(~cross)
#Why are all 28 not showing up...

#I think it's the limit. Why in the world did I put that in....
ggplot(ld_df, aes(reorder(locus_pair, sample_by_strain), sample_by_strain)) + geom_col(aes(fill=locus_pair), position = "dodge") + geom_hline(yintercept = 0, linetype='dotted') + theme(axis.text.x = element_text(angle = 90), legend.position = "none") + coord_flip() + facet_wrap(~cross)

#So cross18 is the highest reassortment cross CA09xSI86
#All values positive which means those combinations are over-represented

#Now let's try with HK68xSI86, cross 8, reference strain HK68
ld_df <- NULL
d_list <- NULL
#Run for loop
for (i in 1:nrow(loci_pairings)) {
  #print(unlist(loci_pairings[i, ], use.names=FALSE))
  d_list <- linkage_disequilibrium(locus_strain_full, "cross8", unlist(loci_pairings[i, ], use.names=FALSE), "HK68")
  #print(loci_pairings[i, ])
  #print(d_list)
  print(cbind(loci_pairings[i, ], d_list))
  #ld_df <- cbind(loci_pairings[i, ], d_list)
  ld_df <- rbind(ld_df, cbind(loci_pairings[i, ], d_list, cross = "cross8"))
}  

#make a locus pair column
ld_df <- unite(ld_df, locus_pair, c("locusA", "locusB"), sep = "_", remove = FALSE)

#Graph
ggplot(ld_df, aes(reorder(locus_pair, sample_by_strain), sample_by_strain)) + geom_col(aes(fill=locus_pair), position = "dodge") + geom_hline(yintercept = 0, linetype='dotted') + theme(axis.text.x = element_text(angle = 90), legend.position = "none") + coord_flip() + facet_wrap(~cross)

#Now let's try with PAN99xTX12, cross 7, reference strain PAN99
ld_df <- NULL
d_list <- NULL
#Run for loop
for (i in 1:nrow(loci_pairings)) {
  #print(unlist(loci_pairings[i, ], use.names=FALSE))
  d_list <- linkage_disequilibrium(locus_strain_full, "cross7", unlist(loci_pairings[i, ], use.names=FALSE), "PAN99")
  #print(loci_pairings[i, ])
  #print(d_list)
  print(cbind(loci_pairings[i, ], d_list))
  #ld_df <- cbind(loci_pairings[i, ], d_list)
  ld_df <- rbind(ld_df, cbind(loci_pairings[i, ], d_list, cross = "cross7"))
}  

#make a locus pair column
ld_df <- unite(ld_df, locus_pair, c("locusA", "locusB"), sep = "_", remove = FALSE)

#Graph
ggplot(ld_df, aes(reorder(locus_pair, sample_by_strain), sample_by_strain)) + geom_col(aes(fill=locus_pair), position = "dodge") + geom_hline(yintercept = 0, linetype='dotted') + theme(axis.text.x = element_text(angle = 90), legend.position = "none") + coord_flip() + facet_wrap(~cross)

#Ok, how do we scale this up so we can get LD values for each segment pair in each cross
#Make a second loop over that will go over the  

#Test loop
for (c in 1:nrow(cross_stats)) {
  print(cross_stats$cross[c])
  print(cross_stats$strainA[c])
}

ld_df <- NULL
d_list <- NULL
for (c in 1:nrow(cross_stats)) {
  print(cross_stats$cross[c])
  print(as.character(cross_stats$strainA[c]))
  #Run for loop
    for (i in 1:nrow(loci_pairings)) {
      #print(unlist(loci_pairings[i, ], use.names=FALSE))
      d_list <- linkage_disequilibrium(locus_strain_full, cross_stats$cross[c], unlist(loci_pairings[i, ], use.names=FALSE), as.character(cross_stats$strainA[c]))
      #print(loci_pairings[i, ])
      #print(d_list)
      #print(cbind(loci_pairings[i, ], d_list))
      #ld_df <- cbind(loci_pairings[i, ], d_list)
      ld_df <- rbind(ld_df, cbind(loci_pairings[i, ], d_list, cross_stats$cross[c]))
      print(ld_df)
    }
}

#START HERE: Aug 26, 2022
#Could just make each one individually (just ten) and then join them.... 
#Quicker probably

#I have evidence of a biased association through stats, now I want to know what that association is (which strain)
#Start with PB2f and PB1c, wanto to figure out how many samples are CA09, HK68, or mixed

#Now question is how do I scale up to arbitrary number of crosses and loci combinations
#The key is loci_by_sample

scale_up_test <- controls_df_full %>% 
  group_by(sample) %>% 
  mutate(loci_by_sample = paste0(majority_strain_locus, collapse = ""))

#Huh, this actually looks a lot like the genotypes table, maybe this is what I need. Let's take a peek
#genotypes_table  %>% group_by(cross) %>% 

#Let's actually look at the ld_df because going to want to do something like that 
ld_df <- NULL
d_list <- NULL
#Run for loop
for (i in 1:nrow(loci_pairings)) {
  #print(unlist(loci_pairings[i, ], use.names=FALSE))
  d_list <- linkage_disequilibrium(locus_strain_full, "cross18", unlist(loci_pairings[i, ], use.names=FALSE), "CA09")
  #print(loci_pairings[i, ])
  #print(d_list)
  print(cbind(loci_pairings[i, ], d_list))
  #ld_df <- cbind(loci_pairings[i, ], d_list)
  ld_df <- rbind(ld_df, cbind(loci_pairings[i, ], d_list, cross = "cross18"))
}  


#make a locus pair column
ld_df <- unite(ld_df, locus_pair, c("locusA", "locusB"), sep = "_", remove = FALSE)

#Wait, why isn't it just locus_strain_full??  
locus_strain_full 

locus_strain_full <- unite(locus_strain_full, locus_strain, c("locus", "majority_strain"), sep = "_", remove = F)

#Plot!
ggplot(locus_strain_full, aes(x = locus_strain, y= sample_by_strain)) + geom_col(aes(colour = locus)) + facet_wrap(~ cross, scales = "free_x")

#Ok, let's zoom in to get an idea of how we're doing this
ggplot(subset(locus_strain_full, cross == "cross15"), aes(x = locus_strain, y= sample_by_strain, fill = majority_strain)) + geom_col() + coord_flip() + facet_grid(~ locus)

#Ok, no I was wrong. Need to make a pairwise genotype for each strain.
# This is analogous to the 8 segment genotype, expect now there will be 28 loci per sample (for pairwise comparisons of loci)
# So need to unite in genotypes table across *only* the loci I'm pairing (say PB2f, PB1c) and then make that the title for the column
# Therefore, let's go back to the genotypes table 

#genotypes_table <- unite(genotypes_table, "HA_Mg",  c("HA", "Mg", "NA", "NPd", "NS1d", "PAc", "PB1c", "PB2f"), remove = F)

genotypes_table <- group_by(controls_df, cross_sample, strainAB, locus) %>%
  summarise(
    majority_strain = majority_strain,
  )

#Looking good!
genotypes_table <- spread(genotypes_table, locus, majority_strain)

#Merge all genotypes into one and make a code for each genotype
genotypes_table <- unite(genotypes_table, "genotype",  c("HA", "Mg", "NA", "NPd", "NS1d", "PAc", "PB1c", "PB2f"), remove = F)

genotypes_table <- separate(genotypes_table, cross_sample, into = c("cross", "sample"), sep = "_", remove = F)

genotypes_table_test <- unite(genotypes_table, "HA_Mg",  c("HA", "Mg"), remove = F)

ggplot(genotypes_table_test, aes(x = HA_Mg)) + geom_bar() + facet_wrap(~ cross, scales = "free_x")
#Wow, cool!!!

#Need to restrict to only complete genotypes
#genotypes_table_test <- na.omit(genotypes_table_test)

#Again
ggplot(genotypes_table_test, aes(x = HA_Mg)) + geom_bar() + facet_wrap(~ cross, scales = "free_x")

#With proportions
ggplot(genotypes_table_test, aes(x = HA_Mg)) + geom_bar(aes(y = ..prop.., group = 1)) + facet_wrap(~ cross, scales = "free_x")
#Much easier to see the differences 

#NEXT STEPS: April 4, 2022
#Now from here need to:
# 1. Make all the other pairings
# 2. Potentially summarize the pairings into a new tidy data frame
# 3. Then label by heterologous and homologous for coloring bars
# 4. Potentially can graph facets by cross or locus

#First let's see if we can automate the process. This is the key command
genotypes_table_test <- unite(genotypes_table, "HA_Mg",  c("HA", "Mg"), remove = F)


#Autmating below for each locus pairing
#Remeber to remove NA's at this step because 
genotypes_table_two_locus <- na.omit(genotypes_table)
for (i in 1:nrow(loci_pairings)) {
  a <- as.vector(loci_pairings$locusA[i])
  b <- as.vector(loci_pairings$locusB[i])
  ab <- paste(a, b, sep = "_")
  print(a)
  print(b)
  print(ab)
  
  #The !! (bang, bang) is crucial for ab variable to be evaluated
  genotypes_table_two_locus <- unite(genotypes_table_two_locus, !!ab, c(a, b), remove = F)
}

#Now we can use others!
ggplot(genotypes_table_two_locus, aes(x = HA_NA)) + geom_bar(aes(y = ..prop.., group = 1)) + facet_wrap(~ cross, scales = "free_x")

#Ok, all done with two-locus. Now need to summarize to be able to plot and make a table of all values
# Summary has to be at the cross level. Start with the big data frame 

genotypes_table <- group_by(controls_df, cross, strainA, strainB, strainAB, cross_sample, locus) %>%
  summarise(
    majority_strain = majority_strain,
  )

#Changing to pivot_wider() function, as spread() is deprecated
genotypes_table <- pivot_wider(genotypes_table, names_from = locus, values_from =  majority_strain)

#Merge all genotypes into one and make a code for each genotype
genotypes_table <- unite(genotypes_table, "genotype",  c("HA", "Mg", "NA", "NPd", "NS1d", "PAc", "PB1c", "PB2f"), remove = F)

#Drop NA's
genotypes_table_two_locus <- na.omit(genotypes_table)
#Now remove the 8-segment genotype column
genotypes_table_two_locus <-select(genotypes_table_two_locus, -genotype) 
for (i in 1:nrow(loci_pairings)) {
  a <- as.vector(loci_pairings$locusA[i])
  b <- as.vector(loci_pairings$locusB[i])
  ab <- paste(a, b, sep = "_")
  print(a)
  print(b)
  print(ab)
  
  #The !! (bang, bang) is crucial for ab variable to be evaluated
  genotypes_table_two_locus <- unite(genotypes_table_two_locus, !!ab, c(a, b), remove = F)
}

#Now let's get rid of the single loci
genotypes_table_two_locus <- select(genotypes_table_two_locus, -loci)

#NOw make long, count, and calculate proportions
genotypes_table_two_locus_long <- pivot_longer(genotypes_table_two_locus, 6:33, names_to = "multilocus", values_to = "genotype") %>%
  group_by(cross, strainA, strainB, strainAB, multilocus, genotype) %>%
  summarise(
    count = n(),
    ) %>% 
  mutate(freq = count / sum(count))

#Looks good! Now with frequency

#Now let's recreate the plot from above
ggplot(subset(genotypes_table_two_locus_long, multilocus == "HA_NA"), aes(x = genotype, y = count)) + geom_col() + facet_wrap(~ cross, scales = "free_x")

#Try with proportion 
ggplot(subset(genotypes_table_two_locus_long, multilocus == "HA_NA"), aes(x = genotype, y = freq)) + geom_col() + facet_wrap(~ cross, scales = "free_x")

#Now let's label with strainAB instead
ggplot(subset(genotypes_table_two_locus_long, multilocus == "HA_NA"), aes(x = genotype, y = freq)) + geom_col() + facet_wrap(~ strainAB, scales = "free_x")

#Now by locus if that makes sense??
ggplot(subset(genotypes_table_two_locus_long, cross == "cross8"), aes(x = genotype, y = freq)) + geom_col() + facet_wrap(~ multilocus, scales = "free_x")
#Actually this makes sense and may be the best way to present. We'll go with this.

#Last thing is to color by the locus
genotypes_table_two_locus_long <- separate(genotypes_table_two_locus_long, genotype, c("genotype_locusA", "genotype_locusB"), remove = F)

homologous <- ifelse(genotypes_table_two_locus_long$genotype_locusA == genotypes_table_two_locus_long$genotype_locusB, "Homologous", "Heterologous" )

genotypes_table_two_locus_long <- data.frame(genotypes_table_two_locus_long, homologous)

#Now let's plot again
ggplot(subset(genotypes_table_two_locus_long, cross == "cross1"), aes(x = genotype, y = freq, fill = homologous)) + geom_col() + scale_fill_manual(values=c("orange", "darkblue")) + facet_wrap(~ multilocus, scales = "free_x")

#This is really awesome! Now that I have the homologous/heterologous tag, I can make a summary plot, 
# then also recreate the target plots. 

#Summary plot across loci, per cross
ggplot(genotypes_table_two_locus_long, aes(x = homologous, fill = homologous)) + geom_bar() + scale_fill_manual(values=c("orange", "darkblue")) + facet_wrap(~ strainAB, scales = "free_x")

#Add colors!
ggplot(genotypes_table_two_locus_long, aes(x = homologous, fill = homologous)) + geom_bar() + scale_fill_manual(values=c("orange", "darkblue")) + facet_wrap(~ strainAB, scales = "free_x")
#Wow, this is really cool. Looks like there are actually few pairings where on average there is a bias towards homologous
#Actually this just says how many have heterologous and homologous pairings, but does not quanify by plaques! It's presence/absence
#So the way to phrase it is, there are few strain pairings where only homologous combinations 

#Last one, with the matrix to see strain-specific effects
ggplot(genotypes_table_two_locus_long, aes(x = homologous, fill = homologous)) + geom_bar() + scale_fill_manual(values=c("orange", "darkblue")) + facet_grid(strainA ~ strainB, scales = "free_x")
#Nah

#Let's settle of this one
ggplot(genotypes_table_two_locus_long, aes(x = homologous, y = count, fill = homologous)) + geom_col() + scale_fill_manual(values=c("orange", "darkblue")) + facet_wrap(~ strainAB, scales = "free_y")

#What about frequencies instead?


#Let's axis labels Figure 7A
ggplot(genotypes_table_two_locus_long, aes(x = homologous, y = count, fill = homologous)) + geom_col() + theme(axis.text.x = element_blank()) + scale_fill_manual(values=c("orange", "darkblue")) + ylab("Frequency of progeny plaque isolates") + xlab("Strain combinations of pairwise segment combinations") + facet_wrap(~ strainAB, scales = "free_y", ncol = 5)


#Now target plot
#First need to reshape the data frame so I have two columns
#Actually this giving me a bit of a headache, going to skip this for now.

#Huge plot!!!!!
#ggplot(genotypes_table_two_locus_long, aes(x = genotype, y = freq, fill = homologous)) + geom_col() + scale_fill_manual(values=c("orange", "darkblue")) + facet_wrap(strainAB ~ multilocus, scales = "free_x")
#ggplot(subset(genotypes_table_two_locus_long), aes(x = genotype, y = count)) + geom_col() + facet_wrap(multilocus ~ cross, scales = "free_x")
#R did it, but I couldn't see anything 
#Need to remove the single locus and genotype data, going back from top

#Ordered Grid #Huge plot!!!!!
#ggplot(genotypes_table_two_locus_long, aes(x = genotype, y = freq, fill = homologous)) + geom_col() + scale_fill_manual(values=c("orange", "darkblue")) + facet_grid(multilocus ~ strainAB, scales = "free_x")

#Now let's just subset just the most extreme crosses in terms of least evenly distributed, to most distributed
ggplot(subset(genotypes_table_two_locus_long, strainAB == "CA09_HK68"), aes(x = genotype, y = freq, fill = homologous)) + geom_col() + theme(axis.text.x = element_text(angle=90)) + scale_fill_manual(values=c("darkblue", "orange")) + facet_grid(strainAB ~ multilocus, scales = "free_x")

ggplot(subset(genotypes_table_two_locus_long, strainAB == "CA09_SI86"), aes(x = genotype, y = freq, fill = homologous)) + geom_col() + theme(axis.text.x = element_text(angle=90)) + scale_fill_manual(values=c("orange", "darkblue")) + facet_grid(strainAB ~ multilocus, scales = "free_x")

ggplot(subset(genotypes_table_two_locus_long, strainAB == "CA09_SI86"), aes(x = genotype, y = count, fill = homologous)) + geom_col() + theme(axis.text.x = element_text(angle=90)) + scale_fill_manual(values=c("orange", "darkblue")) + facet_wrap( ~ multilocus, scales ="free")

ggplot(subset(genotypes_table_two_locus_long, strainAB == "CA09_PAN99"), aes(x = genotype, y = count, fill = homologous)) + geom_col() + theme(axis.text.x = element_text(angle=90)) + scale_fill_manual(values=c("orange", "darkblue")) + facet_wrap( ~ multilocus, scales ="free")

ggplot(subset(genotypes_table_two_locus_long, strainAB == "CA09_PAN99"), aes(x = genotype, y = freq, fill = homologous)) + geom_col() + theme(axis.text.x = element_text(angle=90)) + scale_fill_manual(values=c("orange", "darkblue")) + facet_wrap( ~ multilocus, scales = "free_x")

#HA and NA for all coinfections 
ggplot(subset(genotypes_table_two_locus_long, multilocus == "HA_NA"), aes(x = genotype, y = freq, fill = homologous)) + geom_col() + scale_fill_manual(values=c("orange", "darkblue")) + facet_wrap(strainA ~ strainB, scales = "free")

#PB2 and PA for all coinfections
ggplot(subset(genotypes_table_two_locus_long, multilocus == "PAc_PB2f"), aes(x = genotype, y = freq, fill = homologous)) + geom_col() + scale_fill_manual(values=c("orange", "darkblue")) + facet_wrap(strainA ~ strainB, scales = "free")

#PB2 and PB1 for all coinfections
ggplot(subset(genotypes_table_two_locus_long, multilocus == "PB1c_PB2f"), aes(x = genotype, y = freq, fill = homologous)) + geom_col() + scale_fill_manual(values=c("orange", "darkblue")) + facet_wrap(strainA ~ strainB, scales = "free")

#PA and PB1 for all coinfections
ggplot(subset(genotypes_table_two_locus_long, multilocus == "PAc_PB1c"), aes(x = genotype, y = freq, fill = homologous)) + geom_col() + scale_fill_manual(values=c("orange", "darkblue")) + facet_wrap(strainA ~ strainB, scales = "free")



#So the way to make the point is to look at whether there are loci that have majority heterologous plaques
#The real question is how many loci per cross generated 4, 3, 2, 1 loci. The majority of loci across crosses generated all segment combinations
pairwise_homologous_by_locus <- genotypes_table_two_locus_long %>% group_by(cross, strainAB, multilocus, homologous) %>%
  summarise(
    count_homologous = n(),
  )  %>%  pivot_wider(names_from = homologous, values_from = count_homologous, values_fill = 0) %>% unite("pairwise_code", c("Homologous", "Heterologous"), sep = "_")

pairwise_homologous_by_locus_count  <- pairwise_homologous_by_locus %>% group_by(cross, strainAB, pairwise_code) %>%
  summarise(
    number_pairwise_loci = length(multilocus),
  )

pairwise_homologous_by_locus_count$pairwise_code <- gsub("1_0", "One Homologous", pairwise_homologous_by_locus_count$pairwise_code)
pairwise_homologous_by_locus_count$pairwise_code <- gsub("1_1", "One Homologous, One Heterologous", pairwise_homologous_by_locus_count$pairwise_code)
pairwise_homologous_by_locus_count$pairwise_code <- gsub("1_2", "One Homologous, Two Heterologous", pairwise_homologous_by_locus_count$pairwise_code)
pairwise_homologous_by_locus_count$pairwise_code <- gsub("2_0", "Two Homologous", pairwise_homologous_by_locus_count$pairwise_code)
pairwise_homologous_by_locus_count$pairwise_code <- gsub("2_1", "Two Homologous, One Heterologous", pairwise_homologous_by_locus_count$pairwise_code)
pairwise_homologous_by_locus_count$pairwise_code <- gsub("2_2", "Two Homologous, Two Heterologous", pairwise_homologous_by_locus_count$pairwise_code)
pairwise_homologous_by_locus_count$pairwise_code <- gsub("0_1", "One Heterologous", pairwise_homologous_by_locus_count$pairwise_code)
pairwise_homologous_by_locus_count$pairwise_code <- gsub("0_1", "One Heterologous", pairwise_homologous_by_locus_count$pairwise_code)

#Figure 7C
ggplot(pairwise_homologous_by_locus_count, aes(x = reorder(pairwise_code, number_pairwise_loci), y = number_pairwise_loci)) + ylab("Number of Pairwise Loci Combinations") + xlab("Genotypes represented at least one in plaque isolates") + geom_col() + coord_flip() 

#Calculate percent of loci that had at least one heterologous pairing  
aggregate(pairwise_homologous_by_locus_count$number_pairwise_loci, by=list(Category=pairwise_homologous_by_locus_count$pairwise_code), FUN=sum)
#Category   x
#1                   One Homologous  63
#2 One Homologous, One Heterologous  56
#3 One Homologous, Two Heterologous  22
#4                   Two Homologous   2
#5 Two Homologous, One Heterologous  29
#6 Two Homologous, Two Heterologous 108

#Calculate by hand
sum(aggregate(pairwise_homologous_by_locus_count$number_pairwise_loci, by=list(Category=pairwise_homologous_by_locus_count$pairwise_code), FUN=sum)$x)
#280

#Number with at least one heterologous
((280-63-2) / 280)*100
#76.78571

#Maybe above analysis is not the most important or best. Let's go more basic and generate numbers

#Overall
aggregate(count ~ homologous, data=genotypes_table_two_locus_long, sum, na.rm=TRUE)
#  homologous   count
#1 Heterologous 2602
#2 Homologous   9074
#This is meaningless because different numbers of plaques isolated per cross. So.. 

#Lets look at proportions on average for homologous 
overall_homologous_proportions <- genotypes_table_two_locus_long %>% group_by(homologous) %>%
  summarise(
    freq = freq,
  )
aggregate(freq ~ homologous, data=overall_homologous_proportions, mean, na.rm=TRUE)
#  homologous   freq
#1 Heterologous 0.1596215
#2 Homologous   0.5368272
#On average, pairwise loci combinations were represented by 53.68% homolgous plaques and 15.96% heterologous plauqes

#Super Plot
#ggplot(genotypes_table_two_locus_long, aes(x = homologous, y = count, fill = homologous)) + geom_col() + theme(axis.text.x = element_blank()) + scale_fill_manual(values=c("orange", "darkblue")) + ylab("Frequency of progeny plaque isolates") + xlab("Strain combinations of pairwise segment combinations") + facet_grid(strainAB ~ multilocus, scales = "free_y")
#Need to export to large PDF to see

#The strain patterns seem much more consistent than the segment patterns, so going to hit avergae heterologous

subset(genotypes_table_two_locus_long, homologous == "Heterologous")

genotypes_table_two_locus_long %>% filter(homologous == "Heterologous") %>% group_by(cross, strainAB) %>%
  summarise(
    mean_frequency = mean(freq),
    max_frequency = max(freq),
    sd_frequency = sd(freq),
  )
# A tibble: 9 × 5
# Groups:   cross [9]
#cross   strainAB   mean_frequency max_frequency sd_frequency
#<chr>   <chr>               <dbl>         <dbl>        <dbl>
#1 cross1  HK68_TX12          0.0909        0.0909       0     
#2 cross11 CA09_PAN99         0.263         0.636        0.198 
#3 cross14 CA09_TX12          0.165         0.489        0.133 
#4 cross15 PAN99_SI86         0.122         0.345        0.0842
#5 cross18 CA09_SI86          0.275         0.844        0.222 
#6 cross19 SI86_TX12          0.204         0.375        0.112 
#7 cross3  HK68_PAN99         0.138         0.333        0.0897
#8 cross7  PAN99_TX12         0.0258        0.0610       0.0175
#9 cross8  HK68_SI86          0.158         0.348        0.0912

multilocus_average_heterologous <- genotypes_table_two_locus_long %>% filter(homologous == "Heterologous") %>% group_by(cross, strainAB, multilocus) %>%
  summarise(
    mean_frequency = mean(freq),
  ) 

#Add subtype information to data frame
multilocus_average_heterologous <- right_join(multilocus_average_heterologous, cross_stats[, c(1, 14)])
#Remove NA
multilocus_average_heterologous <- na.omit(multilocus_average_heterologous)

#New Figure 7A,  
ggplot(multilocus_average_heterologous, aes(x = reorder(strainAB, -mean_frequency), y = mean_frequency)) + geom_boxplot() + geom_point() + ylab("Proportion of plaques with heterologous pairwise segments") + xlab("Strains in experimental coinfection") + theme(legend.position = "none") + theme_tufte() + theme(text = element_text(size = 20,  family="Helvetica")) + theme(axis.text.x = element_text(size = 12)) + theme(legend.position = "none") + theme(panel.grid.major.y = element_line(color = "lightgray",size = 0.5)) + facet_wrap(~ same_subtype, scales = "free_x")

model9 <- lm(mean_frequency ~ strainAB, multilocus_average_heterologous)
summary(model9)
anova(model9)
#Analysis of Variance Table

#Response: mean_frequency
#Df Sum Sq  Mean Sq F value    Pr(>F)    
#strainAB    8 1.5378 0.192230  14.698 < 2.2e-16 ***
#Residuals 206 2.6943 0.013079   
#So pairwise associations between segments are different among experimental coinfections

#Is there a correlation between proportion heterologous versus propostion reassortants 
#Need to bind to cross_stats dataframe
multilocus_average_heterologous[c(1,2,4)]

cross_average_heterologous <- genotypes_table_two_locus_long %>% filter(homologous == "Heterologous") %>% group_by(cross, strainAB) %>%
  summarise(
    pairwise_heterologous_mean_frequency = mean(freq),
  )

#Join with cross stats
cross_stats <- left_join(cross_stats, cross_average_heterologous, by = c("cross", "strainAB"))

#Have an NA, which is a zero, because CA_HK did not have reassortants (it did but those were not 8 segment genotypes)
cross_stats$pairwise_heterologous_mean_frequency[is.na(cross_stats$pairwise_heterologous_mean_frequency)] <- 0

#This is Figure 7C
ggplot(cross_stats, aes(x = prop_reassortant, y = pairwise_heterologous_mean_frequency)) + geom_point(aes(colour = strainAB), size = 6) + geom_smooth(method = "lm") + xlab("Proportion of Reassortant Plaque Isolates") + ylab("Proportion of plaques with heterologous pairwise segments") + theme(legend.position = "none") + theme_tufte() + theme(text = element_text(size = 20,  family="Helvetica")) + theme(axis.text.x = element_text(size = 12)) + theme(legend.position = "none") + theme(panel.grid.major.y = element_line(color = "lightgray",size = 0.5))  

reassortment_pairwise_model <- lm(pairwise_heterologous_mean_frequency ~ prop_reassortant, cross_stats)
summary(reassortment_pairwise_model)
#Call:
#lm(formula = pairwise_heterologous_mean_frequency ~ prop_reassortant, 
#   data = cross_stats)

#Residuals:
#  Min        1Q    Median        3Q       Max 
#-0.072181 -0.042234 -0.007869  0.029267  0.130467 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)  
#(Intercept)       0.03404    0.03996   0.852   0.4191  
#prop_reassortant  0.22510    0.07052   3.192   0.0128 *
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.06364 on 8 degrees of freedom
#Multiple R-squared:  0.5601,	Adjusted R-squared:  0.5052 
#F-statistic: 10.19 on 1 and 8 DF,  p-value: 0.01277


#### 3. Coinfection Supernatant Titers ####
#AFTER SENDING DRAFT - This analysis is in Excel right now, need to move here

#Import file: NEED TO PUT RAW COUNTS IN HERE! USE PHAGE REPO AS TEMPLATE FOR CALCULATION
supernatant_titers <- read.csv("~/Dropbox/mixtup/Documentos/ucdavis/papers/influenza_GbBSeq_human/influenza_GbBSeq/data/supernatant_titers.csv")

#Plot raw titers
ggplot(supernatant_titers, aes(x = reorder(cross, titer), y = titer)) + geom_col() 

#Remove Chile83 strains
supernatant_titers <- filter(supernatant_titers, strainA != "Chile83") %>% filter(strainB != "Chile83")

#Plot raw titers again
ggplot(supernatant_titers, aes(x = reorder(cross, titer), y = titer)) + geom_col() + scale_y_continuous(labels = function(x) format(x, scientific = TRUE))

ggplot(supernatant_titers, aes(x = reorder(cross, titer), y = titer)) + geom_col() + scale_y_continuous(labels = scientific, breaks=seq(0, 10, 1))

library(scales)
ggplot(supernatant_titers, aes(x = reorder(cross, titer), y = titer)) + geom_col() + scale_y_log10()

#Heatmap
ggplot(supernatant_titers, aes(x = strainA, y = strainB, fill= titer)) + geom_tile() + geom_text(aes(label = round(titer, 4))) + scale_fill_gradient(low = "white", high = "red") 

#### 4. Population Genetics ####
#Making a table for the GenAlEx format for population genetics analyses

genotypes_table <- group_by(controls_df, cross_sample, strainAB, locus) %>%
  summarise(
    majority_strain = majority_strain,
  )

#Looking good!
genotypes_table <- spread(genotypes_table, locus, majority_strain)

#Merge all genotypes into one and make a code for each genotype
genotypes_table <- unite(genotypes_table, "genotype",  c("HA", "Mg", "NA", "NPd", "NS1d", "PAc", "PB1c", "PB2f"), remove = F)

genotypes_table <- separate(genotypes_table, cross_sample, into = c("cross", "sample"), sep = "_", remove = F)

#Remove NA's
genotypes_table <- na.omit(genotypes_table)

######################## These Seeem Redundant now ########################

#Graph by cross first
ggplot(genotypes_table, aes(x = strainAB, fill = genotype)) + geom_bar() + theme(legend.position = "none")

genotype_counts <- group_by(genotypes_table, cross, strainAB, genotype) %>%
  summarise(
    number_plaques = n()
  )

#Graph Counts per cross
ggplot(genotype_counts, aes(x = genotype, y = number_plaques, fill = genotype)) + geom_col() + theme(legend.position = "none", axis.text.x = element_blank()) + facet_wrap(~ strainAB, ncol = 1, scales = "free") 

#Now need to count the number of unique genotypes per cross
genotype_counts_cross <- group_by(genotype_counts, cross, strainAB) %>%
  summarise(
    total_genotypes = length(unique(genotype)),
    total_plaques = sum(number_plaques)
  )

genotype_counts_cross <- mutate(genotype_counts_cross, proportion_unique_genotypes = total_genotypes/total_plaques)

genotypes_df <- right_join(genotype_counts, genotype_counts_cross, by = c("cross", "strainAB"))

#Same graph as above Figure 8A
ggplot(genotypes_df, aes(x = reorder(genotype, number_plaques), y = number_plaques, fill = strainAB)) + geom_col() + theme_tufte() + theme(text = element_text(size = 20,  family="Helvetica")) + theme(legend.position = "none", axis.text.x = element_blank()) + theme(strip.text = element_text(face = "bold")) + xlab("8 Segment Genotype") + ylab("Number of Progeny Plaque Isolates") + facet_wrap(~ strainAB, ncol = 1, scales = "free") 

#Now genotypes per cross
ggplot(genotype_counts_cross, aes(x = reorder(strainAB, total_genotypes), y = total_genotypes)) + geom_col()

#Old reasortment graph from above
ggplot(cross_stats, aes(x = reorder(cross, prop_reassortant), y = prop_reassortant)) + geom_col() + geom_hline(yintercept = 0.40, linetype = 2, color = "grey", alpha = 0.75) + geom_hline(yintercept = 0.9921875, linetype = 2, color = "red", alpha = 0.75)

#Now let's merge these two into one data frame to compare
cross_stats <- right_join(cross_stats, genotype_counts_cross, by = c("cross", "strainAB"))

#Side by side
ggplot(cross_stats, aes(x = reorder(strainAB, prop_reassortant), y = prop_reassortant)) + geom_col() + geom_hline(yintercept = 0.40, linetype = 2, color = "red", alpha = 0.75)
#Genotypes per cross
#Leaving out for now....
ggplot(genotype_counts_cross, aes(x = reorder(strainAB, total_genotypes), y = total_genotypes)) +  geom_col()

#Looks like some of them actually are in different order. Makes sense, reassortant count doesn't count unique reassorants 
#Is there a correlation between prop of reassortants and unique genotypes?
  #May be good to point out strain combinations that produce the most diversity:
    #So there's: number of reassortants, productivity, and diversity generated
      #Different implications for say RNA-RNA interactions during packaging

#Clearly correlated
ggplot(cross_stats, aes(x = prop_reassortant, y = total_genotypes)) + geom_point(aes(colour = strainAB)) + geom_smooth(method = "lm") 
#But also need to control for the number of plaques used to calculate the number of genotypes

#Quick stats - With raw number of genotypes - this is not right becase didn't use same number of plaques to 
# calculate genotypes, which is also different than calculation of reassortment rate
reassortment_genotypes_model <- lm(total_genotypes ~ reassortants, cross_stats)
summary(reassortment_genotypes_model)
#summary(reassortment_genotypes_model)

#Call:
#  summary(reassortment_genotypes_model)

#Call:
#  lm(formula = total_genotypes ~ reassortants, data = cross_stats)

#Residuals:
#  Min       1Q   Median       3Q      Max 
#-15.8836  -6.9707  -0.9577   5.8789  14.7998 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)  
#(Intercept)   -0.9592     6.5370  -0.147   0.8870  
#reassortants   0.3417     0.1202   2.843   0.0217 *
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 10.41 on 8 degrees of freedom
#Multiple R-squared:  0.5027,	Adjusted R-squared:  0.4405 
#F-statistic: 8.085 on 1 and 8 DF,  p-value: 0.0217

#Now using the correct measure, the proportion of unique genotypes per plaque
#Figure 8B
ggplot(cross_stats, aes(x = prop_reassortant, y = proportion_unique_genotypes)) + geom_point(aes(colour = strainAB), size = 6) + geom_smooth(method = "lm") + xlab("Proportion of Reassortant Plaque Isolates") + ylab("Proportion of Plaque Isolates wih Unique Genotypes") + theme_tufte() + theme(text = element_text(size = 20,  family="Helvetica")) + theme(axis.text.x = element_text(size = 12)) + theme(legend.position = "none") + theme(panel.grid.major.y = element_line(color = "lightgray",size = 0.5))  
#Much tighter correlation

reassortment_genotypes_model <- lm(proportion_unique_genotypes ~ reassortants, cross_stats)
summary(reassortment_genotypes_model)
#summary(reassortment_genotypes_model)

#Call:
#  lm(formula = proportion_unique_genotypes ~ reassortants, data = cross_stats)

#Residuals:
#  Min       1Q   Median       3Q      Max 
#-0.09521 -0.05199 -0.01805  0.04423  0.11729 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)   0.09533    0.04896   1.947 0.087383 .  
#reassortants  0.00576    0.00090   6.400 0.000209 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.07797 on 8 degrees of freedom
#Multiple R-squared:  0.8366,	Adjusted R-squared:  0.8162 
#F-statistic: 40.96 on 1 and 8 DF,  p-value: 0.0002092


#Can we come up with a risk assesement measure? Based on three criteria (reassortment, genotypes, productivity)

#Let's bring in supernatant titers into the cross_stats data frame
#supernatant_titers <- unite(supernatant_titers, strainAB, c("strainA", "strainB"))

#supernatant_titers_small <- group_by(supernatant_titers, strainAB) %>%
  summarise(
    average_titer = mean(titer)
  )

#CHECK!! #Need to fix the order of names in StrainAB
View(right_join(supernatant_titers_small, cross_stats))
#CHECK!! NOT ALL TITERS AVAILABLE
ggplot(cross_stats, aes(x = prop_reassortant, y = total_genotypes)) + geom_point(aes(colour = strainAB, size = average_titer)) + geom_smooth(method = "lm") 
#Statistically different from free reassortment. So parentals still over-represented statistically, but close!
#Maybe let's add threshold to the reassortment graph (see Line 264)

#Let's scale up
p_values_free_reasortment <- sapply(cross_stats$reassortants, function(x) binom.test(x, 96, 254/256)$p.value)
#None non-significant, so all frequencies statistically different from free reassortment

#Incorporate into cross stats table 
cross_stats <- data.frame(cross_stats, p_values_free_reasortment)

#Now let's try to crack rarefaction analysis

#Packages are a hassle so let's come up with a function:

#First sample loci. I will sample as a coin flip (bernoulli) eight times. Then loci will be 0/1


#This loop takes a fixed number of attempts and figures out how many unique genotypes you get
#Set a seed for reproducibility
set.seed(68)
n <- 1
simulation_df <- NULL

for (n in 1:1000) {
  uniques <- NULL
  simul_genotypes <- NULL
  tries <- NULL
  current_try <- NULL
  i <- 1
  while(i < 97) {
    simul_genotypes <- rbind(simul_genotypes, rbinom(8, 1,.5))
    nrow(simul_genotypes)
    uniques <- nrow(unique(simul_genotypes))
    #print(uniques)
    current_try <- data.frame(attempt = i, unique_genotypes = as.numeric(uniques), trial_number = n)
    tries <- rbind(tries, current_try)
    #trial_number <- rep(n, times = nrow(tries))
    #tries <- cbind(tries, trial_number)
    i <- i+1
  }
  n <- n + 1
  simulation_df <- rbind(simulation_df, tries)
}

#Figure supplement?
ggplot(simulation_df, aes(x = attempt, y = unique_genotypes, colour = trial_number)) + geom_point(aes(alpha=0.3)) + geom_point(aes(x=32, y=20), colour="red", size = 3.5) + theme(legend.position = "none")

simulation_df %>% filter(attempt == 96) %>% group_by(attempt) %>% 
  summarise(
    min_unique = min(unique_genotypes),
    max_unique = max(unique_genotypes),
    mean_unique = mean(unique_genotypes),
    sd_unique = sd(unique_genotypes),
    n_trials = n(),
  )
# A tibble: 1 x 6
#attempt min_unique max_unique mean_unique sd_unique n_trials
#<dbl>      <dbl>      <dbl>       <dbl>     <dbl>    <int>
# 96         72         88        80.3      3.05     1000


#This loop takes a fixed number of unique genotypes and figures out how many attempts you need
set.seed(12)
n <- 1
simulation_df <- NULL

for (n in 1:100) {
  uniques <- 1
  simul_genotypes <- NULL
  tries <- NULL
  current_try <- NULL
  i <- 1
  while(uniques < 256) {
    simul_genotypes <- rbind(simul_genotypes, rbinom(8, 1,.5))
    nrow(simul_genotypes)
    uniques <- nrow(unique(simul_genotypes))
    #print(uniques)
    current_try <- data.frame(attempt = i, unique_genotypes = as.numeric(uniques), trial_number = n)
    tries <- rbind(tries, current_try)
    #trial_number <- rep(n, times = nrow(tries))
    #tries <- cbind(tries, trial_number)
    i <- i+1
  }
  n <- n + 1
  simulation_df <- rbind(simulation_df, tries)
}

ggplot(simulation_df, aes(x = attempt, y = unique_genotypes, colour = trial_number)) + geom_point(aes(alpha=0.3)) + theme(legend.position = "none")

#Add some lines for context and point of highest reassortment cross
ggplot(simulation_df, aes(x = attempt, y = unique_genotypes, colour = trial_number)) + geom_point(aes(alpha=0.3)) + ylim(0, 275) + scale_y_continuous(breaks = seq(0, 275, by = 25)) + theme(legend.position = "none") + geom_hline(yintercept = 256, linetype = 2, color = "red", alpha = 0.75) + geom_vline(xintercept = c(24, 96, 384), linetype = 1, color = "black") + geom_point(aes(x=32, y=20), colour="red", size = 3.5) + xlab("Simulated Number of Plaque Isolates with Complete Genotypes") + ylab("Number of Total Unique Genotypes")

#Now calculate the number of genotypes you would see for each sample number
simulation_df %>% filter(attempt == 24 | attempt == 96 | attempt == 384) %>% group_by(attempt) %>% 
  summarise(
    min_unique = min(unique_genotypes),
    max_unique = max(unique_genotypes),
    mean_unique = mean(unique_genotypes),
    sd_unique = sd(unique_genotypes),
    n_trials = n(),
  )
#attempt min_unique max_unique mean_unique sd_unique n_trials
#<dbl>      <dbl>      <dbl>       <dbl>     <dbl>    <int>
#24         20         24        23.0      1.03      100
#96         71         88        79.9      3.29      100
#384        185        212       198.       4.58      100


#Now calculate number of progeny you would have to sample 
simulation_df %>% filter(unique_genotypes == 256) %>% group_by(unique_genotypes) %>% 
  summarise(
    min_attempt = min(attempt),
    max_attempt = max(attempt),
    mean_attempt = mean(attempt),
    sd_attempt = sd(attempt),
    n_trials = n(),
  )
## A tibble: 1 x 6
#unique_genotypes min_attempt max_attempt mean_attempt sd_attempt n_trials
#<dbl>            <dbl>       <dbl>        <dbl>      <dbl>    <int>
#256               942        2716         1562       311.      1000

write.csv(genotypes_table, file = "/Users/mixtup/Dropbox/mixtup/Documentos/ucdavis/papers/influenza_GbBSeq_human/influenza_GbBSeq/outputs/genotypes_table.csv")

#Some manual editing in Excel to make GenAlEx format file
#library("poppr")
pairwise_reassortment <- read.genalex("/Users/mixtup/Dropbox/mixtup/Documentos/ucdavis/papers/influenza_GbBSeq_human/influenza_GbBSeq/outputs/genotypes_table.csv")
pairwise_reassortment

#Never got the software to work, so just going to input here


#### 5. Testing Strain Specificity ####

#Import matrix with crosses and strains for an ANOVA
#cross_data_runA_matrix <- 
cross_data_runA_matrix <- read.csv("~/Dropbox/mixtup/Documentos/ucdavis/papers/influenza_GbBSeq_human/influenza_GbBSeq/data/cross_data_runA_matrix.csv")

cross_matrix <- as.matrix(cross_data_runA_matrix[2:6])

#Now add strain information to cross_stats
cross_stats <- right_join(cross_stats, unique(subset(controls_df, select = c(cross, strainA, strainB))))

table(cross_stats$strainA, cross_stats$prop_reassortant)

model5 = lm(prop_reassortant ~ cross_matrix, cross_stats)

#What about 
model15 <- lm(reassortants ~ cross_matrix, cross_stats)

#Make a point plot with reasosrtment rates for each strain 
cross_stats$prop_reassortant*cross_matrix
#Works, but they are not in the right order!!

cross_stats_ordered <- arrange(cross_stats, cross)
cross_stats_ordered$prop_reassortant

cross_data_runA_matrix_ordered <- arrange(cross_data_runA_matrix, cross)

cross_matrix_ordered <- as.matrix(cross_data_runA_matrix_ordered[2:6])

strain_individual_reassortment <- cross_stats_ordered$prop_reassortant*cross_matrix_ordered

ggplot(strain_individual_reassortment, aes(x = strainA, y = strainB, fill= titer)) + geom_tile() + geom_text(aes(label = round(titer, 4))) + scale_fill_gradient(low = "white", high = "red") 

model6 = lm(prop_reassortant ~ cross_matrix_ordered, cross_stats_ordered)
summary(model6)
#Call:
#  lm(formula = prop_reassortant ~ cross_matrix_ordered, data = cross_stats_ordered)

#Residuals:
#  1        2        3        4        5        6        7        8        9       10 
#-0.11285 -0.20660  0.24826  0.01910  0.07465 -0.02604 -0.11632  0.29688 -0.10938 -0.06771 

#Coefficients: (1 not defined because of singularities)
#Estimate Std. Error t value Pr(>|t|)
#(Intercept)                0.23785    0.23757   1.001    0.363
#cross_matrix_orderedCA09   0.25347    0.17958   1.411    0.217
#cross_matrix_orderedHK68   0.05208    0.17958   0.290    0.783
#cross_matrix_orderedPAN99  0.19444    0.17958   1.083    0.328
#cross_matrix_orderedSI86   0.36111    0.17958   2.011    0.101
#cross_matrix_orderedTX12        NA         NA      NA       NA

#Residual standard error: 0.2199 on 5 degrees of freedom
#Multiple R-squared:  0.5186,	Adjusted R-squared:  0.1335 
#F-statistic: 1.347 on 4 and 5 DF,  p-value: 0.3692

coef(model6)

#### Daniel Runcie Meeting Dec 8 2021

#Key for model, where X is matrix on board that is *not* a column name in the data frame
model1 = lm(prop_reassortant ~ X,cross_stats)
anova(model1)

### Comand line dump from Daniel
View(cross_stats)
View(cross_data_runA)
View(cross_data_runA)
xA = model.matrix(~0+strainA,cross_data_runA)
xA
xB = model.matrix(~0+strainB,cross_data_runA)
xB
colnames(xA)
colnames(xB)
X = xA
X
cross_stats$class = rep(A:D,2)
cross_stats$class = rep(c('A','B','C','D','E'),2)
cross_stats
View(cross_stats)
X = model.matrix(~0+class,cross_stats)
X
X2 = X[c(2:10,1),]
X2
X3 = X+X2
X3
anova(lm(prop_reassortant~X3,cross_stats)
)
colnames(cross_stats)[5]
colnames(cross_stats)[5] = 'strainA'
cross_stats$strainB = cross_stats$strainA[c(2:10,1)]
View(cross_stats)
xA = model.matrix(~0+strainA,cross_stats)
xA
xB = model.matrix(~0+strainB,cross_stats)
X = xA+xB
X
model = lm(prop_reassortant ~ X,cross_stats)
model.matrix(model)
summary(model)
model = lm(prop_reassortant ~ X-1,cross_stats)
summary(model)
model.matrix(model)
model2 = lm(prop_reassortant ~ X-1,cross_stats)
model1 = lm(prop_reassortant ~ X,cross_stats)
coef(model2)
coef(model1)
coef(model1) + coef(model1)[1]
coef(model2)*2
coef(model2)-coef(model2)[1]
model.matrix(model)
model.matrix(model2)
rowSums(model.matrix(model2))
fitted(model2)
fitted(model1)
coef(model2)
model.matrix(model2)
coef(model1)
coef(model2)
coef(model1)*2+coef(model1)[1]
coef(model2)*2
summary(model1)
anova(model1)
anova(model2)
anova(model1)
model1 = lm(prop_reassortant ~ X,cross_stats)
model1 = lm(prop_reassortant ~ X,cross_stats)
anova(model1)





uniques <- 1
simul_genotypes <- NULL
tries <- NULL
while(uniques < 256) {
print(uniques)
tries <- rbind(tries, uniques)
simul_genotypes <- rbind(simul_genotypes, rbinom(8, 1,.5))
nrow(simul_genotypes)
uniques <- nrow(unique(simul_genotypes))
}
nrow(simul_genotypes)
plot(tries)

####### SLATED FOR DELETION #######
#What about 
model15 <- lm(reassortants ~ cross_matrix, cross_stats)

#Make a point plot with reasosrtment rates for each strain 
cross_stats$prop_reassortant*cross_matrix
#Works, but they are not in the right order!!

cross_stats_ordered <- arrange(cross_stats, cross)
cross_stats_ordered$prop_reassortant

cross_data_runA_matrix_ordered <- arrange(cross_data_runA_matrix, cross)

cross_matrix_ordered <- as.matrix(cross_data_runA_matrix_ordered[2:6])

strain_individual_reassortment <- cross_stats_ordered$prop_reassortant*cross_matrix_ordered

model6 = lm(prop_reassortant ~ cross_matrix_ordered, cross_stats_ordered)
summary(model6)
#Call:
#  lm(formula = prop_reassortant ~ cross_matrix_ordered, data = cross_stats_ordered)

#Residuals:
#  1        2        3        4        5        6        7        8        9       10 
#-0.11285 -0.20660  0.24826  0.01910  0.07465 -0.02604 -0.11632  0.29688 -0.10938 -0.06771 

#Coefficients: (1 not defined because of singularities)
#Estimate Std. Error t value Pr(>|t|)
#(Intercept)                0.23785    0.23757   1.001    0.363
#cross_matrix_orderedCA09   0.25347    0.17958   1.411    0.217
#cross_matrix_orderedHK68   0.05208    0.17958   0.290    0.783
#cross_matrix_orderedPAN99  0.19444    0.17958   1.083    0.328
#cross_matrix_orderedSI86   0.36111    0.17958   2.011    0.101
#cross_matrix_orderedTX12        NA         NA      NA       NA

#Residual standard error: 0.2199 on 5 degrees of freedom
#Multiple R-squared:  0.5186,	Adjusted R-squared:  0.1335 
#F-statistic: 1.347 on 4 and 5 DF,  p-value: 0.3692
coef(model6)
####### SLATED FOR DELETION #######


################################# CHECK APRIL 4 - BUT POTENTIALLY DELETE ###############################
#Now need some stats from these
locus_strain_cross <- group_by(locus_strain, cross, locus) %>%
  summarise(
    top_samples_cross = max(total_samples),
    total_samples_cross = sum(total_samples),
    proportion = max(total_samples)/sum(total_samples)
  )

#Let's try stats just on the raw number, but subset to loci with more than 50 total
summary.lm(aov(proportion ~ locus, locus_strain_cross))

#Get top proprtion
#bias_cross_locus <- group_by(locus_strain_cross, cross, locus) %>%
#  summarise(
#    top_proportion = max(proportion),
#  )

#Graph
ggplot(locus_strain_cross, aes(x = locus, y = proportion)) + geom_point(aes(colour = cross)) + geom_boxplot(alpha = 0.25) + geom_hline(yintercept = 0.50)

#Stats add 
bias_locus <- group_by(locus_strain_cross, locus) %>%
  summarise(
    mean_proportion = mean(proportion),
    sd_proportion = sd(proportion)
  )

#Add to big DF
locus_strain_cross <- right_join(bias_locus, locus_strain_cross)

#Proportions
ggplot(locus_strain_cross, aes(x = reorder(locus, mean_proportion), y = proportion)) + geom_point()

#Means only
ggplot(locus_strain_cross, aes(x = reorder(locus, mean_proportion), y = mean_proportion)) + geom_point()

#
summary.lm(aov(proportion ~ locus, locus_strain_cross))
#Call:
#  aov(formula = proportion ~ locus, data = locus_strain_cross)

#Residuals:
#  Min       1Q   Median       3Q      Max 
#-0.34518 -0.12016  0.00415  0.13074  0.22700 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.81377    0.04137  19.672   <2e-16 ***
#locusMg      0.03142    0.05850   0.537    0.592    
#locusNA     -0.02956    0.05850  -0.505    0.614    
#locusNPd     0.04379    0.05850   0.748    0.456    
#locusNS1d   -0.01578    0.05850  -0.270    0.788    
#locusPAc    -0.04077    0.05850  -0.697    0.487    
#locusPB1c    0.03537    0.05850   0.605    0.547    
#locusPB2f    0.04107    0.05962   0.689    0.492    
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.1548 on 103 degrees of freedom
#Multiple R-squared:  0.04348,	Adjusted R-squared:  -0.02153 
#F-statistic: 0.6689 on 7 and 103 DF,  p-value: 0.6979

#Proportion against cross
summary.lm(aov(proportion ~ cross, locus_strain_cross))
#Residual standard error: 0.1285 on 97 degrees of freedom
#Multiple R-squared:  0.3791,	Adjusted R-squared:  0.2959 
#F-statistic: 4.556 on 13 and 97 DF,  p-value: 4.938e-06

#Graph to show this:
ggplot(locus_strain_cross, aes(x = cross, y = proportion)) + geom_point(aes(colour = locus)) + geom_boxplot(alpha = 0.25)
#Reorder

################################# CHECK APRIL 4 - BUT POTENTIALLY DELETE ###############################
