# make the sequoia inferred pedigree ready to be uploaded to the database

# for inference & sleuthing, see Soay_sheep_pedigree_2018_cohort.Rmd

# A template can be found in the database: Menu > Genetics > Parentage data > Import new data

# TODO: compare R_PED with R_GRM

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


setwd("E:/Soay sheep/Sheep 2019")
load("Sequoia_out_2019-06-22.RData")   # SeqSheep
SNPd <- SeqSheep$PedigreePar$id
DBped <- read.table("20190604_Soay_Plates_1to87_for_Jisca/20190604_Genetic_Field_Mothers.txt",
                    header=T, stringsAsFactors=F)  # n=9284
IDstatus <- read.table("20190604_Soay_Plates_1to87_for_Jisca/20190604_ID_Status.txt",
                       header=T, stringsAsFactors=F)  # n=12035

library(sequoia)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# step 1: replace dummy mothers by non-genotyped field mothers
# (minimise number of dummy - non-genotyped individual duplicate individuals in database )


PC <- PedCompare(DBped[, c("ID", "ConsensusMumID", "DadID")], SeqSheep$Pedigree)

PedTMP <- PC$MergedPed[!is.na(PC$MergedPed$id),    # not non-genotyped individuals in DBped only
                       c("id", "dam.1", "dam.2", "dam.r")]
# dam.r: real, replacement dam
# NOTE: matching works based on OFFSPRING ONLY, and currently ignores sibship-grandparents
#  it is considered a match if the inferred sibship which contains the most offspring 
#  of a non-genotyped parent, consists for more than half of this individual's offspring.
# for more elaborate matching, see Rum deer pedigree inference tutorial. 

View(PedTMP[!is.na(PedTMP$dam.r), ])
# nomatch: no match found
# negative number, smaller than -1000: dummy in database, will be matched automatically 
#      when uploading to database

PedTMP$Mother <- with(PedTMP, ifelse(is.na(dam.2), NA,
                                     ifelse(substr(dam.2,1,2) != "F0", dam.2,
                                            ifelse(dam.r == "nomatch", dam.2,
                                                   ifelse(as.numeric(dam.r) < -1000, dam.2,
                                                          dam.r)))))
table(PedTMP$Mother == PedTMP$dam.r, useNA="ifany")
# FALSE  TRUE  <NA> 
#   336    49  7786   # this exercise hardly necessary in soay sheep, v. high sampling rate

chk.r <- unique(PedTMP[which(PedTMP$Mother == PedTMP$dam.r), c("dam.2", "dam.r")])
chk.r <- merge(chk.r, setNames(PedTMP[,1:3], c("dam.2", "MGM", "MGF")), all.x=TRUE)
chk.r <- merge(chk.r, DBped[, 1:3], by.x="dam.r", by.y="ID", all.x=TRUE)
# dam.r dam.2   MGM   MGF ConsensusMumID   MumIDSource
#  -143 F0045  <NA>  <NA>             NA          <NA>
#  1103 F0001 F0023  <NA>          -1045       Genetic  # 1103 MHS with 54,1600,2207. plausible.
# 11395 F0185 F0009 11085             NA          <NA>  # see sleuthing, F0009 = 11396 = 8955
#  1787 F0087  <NA>  1201             NA          <NA>
#  2847 F0015  2165  2282           2165 Field/Genetic  # OK
#  4106 F0008  3088    47           3088 Field/Genetic  # OK
#   416 F0044  <NA>  <NA>             NA          <NA>
#     5 F0062  <NA>  <NA>             NA          <NA>
#  8689 F0160 10721  6721          10721       Genetic  # plausible
#  8751 F0006  6929  5683           6929 Field/Genetic  # OK
#  8763 F0196  6435  6721           6435         Field  # OK
#  8774 F0048  7228  8997           7228 Field/Genetic  # OK
#  8785 F0012  7252  7742           7252 Field/Genetic  # OK


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Mother.cat & Father.cat

# For historical reasons:
# PB : parentage (GG in pedcompare) + behavioural (MumIDSource = Field/Genetic)
# P : parentage (GG), no observation-based confirmation ( = )
# SB : sibship + behavioural (i.e., dummy replaced)
# S : sibship (dummy parent to database)
# G : grandparent; for (replaced) dummies: their parent is assigned as their offspring's grandparent
# GS : grandparent + sibship; individual & parent are both (replaced) dummies
# x : none
# ! : mismatch

# source("Pedigree_Rum_utils.R")  # . See end of this file

PedTMP$Mother.Cat <- MakePCat(PedTMP[, c("id", "dam.2", "dam.1")], SNPd = SeqSheep$PedigreePar$id)
table(PedTMP$Mother.Cat, useNA="ifany")
# !B    B    G   GS    P  P!B  P?B   PB    S  S!B   SB    x 
#  3   44   52   34   28    8    3 6703   44    1  306  945 

PedTMP[which(PedTMP$Mother.Cat=="!B"), ]  # Behavioural mum not genetic match
#        id dam.1 dam.2 dam.r  dam Mother.cat
# 6451 1078    23  <NA>  <NA> <NA>         !B
# 6561 1639    21  <NA>  <NA> <NA>         !B
# 7443 1988    24  <NA>  <NA> <NA>         !B
# all known, unresolved issues.

PedTMP[which(PedTMP$Mother.Cat=="P!B"), ]   # mismatches, all resolved

PedTMP[which(PedTMP$Mother.Cat=="S!B"), ]
#      id dam.1 dam.2   dam.r   dam Mother.cat
# 76 9863  6866 F0049 nomatch F0049        S!B  # resolved


PedOUT <- merge(PedTMP[, c("id", "Mother", "Mother.Cat")], 
                SeqSheep$Pedigree[ ,-2])
PedOUT$Mother.Cat[is.na(PedOUT$Mother)] <- "x"

# dummy dams have been replaced when occuring as mum, but not yet as focal indiv
table(chk.r$dam.r %in% PedOUT$id)  # all FALSE, OK.
for (i in 1:nrow(chk.r)) {
        PedOUT$id[PedOUT$id == chk.r$dam.2[i]] <- chk.r$dam.r[i] 
}


PedOUT$Father.Cat <- MakePCat(cbind(PedOUT[, c("id", "sire")], sire.obs=NA),
                              SNPd = SeqSheep$PedigreePar$id)
table(PedOUT$Father.Cat, useNA="ifany")
#   G   GS    P    S    x 
# 266   25 6328  802  750 

names(PedOUT)[names(PedOUT)=="sire"] <- "Father"
names(PedOUT)[names(PedOUT)=="id"] <- "Code"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Confidence probabilities

# Confidence probabilities are estimated by repeatedly simulating genotype data according to the 
# empirical pedigree, running sequoia, and counting mismatches between the empirical pedigree and 
# based-on-simulated-data pedigree. 

# Done for just the parentage assignment...
# Doing so for full pedigree reconstruction is on to-do list, but is going to take few days
# of computing time. Maybe possible to run on Eddie?

# NOTE: sequoia's SimGeno assumes completely independent SNPs, maybe try with AlphaSim?


LH <- read.table("20190604_Soay_Plates_1to87_for_Jisca/20190604_Plates_1_to_87_parentage.phenos",
                 header=TRUE)
ExtraBY <- read.csv("SoaySheep_BirthYear_guestimates.csv", header=T, stringsAsFactors=FALSE)
these <- match(ExtraBY$ID, LH$ID)
LH2 <- LH
LH2$BirthYear[these] <- ExtraBY$BirthYear


ConfPr <- EstConf(Ped = SeqSheep$Pedigree, 
                  LifeHistData = LH2[, c("ID", "Sex", "BirthYear")],
                  args.sim= list(nSnp = 431, SnpError = 5e-4, 
                                 ParMis = c(0.15, 0.2)),
                  args.seq = list(Err = 1e-4, MaxMismatch = 3, 
                                  MaxSibshipSize=250, MaxSibIter=0),
                  nSim = 20)

# 0 assignment errors in 7 rounds
nrow(SeqSheep$Pedigree)   # 8171
1/(7*2*8171)   # 8.7E-6

# so, parentage confidence *if SNPs were completely independent* is >0.9999 

# For now, assume confidence probabilities are similar to those in the Rum deer:
ConfProb <- c(GG=0.999, GD=0.996, DG=0.998, DD=0.994)
# sham values for parent pairs:
ParProb3 <- c("GGG" = 0.9999,
              "GGD" = 0.998,
              "GDD" = 0.99,
              "DGG" = 0.999,
              "DGD" = 0.99,
              "DDD" = 0.98)


f.GProb <- function(id, par.G, SNPd, PG) {
        Gprob <- ifelse(is.na(par.G), NA,
                        ifelse(id %in% SNPd,
                               ifelse(par.G %in% SNPd, PG["GG"], PG["GD"]),
                               ifelse(par.G %in% SNPd, PG["DG"], PG["DD"])))
}


PedOUT$Mother.Prob <- with(PedOUT, sapply(seq_along(Code), function(i) 
        f.GProb(Code[i], Mother[i], SNPd=SNPd, PG=ConfProb)))
PedOUT$Father.Prob <- with(PedOUT, sapply(seq_along(Code), function(i) 
        f.GProb(Code[i], Father[i], SNPd=SNPd, PG=ConfProb)))
PedOUT$Trio.Prob <- with(PedOUT, ifelse(is.na(Mother) | is.na(Father), NA,
                                        ifelse(Code %in% SNPd,
                                               ifelse(Mother %in% SNPd & Father %in% SNPd, ParProb3["GGG"],
                                                      ifelse(Mother %in% SNPd | Father %in% SNPd, 
                                                             ParProb3["GGD"], ParProb3["GDD"])),
                                               ifelse(Mother %in% SNPd & Father %in% SNPd, ParProb3["DGG"],
                                                      ifelse(Mother %in% SNPd | Father %in% SNPd, 
                                                             ParProb3["DGD"], ParProb3["DDD"])))))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ExclPrevFather

# This column can be used to flag cases where a previous (microsatellite-based) assigned father
# is not a match according to the newly inferred pedigree, but no new father has been assigned.
# If not flagged, by default the database uses the last-known father in the consensus pedigree.

PedOUT$ExclPrevFather <- FALSE
PedOUT$ExclPrevFather[PedOUT$Code %in% c(1988, 2676)] <- TRUE   # 2x microsat father; see sleuthing.



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# write to file

names(PedOUT)[names(PedOUT)=="LLRdam"] <- "Mother.LLR"
names(PedOUT)[names(PedOUT)=="LLRsire"] <- "Father.LLR"
names(PedOUT)[names(PedOUT)=="LLRpair"] <- "Trio.LLR"
PedOUT <- PedOUT[, c("Code", "Mother", "Father", "Mother.Cat", "Father.Cat",
                     "Mother.LLR", "Father.LLR", "Trio.LLR",
                     "Mother.Prob", "Father.Prob", "Trio.Prob")]

# double check all IDs are valid
table(PedOUT$Code %in% IDstatus$ID, substr(PedOUT$Code,1,2) %in% c("F0", "M0"))
#       FALSE TRUE
# FALSE     0  455
# TRUE   7716    0

write.csv(PedOUT, "Pedigree_SoaySheep_2019-07-03_toDB.csv", row.names=F, quote=FALSE)

write.table(SeqSheep$Pedigree, "Pedigree_SoaySheep_2019-07-03_SNPonly.txt",
            row.names=F, quote=F)







#############################################################################
#############################################################################
# TODO: shorten & annotate
MakePCat <- function(DF,  # ID, GPar, Bpar
                     SNPd)  # vector with IDs of SNPd individuals 
{
        names(DF) <- c("ID", "Gpar", "Bpar")
        Pcat <- with(DF, ifelse(ID %in% SNPd,
                                ifelse(!is.na(Gpar),
                                       ifelse(Gpar %in% SNPd,      
                                              ifelse(!is.na(Bpar) & Bpar!="none",
                                                     ifelse(Gpar == Bpar,
                                                            "PB",
                                                            ifelse(Bpar %in% SNPd,
                                                                   "P!B",
                                                                   "P?B")),
                                                     "P"),
                                              ifelse(!is.na(Bpar) & Bpar!="none",
                                                     ifelse(Bpar %in% SNPd,
                                                            "S!B",   # should have been picked up
                                                            "SB"),
                                                     "S")),
                                       ifelse(!is.na(Bpar) & Bpar!="none",
                                              ifelse(Bpar %in% SNPd,
                                                     "!B",
                                                     "B"),
                                              "x")),
                                ifelse(ID %in% Bpar,
                                       ifelse(!is.na(Gpar),
                                              ifelse(Gpar %in% SNPd,
                                                     ifelse(!is.na(Bpar) & Bpar!="none",
                                                            ifelse(Gpar == Bpar,
                                                                   "GB",
                                                                   ifelse(Bpar %in% SNPd,
                                                                          "G!B",
                                                                          "G?B")),
                                                            "G"),
                                                     ifelse(!is.na(Bpar) & Bpar!="none",
                                                            ifelse(Bpar %in% SNPd,
                                                                   "GS!B",   # should have been picked up
                                                                   "GSB"),
                                                            "GS")),
                                              ifelse(!is.na(Bpar) & Bpar!="none",
                                                     ifelse(Bpar %in% SNPd,
                                                            "!B",
                                                            "B"),
                                                     "x")),
                                       ifelse(!is.na(Gpar),
                                              ifelse(Gpar %in% SNPd,
                                                     "G",
                                                     "GS"),
                                              "x"))))
        Pcat                                      
}