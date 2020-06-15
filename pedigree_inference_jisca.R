# sheep pedigree reconstruction

library(sequoia)
library(plyr)
library(knitr)

# Genotype Matrix
GenoM <- GenoConvert("data/20190604_Plates_1_to_87_parentage.ped", InFormat = "ped")
#GenoM <- as.matrix(read.table("data/Soay_geno_20190604_Plates_1_to_87.txt", row.names=1, header=FALSE))

# life history
LH <- read.table("data/20190604_Plates_1_to_87_parentage.phenos", header=TRUE)

# extra birth years
ExtraBY <- read.csv("data/SoaySheep_BirthYear_guestimates.csv", header=T, stringsAsFactors=FALSE)
these <- match(ExtraBY$ID, LH$ID)
LH2 <- LH
LH2$BirthYear[these] <- ExtraBY$BirthYear

DBped <- read.table("data/20190604_Genetic_Field_Mothers.txt", 
                    header=T, stringsAsFactors=F)

IDstatus <- read.table("data/20190604_ID_Status.txt",
                       header=T, stringsAsFactors=F)

DBdata <- merge(DBped[,1:4], setNames(IDstatus, c("ID", "id.status")), all=TRUE)
DBdata <- merge(DBdata, setNames(IDstatus, c("ConsensusMumID", "mum.status")), all.x=TRUE)
DBdata <- merge(DBdata, setNames(IDstatus, c("DadID", "dad.status")), all.x=TRUE)
DBdata <- merge(DBdata, LH[,1:4], all=TRUE)

# sequoia output
load("data/Sequoia_out_2019-06-22.RData")
SeqOUT <- SeqSheep


# Check for Duplicates
SeqDup <- sequoia(GenoM = GenoM,
                  LifeHistData = LH[, c("ID", "Sex", "BirthYear")],
                  Err = 5e-4,          # guestimate of genotyping error rate
                  MaxMismatch = 5,     # max. number of mismatching SNPs
                  MaxSibIter = -1)     # no parentage assignment or sibship clustering yet

SeqDup$DupGenotype
# need to check whether these individuals are really matched across all SNPs (not only the 400)


# There are `r nrow(SeqDup$DupGenotype)` sample pairs that, based on the `r ncol(GenoM)` 
# SNPs used, are more likely to come from the same individual than from two close relatives. 
# Duplicates need to be removed, as they are not handled by the rest of the program. 
# Where the sex or birth year differs between the two putative IDs, it is set to missing for the non-removed ID.

# drop duplicates
GenoM2 <- GenoM[-c(SeqDup$DupGenotype$row2), ]

DUP <- merge(SeqDup$DupGenotype[, c("ID1", "ID2")], LH[,c("ID","Sex", "BirthYear")], 
             by.x="ID1", by.y="ID", all.x=TRUE)
DUP <- merge(DUP, LH[,c("ID","Sex", "BirthYear")], by.x="ID2", by.y="ID", all.x=TRUE,
             suffixes=c(".1", ".2"))
DUP[, c("ID1", "ID2", "Sex.1", "Sex.2", "BirthYear.1", "BirthYear.2")]

# drop duplicates2
LH2$Sex[LH2$ID %in% c(1741, 11554)] <- NA


# Parentage assignment
ParSheep <- sequoia(GenoM = GenoM2,
                    LifeHistData = LH2[, c("ID", "Sex", "BirthYear")],
                    Err = 5e-4,        # guestimate of genotyping error rate
                    MaxSibshipSize = 150,
                    MaxMismatch = 3,   # max. number of parent-offspring mismatches
                    MaxSibIter = 0)    # no sibship clustering yet
# assigned 6741 dams and 6322 sires to 7703 individuals
#save(ParSheep, file="data/Sequoia-out_parentage.RData")

SummarySeq(ParSheep)


# compare new to old pedigree
PC.par <- PedCompare(DBdata[, c("ID", "ConsensusMumID", "DadID")], ParSheep$PedigreePar)
# parent-offspring (GG = Genotyped parent, genotyped offspring)
# T = total number
# check mismatches
PC.par$Counts["GG",,]

### Mismatches in genotyped parents 
MisMatch <- merge(PC.par$Mismatch[,c("id", "dam.2", "sire.2", "Parent")], DBdata, by.x="id", by.y="ID", all.x=TRUE)
MisMatch[MisMatch$Parent=="dam", ]

MisMatch[MisMatch$Parent=="sire", ]

PC.par$MergedPed[which(PC.par$MergedPed$sire.1==-1382), ]

# find maybe relatives
MaybePar <- GetMaybeRel(GenoM = GenoM2, SeqList = ParSheep, ParSib = "par")

PC.par$P1only[PC.par$P1only$Cat == "GG", ]

# there are  10  likely parent-offspring pairs
save(MaybePar, file="data/MaybePar.RData")
MaybePar

SNPstats <- SnpStats(GenoM = GenoM2, Ped = ParSheep$PedigreePar)

ExtraBY <- read.csv("data/SoaySheep_BirthYear_guestimates.csv", header=T, stringsAsFactors=FALSE)
these <- match(ExtraBY$ID, LH$ID)
LH3 <- LH
LH3$BirthYear[these] <- ExtraBY$BirthYear

SeqSheep <- sequoia(GenoM = GenoM2,
                    LifeHistData = LH3[, c("ID", "Sex", "BirthYear")],
                    MaxSibshipSize = 250,
                    Err = 1e-4,        # guestimate of genotyping error rate
                    MaxMismatch = 3,   # max. number of parent-offspring mismatches
                    MaxSibIter = 20)
save(SeqSheep, file = "data/Sequoia_out_2019_full.RData")
load("data/Sequoia_out_2019-06-22.RData")

SummarySeq(SeqSheep)

PC <- PedCompare(DBdata[, c("ID", "ConsensusMumID", "DadID")], SeqSheep$Pedigree)
PC$Counts["TT",,]

MisMatch <- merge(PC$Mismatch[,c("id", "dam.2", "sire.2", "dam.r", "Parent")], DBdata, by.x="id", by.y="ID", all.x=TRUE)
with(MisMatch[MisMatch$Parent=="dam", ], kable(table(mum.status, substr(dam.2,1,2)=="F0")))
MisMatch[which(MisMatch$Parent=="dam" & MisMatch$mum.status!="Dummy"), ]

PC$MergedPed[which(PC$MergedPed$dam.2=="F0009"), ]

11396 %in% rownames(GenoM)

PC$MergedPed[which(PC$MergedPed$dam.1 %in% c(11396, -1378)), ]

NewMums <- PC$P2only[PC$P2only$Parent=="dam", ]
NewMums <- merge(NewMums[, c("id", "dam.2", "dam.r", "sire.2", "Cat")], 
                 DBdata[, c("ID", "id.status", "BirthYear", "Sex", "DeathYear")], by.x="id", by.y="ID", all.x=TRUE)

MaybeRel <- GetMaybeRel2(GenoM = GenoM2, SeqList = SeqSheep, ParSib = "sib") 
save(MaybeRel , file = "data/MaybeRel_2019_full.RData")
#load("data/MaybeRel_2019-06-23.RData")
MaybeRel$MaybeRel

MaybeDF <- MaybeRel$MaybeRel

ONLYP1 <- merge(PC$P1only[PC$P1only$Cat %in% c("GG","GD"), c("id", "dam.2", "sire.2", "Parent", "Cat")], 
                DBdata, by.x="id", by.y="ID", all.x=TRUE)
ONLYP1$id.in.MaybeDF <- ONLYP1$id %in% unlist(MaybeDF[,c("ID1","ID2")])


ONLYP1 <- ONLYP1[order(as.numeric(ONLYP1$id)), ]
