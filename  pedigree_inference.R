# pedigree reconstruction

library(sequoia)
library(plyr)
library(knitr)
library(tidyverse)
library(RJDBC)
library(dbplyr)

# remove problematic snps (From Susies checks)
snps_prob <- c("OAR2_129409414.1", "OAR23_39886336.1")

# get pedigree snps and remove problematic
ped_snps <- read_lines("data/2019/ParentageSNPs_431.txt")
ped_snps <- ped_snps[!(ped_snps %in% snps_prob)]
write_lines(ped_snps, "data/2019/ParentageSNPs_429.txt")

# filter plink pedigree snps
system(paste0("/usr/local/bin/plink --bfile data/2019/20200323_Plates_1to91 ",
              "--sheep --extract data/2019/ParentageSNPs_429.txt ",
              "--recode --out data/2019/20200323_Plates_1to91_parentage"))

# prepare data -----------------------------------------------------------------
# Genotype Matrix
GenoM <- GenoConvert("data/2019/20200323_Plates_1to91_parentage.ped", 
                     InFormat = "ped")

# get data from database
dbname <- "../sheep/data/db/StKilda_Data.accdb"
driver <- "net.ucanaccess.jdbc.UcanloadDriver"
driverpath <- "../sheep/data/db/UCanAccess/loader/ucanload.jar"
options <- paste0("jdbc:ucanaccess://", dbname, ";memory=false")

con <- DBI::dbConnect(JDBC(driver, driverpath), options)
tbls <- dbGetTables(con) %>% as_tibble()
sheep <- dbGetQuery(con, "Select * from Sheep") %>% as_tibble()
pregs <- dbGetQuery(con, "Select * from tblPregnancies") %>% as_tibble()
dbDisconnect(con)

# life history (ID, BirthYear, Sex, DeathYear, DeathMonth)
LH <- sheep %>% 
        left_join(pregs, by = "BirthRef") %>% 
        filter(Status != "Dummy") %>% 
        dplyr::select(ID, BirthYear, Sex, DeathYear, DeathMonth) %>% 
        mutate_all(as.integer)

# Pedigree until now (sys_PedigreeCombined query)
DBped <- read_delim("data/2019/sys_PedigreeCombined.txt", "\t")
IDstatus <- sheep %>% dplyr::select(ID, Status) %>% setNames(c("ID", "id.status"))

DBdata <- DBped %>% 
                dplyr::select(1:4) %>% 
                dplyr::full_join(IDstatus) %>% 
                dplyr::left_join(setNames(IDstatus, c("ConsensusMumID", "mum.status"))) %>% 
                dplyr::left_join(setNames(IDstatus,  c("DadID", "dad.status"))) %>% 
                dplyr::full_join(LH[,1:4]) %>% 
                dplyr::select(ID, DadID, ConsensusMumID, MumIDSource, id.status, 
                              mum.status, dad.status, BirthYear, Sex, DeathYear)

# List the maximum & minimum cohorts estimated for sheep of unknown age
# from sys_OPminimumAge
BYs <- read_delim("data/sys_OPminimumAge.txt", delim = ",", 
                  col_names = c("ID", "BY.min", "BY.max"))

# join minimum and maximum years with LH data
LH2 <- LH %>% 
        left_join(BYs) %>% 
        as.data.frame()
# check
LH2 %>% filter(ID %in% BYs$ID)
# oddities
# ID 8486 has birth year 2011, but 2014/2014 in cohort estimates, why?

# Duplicates -------------------------------------------------------------------
# Check for Duplicates
SeqDup <- sequoia(GenoM = GenoM,
                  LifeHistData = LH2[, c("ID", "Sex", "BirthYear")],
                  Err = 5e-4,          # guestimate of genotyping error rate
                  MaxSibIter = -1)     # no parentage assignment or sibship clustering yet
SeqDup$DupGenotype
write_delim(SeqDup$DupGenotype, "data/2019/output/duplicates.txt")
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

# drop duplicates2 \\ set duplicated IDs with different sex NA
# onle one duplicate had different sex
LH2$Sex[LH2$ID %in% c(1741)] <- NA

# Parentage assignment ---------------------------------------------------------
ParSheep <- sequoia(GenoM = GenoM2,
                    LifeHistData = LH2[,  c("ID", "Sex", "BirthYear", "BY.min", "BY.max")],
                    Err = 5e-4,        # guestimate of genotyping error rate
                    MaxSibshipSize = 150,
                    MaxSibIter = 0)    # no sibship clustering yet

# assigned 6974 dams and 6627 sires to 8014 individuals
save(ParSheep, file="data/2019/output/Sequoia-out_parentage.RData")
load("data/2019/output/Sequoia-out_parentage.RData")
SummarySeq(ParSheep)


## Comparison with database pedigree for genotyped parents

# The differences between the pedigree as currently in the database and the newly 
# inferred pedigree are split into three categories:
#         
# Mismatches: a parent assigned is assigned in each pedigree, but they differ
# Pedigree 1 (P1) only: a parent is assigned in the first (database) pedigree, but 
# not in the second (newly inferred) pedigree
# Pedigree 2 (P2) only: no parent is assigned in the first pedigree, but there 
# is one assiged in the second pedigree.
# 
# For now we are not yet looking at sibship clusters with dummy parents, and 
# therefore only considering those cases were both the focal individual and the 
# assigned parent (in the first and/or second pedigree) are genotyped 
# (category **GG**, for **G**enotyped individual **G**enotyped parent).


# compare new to old pedigree
PC.par <- PedCompare(Ped1 = DBdata[, c("ID", "ConsensusMumID", "DadID")], Ped2 = ParSheep$PedigreePar)
# parent-offspring (GG = Genotyped parent, genotyped offspring)
# T = total number
# check mismatches
PC.par$Counts["GG",,]
# PC.par$Counts["TT",,]
# parent
# class       dam sire
# Total    6971 6630
# Match    6938 6317
# Mismatch    4   51
# P1only     11    3
# P2only     18  259

### Mismatches in genotyped parents
MisMatchPar <- merge(PC.par$Mismatch[,c("id","dam.1", "sire.1", "dam.2", "sire.2", "dam.class", "sire.class")], 
                  DBdata, by.x="id", by.y="ID", all.x=TRUE)

# ask Jisca here whether things are alright

# mother mismatches #
MisMatchPar %>% filter(dam.class != "Match") %>% 
        kable(row.names=FALSE)
# duplicates
DUP[, c("ID1", "ID2", "Sex.1", "Sex.2", "BirthYear.1", "BirthYear.2")]

# These 6 mismatches include
# 2 mismatches where -1738 has been replaced by 8955 (which is now genotyped)
# 1 mismatch where -1789 has been replaced by NA
# 1 mismatch where NA is still NA (i.e. not really a mismatch)
# 1 mismatch where 10527 is now 11597. 10527 was previously assigned only as field mother
# 1 mismatch where 7717 is now 7605. These two were identified as duplicates.

# father mismatch # 
MisMatch %>% filter(sire.class != "Match") %>% 
        kable(row.names=FALSE)

# lots of stuff here, 52 mismatches
MisMatch %>% filter(sire.class != "Match") %>% nrow()

# lets dive in it:
# ID 8773 has been newly genotyped, accounting for 28 mismatches, i.e. 28 offspring!!!
MM1 <- MisMatch %>% filter(sire.2 == 8773)

# check that all his offspring are now assigned the correct dad
PC.par$MergedPed[which(PC.par$MergedPed$sire.1==-1957), ]

# which other IDs have now been genotyped? 
# -> 22 offspring from 9 sires
MM2 <- MisMatch %>% 
        filter((sire.class != "Match") & (sire.2 != 8773) & ((as.numeric(sire.1) < -1000) & (as.numeric(sire.2) > 0))) #%>% 
        #.$sire.2 %>% 
        #table() %>% sum()

# check that all of their offspring are now assigned correct dad
# looks good
map_df(o_sires, function(sire) PC.par$MergedPed[which(PC.par$MergedPed$sire.1==sire), c("id", "sire.1", "sire.2")])

# 2 mismatches left
MisMatch %>% 
        filter(sire.class != "Match") %>% 
        anti_join(rbind(MM1, MM2), by = c("id", "sire.1", "sire.2"))
# id -145 now assigned, previously 9245 (CO249) // duplicates
# id 11958, with previously NA sire has now found its father, with ID 9063


### Genotyped parents in database only
# Database parents may be non-assigned because they are not a genetic match, 
# which may either mean that the database is incorrect, or that either sample got 
# mislabeled in the field or the lab. However, parents may also be not-yet-assigned 
# when they are a good genetic match but the birthyear of one or both individuals 
# is unknown. Then, it is impossible to tell which is the two is the parent and 
# which is the offspring, unless a parent-parent-offspring trio can be found. 
# 
# Note that the lifehistory data used includes birth years guestimated during 
# previous rounds ('SoaySheep_BirthYear_guestimates.csv'), by iteratively 
# 
# * running parentage assignment
# * checking for non-assigned parent-offspring pairs
# * see if cohort estimate is compatible with all parent-offspring pairs an individual 
# is involved in (cohort not identical to known birth year of putative parent(s) or offspring)
# * see if any parent-offspring pairs match field mums, and pick a plausible birth 
# year (i.e., if an individual *is* a field mum, assign it a birth year some time before its eldest offspring)
# 
# See also Section 'FAQ - Is age information really needed?' of the sequoia user 
# guide ( https://cran.r-project.org/web/packages/sequoia/vignettes/sequoia.pdf )

MaybePar <- GetMaybeRel(GenoM = GenoM2, SeqList = ParSheep, ParSib = "par")
MaybePar$MaybePar

# Found 9 likely parent-offspring pairs, and 0 other non-assigned pairs of possible relatives
save(MaybePar, file="data/2019/output/MaybePar.RData")
#load("data/2019/output/MaybePar.RData")
MaybePar

# Some of those issues are resolved during full pedigree reconstruction (see further), 
# when parent - dummy parent - offspring trios help to orient the parent-offspring (PO) pairs correctly:
        
ONLYP1 <- PC.par$P1only %>% 
                filter((id.dam.cat == "GG")) %>% 
                left_join(rename(DBdata, id = ID) %>% mutate(DBdata, id = as.character(id))) %>% 
                arrange(as.numeric(id)) %>% 
                select(c("id", "dam.2", "ConsensusMumID", "MumIDSource", "id.status",
                 "BirthYear", "Sex", "DeathYear")) %>% 
                as_tibble()
print(ONLYP1, n = 100)


### Genotyped parents newly assigned only
#These are not (yet) checked, as the expectation is that (nearly) all of these will involve newly genotyped individuals. 

## Estimate genotyping error
# Sequoia isn't very sensitive to the value of the presumed genotyping error rate, but it is good to check we are at least in the right ballpark. 
SNPstats <- SnpStats(GenoM = GenoM2, Ped = ParSheep$PedigreePar)
table(SNPstats[,"Err.hat"])


## Sibship clustering
SeqSheep <- sequoia(GenoM = GenoM2,
                    LifeHistData =  LH2[,  c("ID", "Sex", "BirthYear", "BY.min", "BY.max")],
                    MaxSibshipSize = 250,
                    Err = 1e-4,        # guestimate of genotyping error rate
                    MaxSibIter = 20)
# save(SeqSheep, GenoM2, file="data/2019/Sequoia_out_2019-06-12.RData")
load(file="data/2019/Sequoia_out_2019-06-12.RData")

SummarySeq(SeqSheep)

## Comparison between database pedigree and full reconstructed pedigree
DBdata2 <- DBdata %>% mutate_at(vars(1:3), as.character) %>% as.data.frame()
PC <- PedCompare(DBdata2[, c("ID", "ConsensusMumID", "DadID")], SeqSheep$Pedigree)
PC$Counts["TT",,]
PC$P1only

# The number of mismatches is considerably higher now, largely because all used-to-be-dummy parents 
# that are now replaced by newly genotyped individuals are counted as mismatches. 
nrow(PC$Mismatch)
MisMatch <- PC$Mismatch %>% 
                select(id, dam.1, sire.1, dam.2, sire.2, dam.r, dam.class, sire.class) %>% 
                left_join(DBdata %>% rename(id = ID) %>% mutate(id = as.character(id)), by = "id")
MisMatch %>% 
        filter(sire.class != "Match")

MisMatch %>% filter(grepl("F0", dam.2))

# check out mismatches for mothers
MisMatch %>% filter(dam.class != "Match") %>% filter(mum.status != "Dummy")


# 
MaybeRel <- GetMaybeRel2(GenoM = GenoM2, SeqList = SeqSheep, ParSib = "sib") 
MaybeRel$MaybeRel


