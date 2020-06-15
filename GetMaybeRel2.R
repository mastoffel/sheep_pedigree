# use updated function as replacement in CRAN version of sequoia.

GetMaybeRel2 <- function (GenoM = NULL, SeqList = NULL, Pedigree = NULL, LifeHistData = NULL, 
                          ParSib = "par", Complex = "full", Err = 1e-04, ErrFlavour = "version2.0", 
                          MaxMismatch = NA, Tassign = 0.5, MaxPairs = 7 * nrow(GenoM), 
                          DumPrefix = c("F0", "M0"), quiet = FALSE) 
{
        #on.exit(.Fortran("deallocall"), add = TRUE)
        Excl <- CheckGeno(GenoM, quiet = quiet)
        if (ParSib == "full") 
                ParSib <- "sib"
        if (!ParSib %in% c("par", "sib")) 
                stop("Invalid value for 'ParSib', choose 'par' or 'sib'")
        if (!is.null(SeqList) && (!is.list(SeqList) | is.data.frame(SeqList))) {
                stop("SeqList must be a list or NULL")
        }
        if (!is.na(MaxMismatch)) {
                warning("NOTE: 'MaxMismatch' is deprecated & ignored; now calculated internally by 'CalcMaxMismatch'", 
                        immediate. = TRUE)
        }
        if ("Pedigree" %in% names(SeqList)) {
                Pedigree <- SeqList$Pedigree
                if (!all(rownames(GenoM) %in% Pedigree$id)) 
                        stop("SeqList Pedigree does not match GenoM")
        }
        else if ("PedigreePar" %in% names(SeqList)) {
                Pedigree <- SeqList$PedigreePar
                if (!all(Pedigree$id == rownames(GenoM))) 
                        stop("SeqList PedigreePar does not match GenoM")
        }
        else if (!is.null(Pedigree)) {
                if (!all(rownames(GenoM) %in% Pedigree[, 1])) 
                        stop("Pedigree must include all individuals in GenoM, or be NULL")
        }
        PrSb <- switch(ParSib, dup = 0, par = 1, sib = 2)
        if (is.null(MaxPairs) || MaxPairs < 0 || !sequoia:::is.wholenumber(MaxPairs)) {
                stop("'MaxPairs' but be a positive whole number")
        }
        if ("LifeHist" %in% names(SeqList)) {
                LhIN <- SeqList$LifeHist
        }
        else if (!is.null(LifeHistData)) {
                LhIN <- sequoia:::CheckLH(LifeHistData)
        }
        else {
                LhIN <- sequoia:::orderLH(LH = NULL, gID = rownames(GenoM))
        }
        if (any(LhIN$Sex == 4)) {
                if (!quiet) 
                        message("detected hermaphrodites (sex=4), changing Complex to 'herm'")
                GenoM <- herm_clone_Geno(GenoM, LhIN, herm.suf = c("f", 
                                                                   "m"))
                if (!"LifeHist" %in% names(SeqList)) {
                        LhIN <- herm_clone_LH(LhIN, herm.suf = c("f", "m"))
                }
                if (!is.null(Pedigree)) {
                        Pedigree <- herm_clone_Ped(Ped = Pedigree, LH = LhIN[, 
                                                                             1:3], herm.suf = c("f", "m"))
                        Pedigree <- Pedigree[match(rownames(GenoM), Pedigree[, 
                                                                             1]), ]
                }
        }
        gID <- rownames(GenoM)
        GenoV <- as.integer(GenoM)
        if (!is.null(Pedigree)) {
                Ped <- PedPolish(Pedigree, GenoNames = rownames(GenoM), 
                                 NAToZero = TRUE, DropNonSNPd = FALSE)
                DPnc <- nchar(DumPrefix)
                PedNum <- sequoia:::IDToNum(Ped[, 1:3], gID, DumPrefix)
                PedPar <- as.matrix(PedNum[match(gID, Ped$id), 2:3])
                PedPar[is.na(PedPar)] <- 0
        }
        else {
                PedPar <- rep(0, 2 * nrow(GenoM))
        }
        if ("Specs" %in% names(SeqList)) {
                Specs <- SeqList$Specs
                Ng <- sequoia:::FacToNum(Specs[, "NumberIndivGenotyped"])
                SMax <- sequoia:::FacToNum(Specs[, "MaxSibshipSize"])
                ErrM <- sequoia:::ErrToM(sequoia:::FacToNum(SeqList$Specs["GenotypingErrorRate"]), 
                               flavour = ErrFlavour, Return = "matrix")
                if (!"MaxMismatchOH" %in% names(SeqList$Specs)) {
                        sts <- sequoia:::SnpStats(GenoM, Plot = FALSE)
                        MaxMismatchV <- CalcMaxMismatch(Err = ErrM, MAF = sts[, 
                                                                              "AF"], ErrFlavour = ErrFlavour, qntl = 0.999^(1/Ng))
                }
                else {
                        MaxMismatchV <- setNames(sequoia:::FacToNum(SeqList$Specs[c("MaxMismatchDUP", 
                                                                          "MaxMismatchOH", "MaxMismatchME")]), c("DUP", 
                                                                                                                 "OH", "ME"))
                }
                Cmplx <- switch(Specs[, "Complexity"], full = 2, simp = 1, 
                                mono = 0, herm = 4)
                AP <- SeqList$AgePriors[, c("M", "P", "FS", "MS", "PS")]
                SpecsInt <- c(nSnp = sequoia:::FacToNum(Specs[, "NumberSnps"]), 
                              MaxMisDUP = sequoia:::FacToNum(Specs[, "MaxMismatchDUP"]), 
                              MaxMisOH = sequoia:::FacToNum(Specs[, "MaxMismatchOH"]), MaxMisME = sequoia:::FacToNum(Specs[, 
                                                                                                       "MaxMismatchME"]), SMax = as.integer(sequoia:::FacToNum(Specs[, 
                                                                                                                                                           "MaxSibshipSize"])), Complx = as.integer(Cmplx), 
                              quiet = as.integer(quiet), nAgeCl = sequoia:::FacToNum(Specs[, 
                                                                                 "nAgeClasses"]))
                SpecsDbl <- c(TF = sequoia:::FacToNum(Specs[, "Tfilter"]), TA = sequoia:::FacToNum(Specs[, 
                                                                                     "Tassign"]))
        }
        else {
                Ng <- nrow(GenoM)
                if (is.null(Pedigree)) {
                        SMax <- 100
                }
                else {
                        SMax <- max(table(Pedigree$dam), table(Pedigree$sire)) + 
                                1
                }
                Complx <- switch(Complex, full = 2, simp = 1, mono = 0, 
                                 herm = 4)
                AP <- MakeAgePrior(Pedigree, LifeHistData, Plot = FALSE, 
                                   quiet = TRUE)
                ErrM <- ErrToM(Err, flavour = ErrFlavour, Return = "matrix")
                sts <- SnpStats(GenoM, Plot = FALSE)
                MaxMismatchV <- CalcMaxMismatch(Err = ErrM, MAF = sts[, 
                                                                      "AF"], ErrFlavour = ErrFlavour, qntl = 0.999^(1/Ng))
                SpecsInt <- c(nSnp = as.integer(ncol(GenoM)), MaxMisDUP = as.integer(MaxMismatchV["DUP"]), 
                              MaxMisOH = as.integer(MaxMismatchV["OH"]), MaxMisME = as.integer(MaxMismatchV["ME"]), 
                              SMax = as.integer(SMax), Complx = as.integer(Complx), 
                              quiet = as.integer(quiet), nAgeCl = as.integer(nrow(AP)))
                SpecsDbl <- c(TF = as.double(-2), TA = as.double(Tassign))
        }
        SpecsIntAmb <- c(ParSib = as.integer(PrSb), nAmbMax = as.integer(MaxPairs))
        Nd <- 0
        DumParRF <- rep(0, 4 * as.integer(Ng/2))
        dID <- NULL
        if (!is.null(Pedigree)) {
                Nd <- c(sum(substr(Pedigree$id, 1, DPnc[1]) == DumPrefix[1]), 
                        sum(substr(Pedigree$id, 1, DPnc[2]) == DumPrefix[2]))
                if (max(Nd) > 0) {
                        SibshipGPs <- array(0, dim = c(2, max(Nd), 2), dimnames = list(c("grandma", 
                                                                                         "granddad"), 1:max(Nd), c("mat", "pat")))
                        for (k in 1:2) {
                                if (Nd[k] > 0) {
                                        SibshipGPs[, 1:Nd[k], k] <- t(as.matrix(PedNum[substr(Ped$id, 
                                                                                              1, DPnc[k]) == DumPrefix[k], 2:3]))
                                        for (s in 1:Nd[k]) {
                                                for (g in 1:2) {
                                                        x <- (k - 1) * 2 * as.integer(Ng/2) + (s - 
                                                                                                       1) * 2 + g
                                                        DumParRF[x] <- SibshipGPs[g, s, k]
                                                }
                                        }
                                }
                        }
                        dID <- c(Ped$id[substr(Ped$id, 1, DPnc[1]) == DumPrefix[1]], 
                                 Ped$id[substr(Ped$id, 1, DPnc[2]) == DumPrefix[2]])
                }
        }
        LHF <- sequoia:::orderLH(LhIN[LhIN$Sex %in% c(1:3), ], gID)
        TMP <- .Fortran(sequoia:::findambig, ng = as.integer(Ng), specsint = as.integer(SpecsInt), 
                        specsintamb = as.integer(SpecsIntAmb), specsdbl = as.double(SpecsDbl), 
                        errv = as.double(ErrM), genofr = as.integer(GenoV), sexrf = as.integer(LHF$Sex), 
                        byrf = as.integer(c(LHF$BirthYear, LHF$BY.min, LHF$BY.max)), 
                        aprf = as.double(AP), parentsrf = as.integer(PedPar), 
                        dumparrf = as.integer(DumParRF), namb = as.integer(0), 
                        ambigid = integer(2 * MaxPairs), ambigrel = integer(2 * 
                                                                                    MaxPairs), ambiglr = double(2 * MaxPairs), ambigoh = integer(MaxPairs), 
                        ntrio = as.integer(0), trioids = integer(3 * Ng), triolr = double(3 * 
                                                                                                  Ng), triooh = integer(3 * Ng))
        TMP$ambiglr[abs(TMP$ambiglr - 999) < 0.1] <- NA
        TMP$ambiglr <- round(TMP$ambiglr, 2)
        TMP$ambigoh[TMP$ambigoh < 0] <- NA
        TMP$triolr[abs(TMP$triolr - 999) < 0.1] <- NA
        TMP$triolr <- round(TMP$triolr, 2)
        TMP$triooh[TMP$triooh < 0] <- NA
        if (TMP$namb > 0) {
                RelName <- c("PO", "FS", "HS", "GP", "FA", "HA", "U ", 
                             "Q", "2nd")
                Na <- TMP$namb
                TMP$ambigid <- sequoia:::NumToID(TMP$ambigid, 0, gID, dID)
                AmbigRel <- factor(TMP$ambigrel, levels = 1:9, labels = RelName)
                MaybeRel <- data.frame(sequoia:::VtoM(TMP$ambigid, Na), sequoia:::VtoM(AmbigRel, 
                                                                   Na), sequoia:::VtoM(TMP$ambiglr, Na), stringsAsFactors = FALSE)
                names(MaybeRel) <- c("ID1", "ID2", "Relx", "TopRel", 
                                     "LLR_Rx_U", "LLR")
                MaybeRel <- MaybeRel[, -which(names(MaybeRel) %in% c("Relx", 
                                                                     "LLR_Rx_U"))]
                MaybeRel$OH <- TMP$ambigoh[sequoia:::s(Na)]
                if (!is.null(LhIN) & nrow(MaybeRel) > 0) {
                        LhIN$BirthYear[LhIN$BirthYear < 0] <- NA
                        MaybeRel <- merge(MaybeRel, setNames(LhIN[, 1:3], 
                                                             c("ID1", "Sex1", "BirthYear1")), all.x = TRUE)
                        MaybeRel <- merge(MaybeRel, setNames(LhIN[, 1:3], 
                                                             c("ID2", "Sex2", "BirthYear2")), all.x = TRUE)
                        MaybeRel$AgeDif <- with(MaybeRel, BirthYear1 - BirthYear2)
                        MaybeRel <- MaybeRel[, c("ID1", "ID2", "TopRel", 
                                                 "LLR", "OH", "BirthYear1", "BirthYear2", "AgeDif", 
                                                 "Sex1", "Sex2")]
                        for (i in 1:Na) {
                                if (is.na(MaybeRel$AgeDif[i])) 
                                        next
                                if (MaybeRel$AgeDif[i] < 0) {
                                        tmpRel <- MaybeRel[i, ]
                                        tmpRel$AgeDif <- abs(tmpRel$AgeDif)
                                        MaybeRel[i, ] <- tmpRel[, c("ID2", "ID1", "TopRel", 
                                                                    "LLR", "OH", "BirthYear2", "BirthYear1", 
                                                                    "AgeDif", "Sex2", "Sex1")]
                                }
                        }
                }
                if (grepl("sib", ParSib)) {
                        MaybeRel <- with(MaybeRel, MaybeRel[TopRel %in% c("PO", 
                                                                          "FS", "HS", "GP", "FA", "2nd", "Q"), ])
                }
                MaybeRel <- MaybeRel[order(ordered(MaybeRel$TopRel, levels = RelName), 
                                           -MaybeRel$LLR), ]
                if (nrow(MaybeRel) == 0) {
                        MaybeRel <- NULL
                }
                else {
                        rownames(MaybeRel) <- 1:nrow(MaybeRel)
                        MaybeRel$SNPdBoth <- sequoia:::CalcSnpdBoth(MaybeRel[, c("ID1", 
                                                                       "ID2")], GenoM)
                        if (any(LhIN$Sex == 4)) {
                                MaybeRel <- sequoia:::herm_unclone_MaybeRel(MaybeRel, Ped = NULL, 
                                                                  LH = LhIN, herm.suf = c("f", "m"))
                        }
                }
        }
        else MaybeRel <- NULL
        if (quiet < 1) {
                if (!is.null(MaybeRel)) {
                        nRel <- nrow(MaybeRel)
                        nPO <- sum(MaybeRel$TopRel == "PO")
                }
                else {
                        nRel <- 0
                        nPO <- 0
                }
                message("Found ", nPO, " likely parent-offspring pairs, and ", 
                        nRel - nPO, " other non-assigned pairs of possible relatives")
        }
        if (TMP$ntrio > 0) {
                trios <- data.frame(sequoia:::VtoM(TMP$trioids, nr = TMP$ntrio, 
                                         nc = 3), sequoia:::VtoM(TMP$triolr, nr = TMP$ntrio, nc = 3), 
                                    sequoia:::VtoM(TMP$triooh, nr = TMP$ntrio, nc = 3), stringsAsFactors = FALSE)
                names(trios) <- c("id", "parent1", "parent2", "LLRparent1", 
                                  "LLRparent2", "LLRpair", "OHparent1", "OHparent2", 
                                  "MEpair")
                for (k in 1:3) trios[, k] <- sequoia:::NumToID(trios[, k], k - 
                                                             1, gID, dID)
                trios$SNPd.id.parent1 <- sequoia:::CalcSnpdBoth(trios[, c("id", 
                                                                "parent1")], GenoM)
                trios$SNPd.id.parent2 <- sequoia:::CalcSnpdBoth(trios[, c("id", 
                                                                "parent2")], GenoM)
                if (any(LhIN$Sex == 4)) {
                        trios <- sequoia:::herm_unclone_Trios(trios, LH = LhIN, herm.suf = c("f", 
                                                                                   "m"))
                }
                if (quiet < 1) {
                        message("Found ", nrow(trios), " parent-parent-offspring trios")
                }
        }
        else trios <- NULL
        if (grepl("par", ParSib)) {
                return(list(MaybePar = MaybeRel, MaybeTrio = trios))
        }
        else {
                return(list(MaybeRel = MaybeRel, MaybeTrio = trios))
        }
}