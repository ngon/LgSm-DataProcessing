traitcovs <- vector("list", length=24)

names(traitcovs) <- c("ppi3.logit", "ppi6.logit", "ppi12.logit", "habituation", "startle",
                      "cpp1.t", "cpp1.1",
                      "cpp8.t", "cpp8.1",
                      "cpp.diff","cpp.diff.p","cpp.diff1",
                      "act1.t", "act1.1",
                      "act2.t",
                      "act3.t",
                      "act4.t",
                      "act5.t",
                      "act8.t",
                      "sens",
                      "wild.binary", "tail","glucose","is.coatW")


traitcovs[["ppi3.logit"]] <- list("one", "sex", "is.ppi.box3", "ppi.weight", "is.batch4")
traitcovs[["ppi6.logit"]] <- list("one", "sex", "is.ppi.box3", "is.ppi.box4", "ppi.weight")
traitcovs[["ppi12.logit"]] <- list("one", "sex", "is.ppi.box3", "is.ppi.box4","ppi.weight",
                                   "is.batch3", "is.batch4", "is.batch7", "is.batch9")
traitcovs[["habituation"]] <- list("one", "sex", "ppi.weight")
traitcovs[["startle"]] <- list("one", "sex", "is.ppi.box3", "is.ppi.box4", "ppi.weight",
                               "is.batch2", "is.batch16", "is.batch17")
traitcovs[["cpp1.t"]]  <- list("one", "sex", "is.gen53", "is.gen54", "is.gen55","is.gen56", "is.batch11", "is.batch13", "is.batch14","is.batch17", "is.batch19","is.batch21", "is.batch22", "is.cpp.box3", "is.cpp.box8", "is.cpp.box9", "is.cpp.box10", "is.cpp.box11")
traitcovs[["cpp1.1"]]  <- list("one", "sex", "is.gen55","is.gen56", "is.batch2", "is.batch7", "is.batch14","is.batch17", "is.batch19","is.batch21", "is.batch22", "is.cpp.box5", "is.cpp.box8")
traitcovs[["cpp8.t"]]  <- list("one", "sex", "is.gen53", "is.gen56", "is.batch2", "is.batch10", "is.batch13","is.batch15", "is.batch21", "is.cpp.box2", "is.cpp.box3", "is.cpp.box4", "is.cpp.box6", "is.cpp.box7", "is.cpp.box8", "is.cpp.box9", "is.cpp.box10", "is.cpp.box11", "is.cpp.box12", "comerr8")
traitcovs[["cpp8.1"]]  <- list("one", "sex", "is.cpp.box2", "is.cpp.box3", "is.cpp.box4", "is.cpp.box8", "is.cpp.box10", "is.cpp.box11", "is.cpp.box12")
traitcovs[["cpp.diff"]]  <- list("one", "sex", "is.cpp.box2", "is.cpp.box3", "is.cpp.box4","is.cpp.box6", "is.cpp.box7", "is.cpp.box11", "is.cpp.box12", "is.batch2", "is.batch3","is.batch9", "is.batch14", "is.batch16", "is.batch17", "is.batch18", "is.batch19", "is.batch21", "is.batch22")
traitcovs[["cpp.diff.p"]]  <- list("one", "sex", "is.cpp.box2", "is.cpp.box3", "is.cpp.box4","is.cpp.box6", "is.cpp.box7", "is.cpp.box11", "is.cpp.box12", "is.batch2", "is.batch3","is.batch9", "is.batch14", "is.batch16", "is.batch17", "is.batch18", "is.batch19", "is.batch21", "is.batch22")
traitcovs[["cpp.diff1"]]  <- list("one", "sex","is.gen55", "is.gen56","is.cpp.box2", "is.cpp.box3", "is.cpp.box4","is.cpp.box6", "is.cpp.box7", "is.cpp.box11", "is.cpp.box12", "is.batch2", "is.batch3","is.batch9","is.batch14","is.batch16","is.batch17","is.batch18","is.batch19","is.batch21", "is.batch22")
traitcovs[["act1.t"]]  <- list("one", "sex", "is.batch14", "is.cpp.box3", "is.cpp.box4", "is.cpp.box5", "is.cpp.box6", "is.cpp.box7", "is.cpp.box8", "is.cpp.box9", "is.cpp.box11", "is.cpp.box12", "comerr1")
traitcovs[["act1.1"]]  <- list("one", "sex", "is.gen51", "is.batch14", "is.batch20", "is.cpp.box2", "is.cpp.box3", "is.cpp.box4", "is.cpp.box5", "is.cpp.box6", "is.cpp.box7", "is.cpp.box8", "is.cpp.box9", "is.cpp.box11", "is.cpp.box12", "comerr1")
traitcovs[["act2.t"]]  <- list("one", "sex", "is.gen51", "is.gen52","is.gen56", "is.cpp.box5", "is.cpp.box7", "is.cpp.box8", "is.cpp.box10", "is.cpp.box11","is.cpp.box12", "is.batch8","is.batch15",
traitcovs[["act3.t"]]  <- list("one", "sex", "is.gen51", "is.gen52", "is.gen53", "is.gen54", "is.gen55", "is.batch2", "is.batch3", "is.batch4", "is.batch5", "is.batch6", "is.batch7", "is.batch8", "is.batch9", "is.batch10", "is.batch11", "is.batch12", "is.batch13", "is.batch14", "is.batch15", "is.batch16", "is.batch17", "is.batch18", "is.batch19", "is.batch20", "is.batch21", "is.batch22", "is.cpp.box7", "is.cpp.box8", "is.cpp.box10", "is.cpp.box11", "is.cpp.box12")
traitcovs[["act4.t"]]  <- list("one", "sex", "is.cpp.box7", "is.cpp.box8", "is.cpp.box11", "is.cpp.box12")
traitcovs[["act5.t"]]  <- list("one", "sex", "is.gen51", "is.gen53", "is.gen54", "is.gen55", "is.gen56", "is.cpp.box7", "is.cpp.box8", "is.cpp.box10", "is.cpp.box11","is.batch6", "is.batch7", "is.batch13", "is.batch15", "is.batch16", "is.batch19", "is.batch20", "is.batch21", "is.batch22")
traitcovs[["act8.t"]]  <- list("one", "sex","is.gen51","is.gen55", "is.gen56", "is.cpp.box2", "is.cpp.box3", "is.cpp.box4","is.cpp.box5","is.cpp.box6", "is.cpp.box7", "is.cpp.box8","is.cpp.box9","is.cpp.box11", "is.cpp.box12", "is.batch3", "is.batch7")
traitcovs[["sens"]]    <- list("one", "sex", "is.gen52", "is.cpp.box7","is.cpp.box10", "is.batch8")
traitcovs[["wild.binary"]]    <- list("one", "sex", "cpp.age")
traitcovs[["tail"]] <-list("one", "sex", "is.gen52", "is.gen53", "is.gen51", "is.gen56", "rip.weight")
traitcovs[["glucose"]] <- list("one", "sex", "glu.weight", "glu.age")
traitcovs[["is.coatW"]] <- list("one", "sex")


cmds <- c()
for (trait in names(traitcovs)) {

    index.pheno     <- which(pheno.names == trait)


        cmds <- c(cmds, paste0("gemma -g /group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/chrAll.filtered.dosage -p /group/palmer-lab/AIL/LgSm-DataProcessing/phenos.allgeno.txt -k /group/palmer-lab/AIL/qtlmapping/kinship/onlyEmpirical/chrAll.cXX.txt -a /group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/chrAll.filtered.snpinfo -c /group/palmer-lab/AIL/qtlmapping/covariates/", trait, ".emp.covs -lmm 2 -maf ", MAF, " -o ", trait, ".allChr -n ", index.pheno))

}
write.table(cmds, file=paste0("./gemma.herit.emp.cmds"),
            row.names=F, col.names=F, quote=F)
