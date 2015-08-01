# Parse GigaMUGA files (in proprietary Illumina iDAT format) using Bioconductor packages crlmm and illuminaio

biocLite("crlmm")
library("crlmm")
library("illuminaio")
setwd("/group/palmer-lab/AIL/knownSNPs/gigaMUGA")


bpmfile <- readBPM(file="./dataFiles/gigaMUGA/GigaMuga_11769261_A.bpm")

idats <- c()
for (file in filenames){
    idats[[file]] <- readIDAT(paste0("./idats/",file)
}

# There is one iDAT file for each channel (red and green) and each individual (48 files in total)

# Each iDAT file includes:
# fileSize
# versionNumver
# nFields
# fields
    # field | fieldCode | byteOffset | Bytes
# nSNPsRead
# IlluminaID
# SD
# Mean
# NBeads
# MidBlock
# RunInfo
# RedGreen
# MostlyNull
# Barcode
# ChipType
# MostlyA
# Unknown.1, 6, 2, 3, 4, 5, 7

# The most important list entry is Quants, which contains average intensity (Mean), number of beads (NBeads), and a measure of variability (SD) for each bead type on the array.




