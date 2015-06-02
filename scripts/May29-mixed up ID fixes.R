### TO FIX EVERYWHERE

### 5-29-15 : fixed in ./covariates.orm, dataFiles/phenotypes/alldata.txt
### removed pheno.allgeno and phenotypes.orm from top folder, to be regenerated
### from the new alldata txt once i'm sure it's fixed.

### fucked up fams
### 46023 & 46024, different parents (but both listed as fam BrF49-24).
id gen      fam   dam  sire
46023.1  50 BrF49-24 45383 45318
46024.1  50 BrF49-24 45392 45446 # should be BrF49-28
### 50210 & 50377. same as above, both listed as BrF51-34
id gen      fam   dam  sire
50210.1  52 BrF51-34 49017 48404 # change parents to 48442 and 49096, as below.
50377.1  52 BrF51-34 48442 49096 # these two mice are sibs from different litters

### fucked up sires
id gen      fam   dam  sire
127 46490.1  50 BrF49-74 45347 45470 # ok
128 46491.1  50 BrF49-74 45347 45470 # ok
131 46640.1  50 BrF49-70 45437 45470 # dam, sire = 45424, 45394
132 46641.1  50 BrF49-70 45437 45470 # "
id gen      fam   dam  sire
390 50225.1  52 BrF51-56 49049 49500 # sire = 49000
392 50230.1  52 BrF51-56 49049 49500 # sire = 49000
491 50716.1  53 BrF52-28 50219 49500 # this should be 50714. 50716 is his sib, a breeder mouse. since they're brothers it really doesn't matter and i'm leaving it as is for now.
494 50723.1  53 BrF52-28 50219 49500 # 49500 is an F52 mouse.
id gen      fam   dam  sire
536 51211.1  53 BrF52-60 50393 50138 # sire = 50124
539 51220.1  53 BrF52-60 50393 50138 # sire = 50124
572 51289.1  53 BrF52-61 50389 50138 # ok
574 51294.1  53 BrF52-61 50389 50138 # ok
id gen      fam   dam  sire
981  54282.1  55 BrF54-77 51957 52184 # ok
986  54291.1  55 BrF54-77 51957 52184 # ok
1060 54853.1  56 BrF55-12 52097 52184 # sire = 52084
1061 54854.1  56 BrF55-12 52097 52184 # sire = 52084

### add rows for this sample that was not included in pheno list because its cpp
### data was messed up. the ppi and glucose data was fine so i included it while
### keeping cpp data as NA.
which(!info$id %in% row.names(ail.kinship)) # 89, 46368, gen 50
info <- info[c(1:88, 90:1124),]

which(grm['49290.1',]>0.05) # sib 49289. gen 51. AIL lib 42. this one looks like a typo in the testing/breeding spreadsheet. it's entered as 48290 even though all the surrounding ids are in the 49000s. so althought 49290 is probably the correct id, i'm leaving it as 48290.

"54386" # flowcell 28 lib 84 - gen 56
which(grm['54386.1',] > 0.05) # 54385 is definitely a sib. 52197 could be a sib, or is at least related. gr between 85,86 = 0.097; with 52197 = 0.051

"57801" # flowcell 25 lib 65 - allegedly generation 54 - should be 51801.
which(grm['57801.1',] > 0.05)
# kinship from the GRM adds up - if 51801, i'd expect her to be most related to the potential sibi=ling, 51794, and that is the case.





