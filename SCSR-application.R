library(glue)
library(devtools)
library(dplyr)
library(stringr)
library(foreach)
library(doParallel)

# Create 4 dirs
baseDir <- '/data/villemin'
package <- '/data2/villemin/SingleCellSignalR/SingleCellSignalR'
Output  <- '/data2/villemin/SingleCellSignalR/Output'
Input   <- '/data2/villemin/SingleCellSignalR/Input'

devtools::install(glue("{baseDir}{package}"))
library(SingleCellSignalR)
#suppressPackageStartupMessages(library(SingleCellSignalR))

cacheInfo()
cacheClear()
cacheInfo()
createRessources(verbose=FALSE)
cacheInfo()

print("reactome")
reactome <-  getRessource(ressource="Reactome")

print("gobp")
GObp     <-  getRessource(ressource="GO-BP")

print("pwc")
pwc     <-  getRessource(ressource="PwC")

stop()
#data <- updatePathwaysFromFile(file=glue("{baseDir}{Input}/c2.cp.reactome.v2022.1.Hs.json"),
						#  pathwaySource="REACTOME")
#head(data)
#tail(data)

#data(example_dataset, package = "SingleCellSignalR")
# data <- as.data.frame(data)
#data <- example_dataset
#object <- dataPrepare(file = data)
#object <- cellClustering(obj = object, n = 10, method = "kmeans")


print("CreateDatabase on Request")
#createDatabase()

print("checkLastVersion FALSE")
#checkDatabaseLastVersion(update=FALSE)

print("checkLastVersion TRUE")
#checkDatabaseLastVersion(update=TRUE)

print("getInteractions 2")

interactions <- getInteractions()
head(interactions)
dim(interactions)

print("getInteractions 1")

interactions <- getInteractions(idRelease=1)
head(interactions)
dim(interactions)

print("intersect - 20  ")
intersect(interactions$Ligand,interactions$Receptor)
# Both Ligand and receptor
# "ADAM17"  "EFNA3"   "NECTIN3" "PLG"    
# "PVR"     "TYROBP"  "ITGB2"   "HLA-F"   "CEACAM1" "PECAM1"    
# "NECTIN1" "CEACAM8"     "CLEC2B"   
# "CD36"  "CD48"  "CD80" "CD177" "CD47" "CD22" "CD1D" 

#getInteractions(idRelease=65)
print("getComplexes")

complexes <- getComplexes()
head(complexes)
dim(complexes)
dim(complexes[complexes$type=="L",])
dim(complexes[complexes$type=="R",])

write.table(interactions ,file = glue("{baseDir}{Output}/interactions.csv"), col.names= TRUE,row.names =  FALSE, quote = F,sep = "\t")
write.table(complexes ,file = glue("{baseDir}{Output}/complexes.csv"), col.names= TRUE,row.names =  FALSE, quote = F,sep = "\t")

#Notes :
########
# Add diff for interactions for LRDB vs OmnipathR::curated_ligand_receptor_interactions() 
# Complexes add from (https://omnipathdb.org/complexes)
#TODO :
########
# Add references (Complexes & LRs)
# Modify source to sources
# Remove Ligand in Receptor
#
name <- "ARRB1_BCL9_BTRC_CBY1_CDC73_DAB2_DKK1_FBXW11_GATA3_GPC3_KREMEN1_KREMEN2_MACF1_MAP3K7_MARK1_MED12_MED23_NKD1_NLK_NPR1_PITX2_PPP2CA_RANBP3_ROR2_SFRP1_SFRP2_SMAD3_SMURF2_STK11_TLE1_WIF1_YWHAB"
name.split <- str_split(name,"_")
print(length(name.split[[1]]))

stoc <- "1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1"
stoc.split <- str_split(stoc,":")
print(length(stoc.split[[1]]))