
library(glue)
library(devtools)
library(dplyr)
library(stringr)
library(foreach)
library(doParallel)

# Create 4 dirs - Server 
baseDir <- '/data/villemin'
package <- '/data2/villemin/SingleCellSignalR/SingleCellSignalR'
Output  <- '/data2/villemin/SingleCellSignalR/Output'
Input   <- '/data2/villemin/SingleCellSignalR/Input'
# Create 4 dirs - Local 
baseDir <- '/home/jp'
package <- '/Documents/workingDir/SingleCellSignalR'
Output  <- '/Documents/workingDir/Output'
Input   <- '/Documents/workingDir/Input'

devtools::install(glue("{baseDir}{package}"))
library(SingleCellSignalR)
#suppressPackageStartupMessages(library(SingleCellSignalR))

# User can set a different repository for cache
#Sys.setenv("SingleCellSignalR_CACHEDIR"
# = "/data/villemin/data2/villemin/SingleCellSignalR/Input")
#cacheClear()
cacheInfo()
#cacheVersion()
cacheClear()
cacheInfo()
createDatabase(verbose=FALSE)
createResources(verbose=FALSE)
#cacheInfo()

print("getInteractions")
interactions <- getInteractions()
head(interactions)
dim(interactions)

print("getComplexes")
complexes <- getComplexes()
head(complexes)

stop()

# 1 ##############################



#createResources()


#cacheInfo()

#checkCacheLastVersion()

# 2 ##############################

reactome <-  getResource(resourceName="Reactome")
GObp     <-  getResource(resourceName="GO-BP")
pwc      <-  getResource(resourceName="PwC")

head(pwc)

dataReactomeJson <- updatePathwaysFromFile(
			file=glue("{baseDir}{Input}/c2.cp.reactome.v2022.1.Hs.json"),
						  resourceName="Reactome",fileType="json")

dataReactomeGmt <- updatePathwaysFromFile(
				file=glue("{baseDir}{Input}/ReactomePathways.gmt"),
						  resourceName="Reactome",fileType="gmt")

head(dataReactomeGmt)


print("getInteractions")
interactions <- getInteractions()
head(interactions)
dim(interactions)

print("getInteractions for older release v1.")
interactions <- getInteractions(idRelease=1)
head(interactions)
dim(interactions)

print("getComplexes")
complexes <- getComplexes()
head(complexes)

stop()
##############################
###   Vrac                 ###
##############################

print("intersect - 20  ")
intersect(interactions$Ligand,interactions$Receptor)
# Both Ligand and receptor
# "ADAM17"  "EFNA3"   "NECTIN3" "PLG"    
# "PVR"     "TYROBP"  "ITGB2"   "HLA-F"   "CEACAM1" "PECAM1"    
# "NECTIN1" "CEACAM8"     "CLEC2B"   
# "CD36"  "CD48"  "CD80" "CD177" "CD47" "CD22" "CD1D" 

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