######################################################
## Mapping of Selleck and LATCA Compounds to ChEMBL ##
######################################################
## Author: Thomas Girke
## Last update: 16-Feb-2024

## Get for Selleck CAS IDs annotations from DrugBank 

## Set run environment from command-line
srun --x11 --partition=girkelab --mem=20gb --cpus-per-task 8 --ntasks 1 --time 20:00:00 --pty bash -l
module load openbabel
cd /rhome/tgirke/Projects/SeanCutler/sean_cutler_Jun2022
vim retrieveCMPids.R # -> start R with \rf
system("hostname")

## (a) Import Selleck table from Sean
library("tidyverse")
slk <- readr::read_tsv("./data/20220301-L1300-FDA-approved-Drug-Library-96-well.tsv")

## (b) DrugBank downloads are here: https://go.drugbank.com/releases/latest
## They require username and password which are provided in wget command below.
## To re-execute download, uncomment next lines
# system("wget https://go.drugbank.com/releases/5-1-9/downloads/all-drug-links --user=thomas.girke@ucr.edu --password=1drug_access2 -O ./data/all-drug-links.zip")
# system("unzip ./data/all-drug-links.zip -d ./data && mv ./data/drug\\ links.csv ./data/drugbank_annot.csv")
db_annot <- readr::read_csv("./data/drugbank_annot.csv")
colnames(db_annot) <- paste("DB", colnames(db_annot), sep="_") # To trace data sources in result table, here DrugBank wiht 'DB_'

## (c) DrugBank data is added using first CAS ID. If not available, match on compound Name column will used instead. 
slk_index <- bind_cols(Row_Index=1:nrow(slk), slk[, c("CAS Number", "Name")])
cas <- setNames(1:length(db_annot$'DB_CAS Number'), db_annot$'DB_CAS Number')
name <- setNames(1:length(db_annot$'DB_Name'), gsub(" .*", "", db_annot$'DB_Name'))
slk_index <- bind_cols(slk_index, 
                       CAS_DrugBank=names(cas[slk_index$'CAS Number']), 
                       CAS_DrugBank_RowNo=cas[slk_index$'CAS Number'], 
                       Name_DrugBank=names(name[gsub(" .*", "", slk_index$'Name')]),
                       Name_DrugBank_RowNo=name[gsub(" .*", "", slk_index$'Name')])
Row_Index <- ifelse(!is.na(slk_index$CAS_DrugBank), slk_index$CAS_DrugBank_RowNo, slk_index$Name_DrugBank_RowNo)
slk_index <- bind_cols(slk_index, RowIndex=Row_Index) 
slk_annot <- bind_cols(slk, db_annot[slk_index$RowIndex,])

## (d) Get CAS IDs from ChEBI that are missing in DrugBank 
# system("wget https://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/database_accession.tsv -O ./data/chebi_annot.tsv")
# system("wget https://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/names.tsv.gz -O ./data/chebi_names.tsv.gz && gunzip ./data/chebi_names.tsv.gz")
chebi_annot <- readr::read_tsv("./data/chebi_annot.tsv")
chebi_annot <- chebi_annot[chebi_annot$TYPE=="CAS Registry Number", ] # keep only CAS entries
chebi_annot <- chebi_annot[!duplicated(chebi_annot$ACCESSION_NUMBER), ] # remove duplicates
chebi_names <- readr::read_tsv("./data/chebi_names.tsv")
chebi_names <- chebi_names[chebi_names$TYPE=='IUPAC NAME',] # keep only IUPAC name
chebi_id_vec <- tapply(chebi_names$NAME, chebi_names$COMPOUND_ID, paste, collapse=", ")
cas_chebi_vec <- tapply(paste0("CHEBI:", chebi_annot$COMPOUND_ID), chebi_annot$ACCESSION_NUMBER, paste, collapse=", ")
chebi_annot_final <- bind_cols(chebi_annot, 
                        rename(as_tibble(cas_chebi_vec[as.character(chebi_annot$ACCESSION_NUMBER)]), CHEBI_ID=value),
                        rename(as_tibble(chebi_id_vec[as.character(chebi_annot$COMPOUND_ID)]), NAME=value)
)
chebi_annot_final <- as.data.frame(chebi_annot_final)
colnames(chebi_annot_final) <- paste("ChEBI", colnames(chebi_annot_final), sep="_")
row.names(chebi_annot_final) <- chebi_annot_final$ChEBI_ACCESSION_NUMBER

## (e) Add ChEBI columns
slk_annot <- bind_cols(slk_annot, chebi_annot_final[slk$'CAS Number', !colnames(chebi_annot_final) %in% "ChEBI_ACCESSION_NUMBER"])
write.table(slk_annot, "./results/selleck_annot.tsv", row.names=FALSE, quote=FALSE, sep="\t")

## (d) Now remaining CAS to PubChem mappings are obtained with webchem where
## fingerprint similarity is used as secondary query to resove ambiguities...
slk_annot <- readr::read_tsv("./results/selleck_annot.tsv")
idMA <- slk_annot[, c("DB_DrugBank ID", "DB_PubChem Compound ID", "DB_PubChem Substance ID", "ChEBI_CHEBI_ID")]
missing_index <- rowSums(is.na(idMA)) == ncol(idMA)
cas_missing <- slk_annot[missing_index, "CAS Number"]$'CAS Number'
cas_missing <- gsub("(,.*)|(\\(.*)|( \\(.*)", "", cas_missing) # Cleans up 7 mixed CAS IDs that are comma separated strings

## PubChem IDs for CAS IDs without matches
library(webchem); library(ChemmineR); library(ChemmineOB)
caspubids <- get_cid(cas_missing, from = "xref/rn")
readr::write_tsv(caspubids, "./results/missing_caspubids.txt") # Creates intermediate file to avoid redownload which is slow
caspubids <- readr::read_tsv("./results/missing_caspubids.txt")
pcids <- as.numeric(caspubids$cid)
pcids <- unique(pcids[!is.na(pcids)])
sdfset <- pubchemCidToSDF(as.numeric(pcids))
sdfset <- sdfset[validSDF(sdfset)]
write.SDF(sdfset, "./results/missing.SDF")
sdfset <- read.SDFset("./results/missing.SDF")
cid(sdfset) <- sdfid(sdfset)
apset <- sdf2ap(sdfset) 
invalid_index <- cid(apset) %in% names(which(sapply(as(apset, "list"), length)==1))
apset <- apset[!invalid_index]

## Create atom pair db for SMILES string from slk
smiles_slk <- setNames(slk_annot$SMILES, slk_annot$'CAS Number')
smiles_slk <- smiles_slk[!grepl("\\||A", smiles_slk)] # Removes misformated SMILES containing '|' or 'A'
smiles_slk <- smiles_slk[!is.na(smiles_slk)]
smiles_slk <- smiles_slk[-c(1676, 2348, 2682, 2760, 3005)] # contain some bad formatting
sdf_slk <- smiles2sdf(smiles_slk) 
sdf_slk <- sdf_slk[validSDF(sdf_slk)]
sdf_slk <- sdf_slk[-c(1696,2882,3010)] # removes compounds with invalid atoms
ap_slk <- sdf2ap(sdf_slk) 

## Perform atom pair searches
caspubids <- caspubids[caspubids$cid %in% cid(apset),]
caspubids_list <- tapply(caspubids$cid, caspubids$query, as.character)
caspubids_list <- caspubids_list[names(caspubids_list) %in% cid(ap_slk)]
search_results <- sapply(names(caspubids_list), function(x) cmp.search(apset[caspubids_list[[x]]], ap_slk[x], type=2, cutoff=Inf, quiet=TRUE)) 

best_match <- sapply(names(search_results), function(x) search_results[[x]] == 1)
best_match <- sapply(names(best_match), function(x) replace(best_match[[x]], 1, TRUE)) # assures that if best match is less than 1 than keep this one
best_match <- sapply(names(best_match), function(x) paste(names(best_match[[x]][best_match[[x]]]), collapse=", "))
best_match_DF <- data.frame(cas=names(search_results), PCID=best_match)
sdf_annot_ma <- datablock2ma(datablocklist=datablock(sdfset))[, c("PUBCHEM_COMPOUND_CID", "PUBCHEM_IUPAC_NAME", "PUBCHEM_IUPAC_INCHIKEY")] 
cas2pcid <- sapply(rownames(best_match_DF), function(x) strsplit(best_match_DF[x, "PCID"], ", "))
cas2pcidDF <- as.data.frame(t(sapply(names(cas2pcid), function(x) apply(sdf_annot_ma[cas2pcid[[x]], ,drop=FALSE], 2, paste, collapse=", "))))
colnames(cas2pcidDF) <- paste("WC_Sim", colnames(cas2pcidDF), sep="_")

slk_annot <- bind_cols(slk_annot, cas2pcidDF[slk$'CAS Number', ])
write.table(slk_annot, "./results/selleck_annot.tsv", row.names=FALSE, quote=FALSE, sep="\t")
# -> note this file has been uploaded to a Google Sheet here: https://bit.ly/3AhpFNm

## Summary stats of matches
MatchStatsDF <- data.frame(Method=c("DrugBank", "ChEBI", "WC_Sim"), 
                    CAS_ID=c(sum(!is.na(slk_annot$'DB_CAS Number')), sum(!is.na(slk_annot$'ChEBI_ID')), "Remaining"), 
                    Name=c(sum(!is.na(slk_annot$DB_Name)), sum(!is.na(slk_annot$'ChEBI_NAME')), "Remaining"),
                    Total_Match=c(sum(!is.na(slk_annot$'DB_CAS Number') | !is.na(slk_annot$DB_Name)), 
                                c(sum(!is.na(slk_annot$'DB_CAS Number') | !is.na(slk_annot$DB_Name) | !is.na(slk_annot$ChEBI_ID))), 
                                c(sum(!is.na(slk_annot$'DB_CAS Number') | !is.na(slk_annot$DB_Name) | !is.na(slk_annot$ChEBI_ID) | !is.na(slk_annot$WC_Sim_PUBCHEM_COMPOUND_CID)))))
write.table(MatchStatsDF, "./results/MatchStatsDF.tsv", row.names=FALSE, quote=FALSE, sep="\t")

######################
## UniChem Mappings ##
######################

## (1) UniChem downloads
## Unichem FTP page: https://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/
## The number to database name mappings can be looked up here: https://www.ebi.ac.uk/unichem/legacy/ucquery/listSources

## (a) DrugBank to ChEMBL: src1src2.txt.gz
# system("wget https://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id1/src1src2.txt.gz -O ./data/src1src2.txt.gz && gunzip ./data/src1src2.txt.gz")
chembl2db <- readr::read_tsv("./data/src1src2.txt")
chembl2db <- setNames(chembl2db$"From src:'1'", chembl2db$"To src:'2'")

## (b) ChEBI to ChEMBL: src1src7.txt
# system("wget https://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id1/src1src7.txt.gz -O ./data/src1src7.txt.gz && gunzip ./data/src1src7.txt.gz")
chembl2chebi <- readr::read_tsv("./data/src1src7.txt")
chembl2chebi <- setNames(chembl2chebi$"From src:'1'", as.character(chembl2chebi$"To src:'7'")) 

## (c) PubChem to ChEMBL: src1src22.txt.gz
# system("wget https://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id1/src1src22.txt.gz -O ./data/src1src22.txt.gz && gunzip ./data/src1src22.txt.gz")
chembl2pc <- readr::read_tsv("./data/src1src22.txt")
chembl2pc <- setNames(chembl2pc$"From src:'1'", as.character(chembl2pc$"To src:'22'")) 

## (2) Merge selleck_annot.tsv with ChEMBL IDs from UniChem
library(dplyr)
slk_annot <- readr::read_tsv("./results/selleck_annot.tsv")
chembl_ids <- mutate(slk_annot[,c(1:2, 13)], DrugBank_ID=chembl2db[slk_annot$'DB_DrugBank ID'], ChEBI_ID=chembl2chebi[as.character(slk_annot$'ChEBI_COMPOUND_ID')])

## WC PubChem ID column needs special treatment since many are concatenation of several IDs
pc_list <- strsplit(as.character(slk_annot$'WC_Sim_PUBCHEM_COMPOUND_CID'), ", ")
pc_ids <- sapply(seq_along(pc_list), function(x) pc_list[[x]][pc_list[[x]] %in% names(chembl2pc)][1])
chembl_ids <- mutate(chembl_ids, PubChem_ID=chembl2pc[pc_ids])

final_chembl_id <- ifelse(!is.na(chembl_ids$DrugBank_ID), chembl_ids$DrugBank_ID, chembl_ids$ChEBI_ID)
final_chembl_id <- ifelse(!is.na(final_chembl_id), final_chembl_id, chembl_ids$PubChem_ID)
chembl_ids <- mutate(chembl_ids, Final_ChEMBL_ID=final_chembl_id)

#########################################
## ChEMBL Annotations, Properties, etc ##
#########################################
library(RSQLite); library(dplyr)
mydb <- dbConnect(SQLite(), "~/Projects/REU_22/REU2022-Lupe/Project/data/chembl_30/chembl_30_sqlite/chembl_30.db")
dbListTables(mydb)
id_lookup <- dbGetQuery(mydb, 'SELECT * FROM chembl_id_lookup')
id_lookup <- setNames(as.character(id_lookup$entity_id), id_lookup$chembl_id) 
mdict <- dbGetQuery(mydb, 'SELECT * FROM molecule_dictionary')
mtype <- setNames(as.character(mdict$molecule_type), as.character(mdict$molregno)) 
mdict <- setNames(as.character(mdict$pref_name), as.character(mdict$molregno)) 
mech <- dbGetQuery(mydb, 'SELECT * FROM drug_mechanism')
mech <- setNames(as.character(mech$mechanism_of_action), as.character(mech$molregno)) 

molreg <- id_lookup[unique(chembl_ids$Final_ChEMBL_ID)] 
molreg <- molreg[!is.na(molreg)] 
chembl_annot <- data.frame(ChEMBL_Name=mdict[molreg], MOA=mech[molreg])
rownames(chembl_annot) <- names(molreg)
chembl_annot <- mutate(chembl_ids, chembl_annot[chembl_ids$Final_ChEMBL_ID,])
write.table(chembl_annot, "./results/chembl_annot.tsv", row.names=FALSE, quote=FALSE, sep="\t")
# -> note the columns from the chembl_annot table were copied into the Google Sheet here https://bit.ly/3J6pWFd

## Alternative name-based method
mdict <- dbGetQuery(mydb, 'SELECT * FROM molecule_dictionary')
mdict <- mdict[!is.na(mdict$pref_name),]
prefname <- setNames(mdict$chembl_id, mdict$pref_name)
chemblid2name <- setNames(mdict$pref_name, mdict$chembl_id)
index <- gsub(" \\(.*\\)", "", toupper(chembl_annot$Name))
index <- gsub(" HCL", " HYDROCHLORIDE", index, ignore.case=TRUE)
index <- gsub(" DIHCL", " DIHYDROCHLORIDE", index, ignore.case=TRUE)
index <- gsub(" hemifumarate", " FUMARATE", index, ignore.case=TRUE)
index <- gsub(".*elafloxacin.*", "DELAFLOXACIN", index, ignore.case=TRUE)
prefname <- prefname[index]
chembl_annot <- mutate(chembl_annot, Final_ChEMBL_ID2=as.character(prefname), ChEMBL_Name2=names(prefname))
chembl_annot$Final_ChEMBL_ID[is.na(chembl_annot$Final_ChEMBL_ID)] <- chembl_annot$Final_ChEMBL_ID2[is.na(chembl_annot$Final_ChEMBL_ID)]
chembl_annot$Final_ChEMBL_ID2[is.na(chembl_annot$Final_ChEMBL_ID2)] <- chembl_annot$Final_ChEMBL_ID[is.na(chembl_annot$Final_ChEMBL_ID2)]
chembl_annot[,"ChEMBL_Name"] <- chemblid2name[as.character(chembl_annot$Final_ChEMBL_ID)]
chembl_annot[,"ChEMBL_Name2"] <- chemblid2name[as.character(chembl_annot$Final_ChEMBL_ID2)]
write.table(chembl_annot, "./results/chembl_annot.tsv", row.names=FALSE, quote=FALSE, sep="\t")
dbDisconnect(mydb) 

## Create master table with parent molecule hierarchies
library(RSQLite); library(dplyr)
mydb <- dbConnect(SQLite(), "~/Projects/REU_22/REU2022-Lupe/Project/data/chembl_30/chembl_30_sqlite/chembl_30.db")
mhier <- dbGetQuery(mydb, 'SELECT * FROM molecule_hierarchy')
mhier <- mhier[order(mhier$parent_molregno, mhier$molregno),c(2,1)]
rownames(mhier) <- NULL
id_lookup <- dbGetQuery(mydb, 'SELECT * FROM chembl_id_lookup')
id_lookup <- id_lookup[id_lookup$entity_type=="COMPOUND",]
id_lookup <- setNames(as.character(id_lookup$chembl_id), id_lookup$entity_id) 
mdict <- dbGetQuery(mydb, 'SELECT * FROM molecule_dictionary')
mtype <- setNames(as.character(mdict$molecule_type), as.character(mdict$molregno)) 
mdict <- setNames(as.character(mdict$pref_name), as.character(mdict$molregno)) 
msyn <- dbGetQuery(mydb, 'SELECT * FROM molecule_synonyms')
msyn <- tapply(msyn$synonyms, msyn$molregno, function(x) paste(unique(x), collapse=", "))
cmpprop <- dbGetQuery(mydb, 'SELECT * FROM compound_properties')
rownames(cmpprop) <- as.character(cmpprop$molregno)
mycol <- c(c("mw_freebase", "full_mwt"), colnames(cmpprop)[!colnames(cmpprop) %in% c("mw_freebase", "full_mwt")])
cmpprop <- cmpprop[, mycol[-3]]
# cmpprop <- cmpprop[,c("mw_freebase", "full_mwt")]
masterDF <- data.frame(mhier, ChEMBL_ID=id_lookup[as.character(mhier$molregno)],
                             ChEMBL_NAME=mdict[as.character(mhier$molregno)],
                             ChEMBL_SYNONYMS=msyn[as.character(mhier$molregno)],
                             Mol_Type=mtype[as.character(mhier$molregno)],
                             cmpprop[as.character(mhier$molregno),]) 
masterDF[masterDF$parent_molregno=="451591",]
masterDF[masterDF$ChEMBL_ID=="CHEMBL2375172",]
write.table(masterDF, "./results/masterDF.tsv", row.names=FALSE, quote=FALSE, sep="\t")
masterDF <- read.delim("./results/masterDF.tsv", sep="\t", quote="")
rownames(masterDF) <- masterDF$ChEMBL_ID

chembl_annot2 <- chembl_annot[,c(1:3,10:11)]
masterDFsub <- masterDF[chembl_annot2$Final_ChEMBL_ID2,]
masterDFsub2 <- masterDF[masterDF$parent_molregno %in% masterDFsub$parent_molregno, ]
mycounts <- table(masterDF$parent_molregno)
mycounts <- mycounts[rep(names(mycounts), mycounts)]
names(mycounts) <- masterDF$ChEMBL_ID
rowindex <- mycounts[chembl_annot2$Final_ChEMBL_ID2]
names(rowindex) <- seq_along(rowindex)
rowindex[is.na(rowindex)] <- 1
rowindex <- as.integer(rep(names(rowindex), rowindex))
chembl_annot_exp <- chembl_annot2[rowindex,]

chembl2parent <- setNames(masterDF$parent_molregno, masterDF$ChEMBL_ID)
chembl2parent <- chembl2parent[chembl_annot2$Final_ChEMBL_ID2]
line_index <- tapply(seq_along(masterDF$parent_molregno), as.factor(masterDF$parent_molregno), paste, collapse=", ")
line_index <- line_index[as.character(chembl2parent)]
masterDF_sub <- masterDF[as.integer(unlist(strsplit(line_index, ", "))),] 
chembl_annot_hierarchy <- cbind(chembl_annot_exp, masterDF_sub)
struct <- dbGetQuery(mydb, 'SELECT * FROM compound_structures')
rownames(struct) <- as.character(struct$molregno)
struct <- struct[,-c(1,2)]
chembl_annot_hierarchy <- data.frame(chembl_annot_hierarchy, struct[as.character(chembl_annot_hierarchy$molregno),])
write.table(chembl_annot_hierarchy, "./results/chembl_annot_hierarchy.tsv", row.names=FALSE, quote=FALSE, sep="\t")

###############################################
## Compute Molecular Descriptors for Selleck ##
###############################################

## Descriptors from PubChem
chembl_annot_hierarchy <- readr::read_tsv("./results/chembl_annot_hierarchy.tsv", col_type=cols(.default = "c")) # Param for last argument forces to read all columns as character 
## Report ChEMBL to PubChem duplicates
chembl2pc <- readr::read_tsv("./data/src1src22.txt", col_type=cols(.default = "c"))
dups <- tapply(chembl2pc$"To src:'22'", as.factor(chembl2pc$"From src:'1'"), paste, collapse=", ")
chembl_annot_hierarchy <- mutate(chembl_annot_hierarchy, PC_ID_dups=dups[chembl_annot_hierarchy$ChEMBL_ID])
library(ChemmineR)
pcid <- as.numeric(gsub(", .*", "", chembl_annot_hierarchy$PC_ID_dups))
pcid <- unique(pcid[!is.na(pcid)])
chembl2pc <- setNames(as.character(chembl2pc$"From src:'1'"), as.character(chembl2pc$"To src:'22'")) 
chembl2pc_id <- chembl2pc[as.character(pcid)]
sdfset <- pubchemCidToSDF(pcid)
for(i in seq_along(sdfset)) sdfset[[i]][[1]][1] <- as.character(chembl2pc_id)[i]
sdfset <- sdfset[validSDF(sdfset)]
write.SDF(sdfset, "./results/selleck_pubchem.sdf")
sdfset <- read.SDFset("./results/selleck_pubchem.sdf")
allids <- sdfid(sdfset)
cid(sdfset) <- allids
PCproperties <- datablock2ma(datablocklist=datablock(sdfset))
PCproperties <- cbind(ChEMBL_ID=row.names(PCproperties), PCproperties)
write.table(PCproperties, "./results/selleck_desc_PC.tsv", sep="\t", quote=FALSE, row.names=FALSE)

## Descriptors from ChemmineR
library("tidyverse"); library(ChemmineR)
slk <- readr::read_tsv("./data/20220301-L1300-FDA-approved-Drug-Library-96-well.tsv", col_type=cols(.default = "c"))
slk_smi <- setNames(slk$SMILES, slk$Cat)
# the following removes 18 misformed smiles that cannot be converted into SDFs 
slk_smi <- slk_smi[!grepl("\\||A", slk_smi)] # Removes misformated SMILES containing '|' or 'A'
slk_smi <- slk_smi[!is.na(slk_smi)]
slk_smi <- slk_smi[-c(1676, 2348, 2682, 2760, 3005)] # contain some bad formatting
sdfset <- smiles2sdf(slk_smi) 
write.SDF(sdfset, "./results/selleck_smiles.sdf")
sdfset <- read.SDFset("./results/selleck_smiles.sdf")
allids <- sdfid(sdfset)
cid(sdfset) <- allids
valid <- validSDF(sdfset); sdfset <- sdfset[valid]
CMproperties <- data.frame(atomcountMA(sdfset, addH=TRUE),
                    MF=MF(sdfset, addH=TRUE),
                    MW=MW(sdfset, addH=TRUE),
                    ChemmineR::groups(sdfset),
                    rings(sdfset, type="count", upper=6, arom=TRUE))
CMproperties <- cbind(Selleck_ID=row.names(CMproperties), CMproperties)
write.table(CMproperties, "./results/selleck_desc_CM.tsv", sep="\t", quote=FALSE, row.names=FALSE)

## Compute molecular descriptors from rcdk
## Currently, same/similar descriptors as in ChEMB that used CDK
## Description of available CDK descriptors is here: 
##     HTML: https://bit.ly/3UVn2su
##     XML: https://bit.ly/3rzhgzg
##     PDF vignette: http://cran.nexr.com/web/packages/rcdk/vignettes/rcdk.pdf

library("tidyverse"); library(rcdk)
slk <- readr::read_tsv("./data/20220301-L1300-FDA-approved-Drug-Library-96-well.tsv")
slk_smi <- setNames(slk$SMILES, slk$Cat)
slk_smi <- parse.smiles(slk_smi)
slk_smi <- slk_smi[sapply(slk_smi, length)!=0] # remove compounds that failed to convert
## Get all possible descriptor names; note some result in seqfault and will not work on many compounds 
# descNames <- unique(unlist(sapply(get.desc.categories(), get.desc.names))) # retrieves unique set of discriptor names
descNames <- c('SmallRingDescriptor', 'AromaticAtomsCountDescriptor', 'AromaticBondsCountDescriptor',
               'APolDescriptor', 'HBondDonorCountDescriptor', 'HBondAcceptorCountDescriptor',
               'ALOGPDescriptor', 'MannholdLogPDescriptor', 'XLogPDescriptor', 'WeightDescriptor',
               'FractionalPSADescriptor', 'TPSADescriptor', 'RotatableBondsCountDescriptor',
               'RuleOfFiveDescriptor', 'AcidicGroupCountDescriptor')
descNames <- paste0("org.openscience.cdk.qsar.descriptors.molecular.", descNames)
CDKproperties <- eval.desc(slk_smi, descNames) # Computes all descriptors, useful description of them is here: https://bit.ly/3UVn2su
formula <- sapply(names(slk_smi), function(x) get.mol2formula(slk_smi[[x]])@string)
CDKproperties <- data.frame(CDKproperties, Formula=formula[rownames(CDKproperties)])
CDKproperties <- cbind(Selleck_ID=row.names(CDKproperties), CDKproperties)
write.table(CDKproperties, "./results/selleck_desc_CDK.tsv", sep="\t", quote=FALSE, row.names=FALSE)

## Merge all Selleck descriptor results into one table
chembl_annot_hierarchy <- readr::read_tsv("./results/chembl_annot_hierarchy.tsv", col_type=cols(.default = "c")) # Param for last argument forces to read all columns as character 
PCproperties <- readr::read_tsv("./results/selleck_desc_PC.tsv", col_type=cols(.default = "c"))
row_index_chembl <- setNames(seq_along(PCproperties$ChEMBL_ID), PCproperties$ChEMBL_ID)
row_index_chembl <- row_index_chembl[chembl_annot_hierarchy$ChEMBL_ID]
PCproperties <- PCproperties[row_index_chembl,]
PCproperties <- rename(PCproperties, ChEMBL_ID_PubChem=ChEMBL_ID)
CMproperties <- readr::read_tsv("./results/selleck_desc_CM.tsv", col_type=cols(.default = "c"))
row_index_selleck_cm <- setNames(seq_along(CMproperties$Selleck_ID), CMproperties$Selleck_ID)
row_index_selleck_cm <- row_index_selleck_cm[chembl_annot_hierarchy$Cat]
CMproperties <- CMproperties[row_index_selleck_cm,]
CMproperties <- rename(CMproperties, Selleck_ID_ChemmineR=Selleck_ID, MW_CM=MW)
CDKproperties <- readr::read_tsv("./results/selleck_desc_CDK.tsv", col_type=cols(.default = "c"))
row_index_selleck_cdk <- setNames(seq_along(CDKproperties$Selleck_ID), CDKproperties$Selleck_ID)
row_index_selleck_cdk <- row_index_selleck_cdk[chembl_annot_hierarchy$Cat]
CDKproperties <- CDKproperties[row_index_selleck_cdk,]
CDKproperties <- rename(CDKproperties, Selleck_ID_CDK=Selleck_ID)
selleck_all <- bind_cols(chembl_annot_hierarchy, PCproperties, CMproperties, CDKproperties)
write_tsv(selleck_all, "./results/selleck_all.tsv")
mycol_names <- as_tibble(colnames(selleck_all))
mycol_names <- rename(mycol_names, Column_Names=value)
write_tsv(mycol_names, "./results/selleck_all_cols.tsv")

#####################
## LATCA Compounds ##
#####################

## LATCA table downloaded from here: https://bit.ly/3pUGNlR
library("tidyverse")

latca <- readr::read_tsv("./data/LATCA_subset-latca_w_props.tsv", col_type="cccc")
## UniChem Mappings: PubChem to ChEMBL: src1src22.txt.gz
# system("wget https://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id1/src1src22.txt.gz -O ./data/src1src22.txt.gz && gunzip ./data/src1src22.txt.gz")
chembl2pc <- readr::read_tsv("./data/src1src22.txt", col_type="cc")
chembl2pc <- setNames(as.character(chembl2pc$"From src:'1'"), as.character(chembl2pc$"To src:'22'")) 
## Name based mappings
library(RSQLite); library(dplyr)
mydb <- dbConnect(SQLite(), "~/Projects/REU_22/REU2022-Lupe/Project/data/chembl_30/chembl_30_sqlite/chembl_30.db")
mdict <- dbGetQuery(mydb, 'SELECT * FROM molecule_dictionary')
mdict <- mdict[!is.na(mdict$pref_name),]
prefname <- setNames(mdict$chembl_id, mdict$pref_name)
chemblid2name <- setNames(mdict$pref_name, mdict$chembl_id)
## To get most complete matches, the following queries hierarchially in this order cid, cat and name
latca_anno <- mutate(latca, ChEMBL_cid=chembl2pc[latca$CID], ChEMBL_cat=chembl2pc[latca$cat], ChEMBL_name=prefname[toupper(latca$name)])
latca_anno <- mutate(latca_anno, ChEMBL_ID_int=ifelse(!is.na(latca_anno$ChEMBL_cid), latca_anno$ChEMBL_cid, latca_anno$ChEMBL_cat))
latca_anno <- mutate(latca_anno, ChEMBL_ID=ifelse(!is.na(latca_anno$ChEMBL_ID_int), latca_anno$ChEMBL_ID_int, latca_anno$ChEMBL_name))

id_lookup <- dbGetQuery(mydb, 'SELECT * FROM chembl_id_lookup')
id_lookup <- setNames(as.character(id_lookup$entity_id), id_lookup$chembl_id) 
latca_anno <- mutate(latca_anno, molregno=id_lookup[latca_anno$ChEMBL_ID])

mdict <- dbGetQuery(mydb, 'SELECT * FROM molecule_dictionary')
mtype <- setNames(as.character(mdict$molecule_type), as.character(mdict$molregno)) 
latca_anno <- mutate(latca_anno, Mol_Type=mtype[latca_anno$molregno])

cmpprop <- dbGetQuery(mydb, 'SELECT * FROM compound_properties')
rownames(cmpprop) <- as.character(cmpprop$molregno)
mycol <- c(c("mw_freebase", "full_mwt"), colnames(cmpprop)[!colnames(cmpprop) %in% c("mw_freebase", "full_mwt")])
cmpprop <- cmpprop[, mycol[-3]]
latca_anno <- mutate(latca_anno, cmpprop[latca_anno$molregno,])

struct <- dbGetQuery(mydb, 'SELECT * FROM compound_structures')
rownames(struct) <- as.character(struct$molregno)
struct <- struct[,-c(1,2)]
latca_anno <- mutate(latca_anno, struct[as.character(latca_anno$molregno),])
latca_anno <- relocate(latca_anno, LATCA_ID=coord) 

write_tsv(latca_anno, "./results/latca/latca_anno.tsv")

## Attempt to get alternative CIDs for LATCA via InChI 
# library(ChemmineR)
# latca <- readr::read_tsv("./data/LATCA_subset-latca_w_props.tsv", col_type="cccc")
# cid <- as.numeric(latca$CID)
# sdfset <- pubchemCidToSDF(cid)
# write.SDF(sdfset, "./results/latca/latca.sdf")
# inchi <- datablock2ma(datablock(sdfset))[,"PUBCHEM_IUPAC_INCHI"]
# cidvec <- pubchemInchi2cid(inchi) # this is slow since the queries have to run one by one
# df <- data.frame(CID=latca_anno$CID, inchi=cidvec)
# df <- df[df$inchi!=0,]
# chembl2pc[df[df[,1]!=df[,2],"inchi"]]
##-> the above did only return one pcid for each inchi. Some of the pcids were different 
## than the original one from LATCA. Only 3 of them returned additional
## ChEMBL Ids but those were not included here since they were all from unnamed compounds.

#############################################
## Compute Molecular Descriptors for LATCA ##
#############################################

## Descriptors from PubChem and ChemmineR 
library(ChemmineR)
latca <- readr::read_tsv("./data/LATCA_subset-latca_w_props.tsv", col_type="cccc")
pcid <- as.numeric(latca$CID)
sdfset <- pubchemCidToSDF(pcid)
for(i in seq_along(sdfset)) sdfset[[i]][[1]][1] <- latca$coord[i]
write.SDF(sdfset, "./results/latca/latca.sdf")
sdfset <- read.SDFset("./results/latca/latca.sdf")
allids <- sdfid(sdfset)
cid(sdfset) <- allids
PCproperties <- datablock2ma(datablocklist=datablock(sdfset))
PCproperties <- cbind(LATCA_ID=row.names(PCproperties), PCproperties)
write.table(PCproperties, "./results/latca/latca_desc_PC.tsv", sep="\t", quote=FALSE, row.names=FALSE)
valid <- validSDF(sdfset); sdfset <- sdfset[valid]
CMproperties <- data.frame(atomcountMA(sdfset, addH=FALSE),
                    MF=MF(sdfset),
                    MW=MW(sdfset),
                    ChemmineR::groups(sdfset),
                    rings(sdfset, type="count", upper=6, arom=TRUE))
CMproperties <- CMproperties[allids,] # adds rows of invalid CMPs back in
rownames(CMproperties) <- allids
CMproperties <- cbind(LATCA_ID=row.names(CMproperties), CMproperties)
write.table(CMproperties, "./results/latca/latca_desc_CM2.tsv", sep="\t", quote=FALSE, row.names=FALSE)

## Descriptors from rcdk
library("tidyverse"); library(rcdk)
sdfset <- read.SDFset("./results/latca/latca.sdf")
allids <- sdfid(sdfset)
slk_sdf <- load.molecules("./results/latca/latca.sdf")
names(slk_sdf) <- allids
slk_sdf <- slk_sdf[sapply(slk_sdf, length)!=0] # remove compounds that failed to convert
## Get all possible descriptor names; note some result in seqfault and will not work on many compounds 
# descNames <- unique(unlist(sapply(get.desc.categories(), get.desc.names))) # retrieves unique set of discriptor names
descNames <- c('SmallRingDescriptor', 'AromaticAtomsCountDescriptor', 'AromaticBondsCountDescriptor',
               'APolDescriptor', 'HBondDonorCountDescriptor', 'HBondAcceptorCountDescriptor',
               'ALOGPDescriptor', 'MannholdLogPDescriptor', 'XLogPDescriptor', 'WeightDescriptor',
               'FractionalPSADescriptor', 'TPSADescriptor', 'RotatableBondsCountDescriptor',
               'RuleOfFiveDescriptor', 'AcidicGroupCountDescriptor')
descNames <- paste0("org.openscience.cdk.qsar.descriptors.molecular.", descNames)
CDKproperties <- eval.desc(slk_sdf, descNames) # Computes all descriptors, useful description of them is here: https://bit.ly/3UVn2su
formula <- sapply(allids, function(x) get.mol2formula(slk_sdf[[x]])@string)
CDKproperties <- data.frame(CDKproperties, Formula=formula[rownames(CDKproperties)])
CDKproperties <- cbind(LATCA_ID=row.names(CDKproperties), CDKproperties)
write.table(CDKproperties, "./results/latca/latca_desc_CDK2.tsv", sep="\t", quote=FALSE, row.names=FALSE)

## Merge all LATCA descriptor results into one table
latca_anno <- readr::read_tsv("./results/latca/latca_anno.tsv", col_type=cols(.default = "c")) # Param for last argument forces to read all columns as character 
latca_anno <- select(latca_anno, -(ChEMBL_cid : ChEMBL_ID_int)) # Removes some unnecessary internal checking columns
sortbyid <- latca_anno$LATCA_ID
PCproperties <- readr::read_tsv("./results/latca/latca_desc_PC.tsv", col_type=cols(.default = "c"))
PCproperties <- slice(PCproperties, match(PCproperties$LATCA_ID, sortbyid))
PCproperties <- rename(PCproperties, LATCA_ID_PubChem=LATCA_ID)
CMproperties <- readr::read_tsv("./results/latca/latca_desc_CM.tsv", col_type=cols(.default = "c"))
CMproperties <- slice(CMproperties, match(CMproperties$LATCA_ID, sortbyid))
CMproperties <- rename(CMproperties, LATCA_ID_ChemmineR=LATCA_ID, MW_CM=MW)
CDKproperties <- readr::read_tsv("./results/latca/latca_desc_CDK.tsv", col_type=cols(.default = "c"))
CDKproperties <- slice(CDKproperties, match(CDKproperties$LATCA_ID, sortbyid))
CDKproperties <- rename(CDKproperties, LATCA_ID_CDK=LATCA_ID)
latca_all <- bind_cols(latca_anno, PCproperties, CMproperties, CDKproperties)
write_tsv(latca_all, "./results/latca/latca_all.tsv")
mycol_names <- as_tibble(colnames(latca_all))
mycol_names <- rename(mycol_names, Column_Names=value)
write_tsv(mycol_names, "./results/latca/latca_all_cols.tsv")


##################################
## Additional Notes and Testing ##
##################################

## (1b) Retrieve various CMP ID types from CAS ID. This approach avoids the cleanup by downloading the source compounds form other sources
## Note option (1b) was not used since (1a) provided better results. 
library(webchem); library(ChemmineR); library(ChemmineOB)
casid <- c("4474-91-3", "102767-28-2", "103060-53-3") # From Selleck FDA drug xls file here: https://bit.ly/3nskoei

## Get pubchem id for cas id 
caspubids <- get_cid(casid, from = "xref/rn")
sdfset <- pubchemCidToSDF(as.numeric(caspubids$cid))
sdfset

## Get ChEBI id for cas id
caschebiids <- get_chebiid(casid, from = "registry numbers")
caschebiids

## (2) Remove disconnected fragments including salts from FDA approvded drugs downloaded from Selleck here: https://bit.ly/3nskoei
library(webchem); library(ChemmineR); library(ChemmineOB)
convertFormatFile(from="SMI", to="SDF", fromFile="smi.txt", toFile="sdfsplit.txt", options=data.frame(names="separate", args=""))
sdfsplit <- read.SDFset("sdfsplit.txt")
## Function to choose which fragment to keep; currently it uses compound with largest heavy atom count
whichToKeep <- function(x=sdfsplit, mode="max_atoms") {
    indexDF <- data.frame(sdfid_org = sdfid(x),
                            sdfid_base = gsub("#.*", "", sdfid(x)),
                            n_atoms = vapply(atomcount(x), sum, FUN.VALUE=numeric(1)))
    indexDF <- indexDF[order(indexDF$sdfid_base, -indexDF$n_atoms),]
    if(mode=="max_atoms") {
        indexDFnodup <- indexDF[!duplicated(indexDF$sdfid_base),]
        return(as.character(indexDFnodup$sdfid_org))
    } else {
        stop("only 'mode=max_mode' implemented so far!")
    }
}
sdfid_keep <- whichToKeep(x=sdfsplit, mode="max_atoms")
sdfset <- sdfsplit[sdfid(sdfsplit) %in% sdfid_keep]
sdfset

################################
## Analysis done on 16-Feb-24 ##
################################
library(RSQLite); library(dplyr)
mydb <- dbConnect(SQLite(), "~/Projects/REU_22/REU2022-Lupe/Project/data/chembl_30/chembl_30_sqlite/chembl_30.db")
## missing_data.xlsx is corresponding tab from here: https://docs.google.com/spreadsheets/d/1u-tH373OoXSrdSDm57cQ6cz1i7HcLZXP/edit?usp=sharing&ouid=109362027298751400283&rtpof=true&sd=true
missing_data <- read.delim("./results/missing_data.xlsx", sep="\t", quote="")
## Now fill columns with corresponding data from ChEMBL
masterDF <- read.delim("./results/masterDF.tsv", sep="\t", quote="")
rownames(masterDF) <- masterDF$ChEMBL_ID
struct <- dbGetQuery(mydb, 'SELECT * FROM compound_structures')
rownames(struct) <- as.character(struct$molregno)
struct <- struct[,-c(1,2)]
masterDFstruct <- data.frame(masterDF, struct[as.character(masterDF$molregno),])
write.table(masterDFstruct, "./results/masterDFstruct.tsv", row.names=FALSE, quote=FALSE, sep="\t")
masterDFstruct <- read.delim("./results/masterDFstruct.tsv", sep="\t", quote="")
missing_complete <- masterDFstruct[masterDFstruct$ChEMBL_ID %in% unique(missing_data$ChEMBL_parent_ID),]
missing_complete <- missing_complete[unique(as.character(missing_data$ChEMBL_parent_ID)),]
write.table(missing_complete, "./results/missing_complete.tsv", row.names=FALSE, quote=FALSE, sep="\t")
missing_complete_parent <- masterDFstruct[as.character(masterDFstruct$parent_molregno) %in% unique(as.character(missing_complete$parent_molregno)),]
write.table(missing_complete_parent, "./results/missing_complete.tsv", row.names=FALSE, quote=FALSE, sep="\t")

## Append same data to Selleck data set with simple name match
library("tidyverse")
slk <- readr::read_tsv("./data/20220301-L1300-FDA-approved-Drug-Library-96-well.tsv")
masterDFstruct_nodup <- masterDFstruct[!duplicated(masterDFstruct$ChEMBL_NAME),]
slk2chembl <- masterDFstruct_nodup[toupper(masterDFstruct_nodup$ChEMBL_NAME) %in% toupper(slk$Name),]
slk <- as.data.frame(slk)
rownames(slk) <- toupper(as.character(slk$Name))
slk2chemblDF <- data.frame(slk2chembl, slk[toupper(slk2chembl$ChEMBL_NAME),])
write.table(slk2chemblDF, "./results/slk2chemblDF.tsv", row.names=FALSE, quote=FALSE, sep="\t")
slk2chemblDF_parent <- masterDFstruct[as.character(masterDFstruct$parent_molregno) %in% unique(as.character(slk2chemblDF$parent_molregno)),]
slk2chemblDF_parent <- data.frame(slk2chemblDF_parent, slk[toupper(slk2chemblDF_parent$ChEMBL_NAME),])
write.table(slk2chemblDF_parent, "./results/slk2chemblDF_parent.tsv", row.names=FALSE, quote=FALSE, sep="\t")

## Check whether Selleck to ChEMBL (src1src20.txt.gz) gives additional mappings for Selleck
system("wget https://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id1/src1src20.txt.gz -O ./data/src1src20.txt.gz && gunzip ./data/src1src20.txt.gz")
selleck2chembl <- readr::read_tsv("./data/src1src20.txt")
slk2unichem <- slk[toupper(slk$Name) %in% toupper(selleck2chembl$"To src:'20'"),]
sum(!toupper(slk2unichem$Name) %in% toupper(slk2chembl$ChEMBL_NAME)) #-> only 14 additional ones

