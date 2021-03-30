#
# Micronutrient biomarker analysis
# SAMBA R package
#
rm(list = ls())
gc()
source("~/OneDrive - Emory University/samba/samba_osf/package/samba_v1.2.R")


# execute the SAMBA R package
samba(biomarker_dataset_template = "~/git/2billion/demo/bio_analysis/code/biomarker_dataset_template_paper.csv",
      biomarker_cutoff_template = "~/git/2billion/demo/bio_analysis/code/biomarker_cutoff_template.csv",
      no_biodata = 1:2,
      save_directory = "~/git/2billion/demo/bio_analysis/final_data/test",
      results_name = "nhanes_kenya")


#######################################
# testing package
########################################
biomarker_dataset_template = "~/git/2billion/demo/bio_analysis/code/biomarker_dataset_template_paper.csv"
biomarker_cutoff_template = "~/git/2billion/demo/bio_analysis/code/biomarker_cutoff_template.csv"
no_biodata = 2
i = 1
save_directory = "~/git/2billion/demo/bio_analysis/final_data/test"
results_name = "nhanes_kenya"


select_biomarker = "serum_ferritin"
input_data <- final_biodata
template   <- temp_template2
data <- input_data

#################################################################################
# Comparison
#################################################################################

#
# Kenya data
#
kenya_psc <- read_sas   ("~/git/2billion/demo/bio_analysis/raw_data/brinda_ke2007c_20160204.sas7bdat")
kenya_samba <- read.csv("~/git/2billion/demo/bio_analysis/final_data/test/biodata1_clean_brinda_ke2007c_20160204.csv",
                         sep = ",")

cbind(kenya_samba$serum_ferritin, kenya_psc$sf)

identical(as.numeric(kenya_samba$serum_ferritin), as.numeric(kenya_psc$sf), num.eq = T, single.NA = F)

# Apply Yaw's adjustment
kenya_psc2 <- brindarc(res_data = kenya_psc)

kenya_psc2$brcadj_ferritin
kenya_psc2$brcadj_rbp
kenya_psc2$brcadj_tfr

# compare serum ferritin
cbind(as.numeric(kenya_samba$serum_ferritin_adj), as.numeric(kenya_psc2$brcadj_ferritin))
#identical(round(kenya_samba$serum_ferritin_adj, 3), round(kenya_psc2$brcadj_ferritin, 3), num.eq = T, single.NA = F)
all.equal(as.numeric(kenya_samba$serum_ferritin_adj), as.numeric(kenya_psc2$brcadj_ferritin))

# compare rbp 
cbind(as.numeric(kenya_samba$retinol_binding_protein_adj), as.numeric(kenya_psc2$brcadj_rbp))
#identical(round(kenya_samba$serum_ferritin_adj, 3), round(kenya_psc2$brcadj_ferritin, 3), num.eq = T, single.NA = F)
all.equal(as.numeric(kenya_samba$retinol_binding_protein_adj), as.numeric(kenya_psc2$brcadj_rbp))

# compare transferrin receptor
cbind(as.numeric(kenya_samba$transferrin_receptor_adj), as.numeric(kenya_psc2$brcadj_tfr))
#identical(round(kenya_samba$serum_ferritin_adj, 3), round(kenya_psc2$brcadj_ferritin, 3), num.eq = T, single.NA = F)
all.equal(as.numeric(kenya_samba$transferrin_receptor_adj), as.numeric(kenya_psc2$brcadj_tfr))

#
# NHANES
#
library(Hmisc)
# NHANES PSC
nhanes_psc <- read_sas("~/git/2billion/demo/bio_analysis/raw_data/brinda_us2006c_20160204.sas7bdat")
contents(nhanes_psc)


#
# NHANES WRA
#
library(haven)

connector <- ifelse(Sys.info()['sysname'] == "Darwin","/","\\")
a <- "\\"

nhanes_wra <- read_sas(paste0("~/git/2billion/demo/bio_analysis/raw_data", connector, "brinda_us2006w_20151130.sas7bdat"))
test <- try( read_sas(paste0("~/git/2billion/demo/bio_analysis/raw_data",  connector,  "brinda_us2006w_20151130.sas7bdat")), silent = F)


nhanes_wra <- read_sas(paste0("/git/2billion/demo/bio_analysis/raw_data", "//", "brinda_us2006w_20151130.sas7bdat"))




nhanes_wra <- read_sas(paste0("/Users/validinternational/git/2billion/demo/bio_analysis/raw_data", "\\", "brinda_us2006w_20151130.sas7bdat"))
nhanes_wra$sex <- 2

# save data
write_sas(nhanes_wra, "~/git/2billion/demo/bio_analysis/raw_data/brinda_us2006w_20210225.sas7bdat")

# read the samba
nhanes_wra_samba <- read.csv("~/git/2billion/demo/bio_analysis/final_data/test/biodata2_clean_brinda_us2006w_20210225.csv",
                         sep = ",")



cbind(nhanes_wra_samba$transferrin_receptor, nhanes_wra_samba$transferrin_receptor_adj)
all.equal(nhanes_wra_samba$transferrin_receptor, nhanes_wra_samba$transferrin_receptor_adj)

nhanes_wra$agp <- rep(0, nrow(nhanes_wra))
nhanes_wra$rbp <- nhanes_wra$sr

nhanes_wra2 <- brindarw(res_data = nhanes_wra)
#res_data = nhanes_wra

# compare the results
# compare serum ferritin
cbind(as.numeric(nhanes_wra_samba$serum_ferritin_adj), as.numeric(nhanes_wra2$brcadj_ferritin))
#identical(round(kenya_samba$serum_ferritin_adj, 3), round(nhanes_wra2$brcadj_ferritin, 3), num.eq = T, single.NA = F)
all.equal(as.numeric(nhanes_wra_samba$serum_ferritin_adj), as.numeric(nhanes_wra2$brcadj_ferritin))

# compare rbp 
cbind(as.numeric(nhanes_wra_samba$retinol_binding_protein_adj), as.numeric(nhanes_wra2$brcadj_rbp))
#identical(round(kenya_samba$serum_ferritin_adj, 3), round(nhanes_wra2$brcadj_ferritin, 3), num.eq = T, single.NA = F)
all.equal(as.numeric(nhanes_wra_samba$retinol_binding_protein_adj), as.numeric(nhanes_wra2$brcadj_rbp))
# serum retinol not adjusted

# compare transferrin receptor
cbind(as.numeric(kenya_samba$transferrin_receptor_adj), as.numeric(nhanes_wra2$brcadj_tfr))
#identical(round(kenya_samba$serum_ferritin_adj, 3), round(nhanes_wra2$brcadj_ferritin, 3), num.eq = T, single.NA = F)
all.equal(as.numeric(kenya_samba$transferrin_receptor_adj), as.numeric(nhanes_wra2$brcadj_tfr))




