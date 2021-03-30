#
# Micronutrient biomarker analysis
# SAMBA R package
#
rm(list = ls())
gc()
source("~/git/2billion/SAMBA/package/samba_v1.6_df.R")


# execute the SAMBA R package
samba(biomarker_dataset_template = "~/git/2billion/SAMBA/code/biomarker_dataset_template_2billionV2.csv",
      biomarker_cutoff_template = "~/git/2billion/SAMBA/code/biomarker_cutoff_template.csv",
      no_biodata = 1:3,
      save_directory = "~/git/2billion/SAMBA/final_data",
      results_name = "kenya_nhanes_ethiopia")
