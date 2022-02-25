#
# Micronutrient biomarker analysis
# SAMBA R package
#
rm(list = ls())
gc()
source("~/git/SAMBA/package/samba_v1.6_df.R")


# execute the SAMBA R package
samba(biomarker_dataset_template = "~/git/SAMBA/code/biomarker_dataset_template_2billionV2.csv",
      biomarker_cutoff_template = "~/git/SAMBA/code/biomarker_cutoff_template.csv",
      no_biodata = 1:2,
      save_directory = "~/git/SAMBA/final_data",
      results_name = "kenya_nhanes_ethiopia")
