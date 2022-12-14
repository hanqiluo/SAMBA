
biomarker_dataset_template <- "~/git/SAMBA2/code/biomarker_dataset_template_2billionV2.csv"
biomarker_cutoff_template <- "~/git/SAMBA2/code/biomarker_cutoff_template.csv"
no_biodata <- 1:2
save_directory <- "~/git/SAMBA/final_data"
results_name <- "kenya_nhanes_ethiopia"

temp_template2 <- 
  load_templates(biomarker_dataset_template = biomarker_dataset_template, 
                 biomarker_cutoff_template = biomarker_cutoff_template, 
                 save_directory = save_directory)


