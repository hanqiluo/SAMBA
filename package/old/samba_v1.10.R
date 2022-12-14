# Main SAMBA function
samba <- function(biomarker_dataset_template,
                  no_biodata,
                  save_directory,
                  biomarker_cutoff_template,
                  results_name) {
 
  #
  # Install libraries
  # if library does not exist, install libraries
  #
  packages_list <- c("plyr", "dplyr", "tidyr", "data.table", "stringr", "purrr", "survey", "haven",  "readxl", "berryFunctions", "tibble")
  new_packages <- packages_list[!(packages_list %in% installed.packages()[,"Package"])]
  if(length(new_packages)) install.packages(new_packages, dependencies = TRUE)
  sapply(packages_list, require, character.only = TRUE)
  rm(new_packages, packages_list)
  
  #
  # setting
  #
  options(survey.lonely.psu="adjust")
  options(stringsAsFactors = FALSE)
  connector <- ifelse(Sys.info()['sysname'] == "Darwin","/","\\")
  
  #
  # Load the datasets
  # Provide error messages when datasets do not exist 
  #
  two_templates <- load_templates(biomarker_dataset_template, biomarker_cutoff_template, save_directory)
  biomarker_template <- two_templates$temp_template 
  reference <- two_templates$reference
  
 
  # results template before the loop
  results_all <- NULL
  
  #
  # LOOP starts 
  #
  
  # loop - TBD
  for (i in no_biodata) {
    #
    # Check if each parameter is correct -----------------
    #
    # transpose template to an easy-to-read form
      
    current_biodata_template <- check_format(biomarker_template, i)
    
    biodata <- import_biodata(current_biodata_template, connector)
    
    #
    # check variables - Reference ---------
    #
    
    variables_checked <- check_variables(reference, current_biodata_template)
    reference <- variables_checked$reference  
    current_biodata_template <- variables_checked$current_biodata_template 
    
    #
    # function: Create a new dataset with standarized variable name
    #
    
    biodata_transformed <- transform_biodata(current_biodata_template, biodata) 
    current_biodata_template <- biodata_transformed$current_biodata_template
    final_biodata <- biodata_transformed$final_biodata
      
    #
    # function: Change unit to standard units
    #
    converted_datasets <- convert_units(current_biodata_template, final_biodata)
    current_biodata_template <- converted_datasets$current_biodata_template 
    final_biodata <- converted_datasets$final_biodata
    
    #
    # Acutal calculate deficiencies ---------
    #
    final_biodata <- calculate_deficiencies(input_data = final_biodata, 
                                            reference = reference, 
                                            template = current_biodata_template)
    
    
    # function: save clean data ---------
    # 
   
    
    # organize the variable list of the results_all
    data_saved <- save_data(no_biodata, final_biodata, save_directory, i, current_biodata_template, results_all, results, results_name)
    results_all <- data_saved$results_all
  }  
  
  fwrite(results_all, file = paste0(save_directory, "/biodata_results_", results_name, ".csv"), row.names = F)
  print(paste0("-------------Analysis of all biomarker dataset(s) is complete-------------"))
}
# SAMBA main function ends --------------------------

#
# Load templates ----------
#
load_templates <- 
    function(biomarker_dataset_template, biomarker_cutoff_template, save_directory){
      # load the template  
      if (file.exists(biomarker_dataset_template)) {
        temp_template <-  fread(file = biomarker_dataset_template, header = TRUE, 
                                sep = ",")  %>% as.data.frame()
        temp_template$information <- str_replace(temp_template$information, " ", "")
        rownames(temp_template) <- temp_template$information
        temp_template$information <- NULL
      } else {
        stop ("Bio template does not exist! Please select the correct file")
      }
      
      # load reference table
      if (file.exists(biomarker_cutoff_template)) {
        reference <-  fread(file = biomarker_cutoff_template, header = TRUE, 
                            sep = ",")  %>% as.data.frame()
      } else {
        stop ("Reference template does not exist! Please select the correct file")
      }
      
      # check if the save directory exist
      if (!dir.exists(save_directory)) {
        stop ("Save directory does not exist! Please check the save directory")
      } 
      return(list(reference = reference, temp_template = temp_template))
}

#
# check formats -----------
#
check_format <- 
    function(temp_template, no_biodata){
      current_biodata_template <- 
        temp_template[no_biodata] %>%
        t() %>%
        as.data.frame() %>%
        mutate_all(~na_if(., ""))
  
      #
      # add data - from different forms
      #
      current_biodata_template$file_format <- tolower(sapply(strsplit(current_biodata_template$data, "\\."), tail, 1))
      
      current_biodata_template$file_name <- substr(current_biodata_template$data, 1, 
                                         nchar(current_biodata_template$data) - nchar(current_biodata_template$file_format) - 1)
      
      if(!(current_biodata_template$file_format %in% c("csv", "dta", "sas7bdat", "sav", "xls", "xlsx", "xpt"))){
        stop(paste0(current_biodata_template$data, " in ", current_biodata_template$file_format, " format. 
    This is not a supported format. Bio_analysis function reads the CSV, excel, XPT, STATA, and SPSS data files"))
      }
      
      return(current_biodata_template)
}

#
# Import biodata --------------
#
import_biodata <- 
    function(current_biodata_template, connector){
      # check if the biomarker dataset exists
      if (file.exists(paste0(current_biodata_template$directory, connector, current_biodata_template$data))) {
        print(paste0("-------------------------------------------------"))
        print(paste0("-------------Loading ", current_biodata_template$directory, connector, current_biodata_template$data, "-------------"))
        
        # read csv
        if(current_biodata_template$file_format == "csv") {
          biodata <- as.data.table((fread(file = paste0(current_biodata_template$directory, connector, current_biodata_template$data), 
                                          header = T, sep = ",")))
        }
        
        # read SAS file
        if(current_biodata_template$file_format == "sas7bdat"){
          biodata <- as.data.table(read_sas(paste0(current_biodata_template$directory, connector, current_biodata_template$data)))
        }
        
        # read SPSS
        if(current_biodata_template$file_format == "sav"){
          biodata <- as.data.table(read_sav(paste0(current_biodata_template$directory, connector, current_biodata_template$data)))
        }
        
        # read STATA
        if(current_biodata_template$file_format == "dta"){
          biodata <- as.data.table(read_dta(paste0(current_biodata_template$directory, connector, current_biodata_template$data)))
        }
        
        # read excel
        if(current_biodata_template$file_format %in% c("xls", "xlsx")){
          biodata <- as.data.table(read_excel(paste0(current_biodata_template$directory, connector, current_biodata_template$data)))
        }
        
        # read xpt
        if(current_biodata_template$file_format == "xpt"){
          biodata <- as.data.table(read_xpt(paste0(current_biodata_template$directory, connector, current_biodata_template$data)))
        }
        
      } else {
        stop(paste0(current_biodata_template$data, " does not exist under ", current_biodata_template$directory))
      }
      
      # add ID if there is no unique ID
      if(is.na(current_biodata_template$id)){
        biodata$id = seq(nrow(biodata))
        current_biodata_template$id = "id"
      }
      
      return(biodata)
}

#
# check variables ------------------
#
check_variables <- 
    function(reference, current_biodata_template){
      #
      # remove the space in variable name
      #
      char_values <- c("nutrients", "biomarker", "relationship", "cutoff",
                       "unit", "population",
                       "lower_age", "upper_age", "note")
  
      sapply(char_values, function(x) {
        reference[, x] <<- str_replace_all(reference[, x], " ", "")
      }) %>%
        invisible()
      rm(char_values)
  
  #
  # make sure all the numeric values were recorded into numeric values
  #
  numeric_var <- c("cutoff", "lower_age", "upper_age")
  
  sapply(numeric_var, function(x) {
    reference[, x] <<- as.numeric(reference[, x])}) %>%
    invisible()
  rm(numeric_var)
  
  # notification
  print(paste0("-------------Analyzing ", current_biodata_template$data, "-------------"))
  
  #
  # check if necessary inputs are specified
  #
  necessary_values <- c("country", "survey_year", "complex_survey_design_yes_no", 
                        "study_design", "sex_var", "age", "age_unit")
  sapply(necessary_values, function(x) {if(is.na(current_biodata_template[, x])) {
    stop(print(paste0("You must fill a value for ", x, " in the template")))}}) %>%
    invisible()
  rm(necessary_values)
  
  return(list(reference = reference, current_biodata_template = current_biodata_template))
}

#
# transform biodata ---------
#
transform_biodata <- function(current_biodata_template, biodata){
  
  biomarker_var_survey_weight <- grep("_survey_weight.*", names(current_biodata_template), value=TRUE)
  
  biomarker_unit              <- grep("_unit.*", names(current_biodata_template), value=TRUE)[-1] # exclude age variable
  
  biomarker_var               <- substr(biomarker_unit, 1, regexpr("\\_[^\\_]*$", biomarker_unit) - 1)
  
  # remove the space in variable name
  char_values <- c("id", "strata", "PSU", "sex_var",
                   "age", "age_unit",
                   "male_value", "pregnancy_value", "lactating_value",
                   "PSU_ID_nested_in_strataID_yes_no",
                   "complex_survey_design_yes_no", 
                   "sex_var",
                   "pregnancy_var", 
                   "lactating_var",
                   "subgroup",
                   biomarker_unit,
                   biomarker_var_survey_weight,
                   biomarker_var
  )
  
  sapply(char_values, function(x) {
    current_biodata_template[, x] <<- str_replace_all(current_biodata_template[, x], " ", "")
  }) %>%
    invisible()
  rm(char_values)
  
  # make sure all the numeric values were recorded into numeric values
  #numeric_var <- c("male_value", "pregnancy_value", "lactating_value")
  
  #sapply(numeric_var, function(x) {
  #    if(exists(numeric_var, current_biodata_template)) {
  #        current_biodata_template[, x] <<- as.numeric(current_biodata_template[, x])}}) %>%
  #    invisible()
  #rm(numeric_var)
  
  # if a biomarker variable is specified, its unit must be specified; also the age variable
  mapply(function(x, y) {
    if(!is.na(current_biodata_template[, x]) & is.na(current_biodata_template[, y])){
      stop(print(paste0("When you specify ", x, ", you must specify ",  y)))
    }
  },  c("age", biomarker_var), c("age_unit", biomarker_unit)) %>%
    invisible()
  
  # lower case the values of unit variables, and the values of yes/no variables
  sapply(c("age_unit", biomarker_unit, "PSU_ID_nested_in_strataID_yes_no",
           "complex_survey_design_yes_no"), function(x) {
             current_biodata_template[, x] <<- tolower(current_biodata_template[, x])
           }) %>%
    invisible()
  
  # Check if the yes_no questions are with valid input
  if(!current_biodata_template$complex_survey_design_yes_no %in% c("yes", "n", "no", "n")){
    stop(print(paste0("Invalid input for complex_survey_design_yes_no")))
  }
  
  if(!current_biodata_template$PSU_ID_nested_in_strataID_yes_no %in% c("yes", "n", "no", "n", NA)){
    stop(print(paste0("Invalid input for complex_survey_design_yes_no")))
  }
  
  # check if the age variable has the correct input
  if(!(current_biodata_template$age_unit %in% c("year", "month"))){
    stop(print("Age unit can only be 'year' or 'month', please correct your input in the template"))
  }
  
  # if there is complex survey design, there must be weight for 
  if(current_biodata_template$complex_survey_design_yes_no %in% c("Yes", "Y")){
    mapply(function(x, y) {
      if(!is.na(current_biodata_template[, x]) & is.na(current_biodata_template[, y])){
        stop(print(paste0("When you specify ", x, ", you must specify ",  y)))
      }
    },  biomarker_var, biomarker_survey_weight) %>%
      invisible()
  }
  
  # If the pregnancy/lactating variables exist, there must be pregnancy/lactating variable value
  preg_lac_var <- c("pregnancy_var", "lactating_var")
  preg_lac_value <- c("pregnancy_value", "lactating_value")
  mapply(function(x, y) {
    if(!is.na(current_biodata_template[, x]) & is.na(current_biodata_template[, y])){
      stop(print(paste0("When you specify ", x, ", you must specify ",  y)))
    }
  },  preg_lac_var, preg_lac_value) %>%
    invisible()
  
  rm(preg_lac_value, preg_lac_var)
  
  #    
  # check every variables exist
  #
  covariates <- as.vector(str_split(current_biodata_template$covariates, " ", simplify =  T))
  
  variable_list <- c("id", "strata", "PSU", "sex_var", "age", "pregnancy_var", 
                     "lactating_var", "subgroup",
                     #biomarker_unit,
                     biomarker_var_survey_weight,
                     biomarker_var)
  
  #
  # check each variable and create new dataset
  #
  final_biodata <- NULL
  for (var in variable_list) {
    if (!is.na(current_biodata_template[, var])) {
      if (exists(current_biodata_template[, var], biodata)) {
        # we need to select variables one by one, due to two weight variable can be the same variable
        final_biodata <- cbind(final_biodata, biodata[, current_biodata_template[, var], with = F])
        setnames(final_biodata, current_biodata_template[, var], var)
      } else {
        stop (paste0(current_biodata_template[, var], " does not exist in ", current_biodata_template$data))
      }
    }
  }
  
  for(var in biomarker_unit) {
    if(!is.na(current_biodata_template[, var])) {
      final_biodata <- final_biodata %>%
        mutate(unit_var = current_biodata_template[, var])
      setnames(final_biodata, "unit_var", var)
    }
  }
  
  # check covariates exist, if so, merge back to the dataset
  if(!is.na(current_biodata_template$covariates)) {
    for (var in covariates) {
      if (exists(var, biodata)) {
        # we need to select variables one by one, due to two weight variable can be the same variable
        final_biodata <- cbind(final_biodata, biodata[, var, with = F])
      } else {
        stop (paste0(var, "does not exist in ", current_biodata_template$data))
      }
    }
  }
  rm(covariates, var)
  
  # generate variables for male, female, pregnant
  if(exists("pregnancy_var", final_biodata)) {
    final_biodata <- 
      final_biodata %>%
      mutate(male = ifelse((sex_var == as.numeric(current_biodata_template$male_value)), 1, 0),
             np_female = ifelse((sex_var != as.numeric(current_biodata_template$male_value) & (pregnancy_var != as.numeric(current_biodata_template$pregnancy_value) |is.na(pregnancy_var))), 1, 0),
             preg = case_when((pregnancy_var == as.numeric(current_biodata_template$pregnancy_value)) ~ 1,
                              pregnancy_var != as.numeric(current_biodata_template$pregnancy_value) | is.na(pregnancy_var) ~ 0))
  } else {
    final_biodata <- 
      final_biodata %>%
      mutate(male = ifelse((sex_var == as.numeric(current_biodata_template$male_value)), 1, 0),
             np_female = ifelse((sex_var != as.numeric(current_biodata_template$male_value)), 1, 0),
             preg = 0)  
  }
  
  # convert age into year
  if(current_biodata_template$age_unit == "year"){
    final_biodata <- 
      final_biodata %>%
      mutate(age = round(age, 1))
  } else {
    final_biodata <- 
      final_biodata %>%
      mutate(age = round(age/12, 1))    
  }
  
  
  # change PSU/STRATA if only one level
  if(length(unique(final_biodata$PSU)) == 1) {
    current_biodata_template$PSU = NA
  } 
  
  if(length(unique(final_biodata$strata)) == 1) {
    current_biodata_template$strata = NA
  } 
  
  return(list(final_biodata = as.tibble(final_biodata), current_biodata_template = current_biodata_template))
}
#
# transform biodata end ---------
#


#
# convert units ------------------
#
convert_units <- 
    function(current_biodata_template, final_biodata){
      
      #
      # serum b-12
      #
      if(!is.na(current_biodata_template$serum_b12)){
        if(!(current_biodata_template$serum_b12_unit %in% c("pg/ml", "pmol/l"))){
          stop(print("Serum/plasma B12 unit can only be 'pg/ml' or 'pmol/l'"))
        } 
        
        if(current_biodata_template$serum_b12_unit == "pg/ml") {
          print("convert serum b12 unit from pg/ml to pmol/l")
          final_biodata$serum_b12 = final_biodata$serum_b12 * 0.738
          current_biodata_template$serum_b12_unit = "pmol/l"
        }
      }
      
      #
      # serum folate
      #
      if(!is.na(current_biodata_template$serum_folate)) {
        
        if(!(current_biodata_template$serum_folate_unit %in% c("nmol/l", "ng/ml"))){
          stop(print("Serum folate unit can only be 'nmol/l' or 'ng/ml'"))
        } 
        
        if(current_biodata_template$serum_folate_unit == "ng/ml") {
          print("convert serum folate unit from ng/ml to nmol/l")
          final_biodata$serum_folate = final_biodata$serum_folate * 2.2655 
          current_biodata_template$serum_folate_unit = "nmol/l"
        }
      }
      
      #
      # RBC folate
      #
      if(!is.na(current_biodata_template$rbc_folate)){
        if(!(current_biodata_template$rbc_folate_unit %in% c("nmol/l", "ng/ml"))){
          stop(print("RBC folate unit can only be 'nmol/l' or 'ng/ml'"))
        } 
        
        if(current_biodata_template$rbc_folate_unit == "ng/ml") {
          print("convert RBC folate unit from ng/ml to nmol/l")
          final_biodata$rbc_folate = final_biodata$rbc_folate * 2.2655 
          current_biodata_template$rbc_folate_unit = "nmol/l"
        }
      }
      
      #
      # serum retinol
      #
      if(!is.na(current_biodata_template$serum_retinol)) {
        
        if(!(current_biodata_template$serum_retinol_unit %in% c("umol/l", "ug/dl"))){
          stop(print("Serum retinol unit can only be 'umol/l' or 'ug/dl'"))
        } 
        
        if(current_biodata_template$serum_retinol_unit == "ug/dl") {
          print("convert serum retinol unit from ug/dl to umol/l")
          final_biodata$serum_retinol = final_biodata$serum_retinol * 0.03491  
          current_biodata_template$serum_retinol_unit = "nmol/l"
        }
      }
      
      #
      # Retinol binding protein
      #
      if(!is.na(current_biodata_template$retinol_binding_protein)) {
        
        if(!(current_biodata_template$retinol_binding_protein_unit %in% c("umol/l", "ug/dl"))){
          stop(print("Retinol binding protein unit can only be 'umol/l' or 'ug/dl'"))
        } 
        
        if(current_biodata_template$retinol_binding_protein == "ug/dl") {
          print("convert serum retinol unit from ug/dl to umol/l")
          final_biodata$retinol_binding_protein = final_biodata$serum_retinol * 0.03491  
          current_biodata_template$serum_retinol_unit = "nmol/l"
        }
      }
      
      #
      # serum zinc
      #
      if(!is.na(current_biodata_template$serum_zinc)) {    
        # serum zinc
        if(!(current_biodata_template$serum_zinc_unit %in% c("umol/l", "ug/dl", "ug/l"))){
          stop(print("Serum zinc unit can only be 'umol/l' or 'ug/ml'"))
        } 
        
        if(current_biodata_template$serum_zinc_unit == "ug/dl") {
          print("convert serum zinc unit from ug/dl to umol/l")
          final_biodata$serum_zinc = final_biodata$serum_zinc * 0.15291
          
          current_biodata_template$serum_zinc_unit = "nmol/l"
        }
        
        if(current_biodata_template$serum_zinc_unit == "ug/l") {
          print("convert serum zinc unit from ug/l to umol/l")
          final_biodata$serum_zinc = final_biodata$serum_zinc * 0.015291
          
          current_biodata_template$serum_zinc_unit = "nmol/l"
        }
      }
      
      #
      # serum ferritin
      #
      if(!is.na(current_biodata_template$serum_ferritin)) {
        if(!(current_biodata_template$serum_ferritin_unit %in% c("ng/ml", "ug/l"))){
          stop(print("Serum ferritin unit can only be 'ng/ml' or 'ug/l'"))
        } 
        
        if(current_biodata_template$serum_ferritin_unit == "ng/ml") {
          print("convert serum ferritin unit from ng/ml to ug/l")
          final_biodata$serum_ferritin = final_biodata$serum_ferritin * 1 
          current_biodata_template$serum_ferritin_unit = "ug/l"
        }
      }
      
      # Soluable Transferrin receptor
      if(!is.na(current_biodata_template$transferrin_receptor)) {
        if(!(current_biodata_template$transferrin_receptor_unit %in% c("mg/l", "mg/dl"))){
          stop(print("Transferrin receptor unit can only be 'mg/l' or 'mg/dl'"))
        } 
        
        if(current_biodata_template$transferrin_receptor_unit == "mg/dl") {
          print("convert transferrin receptor unit from mg/dl to mg/l")
          final_biodata$transferrin_receptor = final_biodata$transferrin_receptor * 10 
          current_biodata_template$transferrin_receptor_unit = "mg/l"
        }
      }
      
      # #
      # # Iodine
      # #
      # if(!is.na(current_biodata_template$urine_iodine)) {
      #   if(!(current_biodata_template$urine_iodine_unit %in% c("ug/l", "ng/ml"))){
      #     stop(print("Urine iodine unit can only be 'ug/l' or 'ng/ml'"))
      #   } 
      #   
      #   if(current_biodata_template$urine_iodine_unit == "ng/ml") {
      #     print("convert urine iodine unit from ng/ml to ug/l")
      #     final_biodata$urine_iodine = final_biodata$urine_iodine * 1 
      #     current_biodata_template$urine_iodine_unit = "ug/l"
      #   }
      # }
      
      # CRP
      if(!is.na(current_biodata_template$crp)) {
        if(!(current_biodata_template$crp_unit %in% c("mg/l"))){
          stop(print("CRP unit can only be 'mg/l'"))
        } 
      }
      
      # AGP
      if(!is.na(current_biodata_template$agp)) {
        if(!(current_biodata_template$agp_unit %in% c("g/l"))){
          stop(print("AGP unit can only be 'g/l'"))
        } 
      }
      
      # vitamin D
      if(!is.na(current_biodata_template$serum_vitaminD)) {
        if(!(current_biodata_template$serum_vitaminD_unit %in% c("ng/ml", "nmol/l"))){
          stop(print("Serum Vitamin D unit can only be 'ng/ml' or 'nmol/l'"))
        } 
        
        if(current_biodata_template$serum_vitaminD_unit == "ng/ml") {
          print("convert serum vitamin D unit from ng/ml to nmol/l")
          final_biodata$serum_vitaminD = final_biodata$urine_iodine * 2.496 
          current_biodata_template$serum_vitaminD_unit = "nmol/l"
        }
      }
      
      return(list(current_biodata_template = current_biodata_template, final_biodata = final_biodata)) 
}

#
# survey design ---------
#
survey_design <- function(input_data, template) {
  # complex survey design
  # strata and PSU both exist 
  if(!is.na(template$strata) & !is.na(template$PSU)) {
    if (template$PSU_ID_nested_in_strataID_yes_no == "yes") {
      mydesign <- 
        svydesign(
          id = ~PSU,
          data = input_data |>
                     filter(!is.na(biomarker) & !is.na(wt_temp)),
          weight = ~wt_temp,
          strata = ~strata,
          nest = TRUE)
    }
    
    if (template$PSU_ID_nested_in_strataID_yes_no == "no") {
      mydesign <- 
        svydesign(
          id = ~PSU,
          data = input_data |>
                    filter(!is.na(biomarker) & !is.na(wt_temp)),
          weight = ~wt_temp,
          strata = ~strata,
          nest = FALSE)
    }
  }
  
  # only strata, no cluster
  if(!is.na(template$strata) & is.na(template$PSU)) { 
    mydesign <- 
      svydesign(
        id = ~1,
        data = input_data |>
                 filter(!is.na(biomarker) & !is.na(wt_temp)),
        weight = ~wt_temp,
        strata = ~strata
      )
  }
  
  # no strata, only cluster
  if(is.na(template$strata) & !is.na(template$PSU)) { 
    mydesign <- 
      svydesign(
        id = ~PSU,
        data = input_data |>
                    filter(!is.na(biomarker) & !is.na(wt_temp)),
        weight = ~wt_temp
      )
  }
  
  # simple site  or (no strata, no cluster, only survey weight) 
  if(is.na(template$strata) & is.na(template$PSU)) { 
    mydesign <- 
      svydesign(
        id = ~1,
        data = input_data |>
                    filter(!is.na(biomarker) & !is.na(wt_temp)),
        weight = ~wt_temp
      )
  }
  # survey design END
  return(mydesign)
}

#
# calculate summary statistics -----------
#
parameters_analysis <- function(input_data, select_biomarker, template){
  #print("parameter analysis")
  
  #print(input_data)
  
  #input_data <- data 
  
  # prepare survey setting
  mydesign <- survey_design(input_data = input_data, template)
  
  # prepare output
  output_temp <-
    cbind(
      # def
      svymean(~def,  mydesign, na.rm = T, deff = T) %>%
        data.frame() %>%
        mutate(
              #deficiency_risk_df = format(round(ifelse(deff == Inf, NA, deff), digit = 1), nsmall = 1),
                deficiency_risk_se = format(round(def * 100,  digit = 2), nsmall = 2),
                deficiency_risk =    format(round(mean * 100, digit = 2), nsmall = 2)) %>%
        dplyr::select(-mean, -def, -deff),
      
      # mean
      svymean(~biomarker,  mydesign, na.rm = T, deff = T) %>%
        data.frame() %>%
        mutate( #mean_df = format(round(ifelse(deff == Inf, NA, deff), digit = 1), nsmall = 1),
                mean_se = format(round(biomarker, digit = 2), nsmall = 2),
                mean =    format(round(mean,      digit = 2), nsmall = 2)) %>%
        dplyr::select(-biomarker, -deff),
      
      # geometric mean
      svymean(~log_biomarker,  mydesign, na.rm = T) %>%
        data.frame() %>%
        mutate( geomean_lci = format(round(exp(mean - 1.96 * log_biomarker), digit = 2), nsmall = 2),
                geomean_uci = format(round(exp(mean + 1.96 * log_biomarker), digit = 2), nsmall = 2),
                geomean =    format(round(exp(mean), digit = 2), nsmall = 2)) %>%
        dplyr::select(-log_biomarker, -mean),
      
      # percentiles 
      oldsvyquantile(~biomarker,  mydesign, quantiles = c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = T) %>%
        as.data.frame() %>%
        mutate(p10 = format(round(`0.1`,  digits = 2), nsmall = 2), 
               p25 = format(round(`0.25`, digits = 2), nsmall = 2), 
               p50 = format(round(`0.5`,  digits = 2), nsmall = 2),
               p75 = format(round(`0.75`, digits = 2), nsmall = 2),
               p90 = format(round(`0.9`,  digits = 2), nsmall = 2)) %>%
        dplyr::select(-`0.1`, -`0.25`, -`0.5`, -`0.75`, -`0.9`),
      
      oldsvyquantile(~biomarker,  mydesign, quantiles = c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = T, ci = T) %>%
        SE() %>%
        t() %>%
        as.data.frame() %>%
        mutate(p10_se = format(round(`0.1`,  digits = 2), nsmall = 2),
               p25_se = format(round(`0.25`, digits = 2), nsmall = 2), 
               p50_se = format(round(`0.5`,  digits = 2), nsmall = 2), 
               p75_se = format(round(`0.75`, digits = 2), nsmall = 2),
               p90_se = format(round(`0.9`, digits = 2), nsmall = 2)) %>%
        dplyr::select(-`0.1`, -`0.25`, -`0.5`, -`0.75`, -`0.9`)
    ) %>%
    mutate(nutrient = select_biomarker,
           subgroup = "overall",
           subgroup_name = "overall",
           n = sum(!is.na(input_data$biomarker))) %>%
    dplyr::select(nutrient, subgroup, subgroup_name, n, deficiency_risk, deficiency_risk_se, mean, mean_se, geomean, geomean_lci, geomean_uci, p10, p10_se, p25, p25_se, p50, p50_se, p75, p75_se, p90, p90_se)
  
  # subgroup analysis
  if(exists("subgroup", input_data)){
    
    # sample sizes
    sample_size <- 
      table(input_data$subgroup[!is.na(input_data$biomarker)]) %>%
      as.data.frame() %>%
      mutate(subgroup = Var1,
             n = Freq) %>%
      dplyr::select(subgroup, n)
    
    # prevalence of deficiency_risk
    def <- 
      svyby(~def, ~subgroup, mydesign, svymean, na.rm = T, drop.empty.groups = FALSE, deff = T) %>% 
      mutate(
        #deficiency_risk_df = format(round(ifelse(DEff.def == Inf, NA, DEff.def), digits = 1), nsmall = 1),
        deficiency_risk =    format(round(def * 100, digits = 2), nsmall = 2), 
        deficiency_risk_se = format(round(se * 100, digits = 2),  nsmall = 2)) %>%
      dplyr::select(-se, -def, -DEff.def)
    
    # mean
    mean <- 
      svyby(~biomarker, ~subgroup, mydesign, svymean,  na.rm = T, drop.empty.groups = FALSE, deff = T) %>%
      mutate(
            #mean_df = format(round(ifelse(DEff.biomarker == Inf, NA, DEff.biomarker), digits = 1), nmall = 1),
             mean =    format(round(biomarker, digits = 2), nsmall = 2), 
             mean_se = format(round(se, digits = 2),        nsmall = 2)) %>%
      dplyr::select(-biomarker, -se, -DEff.biomarker)
    
    # geometric mean
    geo_mean <-    
      svyby(~log_biomarker, ~subgroup, mydesign, svymean,  na.rm = T, drop.empty.groups = FALSE) %>%
      mutate( geomean =    format(round(exp(log_biomarker), digit = 2), nsmall = 2),
              geomean_lci = format(round(exp(log_biomarker - 1.96 * se), digit = 2), nsmall = 2),
              geomean_uci = format(round(exp(log_biomarker + 1.96 * se), digit = 2), nsmall = 2)) %>%
      dplyr::select(-log_biomarker, -se)
    
    # percentiles
    percentiles <- 
      svyby(~biomarker, ~subgroup, mydesign, oldsvyquantile, quantiles =  c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = T, ci = TRUE,
            drop.empty.groups = FALSE) %>%
      mutate(p10 =    format(round(`0.1`,  digits = 2), nsmall = 2), 
             p10_se = format(round(se.0.1, digits = 2), nsmall = 2),
             p25 =    format(round(`0.25`,  digits = 2), nsmall = 2), 
             p25_se = format(round(se.0.25, digits = 2), nsmall = 2),
             p50 =    format(round(`0.5`,   digits = 2), nsmall = 2),  
             p50_se = format(round(se.0.5,  digits = 2), nsmall = 2), 
             p75 =    format(round(`0.75`,  digits = 2), nsmall = 2), 
             p75_se = format(round(se.0.75, digits = 2), nsmall = 2),
             p90 =    format(round(`0.9`,  digits = 2), nsmall = 2), 
             p90_se = format(round(se.0.9, digits = 2), nsmall = 2)) %>%
      dplyr::select(subgroup, p10, p10_se, p25, p25_se, p50, p50_se, p75, p75_se, p90, p90_se)
    
    
    # merge all information together
    output_temp_subgroup <- 
      Reduce(function(x,y) merge(x = x, y = y, by = "subgroup", all = T), 
             list(sample_size, def, mean, geo_mean, percentiles)) %>%
      mutate(
        nutrient = select_biomarker,
        subgroup = as.character(subgroup),
        subgroup_name = template$subgroup)
    
    output_temp <- rbind.fill(output_temp, output_temp_subgroup) 
    
    #print(output_temp)   
  }
  output_temp <- output_temp %>%
    mutate(unit = template[, paste0(select_biomarker, "_unit")])
  # export data
  return(output_temp)
} 

#
# inflammation adjustment -------
#
inflammation_adj <- function(input_data, select_biomarker, template, operation){
  print(paste0("Inflammation adjustment for ", select_biomarker))
  
  input_data <-
    input_data %>%
    mutate(crp = ifelse(is.na(rep(template$crp, times = nrow(input_data))), NA, crp),
           agp = ifelse(is.na(rep(template$agp, times = nrow(input_data))), NA, agp),
           log_crp = log(ifelse(crp == 0, crp + 0.001, crp)),
           log_agp = log(ifelse(agp == 0, agp + 0.001, agp)),
           biomarker = input_data[, eval(parse(text = select_biomarker),)],
           log_biomarker = log(ifelse(biomarker == 0, biomarker + 0.001, biomarker)))
  
  #input_data <<- input_data
  
  input_data <- 
    input_data %>% 
    mutate(group = case_when( age < 0.5 ~ 1, 
                              (age >= 0.5 & age < 5) ~ 2,
                              (age >= 5 & age < 15) ~ 3, 
                              (age >= 15 & age < 50 & male == 0 & preg == 0) ~ 4,
                              (preg == 1) ~ 5,
                              (age >= 15 & male == 1) ~ 6, 
                              (age >= 50 & male == 0 & preg == 0) ~ 7))
  
  
  # calculate CRP/AGP reference
  input_data <-
    input_data %>% 
    mutate(log_crp_ref = case_when(group == 1 ~ log_crp,
                                   group == 2 ~ -2.26,
                                   group == 3 ~ log_crp,
                                   group == 4 ~ -1.83,
                                   group == 5 ~ log_crp,
                                   group == 6 ~ log_crp,
                                   group == 7 ~ log_crp
    ),
    log_agp_ref = case_when(group == 1 ~ log_agp,
                            group == 2 ~ -0.52,
                            group == 3 ~ log_agp,
                            group == 4 ~ -0.63,
                            group == 5 ~ log_agp,
                            group == 6 ~ log_agp,
                            group == 7 ~ log_agp
    ))
  
  # set log_crp_ref as log_crp for transferrin receptor
  if(select_biomarker %in% c("transferrin_receptor")){
    input_data <- 
      input_data %>% 
      mutate(log_crp_ref = log_crp)
  }
  
  # set log_crp_ref as log_crp for retinol binding protein
  if(select_biomarker %in% c("retinol_biding_protein", "serum_retinol")){
    input_data <- 
      input_data %>% 
      mutate(log_crp_ref = ifelse(group %in% c(3, 4), log_crp, log_crp_ref),
             log_agp_ref = ifelse(group %in% c(3, 4), log_agp, log_agp_ref))
  }
  
  #  set log_crp_ref as log_crp for serum zinc
  if(select_biomarker %in% c("serum_zinc")){
    input_data <- 
      input_data %>% 
      mutate(log_crp_ref = ifelse(group == 4, log_crp, log_crp_ref),
             log_agp_ref = ifelse(group == 4, log_agp, log_agp_ref))
  }
  
  ####
  # CRP/AGP reference values done
  ####
  
  # regression coefficients
  # Transferrin receptor only adjust agp
  if(select_biomarker %in% c("transferrin_receptor")){
    if(!is.error(lm(log_biomarker ~ log_agp, data = input_data[group == 2, ], na.action=na.omit))){
      mysvyglm <- lm(log_biomarker ~ log_agp, data = input_data[group == 2, ], na.action=na.omit)
      psc_beta2 <- coef(mysvyglm)[2]
      psc_beta3 <- 0
    } else {
      psc_beta2 = 0
      psc_beta3 = 0
    }
    
    if(!is.error(lm(log_biomarker ~ log_agp, data = input_data[group == 4, ], na.action=na.omit))){
      mysvyglm <- lm(log_biomarker ~ log_agp, data = input_data[group == 4, ], na.action=na.omit)
      wra_beta2 <- coef(mysvyglm)[2]
      wra_beta3 <- 0
    } else {
      wra_beta2 = 0
      wra_beta3 = 0
    }
  } 
  
  # when both AGP and CRP are available, other nutrients
  if(select_biomarker !="transferrin_receptor" & !is.na(template$agp) & !is.na(template$crp)) {
    if(!is.error(lm(log_biomarker ~ log_agp + log_crp, data = input_data[group == 2, ], na.action=na.omit))){
      mysvyglm <- lm(log_biomarker ~ log_agp + log_crp, data = input_data[group == 2, ], na.action=na.omit)
      psc_beta2 <- coef(mysvyglm)[2]
      psc_beta3 <- coef(mysvyglm)[3]
    } else {
      psc_beta2 = 0
      psc_beta3 = 0
    }
    
    if(!is.error(lm(log_biomarker ~ log_agp + log_crp, data = input_data[group == 4, ], na.action=na.omit))){
      mysvyglm <- lm(log_biomarker ~ log_agp + log_crp, data = input_data[group == 4, ], na.action=na.omit)
      wra_beta2 <- coef(mysvyglm)[2]
      wra_beta3 <- coef(mysvyglm)[3]
    } else {
      wra_beta2 = 0
      wra_beta3 = 0
    }
  }
  
  # Only AGP is available
  if(!(select_biomarker %in% c("transferrin_receptor")) & !is.na(template$agp) & is.na(template$crp)){
    if(!is.error(lm(log_biomarker ~ log_agp, data = input_data[group == 2, ], na.action=na.omit))){
      mysvyglm <- lm(log_biomarker ~ log_agp, data = input_data[group == 2, ], na.action=na.omit)
      psc_beta2 <- coef(mysvyglm)[2]
      psc_beta3 <- 0
    } else {
      psc_beta2 = 0
      psc_beta3 = 0
    }
    
    if(!is.error(lm(log_biomarker ~ log_agp, data = input_data[group == 4, ], na.action=na.omit))){
      mysvyglm <- lm(log_biomarker ~ log_agp, data = input_data[group == 4, ], na.action=na.omit)
      wra_beta2 <- coef(mysvyglm)[2]
      wra_beta3 <- 0
    } else {
      wra_beta2 = 0
      wra_beta3 = 0
    }
  } 
  
  # Only CRP is available
  if(!(select_biomarker %in% c("transferrin_receptor")) & is.na(template$agp) & !is.na(template$crp)){
    if(!is.error(lm(log_biomarker ~ log_crp, data = input_data[group == 2, ], na.action=na.omit))){
      mysvyglm <- lm(log_biomarker ~ log_crp, data = input_data[group == 2, ], na.action=na.omit)
      psc_beta2 <- 0
      psc_beta3 <- coef(mysvyglm)[2]
    } else {
      psc_beta2 = 0
      psc_beta3 = 0
    }
    
    if(!is.error(lm(log_biomarker ~ log_crp, data = input_data[group == 4, ], na.action=na.omit))){
      mysvyglm <- lm(log_biomarker ~ log_crp, data = input_data[group == 4, ], na.action=na.omit)
      wra_beta2 <- 0 
      wra_beta3 <- coef(mysvyglm)[2]
    } else {
      wra_beta2 = 0
      wra_beta3 = 0
    }
  } 
  
  # calculate the difference
  input_data <- 
    input_data %>%
    mutate(log_crp_diff = ifelse((((log_crp - log_crp_ref) <= 0) | is.na(log_crp)), 0, (log_crp - log_crp_ref)),
           log_agp_diff = ifelse((((log_agp - log_agp_ref) <= 0) | is.na(log_agp)), 0, (log_agp - log_agp_ref)))
  
  
  
  # apply adjustment algorithm 
  input_data <- 
    input_data %>%
    mutate(adj_biomarker = case_when(group == 2 ~ exp(log_biomarker - psc_beta2 * log_agp_diff - psc_beta3 * log_crp_diff),
                                     group == 4 ~ exp(log_biomarker - wra_beta2 * log_agp_diff - wra_beta3 * log_crp_diff),
                                     group != 2 & group !=4 ~ biomarker,
                                     is.na(group) ~ biomarker))
  
  
  
  
  #define which observation is being adjusted
  input_data <-
    input_data %>%
    mutate(
      adjustment = case_when(group %in% c(1, 3, 5, 6, 7) | is.na(group) ~ "Outside BRINDA adjustment scope",
                             #(group %in% c(2, 4) & is.na(log_crp) & is.na(log_agp)) ~ "No data on both AGP and CRP values",
                             (group %in% c(2, 4) & (log_crp_diff == 0) & (log_agp_diff == 0))     ~ "No adjustment",
                             (group %in% c(2, 4) & (log_crp_diff != 0  & (log_agp_diff == 0)))    ~ "Adjustment by only CRP",
                             (group %in% c(2, 4) & ((log_crp_diff ==0 | log_agp_diff != 0)))      ~ "Adjustment by only AGP",
                             (group %in% c(2, 4) & (log_crp_diff != 0 & log_agp_diff != 0))       ~ "Adjustment by both AGP & CRP"))
  
  # apply cut off
  input_data <- 
    input_data %>%
    mutate(
      cutoff = input_data[, eval(parse(text = paste0(select_biomarker, "_cutoff")))],
      def = case_when(operation == 1
                      ~ ifelse(adj_biomarker < cutoff, 1, 0),
                      operation == 2
                      ~ ifelse(adj_biomarker <= cutoff, 1, 0),
                      operation == 3
                      ~ ifelse(adj_biomarker >  cutoff, 1, 0),
                      operation == 4
                      ~ ifelse(adj_biomarker >=  cutoff, 1, 0)))  
  
  
  
  # remove the crp/age
  # input_data <- 
  #        input_data %>%
  #            select(-log_agp_ref, -log_crp_ref, -log_agp_diff, -log_crp_diff, -log_crp, -log_agp, - group)
  
  #export datad
  return(input_data)
}

#
# sequence function: look up table -----
#
sequence <- function(value, pop) {
    reference_sub <- 
        reference %>%
        filter(biomarker == select_biomarker & population == pop)
    seq_number <- 
        seq_len(nrow(reference_sub)) %>%
        map_lgl(~between(value, reference_sub$lower_age[.x], reference_sub$upper_age[.x])) %>%
        which()
    reference_value <- reference_sub[seq_number, "cutoff"]
    
    if(length(seq_number) == 0L) {
        reference_value <- NA
    }
    return(reference_value)
}
# sequence function over


#
# calculate deficiencies ------------
#
calculate_deficiencies <- function(input_data, reference, template) {
  # input_data = input_data to be analyzed
  # reference = biomarker reference table
  # template = the template for biomarker
  
  # create data frame for output data
  output <- NULL
  
  # calculate every biomarker
  for (select_biomarker in unique(reference$biomarker)) {
    # check if variable exists 
    
    if(!exists(select_biomarker, input_data)) {
      print(paste0(select_biomarker, " does not exist in ", template$data))
    }else if (sum(!is.na(input_data[[select_biomarker]])) == 0) {
      print(paste0("No valid values of ", select_biomarker, " in ", template$data))
    }else{
      # analyze individual biomarker
      print(paste0("Analyzing ", select_biomarker, " for ", template$data))
      
      # create default values
      input_data$biomarker <- input_data[[select_biomarker]]
      
      input_data <- 
          input_data %>%
            mutate(log_biomarker = 
                       ifelse(biomarker == 0, log(biomarker + 0.001), log(biomarker)))
      
      # assign weight; give dummpy wt_temp = 1 for single site
      if(template$complex_survey_design_yes_no == "yes") {
        input_data$wt_temp <- ifelse(is.na(input_data$biomarker), NA, 
                                  input_data[[paste0(select_biomarker, "_survey_weight")]])
        
      }
      
      if(template$complex_survey_design_yes_no == "no") {
        input_data <- 
          input_data %>%
          mutate(wt_temp = 1)
      }

      
      # !!!!
      # check the operator and create deficiency_risk values
      operation <-  case_when(reference[reference$biomarker == select_biomarker, "relationship"][1] == "<"  ~ 1, 
                              reference[reference$biomarker == select_biomarker, "relationship"][1] == "<=" ~ 2,
                              reference[reference$biomarker == select_biomarker, "relationship"][1] == ">"  ~ 3,
                              reference[reference$biomarker == select_biomarker, "relationship"][1] == ">=" ~ 4)
      
      # create cutoff, binary def variable
      input_data <- 
        input_data %>%
        mutate(
          #unit   = template[, paste0(select_biomarker, "_unit")],
          cutoff = case_when(male == 1 ~ map2_dbl(age, "male", sequence),
                             np_female == 1 ~ map2_dbl(age, "female", sequence),
                             np_female == 0 ~ map2_dbl(age, "pregnant women", sequence)),
          def = case_when(operation == 1
                          ~ ifelse(biomarker < cutoff, 1, 0),
                          operation == 2
                          ~ ifelse(biomarker <= cutoff, 1, 0),
                          operation == 3
                          ~ ifelse(biomarker >  cutoff, 1, 0),
                          operation == 4
                          ~ ifelse(biomarker >=  cutoff, 1, 0)))  
     
      output_temp <- parameters_analysis(input_data = input_data, select_biomarker = select_biomarker, template)

      # combine data
      output <- rbind.fill(output, output_temp)  
      
      # change variable name
      setnames(input_data, "def", paste0(select_biomarker, "_def"))
      #setnames(input_data, "unit", paste0(select_biomarker, "_unit"))
      input_data$biomarker <- NULL
      input_data$log_biomarker <- NULL
      
      #
      # inflammation adjusted values ---- high lighted
      # serum ferritin, RBP, and serum retinol
      if(((select_biomarker == "serum_ferritin") & (!is.na(template$crp) | !is.na(template$agp))) |
         ((select_biomarker == "retinol_binding_protein") & (!is.na(template$crp) | !is.na(template$agp))) |
         ((select_biomarker == "serum_retinol") & (!is.na(template$crp) | !is.na(template$agp))) |
         ((select_biomarker == "transferrin_receptor")) |
         ((select_biomarker == "serum_zinc") & (!is.na(template$crp) | !is.na(template$agp))))
      { # inflammation adjustment
        input_data <- inflammation_adj(input_data = input_data, select_biomarker = select_biomarker, template)
        
        
        input_data <-
          input_data %>%
          dplyr::select(-biomarker, -log_biomarker) %>%
          mutate(biomarker = adj_biomarker,
                 log_biomarker = log(ifelse(biomarker == 0, biomarker + 0.001, biomarker))) %>%
          dplyr::select(-adj_biomarker) |>
          mutate()
                def = case_when(operation == 1
                                ~ ifelse(biomarker < cutoff, 1, 0),
                                operation == 2
                                ~ ifelse(biomarker <= cutoff, 1, 0),
                                operation == 3
                                ~ ifelse(biomarker >  cutoff, 1, 0),
                                operation == 4
                                ~ ifelse(biomarker >=  cutoff, 1, 0)))  
        
        # export and analyze
        #input_data <<- input_data
        output_temp <- parameters_analysis(input_data = input_data, select_biomarker = select_biomarker, template)
        
        # rename the nutrient
        output_temp <- 
          output_temp %>%
          mutate(nutrient = paste0(select_biomarker, "_adj"))
        
        #print(output_temp)
        output <- rbind.fill(output, output_temp)
        
        # change dataset variable names
        setnames(input_data, "cutoff", paste0(select_biomarker, "_adj_cutoff"))
        setnames(input_data, "def", paste0(select_biomarker, "_adj_def"))
        setnames(input_data, "adjustment", paste0(select_biomarker, "_adj_method"))
        setnames(input_data, "biomarker", paste0(select_biomarker, "_adj"))
        
        # clean dataset # see if the values should be saved later
        
        setnames(input_data, "cutoff", paste0(select_biomarker, "_cutoff"))
        input_data <-     
          input_data %>%
          dplyr::select(-log_biomarker, -wt_temp)
        
        
      }
      # finish adjusted biomarker
    }
    # finish individual biomarker
    #print(input_data)
  }
  # analyze all biomarker except iodine
  
  
  # #Calculate iodine deficiency_risk
  # # might delete during publications
  # if(is.na(template$urine_iodine)){
  #   print(paste0("Urinary iodine does not exist in ", template$data))
  # } else {
  #   print(paste0("Analyzing urinary iodine for ", template$data))
  #   
  #   select_biomarker = "urine_iodine"
  #   
  #   # create default values
  #   input_data <-
  #     input_data %>%
  #     mutate(biomarker = input_data[, eval(parse(text = select_biomarker),)],
  #            log_biomarker = log(ifelse(biomarker == 0, biomarker + 0.001, biomarker)))
  #   
  #   # assign weight; give dummpy wt_temp = 1 for single site
  #   if(template$complex_survey_design_yes_no == "yes") {
  #     input_data <-
  #       input_data %>%
  #       mutate(wt_temp = ifelse(is.na(biomarker), 0, 
  #                               input_data[, eval(parse(text = paste0(select_biomarker, "_survey_weight")))]))
  #     
  #   }
  #   
  #   if(template$complex_survey_design_yes_no == "no") {
  #     input_data <- 
  #       input_data %>%
  #       mutate(wt_temp = 1)
  #   }
  #   
  #   # put a fake deficiency_risk variable 
  #   input_data$def = 1
  #   
  #   
  #   # calculate parameters
  #   #input_data <<- input_data
  #   output_temp <- parameters_analysis(input_data = input_data, select_biomarker = select_biomarker, template)
  #   
  #   # calcu;ate deficiency_risk using the regression method
  #   output_temp <-
  #     output_temp %>%
  #     mutate(deficiency_risk = format(round(11049 * as.numeric(p50)^(-1.63), digits = 2), nsmall = 2)) %>%
  #     select(-deficiency_risk_se)
  #   
  #   output <- rbind.fill(output, output_temp)
  #   
  #   input_data <-
  #     input_data %>% 
  #     select(-biomarker,  -log_biomarker, -def)
  # }
  # #Calculate iodine deficiency_risk END
  
  results <<- output %>%
    mutate(country = template$country,
           country_code = template$country_code,
           survey_year = template$survey_year,
           study_design = template$study_design) %>%
    dplyr::select("country_code", "country", "survey_year", "study_design", "nutrient", "unit", "subgroup_name", colnames(output)[2: length(colnames(output))])
  
  input_data <-
    input_data %>% 
    dplyr::select(-wt_temp)
  
  if(is.na(template$crp)){
    input_data$crp <- NULL
  }
  
  if(is.na(template$agp)){
    input_data$agp <- NULL
  }
  
  return(input_data)
}


#
# save data ------------
# 
save_data <- function(no_biodata,
                      final_biodata,
                      save_directory,
                      i,
                      current_biodata_template,
                      results_all,
                      results,
                      results_name) {

  fwrite(final_biodata, file = paste0(save_directory, "/biodata", i, "_clean_", current_biodata_template$file_name, ".csv"), row.names = F)
  print(paste0("-------------save `", i, "_clean_", current_biodata_template$file_name, "` under ", save_directory, " -------------"))
  # add results
  results_all <- rbind.fill(results_all, results)
  
  # remove the old data
  #rm(biodata, final_biodata)
  rm(results,  envir = .GlobalEnv)
  
  #
  # additional features to add
  # some variable must exist, such as ID
  # Need to add more checks, if a biomarker exists then its survey weight must exist
  # change biomarkers unit from Standard unit to International System (IS)
  # simulate the distribution
  # inflammation adjustment
  # To be complete
  #
  print(paste0("-------------Analysis of ", current_biodata_template$data, " is complete-------------"))
  print(paste0("--------------------------------------------------------------"))
  return(list(results_all = results_all, results = results))
} 
  
  
  



