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
  connector <<- ifelse(Sys.info()['sysname'] == "Darwin","/","\\")
  
  #
  # Load the datasets
  # Provide error messages when datasets do not exist 
  #
  two_templates <- load_templates(biomarker_dataset_template , 
                                  biomarker_cutoff_template, 
                                  save_directory)
  biomarker_template <- two_templates$template 
  reference <- two_templates$reference
  
  rm(two_templates)
  
 
  # results template before the loop
  results_all <- NULL
  
  #
  # LOOP starts 
  #
  for (i in no_biodata) {
    #
    # Check if each parameter is correct -----------------
    #
    # transpose template to an easy-to-read form
      
    current_biodata_template <- check_format(template = biomarker_template, 
                                             no_biodata = i)
    
    biodata <- import_biodata(template = current_biodata_template, 
                              connector_type = connector)
    
    #
    # check variables - Reference ---------
    #
    variables_checked <- check_variables(reference_table = reference, 
                                         template = current_biodata_template)
    reference_after_checked_var <- variables_checked$reference  
    current_biodata_template_checked_var <- variables_checked$template 
    
    #
    # function: Create a new dataset with standarized variable name
    #
    biodata_transformed <- 
        transform_biodata(template = current_biodata_template_checked_var, 
                          input_data = biodata) 
    current_biodata_template_after_trans <- biodata_transformed$template
    biodata_after_trans <- biodata_transformed$final_biodata
      
    #
    # function: Change unit to standard units
    #
    converted_datasets <- convert_units(
        template = current_biodata_template_after_trans, 
        input_data = biodata_after_trans)
    current_biodata_template_final <- converted_datasets$template 
    final_biodata <- converted_datasets$output_data
    
    # clean the environment
    rm(biodata, biodata_after_trans, biodata_transformed, 
       converted_datasets, current_biodata_template,
       current_biodata_template_after_trans, 
       reference_after_checked_var, variables_checked, 
       current_biodata_template_checked_var)
    
    #
    # Actual calculate deficiencies ---------
    #
    final_biodata <- 
        calculate_deficiencies(input_data = final_biodata, 
                               reference_table = reference, 
                               template = current_biodata_template_final)

    # function: save clean data ---------
    # 
   
    # organize the variable list of the results_all
    data_saved <- save_data(final_biodata = final_biodata, 
                            save_directory = save_directory, 
                            i = i, 
                            template = current_biodata_template_final, 
                            results_all = results_all, 
                            results = results, 
                            results_name = results_name)
    
    results_all <- data_saved$results_all
  }  
  
  fwrite(results_all, file = paste0(save_directory, "/biodata_results_", results_name, ".csv"), row.names = F)
  message(paste0("-------------Analysis of all biomarker dataset(s) is complete-------------"))
}
# SAMBA main function ends --------------------------


#
# Load templates ----------
#
load_templates <- 
    function(biomarker_dataset_template, biomarker_cutoff_template, save_directory){
        # load the template  
        if (file.exists(biomarker_dataset_template)) {
            template <-  fread(file = biomarker_dataset_template, header = TRUE, 
                               sep = ",")  %>% as.data.frame()
            template$information <- str_replace(template$information, " ", "")
            rownames(template) <- template$information
            template$information <- NULL
        } else {
            stop ("biomarker_dataset_template does not exist! Please select the correct file")
        }
        
        # load reference table
        if (file.exists(biomarker_cutoff_template)) {
            reference <-  fread(file = biomarker_cutoff_template, header = TRUE, 
                                sep = ",")  %>% as.data.frame()
        } else {
            stop ("Biomarker_cutoff_template  template does not exist! Please select the correct file")
        }
        
        # check if the save directory exist
        if (!dir.exists(save_directory)) {
            stop ("Save directory does not exist! Please check the save directory")
        } 
        return(list(reference = reference, template = template))
    }

#
# check formats -----------
#
check_format <- 
    function(template, no_biodata){
        current_template <- 
            template[no_biodata] %>%
            t() %>%
            as.data.frame() %>%
            mutate_all(~na_if(., ""))
        
        #
        # add data - from different forms
        #
        current_template$file_format <- tolower(sapply(strsplit(current_template$data, "\\."), tail, 1))
        
        current_template$file_name <- substr(current_template$data, 1, 
                                             nchar(current_template$data) - nchar(current_template$file_format) - 1)
        
        if(!(current_template$file_format %in% c("csv", "dta", "sas7bdat", "sav", "xls", "xlsx", "xpt"))){
            stop(paste0(current_template$data, " in ", current_template$file_format, " format. 
    This is not a supported format. SAMBA function reads the CSV, excel, XPT, STATA, and SPSS data files"))
        }
        
        return(current_template)
    }

#
# Import biodata --------------
#
import_biodata <- 
    function(template, connector_type){
        # check if the biomarker dataset exists
        if (file.exists(paste0(template$directory, connector, template$data))) {
            message(paste0("-------------Loading ", template$directory, connector, template$data, "-------------"))
            
            # read csv
            if(template$file_format == "csv") {
                biodata <- as.data.table((fread(file = paste0(template$directory, connector, template$data), 
                                                header = T, sep = ",")))
            }
            
            # read SAS file
            if(template$file_format == "sas7bdat"){
                biodata <- as.data.table(read_sas(paste0(template$directory, connector, template$data)))
            }
            
            # read SPSS
            if(template$file_format == "sav"){
                biodata <- as.data.table(read_sav(paste0(template$directory, connector, template$data)))
            }
            
            # read STATA
            if(template$file_format == "dta"){
                biodata <- as.data.table(read_dta(paste0(template$directory, connector, template$data)))
            }
            
            # read excel
            if(template$file_format %in% c("xls", "xlsx")){
                biodata <- as.data.table(read_excel(paste0(template$directory, connector, template$data)))
            }
            
            # read xpt
            if(template$file_format == "xpt"){
                biodata <- as.data.table(read_xpt(paste0(template$directory, connector, template$data)))
            }
            
        } else {
            stop(paste0(template$data, " does not exist under ", template$directory))
        }
        
        # add ID if there is no unique ID
        if(is.na(template$id)){
            biodata$id = seq(nrow(biodata))
            template$id = "id"
        }
        
        return(biodata)
    }

#
# check variables ------------------
#
check_variables <- 
    function(reference_table, template){
        #
        # remove the space in variable name
        #
        char_values <- c("nutrients", "biomarker", "relationship", "cutoff",
                         "unit", "population",
                         "lower_age", "upper_age", "note")
        
        sapply(char_values, function(x) {
            reference_table[, x] <<- str_replace_all(reference_table[, x], " ", "")
        }) %>%
            invisible()
        rm(char_values)
        
        #
        # make sure all the numeric values were recorded into numeric values
        #
        numeric_var <- c("cutoff", "lower_age", "upper_age")
        
        sapply(numeric_var, function(x) {
            reference_table[, x] <<- as.numeric(reference_table[, x])}) %>%
            invisible()
        rm(numeric_var)
        
        # notification
        message(paste0("-------------Analyzing ", template$data, "-------------"))
        
        #
        # check if necessary inputs are specified
        #
        necessary_values <- c("country", "survey_year", "complex_survey_design_yes_no", 
                              "study_design", "sex_var", "age", "age_unit")
        sapply(necessary_values, function(x) {if(is.na(template[, x])) {
            stop(message(paste0("You must fill a value for ", x, " in the template")))}}) %>%
            invisible()
        rm(necessary_values)
        
        return(list(reference = reference_table, template = template))
    }

#
# transform biodata ---------
#
transform_biodata <- function(template, input_data){
    
    biomarker_var_survey_weight <- grep("_survey_weight.*", names(template), value=TRUE)
    
    biomarker_unit              <- grep("_unit.*", names(template), value=TRUE)[-1] # exclude age variable
    
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
        template[, x] <<- str_replace_all(template[, x], " ", "")
    }) %>%
        invisible()
    rm(char_values)
    
    # make sure all the numeric values were recorded into numeric values
    #numeric_var <- c("male_value", "pregnancy_value", "lactating_value")
    
    #sapply(numeric_var, function(x) {
    #    if(exists(numeric_var, template)) {
    #        template[, x] <<- as.numeric(template[, x])}}) %>%
    #    invisible()
    #rm(numeric_var)
    
    # if a biomarker variable is specified, its unit must be specified; also the age variable
    mapply(function(x, y) {
        if(!is.na(template[, x]) & is.na(template[, y])){
            stop(message(paste0("When you specify ", x, ", you must specify ",  y)))
        }
    },  c("age", biomarker_var), c("age_unit", biomarker_unit)) %>%
        invisible()
    
    # lower case the values of unit variables, and the values of yes/no variables
    sapply(c("age_unit", biomarker_unit, "PSU_ID_nested_in_strataID_yes_no",
             "complex_survey_design_yes_no"), function(x) {
                 template[, x] <<- tolower(template[, x])
             }) %>%
        invisible()
    
    # Check if the yes_no questions are with valid input
    if(!template$complex_survey_design_yes_no %in% c("yes", "n", "no", "n")){
        stop(message(paste0("Invalid input for complex_survey_design_yes_no")))
    }
    
    if(!template$PSU_ID_nested_in_strataID_yes_no %in% c("yes", "n", "no", "n", NA)){
        stop(message(paste0("Invalid input for complex_survey_design_yes_no")))
    }
    
    # check if the age variable has the correct input
    if(!(template$age_unit %in% c("year", "month"))){
        stop(message("Age unit can only be 'year' or 'month', please correct your input in the template"))
    }
    
    # if there is complex survey design, there must be weight for 
    if(template$complex_survey_design_yes_no %in% c("Yes", "Y")){
        mapply(function(x, y) {
            if(!is.na(template[, x]) & is.na(template[, y])){
                stop(message(paste0("When you specify ", x, ", you must specify ",  y)))
            }
        },  biomarker_var, biomarker_survey_weight) %>%
            invisible()
    }
    
    # If the pregnancy/lactating variables exist, there must be pregnancy/lactating variable value
    preg_lac_var <- c("pregnancy_var", "lactating_var")
    preg_lac_value <- c("pregnancy_value", "lactating_value")
    mapply(function(x, y) {
        if(!is.na(template[, x]) & is.na(template[, y])){
            stop(message(paste0("When you specify ", x, ", you must specify ",  y)))
        }
    },  preg_lac_var, preg_lac_value) %>%
        invisible()
    
    rm(preg_lac_value, preg_lac_var)
    
    #    
    # check every variables exist
    #
    covariates <- as.vector(str_split(template$covariates, " ", simplify =  T))
    
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
        if (!is.na(template[, var])) {
            if (exists(template[, var], input_data)) {
                # we need to select variables one by one, due to two weight variable can be the same variable
                final_biodata <- cbind(final_biodata, input_data[, template[, var], with = F])
                setnames(final_biodata, template[, var], var)
            } else {
                stop (paste0(template[, var], " does not exist in ", template$data))
            }
        }
    }
    
    for(var in biomarker_unit) {
        if(!is.na(template[, var])) {
            final_biodata <- final_biodata %>%
                mutate(unit_var = template[, var])
            setnames(final_biodata, "unit_var", var)
        }
    }
    
    # check covariates exist, if so, merge back to the dataset
    if(!is.na(template$covariates)) {
        for (var in covariates) {
            if (exists(var, input_data)) {
                # we need to select variables one by one, due to two weight variable can be the same variable
                final_biodata <- cbind(final_biodata, input_data[, var, with = F])
            } else {
                stop (paste0(var, "does not exist in ", template$data))
            }
        }
    }
    rm(covariates, var)
    
    # generate variables for male, female, pregnant
    if(exists("pregnancy_var", final_biodata)) {
        final_biodata <- 
            final_biodata %>%
            mutate(male = ifelse((sex_var == as.numeric(template$male_value)), 1, 0),
                   np_female = ifelse((sex_var != as.numeric(template$male_value) & (pregnancy_var != as.numeric(template$pregnancy_value) |is.na(pregnancy_var))), 1, 0),
                   preg = case_when((pregnancy_var == as.numeric(template$pregnancy_value)) ~ 1,
                                    pregnancy_var != as.numeric(template$pregnancy_value) | is.na(pregnancy_var) ~ 0))
    } else {
        final_biodata <- 
            final_biodata %>%
            mutate(male = ifelse((sex_var == as.numeric(template$male_value)), 1, 0),
                   np_female = ifelse((sex_var != as.numeric(template$male_value)), 1, 0),
                   preg = 0)  
    }
    
    # convert age into year
    if(template$age_unit == "year"){
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
        template$PSU = NA
    } 
    
    if(length(unique(final_biodata$strata)) == 1) {
        template$strata = NA
    } 
    
    return(list(final_biodata = as_tibble(final_biodata), template = template))
}
#
# transform biodata end ---------
#


#
# convert units ------------------
#
convert_units <- 
    function(template, input_data){
        
        #
        # serum b-12
        #
        if(!is.na(template$serum_b12)){
            if(!(template$serum_b12_unit %in% c("pg/ml", "pmol/l"))){
                stop(message("Serum/plasma B12 unit can only be 'pg/ml' or 'pmol/l'"))
            } 
            
            if(template$serum_b12_unit == "pg/ml") {
                message("convert serum b12 unit from pg/ml to pmol/l")
                input_data$serum_b12 = input_data$serum_b12 * 0.738
                template$serum_b12_unit = "pmol/l"
            }
        }
        
        #
        # serum folate
        #
        if(!is.na(template$serum_folate)) {
            
            if(!(template$serum_folate_unit %in% c("nmol/l", "ng/ml"))){
                stop(message("Serum folate unit can only be 'nmol/l' or 'ng/ml'"))
            } 
            
            if(template$serum_folate_unit == "ng/ml") {
                message("convert serum folate unit from ng/ml to nmol/l")
                input_data$serum_folate = input_data$serum_folate * 2.2655 
                template$serum_folate_unit = "nmol/l"
            }
        }
        
        #
        # RBC folate
        #
        if(!is.na(template$rbc_folate)){
            if(!(template$rbc_folate_unit %in% c("nmol/l", "ng/ml"))){
                stop(message("RBC folate unit can only be 'nmol/l' or 'ng/ml'"))
            } 
            
            if(template$rbc_folate_unit == "ng/ml") {
                message("convert RBC folate unit from ng/ml to nmol/l")
                input_data$rbc_folate = input_data$rbc_folate * 2.2655 
                template$rbc_folate_unit = "nmol/l"
            }
        }
        
        #
        # serum retinol
        #
        if(!is.na(template$serum_retinol)) {
            
            if(!(template$serum_retinol_unit %in% c("umol/l", "ug/dl"))){
                stop(message("Serum retinol unit can only be 'umol/l' or 'ug/dl'"))
            } 
            
            if(template$serum_retinol_unit == "ug/dl") {
                message("convert serum retinol unit from ug/dl to umol/l")
                input_data$serum_retinol = input_data$serum_retinol * 0.03491  
                template$serum_retinol_unit = "nmol/l"
            }
        }
        
        #
        # Retinol binding protein
        #
        if(!is.na(template$retinol_binding_protein)) {
            
            if(!(template$retinol_binding_protein_unit %in% c("umol/l", "ug/dl"))){
                stop(message("Retinol binding protein unit can only be 'umol/l' or 'ug/dl'"))
            } 
            
            if(template$retinol_binding_protein == "ug/dl") {
                message("convert serum retinol unit from ug/dl to umol/l")
                input_data$retinol_binding_protein = input_data$serum_retinol * 0.03491  
                template$serum_retinol_unit = "nmol/l"
            }
        }
        
        #
        # serum zinc
        #
        if(!is.na(template$serum_zinc)) {    
            # serum zinc
            if(!(template$serum_zinc_unit %in% c("umol/l", "ug/dl", "ug/l"))){
                stop(message("Serum zinc unit can only be 'umol/l' or 'ug/ml'"))
            } 
            
            if(template$serum_zinc_unit == "ug/dl") {
                message("convert serum zinc unit from ug/dl to umol/l")
                input_data$serum_zinc = input_data$serum_zinc * 0.15291
                
                template$serum_zinc_unit = "nmol/l"
            }
            
            if(template$serum_zinc_unit == "ug/l") {
                message("convert serum zinc unit from ug/l to umol/l")
                input_data$serum_zinc = input_data$serum_zinc * 0.015291
                
                template$serum_zinc_unit = "nmol/l"
            }
        }
        
        #
        # serum ferritin
        #
        if(!is.na(template$serum_ferritin)) {
            if(!(template$serum_ferritin_unit %in% c("ng/ml", "ug/l"))){
                stop(message("Serum ferritin unit can only be 'ng/ml' or 'ug/l'"))
            } 
            
            if(template$serum_ferritin_unit == "ng/ml") {
                message("convert serum ferritin unit from ng/ml to ug/l")
                input_data$serum_ferritin = input_data$serum_ferritin * 1 
                template$serum_ferritin_unit = "ug/l"
            }
        }
        
        # Soluable Transferrin receptor
        if(!is.na(template$transferrin_receptor)) {
            if(!(template$transferrin_receptor_unit %in% c("mg/l", "mg/dl"))){
                stop(message("Transferrin receptor unit can only be 'mg/l' or 'mg/dl'"))
            } 
            
            if(template$transferrin_receptor_unit == "mg/dl") {
                message("convert transferrin receptor unit from mg/dl to mg/l")
                input_data$transferrin_receptor = input_data$transferrin_receptor * 10 
                template$transferrin_receptor_unit = "mg/l"
            }
        }
        
        # #
        # # Iodine
        # #
        # if(!is.na(template$urine_iodine)) {
        #   if(!(template$urine_iodine_unit %in% c("ug/l", "ng/ml"))){
        #     stop(message("Urine iodine unit can only be 'ug/l' or 'ng/ml'"))
        #   } 
        #   
        #   if(template$urine_iodine_unit == "ng/ml") {
        #     message("convert urine iodine unit from ng/ml to ug/l")
        #     input_data$urine_iodine = input_data$urine_iodine * 1 
        #     template$urine_iodine_unit = "ug/l"
        #   }
        # }
        
        # CRP
        if(!is.na(template$crp)) {
            if(!(template$crp_unit %in% c("mg/l"))){
                stop(message("CRP unit can only be 'mg/l'"))
            } 
        }
        
        # AGP
        if(!is.na(template$agp)) {
            if(!(template$agp_unit %in% c("g/l"))){
                stop(message("AGP unit can only be 'g/l'"))
            } 
        }
        
        # vitamin D
        if(!is.na(template$serum_vitaminD)) {
            if(!(template$serum_vitaminD_unit %in% c("ng/ml", "nmol/l"))){
                stop(message("Serum Vitamin D unit can only be 'ng/ml' or 'nmol/l'"))
            } 
            
            if(template$serum_vitaminD_unit == "ng/ml") {
                message("convert serum vitamin D unit from ng/ml to nmol/l")
                input_data$serum_vitaminD = input_data$urine_iodine * 2.496 
                template$serum_vitaminD_unit = "nmol/l"
            }
        }
        
        return(list(template = template, output_data = input_data)) 
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
        
        #message(output_temp)   
    }
    output_temp <- output_temp %>%
        mutate(unit = template[, paste0(select_biomarker, "_unit")])
    # export data
    return(output_temp)
} 

#
# inflammation function -------------
#
inflammation_adj <- function(input_data, select_biomarker, template){
    #message(paste0("Inflammation adjustment for ", select_biomarker))
    
    input_data <-
        input_data %>%
        mutate(crp = ifelse(is.na(rep(template$crp, times = nrow(input_data))), NA, crp),
               agp = ifelse(is.na(rep(template$agp, times = nrow(input_data))), NA, agp),
               log_crp = log(ifelse(crp == 0, crp + 0.001, crp)),
               log_agp = log(ifelse(agp == 0, agp + 0.001, agp)))
    
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
    
    #
    # zinc
    #  set log_crp_ref as log_crp for serum zinc
    # add more criteria
    if(select_biomarker %in% c("serum_zinc")){
        
        agp_correlation <- NA
        crp_correlation <- NA
        
        if(!all(is.na(input_data$agp))){
            spearman_agp_results <- cor.test(input_data$biomarker, 
                                             input_data$agp, 
                                             method=c("spearman"), exact=F)
            zn_agp_cor = spearman_agp_results$estimate
            zn_agp_P_value = spearman_agp_results$p.value
            
            agp_correlation <- ifelse(zn_agp_cor < -0.1 & zn_agp_P_value < 0.1, 1, 0)
            
            #rm(zn_agp_cor, zn_agp_P_value )
        }
        
        if(!all(is.na(input_data$crp))){
            spearman_crp_results <- cor.test(input_data$biomarker, 
                                             input_data$crp, 
                                             method=c("spearman"), exact=F)
            zn_crp_cor = spearman_crp_results$estimate
            zn_crp_P_value = spearman_crp_results$p.value
            
            crp_correlation <- ifelse(zn_crp_cor < -0.1 & zn_crp_P_value < 0.1, 1, 0)
            
            rm(zn_crp_cor, zn_crp_P_value)
        }
        
        
        if(!is.na(agp_correlation) && agp_correlation == 1 | 
           !is.na(crp_correlation) && crp_correlation == 1){
            message("There is correlation between zinc and inflammation markers")
            zinc_correction <- 1
        }
    }
    
    ####
    # CRP/AGP reference values done
    ####
    psc_beta2 <- 0
    psc_beta3 <- 0
    
    wra_beta2 <- 0
    wra_beta3 <- 0
    
    # regression coefficients
    # Transferrin receptor only adjust agp
    if(select_biomarker %in% c("transferrin_receptor") & !is.na(template$agp)){
        if(!is.error(lm(log_biomarker ~ log_agp, data = input_data[input_data$group == 2, ], na.action=na.omit))){
            mysvyglm <- lm(log_biomarker ~ log_agp, data = input_data[input_data$group == 2, ], na.action=na.omit)
            psc_beta2 <- coef(mysvyglm)[2]
            psc_beta3 <- 0
        } else {
            psc_beta2 = 0
            psc_beta3 = 0
        }
        
        if(!is.error(lm(log_biomarker ~ log_agp, data = input_data[input_data$group == 4, ], na.action=na.omit))){
            mysvyglm <- lm(log_biomarker ~ log_agp, data = input_data[input_data$group == 4, ], na.action=na.omit)
            wra_beta2 <- coef(mysvyglm)[2]
            wra_beta3 <- 0
        } else {
            wra_beta2 = 0
            wra_beta3 = 0
        }
    } 
     
    
    # calculate the difference
    input_data <- 
        input_data %>%
        group_by(group) |>
        mutate(all_agp_na = all(is.na(agp)),
               all_crp_na = all(is.na(crp))) |>
        ungroup() |>
        mutate(log_agp_diff = case_when(all_agp_na == 1 ~ 0,
                                        (log_agp - log_agp_ref) <= 0 ~ 0,
                                        (log_agp - log_agp_ref) > 0 ~ (log_agp - log_agp_ref)),
               log_crp_diff = case_when(all_crp_na == 1 ~ 0,
                                        (log_crp - log_crp_ref) <= 0 ~ 0,
                                        (log_crp - log_crp_ref) > 0 ~ (log_crp - log_crp_ref)))
               
                            
    # apply adjustment algorithm 
      if(select_biomarker == "serum_ferritin"){
        input_data <- 
            input_data %>%
            mutate(adj_biomarker = 
                       case_when(group == 2 ~ 
                                   exp(log_biomarker - psc_beta2 * log_agp_diff - psc_beta3 * log_crp_diff),
                                 group == 4 ~ 
                                   exp(log_biomarker - wra_beta2 * log_agp_diff - wra_beta3 * log_crp_diff),
                                 group != 2 & group != 4 ~ NA_real_,
                                 is.na(group) ~ NA_real_))
      }
      
      if(select_biomarker == "transferrin_receptor"){
        input_data <- 
          input_data %>%
          mutate(adj_biomarker = 
                   case_when(group == 2 ~ exp(log_biomarker - psc_beta2 * log_agp_diff),
                             group == 4 ~ exp(log_biomarker - wra_beta2 * log_agp_diff),
                             group != 2 & group != 4 ~ NA_real_,
                             is.na(group) ~ NA_real_))
      }
      
      if(select_biomarker %in% c("retinol_binding_protein", "serum_retinol")){
        input_data <- 
          input_data %>%
          mutate(adj_biomarker = 
                   case_when(group == 2 ~ 
                               exp(log_biomarker - psc_beta2 * log_agp_diff - psc_beta3 * log_crp_diff),
                             group == 4 ~ biomarker,
                             group != 2 & group != 4 ~ NA_real_,
                             is.na(group) ~ NA_real_))
      }
      
      if(select_biomarker %in% c("serum_zinc")){
        input_data <- 
          input_data %>%
          mutate(adj_biomarker = 
                   case_when(group == 2 ~ zinc_correlation * exp(log_biomarker - psc_beta2 * log_agp_diff - psc_beta3 * log_crp_diff),
                             group == 4 ~ biomarker,
                             group != 2 & group != 4 ~ NA_real_,
                             is.na(group) ~ NA_real_))
      }
    }
    
    if(is.na(template$agp) & !is.na(template$crp)){
        input_data <- 
            input_data %>%
            mutate(adj_biomarker = 
                       case_when(group == 2 ~ exp(log_biomarker - psc_beta3 * log_crp_diff),
                                 group == 4 ~ exp(log_biomarker - wra_beta3 * log_crp_diff),
                                 group != 2 & group != 4 ~ NA_real_,
                                 is.na(group) ~ NA_real_))
    }
    
    if(!is.na(template$agp) & is.na(template$crp)){
        input_data <- 
            input_data %>%
            mutate(adj_biomarker = 
                       case_when(group == 2 ~ exp(log_biomarker - psc_beta2 * log_agp_diff),
                                 group == 4 ~ exp(log_biomarker - wra_beta2 * log_agp_diff),
                                 group != 2 & group != 4 ~ NA_real_,
                                 is.na(group) ~ NA_real_))
    }
    
    # remove the crp/age
    input_data <- 
        input_data %>%
        dplyr::select(-any_of(c("log_agp_ref", "log_crp_ref", "log_agp_diff", "log_crp_diff", "log_crp", "log_agp")))
    
    return(input_data)
}

#
# cutoff_assign: look up table -----
#
cutoff_assign <- function(pop, 
                          input_data,
                          reference_table = reference_table,
                          select_biomarker = select_biomarker){
    
    reference_sub <- 
        reference_table %>%
        filter(biomarker == select_biomarker & population == pop) |>
        dplyr::select(cutoff, lower_age, upper_age)
    
    #input_data$cutoff <- NA  
    
    for(k in seq(nrow(reference_sub))){
        
        reference_sub_one_age <- reference_sub[k, ]  
        
        input_data$cutoff_temp <- 
            ifelse(input_data$age >= reference_sub_one_age$lower_age & input_data$age <= reference_sub_one_age$upper_age, 
                   reference_sub_one_age$cutoff, NA)
        
        input_data$cutoff <-
            ifelse(!is.na(input_data$cutoff), 
                   input_data$cutoff, 
                   input_data$cutoff_temp)
    }
    
    return(input_data$cutoff)
}
# sequence function over


#
# calculate deficiencies ------------
#
calculate_deficiencies <- function(input_data, 
                                   reference_table, 
                                   template) {
    # input_data = input_data to be analyzed
    # reference = biomarker reference table
    # template = the template for biomarker
    
    # create data frame for output data
    output <- NULL
    
    # calculate every biomarker
    for (select_biomarker in unique(reference_table$biomarker)) {
        # check if variable exists 
        
        if(!exists(select_biomarker, input_data)) {
            message(paste0(select_biomarker, " does not exist in ", template$data))
        }else if (sum(!is.na(input_data[[select_biomarker]])) == 0) {
            message(paste0("No valid values of ", select_biomarker, " in ", template$data))
        }else{
            # analyze individual biomarker
            message(paste0("Analyzing ", select_biomarker, " for ", template$data))
            
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
            
            #
            # check the operator and create deficiency_risk values
            #
            operation <-  case_when(reference_table[reference_table$biomarker == select_biomarker, "relationship"][1] == "<"  ~ 1, 
                                    reference_table[reference_table$biomarker == select_biomarker, "relationship"][1] == "<=" ~ 2,
                                    reference_table[reference_table$biomarker == select_biomarker, "relationship"][1] == ">"  ~ 3,
                                    reference_table[reference_table$biomarker == select_biomarker, "relationship"][1] == ">=" ~ 4)
            
            #
            # create cutoff, binary def variable
            #
            input_data$cutoff <- NA
            
            input_data$cutoff[which(input_data$male == 1)] <- cutoff_assign(pop = "male", reference_table = reference_table,select_biomarker = select_biomarker, input_data = input_data[which(input_data$male == 1), ])
            
            input_data$cutoff[which(input_data$np_female == 1)] <- cutoff_assign(pop = "female", reference_table = reference_table,select_biomarker = select_biomarker, input_data = input_data[which(input_data$np_female == 1), ])
            
            input_data$cutoff[which(input_data$preg == 1)] <- cutoff_assign(pop = "pregnant women", reference_table = reference_table,select_biomarker = select_biomarker, input_data = input_data[which(input_data$preg == 1), ])
            
            input_data <-
                input_data |>
                mutate(
                    def = case_when(operation == 1 ~ ifelse(biomarker < cutoff, 1, 0),
                                    operation == 2 ~ ifelse(biomarker <= cutoff, 1, 0),
                                    operation == 3 ~ ifelse(biomarker >  cutoff, 1, 0),
                                    operation == 4 ~ ifelse(biomarker >=  cutoff, 1, 0)))  
            
            output_temp <- parameters_analysis(input_data = input_data, 
                                               select_biomarker = select_biomarker, 
                                               template = template)
            
            # combine data
            output <- rbind.fill(output, output_temp)  
            
            # change variable name
            setnames(input_data, "def", paste0(select_biomarker, "_def"))
            #setnames(input_data, "unit", paste0(select_biomarker, "_unit"))
            #input_data$biomarker <- NULL
            #input_data$log_biomarker <- NULL
            
            #
            # inflammation adjusted values ---- high lighted
            # serum ferritin, RBP, and serum retinol
            if(((select_biomarker == "serum_ferritin") & (!is.na(template$crp) | !is.na(template$agp))) |
               ((select_biomarker == "retinol_binding_protein") & (!is.na(template$crp) | !is.na(template$agp))) |
               ((select_biomarker == "serum_retinol") & (!is.na(template$crp) | !is.na(template$agp))) |
               ((select_biomarker == "transferrin_receptor") & !is.na(template$agp)) |
               ((select_biomarker == "serum_zinc") & (!is.na(template$crp) | !is.na(template$agp))))
            { # inflammation adjustment
                input_data <- inflammation_adj(input_data = input_data, select_biomarker = select_biomarker, template = template)
                
                input_data <-
                    input_data %>%
                    dplyr::select(-biomarker, -log_biomarker) %>%
                    rename(biomarker = adj_biomarker) |>
                    mutate(log_biomarker = log(ifelse(biomarker == 0, biomarker + 0.001, biomarker)),
                           def = case_when(operation == 1
                                           ~ ifelse(biomarker < cutoff, 1, 0),
                                           operation == 2
                                           ~ ifelse(biomarker <= cutoff, 1, 0),
                                           operation == 3
                                           ~ ifelse(biomarker >  cutoff, 1, 0),
                                           operation == 4
                                           ~ ifelse(biomarker >=  cutoff, 1, 0))) 
                
                
                #input_data <<- input_data
                output_temp <- parameters_analysis(input_data = input_data, select_biomarker = select_biomarker, template)
                
                # rename the nutrient
                output_temp <- 
                    output_temp %>%
                    mutate(nutrient = paste0(select_biomarker, "_adj"))
                
                #message(output_temp)
                output <- rbind.fill(output, output_temp)
                
                # change dataset variable names
                setnames(input_data, "cutoff", paste0(select_biomarker, "_cutoff"))
                setnames(input_data, "def", paste0(select_biomarker, "_adj_def"))
                setnames(input_data, "biomarker", paste0(select_biomarker, "_adj"))
                
                input_data <-     
                    input_data %>%
                    dplyr::select(-log_biomarker, -wt_temp, -group)
                
                
            }
            # finish adjusted biomarker
        }
        # finish individual biomarker
        #message(input_data)
    }
    # analyze all biomarker except iodine
    
    results <<- output %>%
        mutate(country = template$country,
               country_code = template$country_code,
               survey_year = template$survey_year,
               study_design = template$study_design) %>%
        dplyr::select("country_code", "country", "survey_year", "study_design", "nutrient", "unit", "subgroup_name", everything())
    
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

#save_data(no_biodata, final_biodata, save_directory, i, current_biodata_template, results_all, results, results_name)

save_data <- function(final_biodata,
                      save_directory,
                      i,
                      template = current_biodata_template_final,
                      results_all,
                      results,
                      results_name) {
    
    fwrite(final_biodata, file = paste0(save_directory, "/biodata", i, "_clean_", template$file_name, ".csv"), row.names = F)
    message(paste0("-------------save `", i, "_clean_", template$file_name, "` under ", save_directory, " -------------"))
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
    message(paste0("-------------Analysis of ", template$data, " is complete-------------"))
    message(paste0("--------------------------------------------------------------"))
    return(list(results_all = results_all, results = results))
}   



