
#
# inflammation function -------------
#
inflammation_adj <- function(input_data, select_biomarker, template){
    print(paste0("Inflammation adjustment for ", select_biomarker))
    
    input_data <-
        input_data %>%
        mutate(crp = ifelse(is.na(rep(template$crp, times = nrow(input_data))), NA, crp),
               agp = ifelse(is.na(rep(template$agp, times = nrow(input_data))), NA, agp),
               log_crp = log(ifelse(crp == 0, crp + 0.001, crp)),
               log_agp = log(ifelse(agp == 0, agp + 0.001, agp)))
    
    input_data$biomarker <- input_data[[select_biomarker]]
    
    input_data <- 
        input_data |>
            mutate(log_biomarker = log(ifelse(biomarker == 0, biomarker + 0.001, biomarker)))
    
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
    
    #
    # exception - sTfR
    # set log_crp_ref as log_crp for transferrin receptor
    if(select_biomarker %in% c("transferrin_receptor")){
        input_data <- 
            input_data %>% 
            mutate(log_crp_ref = log_crp)
    }
    
    #
    # Retinol
    # set log_agp_ref as log_agp for rbp and retinol 
    # set log_crp_ref as log_crp for retinol binding protein and retinol
    if(select_biomarker %in% c("retinol_biding_protein", "serum_retinol")){
        input_data <- 
            input_data %>% 
            mutate(log_crp_ref = ifelse(group %in% c(3, 4), log_crp, log_crp_ref),
                   log_agp_ref = ifelse(group %in% c(3, 4), log_agp, log_agp_ref))
    }
    
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
        
        rm(zn_agp_cor, zn_agp_P_value )
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
        
        if(max(agp_correlation, crp_correlation)){
            print("correlation between zinc and inflammation markers")
            input_data <- 
                input_data %>% 
                mutate(log_crp_ref = ifelse(group == 4, log_crp, log_crp_ref),
                        log_agp_ref = ifelse(group == 4, log_agp, log_agp_ref))
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
    
    #
    # when both AGP and CRP are available, other nutrients
    #
    if(select_biomarker !="transferrin_receptor" & !is.na(template$agp) & !is.na(template$crp)) {
        if(!is.error(lm(log_biomarker ~ log_agp + log_crp, data = input_data[input_data$group == 2, ], na.action=na.omit))){
            mysvyglm <- lm(log_biomarker ~ log_agp + log_crp, data = input_data[input_data$group == 2, ], na.action=na.omit)
            psc_beta2 <- coef(mysvyglm)[2]
            psc_beta3 <- coef(mysvyglm)[3]
        } else {
            psc_beta2 = 0
            psc_beta3 = 0
        }
        
        if(!is.error(lm(log_biomarker ~ log_agp + log_crp, data = input_data[input_data$group == 4, ], na.action=na.omit))){
            mysvyglm <- lm(log_biomarker ~ log_agp + log_crp, data = input_data[input_data$group == 4, ], na.action=na.omit)
            wra_beta2 <- coef(mysvyglm)[2]
            wra_beta3 <- coef(mysvyglm)[3]
        } else {
            wra_beta2 = 0
            wra_beta3 = 0
        }
    }
    
    #
    # Only AGP is available
    #
    if(!(select_biomarker %in% c("transferrin_receptor")) & !is.na(template$agp) & is.na(template$crp)){
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
            wra_beta2 <- 0
            wra_beta3 <- 0
        }
    } 
    
    #
    # Only CRP is available
    #
 
        if(!is.error(lm(log_biomarker ~ log_crp, data = input_data[input_data$group == 2, ], na.action=na.omit))){
            mysvyglm <- lm(log_biomarker ~ log_crp, data = input_data[input_data$group == 2, ], na.action=na.omit)
            psc_beta2 <- 0
            psc_beta3 <- coef(mysvyglm)[2]
        } else {
            psc_beta2 <- 0
            psc_beta3 <- 0
        }
        
        if(!is.error(lm(log_biomarker ~ log_crp, data = input_data[input_data$group == 4, ], na.action=na.omit))){
            mysvyglm <- lm(log_biomarker ~ log_crp, data = input_data[input_data$group == 4, ], na.action=na.omit)
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
        mutate(log_crp_diff = ifelse(((log_crp - log_crp_ref) <= 0), 0, (log_crp - log_crp_ref)),
               log_agp_diff = ifelse(((log_agp - log_agp_ref) <= 0), 0, (log_agp - log_agp_ref)))
    
# add sTfR    
    
    # apply adjustment algorithm 
if(!is.na(template$agp) & !is.na(template$crp)){
    input_data <- 
        input_data %>%
        mutate(adj_biomarker = 
                   case_when(group == 2 ~ exp(log_biomarker - psc_beta2 * log_agp_diff - psc_beta3 * log_crp_diff),
                             group == 4 ~ exp(log_biomarker - wra_beta2 * log_agp_diff - wra_beta3 * log_crp_diff),
                            group != 2 & group != 4 ~ NA_real_,
                            is.na(group) ~ NA_real_))
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