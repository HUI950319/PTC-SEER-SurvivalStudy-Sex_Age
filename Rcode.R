# ==============================================================================
# PAPILLARY THYROID CANCER SURVIVAL ANALYSIS
# ==============================================================================
# Title: Sex disparities in papillary thyroid cancer survival: divergent patterns
#        of relative and absolute effects across the age spectrum
# 
# Study Design: Retrospective cohort study using SEER database (2004-2015)
# Sample Size: 77,349 PTC patients (Female: 61,197; Male: 16,152)
# 
# Primary Outcome: Cancer-specific survival (CSS)
# Statistical Methods: 
#   - Cox proportional hazards models
#   - Fine-Gray competing risk models  
#   - Restricted cubic spline (RCS) analysis
#   - Interaction analysis (multiplicative and additive scales)
# 
# Adjustment Methods: Unadjusted, Direct standardization, IPTW, PSM
# ==============================================================================
  
#########--------  Load R packages --------  ################ 
# Package management and installation
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  # Survival analysis packages for adjusted curves and competing risks
  adjustedCurves, riskRegression, prodlim, survival, QHScrnomo, tidycmprsk,
  
  # Visualization packages for plots and forest plots
  patchwork, gridExtra, forestplot, forestploter, 
  
  # Data processing and manipulation packages
  tidyverse, 
  
  # Reporting and table formatting packages
  gtsummary, broom.helpers,
  
  # Interaction analysis packages for RCS and interaction effects
  interactionRCS, interactionR
)

#########--------  Functions    --------  ################

# Main function: creates adjusted survival or cumulative incidence function objects
# This is the primary wrapper function for different causal inference methods
# Parameters:
#   data: input dataset
#   cat_var: categorical variable (exposure/treatment)
#   adj_var: adjustment variables for confounding control
#   method: adjustment method (unadjusted, direct standardization, IPTW, PSM)
#   times: time points for analysis (in months)
#   surv: TRUE for survival analysis, FALSE for competing risks
#   conf_level: confidence level for intervals
.get.obj <- function(data, cat_var, adj_var, 
                     method = c("unadj", "direct", "iptw", "psm"),
                     times = c(60,120,180,240),
                     surv = T,
                     conf_level = 0.95
                     ){
  
  # Determine time points for analysis
  # If times is NULL or all specified times exist in data, use all available times
  if(is.null(times) || all(times %in% data$time) ){
    timepoints = NULL
  } else {
    timepoints = c( times, unique(data$time) )
  } 
  
  
  
  # Nested function: fit statistical models for direct standardization method
  # Fits Cox model for survival analysis or Fine-Gray model for competing risks
  .get.fit_for.direct <- function(data, co_var, surv = T){
    
    # Generate appropriate formula based on analysis type
    if (isTRUE(surv)) {
      # Cox proportional hazards model formula for survival analysis
      formula <- as.formula(sprintf("Surv(time, DSS) ~   %s", paste(co_var, collapse = "+") ))
    } else {
      # Fine-Gray competing risks model formula
      formula <- as.formula(sprintf("Hist(time, event) ~ %s", paste(co_var, collapse = "+") ))
    }
    print(formula)
    
    # Fit the appropriate statistical model
    if (isTRUE(surv)) {
      # Fit Cox proportional hazards model
      fit <- coxph(formula, data = data, x=TRUE)
    } else {
      # Fit Fine-Gray competing risks model
      fit <- FGR(formula, data = data, cause=1)
    }
    
    fit$call$formula <- formula
    
    return(fit)
  }
  
  # Nested function: fit propensity score models for PSM and IPTW methods
  # Creates treatment assignment models for causal inference methods
  .get.fit_for.PSM.IPTW <- function(data, cat_var, adj_var){
    
    # Remove categorical variable from adjustment variables to avoid collinearity
    adj_var <- setdiff(adj_var, cat_var)
    
    # Create propensity score model formula
    formula <- as.formula(sprintf("%s ~ %s", cat_var, paste(adj_var, collapse = "+") ))
    
    # Fit appropriate model based on number of treatment levels
    if(length(levels(data[[cat_var]])) == 2L){
      # Binary treatment: use logistic regression
      fit <- glm(formula, data, family="binomial")
    }else{
      # Multiple treatments: use multinomial regression
      fit <- nnet::multinom(formula, data)
    }
    
    fit$call$formula <- formula
    
    return(fit)
  }
  
  # Nested function: create adjusted survival curves using different causal inference methods
  # Applies various adjustment techniques to estimate causal effects on survival
  .get.sur.obj <- function(data, cat_var, adj_var, method = c("unadj", "direct", "iptw", "psm"),
                           timepoints = NULL, conf_level = 0.95){
    
    # Prepare variables and create treatment group variable
    adj_var <- setdiff(adj_var, cat_var)
    data <- data %>% mutate(group = .data[[cat_var]])
    
    # Apply the specified adjustment method
    if(method == "unadj"){
      # Unadjusted analysis using Kaplan-Meier estimator
      obj <- adjustedsurv(data=data,
                          variable= "group",
                          ev_time="time", event="DSS",
                          times = timepoints,
                          method="km", conf_int=T, conf_level = conf_level,
                          bootstrap=F, n_boot= 1000, n_cores=20)
      
      
    }else if(method == "direct"){
      # Direct standardization using Cox regression model
      cox_mod <- .get.fit_for.direct(data,c("group",adj_var),surv = T)
      obj <- adjustedsurv(data=data,
                          variable= "group",
                          ev_time= "time", event= "DSS",
                          times = timepoints,
                          method= "direct",conf_int= T, conf_level = conf_level,
                          outcome_model = cox_mod,
                          bootstrap=F, n_boot=1000, n_cores=20)
      
    }else if(method == "iptw"){
      # Inverse probability of treatment weighting (IPTW)
      iptw_mod <- .get.fit_for.PSM.IPTW(data, "group", adj_var)
      obj <- adjustedsurv(data=data,
                          variable="group",
                          ev_time="time", event="DSS",
                          times = timepoints,
                          method="iptw_km",conf_int=T, conf_level = conf_level,
                          treatment_model=iptw_mod,
                          bootstrap=F, n_boot=1000, n_cores=20)
      
    }else if(method == "psm"){
      # Propensity score matching (only for binary treatments)
      if(length(levels(data[[cat_var]]))!=2L) stop( "levels of group must be two" )
      psm_mod <- .get.fit_for.PSM.IPTW(data, "group", adj_var)
      obj <- adjustedsurv(data=data,
                          variable="group",
                          ev_time="time", event="DSS",
                          times = timepoints,
                          method="matching",conf_int=T, conf_level = conf_level,
                          treatment_model=psm_mod,
                          bootstrap=T, n_boot=50, n_cores=20)
    }
    
    # Generate survival plots for visualization (when appropriate)
    if(  is.null(timepoints) || length(timepoints) > 10 ){
      
      if(  method == "psm" ){
        # Use bootstrap confidence intervals for PSM
        plot(obj, conf_int = T,use_boot = T, steps = F) %>% print
      }else{
        # Use analytical confidence intervals for other methods
        plot(obj, conf_int = T, steps = F) %>% print
      }
      
    }
    
    
    return(obj)
    
  }
  
  .get.cif.obj <- function(data, cat_var, adj_var, method = c("unadj", "direct", "iptw", "psm"),
                           timepoints = NULL, conf_level = 0.95){
    library(cmprsk)
    library(riskRegression)
    library(prodlim)
    
    
    adj_var <- setdiff(adj_var, cat_var)
    data <- data %>% mutate(group = .data[[cat_var]])
    
    if(method == "unadj"){
      obj <- adjustedcif(data=data,
                         variable="group",
                         ev_time="time", event= "event", cause=1,
                         times = timepoints,
                         method="aalen_johansen",conf_int=T,conf_level = conf_level,
                         bootstrap=F, n_boot=1000, n_cores=20)
      
    }else if(method == "direct"){
      fgr_mod <- .get.fit_for.direct(data,c("group",adj_var),surv = F)
      obj <- adjustedcif(data=data,
                         variable="group",
                         ev_time="time", event="event",cause=1,
                         times = timepoints,
                         method="direct",conf_int=T,conf_level = conf_level,
                         outcome_model= fgr_mod,
                         bootstrap=T, n_boot=50, n_cores=20)
      
    }else if(method == "iptw"){
      iptw_mod <- .get.fit_for.PSM.IPTW(data, "group", adj_var)
      obj <- adjustedcif(data=data,
                         variable="group",
                         ev_time="time", event="event",cause=1,
                         times = timepoints,
                         method="iptw",conf_int=T, conf_level = conf_level,
                         treatment_model=iptw_mod,
                         bootstrap=F, n_boot=1000, n_cores=20)
      
    }else if(method == "psm"){
      if(length(levels(data[[cat_var]]))!=2L) stop( "levels of group must be two" )
      psm_mod <- .get.fit_for.PSM.IPTW(data, "group", adj_var)
      obj <- adjustedcif(data=data,
                         variable="group",
                         ev_time="time", event="event",cause=1,
                         times = timepoints,
                         method="matching",conf_int=T, conf_level = conf_level,
                         treatment_model=psm_mod,
                         bootstrap=T, n_boot=50, n_cores=20)
    }
    
    if(  is.null(timepoints) || length(timepoints) > 10 ){
      
      if(  method %in% c("psm", "direct") ){
        plot(obj, conf_int = T,use_boot = T, steps = F) %>% print
      }else{
        plot(obj, conf_int = T, steps = F) %>% print
      }
      
    }
    
    return(obj)
    
    
  }
  
  
  # get res
  if( isTRUE(surv) ){
    obj <- .get.sur.obj(data, cat_var, adj_var, method,timepoints, conf_level)
  } else {
    obj <- .get.cif.obj(data, cat_var, adj_var, method,timepoints, conf_level)
  }
  
  if (method == "psm" || (!surv && method == "direct")) {
    use_boot = T
    cat("------------------------------------------------------------------\n",
        "--- NOTE: The confidence interval was estimated with bootstrap ---\n",
        "------------------------------------------------------------------\n")
  } else {
    use_boot = F
  } 
  
  attr(obj,"use_boot")<- use_boot
  return(obj)
}

# calculates survival differences between groups across different time points and subgroups
.get.sub.dif.times <- function(data, cat_var, sub_var, adj_var, surv = T, method = "unadj",
                               group_1 = NULL,group_2 = NULL,
                               times1 = c(60,120,180),
                               times2 = c(120,180),
                               conf_level1 = 0.95, 
                               conf_level2 = 0.95){
  # all     = "Age"
  
  
  # function 
  .get.sub.dif0 <- function(data, cat_var, adj_var, surv = T, method = "unadj",
                            group_1 = NULL,group_2 = NULL,
                            times1 = c(60,120,180,240),
                            times2 = c(120),
                            conf_level1 = 0.95, 
                            conf_level2 = 0.95){
    
    obj <-  .get.obj(data, cat_var = cat_var, adj_var = adj_var, 
                     method = method,
                     times = times1,
                     surv = surv,
                     conf_level = conf_level1)
    
    if(is.null(group_1) || is.null(group_2) ) {
      
      group_1 <- levels(data[[cat_var]])[1]
      group_2 <- levels(data[[cat_var]])[2]
      
    }
    
    
    .get.dif.from.obj(obj, 
                      times = times2, 
                      group_1 = group_1,group_2 = group_2, 
                      conf_level = conf_level2)
    
    
  }
  
  
  
  .get.sub.dif_Onevar_Multime <- function(data, cat_var, sub_var, adj_var, surv = T, method = "unadj",
                                          group_1 = NULL,group_2 = NULL,
                                          times1 = c(36, 60,120,180),
                                          times2 = c(60,120),
                                          conf_level1 = 0.95, 
                                          conf_level2 = 0.95){
    
    ### get and fmt res
    # fmt adj_var                        
    adj_var <- setdiff(adj_var, c(cat_var, sub_var))
    dat_lis <- split(data, data[[sub_var]])
    levels <- names(dat_lis)
    
    # res_sub and res_sum
    res_sub <- 
      map_dfr(1:length(dat_lis), ~ .get.sub.dif0(dat_lis[[.x]], cat_var, adj_var, surv, method,
                                                 group_1,group_2,
                                                 times1,times2,
                                                 conf_level1,conf_level2 ) %>% 
                mutate(label = levels[.x] ,.before = 1) %>% 
                add_row(label = sub_var, .before = 1)
      ) 
    
    
    res_sub <- 
      res_sub %>% 
      mutate(time = paste0(format(as.character(time/12), justify = "right" ), "-years")) %>% 
      mutate(across(where(is.character), ~ if_else(grepl("NA", .), "", . )  ))
    
    
    return(res_sub)
    
  } 
  
  if(!is.null(sub_var)){
    res <- 
      .get.sub.dif_Onevar_Multime(data, cat_var, sub_var, adj_var,surv,method,
                                  group_1,group_2,
                                  times1 ,
                                  times2 , 
                                  conf_level1, conf_level2)
    
    
    # res_sum <- .get.sub.summary(data, sub_var, add_all = F)
    res_sum <- .get_summary_for_forest(data, cat_var, sub_var, add_all = F)
    res <- left_join(  res, res_sum, by = "label" )
    
    # if( !all(res_sum$label == res$label) ) {
    #   stop( "the label in summary data  is not equal to the label in forest data"  )
    # } else {
    #   res <- cbind(res_sum, res %>% dplyr::select(-c("label")) ) 
    # }
    
  } else {
    
    res <- 
      .get.sub.dif0(data, cat_var, adj_var, surv, method,
                    group_1,group_2,
                    times1,times2,
                    conf_level1,conf_level2) %>% 
      mutate(label = "All patients" ,.before = 1) %>% 
      add_row(label = "group1-2", .before = 1) %>% 
      mutate(time = paste0(format(as.character(time/12), justify = "right" ), "-years")) %>% 
      mutate(across(where(is.character), ~ if_else(grepl("NA", .), "", . )  ))
    
  }
  
  
  
  
  
    res <- res %>% mutate(label = ifelse(label %in% c(names(data), "All patients" ), label, paste0("    ", label) ))
  
  
  attr(res, "data.names") <- c(names(data), "All patients" )
  attr(res, "surv") <- surv
  attr(res, "legend_name") <- cat_var
  attr(res, "legend_value") <- levels(data[[cat_var]])
  attr(res, "group_1") <- group_1
  attr(res, "group_2") <- group_2
  attr(res, "time") <- sprintf("%s years", paste( times2/12, collapse = "."))
  attr(res, "method") <- method
  attr(res, "model") <- if(surv) "Cox" else "Crr"
  attr(res, "header_title") <- paste(attr(res, "time"), attr(res, "method"), attr(res, "model"), "model")
  attr(res, "measure") <- "sur"
  
  
  res
  
}


# Univariate and Multivariate Analysis Functions
.get.UM.haz <- function(data, vars, adj_var, surv = T, method = "unadj"){
  
  time = "time" 
  outcome = "DSS"
  adj_var <- if(method == "unadj") NULL else adj_var
  
  # crr model
  if(!isTRUE(surv)) outcome = "status"
  if( (dplyr::n_distinct(data[[outcome]]) == 3L) & ("numeric" %in% class(data[[outcome]])) ){ data[[outcome]] <- as.factor(data[[outcome]])}
  
  
  # function
  .get.fit <- function(data, co_var, time = "time", outcome = "DSS"){
    
    # get formula
    if (!is.null(time) && !is.null(outcome)) {
      formula <- as.formula(sprintf("Surv(time = %s, event = %s) ~ %s", time, outcome, paste(co_var, collapse = "+") ))
    } else {
      formula <- as.formula(sprintf("%s ~ %s", outcome, paste(co_var, collapse = "+") ))
    }
    
    print(formula)
    
    
    # get fit
    if (!is.null(time) && !is.null(outcome)) {
      if (length(unique(data[[outcome]])) == 2L) {
        fit <- coxph(formula, data = data)
      } else {
        fit <- tidycmprsk::crr(formula, data = data)
      }
    } else {
      fit <- glm(formula, data = data, family = binomial)
    }
    
    fit$call$formula <- formula
    
    return(fit)
  }
  .fmt.fit <- function(model,include){
    t <- model %>% tbl_regression(exponentiate = TRUE, add_estimate_to_reference_rows = T, include = include)
    
    res <- t$table_body %>% mutate(n_obs.percent = n_obs/N_obs*100, n_event.percent = n_event/N_event*100) %>%
      mutate(n1 = if_else(is.na(n_event) | is.na(n_obs),  "",  paste0(sprintf(fmt = "%.0f",n_event), "/", sprintf(fmt = "%.0f",n_obs)))) %>%
      mutate(n1 = format(n1, justify = "centre")) %>%
      mutate(n_obs = if_else(is.na(n_obs), "", sprintf(fmt = "%.0f",n_obs)),
             n_obs = format(n_obs, justify = "right")) %>%
      mutate(n_obs.percent = if_else(is.na(n_obs.percent), "", sprintf(fmt = "(%.1f)",n_obs.percent)),
             n_obs.percent = format(n_obs.percent, justify = "right")) %>%
      mutate(n = paste(n_obs, n_obs.percent, sep = " "),
             n= format(n, justify = "right")) %>%
      dplyr::select(label,n,n1, n_obs, n_event, n_obs.percent, n_event.percent, estimate, conf.low, conf.high, p.value) %>%
      
      mutate(conf.low = if_else(estimate == 1, 1, conf.low),
             conf.high = if_else(estimate == 1, 1, conf.high)) %>%
      mutate(effect = ifelse(
        is.na(estimate) | is.na(conf.low) | is.na(conf.high), " ",
        sprintf(fmt = fmt_ci(),estimate,conf.low,conf.high))) %>%
      mutate(effect = if_else(is.na(conf.low), "",
                              if_else(conf.low == 1, "Reference", effect))) %>%
      
      dplyr::select(label, n, n1, estimate, conf.low, conf.high, effect, p.value) %>%
      rename(`Number(%)` = n,`No. of event/total` = n1, p_value = p.value ) %>%
      mutate(effect. = effect) %>% 
      mutate(effect = case_when(
        p_value < 0.001 ~ paste0(effect, "***"),  
        p_value < 0.01  ~ paste0(effect, "**"),   
        p_value < 0.05  ~ paste0(effect, "*"),   
        TRUE            ~ effect                 
      )) %>% 
      mutate(p_value. = fmt_signif(p_value, digits = 3) )  
    
    
    return(res)
    
    
    
  }
  .get.Uni <- function(data, vars, time = "time", outcome = "DSS"){
    
    # get Uni res
    
    res <-   map_dfr(vars, ~ .get.fit(data, .x, time, outcome) %>% .fmt.fit(include = .x)  )
    return(res)
    
  }
  .get.Mul <- function(data, vars, adj_var, time = "time", outcome = "DSS"){
    
    # get Mul res
    co_var <- c(vars, setdiff(adj_var, vars))
    res <- 
      .get.fit(data, co_var, time, outcome) %>% 
      .fmt.fit(include = vars)
    
    return(res)
    
  }
  chech_res_to_forest <- function(data){
  
    original_attrs <- attributes(data)
    
    data <- data %>%
      rowwise() %>%
      mutate(across(c(estimate, conf.low, conf.high), 
                    ~ if (any(c_across(c(estimate, conf.low, conf.high)) == 0 | 
                              is.infinite(c_across(c(estimate, conf.low, conf.high))), na.rm = TRUE)) {
             
                      print(label)
                      1
                    } else {
                      .
                    }
      )) %>%
      ungroup()
    
    attributes(data) <- original_attrs
    
    return(data)
  }
  
  
  # get res
  res <- if(method == "unadj") .get.Uni(data, vars, time, outcome) else .get.Mul(data, vars, adj_var, time, outcome)
  
  # # print
  # res[c("label", "effect.", "p_value.")] %>% autoReg::myft(digits = 3)%>% print
  
  # fmt res and change my summary data
  # res_sum <- .get.sub.summary(data, vars, add_all = F)
  res_sum <-.get_summary_for_forest(data, cat_var = NULL, sub_var = vars, add_all = F)
  # res <- left_join( res_sum, res, by = "label" )   
  
  if( !all(res_sum$label == res$label) ) {
    stop( "the label in summary data  is not equal to the label in forest data"  )
  } else {
    res <- cbind(res_sum, res %>% dplyr::select(-c("label", "Number(%)","No. of event/total")) ) 
  }
  

  res <- res %>% mutate(label = ifelse(label %in% names(data), label, paste0("    ", label) ))
  
  # check res
  res <- chech_res_to_forest(res)
  
  
  attr(res, "data.names") <- c(names(data), "All patients" )
  # attr(res, "cat_var") <- cat_var
  attr(res, "model") <- if(surv) "Cox" else "Crr"
  attr(res, "method") <- if(is.null(adj_var)) "Unadj" else "Adjus"
  attr(res, "header_title") <- paste(attr(res, "method"), attr(res, "model"), "model")
  attr(res, "measure") <- "haz"
  
  
  
  return(res)
  
}
.get.UM.sur <- function(data, vars, adj_var, surv = T, method = "unadj"){
  
  times1 = c(60,120,180)
  times2 = c(120)
  conf_level1 = 0.95
  conf_level2 = 0.95
  
  # function 
  # one cat_var
  .get.UM.dif0 <- function(data, cat_var, adj_var, surv = T, method = "unadj",                            
                           times1 = c(60,120,180,240),
                           times2 = c(120),
                           conf_level1 = 0.95, 
                           conf_level2 = 0.95){
    
    obj <-  .get.obj(data, cat_var = cat_var, adj_var = adj_var, 
                     method = method,
                     times = times1,
                     surv = surv,
                     conf_level = conf_level1)
    
    
    
    res <- .get.dif.from.obj(obj,  times = times2, 
                             group_1 = NULL,group_2 = NULL, 
                             conf_level = conf_level2)
    
    
    group_level <- levels(data[[cat_var]])
    ref_level <- group_level[1]
    
    res <- 
      res %>% 
      filter(group2 %in% ref_level) %>% 
      dplyr::select(time, group, group1, estimate, conf.low, conf.high, effect, p_value) %>% 
      add_row(time = times2, group = ref_level, group1 = ref_level,
              estimate = 0,conf.low =0, conf.high=0, effect = "reference",
              .before =1 ) %>% 
      rename(label = group1 )
    
    
    res_sur <- 
      .get.sur.from.obj(obj, times = times2) %>% 
      rename(label = group) %>% 
      mutate(label = factor(label, levels = group_level)) %>% arrange(label) %>% mutate(label = as.character(label)) %>%  
      add_row( time = times2,   label = cat_var, .before = 1) %>% 
      dplyr::select(c("time", "label", "estimate", "conf.low", "conf.high", "effect")) %>% 
      rename_with( ~paste0(., "_",  "sur"), .cols = c("estimate", "conf.low", "conf.high", "effect") ) %>% 
      rename(Effect_sur =  "effect_sur")
    
    
    left_join(res_sur, res, by = c("label", "time") )
    
    
  }
  
  
  res <-  map_dfr( vars  , ~ .get.UM.dif0(data, .x, adj_var, surv, method ,                           
                                          times1,
                                          times2, 
                                          conf_level1, 
                                          conf_level2 ) )
  
  # res_sum <- .get.sub.summary(data, vars, outcome = "event", add_all = F)
  res_sum <-.get_summary_for_forest(data, cat_var = NULL, sub_var = vars, add_all = F)
  
  
  if( !all(res_sum$label == res$label) ) {
    stop( "the label in summary data  is not equal to the label in forest data"  )
  } else {
    res <- cbind(res_sum, res %>% dplyr::select(-c("label")) ) 
  }
  
  res <- res %>% mutate(label = ifelse(label %in% names(data), label, paste0("    ", label) ))
  
  
  # # check res
  # res <- chech_res_to_forest(res)
  
  
  attr(res, "data.names") <- c(names(data), "All patients" )
  attr(res, "surv") <- surv
  # attr(res, "legend_name") <- cat_var
  # attr(res, "legend_value") <- levels(data[[cat_var]])
  attr(res, "group_1") <- "Effect_sur"
  # attr(res, "group_2") <- group_2
  attr(res, "time") <- sprintf("%s years", times2/12)
  attr(res, "method") <- method
  attr(res, "model") <- if(surv) "Cox" else "Crr"
  attr(res, "header_title") <- paste(attr(res, "time"), attr(res, "method"), attr(res, "model"), "model")
  attr(res, "measure") <- "sur"
  
  
  
  return(res)
  
  
  
}
.get.UM.lis <- function(data, cat_var, adj_var,surv = T, 
                        methods = c("unadj","adjus","unadj","direct" )){
  
  res1 <- .get.UM.haz(data, cat_var, adj_var,surv,method = methods[1] )
  res2 <- .get.UM.haz(data, cat_var, adj_var,surv,method = methods[2] )
  res3 <- .get.UM.sur(data, cat_var, adj_var,surv,method = methods[3] )
  res4 <- .get.UM.sur(data, cat_var, adj_var,surv,method = methods[4] )
  
  res <- list(res1, res2, res3, res4)
  
  res <- map(res, ~ slice(.x, -1))
  
  res
}


# Univariate and Multivariate Analysis by Subgroups Functions
.get.SUB.haz <- function(data, cat_var, sub_var, adj_var, surv, method){
  
  #get res all
  res_all <- 
    .get.UM.haz(data,cat_var, adj_var, surv, method) %>%
    mutate(label0 = "All patients", .before=1 )
  
  
  # get res sub
  dat_lis <- group_split(data, .data[[sub_var]])
  sub_levels <- levels(data[[sub_var]])
  adj_var <- setdiff(adj_var, c(cat_var, sub_var))
  
  res_sub <-   
    1:length(sub_levels) %>% 
    map_dfr (   ~.get.UM.haz(dat_lis[[.x]] ,cat_var, adj_var, surv, method) %>% 
                  mutate(label0 = sub_levels[.x], .before=1 ))
  
  
  return(bind_rows(res_all, res_sub))
  
}
.get.SUB.sur <- function(data, cat_var, sub_var, adj_var, surv, method){
  
  #get res all
  res_all <- 
    .get.UM.sur(data,cat_var, adj_var, surv, method) %>% 
    mutate(label0 = "All patients", .before=1 )
  
  #get res sub
  dat_lis <- group_split(data, .data[[sub_var]])
  sub_levels <- levels(data[[sub_var]])
  adj_var <- setdiff(adj_var, c(cat_var, sub_var))
  
  res_sub <-   
    1:length(sub_levels) %>% 
    map_dfr (   ~.get.UM.sur(dat_lis[[.x]] ,cat_var, adj_var, surv, method) %>% 
                  mutate(label0 = sub_levels[.x], .before=1 ))
  
  return(bind_rows(res_all, res_sub))
  
}

.get.SUB.lis1 <- function(data, cat_var, sub_var, adj_var, surv = T,
                          methods = c("unadj","adjus","unadj","direct" )){
  
  group_1 <- levels(data[[cat_var]])[1]
  group_2 <- levels(data[[cat_var]])[2]
  res5 <- .get.sub.haz(data, cat_var, sub_var, adj_var, surv, method = methods[1] ) 
  res6 <- .get.sub.haz(data, cat_var, sub_var, adj_var, surv, method = methods[2] ) 
  res7 <- .get.sub.dif(data, cat_var, sub_var, adj_var, surv, method = methods[3] ,  group_1 = group_1, group_2 = group_2 ) 
  res8 <- .get.sub.dif(data, cat_var, sub_var, adj_var, surv, method = methods[4] , group_1 = group_1, group_2 = group_2 )
  
  list(res5, res6, res7, res8)
}
.get.SUB.lis2 <- function(data, cat_var, sub_var, adj_var,surv = T, 
                          methods = c("unadj","adjus","unadj","direct" )){
  
  res1 <- .get.SUB.haz(data, cat_var, sub_var, adj_var, surv, method = methods[1])
  res2 <- .get.SUB.haz(data, cat_var, sub_var, adj_var, surv, method = methods[2])
  res3 <- .get.SUB.sur(data, cat_var, sub_var, adj_var, surv, method = methods[3])
  res4 <- .get.SUB.sur(data, cat_var, sub_var, adj_var, surv, method = methods[4])
  
  list(res1, res2, res3, res4)
  
}



# Univariate and Multivariate Analysis Functions
# Function: Restricted Cubic Spline (RCS) analysis with Cox model and bootstrap confidence intervals
# Analyzes non-linear age-sex interaction effects on hazard ratios using RCS
# Parameters:
#   data: input dataset
#   var1: categorical variable (e.g., sex)
#   var2: continuous variable for RCS (e.g., age)
#   xlim: age range for analysis
#   dead: outcome variable name
#   n: number of knots for RCS
#   co_var: covariates for adjustment
#   conf: confidence level
#   ci.boot.method: bootstrap CI method
#   R: number of bootstrap replicates
#   parallel: parallel processing method
rcs.cox1_boot <- function(data, var1, var2, xlim=c(0,80), dead = "DSS", n=3, co_var = NULL,
                          conf = 0.95 , ci.boot.method = "perc", R = 100, parallel = "multicore"){
  
  co_var <- setdiff(co_var, c(cat_var, con_var))
  
  data <- data %>% rename(dead = dead) %>%
    mutate(across(all_of(var1), ~as.numeric(.)-1)) %>%
    filter(.data[[var2]] %in% c(xlim[1]:xlim[2]))
  
  
  # Set data to environments
  pos <- 1
  envir = as.environment(pos)
  assign("ddist_", rms::datadist(data), envir = envir)
  options(datadist = "ddist_")
  
  
  rcs_var <- paste(var1, " * rcs(", var2, ",", n, ")", sep = "")
  
  if(is.null(co_var)){
    model <- cph(as.formula(paste('Surv(time, dead)',rcs_var, sep=" ~ ")), data = data, x = TRUE , y=TRUE)
  }else{
    model <- cph(as.formula(paste(paste('Surv(time, dead)', paste(co_var, collapse=" + "), sep=" ~ "),rcs_var,  sep=" + ")), data = data, x = TRUE , y=TRUE)
  }
  
  
  x_Lim <- seq(xlim[1],xlim[2],1)
  
  
  
  .bootrcsHR_cat_on_con <- function(data, idx , x , model , var1 , var2){
    form_str <- deparse(formula(model), width.cutoff = 500)
    has_vector_rcs <- "^.+ ~ .+\\s*\\*?\\s*rcs\\([^,]+, c\\(.*\\)\\).*$"
    matches_vector_rcs <- stringr::str_detect(form_str, has_vector_rcs)
    df <- data[idx,]
    mycall <- model$call
    mycall <- pryr::modify_call(mycall , list(data=quote(df)))
    if("cph" %in% class(model) | "lrm" %in% class(model)){
      myformula <- model$sformula
    } else {
      myformula <- model$formula
    }
    
    mymodel <- eval(mycall)
    coefMod <- coef(mymodel)
    
    if("cph" %in% class(model) | "lrm" %in% class(model)){
      separator <- " * "
      k <- mymodel$Design$parms[[var2]]
    } else {
      if (matches_vector_rcs) {
        matches <- stringr::str_match(form_str, "rcs\\((.*)\\)")
        rcs_content <- matches[, 2]
        split_str <- strsplit(rcs_content, ",")[[1]]
        var <- trimws(split_str[1])
        knots_str <- paste(split_str[-1], collapse = ",")
        knots <- eval(parse(text = paste("c(", knots_str, ")", sep="")))
        k <- knots
        separator <- ":"
        rcsTerm <- grep(":" , grep("rcs\\(" , attributes(mymodel$terms)$term.labels , value = TRUE) , invert = TRUE , value = TRUE)
        names(coefMod) <- gsub(rcsTerm , "" , names(coefMod) , fixed = TRUE)
      }   else {
        separator <- ":"
        # need to recreate the knot sequence for object of class coxph
        # remove the rcs part from the names of the variables in case of a class coxph model
        rcsTerm <- grep(":" , grep("rcs\\(" , attributes(model$terms)$term.labels , value = TRUE) , invert = TRUE , value = TRUE)
        names(coefMod) <- gsub(rcsTerm , "" , names(coefMod) , fixed = TRUE)
        indices <- length(grep(paste0("^", var2, "'*$"), names(coefMod)))+1
        k <- attributes(rms::rcs(df[[var2]], indices))$parms
      }
    }
    
    
    ### Get Iterative Fractions
    k_num = length(k)
    k_num_counter = length(k)
    
    fraction_list = list()
    j = 1
    while (k_num_counter >2) {
      temp_frac = (k[k_num] - k[j]) / (k[k_num]-k[k_num-1])
      temp_frac_2 = (k[k_num-1] - k[j]) / (k[k_num]-k[k_num-1])
      fraction_list = c(fraction_list,temp_frac)
      fraction_list = c(fraction_list,temp_frac_2)
      j = j + 1
      k_num_counter = k_num_counter - 1
    }
    
    ### New Variable Definition
    var_list <- c(var1, var2, paste(var1 , var2 , sep = separator),paste(var2 , var1 , sep = separator))
    
    for(i in 1:(length(k)-2)){
      var2_mod <- paste0(var2, paste(rep("'", i), collapse = ""))
      var_list <- c(var_list,
                    paste0(var2_mod),
                    paste(var1, var2_mod, sep = separator),
                    paste(var2_mod, var1, sep = separator))
    }
    
    
    myvars <- sort(intersect(var_list, names(coefMod)))
    
    mycoef <- coefMod[  myvars ]
    mycoefWhich <- sort(sapply( myvars , function(v) which( names(coefMod) %in% v )))
    
    num_ticks <- length(k)-2
    b <- mycoef[ var1 ]
    a <- NULL
    
    for (i in 0:num_ticks) {
      
      var_name <- paste0(var2, paste(rep("'", i), collapse = ""))
      
      if (var_name %in% names(mycoef)) {
        a <- c(a, mycoef[var_name])
      }
    }
    
    l <- mycoef[setdiff(setdiff(myvars, c(var1, names(a))), c(var1, names(a)))]
    
    ### Get Iterative Terms
    
    k_num = length(k)
    k_num_counter = length(k)
    term_list = list()
    j = 1
    frac_index = 1
    while (k_num_counter >2) {
      numer <- vapply(x , function(i) {
        max(i - k[j],0)^3  - (max(i - k[k_num-1],0)^3)*fraction_list[[frac_index]] + (max(i -             k[k_num],0)^3)*fraction_list[[frac_index+1]]
      } , numeric(1))
      denom <- (k[k_num] - k[1])^2
      numDem <- numer/denom
      
      term_list = c(term_list,numDem)
      
      k_num_counter = k_num_counter - 1
      frac_index = frac_index + 2
      j = j + 1
    }
    
    n <- length(term_list)
    elements_per_iteration <- length(x)
    
    subsets <- list()
    
    for (i in seq(1, n, by = elements_per_iteration)) {
      end <- min(i + elements_per_iteration - 1, n)
      subsets[[length(subsets) + 1]] <- term_list[i:end]
    }
    
    sp2 <- vapply(x , function(i) l[1]*i , numeric(1))
    
    for(i in 2:length(l)) {
      sp2 <- sp2 + l[[i]] * unlist(subsets[i-1])
    }
    
    
    HR <- unname(exp( b + sp2))
  }
  
  
  myBoot <- boot::boot(data = data, statistic = .bootrcsHR_cat_on_con,
                       x = x_Lim , model = model, var1 = var1 , var2 = var2,
                       R = R , parallel = parallel)
  
  SE <- apply(myBoot$t , 2 , sd)
  HRci <- t(vapply( seq_len(length(x_Lim)) , function(idx) {
    bci <- boot::boot.ci(boot.out = myBoot,  index = idx , type = ci.boot.method , conf = conf)
    if(ci.boot.method %in% "norm"){
      c(bci$t0 , bci$normal[2] , bci$normal[3] , SE = SE[idx])
    } else {
      c(bci$t0 , bci[[4]][4] , bci[[4]][5] , SE = SE[idx])
    }
  } , numeric(4)))
  HRci <- cbind(x_Lim , HRci)
  colnames(HRci) <- c("x" , "yhat" , "lower" , "upper" , "SE")
  rownames(HRci) <- x_Lim
  HRci <- as.data.frame(HRci)
  class(HRci) <- c("HR1" , class(HRci))
  
  
  
  ## Get P nonlinear
  res.p_nonlin.my <- function(HRci, var1){
    
    dat <- HRci
    
    # Set data to environments
    pos <- 1
    envir = as.environment(pos)
    assign("ddist_dat", rms::datadist(dat), envir = envir)
    options(datadist = "ddist_dat")
    
    fit <- Glm(yhat ~ rcs(x,3), data = dat)
    
    pdata  <- stats::anova(fit) %>% as.data.frame()
    p.value <- pdata[" Nonlinear","P"]
    p.value <- format_pvalue(p.value, digits = 3)
    if (!regex_detect(p.value, "<", fixed = TRUE)) {
      p.value <- paste0("P for nonlinear", " = ", p.value)
    } else{
      p.value <- regex_replace(p.value, "<", replacement = "", fixed = TRUE)
      p.value <- paste0("P for nonlinear", " < ", p.value)
    }
    
    if(exists("ddist_dat")){ rm("ddist_dat", inherits = TRUE, envir = envir)}
    p.string <- paste0(var1, ": ", p.value,  "\n")
    
  }
  
  p.string <- res.p_nonlin.my(HRci, var1)
  
  
  
  # Delete data from environments
  if(exists("ddist_")){
    rm("ddist_", inherits = TRUE, envir = envir)
  }
  if(exists("m_")){
    rm("m_", inherits = TRUE, envir = envir)
  }
  
  
  attr(HRci , "raw_dat") <- list(myBoot = myBoot,
                                 data = data,
                                 model = model,
                                 p.string = p.string,
                                 var1 = var1,  ## 分组变量
                                 var2 = var2, n= n, xlim = xlim, X_Lim = x_Lim,
                                 dead = dead,  co_var = co_var,
                                 conf = conf , ci.boot.method = ci.boot.method, R = R, parallel = parallel
  )
  return(HRci)
  
  
}

# Function: Restricted Cubic Spline (RCS) analysis with Poisson model and bootstrap confidence intervals
# Analyzes non-linear age-sex interaction effects on absolute mortality rates using RCS
# Parameters:
#   data: input dataset
#   var1: categorical variable (e.g., sex) 
#   var2: continuous variable for RCS (e.g., age)
#   xlim: age range for analysis
#   PA: person-years scaling factor
#   dead: outcome variable name
#   n: number of knots for RCS
#   co_var: covariates for adjustment
#   conf: confidence level
#   ci.boot.method: bootstrap CI method
#   R: number of bootstrap replicates
#   parallel: parallel processing method
rcs.poi_boot <- function(data, var1, var2, xlim=c(0,80), PA=12*1000, dead = "DSS", n=3, co_var = NULL,
                         conf = 0.95 , ci.boot.method = "perc", R = 100, parallel = "multicore"){
  
  
  co_var <- setdiff(co_var, c(cat_var, con_var))
  
  data <- data %>% rename(x = var2, dead = dead) %>% filter(x %in% c(xlim[1]:xlim[2]))
  
  # Set data to environments
  pos <- 1
  envir = as.environment(pos)
  assign("ddist_", rms::datadist(data), envir = envir)
  options(datadist = "ddist_")
  
  
  
  rcs_var <- paste("rcs(", "x", ",", n, ")", sep = "")
  
  if(is.null(co_var)){
    model <- Glm(as.formula(paste0("dead ~ offset(log(time/", PA, "))", "+",  rcs_var)), data = data, family=poisson)
  }else{
    model <-  Glm(as.formula(paste(paste0("dead ~ offset(log(time/", PA, "))", "+",  rcs_var), paste(co_var, collapse=" + "), sep =  "+")), data = data, family=poisson)
  }
  
  x_Lim <- seq(xlim[1],xlim[2],1)
  
  
  
  .bootPY_group <- function(data, idx, x_Lim = x_Lim, var1 = var1, model= model){
    
    
    if(is.null(var1)){
      df <- data[idx, ]
      mycall <- model$call
      mycall <- pryr::modify_call(mycall , list(data=quote(df)))
      mymodel <- eval(mycall)
      PY <- rms::Predict(mymodel, x = x_Lim, fun=exp, ref.zero = F, offset=list(time=PA))%>% pull(yhat)
      PY
    }else{
      if(!is.factor(data[[var1]])){stop("var1 must be a factor")}
      
      df <- data[idx, ]
      df_lis <- split(df, df[[var1]])
      PY_df <- map_dfr(df_lis, function(df){
        mycall <- model$call
        mycall <- pryr::modify_call(mycall , list(data=quote(df)))
        mymodel <- eval(mycall)
        
        PY_df<-rms::Predict(mymodel, x = x_Lim, fun=exp, ref.zero = F, offset=list(time=PA)) %>% as.data.frame() %>%  select(x, yhat)
        PY_df })
      
      PY <- PY_df %>% pull(yhat)
      PY
    }
    
  }
  
  # .bootPY_group(data, idx, x_Lim ,var1,model)
  
  
  myBoot <- boot::boot(data = data, statistic = .bootPY_group,
                       x_Lim = x_Lim, var1 = var1, model = model,
                       R = R , parallel = parallel)
  
  
  
  if(is.null(var1)){
    SE <- apply(myBoot$t , 2 , sd)
    LINci <- t(vapply( seq_len(length(x_Lim)) , function(idx) {
      bci <- boot::boot.ci(boot.out = myBoot,  index = idx , type = ci.boot.method , conf = conf)
      if(ci.boot.method %in% "norm"){
        c(bci$t0 , bci$normal[2] , bci$normal[3] , SE = SE[idx])
      } else {
        c(bci$t0 , bci[[4]][4] , bci[[4]][5] , SE = SE[idx])
      }
    } , numeric(4)))
    
    LINci <- cbind(x_Lim , LINci)
    colnames(LINci) <- c("x" , "yhat" , "lower" , "upper" , "SE")
    rownames(LINci) <- x_Lim
    LINci <- as.data.frame(LINci) #%>% mutate(lower = if_else(lower<0, 0, lower))
    
  }else{
    
    myBoot$t[is.na(myBoot$t)] <- 0  
    SE <- apply(myBoot$t , 2 , sd)
    LINci <- t(vapply( seq_len(length(x_Lim)*length(levels(data[[var1]]))) , function(idx) {
      bci <- boot::boot.ci(boot.out = myBoot,  index = idx , type = ci.boot.method , conf = conf)
      if(ci.boot.method %in% "norm"){
        c(bci$t0 , bci$normal[2] , bci$normal[3] , SE = SE[idx])
      } else {
        c(bci$t0 , bci[[4]][4] , bci[[4]][5] , SE = SE[idx])
      }
    } , numeric(4)))
    LINci <- cbind(x_Lim , LINci)
    colnames(LINci) <- c("x" , "yhat" , "lower" , "upper" , "SE")
    rownames(LINci) <- rep(x_Lim, times = length(levels(data[[var1]])))
    LINci <- as.data.frame(LINci) %>% mutate(group = rep(levels(data[[var1]]), each = length(x_Lim)))
  }
  
  
  ## Get P nonlinear
  res.p_nonlin <- function(data, model, var1){
    if(is.null(var1)){
      
      pdata  <- stats::anova(model) %>% as.data.frame()
      p.value <- pdata[" Nonlinear","P"]
      p.value <- format_pvalue(p.value, digits = 3)
      if (!regex_detect(p.value, "<", fixed = TRUE)) {
        p.value <- paste0("P for total nonlinear", " = ", p.value)
      } else{
        p.value <- regex_replace(p.value, "<", replacement = "", fixed = TRUE)
        p.value <- paste0("P for total nonlinear", " < ", p.value)
      }
      p.string <- paste0( p.value,  "\n")
    }else{
      df_lis <- split(data, data[[var1]])
      P_lis <- map(df_lis, function(df){
        
        pos <- 1
        envir = as.environment(pos)
        assign("ddist_df", rms::datadist(df), envir = envir)
        options(datadist = "ddist_df")
        
        
        mycall <- model$call
        mycall <- pryr::modify_call(mycall , list(data=quote(df)))
        mymodel <- eval(mycall)
        pdata  <- stats::anova(mymodel) %>% as.data.frame()
        p.value <- pdata[" Nonlinear","P"]
        p.value <- format_pvalue(p.value, digits = 3)
        
        if (!regex_detect(p.value, "<", fixed = TRUE)) {
          p.value <- paste0("P for nonlinear", " = ", p.value)
        } else{
          p.value <- regex_replace(p.value, "<", replacement = "", fixed = TRUE)
          p.value <- paste0("P for nonlinear", " < ", p.value)
        }
        
        if(exists("ddist_df")){ rm("ddist_df", inherits = TRUE, envir = envir)}
        p.string <- paste0(unique(df[[var1]]),": ", p.value,  "\n")
      })
      
      p.string <- do.call(paste0, P_lis)
      
    }
    
    
  }
  
  p.string <-  res.p_nonlin(data, model, var1)
  
  
  
  # Delete data from environments
  if(exists("ddist_")){rm("ddist_", inherits = TRUE, envir = envir)}
  

  attr(LINci , "raw_dat") <- list(myBoot = myBoot,
                                  data = data,
                                  model = model,
                                  p.string = p.string,
                                  var1 = var1,  
                                  var2 = var2, n    = n,xlim = xlim, X_Lim = x_Lim, PA   = PA,
                                  dead = dead,
                                  co_var = co_var,
                                  conf = conf , ci.boot.method = ci.boot.method, R = R, parallel = parallel
  )
  return(LINci)
  
}



### Visualization Functions

.plot.obj <- function(obj, color = lancet_palette,
                      xlim = c(0,240),xbre = seq(0,240,24),
                      ylim = NULL, ybre = NULL,
                      conf_int = T, risk_table = F){
  
  xexp <- 0.02; yexp <- 0.02
  xlab <- "Time since diagnosis (month)" 
  
  
  
  
  # function
  .plot.sur <- function(obj,color = lancet_palette, conf_int = T, risk_table = F){
    
    # get max_t & x_breaks & use_boot parameters 
    { 
      data <- obj$data
      if( max(data$time) > 200 ) {
        x_interval <- 24
      } else{
        x_interval <- 12
      }
      max_t <- floor( max(data$time)/x_interval ) * x_interval
      x_breaks <- seq(0, max_t, x_interval)
      
      if( !is.null(attr(obj, "use_boot")) ){
        use_boot <- attr(obj, "use_boot")
      } else {
        use_boot <- if( "boot_adj" %in% names(obj) ) TRUE else FALSE
      }
      
    }
    
    # plot
    plot(obj, conf_int = conf_int, use_boot = use_boot,
         max_t = max_t, 
         x_breaks = x_breaks,
         steps = F,
         custom_colors = color,
         conf_int_alpha = 0.2,
         gg_theme = km_theme,
         risk_table = risk_table,risk_table_stratify=T, risk_table_ylab = NULL,
         legend.title = "group", legend.position = "bottom") 
    
    
  }

  .plot.cif <- function(obj,color = lancet_palette, conf_int = T){
    
    { 
      data <- obj$data
      if( max(data$time) > 200 ) {
        x_interval <- 24
      } else{
        x_interval <- 12
      }
      max_t <- floor( max(data$time)/x_interval ) * x_interval
      x_breaks <- seq(0, max_t, x_interval)
      
      if( !is.null(attr(obj, "use_boot")) ){
        use_boot <- attr(obj, "use_boot")
      } else {
        use_boot <- if( "boot_adj" %in% names(obj) ) TRUE else FALSE
      }
      
    }
    
    # plot
    plot(obj, conf_int = conf_int, use_boot = use_boot,
         # max_t = max_t, x_breaks = x_breaks,
         steps = F,
         custom_colors = color,
         conf_int_alpha = 0.2,
         gg_theme = km_theme,
         risk_table = risk_table,risk_table_stratify=T, risk_table_ylab = NULL,
         legend.title = "group", legend.position = "bottom") +
      scale_x_continuous(limits = c(0, max_t), breaks = x_breaks, expand = c(0.02,0.02))
    
    
    
  }
  
  .fmt.plot <- function(p, 
                        xlim = NULL, xbre = NULL, xexp = 0.02,
                        ylim = NULL, ybre = NULL, yexp = 0.02,
                        conf_int_bound =T ) {
    
  
    if (!is.null(xlim)) {
      p <- p + scale_x_continuous(limits = xlim, breaks = xbre, expand = expansion(mult = c(xexp, xexp)))
    } 
    
    if (!is.null(ylim)) {
      p <- p + scale_y_continuous(limits = ylim, breaks = ybre, expand = expansion(mult = c(yexp, yexp)))
    } 
    
    
    if(   !is.null(ylim) & isTRUE(conf_int_bound) & ("ci_lower" %in% names(p$data))   ){
      
      p$data <- 
        p$data %>% 
        mutate(ci_lower = if_else(ci_lower < ylim[1], ylim[1], ci_lower ),
               ci_upper = if_else(ci_upper > ylim[2], ylim[2], ci_upper ) )
    }
    
    
    return(p)
  }
  
  
  # plot
  if( class(obj2) %in% "adjustedsurv" ){
    p <-  .plot.sur(obj, color, conf_int, risk_table)
  } else {
    p <-  .plot.cif(obj, color, conf_int)
  }
  
  
  p <- .fmt.plot(p, xlim, xbre, xexp, ylim, ybre, yexp, conf_int_bound =T)
  
  p <- p+ labs(x = xlab)
  
  return(p)
}
.plot.dif <- function(obj, color = lancet_palette[1],
                      xlim = c(0,240),xbre = seq(0,240,24),
                      ylim = NULL, ybre = NULL,
                      group_1=NULL, group_2=NULL){
  
  conf_int=T; conf_level=0.95; use_boot = ifelse("boot_adj" %in% names(obj), TRUE, FALSE)
  if(is.null(xbre)) xbre <- waiver()
  if(is.null(ybre)) ybre <- waiver()
  xexp = 0.02; yexp = 0.02
  
  xlab <- "Time since diagnosis (month)"
  ylab <-  if( obj$method %in% c("km","aalen_johansen") )  "Survival Difference" else  "Adjusted Survival Difference"
  
  
  # max_t=240
  
  # function
  .fmt.plot <- function(p, 
                        xlim = NULL, xbre = NULL, xexp = 0.02,
                        ylim = NULL, ybre = NULL, yexp = 0.02,
                        conf_int_bound =T ) {
    
   
    if (!is.null(xlim)) {
      p <- p + scale_x_continuous(limits = xlim, breaks = xbre, expand = expansion(mult = c(xexp, xexp)))
    } 
    
  
    if (!is.null(ylim)) {
      p <- p + scale_y_continuous(limits = ylim, breaks = ybre, expand = expansion(mult = c(yexp, yexp)))
    } 
    
    
    
    if(   !is.null(ylim) & isTRUE(conf_int_bound) & ("ci_lower" %in% names(p$data))   ){
      
      p$data <- 
        p$data %>% 
        mutate(ci_lower = if_else(ci_lower < ylim[1], ylim[1], ci_lower ),
               ci_upper = if_else(ci_upper > ylim[2], ylim[2], ci_upper ) )
    }
    
    
    return(p)
  }
  
  
  
  p <- plot_curve_diff(obj, 
                       group_1=group_1, group_2=group_2,
                       conf_int=conf_int, conf_level=conf_level, use_boot=use_boot,
                       type="lines", color=color, conf_int_alpha=0.2, line_at_ref_color = "red",
                       gg_theme= km_theme) 
  
  p <- .fmt.plot(p, xlim, xbre, xexp, ylim, ybre, yexp, conf_int_bound =T)
  
  
  p <- p+ labs(x = xlab, y = ylab)
  
  return(p)
  
}


.plot.forest_for.curve <- function(dat_lis,
                                   xlim = NULL, ticks_at = NULL,ticks_digits = NULL, 
                                   core_size = c(2,2)){
  
  
  
  if(  "haz" %in%  attr(dat_lis[[1]], "measure")  ) {
    label_list = list(label = "..Groups..", 
                      eve1 = "Event/total",
                      effect1 = "HR (95%CI)",
                      Effect_sur2 = "Survival probability",
                      effect2 = "Survival difference")
  } else {
    label_list = list(label = "..Groups..", 
                      eve1 = "Event/total",
                      Effect_sur1 = "Survival probability",
                      effect1 = "Survival difference")
  }
  
  
  
  
  
  .plot.forest(dat_lis,xlim = xlim,ticks_at = ticks_at,
               ticks_digits = ticks_digits,
               core_size = core_size,
               fill_col ="pal_div[[3]][7:10]", 
               type = "h",
               forest_position = "before",
               show.sur = F,
               summary_vars = c("label", "eve1"), #"time"
               label_list = label_list
  ) %>% add_fotest_color_KM %>% as_ggplot
  
} 

add_fotest_color_KM <- function(p, fill_alpha = 0.4){
  
  row_ranges <- 10
  fill_colors <- ggsci::pal_lancet("lanonc", alpha = fill_alpha)(9)
  
  for (i in seq_len(row_ranges)) {
    
    p <- edit_plot(p, row = i, which = "background", gp = gpar(fill = fill_colors[i]))
  }
  
  return(p)
  
}


.plt.for <- function(data_list, display = c("h","v","c"), core_size = c(2, 1.5)){
  
  if("label0" %in% names(data_list[[1]]) ){
    summary_vars <- c("label0","label", "eve12") 
  } else {
    # summary_vars <- if("eve1_1" %in% names(data_list[[1]]))  c("label", "eve12", "eve1_1", "eve1_2" ) else c("label", "eve12")
    summary_vars <- c("label", "eve12")
    
  }
  
  
  if(display == "c"){
    p <- 
      .plot.forest(data_list[c(2,4,1,3)],
                   # xlim = c(0,4),
                   # ticks_at = c(0,1,2,4),
                   ticks_digits = NULL,
                   core_size = core_size,
                   fill_col ="pal_div[[5]][3:9]",
                   type = c("h"),
                   forest_position = "before",
                   show.sur = F,
                   summary_vars = summary_vars,
                   label_list = list(label = "Characteristics", 
                                     num1 = "Number(%)",
                                     eve1 = "event/total",
                                     eve12 = "DSS/OC/total")) %>% as_ggplot()
  } else {
    p1 <- 
      .plot.forest(data_list[c(2,4)],
                   # xlim = c(0,4),
                   # ticks_at = c(0,1,2,4),
                   ticks_digits = NULL,
                   core_size = core_size,
                   fill_col ="pal_div[[5]][3:9]",
                   type = c("h"),
                   forest_position = "before",
                   show.sur = F,
                   summary_vars = summary_vars,
                   label_list = list(label = "Characteristics", 
                                     num1 = "Number(%)",
                                     eve1 = "event/total",
                                     eve12 = "DSS/OC/total")) %>% as_ggplot()
    
    p2 <- 
      .plot.forest(data_list[c(1,3)],
                   # xlim = c(0,4),
                   # ticks_at = c(0,1,2,4),
                   ticks_digits = NULL,
                   core_size = core_size,
                   fill_col ="pal_div[[5]][3:9]",
                   type = c("h"),
                   forest_position = "before",
                   show.sur = F,
                   summary_vars = summary_vars,
                   label_list = list(label = "Characteristics", 
                                     num1 = "Number(%)",
                                     eve1 = "event/total",
                                     eve12 = "DSS/OC/total")) %>% as_ggplot()
    
    p <- if(display=="h") p1|p2 else p1/p2
  }
  
  
  return(p)
  
  
}

.plot.forest <- function(data_list, 
                         xlim = NULL,
                         ticks_at = NULL,
                         ticks_digits = NULL,
                         core_size = c(4,3),
                         fill_col ="pal_div[[3]][7:10]", 
                         type = c("h", "v", "i","n"),
                         forest_position = "before",
                         show.sur = F,
                         summary_vars = c("label", "eve1", "time"),
                         label_list = list(label = "Characteristics", 
                                           num1 = "Number(%)",
                                           eve1 = "event/total",
                                           eve12 = "DSS/OC/total")
                         ){
  
  
  
  
  
  # summary_vars = c("label", "eve1", "time")  #"label", "num1", "eve1", "eve12"
  
  if(is.null(ticks_digits)) ticks_digits = map_vec(data_list, function(d) if( "haz" %in% attr(d, "measure") ) {0} else {1})
  nudge_y = if(show.sur) 0.3 else 0
  core = list(padding = unit(core_size, "mm"))
  
  
  data_list[[1]] <- data_list[[1]] %>% mutate(label = str_sub(label, 1, 16)) 
  
  
  
  
  # get forest_vars
  if(show.sur){
    forest_vars = c("effect_two") 
  }else{
    if( any(map_vec(data_list, ~ !is.null( attr(.,"time"))))   ){
      index <- which( map_vec(data_list, ~ !is.null( attr(.,"time"))) )[1]
      forest_vars <- c( attr(data_list[[index]], "group_1"),attr(data_list[[index]], "group_2"), "effect") 
    } else {
      forest_vars = c("effect") 
    }
  }
  
  # get ref_line
  if(show.sur){
    ref_line <- map_vec(data_list, function(d) if( attr(d, "model") == "Cox" ) {1} else {0})
  } else {
    ref_line <- map_vec(data_list, function(d) if( "haz" %in% attr(d, "measure") ) {1} else {0})
  }
  
  # get col_list
  col_list_length <-  map_vec(data_list, function(d){ length( forest_vars[ forest_vars %in% names(d) ] ) + 1 })   
  .get_col_list <- function(start_col, col_list_length) {
    
    parts_list <- vector("list", length(col_list_length))
    
    current_start <- start_col  
    
    for (i in seq_along(col_list_length)) {
      length <- col_list_length[i]  
      current_end <- current_start + length - 1  
      
     
      parts_list[[i]] <- current_start:current_end
      
      
      current_start <- current_end + 1
    }
    
    return(parts_list)
  }
  col_list <- .get_col_list(length(summary_vars)+1, col_list_length)
  
  
  
  
  # function
  .get_summary_forest_data <- function(data_list,summary_vars,forest_vars,forest_position){
    
    .get_summary_data <- function(data, summary_vars) {
      
      if("No. of event/total" %in% summary_vars) {
        data[summary_vars] %>% rename(`event/total` = "No. of event/total")
      } else {
        data[summary_vars]
      }
      
      
    }
    .get_forest_data  <- function(data_list, index, forest_vars, forest_position) {
      dt <- data_list[[index]]
      
      forest_vars <- forest_vars[ forest_vars %in% names(dt) ] 
      
      vars_new_names <- paste0(forest_vars, index)
      forest_name <- paste("......forest", index, "......", sep = "")
      
      if(forest_position == "before"){
        forest_position <- grep("effect", forest_vars) 
      } else {
        forest_position <- grep("effect", forest_vars) + 1 
      }
      
      
      
      
      # forest_position <- 1
      
      
      
      dt_forest <- 
        dt %>%
        dplyr::select(all_of(forest_vars)) %>%
        setNames(vars_new_names) %>%
        mutate(across(everything(), ~replace_na(., ""))) %>%
        
        add_column(!!forest_name := "", .before = forest_position)
      
      return(dt_forest)
    }
    
    # get summary data
    dt_sum <- .get_summary_data(data_list[[1]], summary_vars)
    # get forest data 
    dt_forest <- map_dfc(1:length(data_list), ~ .get_forest_data(data_list, .x, forest_vars, forest_position))
    # get res
    res <- cbind(dt_sum, dt_forest)
    
    return(res)
  }
  .get.ci_list <- function(data_list, show.sur){
    
    if(show.sur && !is.null(attr(data_list[[1]], "group_1")) ){
      
      group_1 <- attr(data_list[[1]], "group_1")  
      group_2 <- attr(data_list[[1]], "group_2")
      estimate_1 <- paste0("estimate", "_", group_1)
      estimate_2 <- paste0("estimate", "_", group_2)
      conf.low_1 <- paste0("conf.low", "_", group_1)
      conf.low_2 <- paste0("conf.low", "_", group_2)
      conf.high_1 <- paste0("conf.high", "_", group_1)
      conf.high_2 <- paste0("conf.high", "_", group_2)
      
      est1 <- lapply(data_list, function(d) d[[estimate_1]])
      est2 <- lapply(data_list, function(d) d[[estimate_2]])
      lower1 <- lapply(data_list, function(d) d[[conf.low_1]])
      lower2 <- lapply(data_list, function(d) d[[conf.low_2]])
      upper1 <- lapply(data_list, function(d) d[[conf.high_1]])
      upper2 <- lapply(data_list, function(d) d[[conf.high_2]])
      
      list( est = c(est1, est2 ),
            lower= c(lower1, lower2),
            upper = c(upper1, upper2) )
      
    } else {
      
      list(est =   lapply(data_list, function(d) d$estimate),
           lower = lapply(data_list, function(d) d$conf.low),
           upper = lapply(data_list, function(d) d$conf.high) )
      
    }
    
    
    
    
  }
  .plot_forest <- function(dt, data_list, ci_list, 
                           xlim, ticks_at, ref_line,ticks_digits, nudge_y, core, 
                           fill_col, type 
  ){
    
    
    .get_forest_plot <- function(dt, ci_list, ci_column, xlim, ticks_at,ref_line,ticks_digits,nudge_y, core){
      
    
      tm <- forest_theme(base_size = 10,
                         # refline_lty = "solid",
                         ci_pch = c(16),
                         ci_col = c("#00468BFF"),
                         footnote_gp = gpar(col = "blue"),
                         # legend_name = attr(dt1, "legend_name"),
                         # legend_value=  attr(dt1, "legend_value"),
                         # vertline_lty = c("dashed", "dotted"),
                         # vertline_col = c("#d6604d", "#bababa"),
                         core = core,
                         ci_alpha = 0.8,
                         ci_lty = 1,
                         ci_lwd = 1.5,
                         ci_Theight = 0.2,
                         refline_col = "red",
                         arrow_type = "closed",
                         footnote_cex = 0.6,
                         footnote_fontface = "italic",
                         footnote_col = "blue")
      
      
      
      
      p <- forest(dt,
                  est = ci_list$est,
                  lower = ci_list$lower, 
                  upper = ci_list$upper,
                  ci_column = ci_column,
                  arrow_lab = c("Low risk", "High risk"),
                  ref_line = ref_line,
                  ticks_digits = ticks_digits,
                  xlim = xlim,
                  ticks_at = ticks_at,
                  x_trans = "none",
                  nudge_y = nudge_y,
                  theme = tm)
      
      
      return(p)
    }
    .add_forest_border <- function(p, Num_col,col_list, Num_forest, header_title, var_row.index){
      
      
      p <- edit_plot(p, row = var_row.index, col = 1, gp = gpar(fontface = "bold"))
      
      
      # divide_columns <- function(total_cols, N_data) {
      #   
      #         #   cols_per_part <- ceiling((total_cols - 2) / N_data)
      #   
      #   parts_list <- vector("list", N_data)
      #   
      #   for (i in 1:N_data) {
      #          #     start_index <- 3 + (i - 1) * cols_per_part
      #     end_index <- min(start_index + cols_per_part - 1, total_cols)
      #     
      #    
      #     parts_list[[i]] <- start_index:end_index
      #   }
      #   
      #   return(parts_list)
      # }
      # col_list <- divide_columns(Num_col, Num_forest)
      p <- insert_text(p, text =  header_title[[1]],
                       part = "header",
                       col = col_list[[1]],
                       gp = gpar(fontface = "bold"))
      
           if (Num_forest > 1) {
                for (i in 2:length(col_list)) {
          p <- p %>% 
            forestploter:: add_text(text = header_title[[i]],
                                    part = "header",
                                    row = 1,
                                    col = col_list[[i]],
                                    gp = gpar(fontface = "bold"))
          
        }
      }
      
      
          for (i in 1:Num_forest) {
        p <- p %>%
          add_border(part = "header", row = 1, col = col_list[[i]][2], gp = gpar(lwd = 1.2)) %>%
          add_border(part = "body", col = col_list[[i]][1], where = "left", gp = gpar(lwd = 1, lty = 4)) %>%
          add_border(part = "header", col = col_list[[i]][1], where = "left", gp = gpar(lwd = 1, lty = 4))
      }
      
          p <- p %>%
        add_border(part = "header", row = c(0, 2), gp = gpar(lwd = 1.5)) %>% 
        add_border(part = "header", row = length(p$rownames) + 1, gp = gpar(lwd = 1.5)) %>%
        add_border(part = "body", col = c(0, Num_col), where = "right", gp = gpar(lwd = 1.5)) %>%
        add_border(part = "header", col = c(0, Num_col), where = "right", gp = gpar(lwd = 1.5))
      
      return(p)
    }
    .add_fotest_color <- function(p, fill_col ="pal_div[[5]][3:9]", 
                                  type = c("h", "v", "i","n"),
                                  col_list, Num_forest, var_row.index,
                                  head = F){
      
      
      N_data <-   Num_forest
      
      
      add_fotest_color_horizontal <- function(p, fill_col, var_row.index){
        
        row_ranges <- var_row.index
        
        fill_col <- colorRampPalette(fill_col)(100) %>% rev()
        
        fill_colors <- fill_col[seq(5,95, length.out = length(row_ranges))]
        
        for (i in seq_along(row_ranges)) {
          if (i < length(row_ranges)) {
            row_start <- row_ranges[i]
            row_end <- row_ranges[i + 1] - 1
          } else {
            row_start <- row_ranges[i]
            row_end <- 100 # Assuming you want to go to the last row for the last range
          }
          p <- edit_plot(p, row = row_start:row_end, which = "background", gp = gpar(fill = fill_colors[i]))
        }
        return(p)
      }
      add_fotest_color_vertical <- function(p, fill_col, N_data, col_list){
        # N_data <- attr(p, "N_data")
        # col_list <-attr(p, "col_list")
        
        fill_col <- colorRampPalette(fill_col)(100) %>% rev()
        fill_colors <- fill_col[seq(5,95, length.out = N_data+1)]
        
        p <- edit_plot(p,col = 1:c(col_list[[1]][1]-1), which = "background", gp = gpar(fill = fill_colors[1]))
        
        for (i in seq_along(col_list)) {
               cols_indices <- col_list[[i]]
          p <- edit_plot(p, col = cols_indices, which = "background", gp = gpar(fill = fill_colors[i+1]))
        }
        p
        
      }
      add_fotest_color_interval <- function(p, fill_col, N_data, col_list){
        
        # N_data <- attr(p, "N_data")
        # col_list <- attr(p, "col_list")
        
        fill_col <- colorRampPalette(fill_col)(100) %>% rev()
        fill_colors <- fill_col[seq(5,95, length.out = N_data+1)]
        
        
        p <- edit_plot(p, row = seq(1,100,2),col = 1:c(col_list[[1]][1]-1), which = "background", gp = gpar(fill = fill_colors[1]))
        
        for (i in seq_along(col_list)) {
          cols_indices <- col_list[[i]]
          
       
          p <- edit_plot(p, row = seq(1, 100, 2), col = cols_indices, which = "background", gp = gpar(fill = fill_colors[i+1]))
        }
        p
      }
      Brewer_color_selction <- function(fill_col) {
        library(RColorBrewer)
        library(dplyr)
        library(purrr)
        
        Brewer_color_extract <- function(x){
          max_colors <- brewer.pal.info[x, "maxcolors"]
          brewer.pal(max_colors, x) 
        }
        
    
        pal_seq_names <- brewer.pal.info[] %>% filter(category == "seq") %>% row.names() %>% rev()
        pal_div_names <- brewer.pal.info[] %>% filter(category == "div") %>% row.names() %>% rev()
        pal_qua_names <- brewer.pal.info[] %>% filter(category == "qual") %>% row.names() %>% rev()
        
    
        pal_seq <- pal_seq_names %>% map(Brewer_color_extract) %>% set_names(pal_seq_names)
        pal_div <- pal_div_names %>% map(Brewer_color_extract) %>% set_names(pal_div_names)
        pal_qua <- pal_qua_names %>% map(Brewer_color_extract) %>% set_names(pal_qua_names)
        
       
        selected_color <- eval(parse(text = fill_col))
        
        print(selected_color %>% prismatic::color())
        
        return(selected_color )
        
      }
      
      if(head){ p <-  edit_plot(p, row = -1:3, part ="header", which = "background", gp = gpar(fill = "grey"))}
      type  <- match.arg(type)
      if( is.null(fill_col) ){ 
        fill_col <- colorRampPalette(RColorBrewer::brewer.pal(3, "RdYlBu"))(100) %>% rev() 
      } else {
        fill_col <- Brewer_color_selction(fill_col) 
      }
      
      switch(type,
             "h" = add_fotest_color_horizontal(p, fill_col, var_row.index),
             "v" = add_fotest_color_vertical(p, fill_col, N_data, col_list),
             "i" = add_fotest_color_interval(p, fill_col, N_data, col_list),
             "n" = p
      )
    }
    
    
    # get forest parameters
    Num_col <- ncol(dt)
    ci_column <- grep("forest", names(dt))
    Num_forest <- length(ci_column)
    
    header_title <- lapply(data_list, function(d) attr(d, "header_title"))
    var_row.index <- which(data_list[[1]]$label %in% attr(data_list[[1]], "data.names"))
    
    # get and fmt plot
    for (var in names(dt)) {
      label <- attr(dt[[var]], "label")  
      if (!is.null(label) && label != "") {  
        names(dt)[names(dt) == var] <- label 
      }
    }
    
    p <- 
      .get_forest_plot(dt, ci_list, ci_column, xlim, ticks_at, ref_line,ticks_digits, nudge_y, core) %>% 
      .add_forest_border(Num_col,col_list, Num_forest, header_title, var_row.index) %>% 
      .add_fotest_color(fill_col, type, col_list, Num_forest, var_row.index)
    
  }
  
  # get plotdata
  dt <- .get_summary_forest_data(data_list, summary_vars, forest_vars, forest_position)
  
  # 修改dt 的label 
  .set_labels <- function(data, label_list) {
    for (var in names(label_list)) {
      if (var %in% names(data)) {
        attr(data[[var]], "label") <- label_list[[var]]
        print(paste(var, "set to:", label_list[[var]]))
      }
    }
    return(data)
  }
  dt <- .set_labels(dt, label_list)
  
  # get ci_list  
  ci_list <- .get.ci_list(data_list, show.sur)
  
  # plot
  p <- .plot_forest(dt, data_list, ci_list, 
                    xlim, ticks_at, ref_line, ticks_digits, nudge_y, core, 
                    fill_col, type)
  
  
  return(p)
}

.insert.plot <- function(p, obj, pos = c(0.1,0.1)){
  # p <-  .plot.km(dat[[1]], "Sur", c("Age", "Sex"))
  # p1 <- .plot.fr(dat[[1]], "Sur",c("Age", "Sex"))
  # tb <- tibble(x = c(0,1), y = c(2,2) )
  # class(p); class(p1);class(tb)
  
  
  if("data.frame" %in% class(obj) ){
    p + 
      geom_table_npc( data = tibble(x = pos[1], y = pos[2], tb = list(obj) ),
                      aes(npcx = x, npcy =y, label = tb )) 
  } else if("ggplot" %in% class(obj) ){
    p + 
      geom_plot_npc(data = tibble(x = pos[1], y = pos[2], plot = list(obj) ),
                    aes(npcx = x, npcy =y, label = plot ))
  } else if("grob" %in% class(obj) ) {
    p + 
      geom_grob_npc(data = tibble(x = pos[1], y = pos[2], grob = list(obj) ),
                    aes(npcx = x, npcy =y, label = grob ))
  }
  
} 


res.cox1_as_rcs <- function(res, boot_curves = F, my_conf = T,conf.lim= c(0.15,0.9),
                            lin_color=lancet_palette[1], ylim= c(0,10),
                            histlimit = c(0,5), his_color = lancet_palette[3]){
  
  
  raw_dat <- attr(res,"raw_dat")
  data <- raw_dat$data
  var1 <- raw_dat$var1
  xlab <- var2 <- raw_dat$var2
  co_var <- raw_dat$co_var
  xlim <- raw_dat$xlim
  
  
  yexpand = c(0.01,0.01)
  xexpand = c(0.02,0.02)
  # legend.pos = c(0.2,0.8)
  pvalue.position = c(0.3, 0.95); p.string <- attr(res,"raw_dat")$p.string
  
  
  sec_scale <- c(data %>% count(.data[[var2]]) %>% pull(n) %>% max)/histlimit[2]
  
  
  if(is.null(co_var)){
    ylab <-  sprintf("Unadjusted HR of %s", var1)
  }else{
    ylab <-  sprintf("Adjusted HR of %s", var1)
  }
  
  his_color = lancet_palette[3]
  
  
  
  ## boot curves
  plotdata1 <-
    attr(res,"raw_dat")$myBoot$t %>% t() %>%  as.data.frame() %>% mutate(x = res$x, .before = 1) %>%
    pivot_longer(cols = -1 , names_to = "boot", values_to = "yhat") %>%  arrange(boot)
  
  ## observed curves
  if(isTRUE(my_conf)){
    
    plotdata2 <- plotdata1 %>%
      group_by(x) %>%
      dplyr::summarize(lower = quantile(yhat, probs = conf.lim[1], na.rm = TRUE),
                       upper = quantile(yhat, probs = conf.lim[2], na.rm = TRUE),
                       # SE    = sd(yhat, na.rm = TRUE) / sqrt(n()),
                       # yhat = mean(yhat, na.rm = TRUE),
                       yhat = median(yhat, na.rm = TRUE)
      ) %>%
      mutate(upper = if_else(upper > ylim[2], ylim[2], upper))
    
  }else{
    
    plotdata2 <- res %>% mutate(upper = if_else(upper > ylim[2], ylim[2], upper))
  }
  
  
  if(isTRUE(boot_curves)){
    plot <-
      ggplot(plotdata1,aes(x=x, y=yhat)) +
      # geom_line(aes(group=boot), color="grey",size=0.01,  alpha=0.8) +
      geom_line(aes(group=boot), color=lin_color, size=0.001,  alpha=0.03) +
      geom_line(data = plotdata2, aes(x=x, y=yhat),color = lin_color)+
      geom_ribbon(data = plotdata2, aes(ymin = lower, ymax = upper), fill=lin_color, alpha = 0.3)
    
  }else{
    plot <-
      ggplot(plotdata2,aes(x=x, y=yhat)) +
      geom_line(color = lin_color)+
      geom_ribbon(aes(ymin = lower, ymax = upper), fill=lin_color, alpha = 0.3)
    
  }
  plot <- plot+
    geom_histogram(data = data,
                   aes(.data[[var2]], y = scales::rescale(after_stat(count),histlimit)), fill = his_color, col = his_color,
                   binwidth = 1, alpha = 0.1, position = "stack")
  
  
  plot <- plot +
    km_theme + #theme(legend.position = legend.pos)+
    geom_hline(yintercept = 1, colour =  "black", linetype = 2, size = 1) +
    labs(x = xlab, y = ylab)+
    scale_y_continuous(limits = ylim, sec.axis = sec_axis(~.*sec_scale, name = "Count (n)"), expand = yexpand) +
    scale_x_continuous(limits = xlim, expand = xexpand)
  
  
  ## show p value
  px <-  min(xlim) + (max(xlim) - min(xlim)) * pvalue.position[1]
  py <-  min(ylim) + (max(ylim) - min(ylim)) * pvalue.position[2]
  
  plot <- plot + draw_label(p.string,
                            size = 12,
                            x = px,
                            y = py,
                            hjust = 0,
                            vjust = 1)
  
  
  attr(plot,"plot_dat") <- list(data = data, plotdata1 = plotdata1, plotdata2 = plotdata2)
  attr(plot,"raw_dat") <- attr(res,"raw_dat")
  
  
  
  plot
  
}


res.poi_as_rcs <- function(res, boot_curves = F,my_conf = T,conf.lim= c(0.15,0.9),
                           ylim= c(0,10), yexpand = c(0.01,0.01), histlimit = c(0,5)){
  
  
  raw_dat <- attr(res,"raw_dat")
  data <- raw_dat$model$data
  var1 <- raw_dat$var1
  xlab <- var2 <- raw_dat$var2
  co_var <- raw_dat$co_var
  PA <- raw_dat$PA
  xlim <- raw_dat$xlim
  
  
  yexpand = yexpand
  xexpand = c(0.02,0.02)
  legend.pos = c(0.2,0.8)
  pvalue.position = c(0.3, 0.95); p.string <- attr(res,"raw_dat")$p.string
  
  
  if(is.null(var1)){
    sec_scale <- c(data %>% count(.data$x) %>% pull(n) %>% max)/histlimit[2]
  }else{
    sec_scale <- c(data %>% count(.data$x, .data[[var1]]) %>% pull(n) %>% max)/histlimit[2]
  }
  
  if(is.null(co_var)){
    ylab <- sprintf("Unadjusted death rate\n(per %.0f patient-years)", PA/12)
  }else{
    ylab <- sprintf("Adjusted death rate\n(per %.0f patient-years)", PA/12)
  }
  
  sin.his_color = lancet_palette[3]
  
  
  
  ## boot curves
  plotdata1 <-
    attr(res,"raw_dat")$myBoot$t %>% t() %>%  as.data.frame() %>% mutate(x = res$x, .before = 1)
  
  if(is.null(var1)){
    plotdata1 <- plotdata1 %>%
      pivot_longer(cols = -1 , names_to = "boot", values_to = "yhat") %>%  arrange(boot)
  }else{
    plotdata1 <- plotdata1 %>%  mutate(group = res$group, .before = 2) %>%
      pivot_longer(cols = -c(1,2) , names_to = "boot", values_to = "yhat") %>%
      mutate(line_group = paste0(group,boot)) %>%  arrange(boot)
  }
  
  ## observed curves
  if(isTRUE(my_conf)){
    if(is.null(var1)){
      plotdata2 <- plotdata1 %>%
        group_by(x) %>%
        dplyr::summarize(lower = quantile(yhat, probs = conf.lim[1], na.rm = TRUE),
                         upper = quantile(yhat, probs = conf.lim[2], na.rm = TRUE),
                         SE    = sd(yhat, na.rm = TRUE) / sqrt(n()),
                         yhat = median(yhat, na.rm = TRUE)) %>%
        mutate(upper = if_else(upper > ylim[2], ylim[2], upper))
      
    }else{
      plotdata2 <- plotdata1 %>%
        group_by(x, group) %>%
        dplyr::summarize(lower = quantile(yhat, probs = conf.lim[1], na.rm = TRUE),
                         upper = quantile(yhat, probs = conf.lim[2], na.rm = TRUE),
                         SE    = sd(yhat, na.rm = TRUE) / sqrt(n()),
                         yhat = median(yhat, na.rm = TRUE)) %>%
        mutate(upper = if_else(upper > ylim[2], ylim[2], upper))
    }
    
    
  }else{
    
    plotdata2 <- res %>% mutate(upper = if_else(upper > ylim[2], ylim[2], upper))
    
  }
  
  
  
  if(is.null(var1)){
    if(isTRUE(boot_curves)){
      plot <-
        ggplot(plotdata1,aes(x=x, y=yhat)) +
        geom_line(aes(group=boot), color="grey",size=0.01,  alpha=0.8) +
        geom_line(data = plotdata2, aes(x=x, y=yhat),color = sin.his_color)+
        geom_ribbon(data = plotdata2, aes(ymin = lower, ymax = upper), fill=sin.his_color, alpha = 0.3)
      
    }else{
      plot <-
        ggplot(plotdata2,aes(x=x, y=yhat)) +
        geom_line(color = color)+
        geom_ribbon(aes(ymin = lower, ymax = upper), fill=color, alpha = 0.3)
      
    }
    plot <- plot+
      geom_histogram(data = data,
                     aes(x, y = scales::rescale(after_stat(count),histlimit)), fill = color, col = color,
                     binwidth = 1, alpha = 0.1, position = "stack")
    
  } else {
    if(isTRUE(boot_curves)){
      plot <-
        ggplot(plotdata1,aes(x=x, y=yhat)) +
        geom_line(aes(group=line_group, color = group), size=0.001,  alpha=0.03) +
        geom_line(data = plotdata2, aes(x=x, y=yhat, color = group))+
        geom_ribbon(data = plotdata2, aes(ymin = lower, ymax = upper, fill=group ),  alpha = 0.3)
      
    }else{
      plot <-
        ggplot(plotdata2,aes(x=x, y=yhat)) +
        # geom_line(aes(group=boot), color="grey",size=0.01,  alpha=0.8) +
        geom_line(aes(color = group))+
        geom_ribbon(aes(ymin = lower, ymax = upper, fill=group), alpha = 0.3)
    }
    
    
    plot <- plot+
      geom_histogram(data = data %>% rename(group = var1) ,
                     aes(x, y = scales::rescale(after_stat(count),histlimit), fill = group, col = group),
                     binwidth = 1, alpha = 0.1, position = "stack")+
      scale_color_manual(values = lancet_palette)+
      scale_fill_manual(values = lancet_palette)
    
  }
  
  
  
  plot <- plot +
    km_theme + theme(legend.position = legend.pos)+
    geom_hline(yintercept=0,size=1, linetype="dashed", color = "black") +
    labs(x = xlab, y = ylab)+
    scale_y_continuous(limits = ylim, sec.axis = sec_axis(~.*sec_scale, name = "Count (n)"), expand = yexpand) +
    scale_x_continuous(limits = xlim, expand = xexpand)
  
  
  
  ## show p value
  px <-  min(xlim) + (max(xlim) - min(xlim)) * pvalue.position[1]
  py <-  min(ylim) + (max(ylim) - min(ylim)) * pvalue.position[2]
  
  plot <- plot + draw_label(p.string,
                            size = 12,
                            x = px,
                            y = py,
                            hjust = 0,
                            vjust = 1)
  
  
  
  attr(plot,"plot_dat") <- list(data = data, plotdata1 = plotdata1, plotdata2 = plotdata2)
  attr(plot,"raw_dat") <- attr(res,"raw_dat")
  
  plot
  
}


## Get Table of Interaction Analysis
.get.interTable <- function(data, exposure_names, adj_var, surv = T){
  
  outcome <- if(surv) "DSS" else "status"
  
  ci.type = "delta"
  ci.level = 0.95
  em = F
  recode = F
  p.value = T
  
  ## function get table
  Inter_table <- function(obj, p.value = FALSE, file_path = NA) {
    if (!is(obj, "interactionR")) {
      stop("Argument 'obj' must be an object of class 'interactionR',
             use the interactionR() function to generate such object ")
    }
    
    beta1 <- obj$exp_names[1]
    beta2 <- obj$exp_names[2]
    em <- obj$analysis
    d <- obj$dframe
    d$Estimates <- as.character(round(d$Estimates, 2))
    d$CI.ll <- as.character(round(d$CI.ll, 2))
    d$CI.ul <- as.character(round(d$CI.ul, 2))
    E1.absent <- paste(beta1, "absent", sep = " ")
    E1.present <- paste(beta1, "present", sep = " ")
    E2.absent <- paste(beta2, "absent", sep = " ")
    E2.present <- paste(beta2, "present", sep = " ")
    WithinStrataEffect1 <- paste("Effect of", beta2,   sep = " ")
    WithinStrataEffect2 <- paste("Effect of", beta1, sep = " "
    )
    
    if (grepl("glm", obj$call[1])) {
      effect_measure <- "OR [95% CI]"
    } else {
      effect_measure <- "HR [95% CI]"
    }
    
    
    # effect_measure <- "HR [95% CI]"
    
    
    if (p.value) {
      if (em) {
        t <- data.frame(c(
          NA, NA, E1.absent, E1.present, "Multiplicative scale",
          "RERI"
        ), c(NA, effect_measure, NA, NA, NA, NA),
        c(NA, effect_measure, NA, NA, NA, NA), c(
          NA,
          effect_measure, NA, NA, NA, NA
        ),
        stringsAsFactors = FALSE
        )
        names(t) <- c("*", E2.absent, E2.present, WithinStrataEffect1)
        
        t[3, 2] <- paste("1", "[Reference]", sep = " ")
        t[3, 3] <- paste(d[2, 2], " [", d[2, 3], ", ", d[2, 4], "]", " p=", d[2, 5], sep = "")
        t[3, 4] <- paste(d[5, 2], " [", d[5, 3], ", ", d[5, 4], "]", " p=", d[5, 5], sep = "")
        t[4, 2] <- paste(d[3, 2], " [", d[3, 3], ", ", d[3, 4], "]", " p=", d[3, 5], sep = "")
        t[4, 3] <- paste(d[4, 2], " [", d[4, 3], ", ", d[4, 4], "]", " p=", d[4, 5], sep = "")
        t[4, 4] <- paste(d[6, 2], " [", d[6, 3], ", ", d[6, 4], "]", " p=", d[6, 5], sep = "")
        t[5, 2] <- paste(d[7, 2], " [", d[7, 3], ", ", d[7, 4], "]", " p=", d[7, 5], sep = "")
        t[6, 2] <- paste(d[8, 2], " [", d[8, 3], ", ", d[8, 4], "]", " p=", d[8, 5], sep = "")
        
        
        t2 <- flextable::flextable(t)
        t2 <- flextable::set_caption(t2, paste("Modification of the effect of", beta1, "and", beta2, sep = " "))
      } else {
        t <- data.frame(c(
          NA, NA, E1.absent, E1.present, WithinStrataEffect2,
          "Multiplicative scale", "RERI", "AP", "SI"
        ), c(
          NA,
          effect_measure, NA, NA, NA, NA, NA, NA, NA
        ), c(
          NA,
          effect_measure, NA, NA, NA, NA, NA, NA, NA
        ), c(
          NA,
          effect_measure, NA, NA, NA, NA, NA, NA, NA
        ), stringsAsFactors = FALSE)
        
        names(t) <- c("*", E2.absent, E2.present, WithinStrataEffect1)
        
        t[3, 2] <- paste("1", "[Reference]", sep = " ")
        t[3, 3] <- paste(d[2, 2], " [", d[2, 3], ", ", d[2, 4], "]",  d[2, 5], sep = "")
        t[3, 4] <- paste(d[5, 2], " [", d[5, 3], ", ", d[5, 4], "]",  d[5, 5], sep = "")
        t[4, 2] <- paste(d[3, 2], " [", d[3, 3], ", ", d[3, 4], "]",  d[3, 5], sep = "")
        t[4, 3] <- paste(d[4, 2], " [", d[4, 3], ", ", d[4, 4], "]",  d[4, 5], sep = "")
        t[4, 4] <- paste(d[6, 2], " [", d[6, 3], ", ", d[6, 4], "]",  d[6, 5], sep = "")
        t[5, 2] <- paste(d[7, 2], " [", d[7, 3], ", ", d[7, 4], "]",  d[7, 5], sep = "")
        t[5, 3] <- paste(d[8, 2], " [", d[8, 3], ", ", d[8, 4], "]",  d[8, 5], sep = "")
        t[6, 2] <- paste(d[9, 2], " [", d[9, 3], ", ", d[9, 4], "]",  d[9, 5], sep = "")
        t[7, 2] <- paste(d[10, 2], " [", d[10, 3], ", ", d[10, 4], "]",  d[10, 5], sep = "")
        t[8, 2] <- paste(d[11, 2], " [", d[11, 3], ", ", d[11, 4], "]",  d[11, 5], sep = "")
        t[9, 2] <- paste(d[12, 2], " [", d[12, 3], ", ", d[12, 4], "]",  d[12, 5], sep = "")
        t2 <- flextable::flextable(t)
        t2 <- flextable::set_caption(t2, paste("Interaction of", beta1, "and", beta2,"in", obj$model,"model", sep = " "))
      }
    } else {
      if (em) {
        t <- data.frame(c(
          NA, NA, E1.absent, E1.present, "Multiplicative scale",
          "RERI"
        ), c(NA, effect_measure, NA, NA, NA, NA),
        c(NA, effect_measure, NA, NA, NA, NA), c(
          NA,
          effect_measure, NA, NA, NA, NA
        ),
        stringsAsFactors = FALSE
        )
        names(t) <- c("*", E2.absent, E2.present, WithinStrataEffect1)
        
        t[3, 2] <- paste("1", "[Reference]", sep = " ")
        t[3, 3] <- paste(d[2, 2], " [", d[2, 3], ", ", d[2, 4], "]", sep = "")
        t[3, 4] <- paste(d[5, 2], " [", d[5, 3], ", ", d[5, 4], "]", sep = "")
        t[4, 2] <- paste(d[3, 2], " [", d[3, 3], ", ", d[3, 4], "]", sep = "")
        t[4, 3] <- paste(d[4, 2], " [", d[4, 3], ", ", d[4, 4], "]", sep = "")
        t[4, 4] <- paste(d[6, 2], " [", d[6, 3], ", ", d[6, 4], "]", sep = "")
        t[5, 2] <- paste(d[7, 2], " [", d[7, 3], ", ", d[7, 4], "]", sep = "")
        t[6, 2] <- paste(d[8, 2], " [", d[8, 3], ", ", d[8, 4], "]", sep = "")
        
        
        t2 <- flextable::flextable(t)
        t2 <- flextable::set_caption(t2, paste("Modification of the effect of", beta1, "and", beta2, sep = " "))
      } else {
        
        t <- data.frame(c(
          NA, NA, E1.absent, E1.present, WithinStrataEffect2,
          "MS", "RERI", "AP", "SI"
        ), c(
          NA,
          effect_measure, NA, NA, NA, NA, NA, NA, NA
        ), c(
          NA,
          effect_measure, NA, NA, NA, NA, NA, NA, NA
        ), c(
          NA,
          effect_measure, NA, NA, NA, NA, NA, NA, NA
        ), stringsAsFactors = FALSE)
        
        names(t) <- c("*", E2.absent, E2.present, WithinStrataEffect1)
        
        t[3, 2] <- paste("1", "[Reference]", sep = " ")
        t[3, 3] <- paste(d[2, 2], " [", d[2, 3], ", ", d[2, 4], "]", sep = "")
        t[3, 4] <- paste(d[5, 2], " [", d[5, 3], ", ", d[5, 4], "]", sep = "")
        t[4, 2] <- paste(d[3, 2], " [", d[3, 3], ", ", d[3, 4], "]", sep = "")
        t[4, 3] <- paste(d[4, 2], " [", d[4, 3], ", ", d[4, 4], "]", sep = "")
        t[4, 4] <- paste(d[6, 2], " [", d[6, 3], ", ", d[6, 4], "]", sep = "")
        t[5, 2] <- paste(d[7, 2], " [", d[7, 3], ", ", d[7, 4], "]", sep = "")
        t[5, 3] <- paste(d[8, 2], " [", d[8, 3], ", ", d[8, 4], "]", sep = "")
        t[6, 2] <- paste(d[9, 2], " [", d[9, 3], ", ", d[9, 4], "]", sep = "")
        t[7, 2] <- paste(d[10, 2], " [", d[10, 3], ", ", d[10, 4], "]",
                         sep = ""
        )
        t[8, 2] <- paste(d[11, 2], " [", d[11, 3], ", ", d[11, 4], "]",
                         sep = ""
        )
        t[9, 2] <- paste(d[12, 2], " [", d[12, 3], ", ", d[12, 4], "]",
                         sep = ""
        )
        
        # t2 <- flextable::flextable(t)
        # t2 <- flextable::set_caption(t2, paste("Interaction of", beta1, "and", beta2,"in", obj$model,"model", sep = " "))
      }
      
    }
    
    
    t
    
    
    
  }
  
  
  ## check adj_var
  adj_var <- .check.adj_var(adj_var, exposure_names)
  # ## check exposure_names
  # data <- .check.exposure_names(data, exposure_names, adj_var)
  
  # get model
  inter_var <- paste(exposure_names, collapse = "*") 
  model <- .get.fit(data,c(inter_var, adj_var), time = "time", outcome = outcome)
  
  if(is.null(adj_var)){
    attr(model, "model") <- "Uni_Cox.DSS"
  }else {
    attr(model, "model") <- "Mul_Cox.DSS"
  }
  
  
  # Estimates the critical value from the supplied CI.level for
  # subsequent CI estimations
  alpha <- 1 - ci.level
  z <- qnorm(1 - alpha / 2)
  
  # Extracts the names for the main exposure (beta1), the effect modifier
  # (beta2) and their interaction term
  e1 <- grep(exposure_names[1], names(coef(model)), value = TRUE, ignore.case = TRUE)
  e2 <- grep(exposure_names[2], names(coef(model)), value = TRUE, ignore.case = TRUE)
  e1_e2 <- intersect(e1,e2)
  if (length(e1_e2) != 1) {
    stop("The interaction you specified in your exposure_names argument cannot be found in the model")
  }
  beta1 <- e1[1]
  beta2 <- e2[1]
  beta3 <- e1_e2[1]
  
  varNames <- c(beta1, beta2, beta3)
  
  # estimating coefficients to check for any preventive exposures
  b1 <- coef(model)[beta1]
  b2 <- coef(model)[beta2]
  b3 <- coef(model)[beta3]
  
  #### Recode section code is adapted from Marthur and Vanderweele 2018 (doi: 10.1097/EDE.0000000000000752)
  
  # check if any exposure is preventive
  preventive <- function(OR10, OR01) {
    return((OR10 < 1) || (OR01 < 1))
  }
  if (preventive(OR10 = exp(b1), OR01 = exp(b2))) {
    if (!recode) {
      warning("At least one exposure is preventive. Set argument recode=TRUE for the exposures to be automatically recoded. see Knol et al. (2011) European Journal of Epidemiology, 26(6), 433-438")
    }
    if (recode) {
      if ("coxph" %in% class(model)) {
        stop("Currently, interactionR() cannot automatically recode models fitted with coxph or clogit. Recode your exposure variables following the examples in Knol et al. (2011) European Journal of Epidemiology, 26(6), 433-438, re-fit your model, and re-run interactionR()")
      }
      temp <- data.frame(cat = c("OR10", "OR01", "OR11"), value = c(
        exp(b1),
        exp(b2), exp(b1 + b2 + b3)
      ))
      ref.cat <- temp$cat[which.min(temp$value)]
      
      E1.ref <- substr(ref.cat, 3, 3)
      E2.ref <- substr(ref.cat, 4, 4)
      
      # extract the raw data that was used to fit the model
      dat.ir <- model.frame(model)
      
      # recode based on new reference category
      dat.ir[[beta1]] <- ifelse(dat.ir[[beta1]] == E1.ref, 0, 1)
      dat.ir[[beta2]] <- ifelse(dat.ir[[beta2]] == E2.ref, 0, 1)
      
      # inform the user
      warning(
        "Recoding exposures; new reference category for ",
        beta1, " is ", E1.ref, " and for ", beta2, " is ", E2.ref
      )
      
      # refit model with user's original call
      model <- update(model, . ~ ., data = dat.ir)
      
      # get new coefficients and ORs
      b1 <- coef(model)[beta1]
      b2 <- coef(model)[beta2]
      b3 <- coef(model)[beta3]
    }
  }
  
  
  
  if( "tidycrr" %in% class(model)){
    data <- model[["tidy"]] %>% column_to_rownames(var = "term")
  } else { data <- tidy(model) %>% column_to_rownames(var = "term") %>%
    cbind(confint.default(model)) %>%
    rename(conf.low = `2.5 %`, conf.high = `97.5 %`)
  }
  
  
  
  
  v1 <- data[beta1,"std.error"]^2; names(v1) = beta1
  v2 <- data[beta2,"std.error"]^2; names(v2) = beta2
  v3 <- data[beta3,"std.error"]^2; names(v3) = beta3
  
  #Extracts p-values from the model
  pvals <- data[, "p.value"]; names(pvals) = rownames(data)
  
  ### Extracts the variance-covariance matrix from the model### for use in
  ### the delta and MOVER method CI estimation for RERI and AP###
  v_cov <- vcov(model)
  rownames(v_cov) <- rownames(data)
  colnames(v_cov) <- rownames(data)
  
  v_cov1 <- v_cov[varNames, varNames] # for deltamethod
  cov12 <- v_cov[beta1, beta2]
  cov13 <- v_cov[beta1, beta3]
  cov23 <- v_cov[beta2, beta3]
  v123 <- v1 + v2 + v3 + (2 * (cov12 + cov13 + cov23))
  v12 <- v1 + v2 + (2 * (cov12))
  v13 <- v1 + v3 + (2 * (cov13))
  v23 <- v2 + v3 + (2 * (cov23))
  
  # Estimates individual and joint effects ORs (with CI) from the model
  OR00 <- 1 # reference OR
  OR10 <- as.numeric(exp(b1))
  l1 <- exp(data[beta1, "conf.low"])
  u1 <- exp(data[beta1, "conf.high"]) # This is also OR of X on D (A==0)
  p.OR10 <- pvals[beta1]
  OR01 <- as.numeric(exp(b2))
  l2 <- exp(data[beta2, "conf.low"])
  u2 <- exp(data[beta2, "conf.high"]) # This is also OR of A on D (X==0)
  p.OR01 <- pvals[beta2]
  OR11 <- as.numeric(exp(b1 + b2 + b3))
  l3 <- exp(b1 + b2 + b3 - z * sqrt(v123))
  u3 <- exp(b1 + b2 + b3 + z * sqrt(v123))
  q1 <- abs(log(OR11)/sqrt(v123))
  p.OR11 <- exp(-0.717*q1 - 0.416*q1^2) #see BMJ 2011;343:d2304
  
  
  ### Estimates the effect (and CI) of A on D (X==1) ###
  OR_X1 <- as.numeric(exp(b2 + b3)) # OR of A on D (X==1)
  CI.ll_OR_X1 <- exp(b2 + b3 - z * sqrt(v23))
  CI.ul_OR_X1 <- exp(b2 + b3 + z * sqrt(v23))
  q2 <- abs(log(OR_X1)/sqrt(v23))
  p.OR_X1 <- exp(-0.717*q2 - 0.416*q2^2)
  
  ### Estimates the effect (and CI) of X on D (A==1) ###
  OR_A1 <- as.numeric(exp(b1 + b3)) # OR of X on D (A==1)
  CI.ll_OR_A1 <- exp(b1 + b3 - z * sqrt(v13))
  CI.ul_OR_A1 <- exp(b1 + b3 + z * sqrt(v13))
  q3 <- abs(log(OR_A1)/sqrt(v13))
  p.OR_A1 <- exp(-0.717*q3 - 0.416*q3^2)
  
  # Effect modification on the multiplicative scale and CI
  OR_M <- as.numeric(exp(b3))
  CI.ll_OR_M <- exp(data[beta3, "conf.low"])
  CI.ul_OR_M <- exp(data[beta3, "conf.high"])
  p.OR_M <- pvals[beta3]
  
  
  
  # Estimates measures of effect modification on the additive scale
  if (ci.type == "mover") {
    # Estimates measures of effect modification on the additive scale and
    # calculates their CI with the 'MOVER' method (Zou (2018)
    # https://doi.org/10.1093/aje/kwn104)
    
    # RERI, CI and p-value
    RERI <- OR11 - OR01 - OR10 + 1
    r12 <- (v1 + cov12 + cov13) / sqrt(v1 * v123)
    r13 <- (cov12 + v2 + cov23) / sqrt(v2 * v123)
    r23 <- cov12 / sqrt(v1 * v2)
    
    p1 <- (OR11 - l3)^2 + (u1 - OR10)^2 + (u2 - OR01)^2
    p2 <- 2 * r12 * (OR11 - l3) * (u1 - OR10)
    p3 <- 2 * r13 * (OR11 - l3) * (u2 - OR01)
    p4 <- 2 * r23 * (u1 - OR10) * (u2 - OR01)
    p5 <- p1 - p2 - p3 + p4
    p6 <- p5^0.5
    
    L <- 1 + OR11 - OR10 - OR01 - p6
    
    k1 <- (u3 - OR11)^2 + (OR10 - l1)^2 + (OR01 - l2)^2
    k2 <- 2 * r12 * (u3 - OR11) * (OR10 - l1)
    k3 <- 2 * r13 * (u3 - OR11) * (OR01 - l2)
    k4 <- 2 * r23 * (OR10 - l1) * (OR01 - l2)
    k5 <- (k1 - k2 - k3 + k4)^0.5
    
    U <- 1 + OR11 - OR10 - OR01 + k5
    p.RERI <- NA
    
    
    # AP, CI and p-value
    theta1 <- 1 / exp(b1 + b2 + b3)
    theta2 <- 1 / exp(b2 + b3)
    theta3 <- 1 / exp(b1 + b3)
    AP <- theta1 - theta2 - theta3 + 1
    APr12 <- (cov12 + cov13 + v2 + (2 * cov23) + v3) / sqrt(v23 * v123)
    APr13 <- (v1 + cov12 + (2 * cov13) + cov23 + v3) / sqrt(v13 * v123)
    APr23 <- (cov12 + cov23 + cov13 + v3) / sqrt(v23 * v13)
    
    APl1 <- theta1 * exp(-z * sqrt(v123))
    APu1 <- theta1 * exp(z * sqrt(v123))
    
    APl2 <- theta2 * exp(-z * sqrt(v23))
    APu2 <- theta2 * exp(z * sqrt(v23))
    
    APl3 <- theta3 * exp(-z * sqrt(v13))
    APu3 <- theta3 * exp(z * sqrt(v13))
    
    APp1 <- (theta1 - APl1)^2 + (APu2 - theta2)^2 + (APu3 - theta3)^2
    APp2 <- 2 * APr12 * (theta1 - APl1) * (APu2 - theta2)
    APp3 <- 2 * APr13 * (theta1 - APl1) * (APu3 - theta3)
    APp4 <- 2 * APr23 * (APu2 - theta2) * (APu3 - theta3)
    APp5 <- APp1 - APp2 - APp3 + APp4
    APp6 <- APp5^0.5
    
    APL <- 1 + theta1 - theta2 - theta3 - APp6
    
    APk1 <- (APu1 - theta1)^2 + (theta2 - APl2)^2 + (theta3 - APl3)^2
    APk2 <- 2 * APr12 * (APu1 - theta1) * (theta2 - APl2)
    APk3 <- 2 * APr13 * (APu1 - theta1) * (theta3 - APl3)
    APk4 <- 2 * APr23 * (theta2 - APl2) * (theta3 - APl3)
    APk5 <- (APk1 - APk2 - APk3 + APk4)^0.5
    
    APU <- 1 + theta1 - theta2 - theta3 + APk5
    p.AP <- NA
    
    
    # SI, CI and p-value
    SItheta1 <- log((exp(b1 + b2 + b3) - 1))
    SItheta2 <- log((exp(b1) + exp(b2) - 2))
    lnSI <- SItheta1 - SItheta2
    SI <- exp(lnSI)
    
    vSItheta1 <- (exp(b1 + b2 + b3) / (exp(b1 + b2 + b3) - 1))^2 * v123
    vSItheta2 <- ((exp(2 * b1) * v1) + (exp(2 * b2) * v2) + (2 * exp(b1 +
                                                                       b2) * cov12)) / (exp(b1) + exp(b2) - 2)^2
    SIl1 <- SItheta1 - z * sqrt(vSItheta1)
    SIu1 <- SItheta1 + z * sqrt(vSItheta1)
    SIl2 <- SItheta2 - z * sqrt(vSItheta2)
    SIu2 <- SItheta2 + z * sqrt(vSItheta2)
    
    SIr <- ((exp(b1) * (v1 + cov12 + cov13)) + (exp(b2) * (cov12 + v2 +
                                                             cov23))) / sqrt(v123 * ((exp(2 * b1) * v1) + (exp(2 * b2) * v2) +
                                                                                       (2 * exp(b1 + b2) * cov12)))
    
    lnSIL <- (SItheta1 + (-SItheta2)) - sqrt((SItheta1 - SIl1)^2 + ((-SItheta2) -
                                                                      (-SIl2))^2 + (2 * SIr * (SItheta1 - SIl1) * ((-SItheta2) -
                                                                                                                     (-SIl2))))
    lnSIU <- (SItheta1 + (-SItheta2)) + sqrt((SIu1 - SItheta1)^2 + ((-SIu2) -
                                                                      (-SItheta2))^2 + (2 * SIr * (SIu1 - SItheta1) * ((-SIu2) -
                                                                                                                         (-SItheta2))))
    SIL <- exp(lnSIL)
    SIU <- exp(lnSIU)
    p.SI <- NA
    
  } else if (ci.type == "delta") {
    # Estimates measures of effect modification on the additive scale and
    # calculates their CI and p-value with the delta method implemented in
    # the msm package
    
    # RERI, CI and p-value
    RERI <- OR11 - OR01 - OR10 + 1
    se_RERI <- msm::deltamethod(g = ~ exp(x1 + x2 + x3) - exp(x1) - exp(x2) +
                                  1, mean = c(b1, b2, b3), cov = v_cov1)
    L <- RERI - z * se_RERI
    U <- RERI + z * se_RERI
    p.RERI <- 1 - pnorm(RERI/se_RERI)
    
    
    # AP, CI and p-value
    AP <- RERI / OR11
    se_AP <- msm::deltamethod(g = ~ (exp(x1 + x2 + x3) - exp(x1) - exp(x2) +
                                       1) / exp(x1 + x2 + x3), mean = c(b1, b2, b3), cov = v_cov1)
    APL <- AP - z * se_AP
    APU <- AP + z * se_AP
    p.AP <- 1 - pnorm(abs(AP)/se_AP)
    
    
    
    # SI, CI and p-value
    lnSI <- log((exp(b1 + b2 + b3) - 1)) - log((exp(b1) + exp(b2) -
                                                  2))
    SI <- exp(lnSI)
    se_SI <- msm::deltamethod(g = ~ log((exp(x1 + x2 + x3) - 1)) - log((exp(x1) +
                                                                          exp(x2) - 2)), mean = c(b1, b2, b3), cov = v_cov1)
    
    SIL <- exp(lnSI - z * se_SI)
    SIU <- exp(lnSI + z * se_SI)
    p.SI <- 1 - plnorm(exp(lnSI/se_SI))
    
  } else {
    stop("Argument 'ci.type' must be 'delta' or 'mover' ")
  }
  
  
  
  
  d <- data.frame(
    Measures = c(
      "OR00", "OR01", "OR10", "OR11", paste("OR(",
                                            beta2, " on outcome [", beta1, "==0]",
                                            sep = ""
      ), paste("OR(",
               beta2, " on outcome [", beta1, "==1]",
               sep = ""
      ), "Multiplicative scale",
      "RERI"
    ), Estimates = c(
      OR00, OR01, OR10, OR11, OR01, OR_X1, OR_M,
      RERI
    ), CI.ll = c(NA, l2, l1, l3, l2, CI.ll_OR_X1, CI.ll_OR_M, L),
    CI.ul = c(NA, u2, u1, u3, u2, CI.ul_OR_X1, CI.ul_OR_M, U
    ), p = c(
      NA, p.OR01, p.OR10, p.OR11, p.OR01, p.OR_X1, p.OR_M, p.RERI
    ))
  rownames(d) <- NULL
  
  
  
  if (!em) {
    d <- data.frame(Measures = c(
      "OR00", "OR01", "OR10", "OR11", paste("OR(",
                                            beta2, " on outcome [", beta1, "==0]",
                                            sep = ""
      ), paste("OR(",
               beta2, " on outcome [", beta1, "==1]",
               sep = ""
      ), paste("OR(",
               beta1, " on outcome [", beta2, "==0]",
               sep = ""
      ), paste("OR(",
               beta1, " on outcome [", beta2, "==1]",
               sep = ""
      ), "Multiplicative scale",
      "RERI", "AP", "SI"
    ), Estimates = c(
      OR00, OR01, OR10, OR11,
      OR01, OR_X1, OR10, OR_A1, OR_M, RERI, AP, SI
    ), CI.ll = c(
      NA,
      l2, l1, l3, l2, CI.ll_OR_X1, l1, CI.ll_OR_A1, CI.ll_OR_M, L,
      APL, SIL
    ), CI.ul = c(
      NA, u2, u1, u3, u2, CI.ul_OR_X1, u1, CI.ul_OR_A1,
      CI.ul_OR_M, U, APU, SIU
    ),p = c(
      NA, p.OR01, p.OR10, p.OR11, p.OR01, p.OR_X1, p.OR10, p.OR_A1,
      p.OR_M, p.RERI, p.AP, p.SI
    ))
    rownames(d) <- NULL
  }
  
  raw_data <- model$data
  
  
  if (exists("dat.ir")) {
    raw_data <- dat.ir
  }
  
  
  
  fmt_p.signif <- function(x, digits = 3){
    pvalues <- x
    if(is.numeric(pvalues)){
      sign <- ifelse(pvalues <= 0.001, "***",
                     ifelse(pvalues <= 0.01, "**",
                            ifelse(pvalues <= 0.05, "*",
                                   ifelse(pvalues <= 0.1, ".", ""))))
      sign <- ifelse(is.na(sign), "", format(sign, justify = "left"))
      
      pvalues <- fmt_pvalue(pvalues, digits)
      pvalues <- ifelse(is.na(pvalues), "", format(pvalues, justify = "right"))
      
      # paste(pvalues, sign, sep = "")
      return(sign)
      
    }else{
      pvalues.numeric <- as.numeric(gsub(pattern = "<|>", replacement = "", x = pvalues))
      sign <- ifelse(pvalues.numeric <= 0.001, "***",
                     ifelse(pvalues.numeric <= 0.01, "**",
                            ifelse(pvalues.numeric <= 0.05, "*",
                                   ifelse(pvalues.numeric <= 0.1, ".", ""))))
      sign <- ifelse(is.na(sign), "", format(sign, justify = "left"))
      pvalues <- ifelse(is.na(pvalues), "", format(pvalues, justify = "right"))
      # paste(pvalues, sign, sep = "")
      
      return(sign)
    }
  }
  d$p = round(d$p, 5)
  d$p <- fmt_p.signif(d$p)
  
  
  
  
  
  
  if(is.null(model$call)){model$call = model$formula}
  
  ir <- list(dframe = d,
             # exp_names = c(beta1, beta2),
             exp_names = exposure_names,
             analysis = em,
             call = model$call,
             dat = raw_data,
             model = attr(model, "model"))
  
  attr(ir, "class") <- "interactionR"
  
  
  # invisible(ir)
  Inter_table(ir,p.value = T) %>% autoReg::myft(digits = 2) %>% print
  
  Inter_table(ir, p.value = p.value) %>%
    mutate_all(~replace(., is.na(.), "")) %>%
    slice(-1)
  
}



#########--------  Cox Analysis --------  ################
# This section focuses on Cox Proportional Hazards survival analysis for papillary thyroid cancer
# Analyzes sex disparities in survival using multiple causal inference methods

##--Cox---1. Survival curves analysis -------------
# Primary analysis comparing survival curves between sexes using different adjustment methods

# Define analysis variables
cat_var <- "Sex"  # Primary exposure variable (Male vs Female)
# Comprehensive set of adjustment variables for confounding control
adj_var = c("Age", "Sex", "Race","Eth","Marital_4","Income_3", "County_2","Yea_2","Grade", "Tum","Ext","Mul","N_3", "M", "Sur", "Rad" )

##--------1.1 Generate survival curve objects for different adjustment methods
# Create unadjusted and adjusted survival curve objects for sex comparison
cat.obj1 <- .get.obj(data, cat_var, adj_var,  method = c("unadj"), times = NULL, surv = T, conf_level = 0.95)    # Unadjusted Kaplan-Meier curves
cat.obj2 <- .get.obj(data, cat_var, adj_var,  method = c("direct"), times = NULL, surv = T, conf_level = 0.95)  # Direct standardization using Cox model

##--------1.2 Generate data for univariate/multivariate forest plots
# Create forest plot data comparing different statistical approaches
cat.UM <- .get.UM.lis(data, cat_var, adj_var, methods = c("unadj","adjus","unadj","direct" ) )

##--------1.3 get plotdata for forest Mul
group_levels <- levels(data[[cat_var]])
cat.Mul1 <- .get.sub.dif.times(data, cat_var, NULL, adj_var, surv = T, method = "unadj",
                               group_1 = group_levels[1],group_2 = group_levels[2],
                               times1 = c(60,120,180,240),
                               times2 = c(60,120,180,240),
                               conf_level1 = 0.95, 
                               conf_level2 = 0.95)
cat.Mul2 <- .get.sub.dif.times(data, cat_var, NULL, adj_var, surv = T, method = "direct",
                               group_1 = group_levels[1],group_2 = group_levels[2],
                               times1 = c(60,120,180,240),
                               times2 = c(60,120,180,240),
                               conf_level1 = 0.95, 
                               conf_level2 = 0.95)

## ------- plot
# curves
p1 <- .plot.obj(cat.obj1, ylim = c(0.9,1), ybre = seq(0.9, 1, 0.02)) + theme(legend.position="none")
p2 <- .plot.obj(cat.obj2, ylim = c(0.9,1), ybre = seq(0.9, 1, 0.02)) + theme(legend.position="none")
p3 <- .plot.dif(cat.obj1, ylim = c(-0.01,0.042), ybre = seq(-0.01,0.042, 0.01)) 
p4 <- .plot.dif(cat.obj2, ylim = c(-0.01,0.042), ybre = seq(-0.01,0.042, 0.01),color = lancet_palette[2]) 
# ( p1|p2 ) / ( p3|p4 )

# forest UM
dat_lis <- map(cat.UM, ~ slice(.x, -1))
core_size = c(2,4)
p5 <- .plot.forest_for.curve(dat_lis[1], xlim = c(0.8,2),  ticks_at = c(1,2),  core_size = core_size)
p6 <- .plot.forest_for.curve(dat_lis[2], xlim = c(0.8,2),  ticks_at = c(1,2),  core_size = core_size)
p7 <- .plot.forest_for.curve(dat_lis[3], xlim = c(-0.01, 0.025), ticks_at = c(-0.01,0,0.01,0.02), ticks_digits=2, core_size=core_size)
p8 <- .plot.forest_for.curve(dat_lis[4], xlim = c(-0.01, 0.025), ticks_at = c(-0.01,0,0.01,0.02), ticks_digits=2, core_size=core_size)

## Combine plots
P1 <- .insert.plot(p1, p5, pos = c(0.1,0.05)) 
P2 <- .insert.plot(p2, p6, pos = c(0.1,0.05)) 
P3 <- .insert.plot(p3, p7, pos = c(0.24,0.95)) 
P4 <- .insert.plot(p4, p8, pos = c(0.24,0.95)) 
ggsave(filename = "sex_age/1 Survival curves and difference curves between sexes .pdf", 
       plot = ( P1|P2 ) / ( P3|P4 ) + plot_annotation(tag_levels = c('A')),
       width = 20, height = 12)


# forest Mul
label_list =  list(label = "..Groups..", 
                   time = "Times",
                   effect1 = "Survival difference")
core_size = c(1.5,2)
p9  <-.plot.forest(list(cat.Mul1),fill_col = "pal_div[[5]][7:8]", 
                   xlim = c(-0.03, 0.06), ticks_at = seq(0, 0.06, 0.03), ticks_digits=2,
                   core_size = core_size, type = "v", 
                   summary_vars = c("label", "time"), label_list = label_list)
p10 <-.plot.forest(list(cat.Mul2),fill_col = "pal_seq[[5]][2:3]", 
                   xlim = c(-0.03, 0.06), ticks_at = seq(-0.03, 0.06, 0.03), ticks_digits=2,
                   core_size = core_size, type = "v", 
                   summary_vars = c("label", "time"), label_list = label_list)

## Combine plots
P5 <- .insert.plot(p3, p9,  pos = c(0.4,0.8)) 
P6 <- .insert.plot(p4, p10, pos = c(0.4,0.8)) 

ggsave(filename = "sex_age/1.1 Survival curves and difference curves between sexes .pdf", 
       plot = ( P1|P2 ) / ( P5|P6 ) + plot_annotation(tag_levels = c('A')),
       width = 20, height = 12)


##--Cox---2. Sex Effect by Age using Restricted Cubic Splines (RCS) ---------
## Analyze how sex effects on survival vary across age using non-linear modeling

# Generate RCS models for hazard ratios (relative effects)
# Unadjusted RCS model: Sex effect on hazard ratio across age spectrum
res.cat1 <- rcs.cox1_boot(data, var1=cat_var, var2=con_var, xlim=xlim, dead = "DSS", n=3, co_var = NULL,
                          conf = 0.95 , ci.boot.method = "perc", R = 1000, parallel = "multicore")
# Adjusted RCS model: Sex effect on hazard ratio controlling for confounders
res.cat2 <- rcs.cox1_boot(data, var1=cat_var, var2=con_var, xlim=xlim, dead = "DSS", n=3, co_var = adj_var,
                          conf = 0.95 , ci.boot.method = "perc", R = 1000, parallel = "multicore")
# Generate RCS models for absolute mortality rates (absolute effects using Poisson regression)
# Unadjusted Poisson RCS model: Sex effect on absolute mortality rates across age
res.poi1 <-rcs.poi_boot(data, var1=cat_var, var2=con_var, xlim=xlim, dead="DSS", n=3, co_var=NULL,
                        conf=0.95, ci.boot.method="perc", R=1000, parallel="multicore")
# Adjusted Poisson RCS model: Sex effect on absolute mortality rates controlling for confounders  
res.poi2 <-rcs.poi_boot(data, var1=cat_var, var2=con_var, xlim=xlim, dead="DSS", n=3, co_var=adj_var,
                        conf=0.95, ci.boot.method="perc", R=1000, parallel="multicore")


## Create RCS plots showing age-dependent sex effects
# Visualize how hazard ratios for sex vary across the age spectrum
# Plot 1: Unadjusted sex hazard ratios across age using RCS
p1 <- res.cox1_as_rcs(res.cat1, boot_curves = F, my_conf = T,conf.lim= c(0.15,0.85),
                      lin_color=lancet_palette[1], ylim= c(0,3),
                      histlimit = c(0,1.5), his_color = lancet_palette[3])
# Plot 2: Adjusted sex hazard ratios across age using RCS
p2 <- res.cox1_as_rcs(res.cat2, boot_curves = F, my_conf = T,conf.lim= c(0.15,0.85),
                      lin_color=lancet_palette[2], ylim= c(0,3),
                      histlimit = c(0,1.5), his_color = lancet_palette[3])

( p1|p2 )  + plot_annotation(tag_levels = c('A'))

ggsave(filename = "sex_age/2.1 Survival curves and difference curves between sexes .pdf", 
       plot = ( p1|p2 )  + plot_annotation(tag_levels = c('A')),
       width = 16, height = 6)

p3 <- res.poi_as_rcs(res.poi1, boot_curves = F, my_conf = T,conf.lim= c(0.15,0.9),
                     ylim= c(0,15), yexpand = c(0.01,0.01), histlimit = c(0,5))
p4 <- res.poi_as_rcs(res.poi2, boot_curves = F, my_conf = T,conf.lim= c(0.15,0.9),
                     ylim= c(0,15), yexpand = c(0.01,0.01), histlimit = c(0,5))
p5 <- res.poi_as_rcs(res.poi1) %>% 
  diff_poi(group1 = group_levels[2],group2 = group_levels[1],
           boot_curves = F, my_conf = T,conf.lim=c(0.1,0.9),
           ylim= c(-3,3),histlimit = c(0,1.5), lin_color = lancet_palette[1])
p6 <- res.poi_as_rcs(res.poi2) %>% 
  diff_poi(group1 = group_levels[2],group2 = group_levels[1],
           boot_curves = F, my_conf = T,conf.lim=c(0.1,0.9),
           ylim= c(-3,3),histlimit = c(0,1.5), lin_color = lancet_palette[2])

( p3|p4 ) / ( p5|p6 ) + plot_annotation(tag_levels = c('A'))

ggsave(filename = "sex_age/2.2 Survival curves and difference curves between sexes .pdf", 
       plot = ( p3|p4 ) / ( p5|p6 ) + plot_annotation(tag_levels = c('A')),
       width = 16, height = 11)


##--Cox---3. Sex effect by Age -- subgroup analysis -----------
# Analyze sex effects within predefined age subgroups using conventional stratified analysis

# Generate subgroup analysis data for sex effects across age categories
cat.sub <- .get.SUB.lis1( data, cat_var, sub_var, adj_var )

# .plt.for(cat.sub, "v")
core_size <- c(2,2)
p1 <- .plot.forest(cat.sub[c(2,4)],
                   xlim = list(c(0,4), c(-0.04,0.05)),
                   ticks_at =list(c(0,1,2,3,4), c(-0.04,0,0.04)),
                   ticks_digits = c(0, 2),
                   core_size = core_size,
                   fill_col ="pal_div[[5]][3:9]",
                   type = c("h"),
                   summary_vars = c("label", "eve1"),
                   label_list = list(label = "Age groups", 
                                     eve1 = "event/total",
                                     p_inter1 = "p interaction",
                                     effect1 = "HR (95%)",
                                     effect2 = "Survival difference"))

p2 <- .plot.forest(cat.sub[c(1,3)],
                   xlim = list(c(0,4), c(-0.04,0.05)),
                   ticks_at =list(c(0,1,2,3,4), c(-0.04,0,0.04)),
                   ticks_digits = c(0, 2),
                   core_size = core_size,
                   fill_col ="pal_div[[5]][3:9]",
                   type = c("h"),
                   summary_vars = c("label", "eve1"),
                   label_list = list(label = "Age groups", 
                                     eve1 = "event/total",
                                     p_inter1 = "p interaction",
                                     effect1 = "HR (95%)",
                                     effect2 = "Survival difference"))


ggsave(filename = "sex_age/3.1.Adjusted hazard ratio and 10-year survival probability by Sex across different age subgroups based on Cox model.pdf", 
       plot = p1,
       width = 13, height = 5.5)

ggsave(filename = "sex_age/3.2.Unadjusted hazard ratio and 10-year survival probability by Sex across different age subgroups based on Cox model.pdf", 
       plot = p2,
       width = 13, height = 5.5)



##--Cox---4. Formal interaction effect analysis between sex and age -----------
# Statistical testing for multiplicative and additive interaction effects

##--------4.1 Survival curves analysis for sex-age interaction groups
# Create combined exposure variable for joint sex-age analysis
data <- data %>% com_factor(exposure_names[1], exposure_names[2], group.name = "group")
cat_var <- "group"  # Combined sex-age group variable
# Adjustment variables excluding sex and age (now combined in group variable)
adj_var = c("Race","Eth","Marital_4","Income_3", "County_2","Yea_2","Grade", "Tum","Ext","Mul","N_3", "M", "Sur", "Rad" )

# get plotdata
group.obj1 <- .get.obj(data, cat_var, adj_var,  method = c("unadj"), times = NULL, surv = T, conf_level = 0.90)
group.obj2 <- .get.obj(data, cat_var, adj_var,  method = c("direct"), times = NULL, surv = T, conf_level = 0.90)

group.UM <- .get.UM.lis(data, cat_var, adj_var, methods = c("unadj","adjus","unadj","direct" ) )

# plot
dat_lis <- map(group.UM, ~ slice(.x, -1))
core_size = c(2,4)
p1 <- .plot.obj(group.obj1, ylim = c(0.85,1), ybre = seq(0.85, 1, 0.05)) + theme(legend.position="none")
p2 <- .plot.obj(group.obj2, ylim = c(0.85,1), ybre = seq(0.85, 1, 0.05)) + theme(legend.position="none")
p3 <- .plot.forest_for.curve(dat_lis[c(1,3)], core_size = core_size, 
                             xlim = list(c(1,8),c(0, 0.075)), ticks_at = list(c(1,2,4,6,8),c(0,0.06,0.03)), ticks_digits=c(0,2))
p4 <- .plot.forest_for.curve(dat_lis[c(2,4)], core_size = core_size,
                             xlim = list(c(1,8),c(0, 0.075)), ticks_at = list(c(1,2,4,6,8),c(0,0.06,0.03)), ticks_digits=c(0,2))

## combine plot
P1 <- .insert.plot(p1, p3, pos = c(0.25,0.03)) 
P2 <- .insert.plot(p2, p4, pos = c(0.25,0.03)) 

ggsave(filename = "sex_age/4.1 Survival curves for different sex and age groups based on unadjusted (A) and adjusted (B) Cox model.pdf", 
       plot = ( P1/P2 ) + plot_annotation(tag_levels = c('A')),
       width = 14, height = 14)


##--------4.2 interaction table between sex and age
res.inter <- .get.interTable(data, exposure_names, adj_var) 

save.tb(res.inter,
        path = "sex_age/table.interaction.docx",
        title = "Table 2. Interaction analysis between sex and age")



#########--------  Competing Risks Regression (CRR) Analysis --------  ################ 
# This section focuses on Competing Risks Regression using the Fine-Gray model
# Analyzes cancer-specific mortality while accounting for other competing causes of death
##--Crr---1. Cumulative incidence function analysis -------------
##--Crr---1.1 Generate competing risks objects for different adjustment methods
# Create cumulative incidence function objects for cancer-specific death (competing risks analysis)
cat.obj1 <- .get.obj(data, cat_var, adj_var,  method = c("unadj"),  times = NULL, surv = F, conf_level = 0.95)    # Unadjusted Aalen-Johansen estimator
cat.obj2 <- .get.obj(data, cat_var, adj_var,  method = c("direct"), times = NULL, surv = F, conf_level = 0.95)  # Direct standardization using Fine-Gray model

##--------1.2 Generate subdistribution hazard ratios (competing risks forest plots)
# Create forest plot data for competing risks regression models
cat.UM <- .get.UM.lis(data, cat_var, adj_var,surv = F, methods = c("unadj","adjus","unadj","direct" ) )
dat_lis <- map(cat.UM, ~ slice(.x, -1))

##--------1.3 get survival difference 
group_levels <- levels(data[[cat_var]])
cat.Mul1 <- .get.sub.dif.times(data, cat_var, NULL, adj_var, surv = F, method = "unadj",
                               group_1 = group_levels[1],group_2 = group_levels[2],
                               times1 = c(60,120,180,240),
                               times2 = c(60,120,180,240),
                               conf_level1 = 0.95, 
                               conf_level2 = 0.95)
cat.Mul2 <- .get.sub.dif.times(data, cat_var, NULL, adj_var, surv = F, method = "direct",
                               group_1 = group_levels[1],group_2 = group_levels[2],
                               times1 = c(60,120,180,240),
                               times2 = c(60,120,180,240),
                               conf_level1 = 0.95, 
                               conf_level2 = 0.95)

## plot
# curves
p1 <- .plot.obj(cat.obj1, ylim = c(0,0.1), ybre = seq(0, 0.1, 0.02)) + theme(legend.position="none")
p2 <- .plot.obj(cat.obj2, ylim = c(0,0.1), ybre = seq(0, 0.1, 0.02)) + theme(legend.position="none")
p3 <- .plot.dif(cat.obj1, ylim = c(-0.04,0.01), ybre = seq(-0.04,0.01,0.01)) 
p4 <- .plot.dif(cat.obj2, ylim = c(-0.04,0.01), ybre = seq(-0.04,0.01,0.01),color = lancet_palette[2]) 
# ( p1|p2 ) / ( p3|p4 )
# UM forest
core_size = c(2,4)
p5 <- .plot.forest_for.curve(dat_lis[1], xlim = c(0.7,2),  ticks_at = c(1,2),  core_size = core_size)
p6 <- .plot.forest_for.curve(dat_lis[2], xlim = c(0.7,2),  ticks_at = c(1,2),  core_size = core_size)
p7 <- .plot.forest_for.curve(dat_lis[3], xlim = c( -0.022,0.01), ticks_at = c(-0.02,-0.01,0,0.01), ticks_digits=2, core_size=core_size)
p8 <- .plot.forest_for.curve(dat_lis[4], xlim = c( -0.022,0.01), ticks_at = c(-0.02,-0.01,0,0.01), ticks_digits=2, core_size=core_size)

# mul forests
label_list =  list(label = "..Groups..", 
                   time = "Times",
                   effect1 = "Survival difference")
core_size = c(1.5,2)
p9  <-.plot.forest(list(cat.Mul1),fill_col = "pal_div[[5]][7:8]", 
                   xlim = c(-0.045, 0.02), ticks_at = c(-0.04, -0.02, 0, 0.02), ticks_digits=2,
                   core_size = core_size, type = "v", 
                   summary_vars = c("label", "time"), label_list = label_list)
p10 <-.plot.forest(list(cat.Mul2),fill_col = "pal_seq[[5]][2:3]", 
                   xlim = c(-0.045, 0.02), ticks_at = c(-0.04, -0.02, 0, 0.02), ticks_digits=2,
                   core_size = core_size, type = "v", 
                   summary_vars = c("label", "time"), label_list = label_list)


## combine plot
P1 <- .insert.plot(p1, p5, pos = c(0.1,0.95)) 
P2 <- .insert.plot(p2, p6, pos = c(0.1,0.95)) 

P3 <- .insert.plot(p3, p7, pos = c(0.24,0.05)) 
P4 <- .insert.plot(p4, p8, pos = c(0.24,0.05)) 

ggsave(filename = "sex_age/Sup1. Survival curves and difference curves between sexes.pdf", 
       plot = ( P1|P2 ) / ( P3|P4 ) + plot_annotation(tag_levels = c('A')),
       width = 20, height = 12)

P5 <- .insert.plot(p3, p9,  pos = c(0.4,0.14)) 
P6 <- .insert.plot(p4, p10, pos = c(0.4,0.14)) 

ggsave(filename = "sex_age/Sup1.1 Survival curves and difference curves between sexes .pdf", 
       plot = ( P1|P2 ) / ( P5|P6 ) + plot_annotation(tag_levels = c('A')),
       width = 20, height = 12)



##--Crr---3. Sex effect by age -- subgroup analysis (competing risks) -----------
# Analyze sex effects within age subgroups using competing risks framework

# Generate subgroup analysis data for sex effects across age categories in competing risks setting
cat.sub <- .get.SUB.lis1( data, cat_var, sub_var, adj_var, surv = F,  
                          methods = c("unadj","adjus","unadj","direct" ) )
# plot
core_size <- c(2,2)
p1 <- .plot.forest(cat.sub[c(2,4)],
                   xlim = list(c(0,4), c(-0.04,0.05)),
                   ticks_at =list(c(0,1,2,3,4), c(-0.04,0,0.04)),
                   ticks_digits = c(0, 2),
                   core_size = core_size,
                   fill_col ="pal_div[[5]][3:9]",
                   type = c("h"),
                   summary_vars = c("label", "eve1"),
                   label_list = list(label = "Age groups", 
                                     eve1 = "event/total",
                                     p_inter1 = "p interaction",
                                     effect1 = "HR (95%)",
                                     effect2 = "Survival difference"))

p2 <- .plot.forest(cat.sub[c(1,3)],
                   xlim = list(c(0,4), c(-0.04,0.05)),
                   ticks_at =list(c(0,1,2,3,4), c(-0.04,0,0.04)),
                   ticks_digits = c(0, 2),
                   core_size = core_size,
                   fill_col ="pal_div[[5]][3:9]",
                   type = c("h"),
                   summary_vars = c("label", "eve1"),
                   label_list = list(label = "Age groups", 
                                     eve1 = "event/total",
                                     p_inter1 = "p interaction",
                                     effect1 = "HR (95%)",
                                     effect2 = "Survival difference"))


ggsave(filename = "sex_age/Sup3.1.Adjusted hazard ratio and 10-year survival probability by Sex across different age subgroups based on Fine-Gray model.pdf", 
       plot = p1,
       width = 13, height = 5.5)

ggsave(filename = "sex_age/Sup3.2.Unadjusted hazard ratio and 10-year survival probability by Sex across different age subgroups based on Fine-Gray model.pdf", 
       plot = p2,
       width = 13, height = 5.5)



##--Crr---4. Interaction effect analysis between sex and age (competing risks) -----------
# Statistical testing for interaction effects using competing risks framework

##--------4.1 Cumulative incidence curves for sex-age interaction groups
# Create combined exposure variable for joint sex-age analysis in competing risks setting
data <- data %>% com_factor(exposure_names[1], exposure_names[2], group.name = "group")
cat_var <- "group"  # Combined sex-age group variable for competing risks
# Adjustment variables for competing risks models
adj_var = c("Race","Eth","Marital_4","Income_3", "County_2","Yea_2","Grade", "Tum","Ext","Mul","N_3", "M", "Sur", "Rad" )

# get plotdata for curves
group.obj1 <- .get.obj(data, cat_var, adj_var,  method = c("unadj"), times = NULL, surv = F, conf_level = 0.90)
group.obj2 <- .get.obj(data, cat_var, adj_var,  method = c("direct"), times = NULL, surv = F, conf_level = 0.90)
# get plotdata for forest
group.UM <- .get.UM.lis(data, cat_var, adj_var,surv = F, methods = c("unadj","adjus","unadj","direct" ) )
dat_lis <- map(group.UM, ~ slice(.x, -1))


# plot
core_size = c(2,4)
p1 <- .plot.obj(group.obj1, ylim = c(0,0.15), ybre = seq(0, 0.15, 0.05)) + theme(legend.position="none")
p2 <- .plot.obj(group.obj2, ylim = c(0,0.15), ybre = seq(0, 0.15, 0.05)) + theme(legend.position="none")
p3 <- .plot.forest_for.curve(dat_lis[c(1,3)], core_size = core_size, 
                             xlim = list(c(1,8),c(-0.065, 0)), ticks_at = list(c(1,2,4,6,8),c(-0.06,-0.03,0)), ticks_digits=c(0,2))
p4 <- .plot.forest_for.curve(dat_lis[c(2,4)], core_size = core_size,
                             xlim = list(c(1,8),c(-0.065, 0)), ticks_at = list(c(1,2,4,6,8),c(-0.06,-0.03,0)), ticks_digits=c(0,2))

## combine
P1 <- .insert.plot(p1, p3, pos = c(0.25,0.92)) 
P2 <- .insert.plot(p2, p4, pos = c(0.25,0.92)) 

ggsave(filename = "sex_age/Sup4.1 Survival curves for different sex and age groups based on unadjusted (A) and adjusted (B) Fine-Gray model.pdf", 
       plot = ( P1/P2 ) + plot_annotation(tag_levels = c('A')),
       width = 14, height = 14)


##--------4.2 Generate formal interaction table between sex and age (competing risks)
# Create comprehensive interaction analysis table using Fine-Gray competing risks model
res.inter <- .get.interTable(data, exposure_names, adj_var, surv = F) 

# Save interaction analysis results as formatted table
save.tb(res.inter,
        path = "sex_age/Supplementary Table 1.Interaction analysis between sex and age based on Fine-Gray model.docx",
        title = "Supplementary Table 1. Interaction analysis between sex and age based on Fine-Gray model")

# ==============================================================================
# END OF ANALYSIS SCRIPT
# ==============================================================================
# 
# This comprehensive analysis script provides:
# 1. Cox proportional hazards models for relative effects (hazard ratios)
# 2. Fine-Gray competing risks models for absolute effects (cumulative incidence)
# 3. Restricted cubic spline analysis for non-linear age effects
# 4. Formal interaction testing on multiplicative and additive scales
# 5. Multiple causal inference methods (unadjusted, direct standardization, IPTW, PSM)
# 6. Comprehensive visualization and forest plot generation
# 
# Key findings focus on sex disparities in papillary thyroid cancer survival
# across the age spectrum, demonstrating divergent patterns of relative vs
# absolute effects with important clinical and public health implications.
# ==============================================================================



