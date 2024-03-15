# Finn McQuaid 25/08/2023
# Pan-TB impact on amplification of resistance
setwd("C:/Users/eidefmcq/Documents/Simulations/PanTB")
rm(list=ls())
## ================ Load packages ==============================================
library(ggforce)
library(dplyr)
library(ggplot2)
library(ggalt)
library(dampack)
library(rriskDistributions)
library(patchwork)
library("data.table")
library(stringr)
library(tidyverse)

## ================ Function for parameters  ===================================
param_table <- data.frame(E_R    = get.beta.par(p=c(0.025,0.50,0.975), q=c(0.585, 0.709, 0.793), show.output = FALSE),   # Treatment efficacy for rifamycin-based regimen,	Aims 1-2
            E_BX              = get.beta.par(p=c(0.025,0.50,0.975), q=c(0.688, 0.763, 0.831), show.output = FALSE), # Treatment efficacy for pan-TB regimen,	Aims 1-2
            E_ind             = get.beta.par(p=c(0.025,0.50,0.975), q=c(0.337, 0.439, 0.536), show.output = FALSE), # Treatment efficacy for individualised regimen,	Aims 1-2
            CFR               = get.beta.par(p=c(0.025,0.50,0.975), q=c(0.403, 0.483, 0.525), show.output = FALSE), # Case fatality rate given no durable cure, GTB report
            P_R               = get.beta.par(p=c(0.025,0.50,0.975), q=c(0.220, 0.350, 0.500), show.output = FALSE), # Risk ratio of cure for rifamycin-based regimen given RR-TB,	Aims 1-2
            P_B               = get.beta.par(p=c(0.025,0.50,0.975), q=c(0.540, 0.750, 0.910), show.output = FALSE), # Risk ratio of cure for pan-TB regimen given BR-TB,	Ismail 2021 UI  (0.3-1) https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(21)00470-9/fulltext
            #P_X               = get.beta.par(p=c(0.025,0.50,0.975), q=c(0.814, 0.875, 0.957), show.output = FALSE), # Risk ratio of cure for pan-TB regimen given XR-TB,	Gegia 2017 for those treated with WHO-New only, comparing failure/relapse (16% [10–21] vs 4% [3–6]) https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(16)30407-8/fulltext
            P_X               = get.beta.par(p=c(0.025,0.50,0.975), q=c(0.540, 0.750, 0.910), show.output = FALSE), # Risk ratio of cure for pan-TB regimen given XR-TB,	Gegia 2017 for those treated with WHO-New only, comparing failure/relapse (16% [10–21] vs 4% [3–6]) https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(16)30407-8/fulltext
            P_BX              = get.beta.par(p=c(0.025,0.50,0.975), q=c(0.220, 0.350, 0.500), show.output = FALSE), # Risk ratio of cure for pan-TB regimen given BR- and XR-TB,	Aims 1-2, presumably less than rows above if not rethink
            S_R               = get.beta.par(p=c(0.025,0.50,0.975), q=c(0.003, 0.006, 0.012), show.output = FALSE), # Likelihood of RR acquisition after rifamycin-based treatment, Menzies 2009+Kendall 2017 UI (0.003–0.012) https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1000146 https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0172748
            S_B               = get.beta.par(p=c(0.025,0.50,0.975), q=c(0.003, 0.01, 0.023), show.output = FALSE), # Likelihood of BR acquisition after pan-TB treatment,	Ismail 2021, Mallick 2022, The former has 2.3%, the latter review (which doesn't include the former paper) has 2.2% (IQR 1.1%–4.6%) for phenotypic resistance and 4.4% (IQR 1.8%–5.8%) for genotypic https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(21)00470-9/fulltext https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8963286/ 
            #S_X               = get.beta.par(p=c(0.025,0.50,0.975), q=c(0.003, 0.006, 0.009), show.output = FALSE), # Likelihood of XR acquisition after pan-TB treatment,	Ismail 2021, Mallick 2022. We could e.g. use Gegia 2017 to update this to reflect INH instead if better, as these are low (0·6% [0·3–0·9], essentially the same as Menzies et al although this may be for all drugs not just INH)
            S_X               = get.beta.par(p=c(0.025,0.50,0.975), q=c(0.003, 0.01, 0.023), show.output = FALSE), # Likelihood of XR acquisition after pan-TB treatment,	Ismail 2021, Mallick 2022. We could e.g. use Gegia 2017 to update this to reflect INH instead if better, as these are low (0·6% [0·3–0·9], essentially the same as Menzies et al although this may be for all drugs not just INH)
            Q                 = c(shape1=4.00, shape2=16.00), # Risk ratio for resistance acquisition given existing resistance to one component of pan-TB	If we assumed this was ~INH for X, then we could use Gegia 2017 for relative rates given existing resistance (8% acquired resistance for WHO_New when there is existing INH-R, only 1% acquired resistance when there is no INH-R). Kunkel 2016 supp material pp5 may not work as too many unknowns
            rho_soc_new       = get.beta.par(p=c(0.025,0.50,0.975), q=c(0.441*0.75, 0.441, 0.441*1.25), show.output = FALSE), # SoC new patient Rho DST coverage	GTB report, weighted by proportion with bact conf https://www.who.int/teams/global-tuberculosis-programme/tb-reports/global-tuberculosis-report-2022/tb-diagnosis-treatment/3-2-diagnostic-testing-for-tb--hiv-associated-tb-and-drug-resistant-tb
            betachi_soc_new   = get.beta.par(p=c(0.025,0.50,0.975), q=c(0.490*0.75, 0.490, 0.490*1.25), show.output = FALSE), # SoC new patient 2nd line (Beta+Chi) DST coverage given RR-TB,	GTB report for fqr testing coverage, https://www.who.int/teams/global-tuberculosis-programme/tb-reports/global-tuberculosis-report-2022/tb-diagnosis-treatment/3-2-diagnostic-testing-for-tb--hiv-associated-tb-and-drug-resistant-tb
            rho_pan           = get.beta.par(p=c(0.025,0.50,0.975), q=c(0.441*0.75, 0.441, 0.441*1.25), show.output = FALSE), # Pan-TB retreatment Rho DST coverage,	Assumption
            betachi_pan       = get.beta.par(p=c(0.025,0.50,0.975), q=c(0.490*0.75, 0.490, 0.490*1.25), show.output = FALSE), # Pan-TB retreatment beta/Chi DST coverage,	Assumption
            rho_soc_retR      = get.beta.par(p=c(0.025,0.50,0.975), q=c(0.800*0.75, 0.800, 0.800*1.25), show.output = FALSE), # SoC Rho DST coverage for those previously treated with R,	Aims 1-2
            rho_soc_retBX     = c(shape1=1, shape2=1), # SoC Rho DST coverage for those previously treated with BX,	Assumption
            betachi_soc_retR  = get.beta.par(p=c(0.025,0.50,0.975), q=c(0.490*0.75, 0.490, 0.490*1.25), show.output = FALSE), # SoC 2nd line (Beta+Chi) DST coverage given RR-TB for those previously treated with R,	GTB report for fqr testing coverage
            betachi_soc_retBX = get.beta.par(p=c(0.025,0.50,0.975), q=c(0.600*0.75, 0.600, 0.600*1.25), show.output = FALSE)) # SoC 2nd line (Beta+Chi) DST coverage given RR-TB for those previously treated with BX,	Assumtion higher than the row above
  
params_sample <- function(param_table){
    params <- data.frame(E_R               = rbeta(1, shape1=param_table[1,"E_R"], shape2=param_table[2,"E_R"]),
                      E_BX                 = rbeta(1, shape1=param_table[1,"E_BX"], shape2=param_table[2,"E_BX"]),  
                      E_ind                = rbeta(1, shape1=param_table[1,"E_ind"], shape2=param_table[2,"E_ind"]),  
                      CFR                  = rbeta(1, shape1=param_table[1,"CFR"], shape2=param_table[2,"CFR"]), 
                      P_R                  = rbeta(1, shape1=param_table[1,"P_R"], shape2=param_table[2,"P_R"]),
                      P_B                  = rbeta(1, shape1=param_table[1,"P_B"], shape2=param_table[2,"P_B"]), 
                      P_X                  = rbeta(1, shape1=param_table[1,"P_X"], shape2=param_table[2,"P_X"]),
                      P_BX                 = rbeta(1, shape1=param_table[1,"P_BX"], shape2=param_table[2,"P_BX"]),
                      S_R                  = rbeta(1, shape1=param_table[1,"S_R"], shape2=param_table[2,"S_R"]), 
                      S_B                  = rbeta(1, shape1=param_table[1,"S_B"], shape2=param_table[2,"S_B"]),
                      S_X                  = rbeta(1, shape1=param_table[1,"S_X"], shape2=param_table[2,"S_X"]),
                      Q                    = runif(1, min=param_table[1,"Q"], max=param_table[2,"Q"]),
                      rho_soc_new          = rbeta(1, shape1=param_table[1,"rho_soc_new"], shape2=param_table[2,"rho_soc_new"]),
                      betachi_soc_new      = rbeta(1, shape1=param_table[1,"betachi_soc_new"], shape2=param_table[2,"betachi_soc_new"]), 
                      rho_pan              = rbeta(1, shape1=param_table[1,"rho_pan"], shape2=param_table[2,"rho_pan"]), 
                      # TODO ASSUMING THIS IS ZERO IN TEH FIRST INSTANCE, VARY FOR SENSITIVITY ANALYSIS 
                      betachi_pan          = 0,#rbeta(1, shape1=param_table[1,"betachi_pan"], shape2=param_table[2,"betachi_pan"]),
                      rho_soc_retR         = rbeta(1, shape1=param_table[1,"rho_soc_retR"], shape2=param_table[2,"rho_soc_retR"]),
                      rho_soc_retBX        = 1, 
                      betachi_soc_retR     = rbeta(1, shape1=param_table[1,"betachi_soc_retR"], shape2=param_table[2,"betachi_soc_retR"]), 
                      betachi_soc_retBX    = rbeta(1, shape1=param_table[1,"betachi_soc_retBX"], shape2=param_table[2,"betachi_soc_retBX"]))   
  return(params)
}

# Baseline parameter set
param_ext <- data.frame(E_R               = c(0.581, 0.709, 0.8),   
                        E_BX              = c(0.688, 0.763, 0.831), 
                        E_ind             = c(0.337, 0.439, 0.536), 
                        CFR               = c(0.403, 0.483, 0.525), 
                        P_R               = c(0.220, 0.350, 0.500), 
                        P_B               = c(0.540, 0.750, 0.910), 
                        # P_X value to align with bdq    
                        P_X               = c(0.540, 0.750, 0.910), 
                        P_BX              = c(0.220, 0.350, 0.500), 
                        S_R               = c(0.003, 0.006, 0.012), 
                        #S_B               = c(0.003, 0.023, 0.08), # FOR THINKING ABOUT VERY HIGH ACQUISITION RATES
                        S_B               = c(0.003, 0.010, 0.023),
                        # S_X value to align with bdq
                        S_X               = c(0.003, 0.010, 0.023), 
                        Q                 = c(4.000, 7.500, 16.000),
                        rho_soc_new       = c(0.441*0.75, 0.441, 0.441*1.25), 
                        betachi_soc_new   = c(0.490*0.75, 0.490, 0.490*1.25), 
                        rho_pan           = c(0.441*0.75, 0.441, 0.441*1.25),#c(0, 0, 0),
                        betachi_pan       = c(0, 0, 0), #c(0.490*0.75, 0.490, 0.490*1.25),
                        rho_soc_retR      = c(0.800*0.75, 0.800, 0.800*1.25), 
                        rho_soc_retBX     = c(1, 1, 1),
                        betachi_soc_retR  = c(0.490*0.75, 0.490, 0.490*1.25), 
                        betachi_soc_retBX = c(0.600*0.75, 0.600, 0.600*1.25)) 
## ================ Define  prevalence =========================================           
 prev <- data.frame(prev_DS = 1-0.042-0.002-0.009,	# Prevalence of DS-TB,	GTB report https://www.who.int/teams/global-tuberculosis-programme/tb-reports/global-tuberculosis-report-2022
            prev_RR         = 0.042,	# Prevalence of RR-TB,	GTB report 3.6% (95% UI: 2.7–4.4%) new, 18% (95% UI: 11–26%) previous, https://www.who.int/teams/global-tuberculosis-programme/tb-reports/global-tuberculosis-report-2022
            prev_BR         = 0.002,	# Prevalence of BR-TB,	Assumption
            prev_XR         = 0.000,	# Prevalence of XR-TB,	Assumption
            prev_RRBR       = 0.009,	# Prevalence of RR- and BR-TB,	Ismail 2020: 2.2% (1.4-3.4) of RR patients, see also Liu 2021 https://journals.asm.org/doi/10.1128/aac.00479-20 https://pubmed.ncbi.nlm.nih.gov/32667984/ 
            prev_RRXR       = 0.000,	# Prevalence of RR- and XR-TB,	Assumption
            prev_BRXR       = 0.000,	# Prevalence of BR- and XR-TB,	Assumption
            prev_RRBRXR     = 0.000)	# Prevalence of RR-, BR- and XR-TB,	Assumption)

## ================ Function for treatment outcomes ============================
 # Define outcomes for regimen R across different resistance phenotypes
outcomes <- function(params){
  R_cure    <- params$E_R * c(1, params$P_R, 1, 1, params$P_R, params$P_R, 1, params$P_R)
  R_death      <- params$CFR * (1 - R_cure)
  R_failure    <- 1 - R_cure - R_death
  R_acquiredRR <- c(params$S_R, 0, params$S_R, params$S_R, 0, 0, params$S_R, 0)
  R_acquiredRR <- pmin(R_failure,R_acquiredRR)
  R_acquiredBR <- rep(0, 8)
  R_acquiredXR <- rep(0, 8)
  R_outcomes   <- c(R_cure, R_failure, R_death, R_acquiredRR, R_acquiredBR, R_acquiredXR)
  # Define outcomes for regimen BX across different resistance phenotypes
  BX_cure    <- params$E_BX * c(1, 1, params$P_B, params$P_X, params$P_B, params$P_X, params$P_BX, params$P_BX)
  BX_death      <- params$CFR * (1 - BX_cure)
  BX_failure    <- 1 - BX_cure - BX_death
  BX_acquiredRR <- rep(0, 8)
  BX_acquiredBR <- c(params$S_B, params$S_B, 0, params$S_B*params$Q, 0, params$S_B*params$Q, 0, 0)
  BX_acquiredBR <- pmin(BX_failure,BX_acquiredBR)
  BX_acquiredXR <- c(params$S_X, params$S_X, params$S_X*params$Q, 0, params$S_X*params$Q, 0, 0, 0)
  BX_acquiredXR <- pmin(BX_failure,BX_acquiredXR)
  BX_outcomes   <- c(BX_cure, BX_failure, BX_death, BX_acquiredRR, BX_acquiredBR, BX_acquiredXR)
  # Define outcomes for individualised regimen across different resistance phenotypes
  Ind_cure    <- rep(params$E_ind, 8)
  Ind_death      <- params$CFR * (1 - Ind_cure)
  Ind_failure    <- 1 - Ind_cure - Ind_death
  Ind_acquiredRR <- rep(0, 8)
  Ind_acquiredBR <- rep(0, 8)
  Ind_acquiredXR <- rep(0, 8)
  Ind_outcomes   <- c(Ind_cure, Ind_failure, Ind_death, Ind_acquiredRR, Ind_acquiredBR, Ind_acquiredXR)
  # Define outcomes for B-based regimen  across different resistance phenotypes, currently assumed to be an average of BX used for DS-TB and an individualised regimen. Note only used for those with XR or RR+XR-TB
  #Bbased_cure    <- (BX_cure[1] + Ind_cure)/2 * c(NA, NA, NA, 1, NA, 1, NA, NA)
  Bbased_cure    <- (BX_cure[1] + BX_cure[4])/2 * c(NA, NA, NA, 1, NA, 1, NA, NA)
  #Bbased_death      <- (BX_death[1] + Ind_death)/2 * c(NA, NA, NA, 1, NA, 1, NA, NA)
  Bbased_death      <- (BX_death[1] + BX_death[4])/2 * c(NA, NA, NA, 1, NA, 1, NA, NA)
  #Bbased_failure    <- (BX_failure[1] + Ind_failure)/2 * c(NA, NA, NA, 1, NA, 1, NA, NA)
  Bbased_failure    <- (BX_failure[1] + BX_failure[4])/2 * c(NA, NA, NA, 1, NA, 1, NA, NA)
  Bbased_acquiredRR <- rep(0, 8)
  Bbased_acquiredBR <- (BX_acquiredBR[1] + BX_acquiredBR[4])/2 * c(NA, NA, NA, 1, NA, 1, NA, NA)
  Bbased_acquiredXR <- rep(0, 8)
  Bbased_outcomes   <- c(Bbased_cure, Bbased_failure, Bbased_death, Bbased_acquiredRR, Bbased_acquiredBR, Bbased_acquiredXR)
  # Define outcomes for X-based regimen  across different resistance phenotypes, currently assumed to be an average of BX used for DS-TB and an individualised regimen. Note only used for those with BR or RR+BR-TB
  #Xbased_cure    <- (BX_cure[1] + Ind_cure)/2 * c(NA, NA, 1, NA, 1, NA, NA, NA)
  #Xbased_death      <- (BX_death[1] + Ind_death)/2 * c(NA, NA, 1, NA, 1, NA, NA, NA)
  #Xbased_failure    <- (BX_failure[1] + Ind_failure)/2 * c(NA, NA, 1, NA, 1, NA, NA, NA)
  Xbased_cure    <- (BX_cure[1] + BX_cure[3])/2 * c(NA, NA, 1, NA, 1, NA, NA, NA)
  Xbased_death      <- (BX_death[1] + BX_death[3])/2 * c(NA, NA, 1, NA, 1, NA, NA, NA)
  Xbased_failure    <- (BX_failure[1] + BX_failure[3])/2 * c(NA, NA, 1, NA, 1, NA, NA, NA)
  Xbased_acquiredRR <- rep(0, 8)
  Xbased_acquiredBR <- rep(0, 8)
  Xbased_acquiredXR <- (BX_acquiredXR[1] + BX_acquiredXR[3])/2 * c(NA, NA, 1, NA, 1, NA, NA, NA)
  Xbased_outcomes   <- c(Xbased_cure, Xbased_failure, Xbased_death, Xbased_acquiredRR, Xbased_acquiredBR, Xbased_acquiredXR)
  # Compile into an array
  row.names <- c("DS", "RR", "BR", "XR", "RRBR", "RRXR", "BRXR", "RRBRXR")
  column.names <- c("Cure", "Failure", "Death", "Acquired RR", "Acquired BR", "Acquired XR")
  matrix.names <- c("R", "BX", "Bbased", "Xbased", "Ind")
  outcome <- array(c(R_outcomes, BX_outcomes, Bbased_outcomes, Xbased_outcomes, Ind_outcomes),
                   dim = c(8,6,5), dimnames = list(row.names,column.names, matrix.names))
  return(outcome)
  } 

## ================ Function for sankey plots ==================================
staged_sankey <- function(sankey_data) {
# Defining colour palette
  pal <- c("      Cure" = '#A9A9A9',
           "     DS"    = '#AE0D0A',
           "     RR"    = '#177E89',
           "    BR"     = '#CBA715',
           "    XR"     = '#A5C994', 
           "   RR/BR"   = '#CBA715',
           "   RR/XR"   = '#A5C994',
           "  BRXR"    = '#754668',
           " RR/BR/XR"  = '#754668',
           "Death"         = '#030303')
# Defining states and data frame
  start_states <- c("     DS", "     RR", "    BR", "    XR", "   RR/BR", "   RR/XR", "  BR/XR", " RR/BR/XR")
  mid_states <- c("      Cure", "     DS", "     RR", "    BR", "    XR", "   RR/BR", "   RR/XR", "  BR/XR", " RR/BR/XR","Death")
  end_states <- c("      Cure", "     DS", "     RR", "    BR", "    XR", "   RR/BR", "   RR/XR", "  BR/XR", " RR/BR/XR","Death")
  I <- data.frame(Start = rep(c("     DS", "     RR", "    BR", "    XR", "   RR/BR", "   RR/XR", "  BR/XR", " RR/BR/XR"),each=100),
                  Mid = rep(c("      Cure", "     DS", "     RR", "    BR", "    XR", "   RR/BR", "   RR/XR", "  BR/XR", " RR/BR/XR","Death"),each = 10),
                  End = rep(c("      Cure", "     DS", "     RR", "    BR", "    XR", "   RR/BR", "   RR/XR", "  BR/XR", " RR/BR/XR","Death"), 80),
                  value = rep(0,800))
# Populating frequency table  
  exc <- c(sankey_data$DS,sankey_data$RR,sankey_data$BR,sankey_data$XR,sankey_data$RRBR,sankey_data$RRXR,sankey_data$BRXR,sankey_data$RRBRXR)
  I$value <- exc
# Add in time points and colours
I_new <- I %>%
  gather_set_data(1:3) %>%
  mutate(y = factor(y, levels = end_states)) %>%
  mutate(x = factor(x, levels = c(1,2,3)))
  #I_new <- I %>% 
    #gather_set_data(1:2) %>%
    #mutate(y = factor(y, levels = end_states)) %>%
   # mutate(x = factor(x, levels = c(1,2)))
# Defining the sankey diagram
  p <- ggplot(I_new, aes(x, id = id, split = y, value = value)) +
    geom_parallel_sets(aes(fill=End), alpha = 0.7, axis.width = 0.2) +
    geom_parallel_sets_axes(fill = "white", colour = "light grey", axis.width = 0.2) +
    geom_parallel_sets_labels( angle = 0) +
    scale_fill_manual(values = pal, na.translate = FALSE) +
    scale_colour_manual(values = pal, na.translate = FALSE) +
    theme_minimal(base_size = 14) +
    theme(axis.title.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none") +
    coord_cartesian(expand = FALSE, clip = "off") + 
    scale_x_discrete( labels = c("New patients", "Initial outcome", "Final outcome")
    )
  return(p)
}

staged_sankey_inset <- function(sankey_data)
{
  # Defining colour palette
  pal <- c("      Cure" = '#A9A9A9',
           "     DS"    = '#AE0D0A',
           "     RR"    = '#177E89',
           "    BR"     = '#CBA715',
           "    XR"     = '#A5C994', 
           "   RR/BR"   = '#CBA715',
           "   RR/XR"   = '#A5C994',
           "  BR/XR"    = '#754668',
           " RR/BR/XR"  = '#754668',
           "Death"         = '#030303')
  # Defining states and data frame
  start_states <- c("     DS", "     RR", "    BR", "    XR", "   RR/BR", "   RR/XR", "  BR/XR", " RR/BR/XR")
  mid_states <- c("      Cure", "     DS", "     RR", "    BR", "    XR", "   RR/BR", "   RR/XR", "  BR/XR", " RR/BR/XR","Death")
  end_states <- c("      Cure", "     DS", "     RR", "    BR", "    XR", "   RR/BR", "   RR/XR", "  BR/XR", " RR/BR/XR","Death")
  I <- data.frame(Start = rep(c("     DS", "     RR", "    BR", "    XR", "   RR/BR", "   RR/XR", "  BR/XR", " RR/BR/XR"),each=100),
                  Mid = rep(c("      Cure", "     DS", "     RR", "    BR", "    XR", "   RR/BR", "   RR/XR", "  BR/XR", " RR/BR/XR","Death"),each = 10),
                  End = rep(c("      Cure", "     DS", "     RR", "    BR", "    XR", "   RR/BR", "   RR/XR", "  BR/XR", " RR/BR/XR","Death"), 80),
                  value = rep(0,800))
  # Populating frequency table  
  exc <- c(sankey_data$DS,sankey_data$RR,sankey_data$BR,sankey_data$XR,sankey_data$RRBR,sankey_data$RRXR,sankey_data$BRXR,sankey_data$RRBRXR)
  I$value <- exc
  # Add in time points and colours
  I_new <- I %>%
    gather_set_data(1:3) %>%
    mutate(y = factor(y, levels = end_states)) %>%
    mutate(x = factor(x, levels = c(1,2,3)))
  (inset <- ggplot(I_new %>% filter(!str_detect(Mid,"Cure"), !str_detect(Mid,"Death"), x!=1) %>% 
                     select(-Start), aes(x, id = id, split = y, value = value)) +
      geom_parallel_sets(aes(fill=End), alpha = 0.7, axis.width = 0.2) +
      geom_parallel_sets_axes(fill = "white", colour = "light grey", axis.width = 0.2) +
      geom_parallel_sets_labels( angle = 0) +
      scale_fill_manual(values = pal, na.translate = FALSE) +
      scale_colour_manual(values = pal, na.translate = FALSE) +
      theme_minimal(base_size = 14) +
      theme(axis.title.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "none") +
      coord_cartesian(expand = FALSE, clip = "off") + 
      scale_x_discrete( labels = c("Initial outcome", "Final outcome")) 
  )
  return(list("figure"=inset, 
              "yscale"=sum(I_new$value[!str_detect(I_new$Mid,"Cure") & !str_detect(I_new$Mid,"Death") & I_new$x==2])/
                sum(I_new$value[I_new$x==2]))) # denominator ins't really necessary, though shows rounding effects
}
## ================ Function for SoC pathway ===================================
 cohort_soc <- function(prev_data, param_set) {
   outcome_set <- outcomes(param_set)
   ## NEW PATIENTS
   # Proportion with Rho DST
   rho_y     <- prev_data * param_set$rho_soc_new
   rho_n     <- prev_data - rho_y
   # Proportion with Beta/Chi DST, conditional upon existing RR
   betachi_y <- rho_y * param_set$betachi_soc_new * c(0, 1, 0, 0, 1, 1, 0, 1)
   betachi_n <- rho_y - betachi_y
   # Outcomes for each resistance phenotype, conditional upon treatment pathway
   ds_out     <- rho_n$prev_DS * outcome_set["DS",,"R"]         + betachi_n$prev_DS * outcome_set["DS",,"R"]          + betachi_y$prev_DS * outcome_set["DS",,"R"] 
   rr_out     <- rho_n$prev_RR * outcome_set["RR",,"R"]         + betachi_n$prev_RR * outcome_set["RR",,"BX"]         + betachi_y$prev_RR * outcome_set["RR",,"BX"] 
   br_out     <- rho_n$prev_BR * outcome_set["BR",,"R"]         + betachi_n$prev_BR * outcome_set["BR",,"R"]          + betachi_y$prev_BR * outcome_set["BR",,"R"] 
   xr_out     <- rho_n$prev_XR * outcome_set["XR",,"R"]         + betachi_n$prev_XR * outcome_set["XR",,"R"]          + betachi_y$prev_XR * outcome_set["XR",,"R"] 
   rrbr_out   <- rho_n$prev_RRBR * outcome_set["RRBR",,"R"]     + betachi_n$prev_RRBR * outcome_set["RRBR",,"BX"]     + betachi_y$prev_RRBR * outcome_set["RRBR",,"Xbased"]
   rrxr_out   <- rho_n$prev_RRXR * outcome_set["RRXR",,"R"]     + betachi_n$prev_RRXR * outcome_set["RRXR",,"BX"]     + betachi_y$prev_RRXR * outcome_set["RRXR",,"Bbased"]
   brxr_out   <- rho_n$prev_BRXR * outcome_set["BRXR",,"R"]     + betachi_n$prev_BRXR * outcome_set["BRXR",,"R"]      + betachi_y$prev_BRXR * outcome_set["BRXR",,"R"]
   rrbrxr_out <- rho_n$prev_RRBRXR * outcome_set["RRBRXR",,"R"] + betachi_n$prev_RRBRXR * outcome_set["RRBRXR",,"BX"] + betachi_y$prev_RRBRXR * outcome_set["RRBRXR",,"Ind"]
   total_out_new  <- rbind(ds_out, rr_out, br_out, xr_out, rrbr_out, rrxr_out, brxr_out, rrbrxr_out)
   # Distribution of resistance in those who failed treatment, taking acquisition into account
   prev_out   <- total_out_new[,"Failure"] - total_out_new[,"Acquired RR"] - total_out_new[,"Acquired BR"] - total_out_new[,"Acquired XR"]
   prev_out   <- prev_out + c(0,
                              ds_out["Acquired RR"],
                              ds_out["Acquired BR"],
                              ds_out["Acquired XR"],
                              rr_out["Acquired BR"] + br_out["Acquired RR"],
                              rr_out["Acquired XR"] + xr_out["Acquired RR"],
                              br_out["Acquired XR"] + xr_out["Acquired BR"],
                              rrbr_out["Acquired XR"] + rrxr_out["Acquired BR"] + brxr_out["Acquired RR"]) 
   prev_out   <- as.data.frame.list(prev_out)
   prev_out_new <- prev_out
   colnames(prev_out) <- c('prev_DS', 'prev_RR', 'prev_BR', 'prev_XR', 'prev_RRBR', 'prev_RRXR', 'prev_BRXR', 'prev_RRBRXR')
   # Defining the sankey table
   soc_new_sankey <- data.frame(DS     = c(ds_out["Cure"], ds_out["Failure"]-ds_out["Acquired RR"]-ds_out["Acquired BR"]-ds_out["Acquired XR"], ds_out["Acquired RR"], ds_out["Acquired BR"], ds_out["Acquired XR"], 0, 0, 0, 0, ds_out["Death"]),
                              RR     = c(rr_out["Cure"], 0, rr_out["Failure"]-rr_out["Acquired BR"]-rr_out["Acquired XR"], 0, 0, rr_out["Acquired BR"], rr_out["Acquired XR"], 0, 0, rr_out["Death"]),
                              BR     = c(br_out["Cure"], 0, 0, br_out["Failure"]-br_out["Acquired RR"]-br_out["Acquired XR"], 0, br_out["Acquired RR"], 0, br_out["Acquired XR"], 0, br_out["Death"]),
                              XR     = c(xr_out["Cure"], 0, 0, 0, xr_out["Failure"]-xr_out["Acquired RR"]-xr_out["Acquired BR"], 0, xr_out["Acquired RR"], xr_out["Acquired BR"], 0, xr_out["Death"]),
                              RRBR   = c(rrbr_out["Cure"], 0, 0, 0, 0, rrbr_out["Failure"]-rrbr_out["Acquired XR"], 0, 0, rrbr_out["Acquired XR"], rrbr_out["Death"]),
                              RRXR   = c(rrxr_out["Cure"], 0, 0, 0, 0, 0, rrxr_out["Failure"]-rrxr_out["Acquired BR"], 0, rrxr_out["Acquired BR"], rrxr_out["Death"]),
                              BRXR   = c(brxr_out["Cure"], 0, 0, 0, 0, 0, 0, brxr_out["Failure"]-brxr_out["Acquired RR"], brxr_out["Acquired RR"], brxr_out["Death"]),
                              RRBRXR = c(rrbrxr_out["Cure"], 0, 0, 0, 0, 0, 0, 0, rrbrxr_out["Failure"], rrbrxr_out["Death"]))
   ## PREVIOUSLY TREATED PATIENTS
   # Splitting out those previously treated with R vs those previously treated with BX
   prev_R <- data.frame(prev_DS = ds_out["Failure"] - ds_out["Acquired RR"], 
                        prev_RR = ds_out["Acquired RR"] + rho_n$prev_RR * outcome_set["RR","Failure","R"],
                        prev_BR = br_out["Failure"] - br_out["Acquired RR"],
                        prev_XR = xr_out["Failure"] - xr_out["Acquired RR"],
                        prev_RRBR = rho_n$prev_RRBR * outcome_set["RRBR","Failure","R"]+br_out["Acquired RR"],
                        prev_RRXR = rho_n$prev_RRXR * outcome_set["RRXR","Failure","R"]+xr_out["Acquired RR"],
                        prev_BRXR = brxr_out["Failure"] - brxr_out["Acquired RR"],
                        prev_RRBRXR = rrbrxr_out["Failure"]+brxr_out["Acquired RR"])
   prev_BX <- data.frame(prev_DS = 0, 
                        prev_RR = betachi_n$prev_RR * outcome_set["RR","Failure","BX"] + betachi_y$prev_RR * outcome_set["RR","Failure","BX"] - rr_out["Acquired BR"] - rr_out["Acquired XR"],
                        prev_BR = 0,
                        prev_XR = 0,
                        prev_RRBR = betachi_n$prev_RRBR * outcome_set["RRBR","Failure","BX"] + betachi_y$prev_RRBR * outcome_set["RRBR","Failure","Xbased"] - rrbr_out["Acquired XR"] + betachi_n$prev_RR * outcome_set["RR","Acquired BR","BX"] + betachi_y$prev_RR * outcome_set["RR","Acquired BR","BX"],
                        prev_RRXR = betachi_n$prev_RRXR * outcome_set["RRXR","Failure","BX"] + betachi_y$prev_RRXR * outcome_set["RRXR","Failure","Bbased"] - rrxr_out["Acquired BR"] + betachi_n$prev_RR * outcome_set["RR","Acquired XR","BX"] + betachi_y$prev_RR * outcome_set["RR","Acquired XR","BX"],
                        prev_BRXR = 0,
                        prev_RRBRXR = betachi_n$prev_RRBRXR * outcome_set["RRBRXR","Failure","BX"] + betachi_y$prev_RRBRXR * outcome_set["RRBRXR","Failure","Ind"] + rrbr_out["Acquired XR"] + rrxr_out["Acquired BR"])
   # Proportion with Rho DST
   rho_y     <- prev_R * param_set$rho_soc_retR + prev_BX * param_set$rho_soc_retBX
   rho_n     <- prev_out - rho_y
   # Proportion with Beta/Chi DST, conditional upon existing RR
   betachi_y <- (prev_R * param_set$rho_soc_retR * param_set$betachi_soc_retR  + prev_BX * param_set$rho_soc_retBX * param_set$betachi_soc_retBX ) * c(0, 1, 0, 0, 1, 1, 0, 1)
   betachi_n <- rho_y - betachi_y
   # Outcomes for each resistance phenotype, conditional upon treatment pathway
   ds_out     <- rho_n$prev_DS * outcome_set["DS",,"R"]         + betachi_n$prev_DS * outcome_set["DS",,"R"]          + betachi_y$prev_DS * outcome_set["DS",,"R"] 
   rr_out     <- rho_n$prev_RR * outcome_set["RR",,"R"]         + betachi_n$prev_RR * outcome_set["RR",,"BX"]         + betachi_y$prev_RR * outcome_set["RR",,"BX"] 
   br_out     <- rho_n$prev_BR * outcome_set["BR",,"R"]         + betachi_n$prev_BR * outcome_set["BR",,"R"]          + betachi_y$prev_BR * outcome_set["BR",,"R"] 
   xr_out     <- rho_n$prev_XR * outcome_set["XR",,"R"]         + betachi_n$prev_XR * outcome_set["XR",,"R"]          + betachi_y$prev_XR * outcome_set["XR",,"R"] 
   rrbr_out   <- rho_n$prev_RRBR * outcome_set["RRBR",,"R"]     + betachi_n$prev_RRBR * outcome_set["RRBR",,"BX"]     + betachi_y$prev_RRBR * outcome_set["RRBR",,"Xbased"]
   rrxr_out   <- rho_n$prev_RRXR * outcome_set["RRXR",,"R"]     + betachi_n$prev_RRXR * outcome_set["RRXR",,"BX"]     + betachi_y$prev_RRXR * outcome_set["RRXR",,"Bbased"]
   brxr_out   <- rho_n$prev_BRXR * outcome_set["BRXR",,"R"]     + betachi_n$prev_BRXR * outcome_set["BRXR",,"R"]      + betachi_y$prev_BRXR * outcome_set["BRXR",,"R"]
   rrbrxr_out <- rho_n$prev_RRBRXR * outcome_set["RRBRXR",,"R"] + betachi_n$prev_RRBRXR * outcome_set["RRBRXR",,"BX"] + betachi_y$prev_RRBRXR * outcome_set["RRBRXR",,"Ind"]
   total_out_ret  <- rbind(ds_out, rr_out, br_out, xr_out, rrbr_out, rrxr_out, brxr_out, rrbrxr_out)
   # Distribution of resistance in those who failed treatment, taking acquisition into account
   prev_out   <- total_out_ret[,"Failure"] - total_out_ret[,"Acquired RR"] - total_out_ret[,"Acquired BR"] - total_out_ret[,"Acquired XR"]
   prev_out   <- prev_out + c(0,
                              ds_out["Acquired RR"],
                              ds_out["Acquired BR"],
                              ds_out["Acquired XR"],
                              rr_out["Acquired BR"] + br_out["Acquired RR"],
                              rr_out["Acquired XR"] + xr_out["Acquired RR"],
                              br_out["Acquired XR"] + xr_out["Acquired BR"],
                              rrbr_out["Acquired XR"] + rrxr_out["Acquired BR"] + brxr_out["Acquired RR"]) 
   # Removing treatment failures in new patients as they go through the retreatment pathway
   total_out_new[,"Failure"] <- 0
   total_out  <- total_out_new + total_out_ret
   # Defining the sankey outcome table
   soc_ret_sankey <- data.frame(DS     = c(ds_out["Cure"], ds_out["Failure"]-ds_out["Acquired RR"]-ds_out["Acquired BR"]-ds_out["Acquired XR"], ds_out["Acquired RR"], ds_out["Acquired BR"], ds_out["Acquired XR"], 0, 0, 0, 0, ds_out["Death"]),
                                RR     = c(rr_out["Cure"], 0, rr_out["Failure"]-rr_out["Acquired BR"]-rr_out["Acquired XR"], 0, 0, rr_out["Acquired BR"], rr_out["Acquired XR"], 0, 0, rr_out["Death"]),
                                BR     = c(br_out["Cure"], 0, 0, br_out["Failure"]-br_out["Acquired RR"]-br_out["Acquired XR"], 0, br_out["Acquired RR"], 0, br_out["Acquired XR"], 0, br_out["Death"]),
                                XR     = c(xr_out["Cure"], 0, 0, 0, xr_out["Failure"]-xr_out["Acquired RR"]-xr_out["Acquired BR"], 0, xr_out["Acquired RR"], xr_out["Acquired BR"], 0, xr_out["Death"]),
                                RRBR   = c(rrbr_out["Cure"], 0, 0, 0, 0, rrbr_out["Failure"]-rrbr_out["Acquired XR"], 0, 0, rrbr_out["Acquired XR"], rrbr_out["Death"]),
                                RRXR   = c(rrxr_out["Cure"], 0, 0, 0, 0, 0, rrxr_out["Failure"]-rrxr_out["Acquired BR"], 0, rrxr_out["Acquired BR"], rrxr_out["Death"]),
                                BRXR   = c(brxr_out["Cure"], 0, 0, 0, 0, 0, 0, brxr_out["Failure"]-brxr_out["Acquired RR"], brxr_out["Acquired RR"], brxr_out["Death"]),
                                RRBRXR = c(rrbrxr_out["Cure"], 0, 0, 0, 0, 0, 0, 0, rrbrxr_out["Failure"], rrbrxr_out["Death"]))
   # Combining sankey tables, each row is a set out final outcomes (Cure, DS, RR, BR, XR, RRBR, RRXR, BRXR, RRBRXR, death)
   soc_tot_sankey <- data.frame(DS = c(soc_new_sankey[1,"DS"], 0, 0, 0, 0, 0, 0, 0, 0, 0, # Starting=DS, after 1 round of treatment=treatment cure
                                       soc_ret_sankey[,"DS"], # Starting=DS, after 1 round of treatment=DS
                                       soc_new_sankey[3,"DS"]/sum(soc_new_sankey[3,]) * soc_ret_sankey[,"RR"], # Starting=DS, after 1 round of treatment=RR
                                       soc_new_sankey[4,"DS"]/sum(soc_new_sankey[4,]) * soc_ret_sankey[,"BR"], # Starting=DS, after 1 round of treatment=BR
                                       soc_new_sankey[5,"DS"]/sum(soc_new_sankey[5,]) * soc_ret_sankey[,"XR"], # # Starting=DS, after 1 round of treatment=XR
                                       rep(0,10), # Starting=DS, after 1 round of treatment=RRBR
                                       rep(0,10), # Starting=DS, after 1 round of treatment=RRXR
                                       rep(0,10), # Starting=DS, after 1 round of treatment=BRXR
                                       rep(0,10), # Starting=DS, after 1 round of treatment=RRBRXR
                                       0, 0, 0, 0, 0, 0, 0, 0, 0, soc_new_sankey[10,"DS"]), # Starting=DS, after 1 round of treatment=death
                                # Starting RR
                                RR = c(soc_new_sankey[1,"RR"], 0, 0, 0, 0, 0, 0, 0, 0, 0, # Starting=RR, after 1 round of treatment=treatment cure
                                       rep(0,10), # Starting=RR, after 1 round of treatment=DS
                                       soc_new_sankey[3,"RR"]/sum(soc_new_sankey[3,])*soc_ret_sankey[,"RR"], # Starting=RR, after 1 round of treatment=RR
                                       rep(0,10), # Starting=RR, after 1 round of treatment=BR
                                       rep(0,10), # Starting=RR, after 1 round of treatment=XR
                                       soc_new_sankey[6,"RR"]/sum(soc_new_sankey[6,])*soc_ret_sankey[,"RRBR"], # Starting=RR, after 1 round of treatment=RRBR
                                       soc_new_sankey[7,"RR"]/sum(soc_new_sankey[7,])*soc_ret_sankey[,"RRXR"], # Starting=RR, after 1 round of treatment=RRXR
                                       rep(0,10), # Starting=RR, after 1 round of treatment=BRXR
                                       rep(0,10), # Starting=RR, after 1 round of treatment=RRBRXR
                                       0, 0, 0, 0, 0, 0, 0, 0, 0, soc_new_sankey[10,"RR"]), # Starting=RR, after 1 round of treatment=Death
                                #Starting BR
                                BR = c(soc_new_sankey[1,"BR"], 0, 0, 0, 0, 0, 0, 0, 0, 0, # Starting=BR, after 1 round of treatment=treatment cure
                                       rep(0,10), # Starting=BR, after 1 round of treatment=DS
                                       rep(0,10), # Starting=BR, after 1 round of treatment=RR
                                       soc_new_sankey[4,"BR"]/sum(soc_new_sankey[4,])*soc_ret_sankey[,"BR"], # Starting=BR, after 1 round of treatment=BR
                                       rep(0,10), # Starting=BR, after 1 round of treatment=RRXR
                                       soc_new_sankey[6,"BR"]/sum(soc_new_sankey[6,])*soc_ret_sankey[,"RRBR"], # Starting=BR, after 1 round of treatment=RRBR
                                       rep(0,10), # Starting=BR, after 1 round of treatment=RRXR
                                       soc_new_sankey[8,"BR"]/sum(soc_new_sankey[8,])*soc_ret_sankey[,"BRXR"], # Starting=BR, after 1 round of treatment=BRXR
                                       rep(0,10), # Starting=BR, after 1 round of treatment=RRBRXR
                                       0, 0, 0, 0, 0, 0, 0, 0, 0, soc_new_sankey[10,"BR"]), # Starting=BR, after 1 round of treatment=Death)
                                #Starting XR
                                XR = c(soc_new_sankey[1,"XR"], 0, 0, 0, 0, 0, 0, 0, 0, 0, # Starting=BR, after 1 round of treatment=treatment cure
                                       rep(0,10), # Starting=XR, after 1 round of treatment=DS
                                       rep(0,10), # Starting=XR, after 1 round of treatment=RR
                                       rep(0,10), # Starting=XR, after 1 round of treatment=BR
                                       soc_new_sankey[5,"XR"]/sum(soc_new_sankey[5,])*soc_ret_sankey[,"XR"], # Starting=XR, after 1 round of treatment=XR
                                       rep(0,10), # Starting=XR, after 1 round of treatment=RRBR
                                       soc_new_sankey[7,"XR"]/sum(soc_new_sankey[7,])*soc_ret_sankey[,"RRXR"], # Starting=XR, after 1 round of treatment=RRXR
                                       soc_new_sankey[8,"XR"]/sum(soc_new_sankey[8,])*soc_ret_sankey[,"BRXR"], # Starting=XR, after 1 round of treatment=BRXR
                                       rep(0,10), # Starting=XR, after 1 round of treatment=RRBRXR
                                       0, 0, 0, 0, 0, 0, 0, 0, 0, soc_new_sankey[10,"XR"]), # Starting=XR, after 1 round of treatment=Death)   
                                # Starting RRBR
                                RRBR = c(soc_new_sankey[1,"RRBR"], 0, 0, 0, 0, 0, 0, 0, 0, 0, # Starting=RRBR, after 1 round of treatment=treatment cure
                                         rep(0,10), # Starting=RRBR, after 1 round of treatment=DS
                                         rep(0,10), # Starting=RRBR, after 1 round of treatment=RR
                                         rep(0,10), # Starting=RRBR, after 1 round of treatment=BR
                                         rep(0,10), # Starting=RRBR, after 1 round of treatment=XR
                                         soc_new_sankey[6,"RRBR"]/sum(soc_new_sankey[6,])*soc_ret_sankey[,"RRBR"], # Starting=RRBR, after 1 round of treatment=RRBR
                                         rep(0,10), # Starting=RRBR, after 1 round of treatment=RRXR
                                         rep(0,10), # Starting=RRBR, after 1 round of treatment=BRXR
                                         soc_new_sankey[9,"RRBR"]/sum(soc_new_sankey[9,])*soc_ret_sankey[,"RRBRXR"], # Starting=RRBR, after 1 round of treatment=RRBRXR
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, soc_new_sankey[10,"BR"]), # Starting=RRBR, after 1 round of treatment=Death
                                # Starting RRXR
                                RRXR = c(soc_new_sankey[1,"RRXR"], 0, 0, 0, 0, 0, 0, 0, 0, 0, # Starting=XRBR, after 1 round of treatment=treatment cure
                                         rep(0,10), # Starting=RRXR, after 1 round of treatment=DS
                                         rep(0,10), # Starting=RRXR, after 1 round of treatment=RR
                                         rep(0,10), # Starting=RRXR, after 1 round of treatment=BR
                                         rep(0,10), # Starting=RRXR, after 1 round of treatment=XR
                                         rep(0,10), # Starting=RRXR, after 1 round of treatment=RRBR
                                         soc_new_sankey[7,"RRXR"]/sum(soc_new_sankey[7,])*soc_ret_sankey[,"RRXR"], # Starting=RRXR, after 1 round of treatment=RRXR
                                         rep(0,10), # Starting=RRXR, after 1 round of treatment=BRXR
                                         soc_new_sankey[9,"RRXR"]/sum(soc_new_sankey[9,])*soc_ret_sankey[,"RRBRXR"], # Starting=RRXR, after 1 round of treatment=RRBRXR
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, soc_new_sankey[10,"XR"]), # Starting=RRXR, after 1 round of treatment=Death
                                # Starting BRXR
                                BRXR = c(soc_new_sankey[1,"BRXR"], 0, 0, 0, 0, 0, 0, 0, 0, 0, # Starting=BRXR, after 1 round of treatment=treatment cure
                                         rep(0,10), # Starting=BRXR, after 1 round of treatment=DS
                                         rep(0,10), # Starting=BRXR, after 1 round of treatment=RR
                                         rep(0,10), # Starting=BRXR, after 1 round of treatment=BR
                                         rep(0,10), # Starting=BRXR, after 1 round of treatment=XR
                                         rep(0,10), # Starting=BRXR, after 1 round of treatment=RRBR
                                         rep(0,10), # Starting=BRXR, after 1 round of treatment=RRXR
                                         soc_new_sankey[8,"BRXR"]/sum(soc_new_sankey[8,])*soc_ret_sankey[,"RRBRXR"], # Starting=BRXR, after 1 round of treatment=BRXR
                                         soc_new_sankey[9,"BRXR"]/sum(soc_new_sankey[9,])*soc_ret_sankey[,"RRBRXR"], # Starting=BRXR, after 1 round of treatment=RRBRXR
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, soc_new_sankey[10,"BRXR"]), # Starting=BRXR, after 1 round of treatment=Death
                                # Starting RRBRXR
                                RRBRXR = c(soc_new_sankey[1,"RRBRXR"], 0, 0, 0, 0, 0, 0, 0, 0, 0, # Starting=RRBRXR, after 1 round of treatment=treatment cure
                                           rep(0,10), # Starting=RRBRXR, after 1 round of treatment=DS
                                           rep(0,10), # Starting=RRBRXR, after 1 round of treatment=RR
                                           rep(0,10), # Starting=RRBRXR, after 1 round of treatment=BR
                                           rep(0,10), # Starting=RRBRXR, after 1 round of treatment=XR
                                           rep(0,10), # Starting=RRBRXR, after 1 round of treatment=RRBR
                                           rep(0,10), # Starting=RRBRXR, after 1 round of treatment=RRXR
                                           rep(0,10), # Starting=RRBRXR, after 1 round of treatment=BRXR
                                           soc_new_sankey[9,"RRBRXR"]/sum(soc_new_sankey[9,])*soc_ret_sankey[,"RRBRXR"], # Starting=RRBRXR, after 1 round of treatment=RRBRXR
                                           0, 0, 0, 0, 0, 0, 0, 0, 0, soc_new_sankey[10,"RRBRXR"])) # Starting=RRBRXR, after 1 round of treatment=Death                                
   soc_tot_sankey[is.na(soc_tot_sankey)] <- 0
   tmp_out <-prev_out_new+prev_out
   tmp_out$cureN<-sum(total_out_new[,"Cure"])
   tmp_out$cureF<-sum(total_out[,"Cure"])
   tmp_out$death<-sum(total_out[,"Death"])
   
   return(soc_tot_sankey)
   #HEATMAP OUTPUTS: c(new=sum(total_out_new[,"Cure"]), final=sum(total_out[,"Cure"]), death=sum(total_out[,"Death"]), resist=prev_out[5]+prev_out[6]+prev_out[8])
   #SANKEY: soc_tot_sankey
   #PREVALENCE OVER TIME: prev_out_new+prev_out 
   #SENSITIVITY OVER TIME: tmp_out
   }
 
## ================ Function for Pan-TB pathway ================================
 cohort_pan <- function(prev_data, param_set) {
   outcome_set <- outcomes(param_set)
   ## NEW PATIENTS
   # Outcomes for each resistance phenotype
   ds_out     <- prev_data$prev_DS * outcome_set["DS",,"BX"]
   rr_out     <- prev_data$prev_RR * outcome_set["RR",,"BX"]
   br_out     <- prev_data$prev_BR * outcome_set["BR",,"BX"]
   xr_out     <- prev_data$prev_XR * outcome_set["XR",,"BX"]
   rrbr_out   <- prev_data$prev_RRBR * outcome_set["RRBR",,"BX"] 
   rrxr_out   <- prev_data$prev_RRXR * outcome_set["RRXR",,"BX"] 
   brxr_out   <- prev_data$prev_BRXR * outcome_set["BRXR",,"BX"]  
   rrbrxr_out <- prev_data$prev_RRBRXR * outcome_set["RRBRXR",,"BX"] 
   total_out_new  <- rbind(ds_out, rr_out, br_out, xr_out, rrbr_out, rrxr_out, brxr_out, rrbrxr_out)
   # Distribution of resistance in those who failed treatment, taking acquisition into account
   prev_out   <- total_out_new[,"Failure"] - total_out_new[,"Acquired RR"] - total_out_new[,"Acquired BR"] - total_out_new[,"Acquired XR"]
   prev_out   <- prev_out + c(0,
                              ds_out["Acquired RR"],
                              ds_out["Acquired BR"],
                              ds_out["Acquired XR"],
                              rr_out["Acquired BR"] + br_out["Acquired RR"],
                              rr_out["Acquired XR"] + xr_out["Acquired RR"],
                              br_out["Acquired XR"] + xr_out["Acquired BR"],
                              rrbr_out["Acquired XR"] + rrxr_out["Acquired BR"] + brxr_out["Acquired RR"])
   prev_out   <- as.data.frame.list(prev_out)
   prev_out_new <- prev_out
   colnames(prev_out) <- c('prev_DS', 'prev_RR', 'prev_BR', 'prev_XR', 'prev_RRBR', 'prev_RRXR', 'prev_BRXR', 'prev_RRBRXR')
   # Defining the sankey table
   pan_new_sankey <- data.frame(DS   = c(ds_out["Cure"], ds_out["Failure"]-ds_out["Acquired RR"]-ds_out["Acquired BR"]-ds_out["Acquired XR"], ds_out["Acquired RR"], ds_out["Acquired BR"], ds_out["Acquired XR"], 0, 0, 0, 0, ds_out["Death"]),
                              RR     = c(rr_out["Cure"], 0, rr_out["Failure"]-rr_out["Acquired BR"]-rr_out["Acquired XR"], 0, 0, rr_out["Acquired BR"], rr_out["Acquired XR"], 0, 0, rr_out["Death"]),
                              BR     = c(br_out["Cure"], 0, 0, br_out["Failure"]-br_out["Acquired RR"]-br_out["Acquired XR"], 0, br_out["Acquired RR"], 0, br_out["Acquired XR"], 0, br_out["Death"]),
                              XR     = c(xr_out["Cure"], 0, 0, 0, xr_out["Failure"]-xr_out["Acquired RR"]-xr_out["Acquired BR"], 0, xr_out["Acquired RR"], xr_out["Acquired BR"], 0, xr_out["Death"]),
                              RRBR   = c(rrbr_out["Cure"], 0, 0, 0, 0, rrbr_out["Failure"]-rrbr_out["Acquired XR"], 0, 0, rrbr_out["Acquired XR"], rrbr_out["Death"]),
                              RRXR   = c(rrxr_out["Cure"], 0, 0, 0, 0, 0, rrxr_out["Failure"]-rrxr_out["Acquired BR"], 0, rrxr_out["Acquired BR"], rrxr_out["Death"]),
                              BRXR   = c(brxr_out["Cure"], 0, 0, 0, 0, 0, 0, brxr_out["Failure"]-brxr_out["Acquired RR"], brxr_out["Acquired RR"], brxr_out["Death"]),
                              RRBRXR = c(rrbrxr_out["Cure"], 0, 0, 0, 0, 0, 0, 0, rrbrxr_out["Failure"], rrbrxr_out["Death"]))
   ## PREVIOUSLY TREATED PATIENTS
   # Proportion with Rho DST
   rho_y     <- prev_out * param_set$rho_pan
   rho_n     <- prev_out - rho_y
   # Proportion with Beta/Chi DST, conditional upon existing RR
   betachi_y <- (prev_out * param_set$rho_pan * param_set$betachi_pan) * c(0, 1, 0, 0, 1, 1, 0, 1)
   betachi_n <- rho_y - betachi_y
   # Outcomes for each resistance phenotype, conditional upon treatment pathway
   ds_out     <- rho_n$prev_DS * outcome_set["DS",,"BX"]         + betachi_n$prev_DS * outcome_set["DS",,"R"]          + betachi_y$prev_DS * outcome_set["DS",,"R"]
   rr_out     <- rho_n$prev_RR * outcome_set["RR",,"BX"]         + betachi_n$prev_RR * outcome_set["RR",,"BX"]         + betachi_y$prev_RR * outcome_set["RR",,"BX"]
   br_out     <- rho_n$prev_BR * outcome_set["BR",,"BX"]         + betachi_n$prev_BR * outcome_set["BR",,"R"]          + betachi_y$prev_BR * outcome_set["BR",,"R"]
   xr_out     <- rho_n$prev_XR * outcome_set["XR",,"BX"]         + betachi_n$prev_XR * outcome_set["XR",,"R"]          + betachi_y$prev_XR * outcome_set["XR",,"R"]
   rrbr_out   <- rho_n$prev_RRBR * outcome_set["RRBR",,"BX"]     + betachi_n$prev_RRBR * outcome_set["RRBR",,"BX"]     + betachi_y$prev_RRBR * outcome_set["RRBR",,"Xbased"]
   rrxr_out   <- rho_n$prev_RRXR * outcome_set["RRXR",,"BX"]     + betachi_n$prev_RRXR * outcome_set["RRXR",,"BX"]     + betachi_y$prev_RRXR * outcome_set["RRXR",,"Bbased"]
   brxr_out   <- rho_n$prev_BRXR * outcome_set["BRXR",,"BX"]     + betachi_n$prev_BRXR * outcome_set["BRXR",,"R"]      + betachi_y$prev_BRXR * outcome_set["BRXR",,"R"]
   rrbrxr_out <- rho_n$prev_RRBRXR * outcome_set["RRBRXR",,"BX"] + betachi_n$prev_RRBRXR * outcome_set["RRBRXR",,"BX"] + betachi_y$prev_RRBRXR * outcome_set["RRBRXR",,"Ind"]
   total_out_ret  <- rbind(ds_out, rr_out, br_out, xr_out, rrbr_out, rrxr_out, brxr_out, rrbrxr_out)
   # Distribution of resistance in those who failed treatment, taking acquisition into account
   prev_out   <- total_out_ret[,"Failure"] - total_out_ret[,"Acquired RR"] - total_out_ret[,"Acquired BR"] - total_out_ret[,"Acquired XR"]
   prev_out   <- prev_out + c(0,
                              ds_out["Acquired RR"],
                              ds_out["Acquired BR"],
                              ds_out["Acquired XR"],
                              rr_out["Acquired BR"] + br_out["Acquired RR"],
                              rr_out["Acquired XR"] + xr_out["Acquired RR"],
                              br_out["Acquired XR"] + xr_out["Acquired BR"],
                              rrbr_out["Acquired XR"] + rrxr_out["Acquired BR"] + brxr_out["Acquired RR"])
   # Removing treatment failures in new patients as they go through the retreatment pathway
   total_out_new[,"Failure"] <- 0
   total_out  <- total_out_new + total_out_ret
   # Defining the sankey outcome table
   pan_ret_sankey <- data.frame(DS     = c(ds_out["Cure"], ds_out["Failure"]-ds_out["Acquired RR"]-ds_out["Acquired BR"]-ds_out["Acquired XR"], ds_out["Acquired RR"], ds_out["Acquired BR"], ds_out["Acquired XR"], 0, 0, 0, 0, ds_out["Death"]),
                                RR     = c(rr_out["Cure"], 0, rr_out["Failure"]-rr_out["Acquired BR"]-rr_out["Acquired XR"], 0, 0, rr_out["Acquired BR"], rr_out["Acquired XR"], 0, 0, rr_out["Death"]),
                                BR     = c(br_out["Cure"], 0, 0, br_out["Failure"]-br_out["Acquired RR"]-br_out["Acquired XR"], 0, br_out["Acquired RR"], 0, br_out["Acquired XR"], 0, br_out["Death"]),
                                XR     = c(xr_out["Cure"], 0, 0, 0, xr_out["Failure"]-xr_out["Acquired RR"]-xr_out["Acquired BR"], 0, xr_out["Acquired RR"], xr_out["Acquired BR"], 0, xr_out["Death"]),
                                RRBR   = c(rrbr_out["Cure"], 0, 0, 0, 0, rrbr_out["Failure"]-rrbr_out["Acquired XR"], 0, 0, rrbr_out["Acquired XR"], rrbr_out["Death"]),
                                RRXR   = c(rrxr_out["Cure"], 0, 0, 0, 0, 0, rrxr_out["Failure"]-rrxr_out["Acquired BR"], 0, rrxr_out["Acquired BR"], rrxr_out["Death"]),
                                BRXR   = c(brxr_out["Cure"], 0, 0, 0, 0, 0, 0, brxr_out["Failure"]-brxr_out["Acquired RR"], brxr_out["Acquired RR"], brxr_out["Death"]),
                                RRBRXR = c(rrbrxr_out["Cure"], 0, 0, 0, 0, 0, 0, 0, rrbrxr_out["Failure"], rrbrxr_out["Death"]))
   # Combining sankey tables, each row is a set out final outcomes (Cure, DS, RR, BR, XR, RRBR, RRXR, BRXR, RRBRXR, death)
   pan_tot_sankey <- data.frame(DS = c(pan_new_sankey[1,"DS"], 0, 0, 0, 0, 0, 0, 0, 0, 0, # Starting=DS, after 1 round of treatment=treatment cure
                                pan_ret_sankey[,"DS"], # Starting=DS, after 1 round of treatment=DS
                                pan_new_sankey[3,"DS"]/sum(pan_new_sankey[3,]) * pan_ret_sankey[,"RR"], # Starting=DS, after 1 round of treatment=RR
                                pan_new_sankey[4,"DS"]/sum(pan_new_sankey[4,]) * pan_ret_sankey[,"BR"], # Starting=DS, after 1 round of treatment=BR
                                pan_new_sankey[5,"DS"]/sum(pan_new_sankey[5,]) * pan_ret_sankey[,"XR"], # # Starting=DS, after 1 round of treatment=XR
                                rep(0,10), # Starting=DS, after 1 round of treatment=RRBR
                                rep(0,10), # Starting=DS, after 1 round of treatment=RRXR
                                rep(0,10), # Starting=DS, after 1 round of treatment=BRXR
                                rep(0,10), # Starting=DS, after 1 round of treatment=RRBRXR
                                0, 0, 0, 0, 0, 0, 0, 0, 0, pan_new_sankey[10,"DS"]), # Starting=DS, after 1 round of treatment=death
                                # Starting RR
                                RR = c(pan_new_sankey[1,"RR"], 0, 0, 0, 0, 0, 0, 0, 0, 0, # Starting=RR, after 1 round of treatment=treatment cure
                                rep(0,10), # Starting=RR, after 1 round of treatment=DS
                                pan_new_sankey[3,"RR"]/sum(pan_new_sankey[3,])*pan_ret_sankey[,"RR"], # Starting=RR, after 1 round of treatment=RR
                                rep(0,10), # Starting=RR, after 1 round of treatment=BR
                                rep(0,10), # Starting=RR, after 1 round of treatment=XR
                                pan_new_sankey[6,"RR"]/sum(pan_new_sankey[6,])*pan_ret_sankey[,"RRBR"], # Starting=RR, after 1 round of treatment=RRBR
                                pan_new_sankey[7,"RR"]/sum(pan_new_sankey[7,])*pan_ret_sankey[,"RRXR"], # Starting=RR, after 1 round of treatment=RRXR
                                rep(0,10), # Starting=RR, after 1 round of treatment=BRXR
                                rep(0,10), # Starting=RR, after 1 round of treatment=RRBRXR
                                0, 0, 0, 0, 0, 0, 0, 0, 0, pan_new_sankey[10,"RR"]), # Starting=RR, after 1 round of treatment=Death
                                #Starting BR
                                BR = c(pan_new_sankey[1,"BR"], 0, 0, 0, 0, 0, 0, 0, 0, 0, # Starting=BR, after 1 round of treatment=treatment cure
                                rep(0,10), # Starting=BR, after 1 round of treatment=DS
                                rep(0,10), # Starting=BR, after 1 round of treatment=RR
                                pan_new_sankey[4,"BR"]/sum(pan_new_sankey[4,])*pan_ret_sankey[,"BR"], # Starting=BR, after 1 round of treatment=BR
                                rep(0,10), # Starting=BR, after 1 round of treatment=RRXR
                                pan_new_sankey[6,"BR"]/sum(pan_new_sankey[6,])*pan_ret_sankey[,"RRBR"], # Starting=BR, after 1 round of treatment=RRBR
                                rep(0,10), # Starting=BR, after 1 round of treatment=RRXR
                                pan_new_sankey[8,"BR"]/sum(pan_new_sankey[8,])*pan_ret_sankey[,"BRXR"], # Starting=BR, after 1 round of treatment=BRXR
                                rep(0,10), # Starting=BR, after 1 round of treatment=RRBRXR
                                0, 0, 0, 0, 0, 0, 0, 0, 0, pan_new_sankey[10,"BR"]), # Starting=BR, after 1 round of treatment=Death)
                                #Starting XR
                                XR = c(pan_new_sankey[1,"XR"], 0, 0, 0, 0, 0, 0, 0, 0, 0, # Starting=BR, after 1 round of treatment=treatment cure
                                rep(0,10), # Starting=XR, after 1 round of treatment=DS
                                rep(0,10), # Starting=XR, after 1 round of treatment=RR
                                rep(0,10), # Starting=XR, after 1 round of treatment=BR
                                pan_new_sankey[5,"XR"]/sum(pan_new_sankey[5,])*pan_ret_sankey[,"XR"], # Starting=XR, after 1 round of treatment=XR
                                rep(0,10), # Starting=XR, after 1 round of treatment=RRBR
                                pan_new_sankey[7,"XR"]/sum(pan_new_sankey[7,])*pan_ret_sankey[,"RRXR"], # Starting=XR, after 1 round of treatment=RRXR
                                pan_new_sankey[8,"XR"]/sum(pan_new_sankey[8,])*pan_ret_sankey[,"BRXR"], # Starting=XR, after 1 round of treatment=BRXR
                                rep(0,10), # Starting=XR, after 1 round of treatment=RRBRXR
                                0, 0, 0, 0, 0, 0, 0, 0, 0, pan_new_sankey[10,"XR"]), # Starting=XR, after 1 round of treatment=Death)   
                                # Starting RRBR
                                RRBR = c(pan_new_sankey[1,"RRBR"], 0, 0, 0, 0, 0, 0, 0, 0, 0, # Starting=RRBR, after 1 round of treatment=treatment cure
                                rep(0,10), # Starting=RRBR, after 1 round of treatment=DS
                                rep(0,10), # Starting=RRBR, after 1 round of treatment=RR
                                rep(0,10), # Starting=RRBR, after 1 round of treatment=BR
                                rep(0,10), # Starting=RRBR, after 1 round of treatment=XR
                                pan_new_sankey[6,"RRBR"]/sum(pan_new_sankey[6,])*pan_ret_sankey[,"RRBR"], # Starting=RRBR, after 1 round of treatment=RRBR
                                rep(0,10), # Starting=RRBR, after 1 round of treatment=RRXR
                                rep(0,10), # Starting=RRBR, after 1 round of treatment=BRXR
                                pan_new_sankey[9,"RRBR"]/sum(pan_new_sankey[9,])*pan_ret_sankey[,"RRBRXR"], # Starting=RRBR, after 1 round of treatment=RRBRXR
                                0, 0, 0, 0, 0, 0, 0, 0, 0, pan_new_sankey[10,"BR"]), # Starting=RRBR, after 1 round of treatment=Death
                                # Starting RRXR
                                RRXR = c(pan_new_sankey[1,"RRXR"], 0, 0, 0, 0, 0, 0, 0, 0, 0, # Starting=XRBR, after 1 round of treatment=treatment cure
                                rep(0,10), # Starting=RRXR, after 1 round of treatment=DS
                                rep(0,10), # Starting=RRXR, after 1 round of treatment=RR
                                rep(0,10), # Starting=RRXR, after 1 round of treatment=BR
                                rep(0,10), # Starting=RRXR, after 1 round of treatment=XR
                                rep(0,10), # Starting=RRXR, after 1 round of treatment=RRBR
                                pan_new_sankey[7,"RRXR"]/sum(pan_new_sankey[7,])*pan_ret_sankey[,"RRXR"], # Starting=RRXR, after 1 round of treatment=RRXR
                                rep(0,10), # Starting=RRXR, after 1 round of treatment=BRXR
                                pan_new_sankey[9,"RRXR"]/sum(pan_new_sankey[9,])*pan_ret_sankey[,"RRBRXR"], # Starting=RRXR, after 1 round of treatment=RRBRXR
                                0, 0, 0, 0, 0, 0, 0, 0, 0, pan_new_sankey[10,"XR"]), # Starting=RRXR, after 1 round of treatment=Death
                                # Starting BRXR
                                BRXR = c(pan_new_sankey[1,"BRXR"], 0, 0, 0, 0, 0, 0, 0, 0, 0, # Starting=BRXR, after 1 round of treatment=treatment cure
                                rep(0,10), # Starting=BRXR, after 1 round of treatment=DS
                                rep(0,10), # Starting=BRXR, after 1 round of treatment=RR
                                rep(0,10), # Starting=BRXR, after 1 round of treatment=BR
                                rep(0,10), # Starting=BRXR, after 1 round of treatment=XR
                                rep(0,10), # Starting=BRXR, after 1 round of treatment=RRBR
                                rep(0,10), # Starting=BRXR, after 1 round of treatment=RRXR
                                pan_new_sankey[8,"BRXR"]/sum(pan_new_sankey[8,])*pan_ret_sankey[,"RRBRXR"], # Starting=BRXR, after 1 round of treatment=BRXR
                                pan_new_sankey[9,"BRXR"]/sum(pan_new_sankey[9,])*pan_ret_sankey[,"RRBRXR"], # Starting=BRXR, after 1 round of treatment=RRBRXR
                                0, 0, 0, 0, 0, 0, 0, 0, 0, pan_new_sankey[10,"BRXR"]), # Starting=BRXR, after 1 round of treatment=Death
                                # Starting RRBRXR
                                RRBRXR = c(pan_new_sankey[1,"RRBRXR"], 0, 0, 0, 0, 0, 0, 0, 0, 0, # Starting=RRBRXR, after 1 round of treatment=treatment cure
                                rep(0,10), # Starting=RRBRXR, after 1 round of treatment=DS
                                rep(0,10), # Starting=RRBRXR, after 1 round of treatment=RR
                                rep(0,10), # Starting=RRBRXR, after 1 round of treatment=BR
                                rep(0,10), # Starting=RRBRXR, after 1 round of treatment=XR
                                rep(0,10), # Starting=RRBRXR, after 1 round of treatment=RRBR
                                rep(0,10), # Starting=RRBRXR, after 1 round of treatment=RRXR
                                rep(0,10), # Starting=RRBRXR, after 1 round of treatment=BRXR
                                pan_new_sankey[9,"RRBRXR"]/sum(pan_new_sankey[9,])*pan_ret_sankey[,"RRBRXR"], # Starting=RRBRXR, after 1 round of treatment=RRBRXR
                                0, 0, 0, 0, 0, 0, 0, 0, 0, pan_new_sankey[10,"RRBRXR"])) # Starting=RRBRXR, after 1 round of treatment=Death                                
       pan_tot_sankey[is.na(pan_tot_sankey)] <- 0
       tmp_out <-prev_out_new+prev_out
       tmp_out$cureN<-sum(total_out_new[,"Cure"])
       tmp_out$cureF<-sum(total_out[,"Cure"])
       tmp_out$death<-sum(total_out[,"Death"])
       
       return(pan_tot_sankey)
       #HEATMAP OUTPUTS: c(new=sum(total_out_new[,"Cure"]), final=sum(total_out[,"Cure"]), death=sum(total_out[,"Death"]), resist=prev_out[5]+prev_out[6]+prev_out[8])
       #SANKEY: pan_tot_sankey
       #PREVALENCE OVER TIME: prev_out_new+prev_out
       #SENSITIVITY OVER TIME: tmp_out
 }

## ================ Heatmap plots ==============================================
 z        <- 1
 rr_split <- seq(0, 0.3, 0.03)
 br_split <- seq(0, 0.3, 0.03)
 samples <- 1000
 #pan_succ <- data.frame(X=rep(rr_split,each=length(br_split)), Y=rep(br_split,length(rr_split)), succ=rep(0,length(rr_split)*length(br_split)))
 #soc_succ <- data.frame(X=rep(rr_split,each=length(br_split)), Y=rep(br_split,length(rr_split)), succ=rep(0,length(rr_split)*length(br_split)))
 tot_succ <- data.frame(X=100*rep(rr_split,each=length(br_split)), Y=100*rep(br_split,length(rr_split)), new=rep(0,length(rr_split)*length(br_split)), final=rep(0,length(rr_split)*length(br_split)), death=rep(0,length(rr_split)*length(br_split)), resist=rep(0,length(rr_split)*length(br_split)))
 comp_new <- rep(0,samples)
 comp_final <- rep(0,samples)
 comp_death <- rep(0,samples)
 comp_resist <- rep(0,samples)
 for (x in rr_split){
   for (y in br_split){
     prev_temp <- data.frame(prev_DS = 1-x-y+x*y, prev_RR = x-x*y, prev_BR = y-x*y, prev_XR = 0, prev_RRBR = x*y, prev_RRXR = 0, prev_BRXR = 0, prev_RRBRXR = 0)
     for (samp in 1:samples){
       params_temp <- params_sample(param_table)
       pan_out <- cohort_pan(prev_temp, params_temp)
       soc_out <- cohort_soc(prev_temp, params_temp)
       comp_new[samp] <-  pan_out["new"] - soc_out["new"]
       comp_final[samp] <-  pan_out["final"] - soc_out["final"]
       comp_death[samp] <-  pan_out["death"] - soc_out["death"]
       comp_resist[samp] <-  pan_out["resist.rrbr_out"] - soc_out["resist.rrbr_out"]
     }
     tot_succ$new[z] <- sum(comp_new > 0)/samples*100
     tot_succ$final[z] <- sum(comp_final > 0)/samples*100
     tot_succ$death[z] <- sum(comp_death > 0)/samples*100
     tot_succ$resist[z] <- sum(comp_resist > 0)/samples*100
     z           <- z+1
 }
 }
 
 # Plotting treatment cure, i.e. where <50% value indicates SoC is better than pan-TB (more likely less cure due to pan-TB than SoC)
 ggplot(tot_succ, aes(X, Y, fill = new)) +
   geom_tile() +
   scale_fill_gradientn(limits= c(0,100), colors = c("#075AFF", "#FFFFCC", "#FF0000")) +
   ylab("BR-TB as a % of all TB") +
   xlab("RR-TB as a % of all TB") +
   guides(fill=guide_legend(title="Prob panTB>SoC (%)"))
 ggsave("New_cure_PB.png",device="png",width=13,height=10,units=c("cm"))
 ggplot(tot_succ, aes(X, Y, fill = final)) +
   geom_tile() +
   scale_fill_gradientn(limits= c(0,100), colors = c("#075AFF", "#FFFFCC", "#FF0000")) +
   ylab("BR-TB as a % of all TB") +
   xlab("RR-TB as a % of all TB") +
   guides(fill=guide_legend(title="Prob panTB>SoC (%)"))
 ggsave("Final_cure_PB.png",device="png",width=13,height=10,units=c("cm"))
 # Plotting deaths or resistance, i.e. where >50% value indicates SoC is better than pan-TB (more likely more deaths/resistance due to pan-TB than SoC)
 ggplot(tot_succ, aes(X, Y, fill = death)) +
   geom_tile() +
   scale_fill_gradientn(limits= c(0,100),colors = c("#FF0000", "#FFFFCC", "#075AFF")) +
   ylab("BR-TB as a % of all TB") +
   xlab("RR-TB as a % of all TB") +
   guides(fill=guide_legend(title="Prob panTB>SoC (%)"))
 ggsave("Final_death_PB.png",device="png",width=13,height=10,units=c("cm"))
 ggplot(tot_succ, aes(X, Y, fill = resist)) +
   geom_tile() +
   scale_fill_gradientn(limits= c(0,100), colors = c("#FF0000", "#FFFFCC", "#075AFF")) +
   ylab("BR-TB as a % of all TB") +
   xlab("RR-TB as a % of all TB") +
   guides(fill=guide_legend(title="Prob panTB>SoC (%)"))
 ggsave("Final_resist_PB.png",device="png",width=13,height=10,units=c("cm"))
 
## ================ Univariate sensitivity analysis ============================
sens <- data.frame(parameter = names(param_ext), final_diff_lo = 0, final_diff_hi = 0, death_diff_lo = 0, death_diff_hi = 0, resist_diff_lo = 0, resist_diff_hi = 0)
final_temp_hi <- c()
final_temp_lo <- c()
death_temp_hi <- c()
death_temp_lo <- c()
resist_temp_hi <- c()
resist_temp_lo <- c()
final_mean <- c()
death_mean <- c()
resist_mean <- c()

params <- data.frame(params_sample(param_table))
runs <- 1000
for (y in 1:runs){
  params[y,] <- params_sample(param_table)
}

for (x in 1:20){
  for (y in 1:runs){
    params_temp_hi <- params[y,]
    params_temp_lo <- params[y,]
    params_temp_hi[,x] <- param_ext[3,x]
    params_temp_lo[,x] <- param_ext[1,x]
    cp_hi <- cohort_pan(prev, params_temp_hi)
    cp_lo <- cohort_pan(prev, params_temp_lo)
    cs_hi <- cohort_soc(prev, params_temp_hi)
    cs_lo <- cohort_soc(prev, params_temp_lo)
    cp_mean <- cohort_pan(prev, params[y,])
    cs_mean <- cohort_soc(prev, params[y,])
    out_hi <- cp_hi - cs_hi
    out_lo <- cp_lo - cs_lo
    out_mean <- cp_mean - cs_mean
    final_temp_hi[y]    <- out_hi["final"]
    final_temp_lo[y]    <- out_lo["final"]
    death_temp_hi[y]    <- out_hi["death"]
    death_temp_lo[y]    <- out_lo["death"]
    resist_temp_hi[y]    <- out_hi["resist.rrbr_out"]
    resist_temp_lo[y]    <- out_lo["resist.rrbr_out"]
    final_mean[y] <- out_mean["final"]
    death_mean[y] <- out_mean["death"]
    resist_mean[y] <- out_mean["resist.rrbr_out"]
  }
  sens[x,2] <- sum(final_temp_hi>0)/runs*100
  sens[x,3] <- sum(final_temp_lo>0)/runs*100
  sens[x,4] <- sum(death_temp_hi>0)/runs*100
  sens[x,5] <- sum(death_temp_lo>0)/runs*100
  sens[x,6] <- sum(resist_temp_hi>0)/runs*100
  sens[x,7] <- sum(resist_temp_lo>0)/runs*100
  
}
# parameter means for comparison
final_comp <- sum(final_mean>0)/runs*100
death_comp <- sum(death_mean>0)/runs*100
resist_comp <- sum(resist_mean>0)/runs*100
# Plot results
 gg1 <- ggplot(sens, aes(x=final_diff_lo, xend=final_diff_hi, y=parameter, group=parameter)) + 
   geom_dumbbell(size_x=1.75, size_xend=1.75, colour_x="#a3c4dc", colour_xend="#e37371") + 
   geom_vline(xintercept = final_comp, linetype="dashed") +
   labs(x="Probability Pan-TB cure>SoC", y="parameter") +
   theme(plot.title = element_text(hjust=0.5, face="bold"),
         plot.background=element_rect(fill="#f7f7f7"),
         panel.background=element_rect(fill="#f7f7f7"),
         panel.grid.minor=element_blank(),
         panel.grid.major.y=element_blank(),
         panel.grid.major.x=element_line(),
         axis.ticks=element_blank(),
         legend.position="top",
         panel.border=element_blank())
 plot(gg1)
 ggsave("Uni_cure.png",device="png",width=10,height=10,units=c("cm"))
 gg2 <- ggplot(sens, aes(x=death_diff_lo, xend=death_diff_hi, y=parameter, group=parameter)) + 
   geom_dumbbell(size_x=1.75, size_xend=1.75, colour_x="#a3c4dc", colour_xend="#e37371") + 
   geom_vline(xintercept = death_comp, linetype="dashed") +
   labs(x="Probability Pan-TB deaths>SoC", y="parameter") +
   theme(plot.title = element_text(hjust=0.5, face="bold"),
         plot.background=element_rect(fill="#f7f7f7"),
         panel.background=element_rect(fill="#f7f7f7"),
         panel.grid.minor=element_blank(),
         panel.grid.major.y=element_blank(),
         panel.grid.major.x=element_line(),
         axis.ticks=element_blank(),
         legend.position="top",
         panel.border=element_blank())
 plot(gg2)
 ggsave("Uni_death.png",device="png",width=10,height=10,units=c("cm"))
 gg3 <- ggplot(sens, aes(x=resist_diff_lo, xend=resist_diff_hi, y=parameter, group=parameter)) + 
   geom_dumbbell(size_x=1.75, size_xend=1.75, colour_x="#a3c4dc", colour_xend="#e37371") + 
   geom_vline(xintercept = resist_comp, linetype="dashed") +
   labs(x="Probability Pan-TB resistance>SoC", y="parameter") +
   theme(plot.title = element_text(hjust=0.5, face="bold"),
         plot.background=element_rect(fill="#f7f7f7"),
         panel.background=element_rect(fill="#f7f7f7"),
         panel.grid.minor=element_blank(),
         panel.grid.major.y=element_blank(),
         panel.grid.major.x=element_line(),
         axis.ticks=element_blank(),
         legend.position="top",
         panel.border=element_blank())
 plot(gg3)
 ggsave("Uni_resist.png",device="png",width=10,height=10,units=c("cm"))
 
 gg1+gg2+gg3
## ================ Prevalence over generations ===============================
time <- 10
runs <- 10
prev_soc <- prev
prev_soc <- prev_soc[rep(seq_len(1), each = time*runs), ]
prev_soc$cohort <- rep(1:10, runs)
prev_soc$run <- rep(1:runs, each=time)
prev_pan <- prev
prev_pan <- prev_pan[rep(seq_len(1), each = time*runs), ]
prev_pan$cohort <- rep(1:10, runs)
prev_pan$run <- rep(1:runs, each=time)
params<-param_ext[2,]  
for (y in 1:runs){
  params[y,] <- params_sample(param_table)
 for (t in 2:time){
 out_soc_t <- cohort_soc(prev_soc[(y-1)*time+t-1,1:8], params[y,]) # Outcomes for previous cohort
 out_soc <- sum(out_soc_t)/(1+sum(out_soc_t)) # proportion of next cohort that is dependent on previous is given by "# who fail"/(1+"# who fail"), assumes as much transmission after failure as before and that incidence ~constant over time
 prev_soc[(y-1)*time+t,1:8] <- (1-out_soc)*prev_soc[(y-1)*time+t-1,1:8]+out_soc*out_soc_t/sum(out_soc_t) # Prevalence is a combination of previous prevalence + outcomes from previous cohort
 out_pan_t <- cohort_pan(prev_pan[(y-1)*time+t-1,1:8], params[y,]) # Outcomes for previous cohort
 out_pan <- sum(out_pan_t)/(1+sum(out_pan_t))
 prev_pan[(y-1)*time+t,1:8] <- (1-out_pan)*prev_pan[(y-1)*time+t-1,1:8]+out_pan*out_pan_t/sum(out_pan_t)  # Prevalence is a combination of previous prevalence + outcomes from previous cohort
 }
}

out_soc <- sum(out_soc_t[1:8])/(1+sum(out_soc_t[1:8]))
prev_soc[(y-1)*time+t,1:8] <- (1-out_soc)*prev_soc[(y-1)*time+t-1,1:8]+out_soc*out_soc_t[1:8]/sum(out_soc_t[1:8]) # Prevalence is a combination of previous prevalence + outcomes from previous cohort
out_pan_t <- cohort_pan(prev_pan[(y-1)*time+t-1,1:8], params[y,]) # Outcomes for previous cohort
out_pan <- sum(out_pan_t[1:8])/(1+sum(out_pan_t[1:8]))
prev_pan[(y-1)*time+t,1:8] <- (1-out_pan)*prev_pan[(y-1)*time+t-1,1:8]+out_pan*out_pan_t[1:8]/sum(out_pan_t[1:8])


soc_mean <- prev
soc_sd <- prev
pan_mean <- prev
pan_sd <- prev
for (t in 1:time){
  soc_mean[t,] <-colMeans(prev_soc[prev_soc$cohort==t,1:8])
  soc_sd[t,] <- sapply(prev_soc[prev_soc$cohort==t,1:8], sd)
  pan_mean[t,] <-colMeans(prev_pan[prev_pan$cohort==t,1:8])
  pan_sd[t,] <- sapply(prev_pan[prev_pan$cohort==t,1:8], sd)
}
# Plot results
data_soc <- data.frame(cohort=rep(1:time,8), pheno=rep(c("DS", "RR", "BR", "XR", "RRBR", "RRXR", "BRXR", "RRBRXR"), each=10), prev_mean=c(soc_mean$prev_DS, soc_mean$prev_RR, soc_mean$prev_BR, soc_mean$prev_XR, soc_mean$prev_RRBR, soc_mean$prev_RRXR, soc_mean$prev_BRXR, soc_mean$prev_RRBRXR),  prev_sd=c(soc_sd$prev_DS, soc_sd$prev_RR, soc_sd$prev_BR, soc_sd$prev_XR, soc_sd$prev_RRBR, soc_sd$prev_RRXR, soc_sd$prev_BRXR, soc_sd$prev_RRBRXR))
ggplot(data=data_soc, aes(x=cohort))+
  geom_line(aes(y=prev_mean, color=pheno), size=1)+
geom_ribbon(aes(y = prev_mean, ymin = prev_mean - 1.96*prev_sd, ymax = prev_mean +  1.96*prev_sd, fill = pheno), alpha = .2)+
  # To look at DR-TB only
  coord_cartesian(ylim = c(0, 1)) #scale_y_continuous(limits = c(0, 1))
ggsave("Time_soc_DS_betachi.png",device="png",width=10,height=10,units=c("cm"))

data_pan <- data.frame(cohort=rep(1:time,8), pheno=rep(c("DS", "RR", "BR", "XR", "RRBR", "RRXR", "BRXR", "RRBRXR"), each=10), prev_mean=c(pan_mean$prev_DS, pan_mean$prev_RR, pan_mean$prev_BR, pan_mean$prev_XR, pan_mean$prev_RRBR, pan_mean$prev_RRXR, pan_mean$prev_BRXR, pan_mean$prev_RRBRXR),  prev_sd=c(pan_sd$prev_DS, pan_sd$prev_RR, pan_sd$prev_BR, pan_sd$prev_XR, pan_sd$prev_RRBR, pan_sd$prev_RRXR, pan_sd$prev_BRXR, pan_sd$prev_RRBRXR))
ggplot(data=data_pan, aes(x=cohort))+
  geom_line(aes(y=prev_mean, color=pheno), size=1)+
  geom_ribbon(aes(y = prev_mean, ymin = prev_mean - 1.96*prev_sd, ymax = prev_mean +  1.96*prev_sd, fill = pheno), alpha = .2)+
  # To look at DR-TB only
  coord_cartesian(ylim = c(0, 1)) #scale_y_continuous(limits = c(0, 1))
ggsave("Time_pan_DS_betachi.png",device="png",width=10,height=10,units=c("cm"))

# Problematic resistance
problem <- data.frame(cohort=c(seq(1,10,by=1),seq(1,10,by=1)))
problem$scenario <- c(rep("SoC",10),rep("Pan",10))
problem$prev_mean <- c(data_soc[data_soc$pheno=="RRBR",3]+data_soc[data_soc$pheno=="RRXR",3]+data_soc[data_soc$pheno=="RRBRXR",3],data_pan[data_pan$pheno=="RRBR",3]+data_pan[data_pan$pheno=="RRXR",3]+data_pan[data_pan$pheno=="RRBRXR",3])
problem$prev_sd <- c(sqrt(data_soc[data_soc$pheno=="RRBR",4]^2+data_soc[data_soc$pheno=="RRXR",4]^2+data_soc[data_soc$pheno=="RRBRXR",4]^2),sqrt(data_pan[data_pan$pheno=="RRBR",4]^2+data_pan[data_pan$pheno=="RRXR",4]^2+data_pan[data_pan$pheno=="RRBRXR",4]^2))

ggplot(data=problem, aes(x=cohort))+
  geom_line(aes(y=prev_mean, color=scenario), size=1)+
  geom_ribbon(aes(y = prev_mean, ymin = prev_mean - 1.96*prev_sd, ymax = prev_mean +  1.96*prev_sd, fill = scenario), alpha = .2)+
  # To look at DR-TB only
  coord_cartesian(ylim = c(0, 0.1)) #scale_y_continuous(limits = c(0, 1))
ggsave("Time_problem_betachi.png",device="png",width=10,height=10,units=c("cm"))

## ================ Heatmap of durable cure ====================================
z        <- 1
eff_split <- seq(0, 0.178, 0.0178)
br_split <- seq(0, 0.3, 0.03)
samples <- 100
#pan_succ <- data.frame(X=rep(rr_split,each=length(br_split)), Y=rep(br_split,length(rr_split)), succ=rep(0,length(rr_split)*length(br_split)))
#soc_succ <- data.frame(X=rep(rr_split,each=length(br_split)), Y=rep(br_split,length(rr_split)), succ=rep(0,length(rr_split)*length(br_split)))
tot_succ <- data.frame(X=100*rep(eff_split,each=length(br_split)), Y=100*rep(br_split,length(eff_split)), new=rep(0,length(eff_split)*length(br_split)), final=rep(0,length(eff_split)*length(br_split)), death=rep(0,length(eff_split)*length(br_split)), resist=rep(0,length(eff_split)*length(br_split)))
comp_new <- rep(0,samples)
comp_final <- rep(0,samples)
comp_death <- rep(0,samples)
comp_resist <- rep(0,samples)
for (x in eff_split){
  for (y in br_split){
    prev_temp <- data.frame(prev_DS = 1-y+0.042*y, prev_RR = 0.042, prev_BR = y-0.042*y, prev_XR = 0, prev_RRBR = 0.042*y, prev_RRXR = 0, prev_BRXR = 0, prev_RRBRXR = 0)
    for (samp in 1:samples){
      params_temp <- params_sample(param_table)
      params_temp$E_BX <- 0.763
      params_temp$E_R <- 0.763-x
      pan_out <- cohort_pan(prev_temp, params_temp)
      soc_out <- cohort_soc(prev_temp, params_temp)
      comp_new[samp] <-  pan_out["new"] - soc_out["new"]
      comp_final[samp] <-  pan_out["final"] - soc_out["final"]
      comp_death[samp] <-  pan_out["death"] - soc_out["death"]
      comp_resist[samp] <-  pan_out["resist.rrbr_out"] - soc_out["resist.rrbr_out"]
    }
    tot_succ$new[z] <- sum(comp_new > 0)/samples*100
    tot_succ$final[z] <- sum(comp_final > 0)/samples*100
    tot_succ$death[z] <- sum(comp_death > 0)/samples*100
    tot_succ$resist[z] <- sum(comp_resist > 0)/samples*100
    z           <- z+1
  }
}

# Plotting treatment cure, i.e. where <50% value indicates SoC is better than pan-TB (more likely less cure due to pan-TB than SoC)
ggplot(tot_succ, aes(X, Y, fill = new)) +
  geom_tile() +
  scale_fill_gradientn(limits= c(0,100), colors = c("#075AFF", "#FFFFCC", "#FF0000")) +
  ylab("BR-TB as a % of all TB") +
  xlab("Pan-TB % improvement in durable cure (E_BX - E_R)") +
  guides(fill=guide_legend(title="Prob panTB>SoC (%)"))
ggsave("Eff_New_cure.png",device="png",width=13,height=10,units=c("cm"))
ggplot(tot_succ, aes(X, Y, fill = final)) +
  geom_tile() +
  scale_fill_gradientn(limits= c(0,100), colors = c("#075AFF", "#FFFFCC", "#FF0000")) +
  ylab("BR-TB as a % of all TB") +
  xlab("Pan-TB % improvement in durable cure (E_BX - E_R)") +
  guides(fill=guide_legend(title="Prob panTB>SoC (%)"))
ggsave("Eff_Final_cure.png",device="png",width=13,height=10,units=c("cm"))
# Plotting deaths or resistance, i.e. where >50% value indicates SoC is better than pan-TB (more likely more deaths/resistance due to pan-TB than SoC)
ggplot(tot_succ, aes(X, Y, fill = death)) +
  geom_tile() +
  scale_fill_gradientn(limits= c(0,100),colors = c("#FF0000", "#FFFFCC", "#075AFF")) +
  ylab("BR-TB as a % of all TB") +
  xlab("Pan-TB % improvement in durable cure (E_BX - E_R)") +
  guides(fill=guide_legend(title="Prob panTB>SoC (%)"))
ggsave("Eff_Final_death.png",device="png",width=13,height=10,units=c("cm"))
ggplot(tot_succ, aes(X, Y, fill = resist)) +
  geom_tile() +
  scale_fill_gradientn(limits= c(0,100), colors = c("#FF0000", "#FFFFCC", "#075AFF")) +
  ylab("BR-TB as a % of all TB") +
  xlab("Pan-TB % improvement in durable cure (E_BX - E_R)") +
  guides(fill=guide_legend(title="Prob panTB>SoC (%)"))
ggsave("Eff_Final_resist.png",device="png",width=13,height=10,units=c("cm"))

## ================ Sensitivity updating prevalence  ===============================
#prev_DS    prev_RR    prev_BR    prev_XR  prev_RRBR  prev_RRXR prev_BRXR prev_RRBRXR
#0.6349823 0.02960824 0.07974126 0.08032826 0.01349257 0.00432169 0.1221082  0.03541745
prev1 <- data.frame(prev_DS = pan_mean$prev_DS[10],	
                    prev_RR         = pan_mean$prev_RR[10],	
                    prev_BR         = pan_mean$prev_BR[10],	
                    prev_XR         = pan_mean$prev_XR[10],	
                    prev_RRBR       = pan_mean$prev_RRBR[10],	
                    prev_RRXR       = pan_mean$prev_RRXR[10],	
                    prev_BRXR       = pan_mean$prev_BRXR[10],	
                    prev_RRBRXR     = pan_mean$prev_RRBRXR[10])
sens1 <- data.frame(parameter = names(param_ext), resist_diff_lo = 0, resist_diff_hi = 0, resist_diff_mean = 0, cureN_diff_lo = 0, cureN_diff_hi = 0, cureN_diff_mean = 0, cureF_diff_lo = 0, cureF_diff_hi = 0, cureF_diff_mean = 0, death_diff_lo = 0, death_diff_hi = 0, death_diff_mean = 0)
sens2 <- data.frame(parameter = names(param_ext), resist_diff_lo = 0, resist_diff_hi = 0, resist_diff_mean = 0, cureN_diff_lo = 0, cureN_diff_hi = 0, cureN_diff_mean = 0, cureF_diff_lo = 0, cureF_diff_hi = 0, cureF_diff_mean = 0, death_diff_lo = 0, death_diff_hi = 0, death_diff_mean = 0)
sens3 <- data.frame(parameter = names(param_ext), resist_diff_lo = 0, resist_diff_hi = 0, resist_diff_mean = 0, cureN_diff_lo = 0, cureN_diff_hi = 0, cureN_diff_mean = 0, cureF_diff_lo = 0, cureF_diff_hi = 0, cureF_diff_mean = 0, death_diff_lo = 0, death_diff_hi = 0, death_diff_mean = 0)
sens4 <- data.frame(parameter = names(param_ext), resist_diff_lo = 0, resist_diff_hi = 0, resist_diff_mean = 0, cureN_diff_lo = 0, cureN_diff_hi = 0, cureN_diff_mean = 0, cureF_diff_lo = 0, cureF_diff_hi = 0, cureF_diff_mean = 0, death_diff_lo = 0, death_diff_hi = 0, death_diff_mean = 0)
sens5 <- data.frame(parameter = names(param_ext), resist_diff_lo = 0, resist_diff_hi = 0, resist_diff_mean = 0, cureN_diff_lo = 0, cureN_diff_hi = 0, cureN_diff_mean = 0, cureF_diff_lo = 0, cureF_diff_hi = 0, cureF_diff_mean = 0, death_diff_lo = 0, death_diff_hi = 0, death_diff_mean = 0)
sens6 <- data.frame(parameter = names(param_ext), resist_diff_lo = 0, resist_diff_hi = 0, resist_diff_mean = 0, cureN_diff_lo = 0, cureN_diff_hi = 0, cureN_diff_mean = 0, cureF_diff_lo = 0, cureF_diff_hi = 0, cureF_diff_mean = 0, death_diff_lo = 0, death_diff_hi = 0, death_diff_mean = 0)
sens7 <- data.frame(parameter = names(param_ext), resist_diff_lo = 0, resist_diff_hi = 0, resist_diff_mean = 0, cureN_diff_lo = 0, cureN_diff_hi = 0, cureN_diff_mean = 0, cureF_diff_lo = 0, cureF_diff_hi = 0, cureF_diff_mean = 0, death_diff_lo = 0, death_diff_hi = 0, death_diff_mean = 0)
sens8 <- data.frame(parameter = names(param_ext), resist_diff_lo = 0, resist_diff_hi = 0, resist_diff_mean = 0, cureN_diff_lo = 0, cureN_diff_hi = 0, cureN_diff_mean = 0, cureF_diff_lo = 0, cureF_diff_hi = 0, cureF_diff_mean = 0, death_diff_lo = 0, death_diff_hi = 0, death_diff_mean = 0)
sens9 <- data.frame(parameter = names(param_ext), resist_diff_lo = 0, resist_diff_hi = 0, resist_diff_mean = 0, cureN_diff_lo = 0, cureN_diff_hi = 0, cureN_diff_mean = 0, cureF_diff_lo = 0, cureF_diff_hi = 0, cureF_diff_mean = 0, death_diff_lo = 0, death_diff_hi = 0, death_diff_mean = 0)
sens10 <- data.frame(parameter = names(param_ext), resist_diff_lo = 0, resist_diff_hi = 0, resist_diff_mean = 0, cureN_diff_lo = 0, cureN_diff_hi = 0, cureN_diff_mean = 0, cureF_diff_lo = 0, cureF_diff_hi = 0, cureF_diff_mean = 0, death_diff_lo = 0, death_diff_hi = 0, death_diff_mean = 0)
resist_temp_hi <- data.frame()
resist_temp_lo <- data.frame()
resist_mean <- data.frame()

cureN_temp_hi <- data.frame()
cureN_temp_lo <- data.frame()
cureN_mean <- data.frame()

cureF_temp_hi <- data.frame()
cureF_temp_lo <- data.frame()
cureF_mean <- data.frame()

death_temp_hi <- data.frame()
death_temp_lo <- data.frame()
death_mean <- data.frame()
time <- 10

params <- data.frame(params_sample(param_table))
runs <- 100
for (y in 1:runs){
  params[y,] <- params_sample(param_table)
}

cp_hi <- c()
cp_lo <- c()
cs_hi <- c()
cs_lo <- c()
cp_mean <- c()
cs_mean <- c()
out_hi <- c()
out_lo <- c()
out_mean <-c()

cureN_hi <-   c()
cureN_lo <-  c()
cureN_t_mean <-   c()
cureF_hi <-   c()
cureF_lo <-   c()
cureF_t_mean <-   c()
death_hi <-   c()
death_lo <-  c()
death_t_mean <-  c()

start_time <- Sys.time()

for (x in 1:20){
  prev_soc_lo <- prev
  prev_soc_lo <- prev_soc_lo[rep(seq_len(1), each = time*runs), ]
  prev_soc_lo$cohort <- rep(1:time, runs)
  prev_soc_lo$run <- rep(1:runs, each=time)
  prev_soc_hi <- prev
  prev_soc_hi <- prev_soc_hi[rep(seq_len(1), each = time*runs), ]
  prev_soc_hi$cohort <- rep(1:time, runs)
  prev_soc_hi$run <- rep(1:runs, each=time)
  prev_soc <- prev
  prev_soc <- prev_soc[rep(seq_len(1), each = time*runs), ]
  prev_soc$cohort <- rep(1:time, runs)
  prev_soc$run <- rep(1:runs, each=time)
  prev_pan_lo <- prev
  prev_pan_lo <- prev_pan_lo[rep(seq_len(1), each = time*runs), ]
  prev_pan_lo$cohort <- rep(1:time, runs)
  prev_pan_lo$run <- rep(1:runs, each=time)
  prev_pan_hi <- prev
  prev_pan_hi <- prev_pan_hi[rep(seq_len(1), each = time*runs), ]
  prev_pan_hi$cohort <- rep(1:time, runs)
  prev_pan_hi$run <- rep(1:runs, each=time)
  prev_pan <- prev
  prev_pan <- prev_pan[rep(seq_len(1), each = time*runs), ]
  prev_pan$cohort <- rep(1:time, runs)
  prev_pan$run <- rep(1:runs, each=time)
  
  for (y in 1:runs){
    params_temp_hi <- params[y,]
    params_temp_lo <- params[y,]
    params_temp_hi[,x] <- param_ext[3,x]
    params_temp_lo[,x] <- param_ext[1,x]
    
    cp_hi[1] <- prev_pan_hi[y*1,5]+prev_pan_hi[y*1,6]+prev_pan_hi[y*1,8]
    cp_lo[1] <- prev_pan_lo[y*1,5]+prev_pan_lo[y*1,6]+prev_pan_lo[y*1,8]
    cs_hi[1] <- prev_soc_hi[y*1,5]+prev_soc_hi[y*1,6]+prev_soc_hi[y*1,8]
    cs_lo[1] <- prev_soc_lo[y*1,5]+prev_soc_lo[y*1,6]+prev_soc_lo[y*1,8]
    cp_mean[1] <- prev_pan[y*1,5]+prev_pan[y*1,6]+prev_pan[y*1,8]
    cs_mean[1] <- prev_soc[y*1,5]+prev_soc[y*1,6]+prev_soc[y*1,8]
    out_hi[1] <- cp_hi[1] - cs_hi[1]
    out_lo[1] <- cp_lo[1] - cs_lo[1]
    out_mean[1] <- cp_mean[1] - cs_mean[1]
    
    out_s_hi <- cohort_soc(prev, params_temp_hi)
    out_s_lo <- cohort_soc(prev, params_temp_lo)
    out_p_hi <- cohort_pan(prev, params_temp_hi)
    out_p_lo <- cohort_pan(prev, params_temp_lo)
    out_s_mean <- cohort_soc(prev, params[y,])
    out_p_mean <- cohort_pan(prev, params[y,])
    
    cureN_hi[1] <-  out_p_hi$cureN - out_s_hi$cureN
    cureN_lo[1] <-  out_p_lo$cureN - out_s_lo$cureN
    cureN_t_mean[1] <-  out_p_mean$cureN - out_s_mean$cureN
    cureF_hi[1] <-  out_p_hi$cureF - out_s_hi$cureF
    cureF_lo[1] <-  out_p_lo$cureF - out_s_lo$cureF
    cureF_t_mean[1] <-  out_p_mean$cureF - out_s_mean$cureF
    death_hi[1] <-  out_p_hi$death - out_s_hi$death
    death_lo[1] <-  out_p_lo$death - out_s_lo$death
    death_t_mean[1] <-  out_p_mean$death - out_s_mean$death
    
    for (t in 2:time){
      #lo values
      out_soc_t_lo <- cohort_soc(prev_soc_lo[(y-1)*time+t-1,1:8], params_temp_lo) # Outcomes for previous cohort
      out_soc_lo <- sum(out_soc_t_lo[1:8])/(1+sum(out_soc_t_lo[1:8]))
      prev_soc_lo[(y-1)*time+t,1:8] <- (1-out_soc_lo)*prev_soc_lo[(y-1)*time+t-1,1:8]+out_soc_lo*out_soc_t_lo[1:8]/sum(out_soc_t_lo[1:8]) # Prevalence is a combination of previous prevalence + outcomes from previous cohort
      out_pan_t_lo <- cohort_pan(prev_pan_lo[(y-1)*time+t-1,1:8], params_temp_lo) # Outcomes for previous cohort
      out_pan_lo <- sum(out_pan_t_lo[1:8])/(1+sum(out_pan_t_lo[1:8]))
      prev_pan_lo[(y-1)*time+t,1:8] <- (1-out_pan_lo)*prev_pan_lo[(y-1)*time+t-1,1:8]+out_pan_lo*out_pan_t_lo[1:8]/sum(out_pan_t_lo[1:8]) # Prevalence is a combination of previous prevalence + outcomes from previous cohort
      #hi values
      out_soc_t_hi <- cohort_soc(prev_soc_hi[(y-1)*time+t-1,1:8], params_temp_hi) # Outcomes for previous cohort
      out_soc_hi <- sum(out_soc_t_hi[1:8])/(1+sum(out_soc_t_hi[1:8]))
      prev_soc_hi[(y-1)*time+t,1:8] <- (1-out_soc_hi)*prev_soc_hi[(y-1)*time+t-1,1:8]+out_soc_hi*out_soc_t_hi[1:8]/sum(out_soc_t_hi[1:8]) # Prevalence is a combination of previous prevalence + outcomes from previous cohort
      out_pan_t_hi <- cohort_pan(prev_pan_hi[(y-1)*time+t-1,1:8], params_temp_hi) # Outcomes for previous cohort
      out_pan_hi <- sum(out_pan_t_hi[1:8])/(1+sum(out_pan_t_hi[1:8]))
      prev_pan_hi[(y-1)*time+t,1:8] <- (1-out_pan_hi)*prev_pan_hi[(y-1)*time+t-1,1:8]+out_pan_hi*out_pan_t_hi[1:8]/sum(out_pan_t_hi[1:8]) # Prevalence is a combination of previous prevalence + outcomes from previous cohort
      #mean values
      out_soc_t <- cohort_soc(prev_soc[(y-1)*time+t-1,1:8], params[y,]) # Outcomes for previous cohort
      out_soc <- sum(out_soc_t[1:8])/(1+sum(out_soc_t[1:8]))
      prev_soc[(y-1)*time+t,1:8] <- (1-out_soc)*prev_soc[(y-1)*time+t-1,1:8]+out_soc*out_soc_t[1:8]/sum(out_soc_t[1:8]) # Prevalence is a combination of previous prevalence + outcomes from previous cohort
      out_pan_t <- cohort_pan(prev_pan[(y-1)*time+t-1,1:8], params[y,]) # Outcomes for previous cohort
      out_pan <- sum(out_pan_t[1:8])/(1+sum(out_pan_t[1:8]))
      prev_pan[(y-1)*time+t,1:8] <- (1-out_pan)*prev_pan[(y-1)*time+t-1,1:8]+out_pan*out_pan_t[1:8]/sum(out_pan_t[1:8]) # Prevalence is a combination of previous prevalence + outcomes from previous cohort
      
    cp_hi[t] <- prev_pan_hi[y*t,5]+prev_pan_hi[y*t,6]+prev_pan_hi[y*t,8]
    cp_lo[t] <- prev_pan_lo[y*t,5]+prev_pan_lo[y*t,6]+prev_pan_lo[y*t,8]
    cs_hi[t] <- prev_soc_hi[y*t,5]+prev_soc_hi[y*t,6]+prev_soc_hi[y*t,8]
    cs_lo[t] <- prev_soc_lo[y*t,5]+prev_soc_lo[y*t,6]+prev_soc_lo[y*t,8]
    cp_mean[t] <- prev_pan[y*t,5]+prev_pan[y*t,6]+prev_pan[y*t,8]
    cs_mean[t] <- prev_soc[y*t,5]+prev_soc[y*t,6]+prev_soc[y*t,8]
    out_hi[t] <- cp_hi[t] - cs_hi[t]
    out_lo[t] <- cp_lo[t] - cs_lo[t]
    out_mean[t] <- cp_mean[t] - cs_mean[t]
    
    cureN_hi[t] <-  out_pan_t_hi$cureN - out_soc_t_hi$cureN
    cureN_lo[t] <-  out_pan_t_lo$cureN - out_soc_t_lo$cureN
    cureN_t_mean[t] <-  out_pan_t$cureN - out_soc_t$cureN
    cureF_hi[t] <-  out_pan_t_hi$cureF - out_soc_t_hi$cureF
    cureF_lo[t] <-  out_pan_t_lo$cureF - out_soc_t_lo$cureF
    cureF_t_mean[t] <-  out_pan_t$cureF - out_soc_t$cureF
    death_hi[t] <-  out_pan_t_hi$death - out_soc_t_hi$death
    death_lo[t] <-  out_pan_t_lo$death - out_soc_t_lo$death
    death_t_mean[t] <-  out_pan_t$death - out_soc_t$death
    }
    
    resist_temp_hi[y,1:time]    <- out_hi
    resist_temp_lo[y,1:time]    <- out_lo
    resist_mean[y,1:time] <- out_mean
    
    cureN_temp_hi[y,1:time]    <- cureN_hi
    cureN_temp_lo[y,1:time]    <- cureN_lo
    cureN_mean[y,1:time] <- cureN_t_mean
    
    cureF_temp_hi[y,1:time]    <- cureF_hi
    cureF_temp_lo[y,1:time]    <- cureF_lo
    cureF_mean[y,1:time] <- cureF_t_mean
    
    death_temp_hi[y,1:time]    <- death_hi
    death_temp_lo[y,1:time]    <- death_lo
    death_mean[y,1:time] <- death_t_mean
  }
  sens1[x,3] <- sum(resist_temp_hi[,1]>0)/runs*100
  sens1[x,2] <- sum(resist_temp_lo[,1]>0)/runs*100
  sens1[x,4] <- sum(resist_mean[,1]>0)/runs*100
  sens1[x,6] <- sum(cureN_temp_hi[,1]>0)/runs*100
  sens1[x,5] <- sum(cureN_temp_lo[,1]>0)/runs*100
  sens1[x,7] <- sum(cureN_mean[,1]>0)/runs*100
  sens1[x,9] <- sum(cureF_temp_hi[,1]>0)/runs*100
  sens1[x,8] <- sum(cureF_temp_lo[,1]>0)/runs*100
  sens1[x,10] <- sum(cureF_mean[,1]>0)/runs*100
  sens1[x,12] <- sum(death_temp_hi[,1]>0)/runs*100
  sens1[x,11] <- sum(death_temp_lo[,1]>0)/runs*100
  sens1[x,13] <- sum(death_mean[,1]>0)/runs*100
  
  sens2[x,3] <- sum(resist_temp_hi[,2]>0)/runs*100
  sens2[x,2] <- sum(resist_temp_lo[,2]>0)/runs*100
  sens2[x,4] <- sum(resist_mean[,2]>0)/runs*100
  sens2[x,6] <- sum(cureN_temp_hi[,2]>0)/runs*100
  sens2[x,5] <- sum(cureN_temp_lo[,2]>0)/runs*100
  sens2[x,7] <- sum(cureN_mean[,2]>0)/runs*100
  sens2[x,9] <- sum(cureF_temp_hi[,2]>0)/runs*100
  sens2[x,8] <- sum(cureF_temp_lo[,2]>0)/runs*100
  sens2[x,10] <- sum(cureF_mean[,2]>0)/runs*100
  sens2[x,12] <- sum(death_temp_hi[,2]>0)/runs*100
  sens2[x,11] <- sum(death_temp_lo[,2]>0)/runs*100
  sens2[x,13] <- sum(death_mean[,2]>0)/runs*100
  
  sens3[x,3] <- sum(resist_temp_hi[,3]>0)/runs*100
  sens3[x,2] <- sum(resist_temp_lo[,3]>0)/runs*100
  sens3[x,4] <- sum(resist_mean[,3]>0)/runs*100
  sens3[x,6] <- sum(cureN_temp_hi[,3]>0)/runs*100
  sens3[x,5] <- sum(cureN_temp_lo[,3]>0)/runs*100
  sens3[x,7] <- sum(cureN_mean[,3]>0)/runs*100
  sens3[x,9] <- sum(cureF_temp_hi[,3]>0)/runs*100
  sens3[x,8] <- sum(cureF_temp_lo[,3]>0)/runs*100
  sens3[x,10] <- sum(cureF_mean[,3]>0)/runs*100
  sens3[x,12] <- sum(death_temp_hi[,3]>0)/runs*100
  sens3[x,11] <- sum(death_temp_lo[,3]>0)/runs*100
  sens3[x,13] <- sum(death_mean[,3]>0)/runs*100
  
  sens4[x,3] <- sum(resist_temp_hi[,4]>0)/runs*100
  sens4[x,2] <- sum(resist_temp_lo[,4]>0)/runs*100
  sens4[x,4] <- sum(resist_mean[,4]>0)/runs*100
  sens4[x,6] <- sum(cureN_temp_hi[,4]>0)/runs*100
  sens4[x,5] <- sum(cureN_temp_lo[,4]>0)/runs*100
  sens4[x,7] <- sum(cureN_mean[,4]>0)/runs*100
  sens4[x,9] <- sum(cureF_temp_hi[,4]>0)/runs*100
  sens4[x,8] <- sum(cureF_temp_lo[,4]>0)/runs*100
  sens4[x,10] <- sum(cureF_mean[,4]>0)/runs*100
  sens4[x,12] <- sum(death_temp_hi[,4]>0)/runs*100
  sens4[x,11] <- sum(death_temp_lo[,4]>0)/runs*100
  sens4[x,13] <- sum(death_mean[,4]>0)/runs*100
  
  sens5[x,3] <- sum(resist_temp_hi[,5]>0)/runs*100
  sens5[x,2] <- sum(resist_temp_lo[,5]>0)/runs*100
  sens5[x,4] <- sum(resist_mean[,5]>0)/runs*100
  sens5[x,6] <- sum(cureN_temp_hi[,5]>0)/runs*100
  sens5[x,5] <- sum(cureN_temp_lo[,5]>0)/runs*100
  sens5[x,7] <- sum(cureN_mean[,5]>0)/runs*100
  sens5[x,9] <- sum(cureF_temp_hi[,5]>0)/runs*100
  sens5[x,8] <- sum(cureF_temp_lo[,5]>0)/runs*100
  sens5[x,10] <- sum(cureF_mean[,5]>0)/runs*100
  sens5[x,12] <- sum(death_temp_hi[,5]>0)/runs*100
  sens5[x,11] <- sum(death_temp_lo[,5]>0)/runs*100
  sens5[x,13] <- sum(death_mean[,5]>0)/runs*100
  
  sens6[x,3] <- sum(resist_temp_hi[,6]>0)/runs*100
  sens6[x,2] <- sum(resist_temp_lo[,6]>0)/runs*100
  sens6[x,4] <- sum(resist_mean[,6]>0)/runs*100
  sens6[x,6] <- sum(cureN_temp_hi[,6]>0)/runs*100
  sens6[x,5] <- sum(cureN_temp_lo[,6]>0)/runs*100
  sens6[x,7] <- sum(cureN_mean[,6]>0)/runs*100
  sens6[x,9] <- sum(cureF_temp_hi[,6]>0)/runs*100
  sens6[x,8] <- sum(cureF_temp_lo[,6]>0)/runs*100
  sens6[x,10] <- sum(cureF_mean[,6]>0)/runs*100
  sens6[x,12] <- sum(death_temp_hi[,6]>0)/runs*100
  sens6[x,11] <- sum(death_temp_lo[,6]>0)/runs*100
  sens6[x,13] <- sum(death_mean[,6]>0)/runs*100
  
  sens7[x,3] <- sum(resist_temp_hi[,7]>0)/runs*100
  sens7[x,2] <- sum(resist_temp_lo[,7]>0)/runs*100
  sens7[x,4] <- sum(resist_mean[,7]>0)/runs*100
  sens7[x,6] <- sum(cureN_temp_hi[,7]>0)/runs*100
  sens7[x,5] <- sum(cureN_temp_lo[,7]>0)/runs*100
  sens7[x,7] <- sum(cureN_mean[,7]>0)/runs*100
  sens7[x,9] <- sum(cureF_temp_hi[,7]>0)/runs*100
  sens7[x,8] <- sum(cureF_temp_lo[,7]>0)/runs*100
  sens7[x,10] <- sum(cureF_mean[,7]>0)/runs*100
  sens7[x,12] <- sum(death_temp_hi[,7]>0)/runs*100
  sens7[x,11] <- sum(death_temp_lo[,7]>0)/runs*100
  sens7[x,13] <- sum(death_mean[,7]>0)/runs*100
  
  sens8[x,3] <- sum(resist_temp_hi[,8]>0)/runs*100
  sens8[x,2] <- sum(resist_temp_lo[,8]>0)/runs*100
  sens8[x,4] <- sum(resist_mean[,8]>0)/runs*100
  sens8[x,6] <- sum(cureN_temp_hi[,8]>0)/runs*100
  sens8[x,5] <- sum(cureN_temp_lo[,8]>0)/runs*100
  sens8[x,7] <- sum(cureN_mean[,8]>0)/runs*100
  sens8[x,9] <- sum(cureF_temp_hi[,8]>0)/runs*100
  sens8[x,8] <- sum(cureF_temp_lo[,8]>0)/runs*100
  sens8[x,10] <- sum(cureF_mean[,8]>0)/runs*100
  sens8[x,12] <- sum(death_temp_hi[,8]>0)/runs*100
  sens8[x,11] <- sum(death_temp_lo[,8]>0)/runs*100
  sens8[x,13] <- sum(death_mean[,8]>0)/runs*100
  
  sens9[x,3] <- sum(resist_temp_hi[,9]>0)/runs*100
  sens9[x,2] <- sum(resist_temp_lo[,9]>0)/runs*100
  sens9[x,4] <- sum(resist_mean[,9]>0)/runs*100
  sens9[x,6] <- sum(cureN_temp_hi[,9]>0)/runs*100
  sens9[x,5] <- sum(cureN_temp_lo[,9]>0)/runs*100
  sens9[x,7] <- sum(cureN_mean[,9]>0)/runs*100
  sens9[x,9] <- sum(cureF_temp_hi[,9]>0)/runs*100
  sens9[x,8] <- sum(cureF_temp_lo[,9]>0)/runs*100
  sens9[x,10] <- sum(cureF_mean[,9]>0)/runs*100
  sens9[x,12] <- sum(death_temp_hi[,9]>0)/runs*100
  sens9[x,11] <- sum(death_temp_lo[,9]>0)/runs*100
  sens9[x,13] <- sum(death_mean[,9]>0)/runs*100
  
  sens10[x,3] <- sum(resist_temp_hi[,10]>0)/runs*100
  sens10[x,2] <- sum(resist_temp_lo[,10]>0)/runs*100
  sens10[x,4] <- sum(resist_mean[,10]>0)/runs*100
  sens10[x,6] <- sum(cureN_temp_hi[,10]>0)/runs*100
  sens10[x,5] <- sum(cureN_temp_lo[,10]>0)/runs*100
  sens10[x,7] <- sum(cureN_mean[,10]>0)/runs*100
  sens10[x,9] <- sum(cureF_temp_hi[,10]>0)/runs*100
  sens10[x,8] <- sum(cureF_temp_lo[,10]>0)/runs*100
  sens10[x,10] <- sum(cureF_mean[,10]>0)/runs*100
  sens10[x,12] <- sum(death_temp_hi[,10]>0)/runs*100
  sens10[x,11] <- sum(death_temp_lo[,10]>0)/runs*100
  sens10[x,13] <- sum(death_mean[,10]>0)/runs*100
}
# save results
write.csv(sens1, "sens1.csv", row.names=FALSE)
write.csv(sens2, "sens2.csv", row.names=FALSE)
write.csv(sens3, "sens3.csv", row.names=FALSE)
write.csv(sens4, "sens4.csv", row.names=FALSE)
write.csv(sens5, "sens5.csv", row.names=FALSE)
write.csv(sens6, "sens6.csv", row.names=FALSE)
write.csv(sens7, "sens7.csv", row.names=FALSE)
write.csv(sens8, "sens8.csv", row.names=FALSE)
write.csv(sens9, "sens9.csv", row.names=FALSE)
write.csv(sens10, "sens10.csv", row.names=FALSE)

end_time <- Sys.time()
end_time - start_time

sens1   <- fread('sens1.csv')
sens10   <- fread('sens10.csv')

gg <- ggplot(sens1, aes(x=cureN_diff_lo, xend=cureN_diff_hi, y=parameter, group=parameter)) + 
  geom_dumbbell(size_x=1.75, size_xend=1.75, colour_x="#e37371", colour_xend="#a3c4dc") + 
  geom_vline(xintercept = mean(sens10$cureN_diff_mean), linetype="dashed") +
  labs(x="Probability Pan-TB initial cure>SoC", y="parameter") +
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        plot.background=element_rect(fill="#f7f7f7"),
        panel.background=element_rect(fill="#f7f7f7"),
        panel.grid.minor=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.major.x=element_line(),
        axis.ticks=element_blank(),
        legend.position="top",
        panel.border=element_blank())+
  scale_x_continuous(limits = c(0, 100))
plot(gg)
ggsave("Uni_cureN_time10.png",device="png",width=10,height=10,units=c("cm"))


gg <- ggplot(sens1, aes(x=cureF_diff_lo, xend=cureF_diff_hi, y=parameter, group=parameter)) + 
  geom_dumbbell(size_x=1.75, size_xend=1.75, colour_x="#e37371", colour_xend="#a3c4dc") + 
  geom_vline(xintercept = mean(sens1$cureF_diff_mean), linetype="dashed") +
  labs(x="Probability Pan-TB final cure>SoC", y="parameter") +
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        plot.background=element_rect(fill="#f7f7f7"),
        panel.background=element_rect(fill="#f7f7f7"),
        panel.grid.minor=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.major.x=element_line(),
        axis.ticks=element_blank(),
        legend.position="top",
        panel.border=element_blank())+
  scale_x_continuous(limits = c(0, 100))
plot(gg)
ggsave("Uni_cureF_time1_hiprev.png",device="png",width=10,height=10,units=c("cm"))


gg <- ggplot(sens10, aes(x=death_diff_lo, xend=death_diff_hi, y=parameter, group=parameter)) + 
  geom_dumbbell(size_x=1.75, size_xend=1.75, colour_x="#e37371", colour_xend="#a3c4dc") + 
  geom_vline(xintercept = mean(sens10$death_diff_mean), linetype="dashed") +
  labs(x="Probability Pan-TB death>SoC", y="parameter") +
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        plot.background=element_rect(fill="#f7f7f7"),
        panel.background=element_rect(fill="#f7f7f7"),
        panel.grid.minor=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.major.x=element_line(),
        axis.ticks=element_blank(),
        legend.position="top",
        panel.border=element_blank())+
  scale_x_continuous(limits = c(0, 100))
plot(gg)
ggsave("Uni_death_time10.png",device="png",width=10,height=10,units=c("cm"))

gg <- ggplot(sens1, aes(x=resist_diff_lo, xend=resist_diff_hi, y=parameter, group=parameter)) + 
  geom_dumbbell(size_x=1.75, size_xend=1.75, colour_x="#e37371", colour_xend="#a3c4dc") + 
  geom_vline(xintercept = mean(sens1$resist_diff_mean), linetype="dashed") +
  labs(x="Probability Pan-TB resistance>SoC", y="parameter") +
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        plot.background=element_rect(fill="#f7f7f7"),
        panel.background=element_rect(fill="#f7f7f7"),
        panel.grid.minor=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.major.x=element_line(),
        axis.ticks=element_blank(),
        legend.position="top",
        panel.border=element_blank())+
  scale_x_continuous(limits = c(0, 100))
plot(gg)
ggsave("Uni_resist_time1.png",device="png",width=10,height=10,units=c("cm"))
## ================ Heatmap of acq vs imp ====================================
z        <- 1
time     <- 10
imp_split <- seq(0.54, 0.91, 0.037)
acq_split <- seq(0.003, 0.023, 0.002)
samples <- 1000
#pan_succ <- data.frame(X=rep(rr_split,each=length(br_split)), Y=rep(br_split,length(rr_split)), succ=rep(0,length(rr_split)*length(br_split)))
#soc_succ <- data.frame(X=rep(rr_split,each=length(br_split)), Y=rep(br_split,length(rr_split)), succ=rep(0,length(rr_split)*length(br_split)))
tot_succ <- data.frame(X=rep(imp_split,each=length(acq_split)), Y=100*rep(acq_split,length(imp_split)), soc_resist=rep(0,length(imp_split)*length(acq_split)), final=rep(0,length(imp_split)*length(acq_split)), pan_resist=rep(0,length(imp_split)*length(acq_split)), resist=rep(0,length(imp_split)*length(acq_split)))
comp_soc <- rep(0,samples)
comp_final <- rep(0,samples)
comp_pan <- rep(0,samples)
comp_resist <- rep(0,samples)
start_time <- Sys.time()
cp_mean <- c()
cs_mean <- c()
out_mean <-c()
cureN_t_mean <-   c()
cureF_t_mean <-   c()
death_t_mean <-  c()

start_time <- Sys.time()
for (x in imp_split){
  for (y in acq_split){
    prev_soc <- prev
    prev_soc <- prev_soc[rep(seq_len(1), each = time*samples), ]
    prev_soc$cohort <- rep(1:time, samples)
    prev_soc$run <- rep(1:samples, each=time)
    prev_pan <- prev
    prev_pan <- prev_pan[rep(seq_len(1), each = time*samples), ]
    prev_pan$cohort <- rep(1:time, samples)
    prev_pan$run <- rep(1:samples, each=time)
    
    for (samp in 1:samples){
      params_temp <- params_sample(param_table)

      params_temp$S_B <- y
      #params_temp$S_X <- y
      params_temp$P_B <- x
      #params_temp$P_X <- x
      
      
      cp_mean[1] <- prev_pan[(samp-1)*time+1,5]+prev_pan[(samp-1)*time+1,6]+prev_pan[(samp-1)*time+1,8]
      cs_mean[1] <- prev_soc[(samp-1)*time+1,5]+prev_soc[(samp-1)*time+1,6]+prev_soc[(samp-1)*time+1,8]
      out_mean[1] <- cp_mean[1] - cs_mean[1]
      out_s_mean <- cohort_soc(prev, params_temp)
      out_p_mean <- cohort_pan(prev, params_temp)
      #ureN_t_mean[1] <-  out_p_mean$cureN - out_s_mean$cureN
      cureF_t_mean[1] <-  out_p_mean$cureF - out_s_mean$cureF
      #death_t_mean[1] <-  out_p_mean$death - out_s_mean$death
      
      for (t in 2:time){
        #mean values
        out_soc_t <- cohort_soc(prev_soc[(samp-1)*time+t-1,1:8], params_temp) # Outcomes for previous cohort
        out_soc <- sum(out_soc_t[1:8])/(1+sum(out_soc_t[1:8]))
        prev_soc[(samp-1)*time+t,1:8] <- (1-out_soc)*prev_soc[(samp-1)*time+t-1,1:8]+out_soc*out_soc_t[1:8]/sum(out_soc_t[1:8]) # Prevalence is a combination of previous prevalence + outcomes from previous cohort
        out_pan_t <- cohort_pan(prev_pan[(samp-1)*time+t-1,1:8], params_temp) # Outcomes for previous cohort
        out_pan <- sum(out_pan_t[1:8])/(1+sum(out_pan_t[1:8]))
        prev_pan[(samp-1)*time+t,1:8] <- (1-out_pan)*prev_pan[(samp-1)*time+t-1,1:8]+out_pan*out_pan_t[1:8]/sum(out_pan_t[1:8]) # Prevalence is a combination of previous prevalence + outcomes from previous cohort
        
  
        cp_mean[t] <- prev_pan[(samp-1)*time+t,5]+prev_pan[(samp-1)*time+t,6]+prev_pan[(samp-1)*time+t,8]
        cs_mean[t] <- prev_soc[(samp-1)*time+t,5]+prev_soc[(samp-1)*time+t,6]+prev_soc[(samp-1)*time+t,8]
        out_mean[t] <- cp_mean[t] - cs_mean[t]
        
        #cureN_t_mean[t] <-  out_pan_t$cureN - out_soc_t$cureN
        cureF_t_mean[t] <-  out_pan_t$cureF - out_soc_t$cureF
        #death_t_mean[t] <-  out_pan_t$death - out_soc_t$death
      }
      
      comp_resist[samp] <- out_mean[time]
      comp_soc[samp] <- cs_mean[time]
      comp_final[samp] <- cureF_t_mean[time]
      comp_pan[samp] <- cp_mean[time]
    }
    tot_succ$soc_resist[z] <- sum(comp_soc)/samples*100
    tot_succ$final[z] <- sum(comp_final > 0)/samples*100
    tot_succ$pan_resist[z] <- sum(comp_pan)/samples*100
    tot_succ$resist[z] <- sum(comp_resist > 0)/samples*100
    z           <- z+1
  }
}

end_time <- Sys.time()
end_time - start_time

write.csv(tot_succ, "acq_vs_imp_Bonly.csv", row.names=FALSE)
# Plotting treatment cure, i.e. where <50% value indicates SoC is better than pan-TB (more likely less cure due to pan-TB than SoC)
ggplot(tot_succ, aes(X, Y, fill = soc_resist)) +
  geom_tile() +
  scale_fill_gradientn(limits= c(4,9), colors = c("#D0F0C0", "#29AB87", "#0B6623")) +
  ylab("Resistance acquisition rate for B") +
  xlab("Risk ratio of cure given resistance for B") +
  guides(fill=guide_legend(title="Problematic resistance (%)"))
ggsave("Eff_soc_resist_Bonly.png",device="png",width=13,height=10,units=c("cm"))
ggplotot_succt(tot_succ, aes(X, Y, fill = final)) +
  geom_tile() +
  scale_fill_gradientn(limits= c(0,100), colors = c("#075AFF", "#FFFFCC", "#FF0000")) +
  ylab("Resistance acquisition rate for B") +
  xlab("Risk ratio of cure given resistance for B") +
  guides(fill=guide_legend(title="Prob panTB>SoC (%)"))
ggsave("Eff_Final_cure_Bonly.png",device="png",width=13,height=10,units=c("cm"))
# Plotting deaths or resistance, i.e. where >50% value indicates SoC is better than pan-TB (more likely more deaths/resistance due to pan-TB than SoC)
ggplot(tot_succ, aes(X, Y, fill = pan_resist)) +
  geom_tile() +
  scale_fill_gradientn(limits= c(4,9),colors = c("#D0F0C0", "#29AB87", "#0B6623")) +
  ylab("Resistance acquisition rate for B") +
  xlab("Risk ratio of cure given resistance for B") +
  guides(fill=guide_legend(title="Problematic resistance (%)"))
ggsave("Eff_pan_resist_Bonly.png",device="png",width=13,height=10,units=c("cm"))
ggplot(tot_succ, aes(X, Y, fill = resist)) +
  geom_tile() +
  scale_fill_gradientn(limits= c(0,100), colors = c("#FF0000", "#FFFFCC", "#075AFF")) +
  ylab("Resistance acquisition rate for B") +
  xlab("Risk ratio of cure given resistance for B") +
  guides(fill=guide_legend(title="Prob panTB>SoC (%)"))
ggsave("Eff_Final_resist_Bonly.png",device="png",width=13,height=10,units=c("cm"))

## ================ Plot Sankeys  ====================================
soc_figure <- staged_sankey(sankey_data = cohort_soc(prev_data = prev, param_set = param_ext[2,]))
pan_figure <- staged_sankey(sankey_data = cohort_pan(prev_data = prev, param_set = param_ext[2,]))
soc_inset <- staged_sankey_inset(sankey_data = cohort_soc(prev_data = prev, param_set = param_ext[2,]))
pan_inset <- staged_sankey_inset(sankey_data = cohort_pan(prev_data = prev, param_set = param_ext[2,]))

soc_inset
soc_ymax <- ggplot_build(soc_inset$figure)$layout$panel_params[[1]]$y.range[2]
pan_ymax <- ggplot_build(pan_inset$figure)$layout$panel_params[[1]]$y.range[2]


pan_scaled <- pan_inset$figure + ylim(0, pan_ymax * soc_inset$yscale / pan_inset$yscale)

(soc_figure + pan_figure) / 
  (soc_inset$figure + pan_scaled)

ggsave("Sankey_full.png",device="png",width=100,height=70,units=c("cm"))
