extrafont::loadfonts()



## ----'Effect Size Calculation'-------------------------------------------

es <- function(SDe.base, SDe.final, Me.base, Me.final, r) {
  ## ----------------------- Thu Apr 04 18:51:29 2019 ------------------------# Standard Deviation of the Mean Change Score
#https://handbook-5-1.cochrane.org/chapter_16/16_1_3_2_imputing_standard_deviations_for_changes_from_baseline.htm
 SDe.change <- sqrt(SDe.base^2 + SDe.final^2 - (2 * r * SDe.base * SDe.final))
d <- abs(Me.final - Me.base) / SDe.change
  return(d)
}

# ----------------------- Sat Oct 05 15:40:44 2019 ------------------------#
# Independent Groups Pre-test/Post-Test Effect Size Calculation
# (Morris & DeShon, 2002) p 108

d_igpp <- function(m_post_e, m_pre_e, m_post_c, m_pre_c, s_pre_e, s_pre_c) {
  d_igpp <- (m_post_e - m_pre_e) / s_pre_e - (m_post_c - m_pre_c) / s_pre_c
  return(d_igpp)
}

## ----'SD Calculation'----------------------------------------------------
# ----------------------- Fri Apr 19 16:16:54 2019 ------------------------#
# SD & Confidence Interval Function from Effect size
# From Lee, D. K. (2016). Alternatives to P value: Confidence interval and effect size. Korean Journal of Anesthesiology, 69(6), 555. https://doi.org/10.4097/kjae.2016.69.6.555 p 559

dCI <- function(n1, n2, d, a = .05, tails = 2) {
  sigma_d <- sqrt((n1 + n2)/(n1 * n2) + d ^ 2 / (2 * (n1 + n2)))
  z <- abs(qnorm(a / tails))
  return(c(cilb = {d - z * sigma_d}, ciub = d + z * sigma_d))
}


Expanded <- list()

# ----'BergenCico'------------------------------------------#
# Bergen-Cico, D., Possemato, K., & Cheon, S. (2013). Examining the Efficacy of a Brief Mindfulness-Based Stress Reduction (Brief MBSR) Program on Psychological Health. Journal of American College Health, 61(6), 348–360. https://doi.org/10.1080/07448481.2013.813853
# p 354
Expanded$CT$`Bergen-Cico`$pre_treatment <- list(mt = 84.61, # treatment mean
                 mc = 82.62, # control mean
                 st = 17, # treatment sd
                 sc = 12, # control sd
                 nt = 72, # treatment pop. size
                 nc = 47) # control pop.size

Expanded$CT$`Bergen-Cico`$post_treatment <- list(mt = 88.57,
                       mc = 81.7,
                       st = 16,
                       sc = 12,
                       nt = 72,
                       nc = 47)
Expanded$CT$`Bergen-Cico`$N <- c(nt = 72, nc = 47)
Expanded$CT$`Bergen-Cico`$f <- 7.3


# ----------------------- Sat Oct 05 16:09:07 2019 ------------------------#
# BergenCico Between Group

# The R value "Covariate outcome correlation or multiple Correlation" at the level of the individual participants was not provided upon request. This value is computed from the correlation between the individual components of the SCS.
# Table 2, p 354
Expanded$TT$`Bergen-Cico`$subscales <- tibble::tribble(
  ~pre_t_M, ~pre_c_M, ~post_t_M, ~post_c_M,
  12, 12, 13, 12, # common humanity
  14, 13, 14, 13, # isolation
  17, 16, 17, 15, # self-judgment
  16, 15, 16, 15, # self-kindness
  13, 13, 13, 12, # mindfulness
  13, 13, 14, 13 # over-identification
) 
## ----'BergenCico Within Group'-------------------------------------------
# We don't have the correlation between the pre/post for the treatment group so according to the Cochrane handbook is must be imputed
#  When baseline and final standard deviations are known, we can impute the missing standard deviation using an imputed value, Corr, for the correlation coefficient. The value Corr might be imputed from another study in the meta-analysis (using the method in (1) above), it might be imputed from elsewhere, or it might be hypothesized based on reasoned argument. In all of these situations, a sensitivity analysis should be undertaken, trying different values of Corr, to determine whether the overall result of the analysis is robust to the use of imputed correlation coefficients.
# ----------------------- Sun Oct 05 08:46:42 2019 ------------------------#
# Attempt to use subscale data to replicate results [DID NOT WORK]
# Transform the data for subscales entered as it appears from the table to how it was entered in the ANCOVA as stated on p 352 in the "Analysis" subsection
# Group by pre - treatment with group as covariate
pre <- Expanded$TT$`Bergen-Cico`$subscales %>% tidyr::gather(key = "group", value = "pre", pre_t_M, pre_c_M) %>% dplyr::mutate_at(dplyr::vars(group), ~ stringr::str_match(., "^[a-z]+\\_(\\w)")[,2]) %>% dplyr::select(group, pre)
# Group by post - treatment with group as covariate
post <- Expanded$TT$`Bergen-Cico`$subscales %>% tidyr::gather(key = "group", value = "post", post_t_M, post_c_M) %>% dplyr::mutate_at(dplyr::vars(group), ~ stringr::str_match(., "^[a-z]+\\_(\\w)")[,2]) %>% dplyr::select(group, post)
# Combine the two (groups are in the same order so cbind is sufficient)
Expanded$TT$`Bergen-Cico`$subscales <- do.call("cbind.data.frame", list(group = pre$group, pre = pre$pre, post = post$post))
# The covariate outcome correlation or multiple corel
psych::setCor(
  data = Expanded$TT$`Bergen-Cico`$subscales, x = c("pre"), y = c("post"), z = c("group"), plot = F
) # Does not match the results in the study likely due to the degrees freedom
# End attempt to use subscale data
#----------------------- Sun Oct 05 08:47:21 2019 ------------------------#

# Here we will use the correlation between the subscale scores to imput the correlation
r_prepost <- cor(Expanded$TT$`Bergen-Cico`$subscales$pre, Expanded$TT$`Bergen-Cico`$subscales$post)

(Expanded$TT$`Bergen-Cico`$d <- eval(rlang::call2(.fn = es, !!!list(
  SDe.base = Expanded$CT$`Bergen-Cico`$pre_treatment$st, 
  SDe.final = Expanded$CT$`Bergen-Cico`$post_treatment$st, 
  Me.base = Expanded$CT$`Bergen-Cico`$pre_treatment$mt,
  Me.final = Expanded$CT$`Bergen-Cico`$post_treatment$mt,
  r = r_prepost))))
Expanded$TT$`Bergen-Cico`$N <- 72
(Expanded$TT$`Bergen-Cico`$CI <- dCI(n1 = 72, n2 = 72, Expanded$TT$`Bergen-Cico`$d))



# ----------------------- Sun Oct 06 08:41:35 2019 ------------------------#
# Attempt to replicate the distribution of results that Bergen-Cico likely had
# A formula to generate a distribution with exact mean and SD
rnorm. <- function(n=10, mean=0, sd=1, fixed=TRUE) { switch(fixed+1, rnorm(n, mean, sd), as.numeric(mean+sd*scale(rnorm(n)))) }

# Generate sample distriubtion that mirrors what Bergen-Cico would have had
mock_data <- list()
mock_data$pre <- c(rnorm.(n = Expanded$CT$`Bergen-Cico`$pre_treatment$nt,
                          mean = Expanded$CT$`Bergen-Cico`$pre_treatment$mt,
                          sd = Expanded$CT$`Bergen-Cico`$pre_treatment$st,
                          fixed = T), # The pre-treatment treatment group distribution
                   rnorm.(n = Expanded$CT$`Bergen-Cico`$pre_treatment$nc,
                          mean = Expanded$CT$`Bergen-Cico`$pre_treatment$mc,
                          sd = Expanded$CT$`Bergen-Cico`$pre_treatment$sc,
                          fixed = T)) # The pre-treatment control group distribution)

mock_data$post <- c(rnorm.(n = Expanded$CT$`Bergen-Cico`$post_treatment$nt,
                           mean = Expanded$CT$`Bergen-Cico`$post_treatment$mt,
                           sd = Expanded$CT$`Bergen-Cico`$post_treatment$st,
                           fixed = T), # the post-treatment treatment group distribution
                    rnorm.(n = Expanded$CT$`Bergen-Cico`$post_treatment$nc,
                           mean = Expanded$CT$`Bergen-Cico`$post_treatment$mc,
                           sd = Expanded$CT$`Bergen-Cico`$post_treatment$sc,
                           fixed = T))
mock_data$group <- c(rep("t", Expanded$CT$`Bergen-Cico`$post_treatment$nt), # the group covariate,, treatment level
                     rep("c", Expanded$CT$`Bergen-Cico`$post_treatment$nc) # the group covariate control 
)
mock_data <- do.call("cbind.data.frame", mock_data)
# Show the means & SDs to ensure the distribution mirrors the actual data
tapply(mock_data$pre, mock_data$group, FUN = mean)
tapply(mock_data$pre, mock_data$group, FUN = sd)
tapply(mock_data$post, mock_data$group, FUN = mean)
tapply(mock_data$post, mock_data$group, FUN = sd)




# Multiple correlation to calculate R as input to compute.es::a.fes. This will not be exact because the correlation is going to be different from the actual distribution of results.
multiple_correlation <- psych::setCor(
y = post ~ group + pre, data = mock_data
, plot = F)

aov(post ~ group + pre, data = mock_data) # The F statistic is close?
ex <- compute.es::mes(m.1 = 88.57, m.2 = 81.7, sd.1 = 16, sd.2 = 12, n.1 = 72, n.2 = 47, verbose = F) # The effect size is relatively close estimate?

# Compute effect size and Confidence interval for input into the fixed-effects model metafor model
BergenCT_tes <- invisible({compute.es::a.fes(f = Expanded$CT$`Bergen-Cico`$f, # F-statistic for composite SCS score
                                n.1 = Expanded$CT$`Bergen-Cico`$pre_treatment$nt, # treatment pop.size
                                n.2 = Expanded$CT$`Bergen-Cico`$post_treatment$nc, # control pop.size
                                R = r_prepost,
                                q = 1,
                                verbose = F) %>% setNames(names(ex))}) # Setnames because the function doesn't name the outputs correctly
# Save the relevant values
Expanded$CT$`Bergen-Cico`$g <- BergenCT_tes$g
Expanded$CT$`Bergen-Cico`$CI <- c(cilb = BergenCT_tes$l.g,
                                  ciub = BergenCT_tes$u.g)




## ----'Dundas CT'---------------------------------------------------------
#Dundas, I., Binder, P.-E., Hansen, T. G. B., & Stige, S. H. (2017). Does a short self-compassion intervention for students increase healthy self-regulation? A randomized control trial. Scandinavian Journal of Psychology, 58(5), 443–450. https://doi.org/10.1111/sjop.12385
#pg 447
Expanded$CT$Dundas$N <-  c(nt = 53, nc = 64)
# Reported
Expanded$CT$Dundas$d <- .94
Expanded$CT$Dundas$CI <- c(cilb = .66,
ciub = 1.23)
Expanded$CT$Dundas$a <- .05 #Alpha level

## ----'Load Dundas PrePost', eval = F-------------------------------------
## (Dundas <- utils::readClipboard() %>% .[!nchar(.) < 1] %>% matrix(nrow = 7, byrow = T) %>% as.data.frame %>% setNames(c("Pair", "Desc", "N","r","sig")) %>% mutate_at(vars(r),funs(gsub, .args = list(pattern = "\\,", replacement = "."))) %>% mutate_at(vars(r,N), funs(as.numeric(as.character(.))))  %>% write.csv(file = "DundasCorrelations.csv"))

## ----'Load Dundas PrePost from CSV'--------------------------------------
# Retreived via email correspondence
Dundas  <- data.frame(
                                                                                                        X = c(1L, 2L, 3L, 4L, 5L, 6L, 7L),
                                                                                                        N = c(53L, 53L, 52L, 53L, 53L, 53L, 53L),
                                                                                                        r = c(0.637, 0.584, 0.656, 0.722, 0.61, 0.728, 0.477),
                                                                                                     Pair = as.factor(c("Pair 1", "Pair 2", "Pair 3", "Pair 4", "Pair 5",
                                                                                                                        "Pair 6", "Pair 7")),
                                                                                                     Desc = as.factor(c("T2SCSum & T3SCSum", "T2JudgeSum & T3JudgeSum",
                                                                                                                        "SumT2HINT & SumT3HINT",
                                                                                                                        "SumT2C & SumT3C", "SumT2PGIS & SumT3PGIS",
                                                                                                                        "SumT2STAI & SumT3STAI",
                                                                                                                        "SumT2MDI & SumT3MDI")),
                                                                                                      sig = as.factor(c(".000", ".000", ".000", ".000", ".000", ".000", ",
                                                                                                                       000"))
                                                                                             )

Expanded$TT$Dundas$N <- c(nt = Dundas[stringr::str_detect(Dundas$Desc, "SC"), "N", drop = T])
## ----'Dundas TT eff size and CI'-------------------------------------------
(Expanded$TT$Dundas$d <- es(SDe.base = .62, # Baseline 2 pg 447
                            Me.base = 2.8, # Baseline 2 pg 447
                            SDe.final = .57, # Post-intervention pg 447
                            Me.final = 3.35, # Post-intervention pg 447
                            r = Dundas[stringr::str_detect(Dundas$Desc, "SC"), "r", drop = T]))
(Expanded$TT$Dundas$CI <- dCI(n1 = Dundas[stringr::str_detect(Dundas$Desc, "SC"), "N", drop = T],
                              n2 = Dundas[stringr::str_detect(Dundas$Desc, "SC"), "N", drop = T],
                              d = Expanded$TT$Dundas$d))


## ----'Erogul CT'---------------------------------------------------------
# Erogul, M., Singer, G., McIntyre, T., & Stefanov, D. G. (2014). Abridged Mindfulness Intervention to Support Wellness in First-Year Medical Students. Teaching and Learning in Medicine, 26(4), 350–356. https://doi.org/10.1080/10401334.2014.945025

Expanded$CT$Erogul$N <- c(nt = 29, # table 1, pg 352
                          nc = 28)
# Reported pg 353
Expanded$CT$Erogul$d <- .58
Expanded$CT$Erogul$CI <- c(cilb = .23,
                 ciub = .92)
Expanded$CT$Erogul$a <- .05


## ----'Mathad TT'-----------------------------------------------------------
# Mathad, M. D. (2017). Effect of Yoga on Psychological Functioning of Nursing Students: A Randomized Wait List Control Trial. JOURNAL OF CLINICAL AND DIAGNOSTIC RESEARCH. https://doi.org/10.7860/JCDR/2017/26517.9833

Expanded$TT$Mathad$N <- 40 # Table 3, pg 3
Expanded$TT$Mathad$r <- .497 # Provided via email correspondence
(Expanded$TT$Mathad$d <- es(SDe.base = 0.46, # Table 4 Self-compassion /Yoga
                            Me.base = 3.03,
                            SDe.final = 0.28,
                            Me.final = 3.19,
                            r = Expanded$TT$Mathad$r))
(Expanded$TT$Mathad$CI <- dCI(n1 = 40, n2 = 40, d = Expanded$TT$Mathad$d))


## ----'Smeets CT'---------------------------------------------------------
# Smeets, E., Neff, K., Alberts, H., & Peters, M. (2014). Meeting Suffering With Kindness: Effects of a Brief Self-Compassion Intervention for Female College Students: Self-Compassion Intervention for Students. Journal of Clinical Psychology, 70(9), 794–807. https://doi.org/10.1002/jclp.22076

Expanded$CT$Smeets$N <- c(nc = 25, # Table 2 pg 8 Self-compassion
                          nt = 27)
Expanded$CT$Smeets$d <- 1.19
(Expanded$CT$Smeets$CI <- dCI(n1 = 25, n2 = 27, d = 1.19))
Expanded$CT$Smeets$a <- .05

## ----'Hedges Correction and Standard Error from CI'----------------------
# The number of studies & Total N

Studies <- list()
Studies$TT <- purrr::imap(.x = Expanded$TT,.f = function(.x, .y){
  
    N <- unlist(purrr::pluck(.x, c("N"))) %>% setNames(nm = NULL)
    d <- purrr::pluck(.x, c("d"))
    
    # Hedge's G correction for pre-post design
    # Hedges, L., & Olkin, I. (1985). Statistical Methods in Meta-Analysis. In Stat Med (Vol. 20). https://doi.org/10.2307/1164953
    # https://stats.stackexchange.com/questions/1850/difference-between-cohens-d-and-hedges-g-for-effect-size-metrics
    df <- 2 * N - 2
    .x$g <- d * (gamma(df / 2) / (sqrt(df / 2) * gamma((df - 1) / 2))) 
    # Standard Error from CI
    CI <- purrr::pluck(.x, c("CI"))
    .x$se <- (CI[["ciub"]] - CI[["cilb"]]) / 3.92 # If alpha is .05
    return(.x)
  })
Studies$CT <- purrr::imap(.x = Expanded$CT,.f = function(.x, .y){
  
  N <- unlist(purrr::pluck(.x, c("N"))) %>% setNames(nm = NULL)
  # print(N)
  d <- purrr::pluck(.x, c("d"))
  
  # Hedge's G correction for pre-post design
  # Hedges, L., & Olkin, I. (1985). Statistical Methods in Meta-Analysis. In Stat Med (Vol. 20). https://doi.org/10.2307/1164953
  df <- sum(N) - 2
  if (!HDA::go(.x$g)) {
  .x$g <- d * (gamma(df / 2) / (sqrt(df / 2) * gamma((df - 1) / 2))) 
  }
  # Standard Error from CI
  CI <- purrr::pluck(.x, c("CI"))
  .x$se <- (CI[["ciub"]] - CI[["cilb"]]) / 3.92 # If alpha is .05
  return(.x)
})




## ----'Get Results Table', results = 'asis'-------------------------------
# googlesheets::gs_auth(token = "~//R//husky_googlesheets_token.rds")
# gsRes %<>% .[["sheet_key"]] %>% googlesheets::gs_key()
# Res_tbl <- gsRes %>% googlesheets::gs_read(.,ws = "CT Table",trim_ws=T,col_names=T)
# The transformed Results table Res_tbl is below
CT_Table  <- tibble::tribble(
  ~ Study, ~ Year, ~ Population, ~ Type, ~`Contact Time (Min)`, ~ `Ind Time (Min)`, ~ `Duration (Wks)`, ~ `Int Type`, ~ `Inperson/Online`, ~ `Treatment Intervention`, ~ `Control Intervention`, ~ Nint, ~ Ncon, ~ `Hedges' g (CT)`,
  "Bergen-Cico, Possemato, Cheon", 2013, "Undergraduate students", "Quasi-experimental pretest/posttest", 600, 0, 5, "MBSR", "I", "brief (5-week) mindfulness-based stress reduction (brief MBSR)", "a didactic lecture course that met once per week for approximately 2.5 hours.", 72, 47, "0.58 [0.19, 0.95]",
  "Dundas, Binder, Hansen, Stige", 2017, "University students", "Randomized controlled trial", 270, 165, 2, "SCC", "I", "two-week self-compassion course - three 90-minute sessions delivered over a period of two weeks. Course was developed based around previous mindfulness interventions", "waitlist control; participants were offered the intervention one week after the intervention-group had completed theirs", 53, 64, "0.94 [0.66, 1.22]",
  "Erogul, Singer, McIntyre, & Stefanov",2014,"First-year medical students", "Randomized controlled trial",900,960,8,"MBSR","I","8-week MBSR intervention. Full day meditation retreat between week 7 and week 8.","no intervention",28,29,"0.57 [0.23,0.91]",
  "Smeets, Neff, Alberts & Peters", 2014, "Undergraduate females", "Randomized controlled trial", 225, 60, 3, "SCC", "I", "3-week Self-compassion course - two 90 minute sessions and a final 45 minute session.", "3-week time management course - two 90 minute sessions and a final 45 minute session.", 27, 25, "1.17 [0.60, 1.78]")
(rma_mods <- CT_Table %>% names %>% .[c(5,6,7,8)])
CT_Table %<>% mutate_at(vars(rma_mods[4]),as.factor)
## ----'Build RMA model for Control Treatment'-----------------------------
library(metafor)
# Sort Studies Alphabetically
Studies$CT %<>% .[sort(Studies$CT %>% names)]
# General model - no moderators
Studies$Results$CT$General <- metafor::rma(yi = sapply(Studies$CT,`[[`,"g"),sei = sapply(Studies$CT,`[[`,"se"), method = "HS", weighted = F, slab = names(Studies$CT), main = "Test") %>% summary
# Intervention Type as moderators
Studies$Results$CT$Type <- metafor::rma(yi = sapply(Studies$CT,`[[`,"g"),sei = sapply(Studies$CT,`[[`,"se"), method = "HS", weighted = F, mods = ~ `Int Type`, data = CT_Table, slab = paste(names(Studies$CT), CT_Table$`Int Type`,sep = ",")) %>% summary
# Contact & Individual Time as moderators
Studies$Results$CT$Time <- metafor::rma(yi = sapply(Studies$CT,`[[`,"g"),sei = sapply(Studies$CT,`[[`,"se"),method = "HS", weighted = F,mods = paste0("~ `", rma_mods[1], "` + `",rma_mods[2],"`") %>% as.formula, data = CT_Table, slab = paste0(names(Studies$CT)," CT:", CT_Table$`Contact Time (Min)`,"m,"," IT:",CT_Table$`Ind Time (Min)`,"m")) %>% summary
Studies$Results$CT$Time$order <- CT_Table %>% rownames_to_column() %>% arrange(CT_Table$`Contact Time (Min)`) %>% .[,"rowname", drop = T] %>% as.numeric
# Duration (Wks) as moderator
Studies$Results$CT$Duration <- metafor::rma(yi = sapply(Studies$CT,`[[`,"g"),sei = sapply(Studies$CT,`[[`,"se"), method = "HS", weighted = F, mods = paste0("~ `", rma_mods[3], "`") %>% as.formula, data = CT_Table, slab = paste(names(Studies$CT), CT_Table$`Duration (Wks)`,"Wks",sep = ","))
Studies$Results$CT$Duration$order <- CT_Table %>% rownames_to_column() %>% arrange(CT_Table$`Duration (Wks)`) %>% .[,"rowname", drop = T] %>% as.numeric
CT_Table[ ,12:13] %>% colSums() %>% Reduce(f=`+`,x=.)


## ----'Moderators CT', results = 'asis'-----------------------------------
purrr::map2(.x = Studies$Results$CT, .y = c("General", "Moderator: Intervention Type", "Moderator: Contact Time & Individual Time", "Moderator: Duration in Weeks"), .f = function(.x, .y){
  tags$strong(.y)
  # Forest Plot
  if(!is.null(.x[["order"]])) o <- .x[["order"]] else  o <- "obs"
  #grDevices::pdf(paste0(gsub("\\:","",.y) %>% gsub("^","CT",.),".pdf"), family="Merriweather", width=16/2, height=9/2)
.out <- metafor::forest(.x,order = o, steps = 5, mlab = "Hunter-Schmidt RE Estimator", main = "Control-Treatment Hunter-Schmidt Random Effects Model")
dev.off()
  # Table
  .x %>% capture.output %>% gsub("$", " \n", .) %>% cat %>% kableExtra::kable("latex",booktabs = T) %>% kableExtra::kable_styling(position = "center")
  as.data.frame(.x$beta) %>% cbind(.x$se, .x$zval, .x$pval %>% HDA::p.txt(), .x$ci.lb, .x$ci.ub) %>% rownames_to_column() %>% setNames(c("Moderator", "Beta", "St.Err", "Z Value","P Value", "CI lower", "CI upper")) %>%  mutate_at(vars(c(-1,-5)), funs(round(.,3))) %>% kableExtra::kable("latex",booktabs = T) %>% kableExtra::kable_styling(position = "center")
  return(.out)
})


## ----'Get Results TT Table'----------------------------------------------
# googlesheets::gs_auth(token = "~//R//husky_googlesheets_token.rds")
# gs <- googlesheets::gs_url("https://docs.google.com/spreadsheets/d/1dq_9WL8PY3B8EGei3GdlKwxPfD1KlcegmC0QVfI70cI/edit#gid=2047158368&range=A1:N5")
# TT_Table <- googlesheets::gs_read(gs, ws = "TT Table", range = "A1:N5")
TT_Table <- tibble::tribble(
  ~ Study, ~ Year, ~ Population, ~ Type, ~`Contact Time (Min)`, ~ `Ind Time (Min)`, ~ `Duration (Wks)`, ~ `Int Type`, ~ `Inperson/Online`, ~ `Treatment Intervention`, ~ `Control Intervention`, ~ Nint, ~ Ncon, ~ `Hedges' g (CT)`,
  "Bergen-Cico, Possemato, Cheon, 2013",  2013,  "Undergraduate students", "Quasi-experimental pretest/posttest",                 600,                      0,               5,    "MBSR",              "I",                                                                                                        "brief (5-week) mindfulness-based stress reduction (brief MBSR)", "a didactic lecture course that met once per week for approximately 2.5 hours.",              72,         47,              "1.23 [0.89,1.60]",
  "Dundas, Binder, Hansen, Stige, 2017",  2017,     "University students",                                 "RCT",                 270,                    165,               2,     "SCC",              "I", "two-week self-compassion course - three 90-minute sessions delivered over a period of two weeks. Course was developed based around previous mindfulness interventions",                                                              "waitlist control",              53,         64,              "1.07 [0.67,1.49]",
  "Mathad, 2017",  2017, "Female nursing students",                                 "RCT",                2400,                      0,               8,    "YOGA",              "I",                                                                                                          "Integrated approach to yoga \ntherapy as designed by S-VYASA",                                                              "waitlist control",              40,         40,             "0.39 [-0.05,0.84]"
)
(rma_mods <- TT_Table %>% names %>% .[c(5,6,7,8)])
TT_Table %<>% mutate_at(vars(rma_mods[4]),as.factor)
## ----'TT RMA Model'------------------------------------------------------
Studies$TT %<>% .[sort(Studies$TT %>% names)]
# General model - no moderators
Studies$Results$TT$General <- metafor::rma(yi = sapply(Studies$TT,`[[`,"g"),sei = sapply(Studies$TT,`[[`,"se"), method = "FE", weighted = F, slab = names(Studies$TT), main = "Test") %>% summary
# Intervention Type as moderators
Studies$Results$TT$Type <- metafor::rma(yi = sapply(Studies$TT,`[[`,"g"),sei = sapply(Studies$TT,`[[`,"se"), method = "FE", weighted = F, mods = ~ `Int Type`, data = TT_Table, slab = paste(names(Studies$TT), TT_Table$`Int Type`,sep = ",")) %>% summary
# Contact & Individual Time as moderators
Studies$Results$TT$Time <- metafor::rma(yi = sapply(Studies$TT,`[[`,"g"),sei = sapply(Studies$TT,`[[`,"se"),method = "FE", weighted = F,mods = paste0("~ `", rma_mods[1], "` + `",rma_mods[2],"`") %>% as.formula, data = TT_Table, slab = paste0(names(Studies$TT), "CT: ", TT_Table$`Contact Time (Min)`,"m,","IT: ",TT_Table$`Individual Time (Min)`,"m")) %>% summary
Studies$Results$TT$Time$order <- TT_Table %>% rownames_to_column() %>% arrange(TT_Table$`Contact Time (Min)`) %>% .[,"rowname", drop = T] %>% as.numeric
# Duration (Wks) as moderator
Studies$Results$TT$Duration <- metafor::rma(yi = sapply(Studies$TT,`[[`,"g"),sei = sapply(Studies$TT,`[[`,"se"), method = "FE", weighted = F, mods = paste0("~ `", rma_mods[3], "`") %>% as.formula, data = TT_Table, slab = paste(names(Studies$TT), TT_Table$`Duration (Wks)`,"Wks",sep = ","))
Studies$Results$TT$Duration$order <- TT_Table %>% rownames_to_column() %>% arrange(TT_Table$`Duration (Wks)`) %>% .[,"rowname", drop = T] %>% as.numeric
TT_Table[ ,12:13] %>% colSums() %>% Reduce(f=`+`,x=.)


## ----'Forest Plot'-------------------------------------------------------
#grDevices::pdf("ForestTT.pdf", family="Merriweather", width=16/2, height=9/2)
metafor::forest(Studies$Results$TT$General, order = "obs", refline = Studies$Results$TT$General$beta, steps = 5, mlab = "Fixed Effects Estimator", main = "Pre/Post Treatment Fixed Effects Model")
dev.off()


## ----'Moderators PrePost', results = 'asis'------------------------------
purrr::map2(.x = Studies$Results$TT, .y = c("General", "Moderator: Intervention Type", "Moderator: Contact Time & Individual Time", "Moderator: Duration in Weeks"), .f = function(.x, .y){
  tags$strong(.y)
  # Forest Plot
  if(!is.null(.x[["order"]])) o <- .x[["order"]] else  o <- "obs"
  #grDevices::pdf(paste0(gsub("\\:","",.y) %>% gsub("^","TT",.),".pdf"), family="Merriweather", width=16/2, height=9/2)
.out <- metafor::forest(.x, order = o, steps = 5, mlab = "Fixed Effects Estimator", main = "Pre/Post Treatment Fixed Effects Model")
dev.off()
  # Table
  .x %>% capture.output %>% gsub("$", " \n", .) %>% cat %>% kableExtra::kable("latex",booktabs = T) %>% kableExtra::kable_styling(position = "center")
  as.data.frame(.x$beta) %>% cbind(.x$se, .x$zval, .x$pval %>% HDA::p.txt(), .x$ci.lb, .x$ci.ub) %>% rownames_to_column() %>% setNames(c("Moderator", "Beta", "St.Err", "Z Value","P Value", "CI lower", "CI upper")) %>%  mutate_at(vars(c(-1,-5)), funs(round(.,3))) %>% kableExtra::kable("html",booktabs = T) %>% kableExtra::kable_styling(position = "center")
return(.out)
  })


