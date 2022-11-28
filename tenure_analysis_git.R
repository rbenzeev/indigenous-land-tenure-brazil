# load libraries
library(tidyverse)
library(broom)
library(plm)
library(rddtools)
library(rdrobust) # for rdplot
library(plyr) # to create 'mu' object in EDA
library(MatchIt)
library(optmatch)
library(lmtest) #coeftest
library(sandwich) #vcovCL
library(did)

library(DIDmultiplegt)
library(panelView)
library(patchwork)
library(kableExtra)

theme_set(theme_light())

# READ IN DATA ====
# read in full dataset
# subset data to exclude 'Indigenous reserves' and only include 'Indigenous territories' 

forest_change <- read.csv("Data/Excel_spreadsheets/data_file.csv") %>% 
  select(-X, -X.1) %>% 
  filter(TI_RI_DI == "TI")

# add additional columns
forest_change <- mutate(forest_change,
                        Year_hom_0 = replace_na(Year_hom, 0),
                        Year_dec_0 = replace_na(Year_dec, 0),
                        Perc_forest = Forest / Size,
                        Perc_forest2 = Forest/ (Size-Forest),
                        Log_forest = log(Forest),
                        Log_size = log(Size),
                        Log_perc_forest = log(Perc_forest),
                        Forest_change_sign = Forest_change / abs(Forest_change),
                        Log_forest_change = Forest_change_sign * log(abs(Forest_change))
)


# DIFFERENCE-IN-DIFFERENCES WITH MULTIPLE TIME PERIODS ====
#
# Some recent advances in econometrics make this easy for us. We follow
# Callaway and Sant'Anna (2020), who developed a convenient package as well,
# called DID (https://bcallaway11.github.io/did/articles/multi-period-did.html)
# that allows us to perform a difference-in-differences with multiple time periods
# and staggered treatment adoptions.
#
# We have two main choices: we can choose to compare treated groups (territories
# that get tenure) vs. never-treated groups (territories that never get tenure)
# or we can compare treated groups to never-treated AND pre-treated groups.
# We'll run both methods, as a robustness check, as well as both forest change
# and forest change percentage as dependent variables.
#
# Since this can run on data with more than two periods, we'll go back
# to using the original dataframe

att_wrapper <- function(data, yname, gname) {
  forest_change_attgt_dec <- att_gt(
    yname=yname,
    tname="Year",
    idname="Codigo",
    gname=gname,
    allow_unbalanced_panel = FALSE,
    control_group="notyettreated",
    xformla = ~ Log_size,
    data=data,
    anticipation=1
  )
}

plot_dynamic_att <- function(attgt_es) {
  attgt_es_df <-data.frame(
    egt=attgt_es$egt,
    att.egt=attgt_es$att.egt,
    se.egt=attgt_es$se.egt,
    att.egt.05=attgt_es$att.egt - 1.98 * attgt_es$se.egt,
    att.egt.95=attgt_es$att.egt + 1.98 * attgt_es$se.egt
  )
  
  ggplot(attgt_es_df) +
    geom_line(mapping=aes(x=egt, y=att.egt), color="#661100",size=1.5) +
    geom_ribbon(mapping=aes(x=egt, ymin=att.egt.05, ymax=att.egt.95), fill="#88CCEE", alpha=0.5) +
    geom_hline(yintercept=0, size=2, alpha=0.1) +
    geom_vline(xintercept=0, size=2, alpha=0.1)
}

forest_change_attgt <- att_wrapper(forest_change, "Forest_change_perc", "Year_hom_0")

summary(forest_change_attgt)
forest_change_attgt_es <- aggte(forest_change_attgt, type = "dynamic", min_e=-10, max_e = 10, na.rm=TRUE)
summary(forest_change_attgt_es)

# original y axis: "Estimated ATT of % deforestation
dynamic_att <- plot_dynamic_att(forest_change_attgt_es) +
    xlab("Number of years to land tenure") +
    ylab("Percent forest change") +
    theme_bw()

forest_change_attgt_overall <- aggte(forest_change_attgt, type = "simple", na.rm=TRUE)
summary(forest_change_attgt_overall)
str(forest_change_attgt_overall)

df.atts <- data.frame(
  estimate=forest_change_attgt_overall$overall.att,
  std.error=forest_change_attgt_overall$overall.se,
  method="CS"
)

forest_change_attgt_group <- aggte(forest_change_attgt, type = "group", na.rm=TRUE)
summary(forest_change_attgt_group)

# SUN AB ESTIMATOR ====
#

#mod.sunab <- feols(Forest_change_perc ~ factor(Codigo) + sunab(cohort=Year_hom_0, period=Year, ref.p=c(-4, -3)), forest_change)
mod.sunab <- feols(Forest_change_perc ~ factor(Codigo) + sunab(cohort=Year_hom_0, period=Year, ref.p=c(-2, -1)), forest_change)
summary(mod.sunab)
df.atts <- rbind(df.atts, tidy(summary(mod.sunab, agg="ATT"))  %>% 
  filter(term == "ATT") %>%
  select(c("estimate", "std.error")) %>%
  mutate(method="sunab")
)
  
summary(mod.sunab, agg="cohort")
iplot(mod.sunab)

mod.sunab.fc <- feols(Forest_change ~ factor(Codigo) + sunab(cohort=Year_hom_0, period=Year, ), forest_change)
summary(mod.sunab.fc)
summary(mod.sunab.fc, agg="ATT")
iplot(mod.sunab.fc)

# DIDMULTIPLEGT ====
#
# The approach of de Chaisemartin and D'Haultfoeuille (2020)
# https://www.aeaweb.org/doi/10.1257/aer.20181169
#
# We run the corresponding placebo test as well, to test the common trend
# assumption

forest_change_didmgt <- did_multiplegt(
  df=forest_change,
  Y="Forest_change_perc",
  G="Codigo",
  T="Year",
  D="Tenure",
  brep=100
)

df.atts <- rbind(df.atts, data.frame(list(
  estimate=forest_change_didmgt$effect,
  std.error=forest_change_didmgt$se_effect,
  method="dC-DH"
)))

forest_change_didmgt_placebo <- did_multiplegt(
  df=forest_change,
  Y="Forest_change_perc",
  G="Codigo",
  T="Year",
  D="Tenure",
  placebo=1,
  brep=100,
)

forest_change_didmgt_placebo

# DECLARED TERRITORIES ONLY ANALYSIS ==== 
#
# We're going to look at the declared territories in treatment and in control
# and only at the years after their declaration

forest_change_declared_only <- forest_change %>%
  filter(Year_dec_0 != 0 & Year_hom_0 == 0) %>%
  mutate(Year_since_dec = Year - Year_dec,
         Declared = as.numeric(Year > Year_dec))

forest_change_didmgt_declared <- did_multiplegt(
  df=forest_change_declared,
  Y="Forest_change_perc",
  G="Codigo",
  T="Year",
  D="Declared",
  brep=100
)

df.atts.dec <- data.frame(list(
  dec.estimate=forest_change_didmgt_declared$effect,
  dec.std.error=forest_change_didmgt_declared$se_effect,
  method="dC-DH"
))

forest_change_didmgt_declared

mod.sunab.declared.only <- feols(
  Forest_change_perc ~ factor(Codigo) + sunab(cohort=Year_dec, period=Year, ref.p=c(-2, -1)), forest_change_declared_only
)
summary(mod.sunab.declared.only)
summary(mod.sunab.declared.only, agg="ATT")
iplot(mod.sunab.declared.only)

tidysun <- tidy(summary(mod.sunab.declared.only, agg="ATT"))  %>% 
  filter(term == "ATT")

df.atts.dec <- rbind(df.atts.dec, data.frame(list(
  dec.estimate = tidysun$estimate,
  dec.std.error = tidysun$std.error,
  method = "sunab"
)))

forest_change_attgt_dec <- att_wrapper(forest_change_declared_only, "Forest_change_perc", "Year_dec")

forest_change_attgt_overall_dec <- aggte(forest_change_attgt_dec, type = "simple", na.rm=TRUE)
forest_change_attgt_dynamic_dec <- aggte(forest_change_attgt_dec, type = "dynamic", na.rm=TRUE, min_e=-10, max_e=10)

plot_dynamic_att(forest_change_attgt_dynamic_dec) +
    xlab("Number of years to land tenure") +
    ylab("Percent forest change") +
  xlim(-10, 10)  +
  theme_bw()

df.atts.dec <- rbind(df.atts.dec, data.frame(
  dec.estimate=forest_change_attgt_overall_dec$overall.att,
  dec.std.error=forest_change_attgt_overall_dec$overall.se,
  method="CS"
))

df.atts.dec

overall.results <- merge(df.atts, df.atts.dec, by="method")
write_csv(overall.results, "./Writing/results.csv")

# EVENT STUDY

eventstudy_df <- forest_change %>%
  mutate(
    Year_since_hom = Year - Year_hom,
    pre = Year_since_hom < 0
  ) %>%
  filter(Year_since_hom >= -10 & Year_since_hom <= 10 & Year_since_hom != -1)

#1
event_study_ttest <- forest_change %>%
  mutate(
    Year_since_hom = Year - Year_hom,
    pre = Year_since_hom < 0
  ) %>%
  filter(Year_since_hom >= -10 & Year_since_hom <= 10 & Year_since_hom != -1) %>%
  lm(Forest_change_perc ~ pre, data=.) %>%
  coeftest(vcov. = vcovHC(., type = 'HC1')) %>%
  tidy %>%
  filter(term == "preTRUE")

#2
event_study_ttest <- forest_change_ri %>%
  mutate(
    Year_since_hom = Year - Year_hom,
    pre = Year_since_hom < 0
  ) %>%
  filter(Year_since_hom >= -10 & Year_since_hom <= 10 & Year_since_hom != -1) %>%
  lm(Forest_change_perc ~ pre, data=.) %>%
  coeftest(vcov. = vcovHC(., type = 'HC1')) %>%
  tidy %>%
  filter(term == "preTRUE")

#3
event_study_ttest <- forest_change %>%
  mutate(
    Year_since_hom = Year - Year_hom,
    pre = Year_since_hom < 0
  ) %>%
  filter(Year_since_hom >= -10 & Year_since_hom <= 10 & Year_since_hom != -1) %>%
  lm(Forest_change ~ pre, data=.) %>%
  coeftest(vcov. = vcovHC(., type = 'HC1')) %>%
  tidy %>%
  filter(term == "preTRUE")

#4
forest_change_no_outliers <- forest_change %>% 
  filter(Codigo != "65101" &
           Codigo != "9301" &
           Codigo != "27501" )

event_study_ttest <- forest_change_no_outliers %>%
  mutate(
    Year_since_hom = Year - Year_hom,
    pre = Year_since_hom < 0
  ) %>%
  filter(Year_since_hom >= -10 & Year_since_hom <= 10 & Year_since_hom != -1) %>%
  lm(Forest_change_perc ~ pre, data=.) %>%
  coeftest(vcov. = vcovHC(., type = 'HC1')) %>%
  tidy %>%
  filter(term == "preTRUE")

# regional analysis
# filter by state
nordeste <- forest_change %>% 
  filter(State == "BA" | State == "AL" | State == "PB" | State == "PB" )
sul <- forest_change %>% 
  filter(State == "SC" | State == "PR" | State == "PR/SC" | State == "RS")
centro_oeste <- forest_change %>% filter(State == "MS")
sudeste <- forest_change %>% filter(State == "SP" | State == "RJ" | State == "MG" | State == "ES")

event_study_ttest <- nordeste %>%
  mutate(
    Year_since_hom = Year - Year_hom,
    pre = Year_since_hom < 0
  ) %>%
  filter(Year_since_hom >= -10 & Year_since_hom <= 10 & Year_since_hom != -1) %>%
  lm(Forest_change_perc ~ pre, data=.) %>%
  coeftest(vcov. = vcovHC(., type = 'HC1')) %>%
  tidy %>%
  filter(term == "preTRUE")

event_study_ttest <- sul %>%
  mutate(
    Year_since_hom = Year - Year_hom,
    pre = Year_since_hom < 0
  ) %>%
  filter(Year_since_hom >= -10 & Year_since_hom <= 10 & Year_since_hom != -1) %>%
  lm(Forest_change_perc ~ pre, data=.) %>%
  coeftest(vcov. = vcovHC(., type = 'HC1')) %>%
  tidy %>%
  filter(term == "preTRUE")

event_study_ttest <- centro_oeste %>%
  mutate(
    Year_since_hom = Year - Year_hom,
    pre = Year_since_hom < 0
  ) %>%
  filter(Year_since_hom >= -10 & Year_since_hom <= 10 & Year_since_hom != -1) %>%
  lm(Forest_change_perc ~ pre, data=.) %>%
  coeftest(vcov. = vcovHC(., type = 'HC1')) %>%
  tidy %>%
  filter(term == "preTRUE")

event_study_ttest <- sudeste %>%
  mutate(
    Year_since_hom = Year - Year_hom,
    pre = Year_since_hom < 0
  ) %>%
  filter(Year_since_hom >= -10 & Year_since_hom <= 10 & Year_since_hom != -1) %>%
  lm(Forest_change_perc ~ pre, data=.) %>%
  coeftest(vcov. = vcovHC(., type = 'HC1')) %>%
  tidy %>%
  filter(term == "preTRUE")

# BANDWIDTHS
event_study_ttest <- forest_change %>%
  mutate(
    Year_since_hom = Year - Year_hom,
    pre = Year_since_hom < 0
  ) %>%
  filter(Year_since_hom >= -20 & Year_since_hom <= 20 & Year_since_hom != -1) %>%
  lm(Forest_change_perc ~ pre, data=.) %>%
  coeftest(vcov. = vcovHC(., type = 'HC1')) %>%
  tidy %>%
  filter(term == "preTRUE")
event_study_ttest

event_study_ttest <- forest_change %>%
  mutate(
    Year_since_hom = Year - Year_hom,
    pre = Year_since_hom < 0
  ) %>%
  filter(Year_since_hom >= -15 & Year_since_hom <= 15 & Year_since_hom != -1) %>%
  lm(Forest_change_perc ~ pre, data=.) %>%
  coeftest(vcov. = vcovHC(., type = 'HC1')) %>%
  tidy %>%
  filter(term == "preTRUE")
event_study_ttest

event_study_ttest <- forest_change %>%
  mutate(
    Year_since_hom = Year - Year_hom,
    pre = Year_since_hom < 0
  ) %>%
  filter(Year_since_hom >= -5 & Year_since_hom <= 5 & Year_since_hom != -1) %>%
  lm(Forest_change_perc ~ pre, data=.) %>%
  coeftest(vcov. = vcovHC(., type = 'HC1')) %>%
  tidy %>%
  filter(term == "preTRUE")
event_study_ttest

# CUTPOINTS
event_study_ttest_test <- forest_change %>%
  mutate(
    Year_since_hom = Year - (Year_hom+1),
    pre = (Year_since_hom+1) < 0
  ) %>%
  filter(Year_since_hom >= -10 & Year_since_hom <= 10 & Year_since_hom != -1) %>%
  lm(Forest_change_perc ~ pre, data=.) %>%
  coeftest(vcov. = vcovHC(., type = 'HC1')) %>%
  tidy %>%
  filter(term == "preTRUE")
event_study_ttest_test

event_study_ttest_test <- forest_change %>%
  mutate(
    Year_since_hom = Year - (Year_hom-1),
    pre = (Year_since_hom-1) < 0
  ) %>%
  filter(Year_since_hom >= -10 & Year_since_hom <= 10 & Year_since_hom != -1) %>%
  lm(Forest_change_perc ~ pre, data=.) %>%
  coeftest(vcov. = vcovHC(., type = 'HC1')) %>%
  tidy %>%
  filter(term == "preTRUE")
event_study_ttest_test

event_study_ttest_test <- forest_change %>%
  mutate(
    Year_since_hom = Year - (Year_hom+5),
    pre = (Year_since_hom+5) < 0
  ) %>%
  filter(Year_since_hom >= -10 & Year_since_hom <= 10 & Year_since_hom != -1) %>%
  lm(Forest_change_perc ~ pre, data=.) %>%
  coeftest(vcov. = vcovHC(., type = 'HC1')) %>%
  tidy %>%
  filter(term == "preTRUE")
event_study_ttest_test

event_study_ttest_test <- forest_change %>%
  mutate(
    Year_since_hom = Year - (Year_hom-5),
    pre = (Year_since_hom-5) < 0
  ) %>%
  filter(Year_since_hom >= -10 & Year_since_hom <= 10 & Year_since_hom != -1) %>%
  lm(Forest_change_perc ~ pre, data=.) %>%
  coeftest(vcov. = vcovHC(., type = 'HC1')) %>%
  tidy %>%
  filter(term == "preTRUE")
event_study_ttest_test

event_study_ttest %>%
  mutate(
    ci_lower = estimate - 1.96*std.error,
    ci_upper = estimate + 1.96*std.error,
  )

# ROBUSTNESS TESTS======
# Going to do every combination of:
# (TI vs. TI/RI) x (With outliers vs. without) x (forest_change vs forest_change_perc)
# for both tenure and declarada dates.
# Note: 
#  - none of the outliers appear in the declared-only data
#  - there are no RI territories in the declared-only data

forest_change_ri <- read.csv("Data/Excel_spreadsheets/forest_change_FINAL.csv") %>% 
  select(-X, -X.1)

outlier.codigos <- c(65101, 9301, 27501)
forest_change_ri_no_outliers <- 
  forest_change_ri %>%
  filter(!(Codigo %in% outlier.codigos))

forest_change_no_outliers_attgt <- att_wrapper(forest_change_no_outliers, "Forest_change_perc", "Year_hom_0")
forest_change_no_outliers_attgt_overall <- aggte(forest_change_no_outliers_attgt, type = "simple", na.rm=TRUE)
forest_change_ri_attgt <- att_wrapper(forest_change_ri, "Forest_change_perc", "Year_hom_0")
forest_change_ri_attgt_overall <- aggte(forest_change_ri_attgt, type = "simple", na.rm=TRUE)
forest_change_ri_no_outliers_attgt <- att_wrapper(forest_change_ri_no_outliers, "Forest_change_perc", "Year_hom_0")
forest_change_ri_no_outliers_attgt_overall <- aggte(forest_change_ri_no_outliers_attgt, type = "simple", na.rm=TRUE)

nonperc_forest_change_attgt <- att_wrapper(forest_change, "Log_forest_change", "Year_hom_0")
nonperc_forest_change_attgt_overall <- aggte(nonperc_forest_change_attgt, type = "simple", na.rm=TRUE)
nonperc_forest_change_attgt_dec <- att_wrapper(forest_change_declared_only, "Log_forest_change", "Year_dec")
nonperc_forest_change_attgt_overall_dec <- aggte(nonperc_forest_change_attgt_dec, type = "simple", na.rm=TRUE)
nonperc_forest_change_no_outliers_attgt <- att_wrapper(forest_change_no_outliers, "Log_forest_change", "Year_hom_0")
nonperc_forest_change_no_outliers_attgt_overall <- aggte(nonperc_forest_change_no_outliers_attgt, type = "simple", na.rm=TRUE)
nonperc_forest_change_ri_attgt <- att_wrapper(forest_change_ri, "Log_forest_change", "Year_hom_0")
nonperc_forest_change_ri_attgt_overall <- aggte(nonperc_forest_change_ri_attgt, type = "simple", na.rm=TRUE)
nonperc_forest_change_ri_no_outliers_attgt <- att_wrapper(forest_change_ri_no_outliers, "Log_forest_change", "Year_hom_0")
nonperc_forest_change_ri_no_outliers_attgt_overall <- aggte(nonperc_forest_change_ri_no_outliers_attgt, type = "simple", na.rm=TRUE)

summary(aggte(att_wrapper(sul, "Forest_change_perc", "Year_hom_0"), type="simple", na.rm=TRUE))
summary(aggte(att_wrapper(nordeste, "Forest_change_perc", "Year_hom_0"), type="simple", na.rm=TRUE))
summary(aggte(att_wrapper(centro_oeste, "Forest_change_perc", "Year_hom_0"), type="simple", na.rm=TRUE))
summary(aggte(att_wrapper(sudeste, "Forest_change_perc", "Year_hom_0"), type="simple", na.rm=TRUE))
summary(aggte(att_wrapper(rbind(nordeste, centro_oeste, sudeste), "Forest_change_perc", "Year_hom_0"), type="simple", na.rm=TRUE))

robustness.table <- rbind(
  data.frame(  with.ri = FALSE, drop.outliers = FALSE, dv = "Forest change percent", stage="hom", ATT = forest_change_attgt_overall$overall.att, SE = forest_change_attgt_overall$overall.se ),
  data.frame(  with.ri = TRUE, drop.outliers = FALSE, dv = "Forest change percent", stage="hom", ATT = forest_change_ri_attgt_overall$overall.att, SE = forest_change_ri_attgt_overall$overall.se ),
  data.frame( with.ri = FALSE, drop.outliers = TRUE, dv = "Forest change percent", stage="hom", ATT = forest_change_no_outliers_attgt_overall$overall.att, SE = forest_change_no_outliers_attgt_overall$overall.se ),
  data.frame(  with.ri = TRUE, drop.outliers = TRUE, dv = "Forest change percent", stage="hom", ATT = forest_change_ri_no_outliers_attgt_overall$overall.att, SE = forest_change_ri_no_outliers_attgt_overall$overall.se ),
  
  data.frame(  with.ri = FALSE, drop.outliers = FALSE, dv = "Forest change", stage="hom", ATT = nonperc_forest_change_attgt_overall$overall.att, SE = nonperc_forest_change_attgt_overall$overall.se ),
  data.frame( with.ri = FALSE, drop.outliers = TRUE, dv = "Forest change", stage="hom", ATT = nonperc_forest_change_no_outliers_attgt_overall$overall.att, SE = nonperc_forest_change_no_outliers_attgt_overall$overall.se ),
  data.frame(  with.ri = TRUE, drop.outliers = FALSE, dv = "Forest change", stage="hom", ATT = nonperc_forest_change_ri_attgt_overall$overall.att, SE = nonperc_forest_change_ri_attgt_overall$overall.se ),
  data.frame(  with.ri = TRUE, drop.outliers = TRUE, dv = "Forest change", stage="hom", ATT = nonperc_forest_change_ri_no_outliers_attgt_overall$overall.att, SE = nonperc_forest_change_ri_no_outliers_attgt_overall$overall.se ),
  
  data.frame( with.ri = FALSE, drop.outliers = FALSE, dv = "Forest change percent", stage="dec", ATT=forest_change_attgt_overall_dec$overall.att, SE=forest_change_attgt_overall_dec$overall.se ),
  data.frame( with.ri = FALSE, drop.outliers = FALSE, dv = "Forest change", stage="dec", ATT=nonperc_forest_change_attgt_overall_dec$overall.att, SE=nonperc_forest_change_attgt_overall_dec$overall.se )
)

robustness.table.display <- robustness.table %>%
  mutate(
    pval = pnorm(ATT / SE, lower.tail=FALSE),
    statstar = ifelse(pval < 0.0005, "***", ifelse(pval < 0.005, "**", ifelse(pval < 0.025, "*", ifelse(pval < 0.05, "\\dagger", "")))),
    ATT.display = paste("\\makecell{", format(round(1000*ATT)/1000), " ", statstar, "\\\\(", round(1000*SE)/1000, ")}", sep=""),
    with.ri = ifelse(with.ri, "Y", "N"),
    drop.outliers = ifelse(drop.outliers, "Y", "N")
  ) %>%
  select(with.ri, drop.outliers, dv, stage, ATT.display) %>%
  dplyr::rename(`With RI` = with.ri,
         `Drop outliers` = drop.outliers,
         `Dependent var.` = dv,
         `Stage` = stage,
         `ATT` = ATT.display)

kbl(robustness.table.display, booktabs=T, format="latex", escape=FALSE)

