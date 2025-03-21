## TCC MBA USP ESALQ ##

# Autor: Gabriel P. Esteves

# load packages
library(tidyverse); library(here); library(survey); library(broom); 
library(writexl); library(dagitty); library(patchwork); library(glue); 
library(gtsummary); library(marginaleffects); library(ggeffects);
library(splines); library(colorspace); library(gridExtra)

# no scientific notation
options(scipen=999)


#custom ggplot theme
theme_avp <- function() {
  theme_bw(base_size=10) +
    theme(axis.title.x = element_text(size = 9),
          axis.title.y = element_text(size = 9),
          panel.grid.minor = element_blank())
}

# function for assessing svyGLM assumptions

svy_check_model <- function(model) {
  dat <- tibble(fitted = fitted(model),
                resid = resid(model),
                std_resid = svydiags::svystdres(model, stvar="SDMVSTRA", clvar="SDMVPSU")[["stdresids"]])
  
  linearity_p_std <- ggplot(dat, aes(x=fitted, y=std_resid)) +
    geom_point(alpha=.05) +
    geom_hline(yintercept=0,color="red",linetype="dashed") +
    theme_minimal() +
    labs(title = "Fitted values e resíduos padronizados",
         x="Fitted values",
         y="Resíduos padronizados")
  
  normality_p <- ggplot(dat, aes(x = resid)) +
    geom_histogram(color = "black", fill = "lightblue") +
    theme_minimal() +
    labs(title = "Histograma de resíduos",
         x = "Resíduos",
         y = "Frequência")
  
  qq_p_std <- ggplot(dat, aes(sample = std_resid)) +
    stat_qq(alpha=.2) +
    stat_qq_line(color = "red") +
    theme_minimal() +
    labs(title = "Q-Q plot de resíduos padronizados",
         x = "Quantis teóricos",
         y = "Quantis amostrais")
  
  vif_table <- car::vif(model) |> as.data.frame() |> 
    rownames_to_column("Variável") |> select('Variável', 'GVIF') |> 
    rename(VIF = GVIF) |> 
    mutate('Variável' = c("PTN", "Atividade Fisica", "CHO", "LIP", "Gênero",
                        "Peso corporal", "Altura", "Idade", "Uso glicocorticoide",
                        "Interação PTN:Atv Fisica"),
           VIF = round(VIF, 2)) |> as.data.frame() |> tableGrob(rows = NULL)
  
  p_full <- (linearity_p_std + normality_p) / (qq_p_std + vif_table)
  
  return(p_full)
}

# function for df for predictions

get_pred_df <- function(design, predictor) {
  if (predictor == "PTN") {
    newdata <- data.frame(   
      PTN = design$variables$PTN,
      total_vig_and_mod_min = mean(design$variables$total_vig_and_mod_min),
      CHO = mean(design$variables$CHO),
      LIP = mean(design$variables$LIP),
      RIAGENDR = "1", 
      BMXWT = mean(design$variables$BMXWT),
      BMXHT = mean(design$variables$BMXHT),
      RIDAGEYR = mean(design$variables$RIDAGEYR),
      days_gc = mean(design$variables$days_gc))
  } else {
    newdata <- data.frame(   
      PTN = mean(design$variables$PTN),
      total_vig_and_mod_min = design$variables$total_vig_and_mod_min,
      CHO = mean(design$variables$CHO),
      LIP = mean(design$variables$LIP),
      RIAGENDR = "1", 
      BMXWT = mean(design$variables$BMXWT),
      BMXHT = mean(design$variables$BMXHT),
      RIDAGEYR = mean(design$variables$RIDAGEYR),
      days_gc = mean(design$variables$days_gc))
  }
  
  return(newdata)
}

# load raw dataset
load("data/nhanes_raw_dataframe.Rdata")

# creating and transforming variables before setting up survey design
# nutrients

rh_dxa_all <- nhanes_raw_dataframe |> 
  mutate(kcal = case_when(
    DRDINT.x == 2 ~ ((DR1TKCAL+DR2TKCAL)/2), # these lines calculate the mean between two dietary recalls, when two are available
    DRDINT.x == 1 ~ DR1TKCAL # or just pick the first one that is available, if only one is.
  )) |> 
  mutate(PTN = case_when(
    DRDINT.x == 2 ~ ((DR1TPROT+DR2TPROT)/2), # same for other nutrients.
    DRDINT.x == 1 ~ DR1TPROT
  )) |> 
  mutate(CHO = case_when(
    DRDINT.x == 2 ~ ((DR1TCARB+DR2TCARB)/2),
    DRDINT.x == 1 ~ DR1TCARB
  )) |> 
  mutate(LIP = case_when(
    DRDINT.x == 2 ~ ((DR1TTFAT+DR2TTFAT)/2),
    DRDINT.x == 1 ~ DR1TTFAT
  ))

# other relevant variables

rh_dxa_all <- rh_dxa_all |> 
  mutate(PTNgkg = PTN/BMXWT) |> 
  mutate(femur_popmean_male = 1.04, # Reference values from NHANES III, Looker AC et al 1997.
         femur_popsd_male = 0.144,
         femurneck_popmean_male = 0.93,
         femurneck_popsd_male = 0.137,
         femur_popmean_female = 0.94,  
         femur_popsd_female = 0.122, 
         femurneck_popmean_female = 0.86,
         femurneck_popsd_female = 0.12,
         femur_t_score = case_when(
           RIAGENDR == "1" ~ ((DXXOFBMD-femur_popmean_male)/femur_popsd_male),
           RIAGENDR == "2" ~ ((DXXOFBMD-femur_popmean_female)/femur_popsd_female))
         ) |> 
  mutate(femurneck_t_score = case_when(
    RIAGENDR == "1" ~ ((DXXNKBMD-femurneck_popmean_male)/femurneck_popsd_male),
    RIAGENDR == "2" ~ ((DXXNKBMD-femurneck_popmean_female)/femurneck_popsd_female))) |> 
  mutate(femur_osteoporosis = as.factor(case_when(
    femur_t_score >= -2.5 ~ "No",
    femur_t_score < -2.5 ~ "Yes"
  )),
  femurneck_osteoporosis = as.factor(case_when(
    femurneck_t_score >= -2.5 ~ "No",
    femurneck_t_score < -2.5 ~ "Yes"
  )),
  ALM = DXDLALE + DXDRALE + DXDLLLE + DXDRLLE,
  ALM_BMI = (ALM/1000)/BMXBMI,
  sarc_class = as.factor(case_when(RIAGENDR == "1" & ALM_BMI <.789 ~ '1',
                                   RIAGENDR == "2" & ALM_BMI <.512 ~ '1',
                                   RIAGENDR == "1" & ALM_BMI >=.789 ~ '0',
                                   RIAGENDR == "2" & ALM_BMI >=.512 ~ '0'))) 

# creating adequate sample weights
# herein we use the 24h recall day 1 weight, and divide by the number of cycles (6)

rh_dxa_all <- rh_dxa_all |> 
  mutate(weights_all = case_when(
    DRDINT.x == 2 ~ (WTDR2D.x * (1/7)),
    DRDINT.x == 1 ~ (WTDRD1.x * (1/7))), 
    weights_all_exam = WTMEC2YR * (1/7))

# calculating physical activity
# obtain METs for each domain and intensity, then sum it up, then classify

rh_dxa_all <- rh_dxa_all |> 
  # cleaning variables first
  mutate(PAQ605 = case_when(PAQ605 %in% c(7, 9) ~ NA_real_, # these are codes for non-responses within the questionnaires, important to remove them
                            TRUE ~ PAQ605),
         PAQ610 = case_when(PAQ610 %in% c(77, 99) ~ NA_real_,
                            TRUE ~ PAQ610),
         PAD615 = case_when(PAD615 %in% c(7777, 9999) ~ NA_real_,
                            TRUE ~ PAD615),
         
         PAQ620 = case_when(PAQ620 %in% c(7, 9) ~ NA_real_,
                            TRUE ~ PAQ620),
         PAQ625 = case_when(PAQ625 %in% c(77, 99) ~ NA_real_,
                            TRUE ~ PAQ625),
         PAD630 = case_when(PAD630 %in% c(7777, 9999) ~ NA_real_,
                            TRUE ~ PAD630),
         
         PAQ635 = case_when(PAQ635 %in% c(7, 9) ~ NA_real_,
                            TRUE ~ PAQ635),
         PAQ640 = case_when(PAQ640 %in% c(77, 99) ~ NA_real_,
                            TRUE ~ PAQ640),
         PAD645 = case_when(PAD645 %in% c(7777, 9999) ~ NA_real_,
                            TRUE ~ PAD645),
         
         PAQ650 = case_when(PAQ650 %in% c(7, 9) ~ NA_real_,
                            TRUE ~ PAQ650),
         PAQ655 = case_when(PAQ655 %in% c(77, 99) ~ NA_real_,
                            TRUE ~ PAQ655),
         PAD660 = case_when(PAD660 %in% c(7777, 9999) ~ NA_real_,
                            TRUE ~ PAD660),
         
         PAQ665 = case_when(PAQ665 %in% c(7, 9) ~ NA_real_,
                            TRUE ~ PAQ665),
         PAQ670 = case_when(PAQ670 %in% c(77, 99) ~ NA_real_,
                            TRUE ~ PAQ670),
         PAD675 = case_when(PAD675 %in% c(7777, 9999) ~ NA_real_,
                            TRUE ~ PAD675),
         # then calculate METs per domain and overall
         work_vig_min = if_else(PAQ605 == 2, 0, PAQ610 * PAD615), # calculating MET-minutes based on GPAQ questionnaire scoring
         work_mod_min = if_else(PAQ620 == 2, 0, PAQ625 * PAD630),
         work_vig_and_mod_min = work_vig_min + work_mod_min,
         transport_MET = if_else(PAQ635 == 2, 0, PAQ640 * PAD645),
         rec_vig_min = if_else(PAQ650 == 2, 0, PAQ655 * PAD660),
         rec_mod_min = if_else(PAQ665 == 2, 0, PAQ670 * PAD675),
         rec_vig_and_mod_min = rec_vig_min + rec_mod_min,
         total_vig_min = work_vig_min + work_mod_min,
         total_vig_and_mod_min = work_vig_and_mod_min + rec_vig_and_mod_min,
         work_MET = (work_vig_min * 8) + (work_mod_min * 4),
         rec_MET = (rec_vig_min * 8) + (rec_mod_min * 4),
         MET = work_MET + rec_MET + transport_MET, # summing it up
         # classification based on METs
         MET_class = case_when(MET < 600 ~ "Very low PA",
                               MET >= 600 & MET < 1200 ~ "Moderate PA",
                               MET >= 1200 & MET < 1800 ~ "High PA",
                               MET >= 1800 ~ "Very high PA")
  )

# medication use, needs to be processed in separate due to different dataframe structure
# calculating gc use and days of gc use

# load previously joined medication dataframe

load("data/medication_df.Rdata")

# create a string with the names of glucocorticoids 
gluc_names <- c('^(DEFLAZACORT|PREDNISOLONE|DEXAMETHASONE|HYDROCORTISONE|METHYLPREDNISOLONE|PREDNISONE|CORTISONE)$')

med_df <- med_df |> 
  mutate(uses_gc = if_else(stringr::str_detect(RXDDRUG, gluc_names), 1, 0), # detect these names
         days_gc = case_when(uses_gc == 0 ~ 0,
                             RXDDAYS %in% c(77777, 99999) ~ 77777,
                             TRUE ~ RXDDAYS * uses_gc)) |> 
  select(SEQN, uses_gc, days_gc, RXDDRUG) |> 
  filter(stringr::str_detect(RXDDRUG, gluc_names)) |> 
  group_by(SEQN) |> # sometimes some participants would report using more than one GC type,
  summarise(days_gc = mean(days_gc)) |> # so we group by participant and take the mean of both medications to have one single number
  ungroup()

# merge into main and fix

rh_dxa_all <- rh_dxa_all |> full_join(med_df) |> 
  mutate(days_gc = if_else(is.na(days_gc), 0, days_gc), # if days_gc is NA, it means the person does not take GC, i.e. a vlaue of 0
         days_gc = if_else(days_gc == 77777, NA_real_, days_gc)) # non response code

# setting up survey design

rh_dxa_all_analyse <- subset(rh_dxa_all, !is.na(SDMVPSU) & !is.na(SDMVSTRA)) |> 
  dplyr::select(SEQN, PTN, CHO, LIP, kcal, PTNgkg, DXDTOLE, DXDTOBMD, DXXOFBMD, 
                DXXOSBMD, RIAGENDR, BMXHT, BMXWT, ALM, ALM_BMI, sarc_class, RIDAGEYR, 
                SDMVPSU, SDMVSTRA, weights_all, RIDRETH1, BMXBMI, femur_t_score, femurneck_t_score, 
                femurneck_osteoporosis, DR1DRSTZ, DXAEXSTS, DXAFMRST, DXASPNST, 
                DXXNKBMD, MET, MET_class,
                work_vig_min, rec_vig_min, total_vig_min,
                days_gc, cycle,
                work_vig_and_mod_min, rec_vig_and_mod_min,
                total_vig_and_mod_min,
                DXDLALE, DXDRALE, DXDLLLE, DXDRLLE, ALM, ALM_BMI)

# set survey
# set individuals without sample weights to 0, to exclude later
rh_dxa_all_analyse <- rh_dxa_all_analyse |> 
  mutate(weights_all = if_else(is.na(weights_all), 0, weights_all)) |> 
  mutate(RIAGENDR = as.factor(RIAGENDR),
         cycle = as.factor(cycle)) # set categorical variables to factor

NHANES_all <- svydesign(data=rh_dxa_all_analyse, 
                        id=~SDMVPSU, 
                        strata=~SDMVSTRA, 
                        weights=~weights_all,
                        nest=TRUE)

nrow(NHANES_all$variables)

NHANES_rh <- subset(NHANES_all, weights_all != 0) # excluding individuals without sample weight

nrow(NHANES_rh$variables)

NHANES_rh <- subset(NHANES_rh, DR1DRSTZ == "1") # selecting adequate r24 data

nrow(NHANES_rh$variables)

NHANES_rh <- subset(NHANES_rh, RIDAGEYR >= 18 & RIDAGEYR <= 86) # selecting adults 

nrow(NHANES_rh$variables) # adults older than 85 are still coded as 85 in NHANES.

NHANES_rh <- subset(NHANES_rh, !is.na(total_vig_and_mod_min)) # sem atv fisica 

nrow(NHANES_rh$variables) 

NHANES_subset_lm <- subset(NHANES_rh, DXAEXSTS == "1" & !is.na(DXDTOLE))

nrow(NHANES_subset_lm$variables)

NHANES_subset_femur <- subset(NHANES_rh, DXAFMRST == "1" & !is.na(DXXOFBMD))

nrow(NHANES_subset_femur$variables)

# creating subsets per outcome

NHANES_subset_lm <- subset(NHANES_subset_lm, DXAEXSTS == "1" & !is.na(DXDTOLE) # make sure subsets
                           & !is.na(PTN) & !is.na(CHO) & !is.na(LIP) # have available data on all covariates
                           & !is.na(RIAGENDR) & !is.na(BMXWT) & !is.na(RIDAGEYR)
                           & !is.na(total_vig_and_mod_min) & !is.na(days_gc) & !is.na(BMXHT))


NHANES_subset_sarc <- subset(NHANES_subset_lm, DXAEXSTS == "1" & !is.na(DXDTOLE) # make sure subsets
                           & !is.na(PTN) & !is.na(CHO) & !is.na(LIP) # have available data on all covariates
                           & !is.na(RIAGENDR) & !is.na(BMXWT) & !is.na(RIDAGEYR)
                           & !is.na(total_vig_and_mod_min) & !is.na(days_gc) & !is.na(BMXHT) & !is.na(sarc_class))

nrow(NHANES_subset_lm$variables) # subset sample size
nrow(NHANES_subset_sarc$variables) # subset sample size

NHANES_subset_femur <- subset(NHANES_subset_femur, DXAFMRST == "1" & !is.na(DXXOFBMD) & !is.na(DXXNKBMD) 
                              & !is.na(PTN) & !is.na(CHO) & !is.na(LIP)
                              & !is.na(RIAGENDR) & !is.na(BMXWT) & !is.na(RIDAGEYR)
                              & !is.na(total_vig_and_mod_min) & !is.na(days_gc) & !is.na(BMXHT) & !is.na(femurneck_osteoporosis))

nrow(NHANES_subset_femur$variables)

# separate dfs for each outcome, for table 1

lm_dat <- as_tibble(NHANES_subset_lm$variables)

fem_dat <- as_tibble(NHANES_subset_femur$variables)

sample_dat <- reduce(list(lm_dat, fem_dat), bind_rows) |> distinct()

nrow(sample_dat) # sample size

## Visualizando distribuições

sample_dat |> select(BMXWT, BMXHT, BMXBMI, days_gc, total_vig_and_mod_min, kcal, 
                     CHO, PTN, LIP, DXDTOLE, ALM, ALM_BMI, DXXOFBMD, femurneck_t_score) |> 
  pivot_longer(everything()) |> 
  ggplot(aes(x=value)) +
  facet_wrap(~name, scales="free") +
  geom_density()

## Modelling ###################################################################

# Modelos Lean mass e Atividade física vigorosa e moderada total (trabalho + lazer)

# simple models

lm_ptn_total_simple <- svyglm(DXDTOLE ~ PTN + total_vig_and_mod_min + 
                                CHO + LIP + RIAGENDR + BMXWT + BMXHT +
                                 + RIDAGEYR + days_gc,
                              design=NHANES_subset_lm) 

car::Anova(lm_ptn_total_simple, "II", "F")

plot(predict_response(lm_ptn_total_simple, "PTN [all]"), 
     show_residuals=TRUE, show_residuals_line = TRUE,
     dot_size=.5) + # Usando partial regressions plots para visualizar fits inadequados
  xlim(c(0,350)) +
  labs(x="Ingestão de proteína", y="Massa magra total (Kg)", 
       title="Valores preditos de proteína vs. resíduos parciais") +
  theme_avp()

lm_ptnPoly2_total_simple <- svyglm(DXDTOLE ~ poly(PTN, 2) + total_vig_and_mod_min + CHO + LIP + 
                                   RIAGENDR + BMXWT + BMXHT + RIDAGEYR + days_gc,
                                   design=NHANES_subset_lm)

lm_ptnPoly3_total_simple <- svyglm(DXDTOLE ~ poly(PTN, 3) + total_vig_and_mod_min + CHO + LIP + RIAGENDR + 
                                     BMXWT + BMXHT
                                    + RIDAGEYR + days_gc,
                                    design=NHANES_subset_lm) 

lm_ptnns3_total_simple <- svyglm(DXDTOLE ~ splines::ns(PTN, df=3) + total_vig_and_mod_min + 
                                   CHO + LIP + RIAGENDR + BMXWT + BMXHT
                                  + RIDAGEYR + days_gc,
                                  design=NHANES_subset_lm) 

AIC(lm_ptn_total_simple, lm_ptnPoly2_total_simple,
    lm_ptnPoly3_total_simple, lm_ptnns3_total_simple) # Modelo linear simples com menor AIC

anova(lm_ptn_total_simple, lm_ptnPoly2_total_simple) 
anova(lm_ptn_total_simple, lm_ptnPoly3_total_simple) 
anova(lm_ptn_total_simple, lm_ptnns3_total_simple) # Ausencia de diferença com qualquer outro fit

car::Anova(lm_ptn_total_simple, "II", "F")

## Mantemos o modelo linear

# Agora olhando para a atv física

plot(predict_response(lm_ptn_total_simple, "total_vig_and_mod_min [all]"), 
     show_residuals=TRUE, show_residuals_line = TRUE,
     dot_size=.5) + # use partial plots to check the fit of protein 
  xlim(c(0, 1000)) +
  labs(x="Minutos de atividade física\nmoderada a vigorosa", y="Massa magra total (Kg)", 
       title="Valores preditos de proteína vs. resíduos parciais") +
  theme_avp()

lm_ptn_total_poly2_simple <- svyglm(DXDTOLE ~ PTN + poly(total_vig_and_mod_min, 2) + CHO + LIP + RIAGENDR + 
                                      BMXWT + BMXHT 
                              + RIDAGEYR + days_gc,
                              design=NHANES_subset_lm) 

lm_ptn_total_poly3_simple <- svyglm(DXDTOLE ~ PTN + poly(total_vig_and_mod_min, 3) + CHO + LIP + 
                                      RIAGENDR + BMXWT + BMXHT
                                    + RIDAGEYR + days_gc,
                                    design=NHANES_subset_lm) 

lm_ptn_total_ns3_simple <- svyglm(DXDTOLE ~ PTN + splines::ns(total_vig_and_mod_min, df=3) + CHO + LIP + 
                                    RIAGENDR + BMXWT + BMXHT
                                    + RIDAGEYR + days_gc,
                                    design=NHANES_subset_lm) 

AIC(lm_ptn_total_simple, lm_ptn_total_poly2_simple,
    lm_ptn_total_poly3_simple, lm_ptn_total_ns3_simple)

## spline natural com 3 graus de liberdade mostrou menores AIC

anova(lm_ptn_total_simple, lm_ptn_total_ns3_simple) ## diferença significativa entre modelos
svy_check_model(lm_ptn_total_ns3_simple) # pressupostos atingidos

plot(predict_response(lm_ptn_total_ns3_simple, "total_vig_and_mod_min [all]"), 
     show_residuals=TRUE, show_residuals_line = TRUE,
     dot_size=.5) +
  labs(x="Minutos de atividade física\nmoderada a vigorosa", y="Massa magra total (Kg)", 
       title="Valores preditos de proteína vs. resíduos parciais") +
  theme_avp() # boa concordância da linha predita com a linha de resíduos

# checando interação

lm_interaction  <- svyglm(DXDTOLE ~ PTN * 
                            splines::ns(total_vig_and_mod_min, df=3) + CHO + LIP + 
                            RIAGENDR + BMXWT + BMXHT +
                            + RIDAGEYR + days_gc,
                          design=NHANES_subset_lm) 

lm_interaction_centered  <- svyglm(DXDTOLE ~ scale(PTN, center=TRUE, scale=FALSE) * 
                            scale(splines::ns(total_vig_and_mod_min, df=3), center=TRUE, scale=FALSE) + CHO + LIP + 
                                          RIAGENDR + BMXWT + BMXHT +
                                        + RIDAGEYR + days_gc,
                                        design=NHANES_subset_lm) 

car::Anova(lm_interaction, type="II", test="F")

summary(lm_interaction)

AIC(lm_interaction, lm_ptn_total_ns3_simple)

anova(lm_ptn_total_ns3_simple, 
      lm_interaction) # O modelo com interação abaixa a AIC ainda mais 
                                    # e apresenta diferença significativa no F test.

svy_check_model(lm_interaction_centered)

ggsave("figuras/lm_finalmodel_assump.png", dpi=600, unit="in", height=8, width=8)

# gráficos de predição

pred_1 <- ggpredict(lm_interaction1, 
                    terms=c("PTN[1:500]", "total_vig_and_mod_min[0:500]")) |> as_tibble() |> # predições nos valores específicos 
          rename(PTN = x, atv_fisica = group)

pred_1_selected <- pred_1 |> filter(atv_fisica %in% c(0, 100, 200, 300, 400, 500)) 

(
lm_pred_p_1 <- ggplot(pred_1_selected, aes(x=PTN, y=predicted, color=atv_fisica, fill=atv_fisica,
                            ymin=conf.low, ymax=conf.high)) +
  facet_wrap(~atv_fisica, labeller = labeller(atv_fisica = ~ paste(.x, "min/sem"))) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.35, color=NA) +
  geom_line() +
  labs(x="Ingestão de proteína (g/dia)", y="Massa magra predita (kg)", 
       color="Atividade física\nmoderada e vigorosa\n(min/sem)",
       fill="Atividade física\nmoderada e vigorosa\n(min/sem)") +
  scale_color_discrete_sequential(labels=paste0(seq(0, 500, 100), " min/sem"), palette = "Reds", l1 = 25, c2 = 150, p1 = 1) + 
  scale_fill_discrete_sequential(labels=paste0(seq(0, 500, 100), " min/sem"), palette = "Reds", l1 = 25, c2 = 150, p1 = 1) + 
  scale_x_continuous(limits=c(25, 250)) +
  scale_y_continuous(labels = scales::label_number(scale = 1e-3),
                     limits=c(51000, 64000), breaks=seq(51000, 64000, 2500)) +
  guides(color='none', fill='none') +
  theme_avp()
)

# gráfico de slope para visualizar interação

pred_1_slopes <- pred_1 |> filter(PTN %in% c(1, 2)) |> 
  select(PTN, predicted, atv_fisica) |> 
  pivot_wider(values_from=predicted, names_from=PTN) |> 
  mutate(slope = `2` - `1`) |> 
  mutate(atv_fisica = as.numeric(atv_fisica))

(
lm_pred_p_2 <- ggplot(pred_1_slopes, aes(x=atv_fisica, y=slope)) +
  geom_line(aes(group=1)) +
  geom_vline(xintercept=230, linetype="dashed", alpha=.25) +
  annotate(geom='label', label='230 min.', x=230, y=19.85, size=2.5) +
  scale_x_continuous(breaks=seq(0, 500, 100)) +
  scale_y_continuous(breaks=c(seq(10, 20, 2.5), 19.1),
                     limits=c(9, 20)) +
  theme_avp() +
  labs(x="Atividade física moderada e vigorosa (min/sem)", y="Slope para ingestão de proteína")
)

(lm_pred_p_1 / lm_pred_p_2) +
  plot_annotation(tag_levels="A") +
  plot_layout(height=c(1, .5))

#ggsave("figuras/plot1.png", dpi=600, unit="in", height=8, width=6)

## análise osso

# simple models

femur_ptn_total_simple <- svyglm(DXXNKBMD ~ PTN + total_vig_and_mod_min + 
                                CHO + LIP + RIAGENDR + BMXWT + BMXHT +
                                + RIDAGEYR + days_gc,
                              design=NHANES_subset_femur) 

summary(femur_ptn_total_simple)

car::Anova(femur_ptn_total_simple, "II", "F")

svy_check_model(femur_ptn_total_simple)

plot(predict_response(femur_ptn_total_simple, "PTN [all]"), 
     show_residuals=TRUE, show_residuals_line = TRUE,
     dot_size=.5) + # Usando partial regressions plots para visualizar fits inadequados
  xlim(c(0,350)) +
  labs(x="Ingestão de proteína", y="DMO colo femural", 
       title="Valores preditos de proteína vs. resíduos parciais") +
  theme_avp()

femur_ptnPoly2_total_simple <- svyglm(DXXNKBMD ~ poly(PTN, 2) + total_vig_and_mod_min + CHO + LIP + 
                                     RIAGENDR + BMXWT + BMXHT + RIDAGEYR + days_gc,
                                   design=NHANES_subset_femur)

femur_ptnPoly3_total_simple <- svyglm(DXXNKBMD ~ poly(PTN, 3) + total_vig_and_mod_min + CHO + LIP + RIAGENDR + 
                                     BMXWT + BMXHT
                                   + RIDAGEYR + days_gc,
                                   design=NHANES_subset_femur) 

femur_ptnns3_total_simple <- svyglm(DXXNKBMD ~ splines::ns(PTN, df=3) + total_vig_and_mod_min + 
                                   CHO + LIP + RIAGENDR + BMXWT + BMXHT
                                 + RIDAGEYR + days_gc,
                                 design=NHANES_subset_femur) 

AIC(femur_ptn_total_simple, femur_ptnPoly2_total_simple,
    femur_ptnPoly3_total_simple, femur_ptnns3_total_simple) # Fit marginalmente melhor no quarto modelo com splines

anova(femur_ptn_total_simple, femur_ptnPoly2_total_simple) 
anova(femur_ptn_total_simple, femur_ptnPoly3_total_simple) 
anova(femur_ptn_total_simple, femur_ptnns3_total_simple) # Ausencia de diferença com qualquer outro fit

## Mantemos o modelo linear

# Agora olhando para a atv física

plot(predict_response(femur_ptn_total_simple, "total_vig_and_mod_min [all]"), 
     show_residuals=TRUE, show_residuals_line = TRUE,
     dot_size=.5) + # use partial plots to check the fit of protein 
  xlim(c(0, 1000)) +
  labs(x="Minutos de atividade física\nmoderada a vigorosa", y="DMO colo femural", 
       title="Valores preditos de proteína vs. resíduos parciais") +
  theme_avp()

femur_ptn_total_poly2_simple <- svyglm(DXXNKBMD ~ PTN + poly(total_vig_and_mod_min, 2) + CHO + LIP + RIAGENDR + 
                                      BMXWT + BMXHT 
                                    + RIDAGEYR + days_gc,
                                    design=NHANES_subset_femur) 

femur_ptn_total_poly3_simple <- svyglm(DXXNKBMD ~ PTN + poly(total_vig_and_mod_min, 3) + CHO + LIP + 
                                      RIAGENDR + BMXWT + BMXHT
                                    + RIDAGEYR + days_gc,
                                    design=NHANES_subset_femur) 

femur_ptn_total_ns3_simple <- svyglm(DXXNKBMD ~ PTN + splines::ns(total_vig_and_mod_min, df=3) + CHO + LIP + 
                                    RIAGENDR + BMXWT + BMXHT
                                  + RIDAGEYR + days_gc,
                                  design=NHANES_subset_femur) 

AIC(femur_ptn_total_simple, femur_ptn_total_poly2_simple,
    femur_ptn_total_poly3_simple, femur_ptn_total_ns3_simple) ## Terceiro modelo com poly 3 mostrou melhor fit

anova(femur_ptn_total_simple, femur_ptn_total_ns3_simple) ## diferença significativa entre modelos
svy_check_model(femur_ptn_total_poly3_simple) # pressupostos atingidos

plot(predict_response(femur_ptn_total_poly3_simple, "total_vig_and_mod_min"), 
     show_residuals=TRUE, show_residuals_line = TRUE,
     dot_size=.5) +
  labs(x="Minutos de atividade física\nmoderada a vigorosa", y="DMO colo femural", 
       title="Valores preditos de proteína vs. resíduos parciais") +
  theme_avp() # boa concordância da linha predita com a linha de resíduos

# checando interação

femur_interaction <- svyglm(DXXNKBMD ~ PTN * poly(total_vig_and_mod_min, 3) + CHO + LIP + 
                            RIAGENDR + BMXWT + BMXHT +
                            + RIDAGEYR + days_gc,
                          design=NHANES_subset_femur) 

femur_ptn_total_ns3_simple <- svyglm(DXXNKBMD ~ PTN + splines::ns(total_vig_and_mod_min, df=3) + CHO + LIP + 
                                       RIAGENDR + BMXWT + BMXHT
                                     + RIDAGEYR + days_gc,
                                     design=NHANES_subset_femur) 

car::Anova(femur_interaction, type="II", test="F")

summary(femur_interaction)

AIC(femur_interaction, femur_ptn_total_poly3_simple)

anova(femur_interaction, femur_ptn_total_poly3_simple) # O modelo com interação abaixa a AIC ainda mais 
                                                       # e apresenta diferença significativa no F test.

femur_interaction_centered <- svyglm(DXXNKBMD ~ scale(PTN, center=TRUE, scale=FALSE) * 
                                       scale(poly(total_vig_and_mod_min, 3), center=TRUE, scale=FALSE)
                                     + CHO + LIP + 
                              RIAGENDR + BMXWT + BMXHT +
                              + RIDAGEYR + days_gc,
                            design=NHANES_subset_femur) 

svy_check_model(femur_interaction_centered)

ggsave("figuras/femur_finalmodel_assump.png", dpi=600, unit="in", height=8, width=8)

# gráficos de predição

pred_1 <- ggpredict(femur_interaction, 
                    terms=c("PTN[1:500]", "total_vig_and_mod_min[0:500]")) |> as_tibble() |> # predições nos valores específicos 
  rename(PTN = x, atv_fisica = group)

pred_1_selected <- pred_1 |> filter(atv_fisica %in% c(0, 100, 200, 300, 400, 500)) 

#colors <- c("#A7D8FF", "#A7D8FF", "#7CBFFF", "#529EFF", "#297EFF", "#005FFF")

(
  femur_pred_p_1 <- ggplot(pred_1_selected, aes(x=PTN, y=predicted, color=atv_fisica, fill=atv_fisica,
                                             ymin=conf.low, ymax=conf.high)) +
    facet_wrap(~atv_fisica, labeller = labeller(atv_fisica = ~ paste(.x, "min/sem"))) +
    geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.35, color=NA) +
    geom_line() +
    labs(x="Ingestão de proteína (g/dia)", y="DMO de colo femoral (g/cm²)", 
         color="Atividade física\nmoderada e vigorosa\n(min/sem)",
         fill="Atividade física\nmoderada e vigorosa\n(min/sem)") +
    scale_color_discrete_sequential(labels=paste0(seq(0, 500, 100), " min/sem"), palette = "Reds", l1 = 25, c2 = 150, p1 = 1) + 
    scale_fill_discrete_sequential(labels=paste0(seq(0, 500, 100), " min/sem"), palette = "Reds", l1 = 25, c2 = 150, p1 = 1) + 
    scale_x_continuous(limits=c(25, 250)) +
    guides(color='none', fill='none') +
    theme_avp()
)

# gráfico de slope para visualizar interação

pred_1_slopes <- pred_1 |> filter(PTN %in% c(1, 2)) |> 
  select(PTN, predicted, atv_fisica) |> 
  pivot_wider(values_from=predicted, names_from=PTN) |> 
  mutate(slope = `2` - `1`) |> 
  mutate(atv_fisica = as.numeric(atv_fisica))

(
  femur_pred_p_2 <- ggplot(pred_1_slopes, aes(x=atv_fisica, y=slope*100)) +
    geom_line(aes(group=1)) +
    #geom_vline(xintercept=230, linetype="dashed", alpha=.25) +
    #annotate(geom='label', label='230 min.', x=230, y=19.85, size=2.5) +
    scale_x_continuous(breaks=seq(0, 500, 100)) +
    scale_y_continuous(breaks=c(seq(0, 0.05, 0.01)),
                       limits=c(0, 0.05)) +
    theme_avp() +
    labs(x="Atividade física moderada e vigorosa (min/sem)", y="Slope para aumento de\n 100 g de proteína")
)

(femur_pred_p_1 / femur_pred_p_2) +
  plot_annotation(tag_levels="A") +
  plot_layout(height=c(1, .5))

#ggsave("figuras/plot2.png", dpi=600, unit="in", height=8, width=6)

## analise binaria sarcopenia e osteoporose

# sarc

sarc_model_ptn_simple <- svyglm(sarc_class ~ PTN + total_vig_and_mod_min + CHO + LIP + 
                               RIAGENDR + BMXWT + BMXHT +
                               + RIDAGEYR + days_gc, family=quasibinomial(),
                             design=NHANES_subset_sarc)

summary(sarc_model_ptn_simple)

car::Anova(sarc_model_ptn_simple, "II", "F")

svy_check_model(sarc_model_ptn_simple)

plot(predict_response(sarc_model_ptn_simple, "PTN [all]"), 
     show_residuals=TRUE, show_residuals_line = TRUE,
     dot_size=.5) + # Usando partial regressions plots para visualizar fits inadequados
  labs(x="Ingestão de proteína", y="Probabilidade de sarcopenia", 
       title="Valores preditos de proteína vs. resíduos parciais") +
  theme_avp()

sarc_model_ptnPoly2_simple <- svyglm(sarc_class ~ poly(PTN, 2) + total_vig_and_mod_min + CHO + LIP + 
                                       RIAGENDR + BMXWT + BMXHT +
                                       + RIDAGEYR + days_gc, family=quasibinomial(),
                                     design=NHANES_subset_sarc)

sarc_model_ptnPoly3_simple <- svyglm(sarc_class ~ poly(PTN, 3) + total_vig_and_mod_min + CHO + LIP + 
                                       RIAGENDR + BMXWT + BMXHT +
                                       + RIDAGEYR + days_gc, family=quasibinomial(),
                                     design=NHANES_subset_sarc) 

sarc_model_ptnNs3_simple <- svyglm(sarc_class ~ ns(PTN, 3) + total_vig_and_mod_min + CHO + LIP + 
                                     RIAGENDR + BMXWT + BMXHT +
                                     + RIDAGEYR + days_gc, family=quasibinomial(),
                                   design=NHANES_subset_sarc) 

AIC(sarc_model_ptn_simple, sarc_model_ptnPoly2_simple,
    sarc_model_ptnPoly3_simple, sarc_model_ptnNs3_simple) # Modelo com polinomial 2 grau parece ter melhor fit

anova(sarc_model_ptn_simple, sarc_model_ptnPoly2_simple) 
anova(sarc_model_ptn_simple, sarc_model_ptnPoly3_simple) 
anova(sarc_model_ptn_simple, sarc_model_ptnNs3_simple) 

car::Anova(sarc_model_ptnPoly2_simple, "II", "F")

predictions(sarc_model_ptnPoly2_simple, by="PTN",
            newdata=get_pred_df(NHANES_subset_sarc, predictor="PTN"),
            wts=NHANES_subset_sarc$variables$weights_all) |> 
  mutate(conf.low = if_else(conf.low < 0, 0, conf.low),
         conf.high = if_else(conf.high > 1, 1, conf.high)) |> 
  ggplot(aes(x=PTN, y=estimate, ymin=conf.low, ymax=conf.high)) +
  geom_line() +
  xlim(c(0, 250)) +
  geom_ribbon(alpha=.1) +
  theme_avp() # Visualização do efeito

## olhando para atv física

plot(predict_response(sarc_model_ptn_simple, "total_vig_and_mod_min [all]"), 
     show_residuals=TRUE, show_residuals_line = TRUE,
     dot_size=.5) + # Usando partial regressions plots para visualizar fits inadequados
  labs(x="Atv física", y="Probabilidade de sarcopenia", 
       title="Valores preditos de proteína vs. resíduos parciais") +
  theme_avp()

sarc_model_atv_simple <- svyglm(sarc_class ~ poly(PTN, 2) + total_vig_and_mod_min + CHO + LIP + 
                                  RIAGENDR + BMXWT + BMXHT +
                                  + RIDAGEYR + days_gc, family=quasibinomial(),
                                design=NHANES_subset_sarc)

summary(sarc_model_atv_simple)

car::Anova(sarc_model_atv_simple, "II", "F")

sarc_model_atvPoly2_simple <- svyglm(sarc_class ~ poly(PTN, 2) + poly(total_vig_and_mod_min, 2) + CHO + LIP + 
                                       RIAGENDR + BMXWT + BMXHT +
                                       + RIDAGEYR + days_gc, family=quasibinomial(),
                                     design=NHANES_subset_sarc) 

sarc_model_atvPoly3_simple <- svyglm(sarc_class ~ poly(PTN, 2) + poly(total_vig_and_mod_min, 3) + CHO + LIP + 
                                     RIAGENDR + BMXWT + BMXHT +
                                     + RIDAGEYR + days_gc, family=quasibinomial(),
                                   design=NHANES_subset_sarc) 

sarc_model_atvNs3_simple <- svyglm(sarc_class ~ poly(PTN, 2) + ns(total_vig_and_mod_min, 3) + CHO + LIP + 
                                     RIAGENDR + BMXWT + BMXHT +
                                     + RIDAGEYR + days_gc, family=quasibinomial(),
                                   design=NHANES_subset_sarc) 

AIC(sarc_model_atv_simple,
    sarc_model_atvPoly2_simple,
    sarc_model_atvPoly3_simple,
    sarc_model_atvNs3_simple) ## Modelo linear parece o mais adequado

# Devido a erros no predict svyglm, precisamos recomputar o modelo com os valores polinomiais já calculados

NHANES_subset_sarc <- update(NHANES_subset_sarc, PTN_poly1 = poly(NHANES_subset_sarc$variables$PTN, 2)[,1])
NHANES_subset_sarc <- update(NHANES_subset_sarc, PTN_poly2 = poly(NHANES_subset_sarc$variables$PTN, 2)[,2])

sarc_model_atvPoly2_simple_refitted <- svyglm(sarc_class ~ PTN_poly1 + PTN_poly2 + poly(total_vig_and_mod_min, 2) + CHO + LIP + 
                                       RIAGENDR + BMXWT + BMXHT +
                                       + RIDAGEYR + days_gc, family=quasibinomial(),
                                     design=NHANES_subset_sarc) 

ggpredict(sarc_model_atvPoly2_simple_refitted, "total_vig_and_mod_min [all]") |> plot()

# interação

sarc_interaction <- svyglm(sarc_class ~ poly(PTN, 2) * poly(total_vig_and_mod_min, 2) + CHO + LIP + 
                             RIAGENDR + BMXWT + BMXHT +
                             + RIDAGEYR + days_gc, family=quasibinomial(),
                           design=NHANES_subset_sarc) 

car::Anova(sarc_interaction, "II", "F")

AIC(sarc_model_atvPoly2_simple, sarc_interaction) ## interação não significativa e AIC pior no modelo 
anova(sarc_model_atvPoly2_simple, sarc_interaction)

sarc_interaction_centered <- svyglm(sarc_class ~ scale(poly(PTN, 2), center=TRUE, scale=FALSE) * 
                                     scale(poly(total_vig_and_mod_min, 2), center=TRUE, scale=FALSE) + 
                                      CHO + LIP + 
                                     RIAGENDR + BMXWT + BMXHT +
                                     + RIDAGEYR + days_gc, family=quasibinomial(),
                                   design=NHANES_subset_sarc)

# gráficos de predição

pred_1 <- ggpredict(sarc_interaction, 
                    terms=c("PTN[1:500]", "total_vig_and_mod_min[0:500]")) |> as_tibble() |> # predições nos valores específicos 
  rename(PTN = x, atv_fisica = group)

pred_1_selected <- pred_1 |> filter(atv_fisica %in% c(0, 100, 200, 300, 400, 500)) 

#colors <- c("#A7D8FF", "#A7D8FF", "#7CBFFF", "#529EFF", "#297EFF", "#005FFF")

(
  sarc_pred_p_1 <- ggplot(pred_1_selected, aes(x=PTN, y=predicted, color=atv_fisica, fill=atv_fisica,
                                                ymin=conf.low, ymax=conf.high)) +
    facet_wrap(~atv_fisica, labeller = labeller(atv_fisica = ~ paste(.x, "min/sem"))) +
    geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.25, color=NA) +
    geom_line() +
    labs(x="Ingestão de proteína (g/dia)", y="Probabilidade de sarcopenia (%)", 
         color="Atividade física\nmoderada e vigorosa\n(min/sem)",
         fill="Atividade física\nmoderada e vigorosa\n(min/sem)") +
    scale_color_discrete_sequential(labels=paste0(seq(0, 500, 100), " min/sem"), palette = "Reds", l1 = 25, c2 = 150, p1 = 1) + 
    scale_fill_discrete_sequential(labels=paste0(seq(0, 500, 100), " min/sem"), palette = "Reds", l1 = 25, c2 = 150, p1 = 1) + 
    scale_x_continuous(limits=c(25, 250)) +
    scale_y_continuous(breaks=seq(0.06, 0.13, 0.01),
                       labels=paste0(seq(6, 13, 1), "%")) +
    guides(color='none', fill='none') +
    theme_avp()
)

ggsave("figuras/plot_sarc_1.png", dpi=600, unit="in", height=6, width=6)

# osteo

osteo_model_ptn_simple <- svyglm(femurneck_osteoporosis ~ PTN + total_vig_and_mod_min + CHO + LIP + 
                                  RIAGENDR + BMXWT + BMXHT +
                                  + RIDAGEYR + days_gc, family=quasibinomial(),
                                design=NHANES_subset_femur)

summary(osteo_model_ptn_simple)

car::Anova(osteo_model_ptn_simple, "II", "F")

plot(predict_response(sarc_model_ptn_simple, "PTN [all]"), 
     show_residuals=TRUE, show_residuals_line = TRUE,
     dot_size=.5) + # Usando partial regressions plots para visualizar fits inadequados
  labs(x="Ingestão de proteína", y="Probabilidade de osteoporose", 
       title="Valores preditos de proteína vs. resíduos parciais") +
  theme_avp()

osteo_model_ptnPoly2_simple <- svyglm(femurneck_osteoporosis ~ poly(PTN, 2) + total_vig_and_mod_min + CHO + LIP + 
                                       RIAGENDR + BMXWT + BMXHT +
                                       + RIDAGEYR + days_gc, family=quasibinomial(),
                                     design=NHANES_subset_femur)

osteo_model_ptnPoly3_simple <- svyglm(femurneck_osteoporosis ~ poly(PTN, 3) + total_vig_and_mod_min + CHO + LIP + 
                                       RIAGENDR + BMXWT + BMXHT +
                                       + RIDAGEYR + days_gc, family=quasibinomial(),
                                     design=NHANES_subset_femur) 

osteo_model_ptnNs3_simple <- svyglm(femurneck_osteoporosis ~ ns(PTN, 3) + total_vig_and_mod_min + CHO + LIP + 
                                     RIAGENDR + BMXWT + BMXHT +
                                     + RIDAGEYR + days_gc, family=quasibinomial(),
                                   design=NHANES_subset_femur) 

AIC(osteo_model_ptn_simple, osteo_model_ptnPoly2_simple,
    osteo_model_ptnPoly3_simple, osteo_model_ptnNs3_simple) 

anova(osteo_model_ptn_simple, osteo_model_ptnPoly3_simple) 
anova(osteo_model_ptn_simple, osteo_model_ptnNs3_simple) 

car::Anova(sarc_model_ptnNs3_simple, "II", "F")

predictions(sarc_model_ptnNs3_simple, by="PTN",
            newdata=get_pred_df(NHANES_subset_femur, predictor="PTN"),
            wts=NHANES_subset_femur$variables$weights_all) |> 
   mutate(conf.low = if_else(conf.low < 0, 0, conf.low),
          conf.high = if_else(conf.high > 1, 1, conf.high)) |> 
  ggplot(aes(x=PTN, y=estimate, ymin=conf.low, ymax=conf.high)) +
  geom_line() +
  xlim(c(0, 250)) +
  geom_ribbon(alpha=.1) +
  theme_avp() # Visualização do efeito

## olhando para atv física

plot(predict_response(osteo_model_ptn_simple, "total_vig_and_mod_min [all]"), 
     show_residuals=TRUE, show_residuals_line = TRUE,
     dot_size=.5) + # Usando partial regressions plots para visualizar fits inadequados
  labs(x="Atv física", y="Probabilidade de sarcopenia", 
       title="Valores preditos de proteína vs. resíduos parciais") +
  ylim(0, .05) +
  theme_avp()

osteo_model_atv_simple <- svyglm(femurneck_osteoporosis ~ ns(PTN, 3) + total_vig_and_mod_min + CHO + LIP + 
                                  RIAGENDR + BMXWT + BMXHT +
                                  + RIDAGEYR + days_gc, family=quasibinomial(),
                                design=NHANES_subset_femur)

summary(osteo_model_atv_simple)

car::Anova(osteo_model_atv_simple, "II", "F")

osteo_model_atvPoly2_simple <- svyglm(femurneck_osteoporosis ~ ns(PTN, 3) + poly(total_vig_and_mod_min, 2) + CHO + LIP + 
                                       RIAGENDR + BMXWT + BMXHT +
                                       + RIDAGEYR + days_gc, family=quasibinomial(),
                                     design=NHANES_subset_femur) 

osteo_model_atvPoly3_simple <- svyglm(femurneck_osteoporosis ~ ns(PTN, 3) + poly(total_vig_and_mod_min, 3) + CHO + LIP + 
                                       RIAGENDR + BMXWT + BMXHT +
                                       + RIDAGEYR + days_gc, family=quasibinomial(),
                                     design=NHANES_subset_femur) 

osteo_model_atvNs3_simple <- svyglm(femurneck_osteoporosis ~ ns(PTN, 3) + ns(total_vig_and_mod_min, 3) + CHO + LIP + 
                                     RIAGENDR + BMXWT + BMXHT +
                                     + RIDAGEYR + days_gc, family=quasibinomial(),
                                   design=NHANES_subset_femur) 

AIC(osteo_model_atv_simple,
    osteo_model_atvPoly2_simple,
    osteo_model_atvPoly3_simple,
    osteo_model_atvNs3_simple) ## Modelo poly 3 parece mais adequado

anova(osteo_model_atv_simple, osteo_model_atvPoly3_simple)

# Devido a erros no predict svyglm, precisamos recomputar o modelo com os valores do spline já calculados

NHANES_subset_femur <- update(NHANES_subset_femur, PTN_ns1 = splines::ns(NHANES_subset_femur$variables$PTN, 3)[,1])
NHANES_subset_femur <- update(NHANES_subset_femur, PTN_ns2 = splines::ns(NHANES_subset_femur$variables$PTN, 3)[,2])
NHANES_subset_femur <- update(NHANES_subset_femur, PTN_ns3 = splines::ns(NHANES_subset_femur$variables$PTN, 3)[,3])

osteo_model_atvPoly3_simple_refitted <- svyglm(femurneck_osteoporosis ~ PTN_ns1 +
                                                 PTN_ns2 + PTN_ns3 + poly(total_vig_and_mod_min, 3) + CHO + LIP + 
                                                 RIAGENDR + BMXWT + BMXHT +
                                                 + RIDAGEYR + days_gc, family=quasibinomial(),
                                               design=NHANES_subset_femur) 

pred_df <- get_pred_df(NHANES_subset_femur, predictor="PA") |> 
  mutate(PTN_ns1 = mean(NHANES_subset_femur$variables$PTN_ns1),
         PTN_ns2 = mean(NHANES_subset_femur$variables$PTN_ns2),
         PTN_ns3 =  mean(NHANES_subset_femur$variables$PTN_ns3))

car::Anova(osteo_model_atvPoly3_simple, type="II", "F")

predictions(osteo_model_atvPoly3_simple_refitted, by="total_vig_and_mod_min",
            newdata=pred_df,
            wts=NHANES_subset_femur$variables$weights_all) |> 
  mutate(conf.low = if_else(conf.low < 0, 0, conf.low),
         conf.high = if_else(conf.high > 1, 1, conf.high)) |> 
  ggplot(aes(x=total_vig_and_mod_min, y=estimate, ymin=conf.low, ymax=conf.high)) +
  geom_line() +
  geom_ribbon(alpha=.1) +
  xlim(c(0, 3000)) +
  theme_avp() # Visualização do efeito 

# interação

osteo_interaction <- svyglm(femurneck_osteoporosis ~ ns(PTN, 3) * poly(total_vig_and_mod_min, 3) + CHO + LIP + 
                             RIAGENDR + BMXWT + BMXHT +
                             + RIDAGEYR + days_gc, family=quasibinomial,
                           design=NHANES_subset_femur) 

car::Anova(osteo_interaction, "II", "F")

AIC(osteo_model_atvPoly3_simple, osteo_interaction) 
anova(osteo_model_atvPoly3_simple, osteo_interaction) 

osteo_interaction_centered <- svyglm(femurneck_osteoporosis ~ scale(ns(PTN, 3), center=TRUE, scale=FALSE) * 
                                       scale(poly(total_vig_and_mod_min, 3), center=TRUE, scale=FALSE) + CHO + LIP + 
                              RIAGENDR + BMXWT + BMXHT +
                              + RIDAGEYR + days_gc, family=quasibinomial,
                            design=NHANES_subset_femur) 

# gráficos de predição

pred_1 <- ggpredict(osteo_interaction, 
                    terms=c("PTN[1:500]", "total_vig_and_mod_min[0:500]")) |> as_tibble() |> # predições nos valores específicos 
  rename(PTN = x, atv_fisica = group)

pred_1_selected <- pred_1 |> filter(atv_fisica %in% c(0, 100, 200, 300, 400, 500)) 

#colors <- c("#A7D8FF", "#A7D8FF", "#7CBFFF", "#529EFF", "#297EFF", "#005FFF")

(
  osteo_pred_p_1 <- ggplot(pred_1_selected, aes(x=PTN, y=predicted, color=atv_fisica, fill=atv_fisica,
                                               ymin=conf.low, ymax=conf.high)) +
    facet_wrap(~atv_fisica, labeller = labeller(atv_fisica = ~ paste(.x, "min/sem"))) +
    geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.35, color=NA) +
    geom_line() +
    labs(x="Ingestão de proteína (g/dia)", y="Probabilidade de osteoporose (%)", 
         color="Atividade física\nmoderada e vigorosa\n(min/sem)",
         fill="Atividade física\nmoderada e vigorosa\n(min/sem)") +
    scale_color_discrete_sequential(labels=paste0(seq(0, 500, 100), " min/sem"), palette = "Reds", l1 = 25, c2 = 150, p1 = 1) + 
    scale_fill_discrete_sequential(labels=paste0(seq(0, 500, 100), " min/sem"), palette = "Reds", l1 = 25, c2 = 150, p1 = 1) + 
    scale_x_continuous(limits=c(25, 250)) +
    # scale_y_continuous(breaks=seq(0, 0.13, 0.01),
    #                    labels=paste0(seq(6, 13, 1), "%")) +
    guides(color='none', fill='none') +
    theme_avp()
)

ggsave("figuras/plot2_int_osteo_full.png", dpi=600, unit="in", height=6, width=6)

##

pred_2_selected <- pred_1 |> filter(atv_fisica %in% c(0, 100, 200, 300)) 

#colors <- c("#A7D8FF", "#A7D8FF", "#7CBFFF", "#529EFF", "#297EFF", "#005FFF")

(
  osteo_pred_p_2 <- ggplot(pred_2_selected, aes(x=PTN, y=predicted, color=atv_fisica, fill=atv_fisica,
                                                ymin=conf.low, ymax=conf.high)) +
    facet_wrap(~atv_fisica, labeller = labeller(atv_fisica = ~ paste(.x, "min/sem"))) +
    geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.35, color=NA) +
    geom_line() +
    labs(x="Ingestão de proteína (g/dia)", y="Probabilidade de sarcopenia (%)", 
         color="Atividade física\nmoderada e vigorosa\n(min/sem)",
         fill="Atividade física\nmoderada e vigorosa\n(min/sem)") +
    scale_color_discrete_sequential(labels=paste0(seq(0, 500, 100), " min/sem"), palette = "Reds", l1 = 25, c2 = 150, p1 = 1) + 
    scale_fill_discrete_sequential(labels=paste0(seq(0, 500, 100), " min/sem"), palette = "Reds", l1 = 25, c2 = 150, p1 = 1) + 
    scale_x_continuous(limits=c(25, 250)) +
    scale_y_continuous(breaks=seq(0, .50, 0.10),
                       labels=paste0(seq(0, 50, 10), "%"),
                       limits=c(0, .5)) +
    guides(color='none', fill='none') +
    theme_avp()
)

ggsave("figuras/plot2_int_osteo_zoom.png", dpi=600, unit="in", height=6, width=6)

# tentando modelo semelhante com grupos em quantis

# calculando quantis com survey

PTN_tile <- svyquantile(~PTN, design=NHANES_subset_femur, quantiles=c(.25, .5, .75))
atv_tile <- svyquantile(~total_vig_and_mod_min, design=NHANES_subset_femur, quantiles=c(.33, .66))

NHANES_subset_femur <- update(NHANES_subset_femur, PTN_tiles = as.factor(case_when(PTN < PTN_tile$PTN[1] ~ "1",
                                                                         PTN >= PTN_tile$PTN[1] & PTN < PTN_tile$PTN[2] ~ "2",
                                                                         PTN >= PTN_tile$PTN[2] & PTN < PTN_tile$PTN[3] ~ "3",
                                                                         PTN >= PTN_tile$PTN[4] ~ "4")))

NHANES_subset_femur <- update(NHANES_subset_femur, mvpa_tiles = as.factor(case_when(total_vig_and_mod_min < atv_tile$total_vig_and_mod_min[1] ~ "1",
                                                                          total_vig_and_mod_min >= atv_tile$total_vig_and_mod_min[1] & total_vig_and_mod_min < atv_tile$total_vig_and_mod_min[2] ~ "2",
                                                                          total_vig_and_mod_min >= atv_tile$total_vig_and_mod_min[3] ~ "3")))

osteo_interaction_cat <- svyglm(femurneck_osteoporosis ~ PTN_tiles * mvpa_tiles + CHO + LIP + 
                              RIAGENDR + BMXWT + BMXHT +
                              + RIDAGEYR + days_gc, family=quasibinomial(),
                            design=NHANES_subset_femur) 

summary(osteo_interaction_cat)
car::Anova(osteo_interaction_cat, type="II", test="F")

osteo_int_pred_df <- ggaverage(osteo_interaction_cat, c("PTN_tiles", "mvpa_tiles")) |> as_tibble() |> 
  rename(ptn = x, mvpa = group) 

ggplot(osteo_int_pred_df, aes(x=ptn, y=predicted, ymin=conf.low, ymax=conf.high, color=mvpa)) +
  #facet_wrap(~mvpa, labeller=as_labeller(c("1" = "1º tercil de atv. física", 
   #                                       "2" = "2º tercil de atv. física", 
   #                                       "3" =  "3º tercil de atv. física"))) +
  geom_pointrange(position=position_dodge(.35)) +
  #geom_line(aes(group=mvpa), position=position_dodge(.35)) +
  scale_y_continuous(limits=c(0, .1),
                    breaks=seq(0, .1, .02),
                    labels=paste0(seq(0, 10, 2), "%")) +
  scale_color_discrete_sequential(palette = "Reds", l1 = 25, c2 = 150, p1 = 1,
                                  labels=c("1º tercil de atv. física", "2º tercil de atv. física", 
                                           "3º tercil de atv. física")) + 
  labs(x="Quartis de consumo de proteína", y="Probabilidade de osteoporose (%)",
       color="Tercis de atividade\nfísica moderada\ne vigorosa") +
  scale_x_discrete(labels=c("1º quartil", "2º quartil", "3º quartil", "4º quartil")) +
  theme_avp()

ggsave("figuras/plot2_int_osteo_cat.png", dpi=600, unit="in", height=4.5, width=7)

## Vifs

car::vif(sarc_interaction)

car::vif(osteo_interaction)

#### Tables ####################################################################

lm_dat <- as_tibble(NHANES_subset_lm$variables)

nrow(lm_dat)

fem_dat <- as_tibble(NHANES_subset_femur$variables)

nrow(fem_dat)

sample_dat <- reduce(list(lm_dat, fem_dat), bind_rows) |> distinct()

nrow(sample_dat) # sample size

# creating table 1 with gtsummary

theme_gtsummary_language(language = "pt",
                         big.mark="", 
                         decimal.mark=",")

table_1_summary <- sample_dat |>
  select(RIDAGEYR, RIAGENDR, RIDRETH1, BMXHT, BMXWT, BMXBMI, days_gc, 
         total_vig_and_mod_min, 
         kcal, CHO, PTN, LIP, DXDTOLE, ALM, ALM_BMI, sarc_class, DXXNKBMD, 
         femurneck_t_score, femurneck_osteoporosis) |>
  mutate(DXDTOLE = DXDTOLE/1000,
         ALM = ALM/1000) |> 
  tbl_summary(by = RIAGENDR, missing = "no", statistic = all_continuous() ~ "{median} ({p25}, {p75})",
              digits = list(RIDAGEYR ~ 0, 
                            BMXWT ~ 1,
                            BMXHT ~ 1,
                            BMXBMI ~ 1,
                            days_gc ~ 0,
                            total_vig_and_mod_min ~ 0,
                            kcal ~ 0,
                            PTN ~ 0,
                            LIP ~ 0,
                            DXDTOLE ~ 1,
                            DXXNKBMD ~ 3,
                            femurneck_t_score ~ 1,
                            ALM ~ 1,
                            ALM_BMI ~ 1)) |>
  add_overall()

table_1_summary

# save
# table_1_summary |> as_hux_table() |> writexl::write_xlsx("tabelas/tabela1.xlsx")
