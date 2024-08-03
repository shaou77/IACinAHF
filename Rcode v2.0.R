library(tidyverse)
library(naniar)
library(twang)
library(Matching)
library(xgboost)
library(ggplot2)
library(viridis)
library(survey)
library(tableone)
library(dplyr)
# library(finalfit)
library(readr)
data <- read_csv("data.csv")
##数据处理，计算生存时间
data <- dplyr::distinct(data) %>% dplyr::select(-1:-3) %>% mutate(dialysis_active = ifelse(is.na(dialysis_active), 0, dialysis_active),
                                                    surtime_from_admit= ifelse(is.na(surtime_from_admit), 9999, surtime_from_admit),
                                                    surtime_from_icu= ifelse(is.na(surtime_from_icu), 9999, surtime_from_icu),
                                                    death_90day = ifelse(surtime_from_icu<=90,1,0)) %>% 
  dplyr::select(-troponin_t,-ntprobnp,-gcs_min,-oasis,-dialysis_active)


##数据处理和队列筛选
data2 <- data %>% mutate(death_28day = ifelse(surtime_from_icu<=28,1,0)) %>% filter(line_before<1) %>% filter(los_icu>=1) %>% filter(vaso_flag==0) %>% 
  mutate(uokg=urineoutput/weight) %>% mutate(diabetes = if_else(diabetes_with_cc == 1 | diabetes_without_cc == 1, 1, 0)) %>% 
  mutate(liver_disease=if_else(mild_liver_disease == 1 | severe_liver_disease == 1, 1, 0)) %>% mutate(cancer=if_else(malignant_cancer == 1 | metastatic_solid_tumor == 1, 1, 0)) %>% 
  dplyr::select(-diabetes_with_cc, -diabetes_without_cc, -mild_liver_disease,-severe_liver_disease,-malignant_cancer,-metastatic_solid_tumor,-urineoutput,-weight, -calcium_mean, -resp_rate_mean,-vaso_flag,-congestive_heart_failure) %>% 
  mutate(icufreeday = case_when(
    (abs(surtime_from_icu - los_icu) <= 1 & hospital_expire_flag == 1 & los_icu <= 28) ~ 0,
    (surtime_from_icu <= 28 & hospital_expire_flag == 1) ~ 0,
    TRUE ~ pmax(0, 28 - los_icu)
  ))


#line_before,-1为没有置管，0为在ICU中置管，1为在ICU前置管
table(data2$line_before)
##查看缺失值和去除缺失值过多的特征
print(miss_var_summary(data2), n=40)
miss_summary <- miss_var_summary(data2)
data2$line_before <- data2$line_before+1
cols_to_keep <- miss_summary %>%  filter(pct_miss <= 3) %>%  pull(variable) #大于3%的均去掉
data1_filter <- data2 %>% dplyr::select(all_of(cols_to_keep)) %>% filter(complete.cases(.)) %>% 
  dplyr::select(-temperature_mean,-heart_rate_mean,-liver_disease,-mbp_mean,-spo2_mean,-peripheral_vascular_disease,-dementia,-rheumatic_disease)
table(data1_filter$line_before)
data2_c <- data1_filter %>%
  mutate(across(
    .cols = c("gender", "hospital_expire_flag", "line_before", 
              "ventilation", "myocardial_infarct", "cerebrovascular_disease",
              "chronic_pulmonary_disease", "renal_disease", "death_90day",
              "death_28day", "diabetes",  "cancer"),
    .fns = as.factor
  )) %>% as.data.frame()
##配对，使用xgboost
fml_w <- line_before~uokg+hemoglobin_mean+wbc_mean+platelets_mean+glucose_mean+aniongap_mean+bicarbonate_mean+bun_mean+chloride_mean+creatinine_mean+sodium_mean+potassium_mean+gender+admission_age+apsiii+ventilation+myocardial_infarct+cerebrovascular_disease+chronic_pulmonary_disease+renal_disease+diabetes+cancer
data2_c$line_before <- as.numeric(data2_c$line_before)-1
aline_ps_ate <- ps(as.formula(fml_w),
                  data = data2_c,
                  interaction.depth = 2,
                  shrinkage = 0.01,
                  perm.test.iters = 0,
                  estimand = "ATE",
                  verbose = FALSE,
                  stop.method = c("es.mean", "es.max", "ks.mean", "ks.max"),
                  n.trees = 10000,
                  train.fraction = 0.8,
                  cv.folds = 3,version = 'xgboost')
plot(aline_ps_ate,plots=5)

pred <- aline_ps_ate$ps$es.mean.ATE
data_c <- data2_c %>% mutate(ps = pred)
label <- data2_c %>% pull(line_before)
ROCR::performance(ROCR::prediction(pred, label), "auc")@y.values %>% first #配对效果
##特征重要度及绘图
ft_importance <- xgb.importance(model = aline_ps_ate$gbm.obj)
total_gain <- sum(ft_importance$Gain)
relative_influence <- ft_importance[, Gain := Gain / total_gain]
relative_influence <- relative_influence %>%
  mutate(Gain = as.numeric(Gain)) %>%
  arrange(desc(Gain))
relative_influence <- relative_influence %>% mutate(Feature = case_when(
    Feature == "ventilation1" ~ "Ventilation in first 24 hours",
    Feature == "admission_age" ~ "Age",
    Feature == "glucose_mean" ~ "Glucose",
    Feature == "apsiii" ~ "APSIII",
    Feature == "platelets_mean" ~ "Platelets",
    Feature == "uokg" ~ "Urine output",
    Feature == "bun_mean" ~ "BUN",
    Feature == "chloride_mean" ~ "Chloride",
    Feature == "wbc_mean" ~ "WBC",
    Feature == "creatinine_mean" ~ "Creatinine",
    Feature == "ntprobnp1" ~ "NT-proBNP",
    Feature == "hemoglobin_mean" ~ "Hemoglobin",
    Feature == "bicarbonate_mean" ~ "Bicarbonate",
    Feature == "sodium_mean" ~ "Sodium",
    Feature == "aniongap_mean" ~ "Anion Gap",
    Feature == "potassium_mean" ~ "Potassium",
    Feature == "myocardial_infarct1" ~ "Myocardial infarction",
    Feature == "gcs_min" ~ "GCS",
    Feature == "cancer1" ~ "Cancer",
    Feature == "chronic_pulmonary_disease1" ~ "Chronic pulmonary disease",
    Feature == "genderF" ~ "Gender",
    Feature == "cerebrovascular_disease1" ~ "Cerebrovascular disease",
    Feature == "diabetes1" ~ "Diabetes",
    Feature == "renal_disease1" ~ "Renal disease"))
relative_influence <- relative_influence %>%
  mutate(Feature = factor(Feature, levels = unique(Feature)))
ggplot(relative_influence, aes(y = Gain, x = reorder(Feature, Gain), fill = Gain)) +
  geom_bar(stat = "identity") +
  xlab("Relative Influence") +
  ylab("Feature") +
  theme(axis.text.y = element_text(hjust = 0)) +
  scale_fill_viridis_c(option = "G", direction = -1)+coord_flip()+theme_classic()
##主要结局,双重稳健
data_c <- data_c %>%
  mutate(ps_weight = get.weights(aline_ps_ate, stop.method = "es.mean"))
data_c$line_before <- as.factor(data_c$line_before)
#ipw
design <- svydesign(ids = ~1, weights = ~ps_weight, data = data_c)
result <- svyglm(death_28day ~ line_before, design, family = binomial())
summary(result)
exp(cbind(OR = coef(result), confint(result)))
#双重稳健，所有变量
data_c$id <- 1:nrow(data_c)
ipw_svydesign <- svydesign(ids = ~ id, weights = ~ ps_weight, data = data_c)
fml<- death_28day~line_before+uokg+hemoglobin_mean+wbc_mean+platelets_mean+glucose_mean+aniongap_mean+bicarbonate_mean+bun_mean+chloride_mean+creatinine_mean+sodium_mean+potassium_mean+gender+admission_age+apsiii+ventilation+myocardial_infarct+cerebrovascular_disease+chronic_pulmonary_disease+renal_disease+diabetes+cancer
logi <- svyglm(as.formula(fml),  family = quasibinomial, design = ipw_svydesign)
summary(logi)
exp(cbind(OR = coef(logi), confint(logi)))

##配对
ps_matches <- Match(Y = data_c$death_28day, Tr = as.integer(data_c$line_before)-1,
                    X = data_c$ps, M = 1, estimand = "ATT", caliper = 0.04,
                    exact = FALSE, replace = FALSE, ties = FALSE)
ps_matches_df <-
  data_c[ps_matches$index.treated, "id", drop = FALSE] %>%
  mutate(match = data_c[ps_matches$index.control, "id"]) %>%
  rename(trtd = id, ctrl = match) %>%
  mutate(match = 1:n()) %>%
  gather("group", "id", trtd, ctrl) %>%
  dplyr::select(id, group, match) %>%
  arrange(group, match)
ids <- data_c[unlist(ps_matches[c("index.control", "index.treated")]), "id"]
length(ids) / 2
head(ids)
summary(ps_matches)
outcome <- data.frame(echo_pt = data_c$death_28day[ps_matches$index.treated],
    match_pt = data_c$death_28day[ps_matches$index.control])
tab_mc <- table(outcome$echo_pt, outcome$match_pt, dnn = c("Aline", "Non-Aline"))
tab_mc[2, 1] / tab_mc[1, 2]
paste("95% Confint",
      round(exp(c(log(tab_mc[2, 1] / tab_mc[1, 2]) - qnorm(0.975) * sqrt(1 / tab_mc[1, 2] + 1 / tab_mc[2, 1]),
                  log(tab_mc[2, 1] / tab_mc[1, 2]) + qnorm(0.975) * sqrt(1 / tab_mc[1, 2] + 1 / tab_mc[2, 1]))), 7))
mcnemar.test(tab_mc) 
mortality <- data.frame(line = c(0, 1),
                        mortality = c((tab_mc[2, 2] + tab_mc[1, 2]) / sum(tab_mc),
                                      (tab_mc[2, 2] + tab_mc[2, 1]) / sum(tab_mc)))
mortality






##table 1
features <- c("line_before", "uokg", "hemoglobin_mean", "wbc_mean", "platelets_mean", "glucose_mean", "aniongap_mean", "bicarbonate_mean", "bun_mean", "chloride_mean", "creatinine_mean", "sodium_mean", "potassium_mean", "gender", "admission_age",  "apsiii", "ventilation", "myocardial_infarct", "cerebrovascular_disease",  "chronic_pulmonary_disease", "renal_disease", "diabetes", "cancer")
features_all <- c("line_before",'death_90day','death_28day','hospital_expire_flag', 'los_hospital',"uokg", "hemoglobin_mean", "wbc_mean", "platelets_mean", "glucose_mean", "aniongap_mean", "bicarbonate_mean", "bun_mean", "chloride_mean", "creatinine_mean", "sodium_mean", "potassium_mean",  "gender", "admission_age", "apsiii", "ventilation", "myocardial_infarct", "cerebrovascular_disease",  "chronic_pulmonary_disease", "renal_disease", "diabetes",  "cancer",'icufreeday')

tab <- CreateTableOne(vars = features_all,
                      strata = "line_before",
                      data = data2_c,
                      argsNormal = list(var.equal = FALSE))
capture.output(tab_df <- tab %>%
                 print(smd = TRUE) %>%
                 as.data.frame(stringsAsFactors = FALSE) %>%
                 dplyr::select(-test)) %>% invisible
tab_wtd <- svyCreateTableOne(vars = features,
                             strata = "line_before",
                             data = ipw_svydesign)
capture.output(tab_wtd_df <- tab_wtd %>%
                 print %>%
                 as.data.frame(stringsAsFactors = FALSE) %>%
                 dplyr::select(-test)) %>% invisible
ps_df <- ps_matches_df %>%
  left_join(data_c, by = "id") %>%
  arrange(desc(group), match)
tab_ps <- CreateTableOne(vars = features,
                         strata = "line_before",
                         data = ps_df)
capture.output(tab_ps_df <- tab_ps %>%
                 print(smd = TRUE) %>%
                 as.data.frame(stringsAsFactors = FALSE) %>%
                 dplyr::select(-test)) %>% invisible     
write.csv(tab_ps_df,file = 'table_m.csv')

##双重稳健分析，回归未均衡的变量
tab_wtd_df
fml_ub <- death_28day~line_before+ventilation
logi_ub <- svyglm(as.formula(fml_ub), family = quasibinomial, design = ipw_svydesign)
summary(logi_ub)
exp(cbind(OR = coef(logi_ub), confint(logi_ub)))
saveRDS(data_c, "data_c.rds")

##单因素+多因素
explanatory <- features
dependent <- 'death_28day'
data2_c %>%  summary_factorlist(dependent, explanatory, p=TRUE, add_dependent_label=TRUE) -> t1
knitr::kable(t1, row.names=FALSE, align=c("l", "l", "r", "r", "r"))
explanatory_multi <-  c('line_before','uokg',
                       'hemoglobin_mean','wbc_mean','aniongap_mean',
                      'bun_mean','chloride_mean','creatinine_mean','potassium_mean','gender',
                      'admission_age','apsiii','cerebrovascular_disease','renal_disease',
                      'cancer')
data2_c %>%  finalfit(dependent, explanatory, explanatory_multi) -> t2
knitr::kable(t2, row.names=FALSE, align=c("l", "l", "r", "r", "r", "r"))
#table1
table1 <- CreateTableOne(vars = features_all, strata = c("line_before"), data = data2_c, 
                         testExact = fisher.test, testNonNormal = kruskal.test) 
p1<-print(table1, nonnormal=c('icufreeday','los_hospital',"uokg", "hemoglobin_mean", "wbc_mean",  "platelets_mean", "glucose_mean", "aniongap_mean", "bicarbonate_mean", "bun_mean", "chloride_mean", "creatinine_mean", "sodium_mean", "potassium_mean", "admission_age",  "apsiii"),
          argsExact=c('death_90day','death_28day','hospital_expire_flag',  "gender",  "ventilation", "myocardial_infarct", "cerebrovascular_disease",  "chronic_pulmonary_disease", "renal_disease", "diabetes", "cancer")) 
write.csv(p1,'table1a.csv')
tableall <- CreateTableOne(vars = features_all, data = data2_c)
p2<-print(tableall,nonnormal=c('icufreeday','los_hospital',"uokg", "hemoglobin_mean", "wbc_mean",  "platelets_mean", "glucose_mean", "aniongap_mean", "bicarbonate_mean", "bun_mean", "chloride_mean", "creatinine_mean", "sodium_mean", "potassium_mean", "admission_age",  "apsiii"),
          argsExact=c('death_90day','death_28day','hospital_expire_flag',  "gender",  "ventilation", "myocardial_infarct", "cerebrovascular_disease",  "chronic_pulmonary_disease", "renal_disease", "diabetes", "cancer")) 
write.csv(p2,'table1b.csv')


max_index <- which.max(data2_c$glucose_mean)
data2_c$glucose_mean[max_index] <- mean(data2_c$glucose_mean)
max_index2 <- which.max(data_c$glucose_mean)
data_c$glucose_mean[max_index2] <- mean(data_c$glucose_mean)
