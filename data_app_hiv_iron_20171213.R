

# By John Ssenkusu
# Purpose: HIV LIVE Data application and iron data application,
# compute bootstrap CI and format them for LATEX

# Date created: December 05, 2016
# Date modified: May 24, 2017

library(sas7bdat)
# 
# hiv_est <- read.sas7bdat("C:/Users/John/Google Drive/PhD/Thesis/Paper 1/Data application/data_hiv_estimates.sas7bdat", debug=FALSE)
# 
# ACME_est_HIV <- rbind(hiv_est[1,1], hiv_est[1,9], hiv_est[1,5])
# ACME_SE_HIV <- rbind(hiv_est[1,2], hiv_est[1,10], hiv_est[1,6])
# ADE_est_HIV <- rbind(hiv_est[1,3], hiv_est[1,11], hiv_est[1,7])
# ADE_SE_HIV <- rbind(hiv_est[1,4], hiv_est[1,12], hiv_est[1,8])
# 
# HIV <- cbind(ACME_est_HIV, ACME_SE_HIV, ADE_est_HIV, ADE_SE_HIV)
# 
# #########################################
# 
# iron_est <- read.sas7bdat("C:/Users/John/Google Drive/PhD/Thesis/Paper 1/Data application/data_app_iron_estimates.sas7bdat", debug=FALSE)
# 
# ACME_est_iron <- rbind(iron_est[1,1], iron_est[1,9], iron_est[1,5])
# ACME_SE_iron <- rbind(iron_est[1,2], iron_est[1,10], iron_est[1,6])
# ADE_est_iron <- rbind(iron_est[1,3], iron_est[1,11], iron_est[1,7])
# ADE_SE_iron <- rbind(iron_est[1,4], iron_est[1,12], iron_est[1,8])
# 
# iron <- cbind(ACME_est_iron, ACME_SE_iron, ADE_est_iron, ADE_SE_iron)
# 
# #########################################
# 
# estimates <- rbind(iron, HIV)
# 
# ############## Printing table results to LATEX  ##################################
# 
# library(xtable)
# 
# colnames(estimates) <- c("Estimate", "SE", "Estimate", "SE")
# rownames(estimates) <- c("Normal calibration", "SNP calibration", "No calibration", 
#                          "1Normal calibration", "1SNP calibration", "1No calibration")
# print(xtable(estimates, digits=c(0,3,3,3,3), align="lcccc",
#              caption="ACME and ADE estimates for the malaria-cognition data and HIV-LIVE 
#                       data"), caption.placement = "top")

##################################################################################

########### COMPUTING BOOTSTRAP CI FOR IRON DATA #######################

iron_med <- read.sas7bdat("C:/Users/John/Google Drive/PhD/Thesis/Paper 1/Data application/Med_longdata_final_20170524.sas7bdat", debug=FALSE)

# Drop variables from data set
iron_med$visit6 <- NULL
iron_med$visit12 <- NULL
iron_med[ is.na(iron_med) ] <- NA # Make NaN become NA

library(nlme)
# Analysis to check interactions in the outcome model for the 
# iron-cognition data

summary(lme(mulout ~ as.factor(grouplog) + hgb + as.factor(visit) + agefin + 
              as.factor(child_educ)+as.factor(grouplog)*as.factor(visit), random =~ 1 | studyid, 
            na.action=na.exclude,data = iron_med))

summary(lme(mulout ~ as.factor(grouplog) + hgb + as.factor(visit) + agefin + 
              as.factor(child_educ)+hgb*as.factor(grouplog), random =~ 1 | studyid, na.action=na.exclude,data = iron_med))

summary(lme(mulout ~ as.factor(grouplog) + hgb + as.factor(visit) + agefin + 
              as.factor(child_educ)+hgb*as.factor(visit), random =~ 1 | studyid, na.action=na.exclude, data = iron_med))

summary(lme(mulout ~ as.factor(grouplog) + hgb + as.factor(visit) + agefin + 
              as.factor(child_educ), random =~ 1 | studyid, na.action=na.exclude,data = iron_med))


######## Checking identifiability assumption 4: No mediator-outcome confounder is
######## affected by the exposure

iron.age_assump4 <- glm(grouplog ~ agefin,family=binomial(link='logit'),data=subset(iron_med, visit==0))
iron.educ_assump4 <- glm(grouplog ~ as.factor(child_educ),family=binomial(link='logit'),data=subset(iron_med, visit==0))

summary(iron.age_assump4)
summary(iron.educ_assump4)

library(lme4)

############# ACME and ADE estimates without resampling  #################

model.hgb <- lmer(hgb ~ as.factor(grouplog) + as.factor(visit) +
                    as.factor(grouplog)*as.factor(visit)+ as.factor(child_educ) + ( 1 | studyid), 
                  na.action=na.exclude, data = iron_med)
# Uncalibrated
outcome_error <- lmer(mulout ~ as.factor(grouplog) + hgb + as.factor(visit) + agefin +
                        as.factor(child_educ) + ( 1 | studyid), na.action=na.exclude, data = iron_med)

ACME.error_BL <- fixef(model.hgb)[2] * fixef(outcome_error)[3]
ACME.error_6 <- (fixef(model.hgb)[2] + fixef(model.hgb)[6])* fixef(outcome_error)[3]
ACME.error_12 <- (fixef(model.hgb)[2] + fixef(model.hgb)[7])* fixef(outcome_error)[3]
ADE.error <- fixef(outcome_error)[2]

# Calibrated assuming normality
outcome_norm <- lmer(mulout ~ grouplog + hgb_Norm_calib + as.factor(visit) + agefin +
                       as.factor(child_educ) + ( 1 | studyid), na.action=na.exclude, 
                     data = iron_med)

ACME.norm_BL <- fixef(model.hgb)[2] * fixef(outcome_norm)[3]
ACME.norm_6 <- (fixef(model.hgb)[2] + fixef(model.hgb)[6])* fixef(outcome_norm)[3]
ACME.norm_12 <- (fixef(model.hgb)[2] + fixef(model.hgb)[7])* fixef(outcome_norm)[3]
ADE.norm <- fixef(outcome_norm)[2]

# Calibrated assuming SNP
outcome_SNP <- lmer(mulout ~ grouplog + hgb_SNP_calib + as.factor(visit) + agefin +
                      as.factor(child_educ)+( 1 | studyid), na.action=na.exclude,
                    data = iron_med)

ACME.SNP_BL <- fixef(model.hgb)[2] * fixef(outcome_SNP)[3]
ACME.SNP_6 <- (fixef(model.hgb)[2] + fixef(model.hgb)[6])* fixef(outcome_SNP)[3]
ACME.SNP_12 <- (fixef(model.hgb)[2] + fixef(model.hgb)[7])* fixef(outcome_SNP)[3]
ADE.SNP <- fixef(outcome_SNP)[2]

iron_est <- rbind(ACME.norm_BL, ACME.norm_6, ACME.norm_12, 
                  ACME.SNP_BL, ACME.SNP_6, ACME.SNP_12,
                  ACME.error_BL, ACME.error_6, ACME.error_12,
                  ADE.norm, ADE.SNP, ADE.error)

############################################################################

iron_dt_wide <- reshape(iron_med, 
             timevar = "visit",
             idvar = c("studyid", "grouplog","child_educ"),
             direction = "wide")

dt.size <- nrow(iron_dt_wide) # Number of obs to resample

nboot = 2000
iron_boot <- matrix(0, ncol=12, nrow=nboot)

for (i in 1: nboot) {
set.seed(i)
iron.sample <- iron_dt_wide[sample(nrow(iron_dt_wide), size=dt.size, 
                                   replace=TRUE), ] # sampling with replacement
iron.sample$id <- seq.int(nrow(iron.sample)) # Add new unique id to sampled records

iron_dt_long <- reshape(iron.sample, idvar="id", timevar="visit", 
                        varying = c("mulout.0","mulout.6","mulout.12",
                                    "agefin.0", "agefin.6", "agefin.12",
                                    "hgb.0","hgb.6","hgb.12",
                                    "hgb_SNP_calib.0","hgb_SNP_calib.6","hgb_SNP_calib.12",
                                    "hgb_Norm_calib.0","hgb_Norm_calib.6","hgb_Norm_calib.12"),
                        direction = "long")


model.hgb <- lmer(hgb ~ as.factor(grouplog) + as.factor(visit) + as.factor(grouplog)*as.factor(visit)
                  +as.factor(child_educ)+ ( 1 | id), na.action=na.omit, data = iron_dt_long)


# When the mediator is uncalibrated 
model.outcome_error <- lmer(mulout ~ as.factor(grouplog) + hgb + as.factor(visit) + agefin +
                              as.factor(child_educ)+( 1 | id), na.action=na.omit, data = iron_dt_long)

ACME.error_BL <- fixef(model.hgb)[2] * fixef(model.outcome_error)[3]
ACME.error_6 <- (fixef(model.hgb)[2] + fixef(model.hgb)[6])* fixef(model.outcome_error)[3]
ACME.error_12 <- (fixef(model.hgb)[2] + fixef(model.hgb)[7])* fixef(model.outcome_error)[3]
ADE_error <- as.numeric(fixef(model.outcome_error)[2])

# When the mediator is calibrated assuming normality
model.outcome_norm <- lmer(mulout ~ as.factor(grouplog) + hgb_Norm_calib + as.factor(visit) + agefin +
                             as.factor(child_educ)+( 1 | id), na.action=na.omit, data = iron_dt_long)

ACME.norm_BL <- fixef(model.hgb)[2] * fixef(model.outcome_norm)[3]
ACME.norm_6 <- (fixef(model.hgb)[2] + fixef(model.hgb)[6])* fixef(model.outcome_norm)[3]
ACME.norm_12 <- (fixef(model.hgb)[2] + fixef(model.hgb)[7])* fixef(model.outcome_norm)[3]
ADE_norm <- as.numeric(fixef(model.outcome_norm)[2])

# When the mediator is calibrated assuming SNP
model.outcome_SNP <- lmer(mulout ~ as.factor(grouplog) + hgb_SNP_calib + as.factor(visit) + agefin +
                            as.factor(child_educ)+( 1 | id), na.action=na.omit, data = iron_dt_long)

ACME.SNP_BL <- fixef(model.hgb)[2] * fixef(model.outcome_SNP)[3]
ACME.SNP_6 <- (fixef(model.hgb)[2] + fixef(model.hgb)[6])* fixef(model.outcome_SNP)[3]
ACME.SNP_12 <- (fixef(model.hgb)[2] + fixef(model.hgb)[7])* fixef(model.outcome_SNP)[3]
ADE_SNP <- as.numeric(fixef(model.outcome_SNP)[2])

iron_boot[i, ] <- cbind(ACME.norm_BL, ACME.norm_6, ACME.norm_12, 
                        ACME.SNP_BL, ACME.SNP_6, ACME.SNP_12,
                        ACME.error_BL, ACME.error_6, ACME.error_12,
                        ADE_norm, ADE_SNP, ADE_error)
print(i)
}

# computing the percentile CI

library(matrixStats)

probs <- c(0.025, 0.975)
iron_boot_CI <- colQuantiles(iron_boot, probs=probs)

iron_results <- cbind(iron_est, iron_boot_CI)

EstBL <- rbind(iron_results[1,], iron_results[4,],iron_results[7,] )
Est6M <- rbind(iron_results[2,], iron_results[5,],iron_results[8,] )
Est12M <- rbind(iron_results[3,], iron_results[6,],iron_results[9,] )
EstADE <- iron_results[10:12,]

Table_iron <- rbind(EstBL, Est6M, Est12M, EstADE)
rownames(Table_iron) <- c("Normal calibration","SNP calibration","No calibration",
                          "Normal calibration6","SNP calibration6","No calibration6",
                          "Normal calibration12","SNP calibration12","No calibration12",
                          "Normal calibration ADE","SNP calibration ADE","No calibration ADE")

colnames(Table_iron) <- c("Estimate","LP","UP")

############## Printing table results to LATEX  ##################################

library(xtable)

print(xtable(Table_iron, digits=c(0,3,3,3), align="lccc",
             caption="ACME and ADE estimates for the malaria-cognition  
             data"), caption.placement = "top")

##################################################################################
####################################################################################

########### COMPUTING BOOTSTRAP CI FOR HIV-LIVE DATA #######################

hiv <- read.sas7bdat("C:/Users/John/Google Drive/PhD/Thesis/Paper 1/Data application/hiv_live_final.sas7bdat", debug=FALSE)
hiv2 <- hiv[,c(1,2,9,11,12,16,17,19,21,22,23)] # Keeping needed variables
hiv2[ is.na(hiv2) ] <- NA # Make NaN become NA

# Analysis to check interactions in the outcome model for the 
# iron-cognition data
model.hiv.nointer <- lme(sqrt_cd4 ~ niaahaz_BL + pct3d_p + age_BL + homeless + 
                           hiv_qol + as.factor(TP), random=~ 1 | ID, 
                         na.action=na.omit,data = hiv2)
summary(model.hiv.nointer)

summary(lme(sqrt_cd4 ~ niaahaz_BL + pct3d_p + age_BL + homeless + hiv_qol + 
       as.factor(TP)+niaahaz_BL*as.factor(TP), random=~ 1 | ID, 
       na.action=na.omit,data = hiv2))

model.hiv.TPcd4 <- lme(sqrt_cd4 ~ niaahaz_BL + pct3d_p + age_BL + homeless + 
                           hiv_qol + as.factor(TP)+pct3d_p*as.factor(TP), 
                         random=~ 1 | ID, na.action=na.omit,data = hiv2)
summary(model.hiv.TPcd4)

# Comment: It is not appropriate to compare two models with different
# fixed effects when REML is used. Why, because the RE are based on the 
# fixed effects, which if different, the models are not comparable.

summary(lme(sqrt_cd4 ~ niaahaz_BL + pct3d_p + age_BL + homeless + hiv_qol + 
              as.factor(TP)+niaahaz_BL*pct3d_p, random=~ 1 | ID, 
            na.action=na.omit,data = hiv2))

######## Checking identifiability assumption 4: No mediator-outcome confounder is
######## affected by the exposure

## Association between exposure and time-varying mediator-outcome confounder
## at different time points

TP = c(1, 3, 5, 6, 7, 8)

p.val_hivqol <- c(rep(0, 6))
p.val_homeless <- c(rep(0, 6))

for(k in 1:length(TP))
{
  hiv.qual_assump4 <- summary(glm(niaahaz_BL ~ hiv_qol,family=binomial(link='logit'),data=subset(hiv2, TP==TP[k])))$coefficients[2,4]
  p.val_hivqol[k] <- round( hiv.qual_assump4, digits = 3)
  
  hiv.homeless_assump4 <- summary(glm(niaahaz_BL ~ homeless,family=binomial(link='logit'),data=subset(hiv2, TP==TP[k])))$coefficients[2,4]
  p.val_homeless[k] <- round(hiv.homeless_assump4, digits = 3)
}


hiv.educ_assump4 <- glm(niaahaz_BL ~ as.factor(child_educ),family=binomial(link='logit'),data=subset(hiv2, visit==0))

summary(iron.age_assump4)
summary(iron.educ_assump4)

############# ACME and ADE estimates without resampling  #################

model.adh <- lmer(pct3d_p ~ niaahaz_BL + age_BL + hiv_qol + realm2 + as.factor(TP) + 
                    as.factor(TP)*niaahaz_BL + ( 1 | ID), na.action=na.omit,data = hiv2)

# Uncalibrated
model.cd4_error <- lmer(sqrt_cd4 ~ niaahaz_BL + pct3d_p + age_BL + homeless + hiv_qol + 
           as.factor(TP) + ( 1 | ID), na.action=na.omit,data = hiv2)

ACME.error_hiv_BL <- fixef(model.adh)[2] * fixef(model.cd4_error)[3]
ACME.error_hiv_2 <- (fixef(model.adh)[2] + fixef(model.adh)[13]) * fixef(model.cd4_error)[3]
ACME.error_hiv_3 <- (fixef(model.adh)[2] + fixef(model.adh)[14]) * fixef(model.cd4_error)[3]
ACME.error_hiv_4 <- (fixef(model.adh)[2] + fixef(model.adh)[15]) * fixef(model.cd4_error)[3]
ACME.error_hiv_5 <- (fixef(model.adh)[2] + fixef(model.adh)[16]) * fixef(model.cd4_error)[3]
ACME.error_hiv_6 <- (fixef(model.adh)[2] + fixef(model.adh)[17]) * fixef(model.cd4_error)[3]
ACME.error_hiv_7 <- (fixef(model.adh)[2] + fixef(model.adh)[18]) * fixef(model.cd4_error)[3]
ACME.error_hiv_8 <- (fixef(model.adh)[2] + fixef(model.adh)[19]) * fixef(model.cd4_error)[3]

ADE.error_hiv <- fixef(model.cd4_error)[2]

ACME.error_hiv <- rbind(ACME.error_hiv_BL, ACME.error_hiv_2, ACME.error_hiv_3,
                        ACME.error_hiv_4, ACME.error_hiv_5, ACME.error_hiv_6,
                        ACME.error_hiv_7, ACME.error_hiv_8)

# Calibrated assuming normality
model.cd4_norm <- lmer(sqrt_cd4 ~ niaahaz_BL + adh_Norm_calib + age_BL + homeless + hiv_qol + 
                          as.factor(TP) + ( 1 | ID), na.action=na.omit,data = hiv2)
ACME.norm_hiv_BL <- fixef(model.adh)[2] * fixef(model.cd4_norm)[3]
ACME.norm_hiv_2 <- (fixef(model.adh)[2] + fixef(model.adh)[13]) * fixef(model.cd4_norm)[3]
ACME.norm_hiv_3 <- (fixef(model.adh)[2] + fixef(model.adh)[14]) * fixef(model.cd4_norm)[3]
ACME.norm_hiv_4 <- (fixef(model.adh)[2] + fixef(model.adh)[15]) * fixef(model.cd4_norm)[3]
ACME.norm_hiv_5 <- (fixef(model.adh)[2] + fixef(model.adh)[16]) * fixef(model.cd4_norm)[3]
ACME.norm_hiv_6 <- (fixef(model.adh)[2] + fixef(model.adh)[17]) * fixef(model.cd4_norm)[3]
ACME.norm_hiv_7 <- (fixef(model.adh)[2] + fixef(model.adh)[18]) * fixef(model.cd4_norm)[3]
ACME.norm_hiv_8 <- (fixef(model.adh)[2] + fixef(model.adh)[19]) * fixef(model.cd4_norm)[3]

ADE.norm_hiv <- fixef(model.cd4_norm)[2]

ACME.norm_hiv <- rbind(ACME.norm_hiv_BL, ACME.norm_hiv_2, ACME.norm_hiv_3,
                        ACME.norm_hiv_4, ACME.norm_hiv_5, ACME.norm_hiv_6,
                        ACME.norm_hiv_7, ACME.norm_hiv_8)

# Calibrated assuming SNP
model.cd4_SNP <- lmer(sqrt_cd4 ~ niaahaz_BL + adh_SNP_calib + age_BL + homeless + hiv_qol + 
                          as.factor(TP) + ( 1 | ID), na.action=na.omit,data = hiv2)

ACME.SNP_hiv_BL <- fixef(model.adh)[2] * fixef(model.cd4_SNP)[3]
ACME.SNP_hiv_2 <- (fixef(model.adh)[2] + fixef(model.adh)[13]) * fixef(model.cd4_SNP)[3]
ACME.SNP_hiv_3 <- (fixef(model.adh)[2] + fixef(model.adh)[14]) * fixef(model.cd4_SNP)[3]
ACME.SNP_hiv_4 <- (fixef(model.adh)[2] + fixef(model.adh)[15]) * fixef(model.cd4_SNP)[3]
ACME.SNP_hiv_5 <- (fixef(model.adh)[2] + fixef(model.adh)[16]) * fixef(model.cd4_SNP)[3]
ACME.SNP_hiv_6 <- (fixef(model.adh)[2] + fixef(model.adh)[17]) * fixef(model.cd4_SNP)[3]
ACME.SNP_hiv_7 <- (fixef(model.adh)[2] + fixef(model.adh)[18]) * fixef(model.cd4_SNP)[3]
ACME.SNP_hiv_8 <- (fixef(model.adh)[2] + fixef(model.adh)[19]) * fixef(model.cd4_SNP)[3]

ADE.SNP_hiv <- fixef(model.cd4_SNP)[2]

ACME.SNP_hiv <- rbind(ACME.SNP_hiv_BL, ACME.SNP_hiv_2, ACME.SNP_hiv_3,
                       ACME.SNP_hiv_4, ACME.SNP_hiv_5, ACME.SNP_hiv_6,
                       ACME.SNP_hiv_7, ACME.SNP_hiv_8)

Table_ACME_hiv <- rbind(ACME.norm_hiv, ACME.SNP_hiv, ACME.error_hiv, 
                        ADE.norm_hiv, ADE.SNP_hiv, ADE.error_hiv)
# colnames(Table_ACME_hiv) <- c("error","normal","SNP")

############ Plot of HIV-LIVE ACME estimates ###############
# 
# x <- c(1:8)
# 
# plot(range(x), range(Table_ACME_hiv[,3]), type="n", xlab="Time points", 
#      ylab="ACME estimates", cex.axis=1.6,cex.lab=1.7 )
# 
# lines(x, Table_ACME_hiv[,3], type="b", lwd=4.5,lty=3, pch=8) 
# lines(x, Table_ACME_hiv[,2], type="b", lwd=4,lty=5, pch=6)
# lines(x, Table_ACME_hiv[,1], type="b", lwd=4,lty=1, pch=0)
# legend(4,-0.2, pch=c(0,6,8),lty=c(1,5,3), lwd=c(2,2,2.5), bty = "n",
#        legend = c("Uncalibrated","Normal calibration","SNP calibration"),
#        cex = 1.6)

############################################################################

hiv_wide <- reshape(hiv2, 
                        timevar = "TP",
                        idvar = c("ID", "niaahaz_BL","age_BL"),
                        direction = "wide")

n_hiv <- nrow(hiv_wide) # Number of obs to resample

nboot = 2000
hiv_boot <- matrix(0, ncol=27, nrow=nboot)

for (i in 1: nboot) {


hiv.sample <- hiv_wide[sample(nrow(hiv_wide), size=n_hiv, 
                                   replace=TRUE), ] # sampling with replacement
hiv.sample$studyid <- seq.int(nrow(hiv.sample)) # Add new unique id to sampled records

hiv_dt_long <- reshape(hiv.sample, idvar="studyid", timevar="TP", 
                        varying = c("sqrt_cd4.1","sqrt_cd4.2","sqrt_cd4.3","sqrt_cd4.4","sqrt_cd4.5","sqrt_cd4.6","sqrt_cd4.7","sqrt_cd4.8",
                                    "hiv_qol.1","hiv_qol.2","hiv_qol.3","hiv_qol.4","hiv_qol.5","hiv_qol.6","hiv_qol.7","hiv_qol.8",
                                    "pct3d_p.1","pct3d_p.2","pct3d_p.3","pct3d_p.4","pct3d_p.5","pct3d_p.6","pct3d_p.7","pct3d_p.8",
                                    "realm2.1","realm2.2","realm2.3","realm2.4","realm2.5","realm2.6","realm2.7","realm2.8",
                                    "adh_SNP_calib.1","adh_SNP_calib.2","adh_SNP_calib.3","adh_SNP_calib.4","adh_SNP_calib.5","adh_SNP_calib.6","adh_SNP_calib.7","adh_SNP_calib.8",
                                    "adh_Norm_calib.1","adh_Norm_calib.2","adh_Norm_calib.3","adh_Norm_calib.4","adh_Norm_calib.5","adh_Norm_calib.6","adh_Norm_calib.7","adh_Norm_calib.8",
                                    "homeless.1","homeless.2","homeless.3","homeless.4","homeless.5","homeless.6","homeless.7","homeless.8"),
                        direction = "long")

model.adh <- lmer(pct3d_p ~ niaahaz_BL + age_BL + hiv_qol + realm2 + as.factor(TP) + 
                    as.factor(TP)*niaahaz_BL + ( 1 | studyid), na.action=na.omit, data = hiv_dt_long )

# Uncalibrated
model.cd4_err <- lmer(sqrt_cd4 ~ niaahaz_BL + pct3d_p + age_BL + homeless + hiv_qol + 
                          as.factor(TP) + ( 1 | studyid), na.action=na.omit,data = hiv_dt_long)

ACME.err_hiv_BL <- fixef(model.adh)[2] * fixef(model.cd4_err)[3]
ACME.err_hiv_2 <- (fixef(model.adh)[2] + fixef(model.adh)[13]) * fixef(model.cd4_err)[3]
ACME.err_hiv_3 <- (fixef(model.adh)[2] + fixef(model.adh)[14]) * fixef(model.cd4_err)[3]
ACME.err_hiv_4 <- (fixef(model.adh)[2] + fixef(model.adh)[15]) * fixef(model.cd4_err)[3]
ACME.err_hiv_5 <- (fixef(model.adh)[2] + fixef(model.adh)[16]) * fixef(model.cd4_err)[3]
ACME.err_hiv_6 <- (fixef(model.adh)[2] + fixef(model.adh)[17]) * fixef(model.cd4_err)[3]
ACME.err_hiv_7 <- (fixef(model.adh)[2] + fixef(model.adh)[18]) * fixef(model.cd4_err)[3]
ACME.err_hiv_8 <- (fixef(model.adh)[2] + fixef(model.adh)[19]) * fixef(model.cd4_err)[3]

ADE.err_hiv <- fixef(model.cd4_err)[2]

# Calibrated assuming normality
model.cd4_norm <- lmer(sqrt_cd4 ~ niaahaz_BL + adh_Norm_calib + age_BL + homeless + hiv_qol + 
                         as.factor(TP) + ( 1 | studyid), na.action=na.omit,data = hiv_dt_long)

ACME.norm_hiv_BL <- fixef(model.adh)[2] * fixef(model.cd4_norm)[3]
ACME.norm_hiv_2 <- (fixef(model.adh)[2] + fixef(model.adh)[13]) * fixef(model.cd4_norm)[3]
ACME.norm_hiv_3 <- (fixef(model.adh)[2] + fixef(model.adh)[14]) * fixef(model.cd4_norm)[3]
ACME.norm_hiv_4 <- (fixef(model.adh)[2] + fixef(model.adh)[15]) * fixef(model.cd4_norm)[3]
ACME.norm_hiv_5 <- (fixef(model.adh)[2] + fixef(model.adh)[16]) * fixef(model.cd4_norm)[3]
ACME.norm_hiv_6 <- (fixef(model.adh)[2] + fixef(model.adh)[17]) * fixef(model.cd4_norm)[3]
ACME.norm_hiv_7 <- (fixef(model.adh)[2] + fixef(model.adh)[18]) * fixef(model.cd4_norm)[3]
ACME.norm_hiv_8 <- (fixef(model.adh)[2] + fixef(model.adh)[19]) * fixef(model.cd4_norm)[3]

ADE.norm_hiv2 <- fixef(model.cd4_norm)[2]

# Calibrated assuming SNP
model.cd4_SNP <- lmer(sqrt_cd4 ~ niaahaz_BL + adh_SNP_calib + age_BL + homeless + hiv_qol + 
                        as.factor(TP) + ( 1 | studyid), na.action=na.omit, data = hiv_dt_long)

ACME.SNP_hiv_BL <- fixef(model.adh)[2] * fixef(model.cd4_SNP)[3]
ACME.SNP_hiv_2 <- (fixef(model.adh)[2] + fixef(model.adh)[13]) * fixef(model.cd4_SNP)[3]
ACME.SNP_hiv_3 <- (fixef(model.adh)[2] + fixef(model.adh)[14]) * fixef(model.cd4_SNP)[3]
ACME.SNP_hiv_4 <- (fixef(model.adh)[2] + fixef(model.adh)[15]) * fixef(model.cd4_SNP)[3]
ACME.SNP_hiv_5 <- (fixef(model.adh)[2] + fixef(model.adh)[16]) * fixef(model.cd4_SNP)[3]
ACME.SNP_hiv_6 <- (fixef(model.adh)[2] + fixef(model.adh)[17]) * fixef(model.cd4_SNP)[3]
ACME.SNP_hiv_7 <- (fixef(model.adh)[2] + fixef(model.adh)[18]) * fixef(model.cd4_SNP)[3]
ACME.SNP_hiv_8 <- (fixef(model.adh)[2] + fixef(model.adh)[19]) * fixef(model.cd4_SNP)[3]

ADE.SNP_hiv2 <- fixef(model.cd4_SNP)[2]

hiv_boot[i, ] <- cbind(ACME.norm_hiv_BL, ACME.norm_hiv_2, ACME.norm_hiv_3,
                       ACME.norm_hiv_4, ACME.norm_hiv_5, ACME.norm_hiv_6,
                       ACME.norm_hiv_7, ACME.norm_hiv_8,
                       ACME.SNP_hiv_BL, ACME.SNP_hiv_2, ACME.SNP_hiv_3,
                       ACME.SNP_hiv_4, ACME.SNP_hiv_5, ACME.SNP_hiv_6,
                       ACME.SNP_hiv_7, ACME.SNP_hiv_8,
                       ACME.err_hiv_BL, ACME.err_hiv_2, ACME.err_hiv_3,
                       ACME.err_hiv_4, ACME.err_hiv_5, ACME.err_hiv_6,
                       ACME.err_hiv_7, ACME.err_hiv_8,
                       ADE.norm_hiv2, ADE.SNP_hiv2, ADE.err_hiv)
print(i)
}


hiv_boot_CI <- colQuantiles(hiv_boot, probs=probs)
# hiv_boot_CI

hiv_results <- cbind(Table_ACME_hiv, hiv_boot_CI)
hiv_results
# colnames(hiv_results) <- c("Estimate", "LP", "UP")
# hiv_results <- as.data.frame(hiv_results)

####### Creating plot for HIV-LIVE ACME estimates #########
library(Hmisc)


m <- c(1:8)
my <- c(-1.5:1.5)
plot(range(m),range(my), type="n", xlab = "Time points", 
     ylab = "ACME estimates",cex.axis=1.3,cex.lab=1.5 )



errbar(m, hiv_results[1:8, 1], hiv_results[1:8, 3], hiv_results[1:8, 2], 
       add=T, pch=1, cap=.03, lwd=3.5,lty = 3)
m1 = m+0.1
errbar(m1, hiv_results[9:16, 1], hiv_results[9:16, 3],hiv_results[9:16, 2], 
       add=T, pch=0, cap=.02, lwd=4,lty = 5)
m2 = m-0.1
errbar(m2, hiv_results[17:24, 1], hiv_results[17:24, 3],hiv_results[17:24, 2], 
       add=T, pch=16, cap=.02, lwd=3,lty = 1)
abline(h=0, lty = 2, lwd = 2, col = 'grey')
legend("top", pch=c(1,0,16),lty=c(3,5,1), lwd=c(2.5,2.5,2.5), bty = "n",
       legend = c("Normal calibration","SNP calibration", "Uncalibrated"),
       cex = 1.3)

########################################################


Table_hiv <- cbind(hiv_results[1:3, ], hiv_results[4:6, ])

rownames(Table_hiv) <- c("1Normal calibration","1SNP calibration", "1No calibration")


############################################################################

Table_hiv_iron <- rbind(Table_iron, Table_hiv)
colnames(Table_hiv_iron) <- c("Estimate","LCI","UCI", "Estimate","LCI","UCI")

############## Printing table results to LATEX  ##################################

library(xtable)

print(xtable(Table_hiv_iron, digits=c(0,4,4,4,4,4,4), align="lcccccc",
             caption="ACME and ADE estimates for the malaria-cognition data and HIV-LIVE 
             data"), caption.placement = "top")

##################################################################################