library(lme4)
library(dplyr)
library(lattice)


#--Read data from other project
load("data/.RData") 
names(df)

#--
df_HLM <- select(df, preference, title, respondent, matches("^combined.*ct$"), matches("^tste_3.*ct"))

model_HLM <- lmer(
  data=df_HLM,
  REML=FALSE,
  formula=preference ~ 1 + (combined_autonomy_ct + combined_competence_ct + combined_relatedness_ct) * (1 + tste_3_0_ct + tste_3_1_ct + tste_3_2_ct) +
                    (1 + combined_autonomy_ct + combined_competence_ct + combined_relatedness_ct + tste_3_0_ct + tste_3_1_ct + tste_3_2_ct | respondent)
)

model_HLMa <- lmer(
  data=df_HLM,
  REML=FALSE,
  formula=preference ~ 1 + (combined_autonomy_ct + combined_competence_ct + combined_relatedness_ct) * (1 + tste_3_0_ct + tste_3_1_ct + tste_3_2_ct) +
    (1 | respondent)
)

model_HLMb <- lmer(
  data=df_HLM,
  REML=FALSE,
  formula=preference ~ 1 + (combined_autonomy_ct + combined_competence_ct + combined_relatedness_ct) * (1 + tste_3_0_ct + tste_3_1_ct + tste_3_2_ct) +
    (0 + combined_autonomy_ct + combined_competence_ct + combined_relatedness_ct + tste_3_0_ct + tste_3_1_ct + tste_3_2_ct | respondent)
)

model_HLMc <- lmer(
  data=df_HLM,
  REML=FALSE,
  formula=preference ~ 1 + (combined_autonomy_ct + combined_competence_ct + combined_relatedness_ct) * (1 + tste_3_0_ct + tste_3_1_ct + tste_3_2_ct) +
    (1 | respondent) + (1 | title)
)

model_lm <- lm(data=df_HLM, formula=preference ~ 1 + (combined_autonomy_ct + combined_competence_ct + combined_relatedness_ct) * (1 + tste_3_0_ct + tste_3_1_ct + tste_3_2_ct))
summary(model_lm)
BIC(model_lm)
AIC(model_lm)

#Chisquare test
anova(model_HLMa, model_HLMb, model_HLMc, model_HLM, model_lm)

summary(model_lm)


#--model fit sensitivity to parameter change
#Parameter profile deviance confidence interval
profile_HLMa <- profile(model_HLMa)

#The vertical lines in the panels delimit the 50%, 80%, 90%, 95% and 99% confidence intervals
xyplot(profile_HLMa, absVal=FALSE)
confint(profile_HLMa, level=0.95)

ranef(model_HLMa)
dotplot(ranef(model_HLMa, condVar=TRUE))
qqmath(ranef(model_HLMa, condVar=TRUE))


#Tabulate the factors
xtabs(~ respondent + title, df_HLM)
