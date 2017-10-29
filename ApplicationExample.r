library(lme4)
library(dplyr)
library(lattice)


#--Read data from other project
load("data/.RData") 
names(df)


#--Examine if nested or fully crossed between factors
#Tabulate the factors
xtabs(~ respondent + title, df_HLM)


#--Model training
df_HLM <- select(df, preference, title, respondent, matches("^combined.*ct$"), matches("^tste_3.*ct"))

as.formula


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


#--Model performance
#Chisquare test
anova(model_HLMa, model_HLMb, model_HLMc, model_HLM, model_lm)

summary(model_lm)


#--Parameter significance
'
Depending on the method specified, confint() computes confidence intervals by:
"profile", "Wald", "boot"
help("pvalues") formore
'
confint(model_HLMa, level=0.95, method="Wald")

#If the diagnal line is streight, the profile interval can be used
#The vertical lines in the panels delimit the 50%, 80%, 90%, 95% and 99% confidence intervals
profile_HLMa <- profile(model_HLMa)
xyplot(profile_HLMa, absVal=FALSE)
confint(model_HLMa, level=0.95, method="profile")


#--Random effects
ranef(model_HLMa)
dotplot(ranef(model_HLMa, condVar=TRUE))
qqmath(ranef(model_HLMa, condVar=TRUE))
