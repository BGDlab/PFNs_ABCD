
library(psych)
library(qgraph)

x <- read.csv("data_ABCD_clinical.csv")

past <- x[,grepl("_P",colnames(x))]
current <- x[,grepl("_C",colnames(x))]

x1 <- data.frame(
past$	ActiveAvoidanceOfPhobicObject_P	,
past$	Anhedonia_P	,
past$	BingeEating_P	,
past$	Compulsions_P	,
past$	DecreasedNeedForSleep_P	,
past$	DepressedMood_P	,
past$	DifficultyRemainingSeated_P	,
past$	DifficultyRemainingSeatedForMoreThanOneSchoolYear_P	,
past$	DifficultySustainingAttention_P	,
past$	DifficultySustainingAttentionForMoreThanOneSchoolYear_P	,
past$	DistressAtInternalRemindersOfTrauma_P	,
past$	DistressDueToFearOrAvoidanceOfPhobicObject_P	,
past$	DistressUponSeparationFromHomeAttachmentFigures_P	,
past$	NegativeImpactRoleObligationsMissWorkDrugUse_P	,
past$	DurationAtLeast6MonthsODD_P	,
past$	NegativeImpactRoleObligationsMissWorkDrugUse_P	,
past$	NegativeImpactRoleObligationsMissWorkDrugUse_P	,
past$	EasilyDistracted_P	,
past$	EasilyDistractedForMoreThanOneSchoolYear_P	,
past$	EffortsToAvoidsThoughtsOfTrauma_P	,
past$	ElevatedMood_P	,
past$	Emaciation_P	,
past$	ExcessiveWorriesMoreDaysThanNot_P	,
past$	ExplosiveIrritability_P	,
past$	FearOfBecomingObese_P	,
past$	FearOfSocialSituations_P	,
past$	Frequencydruguse_P	,
past$	NegativeImpactRoleObligationsMissWorkDrugUse_P	,
past$	FrequencyOfDrinks_P	,
past$	Hallucinations_P	,
past$	NegativeImpactRoleObligationsMissWorkDrugUse_P	,
past$	HomicidalIdeation_P	,
past$	HomicidePlanning_P	,
past$	Hypersexuality_P	,
past$	ImpairmentDueToFearOrAvoidanceOfPhobicObject_P	,
past$	Impulsivity_P	,
past$	ImpulsivityForMoreThanOneSchoolYear_P	,
past$	Insomnia_P	,
past$	Irritability_P	,
past$	MarkedFearOfPhobicObject_P	,
past$	NegativeImpactInterpersonalProblemsAlcohol_P	,
past$	NegativeImpactInterpersonalProblemsDrugUse_P	,
past$	NegativeImpactRoleObligationsMissWorkDrugUse_P	,
past$	Nightmares_P	,
past$	Obsessions_P	,
past$	OftenArguesWithAdultsAuthority_P	,
past$	OftenBulliesOthers_P	,
past$	OftenDisobeysRulesRequests_P	,
past$	OftenHasThreeOrMoreDrinksADay_P	,
past$	OftenInitiatesPhysicalFights_P	,
past$	OftenLies_P	,
past$	OftenLosesTemper_P	,
past$	PanicAttacks_P	,
past$	PersecutoryDelusions_P	,
past$	PhobicObjectSpecialCase_P	,
past$	PoorEyeContact_P	,
past$	SchoolReluctanceRefusal_P	,
past$	Stealing_P	,
past$	StrictRoutines_P	,
past$	Truancy_P	,
past$	UnusualBodyMovements_P	,
past$	WeightControlOtherLaxativesExerciseDietingPills_P	,
past$	WeightControlVomiting_P	,
past$	NegativeImpactRoleObligationsMissWorkDrugUse_P	,
past$	Anhedonia_P_y	,
past$	DecreasedNeedForSleep_P_y	,
past$	DepressedMood_P_y	,
past$	NegativeImpactRoleObligationsMissWorkDrugUse_P	,
past$	ElevatedMood_P_y	,
past$	ExcessiveWorriesMoreDaysThanNot_P_y	,
past$	ExplosiveIrritability_P_y	,
past$	FearOfSocialSituations_P_y	,
past$	Hypersexuality_P_y	,
past$	Insomnia_P_y	,
past$	Irritability_P_y	,
past$	NegativeImpactRoleObligationsMissWorkDrugUse_P	,
past$	SelfInjuriousBehavior_P_y	,
past$	WishesBetterOffDead_P_y	,
past$	SuicidalAttempt_P_y	,
past$	SuicidalIdeation_P_y	)

x2 <- data.frame(
current$	ActiveAvoidanceOfPhobicObject_C	,
current$	Anhedonia_C	,
current$	BingeEating_C	,
current$	Compulsions_C	,
current$	DecreasedNeedForSleep_C	,
current$	DepressedMood_C	,
current$	DifficultyRemainingSeated_C	,
current$	DifficultyRemainingSeatedSinceElementarySchool_C	,
current$	DifficultySustainingAttention_C	,
current$	DifficultySustainingAttentionSinceElementarySchool_C	,
current$	DistressAtInternalRemindersOfTrauma_C	,
current$	DistressDueToFearOrAvoidanceOfPhobicObject_C	,
current$	DistressUponSeparationFromHomeAttachmentFigures_C	,
current$	DrugsTried_C	,
current$	DurationAtLeast6Months_C	,
current$	DurationOfPastPhobiaAtLeast6Months_C	,
current$	DurationOfPhobiaAtLeast6Months_C	,
current$	EasilyDistracted_C	,
current$	EasilyDistractedSinceElementarySchool_C	,
current$	EffortsToAvoidsThoughtsOfTrauma_C	,
current$	ElevatedMood_C	,
current$	Emaciation_C	,
current$	ExcessiveWorriesMoreDaysThanNot_C	,
current$	ExplosiveIrritability_C	,
current$	FearOfBecomingObese_C	,
current$	FearOfSocialSituations_C	,
current$	Frequencydruguse_C	,
current$	FrequencyDrugUseSpecial_C	,
current$	NegativeImpactInterpersonalProblemsDrugUse_C	,
current$	Hallucinations_C	,
current$	HistoryOfTraumaticEvent_C	,
current$	HomicidalIdeation_C	,
current$	HomicidePlanning_C	,
current$	Hypersexuality_C	,
current$	ImpairmentDueToFearOrAvoidanceOfPhobicObject_C	,
current$	Impulsivity_C	,
current$	NegativeImpactInterpersonalProblemsDrugUse_C	,
current$	Insomnia_C	,
current$	Irritability_C	,
current$	MarkedFearOfPhobicObject_C	,
current$	NegativeImpactInterpersonalProblemsAlcohol_C	,
current$	NegativeImpactInterpersonalProblemsDrugUse_C	,
current$	NegativeImpactRoleObligationsMissWorkDrugUse_C	,
current$	Nightmares_C	,
current$	Obsessions_C	,
current$	OftenArguesWithAdultsAuthority_C	,
current$	OftenBulliesOthers_C	,
current$	OftenDisobeysRulesRequests_C	,
current$	OftenHasThreeOrMoreDrinksADay_C	,
current$	OftenInitiatesPhysicalFights_C	,
current$	OftenLies_C	,
current$	OftenLosesTemper_C	,
current$	PanicAttacks_C	,
current$	PersecutoryDelusions_C	,
current$	NegativeImpactInterpersonalProblemsDrugUse_C	,
current$	PoorEyeContact_C	,
current$	SchoolReluctanceRefusal_C	,
current$	Stealing_C	,
current$	StrictRoutines_C	,
current$	Truancy_C	,
current$	UnusualBodyMovements_C	,
current$	WeightControlOtherLaxativesExerciseDietingPills_C	,
current$	WeightControlVomiting_C	,
current$	WorryingHasLastedAtLeast6Months_C	,
current$	Anhedonia_C_y	,
current$	DecreasedNeedForSleep_C_y	,
current$	DepressedMood_C_y	,
current$	DurationAtLeast6Months_C_y	,
current$	ElevatedMood_C_y	,
current$	ExcessiveWorriesMoreDaysThanNot_C_y	,
current$	ExplosiveIrritability_C_y	,
current$	FearOfSocialSituations_C_y	,
current$	Hypersexuality_C_y	,
current$	Insomnia_C_y	,
current$	Irritability_C_y	,
current$	WorryingHasLastedAtLeast6Months_C_y	,
current$	SelfInjuriousBehavior_C_y	,
current$	WishesBetterOffDead_C_y	,
current$	SuicidalAttempt_C_y	,
current$	SuicidalIdeation_C_y	)

x55 <- x1 + x2

colnames(x55) <- c(
"	ActiveAvoidanceOfPhobicObject	",
"	Anhedonia	",
"	BingeEating	",
"	Compulsions	",
"	DecreasedNeedForSleep	",
"	DepressedMood	",
"	DifficultyRemainingSeated	",
"	DifficultyRemainingSeatedForMoreThanOneSchoolYear	",
"	DifficultySustainingAttention	",
"	DifficultySustainingAttentionForMoreThanOneSchoolYear	",
"	DistressAtInternalRemindersOfTrauma	",
"	DistressDueToFearOrAvoidanceOfPhobicObject	",
"	DistressUponSeparationFromHomeAttachmentFigures	",
"	DrugsTried	",
"	DurationAtLeast6MonthsODD	",
"	DurationOfPastPhobiaAtLeast6Months	",
"	DurationOfPhobiaAtLeast6Months	",
"	EasilyDistracted	",
"	EasilyDistractedForMoreThanOneSchoolYear	",
"	EffortsToAvoidsThoughtsOfTrauma	",
"	ElevatedMood	",
"	Emaciation	",
"	ExcessiveWorriesMoreDaysThanNot	",
"	ExplosiveIrritability	",
"	FearOfBecomingObese	",
"	FearOfSocialSituations	",
"	Frequencydruguse	",
"	FrequencyDrugUseSpecial	",
"	FrequencyOfDrinks	",
"	Hallucinations	",
"	HistoryOfTraumaticEvent	",
"	HomicidalIdeation	",
"	HomicidePlanning	",
"	Hypersexuality	",
"	ImpairmentDueToFearOrAvoidanceOfPhobicObject	",
"	Impulsivity	",
"	ImpulsivityForMoreThanOneSchoolYear	",
"	Insomnia	",
"	Irritability	",
"	MarkedFearOfPhobicObject	",
"	NegativeImpactInterpersonalProblemsAlcohol	",
"	NegativeImpactInterpersonalProblemsDrugUse	",
"	NegativeImpactRoleObligationsMissWorkDrugUse	",
"	Nightmares	",
"	Obsessions	",
"	OftenArguesWithAdultsAuthority	",
"	OftenBulliesOthers	",
"	OftenDisobeysRulesRequests	",
"	OftenHasThreeOrMoreDrinksADay	",
"	OftenInitiatesPhysicalFights	",
"	OftenLies	",
"	OftenLosesTemper	",
"	PanicAttacks	",
"	PersecutoryDelusions	",
"	PhobicObjectSpecialCase	",
"	PoorEyeContact	",
"	SchoolReluctanceRefusal	",
"	Stealing	",
"	StrictRoutines	",
"	Truancy	",
"	UnusualBodyMovements	",
"	WeightControlOtherLaxativesExerciseDietingPills	",
"	WeightControlVomiting	",
"	WorryingHasLastedAtLeast6Months	",
"	Anhedonia_y	",
"	DecreasedNeedForSleep_y	",
"	DepressedMood_y	",
"	DurationAtLeast6Months_y	",
"	ElevatedMood_y	",
"	ExcessiveWorriesMoreDaysThanNot_y	",
"	ExplosiveIrritability_y	",
"	FearOfSocialSituations_y	",
"	Hypersexuality_y	",
"	Insomnia_y	",
"	Irritability_y	",
"	WorryingHasLastedAtLeast6Months_y	",
"	SelfInjuriousBehavior_y	",
"	WishesBetterOffDead_y	",
"	SuicidalAttempt_y	",
"	SuicidalIdeation_y	")

x55 <- x55[,describe(x55)$mean > 0.005]

phobic_object <- round((x55$ActiveAvoidanceOfPhobicObject + x55$MarkedFearOfPhobicObject + x55$DistressDueToFearOrAvoidanceOfPhobicObject)/3,0)
impulsivity <- round((x55$Impulsivity + x55$ImpulsivityForMoreThanOneSchoolYear)/2,0)
diff_remain_seated <- round((x55$DifficultyRemainingSeated + x55$DifficultyRemainingSeatedForMoreThanOneSchoolYear)/2,0)
diff_sustain_attention <- round((x55$DifficultySustainingAttention + x55$DifficultySustainingAttentionForMoreThanOneSchoolYear)/2,0)
easily_distracted <- round((x55$EasilyDistracted + x55$EasilyDistractedForMoreThanOneSchoolYear)/2,0)
drug_use <- round((x55$Frequencydruguse + x55$DrugsTried)/2,0)

x55 <- data.frame(x55[,-c(1,7:10,12,14,16:19,27,33,34,29,37,48,61)],phobic_object,impulsivity,diff_remain_seated ,diff_sustain_attention,easily_distracted,drug_use)

x <- data.frame(x[,1:154],x55)

temp <- x[,5:214]
temp[temp>0] <- 1
x[,5:214] <- temp

pfid <- read.csv("p_factor_IDs.csv")

x <- merge(pfid,x,by=1,all=TRUE)

#x <- x[,c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,describe(x[,7:216])$mean > 0.005)]
x <- x[,c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,describe(x[,7:216])$mean > 0.02)]
temp <- x$eventname
temp[temp=="baseline_year_1_arm_1"] <- 1
temp[temp=="2_year_follow_up_y_arm_1"] <- 2
x$eventname <- temp

set.seed(2022)
ids <- data.frame(unique(x[,1]),sample(1:length(unique(x[,1])),length(unique(x[,1])),replace=FALSE))
train_ids <- ids[which(ids$sample.1.length.unique.x...1.....length.unique.x...1.....replace...FALSE. < 5939),] 
test_ids <- ids[which(ids$sample.1.length.unique.x...1.....length.unique.x...1.....replace...FALSE. > 5938),]

train <- x[x$src_subject_id %in% train_ids[,1],]
test <- x[x$src_subject_id %in% test_ids[,1],]
id_train <- train[,1]
id_test <- test[,1]
train[,1] <- rank(train[,1],ties.method="min")
test[,1] <- rank(test[,1],ties.method="min")

# write.csv(train,"ABCD_psychopathology_training_MPLUS_9July2022.csv",na="99999",row.names=FALSE)
# write.csv(test,"ABCD_psychopathology_testing_MPLUS_9July2022.csv",na="99999",row.names=FALSE)

train <- train[,-c(1:6,150,172,171,185,188)]































# pull in Mplus stuff and calculate scores

x1 <- read.csv("ABCD_psychopathology_full_sample.csv")
x2 <- read.csv("ABCD_CFA_loadings_for_scores.csv")
x2[is.na(x2)] <- 0

vars <- x2[,1]

x <- x1[,colnames(x1) %in% vars]
x <- x[,order(names(x))]
x2 <- x2[order(x2[,1]),]

sc <- as.matrix(x) %*% as.matrix(x2[,3:11])


# pull in Mplus scores and manipulate

x <- read.table("ABCD_psychopathology_scores_8f_FULL_SAMPLE.dat")
x2 <- read.csv("ABCD_MPLUS_full_sample.csv", header=FALSE)
colnames(x2)[1:6] <- c("src_subject_id","rel_family_id","site_id_l_br","interview_age","sex_M","eventname")

sc <- x[,c(126,128,130,132,134,136,138,140,142,144:146)]
colnames(sc) <- c("Factor1","Factor2","Factor3","Factor4","Factor5","Factor6","Factor7","Factor8","General_p","event","site","family")
sc[,1:9] <- scale(sc[,1:9])

x <- data.frame(x2[,1:6],sc)

# write.csv(x[,1:15],"ABCD_longitudinal_psychopathology_factors_scores_full_sample.csv",na="",row.names=FALSE)


