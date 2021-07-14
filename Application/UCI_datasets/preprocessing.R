######################## Contraceptive Method Choice 
datCMC <- read.csv("cmc.data.csv", header = FALSE)

head(datCMC)
colnames(datCMC) <- c("wife_age","wife_education","husband_education","nr_children",
                      "wife_religion","wife_is_working","husband_occupation",
                      "standard_of_living","media_exposure","contraceptive_method")

# wife's age: 0 (-Inf,20); 1 [20,30); 2 [30,40); 3 [40,Inf)
datCMC$wife_age <- cut(datCMC$wife_age,breaks = c(min(datCMC$wife_age) - 1,20,30,40,max(datCMC$wife_age) + 1),labels = FALSE) - 1

# number of children: 0 [0,1]; 1 [2,3]; 2 [4,5]; 3 [6,Inf)
datCMC$nr_children <- cut(datCMC$nr_children, breaks = c(-1,2,4,6,16), labels = FALSE) - 1

# wife is working: 0 [not working]; 1 [working]
datCMC$wife_is_working <- 1 - datCMC$wife_is_working

# reduce levels for wife's education, husband's education, husband's occupation, standard-of-living index
datCMC$wife_education <- datCMC$wife_education - 1
datCMC$husband_education <- datCMC$husband_education - 1
datCMC$husband_occupation <- datCMC$husband_occupation - 1
datCMC$standard_of_living <- datCMC$standard_of_living - 1

# media exposure: 0 [not good]; 1 [good]
datCMC$media_exposure <- 1 - datCMC$media_exposure

# contraceptive method: 0 [no use]; 1 [short-term]; 2 [long-term]
datCMC$contraceptive_method[datCMC$contraceptive_method == 1] <- 0
datCMC$contraceptive_method[datCMC$contraceptive_method == 3] <- 1



######################## Congressional Voting Records 
datVotes <- read.csv("house-votes-84.data.csv", header = FALSE)

head(datVotes)
colnames(datVotes) <- c("party","handicapped_infants","water_project_cost_sharing",
                        "adoption_of_the_budget_resolution","physician_fee_freeze","el_salvador_aid",
                        "religious_groups_in_schools","anti_satellite_test_ban","aid_to_nicaraguan_contras",
                        "mx_missile","immigration","synfuels_corporation_cutback","education_spending",
                        "superfund_right_to_sue","crime","duty_free_exports","export_administration_act_south_africa")

# party: republican [0]; democratic [1]
datVotes$party[datVotes$party == "republican"] <- 0
datVotes$party[datVotes$party == "democrat"] <- 1

# replace "n" with 0 and "y" with 1, impute "?"
for (i in c(2:ncol(datVotes))) {
  p_n <- sum(datVotes[,i] == "n")
  p_y <- sum(datVotes[,i] == "y")
  datVotes[datVotes[,i] == "?",i] <- sample(c("n","y"), sum(datVotes[,i] == "?"), replace = TRUE, prob = c(p_n, p_y) / (p_n + p_y))
  datVotes[datVotes[,i] == "n",i] <- 0
  datVotes[datVotes[,i] == "y",i] <- 1
}

for (i in c(1:ncol(datVotes))) {
  datVotes[,i] <- as.integer(datVotes[,i])
}



######################## Primary Tumor 
datTumor <- read.csv("primary-tumor.data.csv", header = FALSE)

head(datTumor)
datTumor$V1 <- NULL
colnames(datTumor) <- c("age","sex","histologic_type","degree_of_diffe","bone","bone_marrow",
                        "lung","pleura","peritoneum","liver","brain","skin","neck","supraclavicular",
                        "axillar","mediastinum","abdominal")
datTumor <- datTumor[-which(datTumor$V34 == "?"),]
# impute missing data
datTumor$sex[datTumor$sex == "?"] <- sample(names(table(datTumor$sex))[-1],
                                            sum(datTumor$sex == "?"),
                                            replace = TRUE,
                                            prob = table(datTumor$sex)[-1] / sum(table(datTumor$sex)[-1]))
datTumor$histologic_type[datTumor$histologic_type == "?"] <- sample(names(table(datTumor$histologic_type))[-1],
                                                                    sum(datTumor$histologic_type == "?"),
                                                                    replace = TRUE,
                                                                    prob = table(datTumor$histologic_type)[-1] / sum(table(datTumor$histologic_type)[-1]))
datTumor$degree_of_diffe[datTumor$degree_of_diffe == "?"] <- sample(names(table(datTumor$degree_of_diffe))[-1],
                                                                    sum(datTumor$degree_of_diffe == "?"),
                                                                    replace = TRUE,
                                                                    prob = table(datTumor$degree_of_diffe)[-1] / sum(table(datTumor$degree_of_diffe)[-1]))
datTumor$skin[datTumor$skin == "?"] <- sample(names(table(datTumor$skin))[-1],
                                              sum(datTumor$skin == "?"),
                                              replace = TRUE,
                                              prob = table(datTumor$skin)[-1] / sum(table(datTumor$skin)[-1]))
datTumor$axillar[datTumor$axillar == "?"] <- sample(names(table(datTumor$axillar))[-1],
                                                    sum(datTumor$axillar == "?"),
                                                    replace = TRUE,
                                                    prob = table(datTumor$axillar)[-1] / sum(table(datTumor$axillar)[-1]))

datTumor <- datTumor - 1
for (i in c(1:ncol(datTumor))) {
  datTumor[,i] <- as.integer(datTumor[,i])
}



######################## SPECT Heart 
datSPECT <- read.csv("SPECT.csv", header = FALSE)

head(datSPECT)

for (i in c(1:ncol(datSPECT))) {
  datSPECT[,i] <- as.integer(datSPECT[,i])
}