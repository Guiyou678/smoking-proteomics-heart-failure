library(car)
library(haven)
library(skimr)
library(tidyverse)
library(nephro)
library(compareGroups)
library(survival)
library(glmnet)
library(survminer)
library(ggrepel)
library(gprofiler2)
library(writexl)
library(STRINGdb)
library(networkD3)
library(corrplot)
library(readxl)
library(forestplot)
library(ggVennDiagram)
library(patchwork)
############################protein annotation file######################################
folder_path <- ""
# Get all Olink CSV files
file_list <- list.files(
  path = folder_path,
  pattern = "^Olink Target 96.*\\.csv$",
  full.names = TRUE
)
# Read and combine all files (semicolon-delimited)
protein_annotation <- file_list %>%
  lapply(function(x) read_delim(x, delim = ";")) %>%
  bind_rows()
# Inspect
glimpse(protein_annotation)
# Save combined file
write_xlsx(
  protein_annotation,
  path = file.path(folder_path, "protein annotation.xlsx")
)

######################################read the dataset for analysis########################
#load sav file
BIOSTAT <- read_sav("BIOSTAT_index.sav")

# Creates a key and annotation file checking what is exactly in the columns
column_keys <- sapply(BIOSTAT, function(x) attributes(x)[["label"]])
COLUMN_ANNOTATION <- data.frame(COLNAME = names(column_keys), ANNOTATION = unlist(unname(column_keys)))

# Calculating the age based on visit and date of birth
# We use 365.25 because of leap year;)
BIOSTAT <- BIOSTAT %>%
  mutate(AGE = round(as.numeric(difftime(as.Date(VIS_DTE11), as.Date(DEM_DBIR), units = "days"))/365.25,0)) %>%
  mutate(SMOKING = case_when(CVRF_SMOKE == 1 ~ "Past/Never",
                             CVRF_SMOKE == 2 ~ "Past/Never",
                             CVRF_SMOKE == 3 ~ "Current"))%>%
  mutate(SMOKING3 = factor(case_when(CVRF_SMOKE == 1 ~ "Never",
                             CVRF_SMOKE == 2 ~ "Past",
                             CVRF_SMOKE == 3 ~ "Current"))) %>% 
  mutate(BMI = as.numeric(PE_WGT1/((PE_HGT1*0.01)^2)))

# Labeling the primary CV disease and create a new variable: ETIO
BIOSTAT <- BIOSTAT %>% 
  mutate(MHHF_DIS1 = factor(MHHF_DIS1, labels = c("Primary", "Contributory","Not Present", "Unknown")),
         MHHF_DIS2 = factor(MHHF_DIS2, labels = c("Primary", "Contributory","Not Present", "Unknown")),
         MHHF_DIS3 = factor(MHHF_DIS3, labels = c("Primary", "Contributory","Not Present", "Unknown")),
         MHHF_DIS4 = factor(MHHF_DIS4, labels = c("Primary", "Contributory","Not Present", "Unknown")),
         MHHF_DIS5 = factor(MHHF_DIS5, labels = c("Primary", "Contributory","Not Present", "Unknown")),
         MHHF_DIS6 = factor(MHHF_DIS6, labels = c("Primary", "Contributory","Not Present", "Unknown"))
  ) %>%
  mutate(ETIO = as.factor(ifelse(MHHF_DIS1 == "Primary", "Ischemic Heart Disease",
                                 ifelse(MHHF_DIS2 == "Primary", "Hypertensive",
                                        ifelse(MHHF_DIS3 == "Primary", "Cardiomyopathy",
                                               ifelse(MHHF_DIS4 == "Primary", "Valvular Disease",
                                                      ifelse(MHHF_DIS6 == "Primary", "Other", "Unknown"))))))) %>%
  mutate(ETIO = as.character(ETIO)) %>%
  mutate(ETIO = ifelse(is.na(ETIO), "Unknown", ETIO))


#calculate the kidney function using CKD epi 2009 formula and create a new variable: EGFRBL
BIOSTAT <- BIOSTAT %>%
  mutate(RACE = factor(DEM_RACE, levels = 1:6, labels = c("Caucasian", 
                                                          "Black", 
                                                          "Asian", 
                                                          "Oriental", 
                                                          "Hispanic", 
                                                          "Other"))) %>%
  mutate(SCREA1_mgDl = SCREA1 / 88.42) %>%
  
  mutate(SEX_creat = ifelse(DEM_SEX == 1, 1,0)) %>%
  mutate(
    RACEcreat = case_when(
      RACE == "Black" ~ 1,
      RACE %in% c("Caucasian", "Asian", "Oriental", "Hispanic", "Other") ~ 0,
      TRUE ~ NA_real_)) %>%
  mutate(EGFRBL = CKDEpi.creat(creatinine = SCREA1_mgDl, sex = SEX_creat, age = AGE, ethnicity = RACEcreat))


# Column names needed to compile baseline characteristics table
# You can find the column names and what they mean in the COLUMN_ANNOTATION variable
# Demographics column names
DEMOGRAPHICS <- c(
  "AGE",
  "DEM_SEX",
  "RACE",
  "BMI",
  "MH.factor",
  "M.factor",
  "H.factor"
)

# Heart failure characteristics
HF_CHARS <- c(
  "SYMP_NYHA1", #NYHA class (Visit 1 (1))
  "ECHO_LVEF1" #
)

# History of diseases
HISTORY <- c(
  "CVRF_NY8",
  "CVRF_NY10",
  "CVRF_NY1",
  "CVRF_NY6",
  "CVRF_NY14"
)

#lifestyle
LIFESTYLE <- c("CVRF_NY9")#alcohol

# Signs and symptoms 
SYMPTOMS_SIGNS_HF <- c(
  "PE_NY21",
  "PE_PULCON1",#oedema
  "PE_NY41" #elevated JVP
)

# Physical parameters
PHYSICAL_EXAM <- c(
  "PE_HGT1",
  "PE_WGT1",
  "PE_HR1",
  "PE_SBP1",
  "PE_DBP1"
)

# Relevant medication use
MEDICATIONS <- c(
  "ACE_BS",
  "BB_BS",
  "MEDCV_ALDO1",
  "MEDCV_DIUR1",
  "MEDCV_LDIUR1"
)

# Some relevant biomarkers often assessed in HF
LAB_VALUES <- c(
  "HB1",
  "NTPBNP1",
  "NPX_195_NT_PRO_BNP",
  "BNP1",
  "SCREA1",
  "UREA1",
  "POT1",
  "SOD1",
  "TBILI1",
  "BIOADM",
  "EGFRBL"
)

# Putting all these column names in one vector --> Use these columns for generating the baseline characteristics
CHARS.4.BASELINE <- c(
  DEMOGRAPHICS,
  HF_CHARS,
  HISTORY,
  "SMOKING3"
)

#generating the baseline characteristics in index cohort
BIOSTAT <- BIOSTAT %>% 
  mutate(M.factor = as.factor(M)) %>% 
  mutate(H.factor = as.factor(H)) %>% 
  mutate(MH.factor = as.factor(MH))
  


baseline_table <- compareGroups(
  formula = SMOKING3 ~ .,  # SMOKING as grouping variable, all other variables as descriptives
  data = BIOSTAT[, CHARS.4.BASELINE]  # subset to only your characteristics
)

skim_without_charts(BIOSTAT)

# Create table with better formatting
baseline_table_output <- createTable(
  baseline_table,
  hide.no = FALSE,           # Show missing values
  digits = 1,                # 1 decimal places for continuous vars
  show.p.overall = TRUE,     # Show p-values
  show.p.trend = FALSE       # Don't show trend p-values
)

# Display with nicer formatting
baseline_table_output

#log rank test for time-to-event outcome (variables: MH, M,H, MHdays, Mdays, Hdays)
event_level <- "1" 
Validation_$event <- as.integer(as.character(Validation_$MH.factor) == event_level)#BIOSTAT or Validation_

# Overall log-rank comparison across SMOKING3 groups
lr <- survdiff(
  Surv(MHdays, event) ~ SMOKING3,
  data = Validation_, #BIOSTAT or Validation_
  rho = 0
)

df <- length(lr$n) - 1
p_value <- pchisq(lr$chisq, df = df, lower.tail = FALSE)

lr
p_value


#######generating the baseline characteristics in validation cohort#########


Validation <- Validation_ %>% 
  mutate(SMOKING = as.factor(case_when(CVRF_SMOKE == 1 ~ "1",#1 never/past; 2 current smoker
                                       CVRF_SMOKE == 2 ~ "1",
                                       CVRF_SMOKE == 3 ~ "2")))

Validation_ <- Validation_ %>% 
  mutate(M.factor = as.factor(M)) %>% 
  mutate(H.factor = as.factor(H)) %>% 
  mutate(MH.factor = as.factor(MH))

skim_without_charts(Validation)

DEMOGRAPHICS <- c(
  "AGE",
  "DEM_SEX",
  "RACE2",
  "PE_BMI",
  "M.factor",
  "H.factor",
  "MH.factor"
)

# Heart failure characteristics
HF_CHARS <- c(
  "MHHF_NYHA", #NYHA class (Visit 1 (1))
  "ECRV_LVEF" #
)

# History of diseases
HISTORY <- c(
  "CVRF_NY8",
  "CVRF_NY10",
  "CVRF_NY1",
  "CVRF_NY6",
  "CVRF_NY14"
)

CHARS.4.BASELINE <- c(
  DEMOGRAPHICS,
  HF_CHARS,
  HISTORY,
  "SMOKING3"
)

baseline_table <- compareGroups(
  formula = SMOKING3 ~ .,  # SMOKING as grouping variable, all other variables as descriptives
  data = Validation_[, CHARS.4.BASELINE]  # subset to only your characteristics
)

# Create table with better formatting
baseline_table_output <- createTable(
  baseline_table,
  hide.no = FALSE,           # Show missing values
  digits = 1,                # 1 decimal places for continuous vars
  show.p.overall = TRUE,     # Show p-values
  show.p.trend = FALSE       # Don't show trend p-values
)

# Display with nicer formatting
baseline_table_output

export2word(baseline_table_output, "baseline table validation.docx")

#the subset of variables I need (use this code)
my_vars <- c("DEM_PT",
             DEMOGRAPHICS,
             HF_CHARS,
             HISTORY,
             LIFESTYLE,
             SYMPTOMS_SIGNS_HF,
             PHYSICAL_EXAM,
             LAB_VALUES,
             MEDICATIONS,
             "SMOKING3",
             "CVRF_SMOKE",
             "SMOKING",
             "M","H","Mdays","Hdays","MH","MHdays")

#create a dataset with only these variables I need
BIOSTAT_ <- BIOSTAT %>% select(all_of(my_vars)) %>% 
  rename(NPX_NT_PRO_BNP = NPX_195_NT_PRO_BNP) #renamed the NTproBNP variable because there is the same protein 
#from CVD3 panel...if the name is the same there will be problems in later analysis
  

BIOSTAT_$DEM_SEX <- factor(
  BIOSTAT_$DEM_SEX,
  levels = c(1, 2),
  labels = c("Male", "Female")
)

BIOSTAT_ <- BIOSTAT_ %>% 
  mutate(SMOKING = as.factor(case_when(CVRF_SMOKE == 1 ~ "1",#1 non/past; 2 current smoker
                             CVRF_SMOKE == 2 ~ "1",
                             CVRF_SMOKE == 3 ~ "2"))) %>% 
  mutate(
    RACE2 = as.factor(case_when(
      RACE == "Caucasian" ~ 1,
      RACE %in% c("Black", "Asian", "Oriental", "Hispanic", "Other") ~ 0,
      TRUE ~ NA_real_)))

skim_without_charts(BIOSTAT)


#############################################Smoking status and HF clinical outcome#############################################
##################in index cohort####################
BIOSTAT_$SMOKING3 <- relevel(factor(BIOSTAT_$SMOKING3), ref = "Never")

fit_km <- survfit(Surv(Mdays, M) ~ SMOKING, data = BIOSTAT_)

#the final model
fit_cox <- coxph(Surv(MHdays, MH) ~ SMOKING3+AGE+DEM_SEX
                 , data = BIOSTAT_)

vif(fit_cox)
ph_test <- cox.zph(fit_cox)
print(ph_test)
fit_cox


##################in validation cohort####################
skim_without_charts(Validation_)

Validation_ <- Validation_ %>% 
  mutate(SMOKING3 = case_when(SMOKING3 == 1 ~ "Never",
                              SMOKING3 == 2 ~ "Past",
                              SMOKING3 == 3 ~ "Current")) 

Validation_$SMOKING3 <- as.factor(Validation_$SMOKING3)

Validation_$SMOKING3 <- relevel(factor(Validation_$SMOKING3), ref = "Never")

Validation_$DEM_SEX <- factor(
  Validation_$DEM_SEX,
  levels = c(1, 2),
  labels = c("Male", "Female")
)

#the final model
fit_cox <- coxph(Surv(MHdays, MH) ~ SMOKING3+AGE+DEM_SEX 
                 , data = Validation_)

ph_test <- cox.zph(fit_cox)
print(ph_test)
fit_cox


################################Volcano plots in index cohort#############################
##########################################################################################
##########################################################################################
##########################################################################################

#load olink omics data

CVD2 <- read_sav("Biostat_IndexCohort_Olink_CVD2_panel.sav")
skim_without_charts(CVD2)
CVD3 <- read_sav("Biostat_IndexCohort_Olink_CVD3_panel.sav")
IMMUNE <- read_sav("Biostat_IndexCohort_Olink_Immune_panel.sav")
ONCO <- read_sav("Biostat_IndexCohort_Olink_Oncology_panel.sav")

##############################################smoking data processing###########################
#Never: 940, Past: 1220, Current: 353
#never/past smoker VS current smoker
BIOSTAT_1 <- BIOSTAT_ %>% 
  mutate(SMOKING = as.factor(case_when(CVRF_SMOKE == 1 ~ "1",#1 never/past; 2 current smoker
                                       CVRF_SMOKE == 2 ~ "1",
                                       CVRF_SMOKE == 3 ~ "2")))

#never smoker VS current smoker
BIOSTAT_1 <- BIOSTAT_ %>% 
  filter(CVRF_SMOKE != 2) %>% 
  mutate(SMOKING = as.factor(case_when(CVRF_SMOKE == 1 ~ "1",#1 never; 2 current smoker
                                       CVRF_SMOKE == 3 ~ "2")))

#past smoker VS current smoker
BIOSTAT_1 <- BIOSTAT_ %>% 
  filter(CVRF_SMOKE != 1) %>% 
  mutate(SMOKING = as.factor(case_when(CVRF_SMOKE == 2 ~ "1",#1 past; 2 current smoker
                                       CVRF_SMOKE == 3 ~ "2")))

#never smoker VS past smoker
BIOSTAT_1 <- BIOSTAT_ %>% 
  filter(CVRF_SMOKE != 3) %>% 
  mutate(SMOKING = as.factor(case_when(CVRF_SMOKE == 1 ~ "1",#1 never; 2 past smoker
                                       CVRF_SMOKE == 2 ~ "2")))

skim_without_charts(BIOSTAT_)
##########################################volcano plot CVD2 panel##############################
#combine CVD2 and BIOSTAT_1 dataset
BIOSTAT.CVD2 <- full_join(BIOSTAT_1, CVD2, by = "DEM_PT") %>% 
  mutate(BMI = as.numeric(PE_WGT1/((PE_HGT1*0.01)^2)))

#protein list for CVD2 panel
protein_list <- grep("^OID", names(BIOSTAT.CVD2), value = TRUE)

#extract protein names and make an annotation file
protein_cols <- grep("^OID", names(BIOSTAT.CVD2), value = TRUE)#Get protein columns (OID only)
labels <- sapply(BIOSTAT.CVD2[protein_cols], function(x) attr(x, "label"))#extract labels
protein_names <- sub(".*:", "", labels)#extract protein names
annotation_CVD2 <- data.frame(
  protein = protein_cols,
  protein_name = protein_names,
  stringsAsFactors = FALSE
)#annotation file


#fit a model for each protein, and extract beta, standard error and p value for smoking, and FDR correction
results <- lapply(protein_list, function(prot) {
  model <- lm(as.formula(paste(prot, "~ SMOKING + AGE + DEM_SEX + BMI + SYMP_NYHA1 +NPX_NT_PRO_BNP")), data = BIOSTAT.CVD2)
  coef <- summary(model)$coefficients
  
  row <- coef[grep("SMOKING2", rownames(coef))[1], ]
  
  data.frame(
    protein = prot,
    beta = row["Estimate"],
    SE = row["Std. Error"],
    pvalue = row["Pr(>|t|)"]
  )
})

results_df <- do.call(rbind, results)
results_df$FDR <- p.adjust(results_df$pvalue, method = "fdr")#FDR correction

#log transform p value
results_df$neglog10p <- -log10(results_df$pvalue)

#add real protein names to results_df
results_df <- results_df %>%
  left_join(annotation_CVD2, by = "protein")

#save a copy of the results_df for making a combined volcano plot later
results_CVD2 <- results_df



##########################################volcano plot CVD3 panel##############################

#combine CVD3 and BIOSTAT_1 dataset
BIOSTAT.CVD3 <- full_join(BIOSTAT_1, CVD3, by = "DEM_PT") %>% 
  mutate(BMI = as.numeric(PE_WGT1/((PE_HGT1*0.01)^2)))


#protein list for CVD3 panel
protein_list <- grep("^NPX", names(BIOSTAT.CVD3), value = TRUE)

#fit a model for each protein, and extract beta and p value for smoking, and FDR correction
results <- lapply(protein_list, function(prot) {
  model <- lm(as.formula(paste(prot, "~ SMOKING + AGE + DEM_SEX + BMI + SYMP_NYHA1+NPX_NT_PRO_BNP")), data = BIOSTAT.CVD3)
  coef <- summary(model)$coefficients
  
  row <- coef[grep("SMOKING2", rownames(coef))[1], ]
  
  data.frame(
    protein = prot,
    beta = row["Estimate"],
    SE = row["Std. Error"],
    pvalue = row["Pr(>|t|)"]
  )
})

results_df <- do.call(rbind, results)
results_df$FDR <- p.adjust(results_df$pvalue, method = "fdr")#FDR correction

#rename protein names
results_df$protein_name <- sub("^NPX_[0-9]+_", "", results_df$protein)

#log transform p value
results_df$neglog10p <- -log10(results_df$pvalue)

#save a copy of the results_df for making a combined volcano plot later
results_CVD3 <- results_df



##########################################volcano plot IMMUNE panel##############################
#combine IMMUNE and BIOSTAT_1 dataset
BIOSTAT.IMM <- full_join(BIOSTAT_1, IMMUNE, by = "DEM_PT") %>% 
  mutate(BMI = as.numeric(PE_WGT1/((PE_HGT1*0.01)^2)))

#protein list for CVD2 panel
protein_list <- grep("^OID", names(BIOSTAT.IMM), value = TRUE)

#extract protein names and make an annotation file
protein_cols <- grep("^OID", names(BIOSTAT.IMM), value = TRUE)#Get protein columns (OID only)
labels <- sapply(BIOSTAT.IMM[protein_cols], function(x) attr(x, "label"))#extract labels
protein_names <- sub(".*:", "", labels)#extract protein names
annotation_IMM <- data.frame(
  protein = protein_cols,
  protein_name = protein_names,
  stringsAsFactors = FALSE
)#annotation file


#fit a model for each protein, and extract beta and p value for smoking, and FDR correction
results <- lapply(protein_list, function(prot) {
  model <- lm(as.formula(paste(prot, "~ SMOKING + AGE + DEM_SEX + BMI + SYMP_NYHA1+NPX_NT_PRO_BNP")), data = BIOSTAT.IMM)
  coef <- summary(model)$coefficients
  
  row <- coef[grep("SMOKING2", rownames(coef))[1], ]
  
  data.frame(
    protein = prot,
    beta = row["Estimate"],
    SE = row["Std. Error"],
    pvalue = row["Pr(>|t|)"]
  )
})

results_df <- do.call(rbind, results)
results_df$FDR <- p.adjust(results_df$pvalue, method = "fdr")#FDR correction

#log transform p value
results_df$neglog10p <- -log10(results_df$pvalue)

#add real protein names to results_df
results_df <- results_df %>%
  left_join(annotation_IMM, by = "protein")

#save a copy of the results_df for making a combined volcano plot later
results_IMM <- results_df


##########################################volcano plot oncology panel##############################
#combine BIOSTAT and oncology dataset
BIOSTAT.ONCO <- full_join(BIOSTAT_1, ONCO, by = "DEM_PT") %>% 
  mutate(BMI = as.numeric(PE_WGT1/((PE_HGT1*0.01)^2)))

#protein list for CVD2 panel
protein_list <- grep("^OID", names(BIOSTAT.ONCO), value = TRUE)

#extract protein names and make an annotation file
protein_cols <- grep("^OID", names(BIOSTAT.ONCO), value = TRUE)#Get protein columns (OID only)
labels <- sapply(BIOSTAT.ONCO[protein_cols], function(x) attr(x, "label"))#extract labels
protein_names <- sub(".*:", "", labels)#extract protein names
annotation_ONCO <- data.frame(
  protein = protein_cols,
  protein_name = protein_names,
  stringsAsFactors = FALSE
)#annotation file


#fit a model for each protein, and extract beta and p value for smoking, and FDR correction
results <- lapply(protein_list, function(prot) {
  model <- lm(as.formula(paste(prot, "~ SMOKING + AGE + DEM_SEX + BMI + SYMP_NYHA1+NPX_NT_PRO_BNP")), data = BIOSTAT.ONCO)
  coef <- summary(model)$coefficients
  
  row <- coef[grep("SMOKING2", rownames(coef))[1], ]
  
  data.frame(
    protein = prot,
    beta = row["Estimate"],
    SE = row["Std. Error"],
    pvalue = row["Pr(>|t|)"]
  )
})

results_df <- do.call(rbind, results)
results_df$FDR <- p.adjust(results_df$pvalue, method = "fdr")#FDR correction

#log transform p value
results_df$neglog10p <- -log10(results_df$pvalue)

#add real protein names to results_df
results_df <- results_df %>%
  left_join(annotation_ONCO, by = "protein")

#save a copy of the results_df for making a combined volcano plot later
results_ONCO <- results_df



################################combine the volcano plot from the index cohort##########################################
#combine the results from each panel
volcano_index <- rbind(results_CVD2, results_CVD3, results_IMM, results_ONCO)

#Average values for other variables,remove "protein" and remove duplicate protein_name 
volcano_index <-volcano_index %>% 
  select(-(protein)) %>% 
  group_by(protein_name) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))

#Perform new FDR correction
volcano_index$FDR <- p.adjust(volcano_index$pvalue, method = "fdr")#FDR correction

#protein_name "CXL7" should be "CXCL17"
volcano_index$protein_name[volcano_index$protein_name == "CXL17"] <- "CXCL17"

#add statistical significance highlight
volcano_index$regulation <- "Not significant"

volcano_index$regulation[volcano_index$FDR < 0.05 & volcano_index$beta > 0] <- "Upregulated"
volcano_index$regulation[volcano_index$FDR < 0.05 & volcano_index$beta < 0] <- "Downregulated"

volcano_index$regulation <- factor(
  volcano_index$regulation,
  levels = c("Downregulated", "Not significant", "Upregulated")
)


#remove "-" from protein names

# This removes underscores completely: A_B becomes AB
volcano_index$protein_name <- gsub("_", "", volcano_index$protein_name)

# volcano plot with annotation
ggplot(volcano_index, aes(x = beta, y = neglog10p, color = regulation)) +
  geom_point(alpha = 0.8, size = 2) +
  # label top proteins 
  geom_text_repel(
    data = subset(volcano_index, FDR < 0.05),
    aes(label = protein_name),
    size = 3,
    max.overlaps = 20,
    show.legend = FALSE
  ) +
  
  scale_color_manual(
    values = c(
      "Downregulated" = "#0072B2",   # blue
      "Not significant" = "grey70",
      "Upregulated" = "#D55E00"      # red/orange
    )
  ) +
  
  # Set y-axis scale (0 to 50)
  scale_y_continuous(limits = c(0, 50), expand = c(0.02, 0)) +
  
  # Set x-axis scale (-0.50 to 0.75)
  scale_x_continuous(limits = c(-0.55, 0.80), expand = c(0.02, 0)) +
  
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  # Add annotation for P = 0.05 line
  annotate("text", 
           x = Inf,  # right edge
           y = -log10(0.05) - 0.3,  # slightly below the line
           label = "p = 0.05",
           hjust = 1.1,  # adjust to position text left of the edge
           size = 3.5,
           color = "black") +
  
  labs(
    x = expression(Log[2]~(fold~change)~"in"~index~cohort),
    y = expression(-Log[10]~(p~value)),
    title = "",
    color = "Regulation"
  ) +
  
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),
    panel.grid.minor = element_blank()
  )

#current vs never
volcano_index.currentvsnever.sens <- volcano_index
write_xlsx(volcano_index, "current smoker vs never smoker sens.xlsx")

#current vs past
volcano_index.currentvspast.sens <- volcano_index
write_xlsx(volcano_index, "current smoker vs past smoker sens.xlsx")

#past vs never
volcano_index.pastvsnever.sens <- volcano_index
write_xlsx(volcano_index, "past smoker vs never smoker sens.xlsx")

################################validation cohort 20260429################################
##########################################################################################
##########################################################################################
##########################################################################################

Validation.outcome <- read_sav("BIOSTAT_validation_outcome_including_safe_haven.sav")
Validation.pre <- read_sav("Biostat_ValidationCohort_Demographics_Clinical_medical_and_Family_History.sav")
Validation.add <- read_sav("Biostat_ValidationCohort_PhysicalExam_SymptomsHF_and_Echo.sav")
Validation.lab <- read_sav("Biostat_ValidationCohort_Laboratory_simple.sav")

Validation.all <- Validation.pre %>%
  full_join(Validation.add, by = "DEM_PT") %>%
  full_join(Validation.outcome, by = "DEM_PT") %>% 
  full_join(Validation.lab, by = "DEM_PT")

#add NPX_195_NT_PRO_BNP froom CVD3.V panel
Validation.all <- CVD3.V %>% 
  select(DEM_PT, NPX_195_NT_PRO_BNP) %>% 
  full_join(Validation.all, by = "DEM_PT") %>% 
  rename(NPX_NT_PRO_BNP = NPX_195_NT_PRO_BNP)

skim_without_charts(Validation.all)

#Creates a key and annotation file checking what is exactly in the columns
column_keys <- vapply(
  Validation.all,
  function(x) {
    lab <- attr(x, "label")
    if (is.null(lab)) NA_character_ else as.character(lab)
  },
  character(1)
)

ANNOTATION.Validation <- data.frame(
  COLNAME = names(Validation.all),
  ANNOTATION = column_keys,
  stringsAsFactors = FALSE
)

# Calculating age based on visit and date of birth
# We use 365.25 because of leap year;)
Validation_ <- Validation.all %>%
  mutate(AGE = round(as.numeric(difftime(as.Date(VIS_DTE), as.Date(DEM_DBIR), units = "days"))/365.25,0)) %>%
  mutate(RACE2 = ifelse(DEM_RACE == 1, 1, 0)) %>%
  mutate(RACE2 = as.factor(RACE2)) %>% 
  mutate(SMOKING3 = as.factor(case_when(CVRF_SMOKE == 1 ~ "1",
                                       CVRF_SMOKE == 2 ~ "2",
                                       CVRF_SMOKE == 3 ~ "3")))

skim_without_charts(Validation_)

skim_without_charts(CVD2.V)

##############################################smoking data processing###########################
#Validation: Never: 884, Past: 605, Current: 237
#Index: Never: 940, Past: 1220, Current: 353
#never/past smoker VS current smoker
Validation <- Validation_ %>% 
  mutate(SMOKING = as.factor(case_when(CVRF_SMOKE == 1 ~ "1",#1 never/past; 2 current smoker
                                       CVRF_SMOKE == 2 ~ "1",
                                       CVRF_SMOKE == 3 ~ "2")))
skim_without_charts(Validation.all)

#never smoker VS current smoker 
Validation <- Validation_ %>% 
  filter(CVRF_SMOKE != 2) %>% 
  mutate(SMOKING = as.factor(case_when(CVRF_SMOKE == 1 ~ "1",#1 never; 2 current smoker
                                       CVRF_SMOKE == 3 ~ "2")))

#past smoker VS current smoker
Validation <- Validation_ %>% 
  filter(CVRF_SMOKE != 1) %>% 
  mutate(SMOKING = as.factor(case_when(CVRF_SMOKE == 2 ~ "1",#1 past; 2 current smoker
                                       CVRF_SMOKE == 3 ~ "2")))

#never smoker VS past smoker
Validation <- Validation_ %>% 
  filter(CVRF_SMOKE != 3) %>% 
  mutate(SMOKING = as.factor(case_when(CVRF_SMOKE == 1 ~ "1",#1 never; 2 past smoker
                                       CVRF_SMOKE == 2 ~ "2")))

#load olink data
CVD2.V <- read_sav("Biostat_ValidationCohort_Olink_CVD2_panel.sav")
CVD3.V <- read_sav("Biostat_ValidationCohort_Olink_CVD3_panel.sav")
IMMUNE.V <- read_sav("Biostat_ValidationCohort_Olink_Immune_panel.sav")
ONCO.V <- read_sav("Biostat_ValidationCohort_Olink_Oncology_panel.sav")


##########################################volcano plot CVD2 panel##############################
#combine CVD2 and BIOSTAT_1 dataset
BIOSTAT.CVD2 <- full_join(Validation, CVD2.V, by = "DEM_PT") 

#protein list for CVD2 panel
protein_list <- grep("^OID", names(BIOSTAT.CVD2), value = TRUE)

#extract protein names and make an annotation file
protein_cols <- grep("^OID", names(BIOSTAT.CVD2), value = TRUE)#Get protein columns (OID only)
labels <- sapply(BIOSTAT.CVD2[protein_cols], function(x) attr(x, "label"))#extract labels
protein_names <- sub(".*:", "", labels)#extract protein names
annotation_CVD2 <- data.frame(
  protein = protein_cols,
  protein_name = protein_names,
  stringsAsFactors = FALSE
)#annotation file


#fit a model for each protein, and extract beta and p value for smoking, and FDR correction
results <- lapply(protein_list, function(prot) {
  model <- lm(as.formula(paste(prot, "~ SMOKING + AGE + DEM_SEX + PE_BMI+ 	
MHHF_NYHA +NPX_NT_PRO_BNP")), data = BIOSTAT.CVD2)
  coef <- summary(model)$coefficients
  
  row <- coef[grep("SMOKING2", rownames(coef))[1], ]
  
  data.frame(
    protein = prot,
    beta = row["Estimate"],
    SE = row["Std. Error"],
    pvalue = row["Pr(>|t|)"]
  )
})

results_df <- do.call(rbind, results)
results_df$FDR <- p.adjust(results_df$pvalue, method = "fdr")#FDR correction

#log transform p value
results_df$neglog10p <- -log10(results_df$pvalue)

#add real protein names to results_df
results_df <- results_df %>%
  left_join(annotation_CVD2, by = "protein")

#save a copy of the results_df for making a combined volcano plot later
results_CVD2 <- results_df


##########################################volcano plot CVD3 panel##############################

#combine CVD3 and Validation dataset
BIOSTAT.CVD3 <- full_join(Validation, CVD3.V, by = "DEM_PT") 


#protein list for CVD3 panel
protein_list <- grep("^NPX", names(BIOSTAT.CVD3), value = TRUE)

#fit a model for each protein, and extract beta and p value for smoking, and FDR correction
results <- lapply(protein_list, function(prot) {
  model <- lm(as.formula(paste(prot, "~ SMOKING + AGE + DEM_SEX + PE_BMI+ 	
MHHF_NYHA+NPX_NT_PRO_BNP")), data = BIOSTAT.CVD3)
  coef <- summary(model)$coefficients
  
  row <- coef[grep("SMOKING2", rownames(coef))[1], ]
  
  data.frame(
    protein = prot,
    beta = row["Estimate"],
    SE = row["Std. Error"],
    pvalue = row["Pr(>|t|)"]
  )
})

results_df <- do.call(rbind, results)
results_df$FDR <- p.adjust(results_df$pvalue, method = "fdr")#FDR correction

#rename protein names
results_df$protein_name <- sub("^NPX_[0-9]+_", "", results_df$protein)

#log transform p value
results_df$neglog10p <- -log10(results_df$pvalue)

#save a copy of the results_df for making a combined volcano plot later
results_CVD3 <- results_df



##########################################volcano plot IMMUNE panel##############################
#combine IMMUNE and Validation dataset
BIOSTAT.IMM <- full_join(Validation, IMMUNE.V, by = "DEM_PT") 

#protein list for CVD2 panel
protein_list <- grep("^OID", names(BIOSTAT.IMM), value = TRUE)

#extract protein names and make an annotation file
protein_cols <- grep("^OID", names(BIOSTAT.IMM), value = TRUE)#Get protein columns (OID only)
labels <- sapply(BIOSTAT.IMM[protein_cols], function(x) attr(x, "label"))#extract labels
protein_names <- sub(".*:", "", labels)#extract protein names
annotation_IMM <- data.frame(
  protein = protein_cols,
  protein_name = protein_names,
  stringsAsFactors = FALSE
)#annotation file


#fit a model for each protein, and extract beta and p value for smoking, and FDR correction
results <- lapply(protein_list, function(prot) {
  model <- lm(as.formula(paste(prot, "~ SMOKING + AGE + DEM_SEX + PE_BMI+ 	
MHHF_NYHA+NPX_NT_PRO_BNP")), data = BIOSTAT.IMM)
  coef <- summary(model)$coefficients
  
  row <- coef[grep("SMOKING2", rownames(coef))[1], ]
  
  data.frame(
    protein = prot,
    beta = row["Estimate"],
    SE = row["Std. Error"],
    pvalue = row["Pr(>|t|)"]
  )
})

results_df <- do.call(rbind, results)
results_df$FDR <- p.adjust(results_df$pvalue, method = "fdr")#FDR correction

#log transform p value
results_df$neglog10p <- -log10(results_df$pvalue)

#add real protein names to results_df
results_df <- results_df %>%
  left_join(annotation_IMM, by = "protein")

#save a copy of the results_df for making a combined volcano plot later
results_IMM <- results_df



##########################################volcano plot oncology panel##############################
#combine BIOSTAT and oncology dataset
BIOSTAT.ONCO <- full_join(Validation, ONCO.V, by = "DEM_PT") 

#protein list for CVD2 panel
protein_list <- grep("^OID", names(BIOSTAT.ONCO), value = TRUE)

#extract protein names and make an annotation file
protein_cols <- grep("^OID", names(BIOSTAT.ONCO), value = TRUE)#Get protein columns (OID only)
labels <- sapply(BIOSTAT.ONCO[protein_cols], function(x) attr(x, "label"))#extract labels
protein_names <- sub(".*:", "", labels)#extract protein names
annotation_ONCO <- data.frame(
  protein = protein_cols,
  protein_name = protein_names,
  stringsAsFactors = FALSE
)#annotation file


#fit a model for each protein, and extract beta and p value for smoking, and FDR correction
results <- lapply(protein_list, function(prot) {
  model <- lm(as.formula(paste(prot, "~ SMOKING + AGE + DEM_SEX + PE_BMI+ 	
MHHF_NYHA+NPX_NT_PRO_BNP")), data = BIOSTAT.ONCO)
  coef <- summary(model)$coefficients
  
  row <- coef[grep("SMOKING2", rownames(coef))[1], ]
  
  data.frame(
    protein = prot,
    beta = row["Estimate"],
    SE = row["Std. Error"],
    pvalue = row["Pr(>|t|)"]
  )
})

results_df <- do.call(rbind, results)
results_df$FDR <- p.adjust(results_df$pvalue, method = "fdr")#FDR correction

#log transform p value
results_df$neglog10p <- -log10(results_df$pvalue)

#add real protein names to results_df
results_df <- results_df %>%
  left_join(annotation_ONCO, by = "protein")

#save a copy of the results_df for making a combined volcano plot later
results_ONCO <- results_df



################################combine the volcano plot from the validation cohort##########################################
#combine the results from each panel
volcano_val <- rbind(results_CVD2, results_CVD3, results_IMM, results_ONCO)

#Average values for other variables,remove "protein" and remove duplicate protein_name 
volcano_val <-volcano_val %>% 
  select(-(protein)) %>% 
  group_by(protein_name) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))

#Perform new FDR correction
volcano_val$FDR <- p.adjust(volcano_val$pvalue, method = "fdr")#FDR correction

#protein_name "CXL7" should be "CXCL17"
volcano_val$protein_name[volcano_val$protein_name == "CXL17"] <- "CXCL17"

#add statistical significance highlight
volcano_val$regulation <- "Not significant"

volcano_val$regulation[volcano_val$FDR < 0.05 & volcano_val$beta > 0] <- "Upregulated"
volcano_val$regulation[volcano_val$FDR < 0.05 & volcano_val$beta < 0] <- "Downregulated"

volcano_val$regulation <- factor(
  volcano_val$regulation,
  levels = c("Downregulated", "Not significant", "Upregulated")
)

volcano_val$protein_name <- gsub("_", "", volcano_val$protein_name)

# volcano plot with annotation
ggplot(volcano_val, aes(x = beta, y = neglog10p, color = regulation)) +
  geom_point(alpha = 0.8, size = 2) +
  # label top proteins 
  geom_text_repel(
    data = subset(volcano_val, FDR < 0.05),
    aes(label = protein_name),
    size = 3,
    max.overlaps = 20,
    show.legend = FALSE
  ) +
  
  scale_color_manual(
    values = c(
      "Downregulated" = "#0072B2",   # blue
      "Not significant" = "grey70",
      "Upregulated" = "#D55E00"      # red/orange
    )
  ) +
  
  # Set y-axis scale (0 to 50)
  scale_y_continuous(limits = c(0, 50), expand = c(0.02, 0)) +
  
  # Set x-axis scale (-0.50 to 0.75)
  scale_x_continuous(limits = c(-0.55, 0.80), expand = c(0.02, 0)) +
  
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  # Add annotation for P = 0.05 line
  annotate("text", 
           x = Inf,  # right edge
           y = -log10(0.05) - 0.5,  # slightly below the line
           label = "p = 0.05",
           hjust = 1.1,  # adjust to position text left of the edge
           size = 3.5,
           color = "black") +
  
  labs(
    x = expression(Log[2]~(fold~change)~"in"~validation~cohort),
    y = expression(-Log[10]~(p~value)),
    title = "",
    color = "Regulation"
  ) +
  
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),
    panel.grid.minor = element_blank()
  )

#current vs never
volcano_val.currentvsnever.sens <- volcano_val
write_xlsx(volcano_val, "current smoker vs never smoker val sens.xlsx")

#current vs past
volcano_val.currentvspast.sens <- volcano_val
write_xlsx(volcano_val, "current smoker vs past smoker val sens.xlsx")

#past vs never
volcano_val.pastvsnever.sens <- volcano_val
write_xlsx(volcano_val, "past smoker vs never smoker val sens.xlsx")

################################Replication results################################
##########################################################################################
##########################################################################################
##########################################################################################
###############use Venn diagram to show the replication 20260918####################
#For the variable smoking, we made three comparisons: current vs never; current vs past,
# and past vs never. For each comparison: we have 4 analyses: index and index sensitivity, validation
# and validation sensitivity. So we will make a Venn diagram for each smoking comparison group

make_count_venn <- function(index, index_sens,
                            validation, validation_sens,
                            title) {
  
  protein_sets <- list(
    "Index" = index %>%
      filter(FDR < 0.05) %>%
      pull(protein_name) %>%
      unique(),
    
    "Index sensitivity" = index_sens %>%
      filter(FDR < 0.05) %>%
      pull(protein_name) %>%
      unique(),
    
    "Validation" = validation %>%
      filter(FDR < 0.05) %>%
      pull(protein_name) %>%
      unique(),
    
    "Validation sensitivity" = validation_sens %>%
      filter(FDR < 0.05) %>%
      pull(protein_name) %>%
      unique()
  )
  
  # Extract set-label positions
  venn_object <- Venn(protein_sets)
  venn_data <- process_data(venn_object)
  
  set_labels <- venn_setlabel(venn_data) %>%
    mutate(
      X = if_else(
        name == "Validation sensitivity",
        X - 0.15,
        X
      ),
      hjust_value = if_else(
        name == "Validation sensitivity",
        0,
        0.5
      )
    )
  
  ggVennDiagram(
    protein_sets,
    label = "count",
    label_size = 5.5,
    label_alpha = 0,
    edge_size = 0.7,
    set_size = 0
  ) +
    geom_label(
      data = set_labels,
      aes(
        x = X,
        y = Y,
        label = name,
        hjust = hjust_value
      ),
      inherit.aes = FALSE,
      size = 3.8,
      fill = "white",
      label.size = NA,
      label.padding = grid::unit(0.08, "lines")
    ) +
    scale_fill_gradient(
      low = "white",
      high = "#4C78A8"
    ) +
    scale_x_continuous(
      expand = expansion(
        mult = c(0.15, 0.30)
      )
    ) +
    scale_y_continuous(
      expand = expansion(mult = 0.12)
    ) +
    labs(title = title) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 13,
        hjust = 0.5
      ),
      plot.margin = margin(
        t = 10,
        r = 20,
        b = 10,
        l = 20
      ),
      legend.position = "none"
    )
}

# Create the three count-based diagrams:
venn_current_never <- make_count_venn(
  volcano_index.currentvsnever,
  volcano_index.currentvsnever.sens,
  volcano_val.currentvsnever,
  volcano_val.currentvsnever.sens,
  "Current versus never smokers"
)

venn_current_past <- make_count_venn(
  volcano_index.currentvspast,
  volcano_index.currentvspast.sens,
  volcano_val.currentvspast,
  volcano_val.currentvspast.sens,
  "Current versus past smokers"
)

venn_past_never <- make_count_venn(
  volcano_index.pastvsnever,
  volcano_index.pastvsnever.sens,
  volcano_val.pastvsnever,
  volcano_val.pastvsnever.sens,
  "Past versus never smokers"
)

venn_current_never
venn_current_past
venn_past_never


#use never/past smokers VS current smokers, and generate the datasets: volcano_index & volcano_val
#make replication table

volcano_val_copy <- volcano_val.currentvspast %>% 
  rename_with(~ paste0(.x, "_val"), .cols = -protein_name)

volcano_rep <- volcano_index.currentvspast %>% 
  filter(regulation != "Not significant") %>% 
  left_join(volcano_val_copy, by = "protein_name")

#make beta-beta effect size plot

beta_beta_plot <- ggplot(volcano_rep, aes(x = beta, y = beta_val, label = protein_name)) +
  geom_point(alpha = 0.5, size = 2, color = "black", fontface = "bold") +  # explicit point color
  geom_text_repel(
    size = 4,
    color = "black",
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
    force = 2,
    min.segment.length = 0
  ) +  # explicit text color + force parameter
  geom_smooth(method = "lm", se = TRUE, color = "orange", fill = "gray80") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue", linewidth = 1) +
  theme_minimal(base_size = 12) +
  labs(x = "Index coefficients", y = "Validation coefficients") +
  theme(
    text = element_text(color = "black", face = "bold"),  # Global text color
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 12)
  )
# Add correlation annotation
cor_value <- cor(volcano_rep$beta, volcano_rep$beta_val, 
                 use = "complete.obs")
r2_value <- cor_value^2

# Get regression line coefficients for annotation
lm_model <- lm(beta_val ~ beta, data = volcano_rep)
slope <- round(coef(lm_model)[2], 3)
intercept <- round(coef(lm_model)[1], 3)

beta_beta_plot <- beta_beta_plot +
  annotate("text", 
           x = min(volcano_rep$beta, na.rm = TRUE), 
           y = max(volcano_rep$beta_val, na.rm = TRUE),
           label = paste("r =", round(cor_value, 3), 
                         "\nR² =", round(r2_value, 3)),
           hjust = 0, vjust = 1, size = 4) +
  
  # Add annotation for the identity line
  annotate("text",
           x = 0.25,
           y = -0.30,
           label = "--- Identity Line",
           hjust = 0, vjust = 0, size = 3.5,
           color = "blue") +
  
  # Add annotation for the regression line
  annotate("text",
           x = 0.25,
           y = -0.25,
           label = "— Regression Line",
           hjust = 0, vjust = 0, size = 3.5,
           color = "orange")

print(beta_beta_plot)



################################proteins and HR outcomes - Index cohort################################
##########################################################################################
##########################################################################################
##########################################################################################

#get the current vs never proteins that are replicated in vaildation cohort

volcano_val_copy <- volcano_val.currentvsnever %>% 
  rename_with(~ paste0(.x, "_val"), .cols = -protein_name)

volcano_rep1 <- volcano_index.currentvsnever %>% 
  filter(regulation != "Not significant") %>% 
  left_join(volcano_val_copy, by = "protein_name") %>% 
  filter(regulation_val != "Not significant")


#match the protein_name in "volcano_rep1" with the protein variable names in the full dataset

temp <- rbind(results_CVD2, results_CVD3, results_IMM, results_ONCO)

protein_rep <- volcano_rep1 %>% 
  select(protein_name) %>% 
  mutate(protein_name = ifelse(protein_name == "CXCL17", "CXL17", protein_name)) %>% 
  left_join(temp, by = "protein_name") %>% 
  mutate(
    protein = case_when(
      protein_name == "TRAP"  ~ "NPX_148_TR_AP", #otherwise cannot find the protein values for these two proteins
      protein_name == "IL2RA" ~ "NPX_109_IL2_RA",
      TRUE ~ protein
    )
  ) %>% 
  select(c(protein_name, protein)) 
  

#get the protein names and corresponding variable names
clipr::write_clip(protein_rep %>% 
                    pull(protein) %>% 
                    paste(collapse = ", "))

clipr::write_clip(protein_rep %>% 
                    select(protein_name))


#OID00970, OID00739, OID00414, OID00725, OID00694, OID00410, OID00400, NPX_109_IL2_RA, 
#OID00666, OID01015, OID00456, OID00660, OID00684, NPX_178_SCGB3A2, NPX_134_TNFRSF10C, 
#OID00991, NPX_148_TR_AP, OID00732

#make a combined dataset with all proteomics data and other variables
BIOSTAT.PRO <- CVD2 %>% 
  full_join(CVD3,by = "DEM_PT") %>% 
  full_join(IMMUNE,by = "DEM_PT") %>% 
  full_join(ONCO,by = "DEM_PT") %>% 
  select(DEM_PT,OID00970, OID00739, OID00414, OID00725, OID00694, OID00410, OID00400, NPX_109_IL2_RA, 
         OID00666, OID01015, OID00456, OID00660, OID00684, NPX_178_SCGB3A2, NPX_134_TNFRSF10C, 
         OID00991, NPX_148_TR_AP, OID00732) %>% 
  full_join(BIOSTAT_,by = "DEM_PT") %>% 
  filter(MHdays >= 0) %>% 
  mutate(BMI = as.numeric(PE_WGT1/((PE_HGT1*0.01)^2))) %>% 
  mutate(SMOKING3 = factor(case_when(CVRF_SMOKE == 1 ~ "Never",
                                     CVRF_SMOKE == 2 ~ "Past",
                                     CVRF_SMOKE == 3 ~ "Current"))) %>% 
  mutate(ETIO = as.factor(ETIO))

BIOSTAT.PRO$SMOKING3 <- relevel(factor(BIOSTAT.PRO$SMOKING3), ref = "Never")

#relationship between proteins and outcome
skim_without_charts(BIOSTAT.PRO)

#fit a model for each protein, and extract coefficients and p value for protein, and FDR correction
protein_vector <- c(
  "OID00970", "OID00739", "OID00414", "OID00725", "OID00694",
  "OID00410", "OID00400", "NPX_109_IL2_RA", "OID00666",
  "OID01015", "OID00456", "OID00660", "OID00684",
  "NPX_178_SCGB3A2", "NPX_134_TNFRSF10C", "OID00991",
  "NPX_148_TR_AP", "OID00732"
)


#covaraites
adjustment_vars <- c(
  "SMOKING3", "AGE", "DEM_SEX", "BMI", "NPX_NT_PRO_BNP",
  "EGFRBL", "SYMP_NYHA1", "ETIO", "ECHO_LVEF1", "CVRF_NY14"
)

run_cox_extract <- function(protein_name, data = BIOSTAT.PRO) {
  
  vars <- unique(c("MHdays", "MH", adjustment_vars, protein_name))
  
  # Skip proteins or covariates that are not present
  if (!all(vars %in% names(data))) return(NULL)
  
  # Use the same complete-case participants to calculate the SD and fit the model
  model_data <- data[, vars, drop = FALSE]
  model_data <- model_data[complete.cases(model_data), , drop = FALSE]
  
  if (nrow(model_data) == 0) return(NULL)
  
  protein_sd <- sd(model_data[[protein_name]])
  if (!is.finite(protein_sd) || protein_sd == 0) return(NULL)
  
  model_data$protein_z <- as.numeric(
    (model_data[[protein_name]] - mean(model_data[[protein_name]])) / protein_sd
  )
  
  formula <- as.formula(paste(
    "Surv(MHdays, MH) ~",
    paste(c(adjustment_vars, "protein_z"), collapse = " + ")
  ))
  
  fit <- tryCatch(
    coxph(formula, data = model_data),
    error = function(e) NULL
  )
  
  if (is.null(fit)) return(NULL)
  
  coef_sum <- summary(fit)$coefficients["protein_z", ]
  beta <- unname(coef_sum["coef"])
  se   <- unname(coef_sum["se(coef)"])
  
  # Test proportional hazards for the whole model
  ph_test <- tryCatch(
    cox.zph(fit, transform = "km"),
    error = function(e) e
  )
  
  if (inherits(ph_test, "error")) {
    ph_protein_p <- NA_real_
  } else {
    ph_table <- ph_test$table
    
    ph_protein_p <- if ("protein_z" %in% rownames(ph_table)) {
      unname(ph_table["protein_z", "p"])
    } else {
      NA_real_
    }
    
    ph_test_error <- NA_character_
  }
  
  data.frame(
    protein = protein_name,
    n = nrow(model_data),
    SD_NPX = protein_sd,
    beta_per_SD = beta,
    HR_per_SD = exp(beta),
    CI_lower = exp(beta - qnorm(0.975) * se),
    CI_upper = exp(beta + qnorm(0.975) * se),
    p_value = unname(coef_sum["Pr(>|z|)"]),
    PH_protein_p = ph_protein_p,
    row.names = NULL
  )
}

results_protein <- lapply(protein_vector, run_cox_extract)

final_results <- bind_rows(results_protein) %>%
  left_join(protein_rep, by = "protein") %>%
  mutate(FDR = p.adjust(p_value, method = "fdr"))%>% 
  mutate(protein_name = ifelse(protein_name == "CXL17", "CXCL17", protein_name))

write_xlsx(final_results, "protein_HR outcome_index.xlsx")


################################proteins and HR outcomes - validation cohort################################
##########################################################################################
##########################################################################################
##########################################################################################

#make a combined dataset with all proteomics data and other variables
BIOSTAT.PRO.V <- CVD2.V %>% 
  full_join(CVD3.V,by = "DEM_PT") %>% 
  full_join(IMMUNE.V,by = "DEM_PT") %>% 
  full_join(ONCO.V,by = "DEM_PT") %>% 
  select(DEM_PT,OID00970, OID00739, OID00414, OID00725, OID00694, OID00410, OID00400, NPX_109_IL2_RA, 
         OID00666, OID01015, OID00456, OID00660, OID00684, NPX_178_SCGB3A2, NPX_134_TNFRSF10C, 
         OID00991, NPX_148_TR_AP, OID00732) %>% 
  full_join(Validation_,by = "DEM_PT") %>% 
  filter(MHdays >= 0) %>% 
  mutate(BMI = PE_BMI) %>% 
  mutate(SMOKING3 = factor(case_when(CVRF_SMOKE == 1 ~ "Never",
                                     CVRF_SMOKE == 2 ~ "Past",
                                     CVRF_SMOKE == 3 ~ "Current"))) %>% 
  mutate(
    ETIO = factor(case_when(
      MHHF_DIS1 == 1 ~ "Ischemic Heart Disease",
      MHHF_DIS2 == 1 ~ "Hypertensive",
      MHHF_DIS3 == 1 ~ "Cardiomyopathy",
      MHHF_DIS4 == 1 ~ "Valvular Disease",
      MHHF_DIS6 == 1 ~ "Other",
      TRUE ~ "Unknown"
    ))
  )

BIOSTAT.PRO.V$SMOKING3 <- relevel(factor(BIOSTAT.PRO.V$SMOKING3), ref = "Never")

#add EGFR
BIOSTAT.PRO.V <- Validation.lab %>% 
  select(DEM_PT, EGFR) %>% 
  full_join(BIOSTAT.PRO.V, by = "DEM_PT")

skim_without_charts(BIOSTAT.PRO.V)


#fit a model for each protein, and extract coefficients and p value for protein, and FDR correction

adjustment_vars <- c(
  "SMOKING3", "AGE", "DEM_SEX", "BMI", "NPX_NT_PRO_BNP",
  "EGFR", "MHHF_NYHA", "ETIO", "ECRV_LVEF", "CVRF_NY14"
)

run_cox_extract <- function(protein_name, data = BIOSTAT.PRO.V) {
  
  vars <- unique(c("MHdays", "MH", adjustment_vars, protein_name))
  
  # Skip proteins or covariates that are not present
  if (!all(vars %in% names(data))) return(NULL)
  
  # Use the same complete-case participants to calculate the SD and fit the model
  model_data <- data[, vars, drop = FALSE]
  model_data <- model_data[complete.cases(model_data), , drop = FALSE]
  
  if (nrow(model_data) == 0) return(NULL)
  
  protein_sd <- sd(model_data[[protein_name]])
  if (!is.finite(protein_sd) || protein_sd == 0) return(NULL)
  
  model_data$protein_z <- as.numeric(
    (model_data[[protein_name]] - mean(model_data[[protein_name]])) / protein_sd
  )
  
  formula <- as.formula(paste(
    "Surv(MHdays, MH) ~",
    paste(c(adjustment_vars, "protein_z"), collapse = " + ")
  ))
  
  fit <- tryCatch(
    coxph(formula, data = model_data),
    error = function(e) NULL
  )
  
  if (is.null(fit)) return(NULL)
  
  coef_sum <- summary(fit)$coefficients["protein_z", ]
  beta <- unname(coef_sum["coef"])
  se   <- unname(coef_sum["se(coef)"])
  
  # Test proportional hazards for the whole model
  ph_test <- tryCatch(
    cox.zph(fit, transform = "km"),
    error = function(e) e
  )
  
  if (inherits(ph_test, "error")) {
    ph_protein_p <- NA_real_
  } else {
    ph_table <- ph_test$table
    
    ph_protein_p <- if ("protein_z" %in% rownames(ph_table)) {
      unname(ph_table["protein_z", "p"])
    } else {
      NA_real_
    }
    
    ph_test_error <- NA_character_
  }
  
  data.frame(
    protein = protein_name,
    n = nrow(model_data),
    SD_NPX = protein_sd,
    beta_per_SD = beta,
    HR_per_SD = exp(beta),
    CI_lower = exp(beta - qnorm(0.975) * se),
    CI_upper = exp(beta + qnorm(0.975) * se),
    p_value = unname(coef_sum["Pr(>|z|)"]),
    PH_protein_p = ph_protein_p,
    row.names = NULL
  )
}


results_protein <- lapply(protein_vector, run_cox_extract)

final_results <- bind_rows(results_protein) %>%
  left_join(protein_rep, by = "protein") %>%
  mutate(FDR = p.adjust(p_value, method = "fdr")) %>% 
  mutate(protein_name = ifelse(protein_name == "CXL17", "CXCL17", protein_name))

write_xlsx(final_results, "protein_HR outcome_validation.xlsx")

#######some protein models violated proportional hazard assumption##########
#######we check the Schoenfeld residual plot#################
#OID00684 (SCF), OID00666 (IL6), OID00970 (CCL11), OID00456 (MMP12)

fit <- coxph(
  Surv(MHdays, MH) ~ OID00456,
  data = BIOSTAT.PRO.V
)
plot(cox.zph(fit), ylab = "Beta for MMP12")
abline(h = 0, lty = 2)

#######################make forest plots for the protein-HR outcome#######################

index.fp <- read_excel("protein_HR outcome_index.xlsx") %>%
  select(protein_name, beta_per_SD, HR_per_SD, CI_lower, CI_upper) %>%
  rename(
    coef = beta_per_SD,
    HR_index = HR_per_SD,
    lower_index = CI_lower,
    upper_index = CI_upper
  )

val.fp <- read_excel("protein_HR outcome_validation.xlsx") %>%
  select(protein_name, beta_per_SD, HR_per_SD, CI_lower, CI_upper) %>%
  rename(
    coef.val = beta_per_SD,
    HR_val = HR_per_SD,
    lower_val = CI_lower,
    upper_val = CI_upper
  )

fp <- index.fp %>%
  left_join(val.fp, by = "protein_name")

# HRs and 95% CIs are already in the updated results
fp_plot <- fp

# Order proteins by HR in the index cohort
protein_order <- fp_plot %>%
  arrange(HR_index) %>%
  pull(protein_name)

# Convert to long format
forest_data <- bind_rows(
  fp_plot %>%
    transmute(
      protein_name,
      cohort = "Index cohort",
      HR = HR_index,
      lower = lower_index,
      upper = upper_index
    ),
  
  fp_plot %>%
    transmute(
      protein_name,
      cohort = "Validation cohort",
      HR = HR_val,
      lower = lower_val,
      upper = upper_val
    )
) %>%
  mutate(
    label = sprintf("%.2f (%.2f–%.2f)", HR, lower, upper)
  )

forest_data$protein_name <- factor(
  forest_data$protein_name,
  levels = protein_order
)

# Forest plot
p <- ggplot(
  forest_data,
  aes(x = HR, y = protein_name)
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    colour = "grey50"
  ) +
  geom_errorbarh(
    aes(xmin = lower, xmax = upper),
    height = 0.2,
    linewidth = 0.5
  ) +
  geom_point(size = 2.5) +
  geom_text(
    aes(x = HR, label = label),
    nudge_y = 0.35,
    size = 3.5
  ) +
  facet_grid(
    . ~ cohort,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_x_log10(
    breaks = scales::breaks_extended(n = 6)
  ) +
  labs(
    x = "Hazard ratio per 1-SD increase (95% CI)",
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(
      size = 10,
      colour = "black",
      face = "bold"
    ),
    axis.title.x = element_text(face = "bold"),
    axis.text.x = element_text(
      size = 10,
      colour = "black",
      face = "bold"
    ),
    plot.title = element_text(face = "bold")
  )

p


################################protein- protein interaction################################
##########################################################################################
##########################################################################################
##########################################################################################

###############explore the correlations among the replicated proteins (index/validation cohort)#######
#correlations among protein residuals after adjustment for key clinical variables 
#first get protein residuals, then perform correlations


#################################in index cohort#############################
pp_vector <- c(
  "OID00732", "OID00725", "OID01015", "OID00991","NPX_148_TR_AP", "OID00660",
  "OID00410", "OID00666", "OID00456","OID00739","OID00684", "OID00970"
)

#index cohort
adjustment_vars <- c(
  "SMOKING3", "AGE", "DEM_SEX", "BMI", "NPX_NT_PRO_BNP",
  "EGFRBL", "SYMP_NYHA1", "ETIO", "ECHO_LVEF1", "CVRF_NY14"
)

#map Olink variable names to protein names for the output
protein_map <- c(
  OID00732 = "WFDC2",
  OID00725 = "CXCL17",
  OID01015 = "LAMP3",
  OID00991 = "TREM1",
  NPX_148_TR_AP = "TRAP",
  OID00660 = "MSLN",
  OID00410 = "FGF-21",
  OID00666 = "IL6",
  OID00456 = "MMP12",
  OID00739 = "CEACAM5",
  OID00684 = "SCF",
  OID00970 = "CCL11"
)

required_vars <- c(pp_vector, adjustment_vars)


# Use the same complete-case participants for all proteins
dat_complete <- BIOSTAT.PRO %>%
  select(all_of(required_vars)) %>%
  mutate(
    SMOKING3 = factor(SMOKING3),
    DEM_SEX  = factor(DEM_SEX)
  ) %>%
  drop_na()


# Regress each protein on the same adjustment variables
# and save its residuals
residual_matrix <- sapply(pp_vector, function(protein_var) {
  fit <- lm(
    reformulate(adjustment_vars, response = protein_var),
    data = dat_complete
  )
  
  residuals(fit)
})

# Replace OID names with protein names
colnames(residual_matrix) <- unname(protein_map[pp_vector])

# Pearson correlations among protein residuals
partial_cor_matrix <- cor(
  residual_matrix,
  method = "pearson"
)

# Extract each unique protein pair from the upper triangle
pair_index <- which(upper.tri(partial_cor_matrix), arr.ind = TRUE)

pair_results <- data.frame(
  Protein_1 = rownames(partial_cor_matrix)[pair_index[, 1]],
  Protein_2 = colnames(partial_cor_matrix)[pair_index[, 2]],
  Pearson_r = round(partial_cor_matrix[pair_index], 3),
  N = nrow(dat_complete)
)


protein_order <- rownames(partial_cor_matrix)

# Convert the correlation matrix to plotting data
cor_plot_data <- as.data.frame(as.table(partial_cor_matrix)) %>%
  setNames(c("Protein_row", "Protein_column", "r")) %>%
  mutate(
    row_index = match(as.character(Protein_row), rownames(partial_cor_matrix)),
    col_index = match(as.character(Protein_column), colnames(partial_cor_matrix))
  ) %>%
  filter(row_index > col_index) %>%  # one triangle; omit diagonal
  mutate(
    Protein_row = factor(
      as.character(Protein_row),
      levels = rev(rownames(partial_cor_matrix))
    ),
    Protein_column = factor(
      as.character(Protein_column),
      levels = colnames(partial_cor_matrix)
    )
  )

p_cor_matrix <- ggplot(
  cor_plot_data,
  aes(x = Protein_column, y = Protein_row)
) +
  geom_tile(fill = "white", colour = "black", linewidth = 0.4) +
  geom_text(
    aes(label = sprintf("%.2f", r)),
    colour = "black",
    size = 4
  ) +
  coord_fixed() +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
    axis.text.y = element_text(colour = "black"),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    plot.title = element_text(face = "bold", colour = "black"),
    plot.subtitle = element_text(colour = "black")
  )

print(p_cor_matrix)

# Save the protein-pair results in Excel
writexl::write_xlsx(
  list("Protein pairs" = pair_results),
  path = "protein_partial_correlations index.xlsx"
)

#############################in validation cohort###############################

pp_vector <- c(
  "OID00732", "OID00725", "OID00991", "OID00660",
  "NPX_178_SCGB3A2", "OID00666", "OID00739"
)

adjustment_vars <- c(
  "SMOKING3", "AGE", "DEM_SEX", "BMI", "NPX_NT_PRO_BNP",
  "EGFR", "MHHF_NYHA", "ETIO", "ECRV_LVEF", "CVRF_NY14"
)

#map Olink variable names to protein names for the output
protein_map <- c(
  OID00732 = "WFDC2",
  OID00725 = "CXCL17",
  OID00991 = "TREM1",
  OID00660 = "MSLN",
  NPX_178_SCGB3A2 = "SCGB3A2",
  OID00666 = "IL6",
  OID00739 = "CEACAM5"
)


required_vars <- c(pp_vector, adjustment_vars)


# Use the same complete-case participants for all proteins
dat_complete <- BIOSTAT.PRO.V %>%
  select(all_of(required_vars)) %>%
  mutate(
    SMOKING3 = factor(SMOKING3),
    DEM_SEX  = factor(DEM_SEX)
  ) %>%
  drop_na()


# Regress each protein on the same adjustment variables
# and save its residuals
residual_matrix <- sapply(pp_vector, function(protein_var) {
  fit <- lm(
    reformulate(adjustment_vars, response = protein_var),
    data = dat_complete
  )
  
  residuals(fit)
})

# Replace OID names with protein names
colnames(residual_matrix) <- unname(protein_map[pp_vector])

# Pearson correlations among protein residuals
partial_cor_matrix <- cor(
  residual_matrix,
  method = "pearson"
)

# Extract each unique protein pair from the upper triangle
pair_index <- which(upper.tri(partial_cor_matrix), arr.ind = TRUE)

pair_results <- data.frame(
  Protein_1 = rownames(partial_cor_matrix)[pair_index[, 1]],
  Protein_2 = colnames(partial_cor_matrix)[pair_index[, 2]],
  Pearson_r = round(partial_cor_matrix[pair_index], 3),
  N = nrow(dat_complete)
)


protein_order <- rownames(partial_cor_matrix)

# Convert the correlation matrix to plotting data
cor_plot_data <- as.data.frame(as.table(partial_cor_matrix)) %>%
  setNames(c("Protein_row", "Protein_column", "r")) %>%
  mutate(
    row_index = match(as.character(Protein_row), rownames(partial_cor_matrix)),
    col_index = match(as.character(Protein_column), colnames(partial_cor_matrix))
  ) %>%
  filter(row_index > col_index) %>%  # one triangle; omit diagonal
  mutate(
    Protein_row = factor(
      as.character(Protein_row),
      levels = rev(rownames(partial_cor_matrix))
    ),
    Protein_column = factor(
      as.character(Protein_column),
      levels = colnames(partial_cor_matrix)
    )
  )

p_cor_matrix <- ggplot(
  cor_plot_data,
  aes(x = Protein_column, y = Protein_row)
) +
  geom_tile(fill = "white", colour = "black", linewidth = 0.4) +
  geom_text(
    aes(label = sprintf("%.2f", r)),
    colour = "black",
    size = 4
  ) +
  coord_fixed() +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
    axis.text.y = element_text(colour = "black"),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    plot.title = element_text(face = "bold", colour = "black"),
    plot.subtitle = element_text(colour = "black")
  )

print(p_cor_matrix)

# Save the protein-pair results in Excel
writexl::write_xlsx(
  list("Protein pairs" = pair_results),
  path = "protein_partial_correlations val.xlsx"
)


# the protein-protein network was then made with CytoSPACE










