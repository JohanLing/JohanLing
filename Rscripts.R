
rm(list=ls())

#Date: 2024/5/20
#Time:11:35
#By: Ling

setwd("C:/JohanLing/Projects/Loasaceae/Mentzelia Section Bartonia/Environmental factors/Analysis/PCA_PGLS")###Load the libraries####
library(terra)
library(dplyr)
library(vegan)
################################STEP 5: pPCA ANALYSES#####################
################################PREPARING REQUIRED FILES##########################
#################################################################################
#Load non correlated environmental variables: file name = non_correlated bioclim layers.tif
file_path <- "non_correlated_bioclim_layers.tif"
#Load the raster date
non_correlated_stack <- rast(file_path)
#Load distribution file containing taxon, staminode,gypsum, and coordinates
dist <- read.csv(file="Distributions.20240523Cleaned.csv")

###plotting and saving species distribution on bioclimate layer
#Open the SVG device to start capturing plot output
svg("species_distribution_on_bioclimate_layer.svg")
#Basic plot of the first layer to use as a background
plot(non_correlated_stack[[1]], main="Species distribution on bioclimate layer")
#Add the point DISTRIBUTION DATA on top of the plot
points(dist$Longitude, dist$Latitude, col="black", pch=16)
#Close the SVG device to save the plot to the file
dev.off()

####Converting point data to a SpatVector and extract climate data####
points_vector <- vect(dist, geom=c("Longitude", "Latitude"), crs=crs(non_correlated_stack))
extracted_data <- terra::extract(non_correlated_stack, points_vector)

#Names of the environmental variables corresponding to the selected non-correlated layers
#If the first column in 'dist' is an unwanted ID column
extracted_data <- extracted_data[, -1]  # Removes the first column
colnames(extracted_data) <- c("bio1", #Annual_Mean_Temp
  "bio2", # Mean_Diurnal_Range
  "bio4", # Temp_Seasonality
  "bio8", # Mean_Temp_Wettest_Quarter
  "bio9", # Mean_Temp_Driest_Quarter
  "bio12", #Annual_Precip
  "bio14", # Precip_Driest_Month
  "bio15", # Precip_Seasonality
  "bio18", # Precip_Warmest_Quarter
  "bio21") # Elevation

#Add Taxon, Staminode type, gypsum type and coordinate columns to the selected non-correlated layers
extracted_data$Taxon <- dist$Taxon
extracted_data$Taxon_ID <- dist$Taxon_ID
extracted_data$Floral_structure <- dist$Floral_structure
extracted_data$Gypsum_type <- dist$Gypsum_type
extracted_data$Latitude <- dist$Latitude
extracted_data$Longitude <- dist$Longitude
extracted_data$Perennial <- dist$Perennial
extracted_data$Country <- dist$Country

#View the first few rows of the extracted data with named columns
head(extracted_data)

############CHECKING NA AND DUPLICATE######################
####Identify rows with any NA values across all columns
rows_with_na <- apply(is.na(extracted_data), 1, any)
#Display row indices that have NA values
na_row_indices <- which(rows_with_na)
#Display rows that contain NA values
na_rows <- extracted_data[na_row_indices, ]
#Print the rows with NAs
print(na_rows)
###Remove rows with NAs from the extracted_data dataframe
extracted_data_na_rem <- na.omit(extracted_data)
#View the cleaned extracted data
head(extracted_data_na_rem)

####Identify and remove duplicate coordinates within each taxon
extracted_data_dup <- extracted_data_na_rem %>%
  group_by(Taxon, Longitude, Latitude, Country, Gypsum_type, Perennial, Floral_structure) %>%
  distinct() %>%
  ungroup()
head(extracted_data_dup)
#verify no duplicates after removal
verify <- extracted_data_dup %>%
  group_by(Taxon, Longitude, Latitude, Country, Gypsum_type, Perennial, Floral_structure) %>%
  filter(n() > 1) %>%
  ungroup()
#Check if there are any duplicates left
if(nrow(verify) > 0) {
  print("Duplicates still exist within taxa after removal.")
} else {
  print("No duplicates exist within taxa after removal.")
}


##############there are two rows of elevation with -36 and -1#################
############################## detete them #################################
#there are negative values in Bi021, check this.use 4 m instead of -1
#Save the duplicate data as a CSV file
#write.csv(extracted_data_dup, file = "extracted_data_dup_withMinus.csv", row.names = FALSE)
###Load the CSV file
extracted_data_dup <- read.csv("extracted_data_dup_withMinus.csv")

################FILTER MENLONLON##############
mlon <- extracted_data_dup %>% filter(Taxon == "MENLONLON")
#GENERATE DENSITY PLOT FOR BIO21
library(ggplot2)
ggplot(mlon, aes(x = bio21)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = "Density Plot of bio21 for MENLONLON", x = "bio21", y = "Density") +
  theme_minimal()

mlon1 <- mlon %>% filter(bio21 != -36)
ggplot(mlon1, aes(x = bio21)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = "Density Plot of bio21 for MENLONLON", x = "bio21", y = "Density") +
  theme_minimal()
###############calculate sample size for each country#############
##################################################################
###OUNTING UNIQUE TAXA AND INDIVIDUALS PER COUNTRY###

#Count the number of individuals for each country after removing duplicates and NAs
individuals_by_country <- extracted_data_dup %>%
  group_by(Country) %>%
  summarise(Num_Individuals = n())

#Count the number of unique taxa for each country after removing duplicates and NAs
taxa_by_country <- extracted_data_dup %>%
  group_by(Country) %>%
  summarise(Unique_Taxa = n_distinct(Taxon))

#Merge the two dataframes to get a complete summary
summary_by_country <- merge(individuals_by_country, taxa_by_country, by = "Country")

#Print the summary
print(summary_by_country)



#Count the total number of individuals after removing duplicates and NAs
total_individuals <- nrow(extracted_data_dup)

#Count the total number of unique taxa across all countries after removing duplicates and NAs
total_taxa <- n_distinct(extracted_data_dup$Taxon)

#Print the total counts
cat("Total number of individuals:", total_individuals, "\n")
cat("Total number of unique taxa:", total_taxa, "\n")

#Count the number of individuals for each species (Taxon)
sample_size_by_species <- extracted_data_dup %>%
  group_by(Taxon) %>%
  summarise(Sample_Size = n()) %>%
  arrange(desc(Sample_Size))  # Arrange by sample size in descending order

#Print the sample size difference across species

#Calculate the sample size for each species (Taxon)
sample_size_by_species <- extracted_data_dup %>%
  group_by(Taxon) %>%
  summarise(Sample_Size = n()) %>%
  arrange(desc(Sample_Size))  #Arrange by sample size in descending order

#Calculate the mean sample size across species
mean_sample_size <- mean(sample_size_by_species$Sample_Size)

#Print the mean sample size
cat("Mean sample size across species:", mean_sample_size, "\n")


#Count the number of individuals and unique taxa for Staminode presence
staminode_presence <- extracted_data_dup %>%
  filter(Floral_structure == "Staminode") %>%
  summarise(Total_Individuals = n(), Unique_Taxa = n_distinct(Taxon))

#Count the number of individuals and unique taxa for Staminode absence
staminode_absence <- extracted_data_dup %>%
  filter(Floral_structure == "Stamen") %>%
  summarise(Total_Individuals = n(), Unique_Taxa = n_distinct(Taxon))

#Print the results
cat("Staminode Presence - Total Individuals:", staminode_presence$Total_Individuals, 
    " | Unique Taxa:", staminode_presence$Unique_Taxa, "\n")

cat("Staminode Absence - Total Individuals:", staminode_absence$Total_Individuals, 
    " | Unique Taxa:", staminode_absence$Unique_Taxa, "\n")


# Count the number of individuals and unique taxa for Gypsum affinity (gypsum-adapted taxa)
gypsum_adapted <- extracted_data_dup %>%
  filter(Gypsum_type == "gypsum") %>%
  summarise(Total_Individuals = n(), Unique_Taxa = n_distinct(Taxon))

# Count the number of individuals and unique taxa for non-Gypsum taxa
non_gypsum_adapted <- extracted_data_dup %>%
  filter(Gypsum_type == "nongypsum") %>%
  summarise(Total_Individuals = n(), Unique_Taxa = n_distinct(Taxon))

# Print the results
cat("Gypsum-Adapted - Total Individuals:", gypsum_adapted$Total_Individuals, 
    " | Unique Taxa:", gypsum_adapted$Unique_Taxa, "\n")

cat("Non-Gypsum-Adapted - Total Individuals:", non_gypsum_adapted$Total_Individuals, 
    " | Unique Taxa:", non_gypsum_adapted$Unique_Taxa, "\n")

# Count the number of individuals and unique taxa for Perennials
perennials <- extracted_data_dup %>%
  filter(Perennial == "Yes") %>%
  summarise(Total_Individuals = n(), Unique_Taxa = n_distinct(Taxon))

# Count the number of individuals and unique taxa for Non-Perennials
non_perennials <- extracted_data_dup %>%
  filter(Perennial == "No") %>%
  summarise(Total_Individuals = n(), Unique_Taxa = n_distinct(Taxon))

# Print the results
cat("Perennials - Total Individuals:", perennials$Total_Individuals, 
    " | Unique Taxa:", perennials$Unique_Taxa, "\n")

cat("Non-Perennials - Total Individuals:", non_perennials$Total_Individuals, 
    " | Unique Taxa:", non_perennials$Unique_Taxa, "\n")

#################################################################################
###CALCULATING LOG GEOMETRIC MEANS FOR EACH ENVIRONMENTAL VARIABLE WITHIN EACH TAXON#####
#################################################################################
library(tidyverse) #For pivot_longer, ggplot2, and dplyr
#Pivot extracted_data_dup
data_long <- pivot_longer(
  extracted_data_dup,
  cols = c("bio1", #Annual_Mean_Temp
            "bio2", # Mean_Diurnal_Range
            "bio4", # Temp_Seasonality
            "bio8", # Mean_Temp_Wettest_Quarter
            "bio9", # Mean_Temp_Driest_Quarter
            "bio12", # Annual_Precip
            "bio14", # Precip_Driest_Month
            "bio15", # Precip_Seasonality
            "bio18", # Precip_Warmest_Quarter
            "bio21"),#Elevation
  names_to = "Env_variable",
  values_to = "Value"
)

#Define the offset
offset <- 37  #We got -36.00, so I will use 37 to get the positive values here
#Adjust values by adding the new offset
data_long <- data_long %>%
  mutate(Adjusted_Value = Value + offset)
#Calculate the geometric means with log transformation and drop offset to make it original data
geometric_mean_values <- data_long %>%
  group_by(Taxon, Env_variable) %>%
  summarise(
    Geometric_Mean = exp(mean(log(Adjusted_Value), na.rm = TRUE)) - offset,
    .groups = 'drop'
  )

###OR######
#no offset
#Define the offset to convert all negative values to positive
#min_value <- min(data_long$Value, na.rm = TRUE)
#offset <- abs(min_value) + 1  # Offset is abs(-36) + 1 = 37

#Convert values to positive by adding the offset
#data_long <- data_long %>%
#  mutate(Adjusted_Value = Value + offset)

#Apply log transformation to the adjusted positive values
#data_long <- data_long %>%
 # mutate(Log_Adjusted_Value = log(Adjusted_Value))

#Calculate the geometric mean using the log-transformed values
#geometric_mean_values <- data_long %>%
#  group_by(Taxon, Env_variable) %>%
 # summarise(
 #   Geometric_Mean = exp(mean(Log_Adjusted_Value, na.rm = TRUE)),  #No offset subtraction
 #   .groups = 'drop'
 # )

#Display the geometric mean values
#print(geometric_mean_values)


#Calculate the geometric mean without explicit log transformation
#geometric_mean_values_nolog <- data_long %>%
#  group_by(Taxon, Env_variable) %>%
#  summarise(
#    Geometric_Mean = (prod(Adjusted_Value, na.rm = TRUE))^(1/n()),  # Using the product and nth root directly
#    .groups = 'drop'
  #)

#Subtract the offset after geometric mean calculation
#geometric_mean_values_nolog <- geometric_mean_values_nolog %>%
#  mutate(Geometric_Mean = Geometric_Mean - offset)
#Display the geometric mean values
#print(geometric_mean_values_nolog)



#Pivot the data wider using geometric_mean_values
geometric_mean_values_wide <- geometric_mean_values %>%
pivot_wider(names_from = Env_variable, values_from = Geometric_Mean)

#Create a summary table of taxa with their respective Staminode_type and Gypsum_type
library(dplyr)
#Convert to tibble and create a summary table
geo_taxon_types <- extracted_data_dup %>%
  as_tibble() %>%
  dplyr::select(Taxon, Floral_structure, Gypsum_type, Perennial) %>%
  distinct()

geo_taxon_types <- dplyr::select(extracted_data_dup, Taxon, Floral_structure, Gypsum_type, Perennial) %>%
  distinct()
#Join this summary table with the geometric_mean_values_wide to add the Staminode_type and Gypsum_type information back
geo_cbind_taxon_env <- geometric_mean_values_wide %>%
  left_join(geo_taxon_types, by = "Taxon")

write.csv(geo_cbind_taxon_env, file = "geo_cbind_taxon_env.csv", row.names = FALSE)
write.csv(geo_cbind_taxon_env, file = "geo_cbind_taxon_env_withoffset.csv", row.names = FALSE)


####################Phylogenetic PCA####################################
####FIRST: PERFORMING PCA ON GEOMETRIC MEAN DATASET#####################

#Load Necessary Libraries
library(dplyr)
library(ggplot2)
library(caper)
library(phytools)
library(ape)
library(geiger)
library(phylolm)

#Load the Data and Phylogenetic Tree
data <- read.csv("geo_cbind_taxon_env_plus_trait.csv", stringsAsFactors = FALSE)

head(data)


##########################################################
######################NPP#################################
################## Extract Average_NPP Data ##############
##########################################################
#Extract the necessary columns for analysis (Average_NPP and Staminode type)
npp_data <- data %>%
  dplyr::select(Taxon, Floral_structure, Showiness, mas_size, petal_size, Habitat_type, Life_history, Average_NPP) %>%
  filter(!is.na(Average_NPP))

#Convert to data frame and set row names
npp_data <- as.data.frame(npp_data)
rownames(npp_data) <- npp_data$Taxon


#################################################################
################### Check Normality for NPP######################
#################################################################
#Q-Q Plot for Visual Inspection
qq_plot <- ggplot(npp_data, aes(sample = Average_NPP)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "Q-Q Plot for Average_NPP", x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal()

print(qq_plot)

#Shapiro-Wilk Test for Normality
shapiro_test <- shapiro.test(npp_data$Average_NPP)
print(shapiro_test)

#################################################################
################### Check Normality for petal size###############
#################################################################
#Q-Q Plot for Visual Inspection
qq_plot_petal <- ggplot(npp_data, aes(sample = petal_size)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "Q-Q Plot for Petal Size", x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal()

print(qq_plot_petal)

#Shapiro-Wilk Test for Normality
shapiro_test <- shapiro.test(npp_data$petal_size)
print(shapiro_test)

#################################################################
################### Check Normality for mas size###############
#################################################################
#Q-Q Plot for Visual Inspection
qq_plot_mas <- ggplot(npp_data, aes(sample = mas_size)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "Q-Q Plot for MAS Size", x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal()

print(qq_plot_mas)

#Shapiro-Wilk Test for Normality
shapiro_test <- shapiro.test(npp_data$mas_size)
print(shapiro_test)

#PETAL, MAS are not normally distributed so do log transformation for checking the model fit
npp_data$log_petal_size <- log(npp_data$petal_size)
npp_data$log_mas_size <- log(npp_data$mas_size)
#Shapiro-Wilk Test for Normality
log_petal_size <- shapiro.test(npp_data$log_petal_size)
print(log_petal_size)
log_mas_size <- shapiro.test(npp_data$log_mas_size)
print(log_mas_size)
#still not improving
################################################################
####################RUN LM PETAL AND MAS MODEL FIT############
################################################################
###Staminode and stamen combined##########
#Define the GLM model with mas_size as the response variable and petal_size as the predictor
glm_model <- glm(mas_size ~ petal_size, data = npp_data,
                family = Gamma(link = "log"))
#Print the summary of the gm model
glm_summary <- summary(glm_model)
print(glm_summary)

####Staminode alone########
glm_staminode <- glm(mas_size ~ petal_size, 
                     data = subset(npp_data, Floral_structure == "Staminode"), 
                     family = Gamma(link = "log"))
summary_staminode <- summary(glm_staminode)
print(summary_staminode)


####Stamen alone########
glm_stamen <- glm(mas_size ~ petal_size, 
                  data = subset(npp_data, Floral_structure == "Stamen"), 
                  family = Gamma(link = "log"))
summary_stamen <- summary(glm_stamen)
print(summary_stamen)




###Find the residuals for both of staminodes and stamens
#Diagnostic plots to check residuals
par(mfrow = c(2, 2))
plot(glm_model)
#Extract residuals
residuals_glm <- residuals(glm_model)
#Create a data frame with Taxon and corresponding residuals
residuals_df <- data.frame(Taxon = rownames(npp_data), showiness_res = residuals_glm)
#Merge the residuals data frame with the original npp_data on the Taxon column
npp_data <- merge(npp_data, residuals_df, by = "Taxon")
#View the first few rows to confirm the new column is added
head(npp_data)
#Calculate correlation between log-transformed petal size and MAS size using Pearson method
correlation <- cor(npp_data$petal_size, npp_data$mas_size, method = "pearson")
#Extract the p-value for the slope coefficient (petal_size) from the linear model summary
p_value <- glm_summary$coefficients[2, 4]  #[2, 4] corresponds to the slope's p-value in the coefficients matrix
#Print correlation and p-value
print(paste("Correlation between petal size and MAS size:", correlation))
print(paste("P-value for the slope coefficient:", p_value))


#Create the density plot for residuals with a vertical line at zero
#################################################################
################### Check Normality for residuals################
#################################################################
#Q-Q Plot for Visual Inspection
qq_plot_res <- ggplot(npp_data, aes(sample = showiness_res)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "Q-Q Plot for Showiness", x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal()

print(qq_plot_res)

#Shapiro-Wilk Test for Normality
shapiro_test <- shapiro.test(npp_data$showiness_res)
print(shapiro_test)

npp_data$log_showiness_res <- log(npp_data$showiness_res+3)
shapiro_test <- shapiro.test(npp_data$log_showiness_res)
print(shapiro_test)

#Load necessary library
library(ggplot2)
residual_density_plot <- ggplot(npp_data, aes(x = showiness_res)) +
  geom_density(fill = "black", alpha = 0.5) +
  geom_vline(xintercept = 0, color = "black", linetype = "solid") +
  labs(x = "Residuals", y = "Density") +
  theme_minimal() +
  theme(
    text = element_text(family = "Times New Roman", size = 22),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.title = element_text(family = "Times New Roman", size = 22, color = "black"),
    axis.text = element_text(family = "Times New Roman", size = 22, color = "black")
  )
#Print the plot
print(residual_density_plot)
#Save the plot as SVG
ggsave("residual_density_plot.svg", plot = residual_density_plot, device = "svg", bg = "white")
#write.csv(npp_data, file = "geo_cbind_taxon_env_with_MASresidual.csv", row.names = FALSE)



#Create the density plot with curve lines
residual_density_plot <- ggplot(npp_data, aes(x = showiness_res, color = Floral_structure, linetype = Floral_structure)) +
  geom_density(size = 1.5) +
  geom_vline(xintercept = 0, color = "black", linetype = "solid") +
  scale_color_manual(values = c("Stamen" = "orange", "Staminode" = "#377EB8")) +  
  scale_linetype_manual(values = c("Stamen" = "solid", "Staminode" = "dashed")) +  
  labs(x = "Residuals", y = "Density", color = "Floral organ", linetype = "Floral organ") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  #Remove major grid lines
    panel.grid.minor = element_blank(),  #Remove minor grid lines
    axis.line.x = element_line(color = "black", size = 0.5),  #Add x-axis line
    axis.line.y = element_line(color = "black", size = 0.5),  #Add y-axis line
    axis.ticks = element_line(color = "black", size = 1),  #Add axis ticks in black
    axis.ticks.length = unit(0.25, "cm"),  # Length of the ticks
    panel.border = element_rect(color = "black", fill = NA, size = 2),  #Add border around the plot
    axis.title = element_text(family = "Times New Roman", size = 22, color = "black"),
    legend.title = element_text(family = "Times New Roman", size = 24, color = "black"),#Axis title style
    axis.text = element_text(family = "Times New Roman", size = 22, color = "black"),  #Axis text style
    legend.text = element_text(family = "Times New Roman", size = 22, color = "black"),
    legend.position = "none"  #Hide the legend
  )
#Print the plot
print(residual_density_plot)
#Save the plot as an SVG file
ggsave(filename = "Residuals with staminode and stamen.svg", plot = residual_density_plot, device = "svg", width = 9, height = 7, units = "in")


################################################################
################Phylogenetic Tree Preparation###################
################################################################
###Load necessary libraries
library(phytools)
library(ape)
library(geiger)
tree <- read.tree("BartoniaIngroupPLnames.tre")
plot(tree)

#Ensure that taxa in npp_data match those in the tree
npp_data <- npp_data[npp_data$Taxon %in% tree$tip.label, ]


######Staminode vs NPP#########
###############################
#npp_data_filtered1 <- npp_data_filtered %>% mutate(Presence = ifelse(Staminode_type == "Presence", 1,0))
#Convert Staminode_type to numeric if not already done in your dataset
#npp_data <-npp_data$Staminode_type <- as.numeric(as.factor(npp_data$Staminode_type))

library(caper)
comp_data <- comparative.data(tree, npp_data, names.col = "Taxon")


#MAS Showiness#both staminode and stamens
pgls_showiness <- pgls(showiness_res ~ Average_NPP + Floral_structure, data = comp_data)
summary(pgls_showiness)



#Staminode#Run pANOVA
library(nlme)
#Perform Phylogenetic ANOVA
#Grouping variable: Floral_organ
pANOVA_result <- phylANOVA(tree, x = npp_data$Floral_structure, y = npp_data$Average_NPP)
#Check the tree tip labels
print(tree$tip.label)

#Check the dataset Taxon column
print(npp_data$Taxon)

#Identify mismatches
mismatched_taxa <- setdiff(tree$tip.label, npp_data$Taxon)
cat("Taxa in tree but not in data:", mismatched_taxa, "\n")

mismatched_taxa_data <- setdiff(npp_data$Taxon, tree$tip.label)
cat("Taxa in data but not in tree:", mismatched_taxa_data, "\n")
#Subset npp_data to include only taxa in the tree
npp_data <- npp_data[npp_data$Taxon %in% tree$tip.label, ]

#Verify that all taxa now match
all(npp_data$Taxon %in% tree$tip.label) # Should return TRUE

#Perform phylogenetic ANOVA
pANOVA_result <- phylANOVA(tree, x = npp_data$Floral_structure, y = npp_data$Average_NPP)

#Display the results
print(pANOVA_result)

#####To avoid warning####
#########################
#Explicitly add labels to y and x
y_named <- setNames(npp_data$Average_NPP, npp_data$Taxon)
x_named <- setNames(npp_data$Floral_structure, npp_data$Taxon)

#Perform phylogenetic ANOVA with labeled variables
pANOVA_result <- phylANOVA(tree, x = x_named, y = y_named)

#Display the results
print(pANOVA_result)


######run PGLS only for taxa with staminode###########
#Extract taxa with staminodes
taxa_with_staminodes <- data$Taxon[data$Floral_structure == "Staminode"]
#Prune the tree to include only taxa with staminodes
pruned_tree <- drop.tip(tree, setdiff(tree$tip.label, taxa_with_staminodes))
#Subset the data to include only taxa with staminodes
pruned_data <- data[data$Taxon %in% taxa_with_staminodes, ]

#Log-transform Average_NPP and add it as a new column
#pruned_data$log_Average_NPP <- log(pruned_data$Average_NPP)

#Ensure row names match tip labels for compatibility with PGLS
rownames(pruned_data) <- pruned_data$Taxon
#Prepare the comparative data object
comp_data <- comparative.data(pruned_tree, pruned_data, names.col = "Taxon")
#Run PGLS
pgls_model <- pgls(showiness_res ~ Average_NPP, data = comp_data)
#Display the summary of the PGLS model
summary(pgls_model)

#Plot diagnostic plots for residuals
par(mfrow = c(2, 2))
plot(pgls_model)




######run PGLS only for taxa without staminode###########
#Extract taxa with staminodes
taxa_without_staminodes <- data$Taxon[data$Floral_structure == "Stamen"]
#Prune the tree to include only taxa without staminodes
pruned_tree1 <- drop.tip(tree, setdiff(tree$tip.label, taxa_without_staminodes))
#Subset the data to include only taxa without staminodes
pruned_data1 <- data[data$Taxon %in% taxa_without_staminodes, ]
#Ensure row names match tip labels for compatibility with PGLS
rownames(pruned_data1) <- pruned_data1$Taxon
#Prepare the comparative data object
comp_data1 <- comparative.data(pruned_tree1, pruned_data1, names.col = "Taxon")
#Run PGLS
pgls_model1 <- pgls(showiness_res ~ Average_NPP, data = comp_data1)
#Display the summary of the PGLS model
summary(pgls_model1)

#Plot diagnostic plots for residuals
par(mfrow = c(2, 2))
plot(pgls_model1)


###########Run PGLMM for Petal size (mm)############
#calculate overdispersion 
#Fit the Poisson model first
poisson_model <- glm(log_petal_size ~ Average_NPP + Floral_structure, 
                     data = npp_data, 
                     family = poisson())

#Calculate the overdispersion statistic
overdispersion_stat <- sum(residuals(poisson_model, type = "pearson")^2) / poisson_model$df.residual
#Print the overdispersion statistic
print(overdispersion_stat)#greater than 1 means overdispersion

library(phyr)
#Run PGLMM model including observation-level random effect
pglmm_model_petal <- pglmm(
  log_petal_size ~ Average_NPP + Floral_structure + (1 | log_mas_size), #Added observation-level random effect
  data = npp_data,
  family = "gaussian",  #Use Gaussian if data is continuous and normally distributed
  cov_ranef = list(Taxon = tree)  #Include phylogenetic tree as a random effect
)
summary(pglmm_model_petal)



#########################################################
##############mean and standard error NPP################
#########################################################
#Load necessary library
library(dplyr)
#Calculate the mean value of Average_NPP for staminode presence and absence
mean_npp_staminode <- npp_data %>%
  group_by(Floral_structure) %>%
  summarise(
    Mean_Average_NPP = mean(Average_NPP, na.rm = TRUE),
    SE_Average_NPP = sd(Average_NPP, na.rm = TRUE) / sqrt(n()) #Standard Error
  )
#Print the results
print(mean_npp_staminode)

#########################################################
##########mean and standard error MAS showiness##########
#########################################################
#Load necessary library
library(dplyr)
#Calculate the mean and standard error of flower showiness for staminode presence and absence
mean_masshowiness_staminode <- npp_data %>%
  group_by(Floral_structure) %>%
  summarise(
    Mean_masShowiness = mean(showiness_res, na.rm = TRUE),
    SE_masShowiness = sd(showiness_res, na.rm = TRUE) / sqrt(n())  # Standard Error
  )

# Print the results
print(mean_masshowiness_staminode)

#########################################################
########mean and standard error Petal size###############
#########################################################
#Load necessary library
library(dplyr)
#Calculate the mean and standard error of flower showiness for staminode presence and absence
mean_petalshowiness_staminode <- npp_data %>%
  group_by(Floral_structure) %>%
  summarise(
    Mean_petalShowiness = mean(petal_size, na.rm = TRUE),
    SE_petalShowiness = sd(petal_size, na.rm = TRUE) / sqrt(n())  # Standard Error
  )

#Print the results
print(mean_petalshowiness_staminode)



###################################################################
##############GGPLOT2: NPP FOR MAS SHOWINESS#######################
###################################################################
#Load necessary library
library(ggplot2)
library(grid)

theme <-  theme(
  panel.grid.major = element_blank(),  #Remove major grid lines
  panel.grid.minor = element_blank(),  #Remove minor grid lines
  axis.line.x = element_line(color = "black", size = 0.5),  #Add x-axis line
  axis.line.y = element_line(color = "black", size = 0.5),  #Add y-axis line
  axis.ticks = element_line(color = "black", size = 1),  #Add axis ticks in black
  axis.ticks.length = unit(0.25, "cm"),  # Length of the ticks
  panel.border = element_rect(color = "black", fill = NA, size = 2),  #Add border around the plot
  axis.title = element_text(family = "Times New Roman", size = 22, color = "black"),  #Axis title style
  legend.title = element_text(family = "Times New Roman", size = 24, color = "black"),#Axis title style
  axis.text = element_text(family = "Times New Roman", size = 22, color = "black"),  #Axis text style
  legend.text = element_text(family = "Times New Roman", size = 22, color = "black"),
  legend.position = "none"  # Hide the legend
)

##########################mas showiness#########################################
#############################ggplot2############################################
#Create the boxplot with jittered points
npp_plot <- ggplot(npp_data, aes(x = Floral_structure, y = Average_NPP, fill = Floral_structure)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, color = "black") +  # Boxplot with no outlier points
  geom_jitter(aes(shape = Floral_structure, color = Floral_structure), size = 8, width = 0.2) +  #Jittered points with specific shape
  scale_fill_manual(name = "Floral organ",  # Rename the legend title for fill
    values = c("Stamen" = "darkorange", "Staminode" = "#377EB8")) +
  scale_color_manual(name = "Floral organ",  # Rename the legend title for color
    values = c("Stamen" = "darkorange", "Staminode" = "#377EB8")) +
  scale_shape_manual(name = "Floral organ",  # Rename the legend title for shape
    values = c("Stamen" = 21, "Staminode" = 17)) +
  labs(x = "Taxa", y = expression("Mean net primary productivity (g C m"^-2~yr^-1*")")) +  #Labels for x and y axes
  theme_minimal() +  # Minimal theme to start with
  theme + scale_x_discrete(labels = c("Without staminode", "With staminode")) +
  annotate("text",x = 0.5,y = max(npp_data$Average_NPP) + 0.05,label = "(a)",hjust = 0,vjust = 1, size = 8, family = "Times New Roman")

#Print the plot
print(npp_plot)
#Save the plot as an SVG file
ggsave(filename = "npp_plot.svg", plot = npp_plot, device = "svg", width = 9, height = 7, units = "in")

############
#Load necessary library
library(ggplot2)
theme1 <-  theme(
  panel.grid.major = element_blank(),  #Remove major grid lines
  panel.grid.minor = element_blank(),  #Remove minor grid lines
  axis.line.x = element_line(color = "black", size = 0.5),  #Add x-axis line
  axis.line.y = element_line(color = "black", size = 0.5),  #Add y-axis line
  axis.ticks = element_line(color = "black", size = 1),  #Add axis ticks in black
  axis.ticks.length = unit(0.25, "cm"),  # Length of the ticks
  panel.border = element_rect(color = "black", fill = NA, size = 2),  #Add border around the plot
  axis.title = element_text(family = "Times New Roman", size = 22, color = "black"),
  legend.title = element_text(family = "Times New Roman", size = 24, color = "black"),#Axis title style
  axis.text = element_text(family = "Times New Roman", size = 22, color = "black"),  #Axis text style
  legend.text = element_text(family = "Times New Roman", size = 22, color = "black"),
  legend.position = "none"  #Hide the legend
)
#Define the shapes for Staminode presence and absence
staminode_shapes <- c("Staminode" = 17, "Stamen" = 16)  #16 for circle, 17 for triangle
#Create a scatter plot with a regression line, differentiated by Floral_structure with specified colors and shapes
npp_showiness_plot <- ggplot(npp_data, aes(x = Average_NPP, y = showiness_res, color = Floral_structure, shape = Floral_structure)) +
  geom_smooth(method = "lm", aes(group = Floral_structure), size = 5, width = 0.2) +  #Add a linear regression line with confidence interval
  geom_point(size = 8) +  #Set the point size
  scale_color_manual(name = "Floral organ", values = c("Staminode" = "#377EB8", "Stamen" = "darkorange")) +  #Custom colors and legend title
  scale_shape_manual(name = "Floral organ", values = staminode_shapes) +  # Set custom shapes
  labs(x = expression("Mean net primary productivity (g C m"^-2~yr^-1*")"), y = "Showiness (residuals)") +  #Label the axes
  ylim(-1.6, 1.6) +  # Set the y-axis limits
  theme_minimal() +  #Apply a minimal theme
  theme1 + 
  annotate("text", x = min(npp_data$Average_NPP), y = max(npp_data$showiness_res) + 0.28,
           label = "(a)", hjust = 0, vjust = 1, size = 8, family = "Times New Roman")
print(npp_showiness_plot)
ggsave(filename = "this is for LEgend.svg", plot = npp_showiness_plot, device = "svg", width = 9, height = 7, units = "in")


##################################Petal showiness##############################
###############################################################################

#Define the shapes for Staminode presence and absence
staminode_shapes <- c("Stamen" = 16, "Staminode" = 17)  #16 for circle, 17 for triangle
#Create a scatter plot with a regression line, differentiated by Floral_structure with specified colors and shapes
npp_petal <- ggplot(npp_data, aes(x = Average_NPP, y = log_petal_size, color = Floral_structure, shape = Floral_structure)) +
  geom_smooth(method = "lm", aes(group = Floral_structure), size = 5, width = 0.2) +  #Add a linear regression line with confidence interval
  geom_point(size = 8) +  #Set the point size
  scale_color_manual(name = "Floral organ", values = c("Staminode" = "#377EB8", "Stamen" = "darkorange")) +  #Custom colors and legend title
  scale_shape_manual(name = "Floral organ", values = staminode_shapes) +  # Set custom shapes
  labs(x = expression("Mean net primary productivity (g C m"^-2~yr^-1*")"), y = "Log petal showiness (mm²)") +  #Label the axes
  theme_minimal() +  #Minimal theme to start with
  theme1 +
  annotate("text", x = min(npp_data$Average_NPP), y = max(npp_data$showiness_res) + 6, label = "(a)", 
           hjust = 0, vjust = 1, size = 8, family = "Times New Roman")  
print(npp_petal)
#ggsave(filename = "this is for LEgend for petal.svg", plot = npp_petal, device = "svg", width = 9, height = 7, units = "in")

##############################################################
########PHYLOGENETIC PCA FOR ENVIRONMENTAL VARIABLES##########
##############################################################
#Select only the environmental variables (columns starting with "bio")
env_data <- dplyr::select(data, starts_with("bio"))
#Standardize the environmental variables (mean = 0, variance = 1)
env_data_scaled1 <- scale(env_data)
#Create a new data frame with scaled data and Taxon as the first column
env_data_scaled1 <- as.data.frame(env_data_scaled1) %>%
  mutate(Taxon = data$Taxon)  # Adding the Taxon column
#Set the row names to be Taxon for matching
rownames(env_data_scaled1) <- env_data_scaled1$Taxon
#Check and match species between the data and the phylogenetic tree
#Get species names from tree and scaled data
tree_species <- tree$tip.label
data_species <- rownames(env_data_scaled1)
pca_result <- phyl.pca(tree, env_data_scaled1[, -ncol(env_data_scaled1)], method = "BM", mode = "corr")
#Get summary of PCA result
pca_summary <- summary(pca_result)
#Print the proportion of variance explained by each PC
print(pca_summary$importance)

###Extract the proportion of variance for the first three PCs
proportion_variance <- pca_summary$importance[2, 1:3] * 100

###Print the percentage of variance explained by PC1, PC2, and PC3
cat("Percentage of Variance Explained by PC1:", round(proportion_variance[1], 2), "%\n")
cat("Percentage of Variance Explained by PC2:", round(proportion_variance[2], 2), "%\n")
cat("Percentage of Variance Explained by PC3:", round(proportion_variance[3], 2), "%\n")

###Get the PCA scores and loadings
pca_scores <- as.data.frame(pca_result$S)
pca_loadings <- pca_result$L

###Add taxon names to the PCA scores
pca_scores$taxon <- rownames(pca_scores)

###Print the PCA scores
print(pca_scores)

###Print the PCA loadings
print(pca_loadings)

###Create a data frame for the variable loadings (first three PCs)
loadings_df <- as.data.frame(pca_loadings[, 1:3])
loadings_df$variable <- rownames(loadings_df)
print(loadings_df)
#("bio1", #Annual_Mean_Temp
          # "bio2", # Mean_Diurnal_Range
          # "bio4", # Temp_Seasonality
         #  "bio8", # Mean_Temp_Wettest_Quarter
         #  "bio9", # Mean_Temp_Driest_Quarter
        #   "bio12", # Annual_Precip
        #   "bio14", # Precip_Driest_Month
        #   "bio15", # Precip_Seasonality
        #   "bio18", # Precip_Warmest_Quarter
        #   "bio21"),#Elevation

###Create a data frame for the scores
scores_df <- as.data.frame(pca_scores)

###Merge taxon, staminode type, and other relevant info
scores_df <- cbind(scores_df, data[, c("Taxon", "Floral_structure", "Habitat_type", "mas_size", "petal_size", "Average_NPP", "Life_history", "Average_NPP")])

###PERMANOVA for BIOCLIMATIC VARIABLES AND STAMINODE TYPE#####
library(RVAideMemoire)
pca_df_filtered <- scores_df[,1:3]
pca_perm <-adonis2(pca_df_filtered ~ Floral_structure,
                   data=scores_df, permutations = 9999, method = "euclidean")
pca_perm
#pairwise.perm.manova(dist(pca_df_filtered,"euclidean"),scores_df$Perennial,nperm=9999, p.method = "none")




#######################################################
####PGLM FOR FLOWER SHOWINESS ON PC1 PC2 AND PC3######
#####################################################
#Combine the data frames of scores df and showiness_res
#check row names of scores_df
head(rownames(scores_df))
#Add taxon column using row names
scores_df$Taxon <- rownames(scores_df)
#check taxon and showiness_res are in npp_data
head(npp_data[, c("Taxon", "showiness_res", "log_petal_size", "log_mas_size")])
#merge scores_df with npp_data to include showiness_res
cbind_data <-merge(scores_df, npp_data[, c("Taxon", "showiness_res","log_mas_size", "log_petal_size")], by = "Taxon")
head(cbind_data)

#MAS showiness: combine staminode presence and absence
##Run PGLS
comp_data_cbind <- comparative.data(tree, cbind_data, names.col = "Taxon")
#Staminodes
pgls_showiness_cbind <- pgls(showiness_res ~ Floral_structure, data = comp_data_cbind)
summary(pgls_showiness_cbind)


######run PGLS only for taxa with staminode###########
#Extract taxa with staminodes
taxa_with_staminodes2 <- cbind_data$Taxon[cbind_data$Floral_structure == "Staminode"]
pruned_tree2 <- drop.tip(tree, setdiff(tree$tip.label, taxa_with_staminodes2))
pruned_data2 <- cbind_data[cbind_data$Taxon %in% taxa_with_staminodes2, ]
rownames(pruned_data2) <- pruned_data2$Taxon
comp_data2 <- comparative.data(pruned_tree2, pruned_data2, names.col = "Taxon")
pgls_model2 <- pgls(showiness_res ~ PC1 + PC2 + PC3, data = comp_data2)
summary(pgls_model2)
par(mfrow = c(2, 2))
plot(pgls_model2)


######run PGLS only for taxa without staminode###########
#Extract taxa without staminodes
taxa_without_staminodes2 <- cbind_data$Taxon[cbind_data$Floral_structure == "Stamen"]
pruned_tree3 <- drop.tip(tree, setdiff(tree$tip.label, taxa_without_staminodes2))
pruned_data3 <- cbind_data[cbind_data$Taxon %in% taxa_without_staminodes2, ]
rownames(pruned_data3) <- pruned_data3$Taxon
comp_data3 <- comparative.data(pruned_tree3, pruned_data3, names.col = "Taxon")
pgls_model3 <- pgls(showiness_res ~ PC1 + PC2 + PC3, data = comp_data3)
summary(pgls_model3)
par(mfrow = c(2, 2))
plot(pgls_model3)








#Petal showiness
#Run PGLMM model including observation-level random effect
pglmm_model_petal <- pglmm(
  log_petal_size ~ PC1 + PC2 + PC3 + Floral_structure + (1 | log_mas_size), #Added observation-level random effect
  data = cbind_data,
  family = "gaussian",  #Use Gaussian if data is continuous and normally distributed
  cov_ranef = list(Taxon = tree)  #Include phylogenetic tree as a random effect
)
summary(pglmm_model_petal)


#########################################FIGURE##################################
############ggplots: BIOCLIMATIC VARIABLES:PLOTTING STAMINODE TYPE###############
################################################
###Define the environmental variable names
#bio_variable_names <- c(
#  "bio1" = "Annual Mean Temperature",
# "bio2" = "Mean Diurnal Range",
#  "bio4" = "Temperature Seasonality",
#  "bio8" = "Mean Temperature of Wettest Quarter",
#  "bio9" = "Mean Temperature of Driest Quarter",
#  "bio12" = "Annual Precipitation",
#  "bio14" = "Precipitation of Driest Month",
##  "bio15" = "Precipitation Seasonality",
#  "bio18" = "Precipitation of Warmest Quarter",
# "bio21" = "Elevation")
###Apply the named vector to the row names of the loadings data frame
#loadings_df$variable <- bio_variable_names[rownames(loadings_df)]

#Scaling factor for the loadings, adjust this as needed
scaling_factor <- 7  #Adjust this scaling factor as appropriate for the plot
scaling_factor <- max(abs(c(scores_df$PC1, scores_df$PC2))) / max(abs(c(loadings_df$PC1, loadings_df$PC2))) * 0.8


#Create a PCA1 and PC2 biplot with ggplot2
staminode_shapes <- c("Staminode" = 17, "Stamen" = 16)  # 16 for circle, 17 for triangle
library(ggplot2)
labs(x = "PC1 (36.71%)", y = "PC2 (27.11%)")
#Adjusted biplot with dashed ellipse for staminode presence
biplot_PC12 <- ggplot(data = scores_df, aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "darkgray", size = 1) +
  geom_vline(xintercept = 0, linetype = "solid", color = "darkgray", size = 1) +
  geom_segment(data = loadings_df, aes(x = 0, y = 0, xend = PC1 * scaling_factor, yend = PC2 * scaling_factor),
               arrow = arrow(length = unit(0.08, "inches")), color = "darkgray", size = 2) +
  geom_point(aes(color = Floral_structure, shape = Floral_structure), size = 8) + 
   # Add ellipses with specific line types for staminode presence and absence
  stat_ellipse(aes(color = Floral_structure, linetype = Floral_structure), 
               type = "norm", level = 0.95, linewidth = 1.5) +
  scale_linetype_manual(name = "Floral organ", 
                        values = c("Stamen" = "solid", "Staminode" = "dashed"),
                        labels = c("Stamen", "Staminode")) +
  #geom_text(data = loadings_df, aes(x = PC1 * scaling_factor, y = PC2 * scaling_factor, label = variable), 
  #          vjust = 0.6, hjust = 0.6, color = "black", size = 9, family = "Times New Roman") +
  scale_color_manual(name = "Floral organ", 
                     values = c("Staminode" = "#377EB8", "Stamen" = "darkorange"),
                     labels = c("Stamen", "Staminode")) +
  scale_shape_manual(name = "Floral organ", values = staminode_shapes, 
                     labels = c("Stamen", "Staminode")) +  
  xlim(-16, 21.5) + ylim(-15, 13) +
  theme_minimal() +
  theme(
    text = element_text(family = "Times New Roman", size = 22, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.5),  # Add x-axis line
    axis.line.y = element_line(color = "black", size = 0.5),  # Add y-axis line
    axis.ticks = element_line(color = "black", size = 1),  # Add axis ticks in black
    axis.ticks.length = unit(0.25, "cm"),  # Length of the ticks
    axis.line = element_line(color = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    legend.position = "none",
    axis.title = element_text(family = "Times New Roman", size = 22, color = "black"),
    legend.title = element_text(family = "Times New Roman", size = 24, color = "black"),#Axis title style
    axis.text = element_text(family = "Times New Roman", size = 22, color = "black"),
    legend.text = element_text(family = "Times New Roman", size = 22, color = "black")
  ) + labs(x = "PC1 (36.71%)", y = "PC2 (27.11%)") +
  annotate("text", x = min(scores_df$PC1) - 8.5, y = max(scores_df$PC2) + 5.6,
           label = "(b)", hjust = 0, vjust = 1, size = 8, family = "Times New Roman")

print(biplot_PC12)
#Save the plot as an SVG file
ggsave(filename = "PCA1_2_Staminode_type_without_taxa_point.svg", plot = biplot_PC12, device = "svg", width = 9, height = 7, units = "in")



########WITH LEGEND#########
# Define the shapes for Stamen and Staminode
#staminode_shapes <- c("Stamen" = 16, "Staminode" = 17)  # 16 for circle, 17 for triangle

# Modify the data labels in the dataset for consistency
#scores_df$Staminode_type <- factor(scores_df$Staminode_type, 
                                   #                     levels = c("Absence", "Presence"), 
                                   #                   labels = c("Stamen", "Staminode"))

# Create the PCA biplot
#biplot_PC12 <- ggplot(data = scores_df, aes(x = PC1, y = PC2)) +
  # geom_hline(yintercept = 0, linetype = "solid", color = "darkgray", size = 1) +
  # geom_vline(xintercept = 0, linetype = "solid", color = "darkgray", size = 1) +
  # geom_segment(data = loadings_df, aes(x = 0, y = 0, xend = PC1 * scaling_factor, yend = PC2 * scaling_factor),
               #              arrow = arrow(length = unit(0.08, "inches")), color = "darkgray", size = 2) +
  # geom_point(aes(color = Staminode_type, shape = Staminode_type), size = 8) + 
  # # Add ellipses with specific line types for Stamen and Staminode
  # stat_ellipse(aes(color = Staminode_type, linetype = Staminode_type), 
               #              type = "norm", level = 0.95, linewidth = 1.5) +
  # scale_linetype_manual(name = "Floral trait", 
                        #                       values = c("Stamen" = "solid", "Staminode" = "dashed"),
                        #                       labels = c("Stamen", "Staminode")) +
  #  scale_color_manual(name = "Floral trait", 
                     #                    values = c("Staminode" = "#377EB8", "Stamen" = "darkorange"),
                     #                   labels = c("Stamen", "Staminode")) +
  #  scale_shape_manual(name = "Floral trait", values = staminode_shapes, 
  #                     labels = c("Stamen", "Staminode")) +  
  # xlim(-16, 21.5) + ylim(-15, 13) +
  # theme_minimal() +
  #  theme(
    #   text = element_text(family = "Times New Roman", size = 22, color = "black"),
    #   panel.grid.major = element_blank(),
    #   panel.grid.minor = element_blank(),
    #   axis.line.x = element_line(color = "black", size = 0.5),  # Add x-axis line
    #   axis.line.y = element_line(color = "black", size = 0.5),  # Add y-axis line
    #    axis.ticks = element_line(color = "black", size = 1),  # Add axis ticks in black
    #    axis.ticks.length = unit(0.25, "cm"),  # Length of the ticks
    #   axis.line = element_line(color = "black"),
    #   panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    #   legend.position = "right",
    #   axis.title = element_text(family = "Times New Roman", size = 22, color = "black"),
    #   axis.text = element_text(family = "Times New Roman", size = 22, color = "black"),
    #   legend.text = element_text(family = "Times New Roman", size = 22, color = "black"),
    #   legend.title = element_text(family = "Times New Roman", size = 22, color = "black")  # Match legend title style to text
    #  ) + 
  #  labs(x = "PC1 (36.71%)", y = "PC2 (27.11%)") +
  #  annotate("text", x = min(scores_df$PC1) - 8.5, y = max(scores_df$PC2) + 5.6,
           #          label = "(b)", hjust = 0, vjust = 1, size = 8, family = "Times New Roman")

# Print the plot
#print(biplot_PC12)
#ggsave(filename = "this is just for legend for pc scores.svg", plot = biplot_PC12, device = "svg", width = 9, height = 7, units = "in")

################GGPLOT2 pPCA scores for MAS showiness###############
####################################################################
#################EACH SEPARATED PC##################################
####################################################################

#Load necessary library
library(ggplot2)
theme1 <-  theme(
  panel.grid.major = element_blank(),  #Remove major grid lines
  panel.grid.minor = element_blank(),  #Remove minor grid lines
  axis.line.x = element_line(color = "black", size = 0.5),  #Add x-axis line
  axis.line.y = element_line(color = "black", size = 0.5),  #Add y-axis line
  axis.ticks = element_line(color = "black", size = 1),  #Add axis ticks in black
  axis.ticks.length = unit(0.25, "cm"),  # Length of the ticks
  panel.border = element_rect(color = "black", fill = NA, size = 2),  #Add border around the plot
  axis.title = element_text(family = "Times New Roman", size = 22, color = "black"),
  legend.title = element_text(family = "Times New Roman", size = 24, color = "black"),#Axis title style
  axis.text = element_text(family = "Times New Roman", size = 22, color = "black"),  #Axis text style
  legend.text = element_text(family = "Times New Roman", size = 22, color = "black"),
  legend.position = "none"  #Hide the legend
)
#Define the shapes for Staminode presence and absence
staminode_shapes <- c("Stamen" = 16, "Staminode" = 17)  #16 for circle, 17 for triangle
#Create a scatter plot with a regression line, differentiated by Floral_structure with specified colors and shapes
pc1_showiness_plot <- ggplot(cbind_data, aes(x = PC1, y = showiness_res, color = Floral_structure, shape = Floral_structure)) +
  geom_smooth(method = "lm", aes(group = Floral_structure), size = 5, width = 0.2) +  #Add a linear regression line with confidence interval
  geom_point(size = 8) +  #Set the point size
  scale_color_manual(name = "Floral organ", values = c("Staminode" = "#377EB8", "Stamen" = "darkorange")) +  #Custom colors and legend title
  scale_shape_manual(name = "Floral organ", values = staminode_shapes) +  # Set custom shapes
  labs(x = "PC1", y = "Showiness (residuals)") +  # Label the axes
  ylim(-1.6, 1.6) +  # Set the y-axis limits
  theme_minimal() +  #Apply a minimal theme
  theme1 + labs(x = "PC1 (36.71%)") +
  annotate("text", x = min(cbind_data$PC1), y = max(cbind_data$showiness_res) + 0.28,
           label = "(b)", hjust = 0, vjust = 1, size = 8, family = "Times New Roman")
#Print the plot
print(pc1_showiness_plot)
#Save the plot as an SVG file
ggsave(filename = "showinesspc1.svg", plot = biplot_PC12, device = "svg", width = 9, height = 7, units = "in")


#Define the shapes for Staminode presence and absence
staminode_shapes <- c("Stamen" = 16, "Staminode" = 17)  #16 for circle, 17 for triangle

#Create a scatter plot with a regression line, differentiated by Floral_structure with specified colors and shapes
pc2_showiness_plot <- ggplot(cbind_data, aes(x = PC2, y = showiness_res, color = Floral_structure, shape = Floral_structure)) +
  geom_smooth(method = "lm", aes(group = Floral_structure), size = 5, width = 0.2) +  #Add a linear regression line with confidence interval
  geom_point(size = 8) +  #Set the point size
  scale_color_manual(name = "Floral organ", values = c("Staminode" = "#377EB8", "Stamen" = "darkorange")) +  #Custom colors and legend title
  scale_shape_manual(name = "Floral organ", values = staminode_shapes) +  #Set custom shapes
  labs(x = "PC2", y = "Showiness (residuals)") +  #Label the axes
  theme_minimal() +  #Apply a minimal theme
  theme1 + labs(x = "PC2 (27.11%)") +
  annotate("text", x = min(cbind_data$PC2), y = max(cbind_data$showiness_res) + 0.05,
           label = "(c)", hjust = 0, vjust = 1, size = 8, family = "Times New Roman")
#Print the plot
print(pc2_showiness_plot)
#Save the plot as an SVG file
ggsave(filename = "showinesspc2.svg", plot = biplot_PC12, device = "svg", width = 9, height = 7, units = "in")


#Define the shapes for Staminode presence and absence
staminode_shapes <- c("Stamen" = 16, "Staminode" = 17)  #16 for circle, 17 for triangle
#Create a scatter plot with a regression line, differentiated by Floral_structure with specified colors and shapes
pc3_showiness_plot <- ggplot(cbind_data, aes(x = PC3, y = showiness_res, color = Floral_structure, shape = Floral_structure)) +
  geom_smooth(method = "lm", aes(group = Floral_structure), size = 5, width = 0.2) +  #Add a linear regression line with confidence interval
  geom_point(size = 8) +  #Set the point size
  scale_color_manual(name = "Floral organ", values = c("Staminode" = "#377EB8", "Stamen" = "darkorange")) +  #Custom colors and legend title
  scale_shape_manual(name = "Floral organ", values = staminode_shapes) +  #Set custom shapes
  labs(x = "PC3", y = "Showiness (residuals)") +  #Label the axes
  theme_minimal() +  #Apply a minimal theme
  theme1 + labs(x = "PC3 (14.79%)") +
  annotate("text", x = min(cbind_data$PC3), y = max(cbind_data$showiness_res) + 0.05,
           label = "(d)", hjust = 0, vjust = 1, size = 8, family = "Times New Roman")
#Print the plot
print(pc3_showiness_plot)
#Save the plot as an SVG file
ggsave(filename = "showinesspc3.svg", plot = biplot_PC12, device = "svg", width = 9, height = 7, units = "in")



#Load necessary libraries
library(ggplot2)
library(patchwork)
#Create the PCA biplot (using the previously created biplot_PC12)
#Combine the plots using patchwork
#with legend#####
#combined_plot_npp <- (npp_plot | biplot_PC12)
#Print the combined plot
#print(combined_plot_npp)
#Save the combined plot as an SVG file
#ggsave(filename = "combined_plot_npp_with legend.svg", plot = combined_plot_npp, device = "svg", width = 14, height = 6, units = "in")


#without legend#####
combined_plot_npp <- (npp_plot | biplot_PC12)
combined_plot_npp1 <- npp_plot / biplot_PC12  # Use "/" to stack plots vertically
#Print the combined plot
print(combined_plot_npp1)
#Save the combined plot as an SVG file
ggsave(filename = "combined_plot_npp.svg", plot = combined_plot_npp, device = "svg", width = 17, height = 8, units = "in")
ggsave(filename = "combined_plot_npp1.svg", plot = combined_plot_npp1, device = "svg", width = 10, height = 15, units = "in")




#Load necessary libraries
library(ggplot2)
library(patchwork)
#Create the PCA biplot (using the previously created biplot_PC12)
#Combine the plots using patchwork
combined_plot <- (npp_showiness_plot | pc1_showiness_plot) / (pc2_showiness_plot | pc3_showiness_plot)

#Print the combined plot
print(combined_plot)

#Save the combined plot as an SVG file
ggsave(filename = "PCA_combined_plot.svg", plot = combined_plot, device = "svg", width = 14, height = 12, units = "in")





###################GGPLOT2 pPCA scores for Petal size###############
####################################################################
#################EACH SEPARATED PC##################################
####################################################################

#Define the shapes for Staminode presence and absence
staminode_shapes <- c("Stamen" = 16, "Staminode" = 17)  #16 for circle, 17 for triangle
#Create a scatter plot with a regression line, differentiated by Floral_structure with specified colors and shapes
pc1_petal <- ggplot(cbind_data, aes(x = PC1, y = log_petal_size, color = Floral_structure, shape = Floral_structure)) +
  geom_smooth(method = "lm", aes(group = Floral_structure), size = 5, width = 0.2) +  #Add a linear regression line with confidence interval
  geom_point(size = 8) +  #Set the point size
  scale_color_manual(name = "Floral organ", values = c("Staminode" = "#377EB8", "Stamen" = "darkorange")) +  #Custom colors and legend title
  scale_shape_manual(name = "Floral organ", values = staminode_shapes) +  # Set custom shapes
  labs(x = "PC1", y = "Log petal showiness (mm²)") +  # Label the axes
  theme_minimal() +  #Apply a minimal theme
  theme1 + labs(x = "PC1 (36.71%)") +
  annotate("text", x = min(cbind_data$PC1), y = max(cbind_data$showiness_res) + 6,
           label = "(b)", hjust = 0, vjust = 1, size = 8, family = "Times New Roman")
#Print the plot
print(pc1_petal)

#Create a scatter plot with a regression line, differentiated by Floral_structure with specified colors and shapes
pc2_petal <- ggplot(cbind_data, aes(x = PC2, y = log_petal_size, color = Floral_structure, shape = Floral_structure)) +
  geom_smooth(method = "lm", aes(group = Floral_structure), size = 5, width = 0.2) +  #Add a linear regression line with confidence interval
  geom_point(size = 8) +  #Set the point size
  scale_color_manual(name = "Floral organ", values = c("Staminode" = "#377EB8", "Stamen" = "darkorange")) +  #Custom colors and legend title
  scale_shape_manual(name = "Floral organ", values = staminode_shapes) +  #Set custom shapes
  labs(x = "PC2", y = "Log petal showiness (mm²)") + #Label the axes
  theme_minimal() +  #Apply a minimal theme
  theme1 + labs(x = "PC2 (27.11%)") +
  annotate("text", x = min(cbind_data$PC2), y = max(cbind_data$showiness_res) + 6,
           label = "(c)", hjust = 0, vjust = 1, size = 8, family = "Times New Roman")
#Print the plot
print(pc2_petal)

#Create a scatter plot with a regression line, differentiated by Floral_structure with specified colors and shapes
pc3_petal <- ggplot(cbind_data, aes(x = PC3, y = log_petal_size, color = Floral_structure, shape = Floral_structure)) +
  geom_smooth(method = "lm", aes(group = Floral_structure), size = 5, width = 0.2) +  #Add a linear regression line with confidence interval
  geom_point(size = 8) +  #Set the point size
  scale_color_manual(name = "Floral organ", values = c("Staminode" = "#377EB8", "Stamen" = "darkorange")) +  #Custom colors and legend title
  scale_shape_manual(name = "Floral organ", values = staminode_shapes) +  #Set custom shapes
  labs(x = "PC3", y = "Log petal showiness (mm²)") + #Label the axes
  theme_minimal() +  #Apply a minimal theme
  theme1 + labs(x = "PC3 (14.79%)") +
  annotate("text", x = min(cbind_data$PC3), y = max(cbind_data$showiness_res) + 6,
           label = "(d)", hjust = 0, vjust = 1, size = 8, family = "Times New Roman")
#Print the plot
print(pc3_petal)

#Load necessary libraries
library(ggplot2)
library(patchwork)
#Create the PCA biplot (using the previously created biplot_PC12)
#Combine the plots using patchwork
combined_petal <- (npp_petal | pc1_petal) / (pc2_petal | pc3_petal)

#Print the combined plot
print(combined_petal)

#Save the combined plot as an SVG file
ggsave(filename = "combined_petal.svg", plot = combined_petal, device = "svg", width = 14, height = 12, units = "in")

#######multiple correction across the model###############
##########################################################
##Bonferroni Correction
#Number of tests I have
#Average_NPP ~ Staminode_type (PGLS)
#showiness_res ~ Average_NPP + Staminode_type (PGLS)
#pca_df_filtered ~ Staminode_type (PERMANOVA)
#showiness_res ~ PC1 + PC2 + PC3 + Staminode_type (PGLS)
## List of p-values fro each test
p_values <- c(0.7976, 0.305, 0.010, 5e-04, 0.293, 0.827, 0.475, 0.025)
#p_values <- c(0.7976, 0.01023, 5e-04, 0.02505)
#Number of tests
n_tests <- length(p_values)
#Adjusted alpha level using Bonferroni correction
adjusted_alpha <- 0.05 / n_tests
#Apply correction and check which tests are significant
significant <- p_values < adjusted_alpha
#Results
adjusted_alpha
significant


#Holm-Bonferroni method
adjusted_p_values <- p.adjust(p_values, method = "holm")
#Print the adjusted p-values
print(adjusted_p_values)
#Optional: Filter significant results based on a chosen significance level (e.g., 0.05)
significant_results <- adjusted_p_values < 0.05
print(significant_results)
#Print the significant p-values
print(p_values[significant_results])


#Benjamini-Hochberg (FDR) Correction#) I will use this as i have multiple tests plus exploratory studies like phylogeney and permanova
adjusted_p_values <- p.adjust(p_values, method = "BH")
#Print the adjusted p-values
print(adjusted_p_values)

######USE for my Study###########
##############NPP###################
##List of p-values fro each test staminode or stamen or petal and their showiness
p_values <- c(0.153, 0.190, 0.580, 0.000, 0.648)#pANOVA for staminode, staminode showiness, stamen showiness, petal showiness, petal of staminode and stamen
adjusted_p_values <- p.adjust(p_values, method = "BH")
print(adjusted_p_values)

###########Bioclimatic variables###########
p_values1 <- c(0.000, 0.547, 0.697, 0.040, 0.207, 0.218, 0.107, 0.007, 0.255, 0.023, 0.334) #permanova, PC1 to PC3 for staminode showiness, PC1 to PC3 for stamen, and PC1 to PC3 for petal showiness, staminode presence/absence
adjusted_p_values1 <- p.adjust(p_values1, method = "BH")
print(adjusted_p_values1)



#Optional: Filter significant results based on a chosen FDR level (e.g., 0.05)
significant_results <- adjusted_p_values < 0.05
print(significant_results)

#Print the significant p-values
print(p_values[significant_results])


######correct across predictors within the model#####
#####################################################
#PGLS Model 1: NPP ~ Staminode Type
p_npp_staminode <- 0.798
#PGLS Model 2: Showiness ~ NPP + Staminode Type
p_showiness_npp <- 0.310
p_showiness_staminode <- 0.010
#PERMANOVA Model: Environmental Differences (PC1, PC2, PC3) ~ Staminode Type
p_perm_staminode <- 0.0005
#PGLS Model 3: Showiness ~ PC1 + PC2 + PC3 + Staminode Type
p_showiness_pc1 <- 0.017
p_showiness_pc2 <- 0.289
p_showiness_pc3 <- 0.710
p_showiness_staminode_pc <- 0.056

#Combine all p-values from individual predictors across models
all_p_values <- c(p_npp_staminode,
                  p_showiness_npp,
                  p_showiness_staminode,
                  p_perm_staminode,
                  p_showiness_pc1,
                  p_showiness_pc2,
                  p_showiness_pc3,
                  p_showiness_staminode_pc)

#Apply Benjamini-Hochberg correction
adjusted_p_values <- p.adjust(all_p_values, method = "BH")

#Print the adjusted p-values
print(adjusted_p_values)


######correct across predictors accross the modelaltogether#####
#####################################################
p_values <- c(
  0.000, #INTERCEPT NPP effect from PGLS (NPP and Staminode Presence/Absence)
  0.798, #Staminode effect on NPP from PGLS
  0.038, #INTERCETP NPP effect from PGLS (MAS Showiness)
  0.305, #NPP effect from PGLS (MAS Showiness)
  0.010, #Staminode effect from PGLS (MAS Showiness)
  0.090, #INTERCEPT effect from PGLS (PC1, PC2, PC3 on MAS Showiness)
  0.293, #PC1 effect on mas from PGLS
  0.827, #PC2 effect on mas from PGLS
  0.475, #PC3 effect on mas from PGLS
  0.025, #Staminode effect on MAS from PGLS
  0.000, #Intercept effect from PGLMM (Petal Size with NPP and Staminode)
  0.000, #NPP effect on petal size from PGLMM
  0.648, #Staminode effect on petal size from PGLMM
  0.000, #Intercept effect from PGLMM (Petal Size with PC1, PC2, PC3, and Staminode)
  0.007, #PC1 effect from PGLMM
  0.255, #PC2 effect from PGLMM
  0.023, #PC3 effect from PGLMM
  0.334, #Staminode effect from PGLMM
  2e-04  #P-value from PERMANOVA (Staminode_type effect on bioclimatic variables)
)

#Apply the Benjamini-Hochberg correction
adjusted_p_values <- p.adjust(p_values, method = "BH")
#Display the adjusted P-values
adjusted_p_values


p_values1 <- c(
  0.000, #INTERCEPT NPP effect from PGLS (NPP and Staminode Presence/Absence)
  0.798, #Staminode effect on NPP from PGLS
  0.038, #INTERCETP NPP effect from PGLS (MAS Showiness)
  0.305, #NPP effect from PGLS (MAS Showiness)
  0.010, #Staminode effect from PGLS (MAS Showiness)
  0.090, #INTERCEPT effect from PGLS (PC1, PC2, PC3 on MAS Showiness)
  0.293, #PC1 effect on mas from PGLS
  0.827, #PC2 effect on mas from PGLS
  0.475, #PC3 effect on mas from PGLS
  0.025, #Staminode effect on MAS from PGLS
  0.000, #Intercept effect from PGLMM (Petal Size with NPP and Staminode)
  0.000, #NPP effect on petal size from PGLMM
  0.648, #Staminode effect on petal size from PGLMM
  0.000, #Intercept effect from PGLMM (Petal Size with PC1, PC2, PC3, and Staminode)
  0.007, #PC1 effect from PGLMM
  0.255, #PC2 effect from PGLMM
  0.023, #PC3 effect from PGLMM
  0.334, #Staminode effect from PGLMM
  2e-04  #P-value from PERMANOVA (Staminode_type effect on bioclimatic variables)
)
#Apply the Benjamini-Hochberg correction
adjusted_p_values1 <- p.adjust(p_values1, method = "BH")
#Display the adjusted P-values
adjusted_p_values1


#########################################################################
####Calculate the mean and standard error for the entire dataset#########
#########################################################################

mean_showiness <- mean(scores_df$Showiness)
se_showiness <- sd(scores_df$Showiness) / sqrt(nrow(scores_df))

#Print the results for the entire dataset
cat("Mean Showiness (Total):", mean_showiness, "\n")
cat("Standard Error (Total):", se_showiness, "\n")

#####Calculate the mean and standard error for each Staminode_type########
showiness_summary <- scores_df %>%
  group_by(Staminode_type) %>%
  summarise(
    Mean_Showiness = mean(Showiness),
    SE_Showiness = sd(Showiness) / sqrt(n())
  )
#Print the results for each Staminode_type
print(showiness_summary)

#####Calculate the sample size for the entire dataset#######
total_sample_size <- nrow(scores_df)
#Print the sample size for the entire dataset
cat("Sample Size (Total):", total_sample_size, "\n")

######Calculate the sample size for each Staminode_type#####
sample_size_by_staminode <- scores_df %>%
  group_by(Staminode_type) %>%
  summarise(Sample_Size = n())

######Print the sample size for each Staminode_type#######
print(sample_size_by_staminode)












##########BELOW SCRIPTS ARE NOT FOR MANUSCRIPT BUT GOOD FOR PRACTICES######
################################NOT FOR MANUSCRIPT#########################
################################SHOWINESS USING pPCA#######################
###########################################################################


#########################################################################
###########################PGLMM FOR SHOWINES VS PCS AND STAMINODE#######
#########################################################################

#Fit a PGLMM with Habitat_type and Life_history as random effects and Gaussian family
pglmm_model <- pglmm(
  Log_Showiness ~ PC1 + PC2 + PC3 + Staminode_type + (1 | Habitat_type) + (1 | Life_history) ,
  data = model_data,
  family = "gaussian",  #Use Gaussian for normally distributed continuous response variable
  cov_ranef = list(Taxon = tree)  #Include phylogenetic tree as a random effect
)
#Print summary of the model
summary(pglmm_model)


#########################################################################
###########################PGLMM FOR PETAL VS PCS AND STAMINODE#######
#########################################################################
#Fit a PGLMM with Habitat_type and Life_history as random effects and Gaussian family
pglmm_model_petal <- pglmm(
  log_petal_size ~ PC1 + PC2 + PC3 + Staminode_type  + (1 | Habitat_type) + (1 | Life_history) ,
  data = model_data,
  family = "gaussian",  #Use Gaussian for normally distributed continuous response variable
  cov_ranef = list(Taxon = tree)  #Include phylogenetic tree as a random effect
)
#Print summary of the model
summary(pglmm_model_petal)

#########################################################################
###########################PGLMM FOR MAS VS PCS AND STAMINODE#######
#########################################################################
#Fit a PGLMM with Habitat_type and Life_history as random effects and Gaussian family
pglmm_model_mas <- pglmm(
  log_mas_size ~ PC1 + PC2 + PC3 + Staminode_type  + (1 | Habitat_type) + (1 | Life_history),
  data = model_data,
  family = "gaussian",  #Use Gaussian for normally distributed continuous response variable
  cov_ranef = list(Taxon = tree)  #Include phylogenetic tree as a random effect
)
#Print summary of the model
summary(pglmm_model_mas)


#########################################################################
###########################PGLM FOR SHOWINES VS PCS AND STAMINODE########
#########################################################################
pglm_model <- phylolm(
  Log_Showiness ~ PC1 + PC2 + PC3 + Staminode_type,  #Model formula
  data = model_data,                                 #Data frame containing the variables
  phy = tree,                                         #Phylogenetic tree
  model = "BM",                                       #Brownian Motion model for phylogenetic signal
  family = gaussian()                                 #Use Gaussian for log-transformed continuous data
)
#Summary of the model
summary(pglm_model)


#########################################################################
###########################PGLM FOR LIFE HISTORY AND HABITAT TYPE########
#########################################################################

#####Convert Life_history column to binary numeric values###

model_data <- model_data %>%
  mutate(Perennial = ifelse(Life_history == "Perennial", 1, 0))  #"Perennial" to 1, "Non_perennial" to 0
#View the updated data frame
head(model_data)

#Convert Habitat_type column to binary numeric values
model_data <- model_data %>%
  mutate(Habitat = ifelse(Habitat_type == "Gypsum", 1, 0))  #"Gypsum" to 1, "Non_gypsum" to 0
#View the updated data frame
head(model_data)


#######################################################################
################################Gypsum affinity########################
#Fit a PGLMM with Habitat_type as random effect and Gaussian family
pglmm_lifehistory <- pglmm(
  Habitat ~ Staminode_type + (1 | Life_history) ,
  data = model_data,
  family = "binomial",  #Use Gaussian for normally distributed continuous response variable
  cov_ranef = list(Taxon = tree)  #Include phylogenetic tree as a random effect
)
#Print summary of the model
summary(pglmm_lifehistory)


###################################RUN PGLMM#########################
##################################Life history#######################
#binary data (presence and absence)
#Fit a PGLMM with Habitat_type and Life_history as random effects and Gaussian family
pglmm_lifeHistory <- pglmm(
  Perennial ~ Staminode_type + (1 | Habitat_type),
  data = model_data,
  family = "binomial",  #Use Gaussian for normally distributed continuous response variable
  cov_ranef = list(Taxon = tree)  #Include phylogenetic tree as a random effect
)
#Print summary of the model
summary(pglmm_lifeHistory)


################################################################
######################NPP Phylogenetic PGLS ####################NOT FOR MANUSCRIPT
################################################################

#Create a comparative data object
library(caper)
comp_data <- comparative.data(phy = tree, data = npp_data, names.col = "Taxon", vcv = TRUE, na.omit = FALSE)

#Define the PGLS model with Average_NPP as the response variable and Staminode type as the predictor
pgls_model <- pgls(Average_NPP ~ Staminode_type, data = comp_data)
#Print the summary of the PGLS model
summary(pgls_model)

#Diagnostic plots to check residuals
par(mfrow = c(2, 2))  #Adjusts the layout for multiple plots
plot(pgls_model)
#Reset plotting layout
par(mfrow = c(1, 1))



################################################################
#####################PGLMM FOR NPP##############################NOT FOR MANUSCRIPT
################################################################
#Load necessary library
library(phyr)
#Fit a PGLMM with Habitat_type and Life_history as random effects and Gaussian family
#pglmm_model <- pglmm(
 # Log_Showiness ~ Average_NPP * Staminode_type + (1 | Habitat_type) + (1 | Life_history),
 # data = npp_data,
#  family = "gaussian",  #Use Gaussian for normally distributed continuous response variable
#  cov_ranef = list(Taxon = tree)  #Include phylogenetic tree as a random effect
#)
#Print summary of the model
#summary(pglmm_model)

########################################
#############PGLMM FOR PETAL############NOT FOR MANUSCRIPT
########################################
###WITH RANDOM
#Log-transform the petal size
npp_data$log_petal_size <- log(npp_data$petal_size + 1) 
pglmm_model_petal <- pglmm(
  log_petal_size ~ Average_NPP * Staminode_type + (1 | Life_history) + (1 | Habitat_type),  
  data = npp_data,  
  family = "gaussian",  
  cov_ranef = list(Taxon = tree)
)
#Print summary of the model
summary(pglmm_model_petal)


###########PGLMM FOR MAS#############NOT FOR MANUSCRIPT
#####################################
#WITH RANDOM
npp_data$log_mas_size <- log(npp_data$mas_size + 1) 
pglmm_model_MAS <- pglmm(
  log_mas_size ~ Average_NPP * Staminode_type + (1 | Life_history) + (1 | Habitat_type), 
  data = npp_data,  
  family = "gaussian", 
  cov_ranef = list(Taxon = tree) 
)
# Print summary of the model
summary(pglmm_model_MAS)

########################################################
####Check phylogenetic signal for showiness#############
########################################################
#Install and load necessary libraries
library(phytools)
library(geiger)

#Prepare the data
model_data <- scores_df
rownames(model_data) <- model_data$Taxon
#Extract the Showiness trait and the corresponding phylogeny
showiness_trait <- setNames(model_data$Showiness, rownames(model_data))
phylo <- tree

# Log-transform the Showiness trait
showiness_trait <- log(setNames(model_data$Showiness, rownames(model_data)))

#Calculate Pagel's Lambda
lambda_result <- phylosig(phylo, showiness_trait, method = "lambda", test = TRUE)
print(lambda_result)
#Calculate Blomberg's K
k_result <- phylosig(phylo, showiness_trait, method = "K", test = TRUE)
print(k_result)


#Load necessary libraries
library(caper)
library(ape)
library(car)
#Find the minimum value in the Showiness column
min_showiness <- min(scores_df$Showiness, na.rm = TRUE)

#Add a small constant to shift the values if the minimum is less than or equal to zero
#Choose a constant that will ensure all values are positive
constant <- if (min_showiness <= 0) { abs(min_showiness) + 1 } else { 0 }

#Apply the log transformation after adding the constant
scores_df$log_Showiness <- log(scores_df$Showiness + constant)

#Check normality again after transformation
ggplot(scores_df, aes(x = log_Showiness)) +
  geom_histogram(binwidth = 0.5, color = "black", fill = "white") +
  theme_minimal()

qqPlot(scores_df$log_Showiness, main = "Q-Q Plot for Log-Transformed Showiness")

shapiro.test(scores_df$log_Showiness)

#Ensure your tree and data are compatible
comparative_data <- comparative.data(tree, scores_df, names.col = "Taxon", vcv = TRUE, na.omit = TRUE)

pgls_pc_model <- pgls(log_Showiness ~ PC1 + PC2 + PC3, data = comparative_data)
summary(pgls_pc_model)

##################################NOT FOR MANUSCRIPT############################
#####################SHOWINESS USING PGLS AND WITHOUT USING pPCA################
################################################################################
#Load necessary libraries
library(phytools)
library(ape)
library(caper)
library(dplyr)

#Load the CSV file
data <- read.csv("geo_cbind_taxon_env_plus_trait.csv")

#Perform normality tests and apply transformations if needed
env_data <- dplyr::select(data, starts_with("bio"))

#Function to decide transformation based on normality
transform_variable <- function(x) {
  p_value <- shapiro.test(x)$p.value
  if (p_value < 0.05) {
    return(log1p(x - min(x, na.rm = TRUE) + 1)) #Log transformation for non-normal distribution
  } else {
    return(x) #No transformation needed
  }
}

#Apply transformation
transformed_env_data <- env_data %>% mutate(across(everything(), transform_variable))

# Perform PCA on the transformed data (Standardized and centered PCA)
pca_result <- prcomp(transformed_env_data, scale. = TRUE, center = TRUE)

# Create a data frame for the PCA scores (first 3 components)
scores_df <- as.data.frame(pca_result$x[, 1:3])
colnames(scores_df) <- c('PC1', 'PC2', 'PC3')

#Log transform 'Showiness'#data$Showiness_log <- log(data$Showiness)

data$Showiness <- data$Showiness##we dont use log transformation as geometric mean is log-transformed

#Merge the PCA scores with the transformed showiness and Taxon information
merged_data <- cbind(scores_df, data[, c("Taxon", "Showiness")])

#Load the phylogenetic tree
tree <- read.tree("BartoniaIngroupPLnames.tre")

#Ensure your tree and data are compatible
comparative_data <- comparative.data(tree, merged_data, names.col = "Taxon", vcv = TRUE, na.omit = TRUE)

#Run PGLS
pgls_model <- pgls(Showiness ~ PC1 + PC2 + PC3, data = comparative_data)

#Summarize the results
summary(pgls_model)#when i used log transformation pc1 is signigicant



##################################################################
####Generate boxplots for each environmental variable#############
############Not for the paper#####################################
bio_variables <- c("bio1", "bio2", "bio4", "bio8", "bio9", "bio12", "bio14", "bio15", "bio18", "bio21")

###Create a list to store the plots
plots <- list()

###Generate and store the plots in the list
for (bio in bio_variables) {
  plot <- ggplot(data, aes_string(x = "Staminode_type", y = bio, fill = "Staminode_type")) +
    geom_boxplot() +
    theme_minimal() +
    labs(title = paste("Boxplot of", bio, "by Staminode Type"), x = "Staminode Type", y = bio) +
    scale_fill_manual(values = c("Absence" = "darkgreen", "Presence" = "orange")) +
    theme(
      text = element_text(family = "Times New Roman", size = 15),
      axis.title = element_text(family = "Times New Roman", size = 15),
      axis.text = element_text(family = "Times New Roman", size = 12),
      legend.position = "none"
    )
  plots[[bio]] <- plot
}
library(gridExtra)
###Combine all the plots into one
combined_plot <- grid.arrange(grobs = plots, ncol = 3)

###Print the combined plot
print(combined_plot)

###Save the combined plot to a file
#ggsave(filename = "combined_boxplots.tif", plot = combined_plot, width = 20, height = 15, units = "in")

######################################################################
####################PCA without phylogeny#############################
####PERFORMING PCA ON GEOMETRIC MEAN DATASET##########################
# Load necessary libraries
library(ggplot2)
library(tidyr)

# Load the CSV file
data <- read.csv("geo_cbind_taxon_env_plus_trait.csv")
#Perform normality tests and apply transformations if needed
env_data <- dplyr::select(data, starts_with("bio"))

#Dont apply log transformation
pca_result <- prcomp(env_data, scale. = TRUE, center = TRUE)# Standardized and centered PCA
#pca_result <- prcomp(env_data, scale. = FALSE, center = TRUE)# for performing the PCA using the covariance matrix, only center the data, do not scale
#Get summary of PCA result
pca_summary <- summary(pca_result)
# Print the proportion of variance explained by each PC
print(pca_summary$importance)
# Extract the proportion of variance for the first three PCs
proportion_variance <- pca_summary$importance[2, 1:3] * 100
# Print the percentage of variance explained by PC1, PC2, and PC3
cat("Percentage of Variance Explained by PC1:", round(proportion_variance[1], 2), "%\n")
cat("Percentage of Variance Explained by PC2:", round(proportion_variance[2], 2), "%\n")
cat("Percentage of Variance Explained by PC3:", round(proportion_variance[3], 2), "%\n")

#Create a data frame for the variable loadings
loadings_df <- as.data.frame(pca_result$rotation[, 1:3])#Select the first three principal components
# Print the loadings
print(loadings_df)

#Optionally, create a more readable data frame
loadings_df1 <- as.data.frame(loadings_df)
loadings_df1$variable <- rownames(loadings_df1)
print(loadings_df1)
loadings_df$variable <- rownames(pca_result$rotation) #Name the loadings
#Create a data frame for the scores
#Directly select all components instead of just the first three for potential full analysis later
scores_df <- as.data.frame(pca_result$x)  # All principal component scores
#Merge taxon, staminode type, and gypsum type info
scores_df <- cbind(scores_df, data[, c("Taxon", "Staminode_type", "Habitat_type", "Life_history", "petal_size", "mas_size", "Average_NPP")])
#write.csv(scores_df, file = "PC_scores_Mentzelia.csv", row.names = FALSE)

###PERMANOVA for PERENNIAL#####
library(RVAideMemoire)
pca_df_filtered <- scores_df[,1:3]
pca_perm <-adonis2(pca_df_filtered~ Staminode_type * Life_history * Habitat_type, 
                   data=scores_df, permutations = 9999, method = "euclidean")
pca_perm
pairwise.perm.manova(dist(pca_df_filtered,"euclidean"),scores_df$Perennial,nperm=9999, p.method = "none")
pairwise.perm.manova(dist(pca_df_filtered,"euclidean"),scores_df$Staminode_type,nperm=9999, p.method = "none")
pairwise.perm.manova(dist(pca_df_filtered,"euclidean"),scores_df$Gypsum_type,nperm=9999, p.method = "none")

################################################
##########PLOTTING STAMINODE TYPE###############
################################################
#Define the environmental variable names
bio_variable_names <- c(
  "bio1" = "Annual Mean Temperature",
  "bio2" = "Mean Diurnal Range",
  "bio3" = "Isothermality",
  "bio4" = "Temperature Seasonality",
  "bio5" = "Max Temperature of Warmest Month",
  "bio6" = "Min Temperature of Coldest Month",
  "bio7" = "Temperature Annual Range",
  "bio8" = "Mean Temperature of Wettest Quarter",
  "bio9" = "Mean Temperature of Driest Quarter",
  "bio10" = "Mean Temperature of Warmest Quarter",
  "bio11" = "Mean Temperature of Coldest Quarter",
  "bio12" = "Annual Precipitation",
  "bio13" = "Precipitation of Wettest Month",
  "bio14" = "Precipitation of Driest Month",
  "bio15" = "Precipitation Seasonality",
  "bio16" = "Precipitation of Wettest Quarter",
  "bio17" = "Precipitation of Driest Quarter",
  "bio18" = "Precipitation of Warmest Quarter",
  "bio19" = "Precipitation of Coldest Quarter",
  "bio20" = "Solar Radiation",
  "bio21" = "Elevation"
)

#Apply the named vector to the row names of the loadings data frame
loadings_df$variable <- bio_variable_names[rownames(loadings_df)]

#Scaling factor for the loadings, adjust this as needed
#scaling_factor <- 7  # Adjust this scaling factor as appropriate for the plot
scaling_factor <- max(abs(c(scores_df$PC1, scores_df$PC2))) / max(abs(c(loadings_df$PC1, loadings_df$PC2))) * 0.8

#Create a PCA1 and PC2 biplot with ggplot2
staminode_shapes <- c("Absence" = 16, "Presence" = 17)  # 16 for circle, 17 for triangle
library(ggplot2)


###Create the first biplot with label "A"
biplot_PC12 <- ggplot(data = scores_df, aes(x = PC1, y = PC2)) +
  geom_segment(data = loadings_df, aes(x = 0, y = 0, xend = PC1 * scaling_factor, yend = PC2 * scaling_factor),
               arrow = arrow(length = unit(0.05, "inches")), color = "darkgray", size = 1) +
  geom_point(aes(color = Staminode_type, shape = Staminode_type), size = 5) + 
  stat_ellipse(aes(color = Staminode_type), type = "norm", level = 0.95, linewidth = 1.5)+
  scale_linetype_manual(name = "Staminode", values = c("Absence" = "dashed", "Presence" = "solid"),
                        labels = c("Absence", "Presence")) +
  geom_text(data = loadings_df, aes(x = PC1 * scaling_factor, y = PC2 * scaling_factor, label = variable), 
            vjust = 0.6, hjust = 0.6, color = "black", size = 4, family = "Times New Roman") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  scale_color_manual(name = "Staminode", values = c("Presence" = "orange", "Absence" = "darkgreen"),
                     labels = c("Absence", "Presence")) +
  scale_shape_manual(name = "Staminode", values = staminode_shapes, 
                     labels = c("Absence", "Presence")) +  
  xlim(-7, 6) + ylim(-5.8, 4) +
  theme_minimal() +
  theme(
    text = element_text(family = "Times New Roman", size = 20, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    legend.position = "right",
    axis.title = element_text(family = "Times New Roman", size = 15, color = "black"),
    axis.text = element_text(family = "Times New Roman", size = 15, color = "black"),
    legend.text = element_text(family = "Times New Roman", size = 15, color = "black")
  ) + labs(x = paste0("PC1 (", round(pca_result$sdev[1]^2 / sum(pca_result$sdev^2) * 100, 2), "%)"),
       y = paste0("PC2 (", round(pca_result$sdev[2]^2 / sum(pca_result$sdev^2) * 100, 2), "%)"))
 # annotate("text", x = -5.5, y = 4, label = "A", size = 6, family = "Times New Roman")

print(biplot_PC12)

#Save the plot as an SVG file
ggsave(filename = "PCA1_2_Staminode_type_without_taxa_point.svg", plot = biplot_PC12, device = "svg", width = 9, height = 7, units = "in")

###########################################################
#####################TAXA on PC1 and PC2####################
############################################################

# Create a PCA1 and PC2 biplot with ggplot2
biplot_PC_taxa <- ggplot(data = scores_df, aes(x = PC1, y = PC2)) +
  geom_segment(data = loadings_df, aes(x = 0, y = 0, xend = PC1 * scaling_factor, yend = PC2 * scaling_factor),
               arrow = arrow(length = unit(0.05, "inches")), color = "darkgray", size = 1) +
  geom_text(aes(label = Taxon, color = Staminode_type), size = 3) + 
  stat_ellipse(aes(color = Staminode_type), type = "norm", level = 0.95, linewidth = 1.5) +
  geom_text(data = loadings_df, aes(x = PC1 * scaling_factor, y = PC2 * scaling_factor, label = variable), 
            vjust = 0.6, hjust = 0.6, color = "black", size = 4, family = "Times New Roman") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  scale_color_manual(name = "Staminode", values = c("Presence" = "red", "Absence" = "darkgreen"),
                     labels = c("Absence", "Presence")) +
  xlim(-7, 6) + ylim(-5.8, 4) +
  theme_minimal() +
  theme(
    text = element_text(family = "Times New Roman", size = 20, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    legend.position = "right",
    axis.title = element_text(family = "Times New Roman", size = 15, color = "black"),
    axis.text = element_text(family = "Times New Roman", size = 15, color = "black"),
    legend.text = element_text(family = "Times New Roman", size = 15, color = "black")
  ) +
  labs(x = paste0("PC1 (", round(pca_result$sdev[1]^2 / sum(pca_result$sdev^2) * 100, 2), "%)"),
       y = paste0("PC2 (", round(pca_result$sdev[2]^2 / sum(pca_result$sdev^2) * 100, 2), "%)"))

print(biplot_PC_taxa)

#Save the plot as an SVG file
ggsave(filename = "PCA1_2_Staminode_type_with_taxa_names.svg", plot = biplot_PC_taxa, device = "svg", width = 9, height = 7, units = "in")




####################################
###########ggplot with no############
# Required libraries
library(ggplot2)
library(grid)

# Example dataset for illustration
scores_df <- data.frame(
  PC1 = rnorm(100),
  PC2 = rnorm(100),
  Staminode_type = sample(c("Presence", "Absence"), 100, replace = TRUE)
)

loadings_df <- data.frame(
  PC1 = runif(10, -1, 1),
  PC2 = runif(10, -1, 1),
  variable = paste("Var", 1:10)
)

scaling_factor <- 1  # Adjust as needed
staminode_shapes <- c("Presence" = 16, "Absence" = 17)

# PCA result placeholder for illustration
pca_result <- list(sdev = c(1.5, 1.0, 0.5))

# Creating the biplot
biplot_PC12 <- ggplot(data = scores_df, aes(x = PC1, y = PC2)) +
  geom_segment(data = loadings_df, aes(x = 0, y = 0, xend = PC1 * scaling_factor, yend = PC2 * scaling_factor),
               arrow = arrow(length = unit(0.05, "inches")), color = "darkgray", linewidth = 1) +
  geom_point(aes(color = Staminode_type, shape = Staminode_type), size = 5) + 
  stat_ellipse(aes(color = Staminode_type, linetype = Staminode_type), type = "norm", level = 0.95, linewidth = 1) +
  scale_linetype_manual(name = "Staminode type", values = c("Absence" = "dashed", "Presence" = "solid"),
                        labels = c("Absence", "Presence")) +
  geom_text(data = loadings_df, aes(x = PC1 * scaling_factor, y = PC2 * scaling_factor, label = variable), 
            vjust = 0.6, hjust = 0.6, color = "black", size = 4, family = "Times New Roman") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  scale_color_manual(name = "Staminode type", values = c("Presence" = "red", "Absence" = "darkgreen"),
                     labels = c("Absence", "Presence")) +
  scale_shape_manual(name = "Staminode type", values = staminode_shapes, 
                     labels = c("Absence", "Presence")) +  
  xlim(-5.8, 5) + ylim(-4.8, 4) +
  theme_minimal() +
  theme(
    text = element_text(family = "Times New Roman", size = 15, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    legend.position = "right",
    axis.title = element_text(family = "Times New Roman", size = 15, color = "black"),
    axis.text = element_text(family = "Times New Roman", size = 15, color = "black"),
    legend.text = element_text(family = "Times New Roman", size = 15, color = "black")
  ) + labs(x = paste0("PC1 (", round(pca_result$sdev[1]^2 / sum(pca_result$sdev^2) * 100, 2), "%)"),
           y = paste0("PC2 (", round(pca_result$sdev[2]^2 / sum(pca_result$sdev^2) * 100, 2), "%)")) +
  annotate("text", x = -5.5, y = 4, label = "A", size = 6, family = "Times New Roman")

print(biplot_PC12)









loadings_df

##########################################################################
#################################PC1 to PC3###############################
##########################################################################
scaling_factor <- max(abs(c(loadings_df$PC1, loadings_df$PC3))) / max(abs(c(loadings_df$PC1, loadings_df$PC3))) * 0.8

#Create a PCA1 and PC2 biplot with ggplot2
staminode_shapes <- c("Absence" = 16, "Presence" = 17)  # 16 for circle, 17 for triangle
# Create the first biplot with label "A"
biplot_PC13 <- ggplot(data = scores_df, aes(x = PC1, y = PC3)) +
  geom_segment(data = loadings_df, aes(x = 0, y = 0, xend = PC1 * scaling_factor, yend = PC3 * scaling_factor),
               arrow = arrow(length = unit(0.05, "inches")), color = "darkgray", size = 1.5) +
  geom_point(aes(color = Staminode_type, shape = Staminode_type), size = 5) + 
  stat_ellipse(aes(color = Staminode_type, linetype = Staminode_type), level = 0.95, size = 1.5) +
  scale_linetype_manual(name = "Staminode type", values = c("Absence" = "dashed", "Presence" = "solid"),
                        labels = c("Absence", "Presence")) +
  geom_text(data = loadings_df, aes(x = PC1 * scaling_factor, y = PC3 * scaling_factor, label = variable), 
            vjust = 0.6, hjust = 0.6, color = "black", size = 4, family = "Times New Roman") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  scale_color_manual(name = "Staminode type", values = c("Presence" = "orange", "Absence" = "darkgreen"),
                     labels = c("Absence", "Presence")) +
  scale_shape_manual(name = "Staminode type", values = staminode_shapes, 
                     labels = c("Absence", "Presence")) +  
  xlim(-5.8, 5) + ylim(-4.8, 4) +
  theme_minimal() +
  theme(
    text = element_text(family = "Times New Roman", size = 20, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    legend.position = "right",
    axis.title = element_text(family = "Times New Roman", size = 20, color = "black"),
    axis.text = element_text(family = "Times New Roman", size = 20, color = "black"),
    legend.text = element_text(family = "Times New Roman", size = 20, color = "black")
  ) + labs(x = paste0("PC1 (", round(pca_result$sdev[1]^2 / sum(pca_result$sdev^2) * 100, 2), "%)"),
           y = paste0("PC3 (", round(pca_result$sdev[3]^2 / sum(pca_result$sdev^2) * 100, 2), "%)")) +
  annotate("text", x = -5.5, y = 4, label = "B", size = 6, family = "Times New Roman")

print(biplot_PC13)

#Save the plot as an SVG file
#ggsave(filename = "PCA1_2_Staminode_type.svg", plot = biplot, device = "svg", width = 8, height = 6, units = "in")
ggsave(filename = "PCA1_3_Staminode_type_without_A.svg", plot = biplot_PC13, device = "svg", width = 9, height = 7, units = "in")


##########################################################################
#################################PC2 to PC3###############################
##########################################################################
scaling_factor <- max(abs(c(scores_df$PC2, scores_df$PC3))) / max(abs(c(loadings_df$PC1, loadings_df$PC3))) * 0.8
staminode_shapes <- c("Absence" = 16, "Presence" = 17)  # 16 for circle, 17 for triangle
# Create the first biplot with label "A"
biplot_PC23 <- ggplot(data = scores_df, aes(x = PC2, y = PC3)) +
  geom_segment(data = loadings_df, aes(x = 0, y = 0, xend = PC2 * scaling_factor, yend = PC3 * scaling_factor),
               arrow = arrow(length = unit(0.05, "inches")), color = "darkgray", size = 1.5) +
  geom_point(aes(color = Staminode_type, shape = Staminode_type), size = 5) + 
  stat_ellipse(aes(color = Staminode_type, linetype = Staminode_type), level = 0.95, size = 1.5) +
  scale_linetype_manual(name = "Staminode type", values = c("Absence" = "dashed", "Presence" = "solid"),
                        labels = c("Absence", "Presence")) +
  geom_text(data = loadings_df, aes(x = PC2 * scaling_factor, y = PC3 * scaling_factor, label = variable), 
            vjust = 0.6, hjust = 0.6, color = "black", size = 4, family = "Times New Roman") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  scale_color_manual(name = "Staminode type", values = c("Presence" = "orange", "Absence" = "darkgreen"),
                     labels = c("Absence", "Presence")) +
  scale_shape_manual(name = "Staminode type", values = staminode_shapes, 
                     labels = c("Absence", "Presence")) +  
  xlim(-5.5, 5) + ylim(-4.8, 4) +
  theme_minimal() +
  theme(
    text = element_text(family = "Times New Roman", size = 20, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    legend.position = "right",
    axis.title = element_text(family = "Times New Roman", size = 20, color = "black"),
    axis.text = element_text(family = "Times New Roman", size = 20, color = "black"),
    legend.text = element_text(family = "Times New Roman", size = 20, color = "black")
  ) + labs(x = paste0("PC2 (", round(pca_result$sdev[2]^2 / sum(pca_result$sdev^2) * 100, 2), "%)"),
           y = paste0("PC3 (", round(pca_result$sdev[3]^2 / sum(pca_result$sdev^2) * 100, 2), "%)")) +
  annotate("text", x = -5.5, y = 4, label = "C", size = 6, family = "Times New Roman")

print(biplot_PC23)

#Save the plot as an SVG file
ggsave(filename = "PCA2_3_Staminode_type_without_A.svg", plot = biplot_PC23, device = "svg", width = 9, height = 7, units = "in")


library(cowplot)
library(ggplot2)

# Combine the two plots with equal heights
staminode_pc123 <- plot_grid(biplot_PC12, biplot_PC13, biplot_PC23, ncol = 1, align = "v", rel_heights = c(1, 1, 1))

# Print the combined plot
print(staminode_pc123)

# Save the combined plot as an SVG file
ggsave(filename = "Staminode_PC123.svg", plot = staminode_pc123, device = "svg", width = 9, height = 19, units = "in")




##########################################
##########PLOTTING LIFE FORM ###############
#################PC1 TO 2#######################
scaling_factor <- max(abs(c(scores_df$PC1, scores_df$PC2))) / max(abs(c(loadings_df$PC1, loadings_df$PC2))) * 0.8
#Create a new column without underscores in Non_perennial
scores_df$Life_history_clean <- gsub("_", " ", scores_df$Life_history)

life_history_shapes <- c("Non perennial" = 16, "Perennial" = 17)  # 16 for circle, 17 for triangle

lbiplot_PC12 <-  ggplot(data = scores_df, aes(x = PC1, y = PC2)) +
  geom_segment(data = loadings_df, aes(x = 0, y = 0, xend = PC1 * scaling_factor, yend = PC2 * scaling_factor),
               arrow = arrow(length = unit(0.05, "inches")), color = "darkgray", size = 1.5) +
  geom_point(aes(color = Life_history_clean, shape = Life_history_clean), size = 6) + 
  stat_ellipse(aes(color = Life_history_clean, linetype = Life_history_clean), level = 0.95, size = 1.5) +
  scale_linetype_manual(name = "Life history", values = c("Non perennial" = "dashed", "Perennial" = "solid"),
                        labels = c("Non perennial", "Perennial")) +
  geom_text(data = loadings_df, aes(x = PC1 * scaling_factor, y = PC2 * scaling_factor, label = variable), 
            vjust = 0.6, hjust = 0.6, color = "black", size = 3.5, family = "Times New Roman") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  scale_color_manual(name = "Life history", values = c("Perennial" = "orange", "Non perennial" = "darkgreen"),
                     labels = c("Non perennial", "Perennial")) +
  scale_shape_manual(name = "Life history", values = life_history_shapes, 
                     labels = c("Non perennial", "Perennial")) +  
  xlim(-7, 5.5) + ylim(-5, 4.5) +
  theme_minimal() +
  theme(
    text = element_text(family = "Times New Roman", size = 22, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    legend.position = "right",
    axis.title = element_text(family = "Times New Roman", size = 20, color = "black"),
    axis.text = element_text(family = "Times New Roman", size = 20, color = "black"),
    legend.text = element_text(family = "Times New Roman", size = 18, color = "black")
  ) +
labs(x = paste0("PC1 (", round(pca_result$sdev[1]^2 / sum(pca_result$sdev^2) * 100, 2), "%)"),
     y = paste0("PC2 (", round(pca_result$sdev[2]^2 / sum(pca_result$sdev^2) * 100, 2), "%)")) +
 annotate("text", x = -6.9, y = 4.5, label = "A", size = 8, family = "Times New Roman")

print(lbiplot_PC12)
#Save the plot as an SVG file
ggsave(filename = "PCA_1_2_Perennial.svg", plot = lbiplot_PC12, device = "svg", width = 8, height = 6, units = "in")

##########PLOTTING LIFE FORM ###############
#################PC1 TO 3#######################
scaling_factor <- max(abs(c(scores_df$PC1, scores_df$PC3))) / max(abs(c(loadings_df$PC1, loadings_df$PC3))) * 0.8
#Create a new column without underscores in Non_perennial
scores_df$Life_history_clean <- gsub("_", " ", scores_df$Life_history)

life_history_shapes <- c("Non perennial" = 16, "Perennial" = 17)  # 16 for circle, 17 for triangle

lbiplot_PC13 <-  ggplot(data = scores_df, aes(x = PC1, y = PC3)) +
  geom_segment(data = loadings_df, aes(x = 0, y = 0, xend = PC1 * scaling_factor, yend = PC3 * scaling_factor),
               arrow = arrow(length = unit(0.05, "inches")), color = "darkgray", size = 1.5) +
  geom_point(aes(color = Life_history_clean, shape = Life_history_clean), size = 6) + 
  stat_ellipse(aes(color = Life_history_clean, linetype = Life_history_clean), level = 0.95, size = 1.5) +
  scale_linetype_manual(name = "Life history", values = c("Non perennial" = "dashed", "Perennial" = "solid"),
                        labels = c("Non perennial", "Perennial")) +
  geom_text(data = loadings_df, aes(x = PC1 * scaling_factor, y = PC3 * scaling_factor, label = variable), 
            vjust = 0.6, hjust = 0.6, color = "black", size = 3.5, family = "Times New Roman") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  scale_color_manual(name = "Life history", values = c("Perennial" = "orange", "Non perennial" = "darkgreen"),
                     labels = c("Non perennial", "Perennial")) +
  scale_shape_manual(name = "Life history", values = life_history_shapes, 
                     labels = c("Non perennial", "Perennial")) +  
  xlim(-7, 5.5) + ylim(-5, 4.5) +
  theme_minimal() +
  theme(
    text = element_text(family = "Times New Roman", size = 22, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    legend.position = "right",
    axis.title = element_text(family = "Times New Roman", size = 20, color = "black"),
    axis.text = element_text(family = "Times New Roman", size = 20, color = "black"),
    legend.text = element_text(family = "Times New Roman", size = 18, color = "black")
  ) +
  labs(x = paste0("PC1 (", round(pca_result$sdev[1]^2 / sum(pca_result$sdev^2) * 100, 2), "%)"),
       y = paste0("PC3 (", round(pca_result$sdev[3]^2 / sum(pca_result$sdev^2) * 100, 2), "%)")) +
  annotate("text", x = -6.9, y = 4.5, label = "B", size = 8, family = "Times New Roman")

print(lbiplot_PC13)
#Save the plot as an SVG file
ggsave(filename = "PCA_1_3_Perennial.svg", plot = lbiplot_PC13, device = "svg", width = 8, height = 6, units = "in")


##########PLOTTING LIFE FORM ###############
#################PC2 TO 3#######################
scaling_factor <- max(abs(c(scores_df$PC2, scores_df$PC3))) / max(abs(c(loadings_df$PC2, loadings_df$PC3))) * 0.8
#Create a new column without underscores in Non_perennial
scores_df$Life_history_clean <- gsub("_", " ", scores_df$Life_history)

life_history_shapes <- c("Non perennial" = 16, "Perennial" = 17)  # 16 for circle, 17 for triangle

lbiplot_PC23 <-  ggplot(data = scores_df, aes(x = PC3, y = PC2)) +
  geom_segment(data = loadings_df, aes(x = 0, y = 0, xend = PC3 * scaling_factor, yend = PC2 * scaling_factor),
               arrow = arrow(length = unit(0.05, "inches")), color = "darkgray", size = 1.5) +
  geom_point(aes(color = Life_history_clean, shape = Life_history_clean), size = 6) + 
  stat_ellipse(aes(color = Life_history_clean, linetype = Life_history_clean), level = 0.95, size = 1.5) +
  scale_linetype_manual(name = "Life history", values = c("Non perennial" = "dashed", "Perennial" = "solid"),
                        labels = c("Non perennial", "Perennial")) +
  geom_text(data = loadings_df, aes(x = PC3 * scaling_factor, y = PC2 * scaling_factor, label = variable), 
            vjust = 0.6, hjust = 0.6, color = "black", size = 3.5, family = "Times New Roman") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  scale_color_manual(name = "Life history", values = c("Perennial" = "orange", "Non perennial" = "darkgreen"),
                     labels = c("Non perennial", "Perennial")) +
  scale_shape_manual(name = "Life history", values = life_history_shapes, 
                     labels = c("Non perennial", "Perennial")) +  
  xlim(-7, 5.5) + ylim(-5, 4.5) +
  theme_minimal() +
  theme(
    text = element_text(family = "Times New Roman", size = 22, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    legend.position = "right",
    axis.title = element_text(family = "Times New Roman", size = 20, color = "black"),
    axis.text = element_text(family = "Times New Roman", size = 20, color = "black"),
    legend.text = element_text(family = "Times New Roman", size = 18, color = "black")
  ) +
  labs(x = paste0("PC2 (", round(pca_result$sdev[2]^2 / sum(pca_result$sdev^2) * 100, 2), "%)"),
       y = paste0("PC3 (", round(pca_result$sdev[3]^2 / sum(pca_result$sdev^2) * 100, 2), "%)")) +
  annotate("text", x = -6.9, y = 4.5, label = "C", size = 8, family = "Times New Roman")

print(lbiplot_PC23)
#Save the plot as an SVG file
ggsave(filename = "PCA_1_2_Perennial.svg", plot = lbiplot_PC23, device = "svg", width = 8, height = 6, units = "in")

##########gypsum type ###############
###############pc1 TO 2#########################
#Create a new column without underscores in Non_perennial
scaling_factor <- max(abs(c(scores_df$PC1, scores_df$PC2))) / max(abs(c(loadings_df$PC1, loadings_df$PC2))) * 0.8
scores_df$Hibat_type_clean <- gsub("_", " ", scores_df$Habitat_type)

habitat_shapes <- c("Gypsum" = 16, "Non gypsum" = 17)  # 16 for circle, 17 for triangle

gbiplot_PC12 <-  ggplot(data = scores_df, aes(x = PC1, y = PC2)) +
  geom_segment(data = loadings_df, aes(x = 0, y = 0, xend = PC1 * scaling_factor, yend = PC2 * scaling_factor),
               arrow = arrow(length = unit(0.05, "inches")), color = "darkgray", size = 1) +
  geom_point(aes(color = Hibat_type_clean, shape = Hibat_type_clean), size = 5) + 
  stat_ellipse(aes(color = Hibat_type_clean, linetype = Hibat_type_clean), level = 0.95, size = 1) +
  scale_linetype_manual(name = "Habitat type", values = c("Non gypsum" = "solid", "Gypsum" = "dashed"),
                        labels = c("Non gypsum", "Gypsum")) +
  geom_text(data = loadings_df, aes(x = PC1 * scaling_factor, y = PC2 * scaling_factor, label = variable), 
            vjust = 0.6, hjust = 0.6, color = "black", size = 3.5, family = "Times New Roman") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(name = "Habitat type", values = c("Gypsum" = "darkgreen", "Non gypsum" = "orange"),
                     labels = c("Non gypsum", "Gypsum")) +
  scale_shape_manual(name = "Habitat type", values = habitat_shapes, 
                     labels = c("Non gypsum", "Gypsum")) +  
  xlim(-5.5, 5.5) + ylim(-4.7, 4.8) +
  theme_minimal() +
  theme(
    text = element_text(family = "Times New Roman", size = 15, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    legend.position = "right",
    axis.title = element_text(family = "Times New Roman", size = 15, color = "black"),
    axis.text = element_text(family = "Times New Roman", size = 15, color = "black"),
    legend.text = element_text(family = "Times New Roman", size = 15, color = "black")
  ) +
  labs(x = paste0("PC1 (", round(pca_result$sdev[1]^2 / sum(pca_result$sdev^2) * 100, 2), "%)"),
       y = paste0("PC2 (", round(pca_result$sdev[2]^2 / sum(pca_result$sdev^2) * 100, 2), "%)")) 
  #annotate("text", x = -6.9, y = 4.8, label = "D", size = 8, family = "Times New Roman")

print(gbiplot_PC12)
#Save the plot as an SVG file
ggsave(filename = "PCA_1_2_gypsum.svg", plot = gbiplot_PC12, device = "svg", width = 8, height = 6, units = "in")

##########gypsum type ###############
###############pc1 TO 3#########################
#Create a new column without underscores in Non_perennial
scaling_factor <- max(abs(c(scores_df$PC1, scores_df$PC3))) / max(abs(c(loadings_df$PC1, loadings_df$PC3))) * 0.8
scores_df$Hibat_type_clean <- gsub("_", " ", scores_df$Habitat_type)

habitat_shapes <- c("Gypsum" = 16, "Non gypsum" = 17)  # 16 for circle, 17 for triangle

gbiplot_PC13 <-  ggplot(data = scores_df, aes(x = PC1, y = PC3)) +
  geom_segment(data = loadings_df, aes(x = 0, y = 0, xend = PC1 * scaling_factor, yend = PC3 * scaling_factor),
               arrow = arrow(length = unit(0.05, "inches")), color = "darkgray", size = 1.5) +
  geom_point(aes(color = Hibat_type_clean, shape = Hibat_type_clean), size = 6) + 
  stat_ellipse(aes(color = Hibat_type_clean, linetype = Hibat_type_clean), level = 0.95, size = 1.5) +
  scale_linetype_manual(name = "Habitat type", values = c("Non gypsum" = "solid", "Gypsum" = "dashed"),
                        labels = c("Non gypsum", "Gypsum")) +
  geom_text(data = loadings_df, aes(x = PC1 * scaling_factor, y = PC3 * scaling_factor, label = variable), 
            vjust = 0.6, hjust = 0.6, color = "black", size = 3.5, family = "Times New Roman") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(name = "Habitat type", values = c("Gypsum" = "darkgreen", "Non gypsum" = "orange"),
                     labels = c("Non gypsum", "Gypsum")) +
  scale_shape_manual(name = "Habitat type", values = habitat_shapes, 
                     labels = c("Non gypsum", "Gypsum")) +  
  xlim(-7, 6) + ylim(-4.7, 4.8) +
  theme_minimal() +
  theme(
    text = element_text(family = "Times New Roman", size = 22, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    legend.position = "right",
    axis.title = element_text(family = "Times New Roman", size = 20, color = "black"),
    axis.text = element_text(family = "Times New Roman", size = 20, color = "black"),
    legend.text = element_text(family = "Times New Roman", size = 18, color = "black")
  ) +
  labs(x = paste0("PC1 (", round(pca_result$sdev[1]^2 / sum(pca_result$sdev^2) * 100, 2), "%)"),
       y = paste0("PC3 (", round(pca_result$sdev[3]^2 / sum(pca_result$sdev^2) * 100, 2), "%)")) +
  annotate("text", x = -6.9, y = 4.8, label = "E", size = 8, family = "Times New Roman")

print(gbiplot_PC13)
#Save the plot as an SVG file
ggsave(filename = "PCA_1_3_gypsum.svg", plot = gbiplot_PC13, device = "svg", width = 8, height = 6, units = "in")


##########gypsum type ###############
###############pc2 TO 3#########################
#Create a new column without underscores in Non_perennial
scaling_factor <- max(abs(c(scores_df$PC2, scores_df$PC3))) / max(abs(c(loadings_df$PC2, loadings_df$PC3))) * 0.8
scores_df$Hibat_type_clean <- gsub("_", " ", scores_df$Habitat_type)

habitat_shapes <- c("Gypsum" = 16, "Non gypsum" = 17)  # 16 for circle, 17 for triangle

gbiplot_PC23 <-  ggplot(data = scores_df, aes(x = PC2, y = PC3)) +
  geom_segment(data = loadings_df, aes(x = 0, y = 0, xend = PC2 * scaling_factor, yend = PC3 * scaling_factor),
               arrow = arrow(length = unit(0.05, "inches")), color = "darkgray", size = 1.5) +
  geom_point(aes(color = Hibat_type_clean, shape = Hibat_type_clean), size = 6) + 
  stat_ellipse(aes(color = Hibat_type_clean, linetype = Hibat_type_clean), level = 0.95, size = 1.5) +
  scale_linetype_manual(name = "Habitat type", values = c("Non gypsum" = "solid", "Gypsum" = "dashed"),
                        labels = c("Non gypsum", "Gypsum")) +
  geom_text(data = loadings_df, aes(x = PC2 * scaling_factor, y = PC3 * scaling_factor, label = variable), 
            vjust = 0.6, hjust = 0.6, color = "black", size = 3.5, family = "Times New Roman") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(name = "Habitat type", values = c("Gypsum" = "darkgreen", "Non gypsum" = "orange"),
                     labels = c("Non gypsum", "Gypsum")) +
  scale_shape_manual(name = "Habitat type", values = habitat_shapes, 
                     labels = c("Non gypsum", "Gypsum")) +  
  xlim(-7, 6) + ylim(-4.7, 4.8) +
  theme_minimal() +
  theme(
    text = element_text(family = "Times New Roman", size = 22, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    legend.position = "right",
    axis.title = element_text(family = "Times New Roman", size = 20, color = "black"),
    axis.text = element_text(family = "Times New Roman", size = 20, color = "black"),
    legend.text = element_text(family = "Times New Roman", size = 18, color = "black")
  ) +
  labs(x = paste0("PC2 (", round(pca_result$sdev[2]^2 / sum(pca_result$sdev^2) * 100, 2), "%)"),
       y = paste0("PC3 (", round(pca_result$sdev[3]^2 / sum(pca_result$sdev^2) * 100, 2), "%)")) +
  annotate("text", x = -6.9, y = 4.8, label = "F", size = 8, family = "Times New Roman")

print(gbiplot_PC23)
#Save the plot as an SVG file
ggsave(filename = "PCA_2_3_gypsum.svg", plot = gbiplot_PC23, device = "svg", width = 8, height = 6, units = "in")



library(cowplot)
library(ggplot2)

# Combine the plots into a 3x2 grid with equal heights
pere_gyp_PC123 <- plot_grid(
  lbiplot_PC12, gbiplot_PC12,
  lbiplot_PC13, gbiplot_PC13,
  lbiplot_PC23, gbiplot_PC23,
  ncol = 2, nrow = 3,
  align = "v", rel_heights = c(1, 1, 1)
)

# Print the combined plot
print(pere_gyp_PC123)

# Save the combined plot as an SVG file
ggsave(filename = "Perennial_gypsumPC123.svg", plot = pere_gyp_PC123, device = "svg", width = 17, height = 17, units = "in")



library(ggplot2)
library(cowplot)

# Adjusting the theme for consistent legend placement
biplot1 <- biplot1 + 
  theme(
    legend.position = "right",
    legend.justification = "center",
    legend.box = "vertical"
  )

biplot2 <- biplot2 + 
  theme(
    legend.position = "right",
    legend.justification = "center",
    legend.box = "vertical"
  )

# Combine the two plots with equal heights
combined_plot <- plot_grid(biplot1, biplot2, ncol = 1, align = "v", rel_heights = c(1, 1))

# Print the combined plot
print(combined_plot)
# Save the combined plot as an SVG file
ggsave(filename = "perennial_gypsum.svg", plot = combined_plot, device = "svg", width = 8, height = 9, units = "in")





# Load necessary libraries
library(ggplot2)
library(tidyr)

# Load the CSV file
data <- read.csv("geo_cbind_taxon_env_plus_trait.csv")

# Define the environmental variable names
bio_variable_names <- c(
  "bio1" = "Annual Mean Temperature",
  "bio2" = "Mean Diurnal Range",
  "bio4" = "Temperature Seasonality",
  "bio8" = "Mean Temperature of Wettest Quarter",
  "bio9" = "Mean Temperature of Driest Quarter",
  "bio12" = "Annual Precipitation",
  "bio14" = "Precipitation of Driest Month",
  "bio15" = "Precipitation Seasonality",
  "bio18" = "Precipitation of Warmest Quarter",
  "bio21" = "Elevation"
)

# Generate boxplots
svg("boxplots_for_staminode_type.svg", width=14, height=10)

par(mfrow=c(3,4), mar=c(4,4,2,1))  # Adjust layout

for (var in names(bio_variable_names)) {
  boxplot(data[[var]] ~ data$Staminode_type,
          main=bio_variable_names[var],
          xlab="Staminode Type",
          ylab="",
          col=c("darkgreen", "orange"))
}

# Close the SVG device
dev.off()

svg("boxplots_for_perennial.svg", width=14, height=10)
# Generate boxplots
par(mfrow=c(3,4), mar=c(4,4,2,1))  # Adjust layout

for (var in names(bio_variable_names)) {
  boxplot(data[[var]] ~ data$Perennial,
          main=bio_variable_names[var],
          xlab="Staminode Type",
          ylab="",
          col=c("darkgreen", "orange"))
}
dev.off()

#########################################
####PGLS for NPP, PETAL, and MAS#########
#########################################
#Load necessary libraries
library(ape)# For phylogenetic analysis
library(caper)# For comparative analysis including PGLS

##Data preparation
pca1 <- scores_df[,1]
#Combine pca scores with taxon and petal size
petal <- cbind(pca1, scores_df[, c("Taxon","Showiness", "Staminode_type", "Average_NPP")])

####Phylogenetic Tree Preparation
#load the phylogenetic tree
tree <- read.tree("BartoniaIngroupPLnames.tre")
plot(tree)#see warning: there are NAs

#Check for NA branch lengths and assign a default value if any NAs are found
if(any(is.na(tree$edge.length))) {
  cat("NA branch lengths found. Assigning default length of 1 to all NA branches.\n")
  tree$edge.length[is.na(tree$edge.length)] <- 1
} else {
  cat("No NA branch lengths found.\n")
}
#Check for matching taxa between the tree and the dataset
library(geiger)

name.check(tree, data$Taxon)

#Prepare the comparative data object for caper
comp_data <- comparative.data(tree, data, "Taxon")#There are duplicated between tips and nodes in phylogeny

####Check and Prepare Taxa and Labels
#Ensure tip labels in the tree are unique and match the taxa in the data
tree$tip.label <- make.unique(as.character(tree$tip.label))
data$Taxon <- make.unique(as.character(data$Taxon))

#Find and report unmatched taxa between the data and the tree
data_not_in_tree <- setdiff(data$Taxon, tree$tip.label)
if (length(data_not_in_tree) > 0) {
  cat("Taxa in the data but not in the tree:\n")
  print(data_not_in_tree)
} else {
  cat("All data taxa have a match in the tree.\n")
}

tree_not_in_data <- setdiff(tree$tip.label, data$Taxon)
if (length(tree_not_in_data) > 0) {
  cat("\nTaxa in the tree but not in the data:\n")
  print(tree_not_in_data)
} else {
  cat("\nAll tree taxa have a match in the data.\n")
}


#Address duplicated labels between tips and nodes
all_labels <- c(tree$tip.label, tree$node.label)
duplicated_labels <- all_labels[duplicated(all_labels)]
if(length(duplicated_labels) > 0) {
  tree$node.label <- make.unique(tree$node.label)
  cat("Duplicated labels in the tree were made unique.\n")
} else {
  cat("No duplicated labels found between tips and nodes.\n")
}

####Comparative Data Object Preparation:PGLS Analysis ####
####PCA1#####
comp_data <- comparative.data(tree, data, "Taxon")

####npp
#check normality###
hist(petal$Average_NPP, probability = TRUE, main = "Histogram for NPP")
hist(log(petal$Average_NPP), probability = TRUE, main = "Log histogram for NPP")

# Shapiro-Wilk test for NPP
shapiro_test_original <- shapiro.test(petal$Average_NPP)
print(shapiro_test_original)

# Shapiro-Wilk test for log-transformed NPP
shapiro_test_log <- shapiro.test(log(petal$Average_NPP))
print(shapiro_test_log)

# Q-Q plot for NPP
qqnorm(petal$Average_NPP, main = "Q-Q Plot for NPP")
qqline(petal$Average_NPP, col = "blue")

# Q-Q plot for log-transformed NPP
qqnorm(log(petal$Average_NPP), main = "Q-Q Plot for Log-transformed NPP")
qqline(log(petal$Average_NPP), col = "blue")

#run pgls
pgls_npp <- pgls(Average_NPP ~ Staminode_type, data = comp_data)
summary(pgls_npp)

# Boxplot for NPP by Staminode type
ggplot(data, aes(x = factor(Staminode_type, levels = c("Absence", "Presence")), y = Average_NPP)) +
  geom_boxplot(fill = "gray") +
  labs(title = "PGLS: P = 0.7976", x = "Staminode type", y = "NPP") +
  theme_minimal()

####petal
#check normality###
hist(petal$petal_size, probability = TRUE, main = "Histogram for Petal Size")
hist(log(petal$petal_size), probability = TRUE, main = "Log histogram for Petal Size")

#run pgls
pgls_petal <- pgls(log(petal_size) ~ pca1, data = comp_data)
summary(pgls_petal)

#add a simple linear regression line
abline(lm(log(petal_size) ~ pca1, data = petal), col = "black")
# Create the plot with the specified customizations
pc1_petal <- ggplot(petal, aes(x = pca1, y = log(petal_size))) +
  geom_smooth(method = "lm", color = "black") +
  geom_point(color = "black", size = 3) +
  labs(x = "PC1", y = "Log mean petal size (mm)") +
  theme_minimal() +
  theme(
    text = element_text(family = "Times New Roman", size = 15),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.title = element_text(family = "Times New Roman", size = 15, color = "black"),
    axis.text = element_text(family = "Times New Roman", size = 15, color = "black")
  ) +
  annotate("text", x = -Inf, y = Inf, label = "A", 
           hjust = -1, vjust = 2, size = 6, family = "Times New Roman")

#Print the plot
print(pc1_petal)

#####mas
#check normality
hist(petal$mas_size, probability = TRUE, main = "Histogram for MAS Size")
hist(log(petal$mas_size), probability = TRUE, main = "Log histogram for MAS Size")

pgls_mas <- pgls(log(mas_size) ~ pca1, data = comp_data)
summary(pgls_mas)

#add a simple linear regression line
abline(lm(log(mas_size) ~ pca1, data = petal), col = "black")
#plot using ggplot
library(ggplot2)
pc1_mas <- ggplot(petal, aes(x = pca1, y = log(mas_size))) +
  geom_smooth(method = "lm", color = "black") +
  geom_point(color = "black", size = 3) +
  labs(x = "PC1", y = "Log mean MAS size (mm)") +
  theme_minimal() +
  theme(
    text = element_text(family = "Times New Roman", size = 15),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.title = element_text(family = "Times New Roman", size = 15, color = "black"),
    axis.text = element_text(family = "Times New Roman", size = 15, color = "black")
  )+
  annotate("text", x = -Inf, y = Inf, label =  "B", 
           hjust = -1, vjust = 2, size = 6, family = "Times New Roman")


print(pc1_mas)

####PCA2#####
##Data preparation
pca2 <- scores_df[,2]
#Combine pca scores with taxon and petal size
petal <- cbind(pca2, scores_df[, c("Taxon","petal_size", "mas_size")])

####petal
#check normality###
hist(petal$petal_size, probability = TRUE, main = "Histogram for Petal Size")
hist(log(petal$petal_size), probability = TRUE, main = "Log histogram for Petal Size")

#run pgls
pgls_petal_pca2 <- pgls(log(petal_size) ~ pca2, data = comp_data)
summary(pgls_petal_pca2)

#add a simple linear regression line
abline(lm(log(petal_size) ~ pca2, data = petal), col = "black")
#plot using ggplot
library(ggplot2)
pc2_petal <- ggplot(petal, aes(x = pca2, y = (log(petal_size)))) +
  geom_smooth(method = "lm", color = "black") +
  geom_point(color = "black", size = 3) +
  labs(x = "PC2", y = "Log mean petal size (mm)") +
  theme_minimal() +
  theme(
    text = element_text(family = "Times New Roman", size = 15),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.title = element_text(family = "Times New Roman", size = 15, color = "black"),
    axis.text = element_text(family = "Times New Roman", size = 15, color = "black")
  )+
  annotate("text", x = -Inf, y = Inf, label = "C", 
           hjust = -1, vjust = 2, size = 6, family = "Times New Roman")

print(pc2_petal)

#####mas
#check normality
hist(petal$mas_size, probability = TRUE, main = "Histogram for MAS Size")
hist(log(petal$mas_size), probability = TRUE, main = "Log histogram for MAS Size")

pgls_mas_pca2 <- pgls(log(mas_size) ~ pca2, data = comp_data)
summary(pgls_mas_pca2)

#add a simple linear regression line
abline(lm(log(mas_size) ~ pca2, data = petal), col = "black")
#plot using ggplot


pc2_mas <- ggplot(petal, aes(x = pca2, y = log(mas_size))) +
  geom_smooth(method = "lm", color = "black") +
  geom_point(color = "black", size = 3) +
  labs(x = "PC2", y = "Log mean MAS size (mm)") +
  theme_minimal() +
  theme(
    text = element_text(family = "Times New Roman", size = 15),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.title = element_text(family = "Times New Roman", size = 15, color = "black"),
    axis.text = element_text(family = "Times New Roman", size = 15, color = "black")
  )+
  annotate("text", x = -Inf, y = Inf, label = "D", 
           hjust = -1, vjust = 2, size = 6, family = "Times New Roman")
print(pc2_mas)





####PCA3#####
##Data preparation
pca3 <- scores_df[,3]
#Combine pca scores with taxon and petal size
petal <- cbind(pca3, scores_df[, c("Taxon","petal_size", "mas_size")])

####petal
#check normality###
hist(petal$petal_size, probability = TRUE, main = "Histogram for Petal Size")
hist(log(petal$petal_size), probability = TRUE, main = "Log histogram for Petal Size")

#run pgls
pgls_petal_pca3 <- pgls(log(petal_size) ~ pca3, data = comp_data)
summary(pgls_petal_pca3)

#add a simple linear regression line
abline(lm(log(petal_size) ~ pca3, data = petal), col = "black")
#plot using ggplot
library(ggplot2)
pc3_petal <- ggplot(petal, aes(x = pca3, y = (log(petal_size)))) +
  geom_smooth(method = "lm", color = "black") +
  geom_point(color = "black", size = 3) +
  labs(x = "PC3", y = "Log mean petal size (mm)") +
  theme_minimal() +
  theme(
    text = element_text(family = "Times New Roman", size = 15),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.title = element_text(family = "Times New Roman", size = 15, color = "black"),
    axis.text = element_text(family = "Times New Roman", size = 15, color = "black")
  )+
  annotate("text", x = -Inf, y = Inf, label = "E", 
           hjust = -1, vjust = 2, size = 6, family = "Times New Roman")

print(pc3_petal)

#####mas
#check normality
hist(petal$mas_size, probability = TRUE, main = "Histogram for MAS Size")
hist(log(petal$mas_size), probability = TRUE, main = "Log histogram for MAS Size")

pgls_mas_pca3 <- pgls(log(mas_size) ~ pca3, data = comp_data)
summary(pgls_mas_pca3)

#add a simple linear regression line
abline(lm(log(mas_size) ~ pca3, data = petal), col = "black")
#plot using ggplot


pc3_mas <- ggplot(petal, aes(x = pca3, y = log(mas_size))) +
  geom_smooth(method = "lm", color = "black") +
  geom_point(color = "black", size = 3) +
  labs(x = "PC3", y = "Log mean MAS size (mm)") +
  theme_minimal() +
  theme(
    text = element_text(family = "Times New Roman", size = 15),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.title = element_text(family = "Times New Roman", size = 15, color = "black"),
    axis.text = element_text(family = "Times New Roman", size = 15, color = "black")
  )+
  annotate("text", x = -Inf, y = Inf, label = "F", 
           hjust = -1, vjust = 2, size = 6, family = "Times New Roman")
print(pc3_mas)

#Combine plots into one figure with labels A, B, C, D
#combined_plot <- grid.arrange(
#  pc1_petal + ggtitle("A"),
#  pc1_mas + ggtitle("B"),
#  pc2_petal + ggtitle("C"),
#  pc2_mas + ggtitle("D"),
#  ncol = 2, nrow = 2
#)
library(gridExtra)
#Combine plots into one figure with labels A, B, C, D
combined <- grid.arrange(
  pc1_petal, pc2_petal,pc3_petal, 
  pc1_mas,
  pc2_mas,
   pc3_mas,
  ncol = 3, nrow = 2
)
print(combined)
ggsave(filename = "combined_plot_petal_mas.svg", plot = combined,
       device = "svg", width = 10, height = 6, bg = "white", units = "in")


#Log transform the petal_size and mas_size
petal$log_petal_size <- log(petal$petal_size)
petal$log_mas_size <- log(petal$mas_size)
# Perform linear regression to obtain correlation and p-value
lm_model <- lm(log_mas_size ~ log_petal_size, data = petal)
lm_summary <- summary(lm_model)
lm_summary
# Extract residuals from the model
residuals_data <- data.frame(residuals = residuals(lm_model))

# Create the density plot for residuals with a vertical line at zero
residual_density_plot <- ggplot(residuals_data, aes(x = residuals)) +
  geom_density(fill = "black", alpha = 0.5) +
  geom_vline(xintercept = 0, color = "black", linetype = "solid") +
  labs(x = "Residuals", y = "Density") +
  theme_minimal()+

theme(
  text = element_text(family = "Times New Roman", size = 15),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(color = "black"),
  panel.border = element_rect(color = "black", fill = NA, size = 0.5),
  axis.title = element_text(family = "Times New Roman", size = 15, color = "black"),
  axis.text = element_text(family = "Times New Roman", size = 15, color = "black")
)
 # annotate("text", x = -Inf, y = Inf, label = "D", 
 #          hjust = -1, vjust = 2, size = 6, family = "Times New Roman")

# Print the plot
print(residual_density_plot)

# Save the plot as SVG
ggsave("residual_density_plot.svg", plot = residual_density_plot, device = "svg", bg = "white")


correlation <- cor(petal$log_petal_size, petal$log_mas_size, method = "pearson")
p_value <- lm_summary$coefficients[2, 4] # p-value for the slope coefficient

print(paste("Correlation between log-transformed petal size and MAS size:", correlation))
print(paste("P-value:", p_value))





########### dont use the one below for manuscript#####

#Calculate the correlation
correlation <- cor(petal$log_petal_size, petal$log_mas_size, method = "pearson")
summary(correlation)
print(paste("Correlation between log-transformed petal size and MAS size:", correlation))

#Create the plot
correlation_plot <- ggplot(petal, aes(x = log_petal_size, y = log_mas_size)) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", color = "red") +
  labs(x = "Log Petal Size", y = "Log MAS Size",
       title = paste("Correlation between Log-Transformed Petal Size and MAS Size\nPearson Correlation:", round(correlation, 2))) +
  theme_minimal()

#Print the plot
print(correlation_plot)

#Save the plot as SVG
ggsave("correlation_plot.svg", plot = correlation_plot, device = "svg", bg = "white")



###plotting#####
####generate a Broken Stick Model plot######
#The Broken Stick model is a method traditionally used to help decide how many principal components to retain in PCA by
#comparing the eigenvalues of the components against a null model where variance is randomly split among components (like breaking a stick into parts).
eigenvalues <- pca_result$sdev^2
proportion_variance <- eigenvalues / sum(eigenvalues)

#Create a Broken Stick Model
broken_stick_model <- function(n) {
  bs <- numeric(n)
  for (i in 1:n) {
    bs[i] <- sum(1/(i:n))
  }
  return(bs / sum(bs))
}

#Calculate the Broken Stick values
broken_stick_values <- broken_stick_model(length(eigenvalues))

#Make a data frame for plotting
scree_data <- data.frame(
  PC = seq_along(eigenvalues),
  Eigenvalue = proportion_variance,
  Broken_Stick = broken_stick_values
)

#Create the scree plot with Broken Stick Model
scree_plot <- ggplot(scree_data, aes(x = PC)) +
  geom_bar(aes(y = Eigenvalue, fill = "Mean of igenvalues"), stat = "identity") +
  geom_line(aes(y = Broken_Stick, color = "Broken stick"), size = 1) +
  geom_point(aes(y = Broken_Stick, color = "Broken stick"), size = 2, shape = 1) +
  scale_fill_manual(values = "grey", labels = "Mean of eigenvalues") +
  scale_color_manual(values = "black", labels = "Broken stick") +
  theme_minimal() +
  theme(
    text = element_text(family = "serif", size = 15, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    legend.position = "right",
    axis.title = element_text(family = "serif", size = 13, color = "black"),  #Adjusting Axis titles
    axis.text = element_text(family = "serif", size = 13, color = "black"),  #Adjusting Axis text (ticks)
  )+
  labs(x = "Principal components",
       y = "Eigenvalue")
#Print the scree plot
print(scree_plot)

#Save the plot
ggsave("PCA_scree_plot_with_broken_stick.svg", scree_plot, width = 10, height = 6, dpi = 300)








################random code###############
#1. Check if the Data is Positive
#Check for negative values
any_negative_values <- any(scores_df$Showiness < 0)
print(any_negative_values)
#Check for zero values
any_zero_values <- any(scores_df$Showiness == 0)
print(any_zero_values)
#Summary statistics
summary(scores_df$Showiness)

#2.Check if the Data is Skewed
#Plot the histogram to visually inspect skewness
hist(scores_df$Showiness, breaks = 20, main = "Histogram of Showiness", xlab = "Showiness")
#Plot the Q-Q plot to inspect skewness
qqnorm(scores_df$Showiness, main = "Q-Q Plot of Showiness")
qqline(scores_df$Showiness, col = "red")
#Calculate skewness using the 'e1071' package
library(e1071)
showiness_skewness <- skewness(scores_df$Showiness)
print(showiness_skewness)

#3Check if the Data is Continuous
#Check for unique values to determine if the data is continuous
unique_values_count <- length(unique(scores_df$Showiness))
print(unique_values_count)

#If the number of unique values is large, the data is likely continuous
if (unique_values_count > 20) {
  print("The Showiness data is likely continuous.")
} else {
  print("The Showiness data may not be continuous.")
}