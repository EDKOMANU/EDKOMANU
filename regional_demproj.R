library(dplyr)

# Global parameters for the regional demographic projection model

time_step <- 1  # For 1-year steps
projection_horizon <- 30 # Default projection horizon in years
age_groups <- c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80-84", "85+") # Standard 5-year age groups
SRB <- 1.05  # Sex Ratio at Birth (males per female)

# Helper functions

# Function to load CSV data and check for required columns
#
# Args:
#   file_path (string): The path to the CSV file.
#   required_cols (vector of strings): A vector of column names that must be present.
#
# Returns:
#   A data frame containing the loaded CSV data.
#   Stops with an error if the file is not found or if any required column is missing.
load_csv_data <- function(file_path, required_cols) {
  if (!file.exists(file_path)) {
    stop(paste("Error: File not found -", file_path))
  }
  df <- read.csv(file_path, stringsAsFactors = FALSE)
  
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(paste("Error: Missing required columns in", file_path, ":", paste(missing_cols, collapse = ", ")))
  }
  
  # TODO: Add more sophisticated validation (e.g., for age groups, regions)
  return(df)
}

# Function to get ASFR rate
#
# Args:
#   asfr_data_df (data frame): Data frame containing ASFR data.
#   target_region (string): The specific region to filter by.
#   target_age_group (string): The specific age group to filter by.
#
# Returns:
#   Numeric ASFR rate for the specified region and age group.
get_asfr_rate <- function(asfr_data_df, target_region, target_age_group) {
  rate <- asfr_data_df %>%
    filter(region == target_region, age_group == target_age_group) %>%
    pull(asfr_rate)
  
  if (length(rate) == 0) {
    # For now, assume found. Could return NA or stop with error.
    # stop(paste("ASFR rate not found for region:", target_region, "age group:", target_age_group))
    return(NA) 
  }
  return(as.numeric(rate))
}

# Function to get survival rate
#
# Args:
#   survival_data_df (data frame): Data frame containing survival rate data.
#   target_region (string): The specific region to filter by.
#   target_age_group (string): The specific age group to filter by.
#   target_sex (string): The specific sex to filter by.
#
# Returns:
#   Numeric survival rate for the specified region, age group, and sex.
get_survival_rate_value <- function(survival_data_df, target_region, target_age_group, target_sex) {
  rate <- survival_data_df %>%
    filter(region == target_region, age_group == target_age_group, sex == target_sex) %>%
    pull(survival_rate)

  if (length(rate) == 0) {
    # For now, assume found. Could return NA or stop with error.
    # stop(paste("Survival rate not found for region:", target_region, "age group:", target_age_group, "sex:", target_sex))
    return(NA)
  }
  return(as.numeric(rate))
}

# Function to get net migration value
#
# Args:
#   net_migration_data_df (data frame): Data frame containing net migration data.
#   target_region (string): The specific region to filter by.
#   target_age_group (string): The specific age group to filter by.
#   target_sex (string): The specific sex to filter by.
#   target_projection_year (numeric): The specific projection year.
#
# Returns:
#   Numeric net migration value. Looks for specific year, then fallback to year 0, then 0 if not found.
get_net_migration_value <- function(net_migration_data_df, target_region, target_age_group, target_sex, target_projection_year) {
  # Try to find migration data for the specific projection year
  migrants <- net_migration_data_df %>%
    filter(region == target_region, 
           age_group == target_age_group, 
           sex == target_sex,
           year == target_projection_year) %>%
    pull(net_migrants)

  if (length(migrants) > 0) {
    return(as.numeric(migrants))
  }
  
  # If not found, try to find data for year 0 (constant rate)
  migrants_year0 <- net_migration_data_df %>%
    filter(region == target_region, 
           age_group == target_age_group, 
           sex == target_sex,
           year == 0) %>%
    pull(net_migrants)
    
  if (length(migrants_year0) > 0) {
    return(as.numeric(migrants_year0))
  }

  # If neither specific year nor year 0 is found, assume 0 net migration
  return(0)
}

# --- Core Demographic Projection Functions ---

# Define age group subsets based on the global 'age_groups'
# These are assumed to be available for the functions below.
# Ensure age_groups is defined before these lines if sourcing the script.
# For example: age_groups <- c("0-4", ..., "85+")

# female_reproductive_ages_vec <- age_groups[4:10] # Typically 15-19 to 45-49
# first_age_group_name <- age_groups[1]
# open_age_group_name <- age_groups[length(age_groups)]

# Note: The actual definitions will be uncommented and potentially adjusted 
# once 'age_groups' is confirmed loaded and has the expected structure.
# For now, these will be passed as arguments or defined locally if needed.


# Function to calculate survivors from a population cohort
#
# Args:
#   current_pop_df (data frame): Population data for current year (cols: age_group, sex, population).
#   survival_rates_df (data frame): Full survival rates data.
#   region_name (string): Name of the region being projected.
#   current_projection_year (numeric): The current projection year (for context, not always used by get_survival_rate_value).
#   age_groups_vec (vector of strings): Vector of all age group names.
#   open_age_group_name (string): Name of the open-ended age group (e.g., "85+").
#
# Returns:
#   A data frame with columns: age_group (original), sex, survived_population.
calculate_survivors <- function(current_pop_df, survival_rates_df, region_name, current_projection_year, age_groups_vec, open_age_group_name) {
  survivors_output <- data.frame(age_group = character(), sex = character(), survived_population = numeric())

  for (s in c("male", "female")) {
    # Loop through all age groups except the open-ended one
    for (i in 1:(length(age_groups_vec) - 1)) {
      current_age_group <- age_groups_vec[i]
      
      pop_cohort_df <- current_pop_df %>% 
        filter(age_group == current_age_group, sex == s)
      
      pop_cohort <- if(nrow(pop_cohort_df) > 0) pop_cohort_df$population else 0
      if (is.na(pop_cohort)) pop_cohort <- 0

      survival_rate <- get_survival_rate_value(survival_rates_df, region_name, current_age_group, s) # current_projection_year not used by current getter
      if(is.na(survival_rate)) survival_rate <- 0 # Default if not found

      survived_pop <- pop_cohort * survival_rate
      
      survivors_output <- rbind(survivors_output, data.frame(age_group = current_age_group, sex = s, survived_population = survived_pop))
    }

    # Handle the open-ended age group
    # Survivors *aging into* the open age group from the last closed age group
    last_closed_age_group <- age_groups_vec[length(age_groups_vec) - 1]
    
    pop_last_closed_df <- current_pop_df %>% 
      filter(age_group == last_closed_age_group, sex == s)
    pop_last_closed <- if(nrow(pop_last_closed_df) > 0) pop_last_closed_df$population else 0
    if (is.na(pop_last_closed)) pop_last_closed <- 0
    
    sr_last_closed <- get_survival_rate_value(survival_rates_df, region_name, last_closed_age_group, s)
    if(is.na(sr_last_closed)) sr_last_closed <- 0

    survivors_into_open <- pop_last_closed * sr_last_closed
    
    # Survivors *already in* the open age group who survive
    pop_open_age_df <- current_pop_df %>% 
      filter(age_group == open_age_group_name, sex == s)
    pop_open_age <- if(nrow(pop_open_age_df) > 0) pop_open_age_df$population else 0
    if (is.na(pop_open_age)) pop_open_age <- 0

    sr_open_age <- get_survival_rate_value(survival_rates_df, region_name, open_age_group_name, s)
    if(is.na(sr_open_age)) sr_open_age <- 0
    
    survivors_in_open <- pop_open_age * sr_open_age
    
    total_open_age_survivors <- survivors_into_open + survivors_in_open 
    # The prompt seems to imply survivors_into_open are those from P[t, 80-84,s] who survive to become P[t+1, 85+,s]
    # and survivors_in_open are those from P[t, 85+,s] who survive to remain P[t+1, 85+,s].
    # However, the current function `calculate_survivors` calculates P[t,age,s] * SR[t,age,s] = S[t,age,s]
    # The aging step comes later. So, the `survivors_into_open` is actually the survivors of the `last_closed_age_group`
    # who are *still in that same age group* after survival is applied, but *before aging*.
    # The prompt's description for `total_open_age_survivors` is more aligned with what happens *after* the aging step.
    # For this function, we should calculate survivors of the open age group *as is*.
    # The aging of `survivors_into_open` into the open group is handled by `age_population_and_add_newborns`.
    # Therefore, for this function, we only need survivors *within* the current open age group.
    # Re-evaluating: The prompt states: "Add data.frame(age_group = open_age_group_name, sex = s, survived_population = total_open_age_survivors)".
    # This means the survived_population for the open_age_group_name *in the output of this function* should be total_open_age_survivors.
    # This is a bit confusing. Let's stick to the prompt's direct instruction for the sum here.
    # The `age_population_and_add_newborns` will then need to be careful about which population it takes for the open age group.

    # Correct interpretation for this function:
    # We need to calculate the survived population for each original age group.
    # For the open age group, this means `pop_open_age * sr_open_age`.
    # The aging of the `last_closed_age_group`'s survivors into the open age group happens in the `age_population_and_add_newborns` function.
    # Let's adjust the logic for the open age group based on this understanding.
    
    # Re-calculating only for those *in* the open age group currently
    pop_open_age_current_df <- current_pop_df %>% 
      filter(age_group == open_age_group_name, sex == s)
    pop_open_age_current <- if(nrow(pop_open_age_current_df) > 0) pop_open_age_current_df$population else 0
    if(is.na(pop_open_age_current)) pop_open_age_current <- 0
    
    sr_open_age_current <- get_survival_rate_value(survival_rates_df, region_name, open_age_group_name, s)
    if(is.na(sr_open_age_current)) sr_open_age_current <- 0
    
    survived_pop_open_age_current <- pop_open_age_current * sr_open_age_current
    
    survivors_output <- rbind(survivors_output, data.frame(age_group = open_age_group_name, sex = s, survived_population = survived_pop_open_age_current))
  }
  return(survivors_output)
}


# Function to calculate births
#
# Args:
#   current_pop_df (data frame): Population data for current year (cols: age_group, sex, population).
#   asfr_df (data frame): Full ASFR data.
#   region_name (string): Name of the region.
#   current_projection_year (numeric): Current projection year (for context).
#   srb_ratio (numeric): Sex ratio at birth (males per female).
#   female_reproductive_ages_vec (vector of strings): Age groups considered reproductive for females.
#
# Returns:
#   A named list: list(male_births = X, female_births = Y).
calculate_births <- function(current_pop_df, asfr_df, region_name, current_projection_year, srb_ratio, female_reproductive_ages_vec) {
  total_births <- 0
  
  fem_pop_reproductive <- current_pop_df %>%
    filter(sex == "female", age_group %in% female_reproductive_ages_vec)
  
  if (nrow(fem_pop_reproductive) > 0) {
    for (r_idx in 1:nrow(fem_pop_reproductive)) {
      row <- fem_pop_reproductive[r_idx, ]
      pop_f_age <- row$population
      if(is.na(pop_f_age)) pop_f_age <- 0

      # current_projection_year is not directly used by get_asfr_rate in current setup, but good for context
      asfr_rate_val <- get_asfr_rate(asfr_df, region_name, row$age_group) 
      if(is.na(asfr_rate_val)) asfr_rate_val <- 0 # Default if not found
      
      total_births <- total_births + (pop_f_age * asfr_rate_val)
    }
  }
  
  male_b <- total_births * srb_ratio / (1 + srb_ratio)
  female_b <- total_births / (1 + srb_ratio)
  
  return(list(male_births = male_b, female_births = female_b))
}


# Function to apply net migration to survived population
#
# Args:
#   pop_df_survived (data frame): Survived population (cols: age_group, sex, survived_population).
#   net_migration_df (data frame): Full net migration data.
#   region_name (string): Name of the region.
#   current_projection_year (numeric): Current projection year.
#   age_groups_vec (vector of strings): Vector of all age group names (unused in current logic but good for context).
#
# Returns:
#   Data frame with columns: age_group, sex, pop_after_migration.
apply_net_migration <- function(pop_df_survived, net_migration_df, region_name, current_projection_year, age_groups_vec) {
  # Ensure pop_after_migration column is initialized correctly from survived_population
  pop_df_after_migration <- pop_df_survived %>%
    mutate(pop_after_migration = survived_population) 

  if(nrow(pop_df_after_migration) == 0) {
    return(pop_df_after_migration %>% select(age_group, sex, pop_after_migration)) # Return empty df with correct structure
  }

  for (i in 1:nrow(pop_df_after_migration)) {
    row <- pop_df_after_migration[i, ]
    
    migrants <- get_net_migration_value(net_migration_df, region_name, row$age_group, row$sex, current_projection_year)
    if(is.na(migrants)) migrants <- 0 # Default if not found
    
    # Ensure survived_population is not NA
    current_survived_pop <- if(is.na(row$survived_population)) 0 else row$survived_population

    pop_df_after_migration[i, "pop_after_migration"] <- current_survived_pop + migrants
  }
  
  return(pop_df_after_migration %>% select(age_group, sex, pop_after_migration))
}


# Function to age the population by one time step and add newborns
#
# Args:
#   pop_after_survival_migration_df (data frame): Population after survival and migration (cols: age_group, sex, pop_after_migration).
#   births_list_result (list): Output from calculate_births (male_births, female_births).
#   survival_rates_df (data frame): Full survival rates data (for newborn survival).
#   region_name (string): Name of the region.
#   current_projection_year (numeric): Current projection year (for context, e.g. for newborn survival rate year).
#   age_groups_vec (vector of strings): Vector of all age group names.
#   first_age_group_name (string): Name of the first age group (e.g., "0-4").
#   open_age_group_name (string): Name of the open-ended age group (e.g., "85+").
#
# Returns:
#   Final population data frame for the next year (cols: age_group (new/aged), sex, population).
age_population_and_add_newborns <- function(pop_after_survival_migration_df, 
                                            births_list_result, 
                                            survival_rates_df, 
                                            region_name, 
                                            current_projection_year, 
                                            age_groups_vec, 
                                            first_age_group_name, 
                                            open_age_group_name) {
                                              
  next_year_pop_df <- data.frame(age_group = character(), sex = character(), population = numeric())

  # Newborns for the next year (become the first age group)
  for (s_newborn in c("male", "female")) {
    # Survival rate for newborns is typically for the 0-4 age group over the first year, 
    # but here it's applied to total births to become the 0-4 population.
    # This assumes the survival rate for "0-4" can approximate survival of newborns into that group.
    # More precise would be a specific "birth to 0-4" survival rate if available.
    newborn_survival_rate <- get_survival_rate_value(survival_rates_df, region_name, first_age_group_name, s_newborn)
    if(is.na(newborn_survival_rate)) newborn_survival_rate <- 0 # Default if not found

    num_births <- if (s_newborn == "male") births_list_result$male_births else births_list_result$female_births
    if(is.na(num_births)) num_births <- 0
    
    surviving_newborns <- num_births * newborn_survival_rate
    
    next_year_pop_df <- rbind(next_year_pop_df, 
                              data.frame(age_group = first_age_group_name, sex = s_newborn, population = surviving_newborns))
  }

  # Aging population from non-open age groups
  # Loop up to the second to last age group (e.g., up to "80-84" if "85+" is open)
  # Survivors from age_groups_vec[i] move to age_groups_vec[i+1]
  for (i in 1:(length(age_groups_vec) - 2)) { # Excludes last closed and open_age_group
    current_original_age_group <- age_groups_vec[i]
    next_actual_age_group <- age_groups_vec[i+1] # Population ages into this group

    for (s in c("male", "female")) {
      pop_to_age_df <- pop_after_survival_migration_df %>%
        filter(age_group == current_original_age_group, sex == s)
      
      pop_to_age <- if(nrow(pop_to_age_df) > 0) pop_to_age_df$pop_after_migration else 0
      if (is.na(pop_to_age)) pop_to_age <- 0
      
      next_year_pop_df <- rbind(next_year_pop_df, 
                                data.frame(age_group = next_actual_age_group, sex = s, population = pop_to_age))
    }
  }
  
  # Handle the population moving into the open-ended age group AND those already in it and surviving
  # Population from the last *closed* age group (e.g., 80-84) who survive and age into the open group (e.g., 85+)
  last_closed_age_group <- age_groups_vec[length(age_groups_vec) - 1]
  
  pop_aging_into_open_male <- pop_after_survival_migration_df %>%
    filter(age_group == last_closed_age_group, sex == "male") %>%
    pull(pop_after_migration)
  if(length(pop_aging_into_open_male) == 0 || is.na(pop_aging_into_open_male)) pop_aging_into_open_male <- 0
    
  pop_aging_into_open_female <- pop_after_survival_migration_df %>%
    filter(age_group == last_closed_age_group, sex == "female") %>%
    pull(pop_after_migration)
  if(length(pop_aging_into_open_female) == 0 || is.na(pop_aging_into_open_female)) pop_aging_into_open_female <- 0

  # Population *already in* the open-ended age group who survive and remain in it
  pop_survived_in_open_male <- pop_after_survival_migration_df %>%
    filter(age_group == open_age_group_name, sex == "male") %>%
    pull(pop_after_migration)
  if(length(pop_survived_in_open_male) == 0 || is.na(pop_survived_in_open_male)) pop_survived_in_open_male <- 0
  
  pop_survived_in_open_female <- pop_after_survival_migration_df %>%
    filter(age_group == open_age_group_name, sex == "female") %>%
    pull(pop_after_migration)
  if(length(pop_survived_in_open_female) == 0 || is.na(pop_survived_in_open_female)) pop_survived_in_open_female <- 0

  # Add to next_year_pop_df, they will be summed later
  next_year_pop_df <- rbind(next_year_pop_df, 
                            data.frame(age_group = open_age_group_name, sex = "male", population = pop_aging_into_open_male + pop_survived_in_open_male),
                            data.frame(age_group = open_age_group_name, sex = "female", population = pop_aging_into_open_female + pop_survived_in_open_female))

  # Group by age_group and sex, then summarise population to ensure distinct rows and sum correctly
  # (e.g. newborns and aged population in the first age group for some models, or open age group contributions)
  final_next_year_pop_df <- next_year_pop_df %>%
    group_by(age_group, sex) %>%
    summarise(population = sum(population, na.rm = TRUE), .groups = 'drop') %>%
    # Order the output consistently by age group (as defined in age_groups_vec) and then by sex
    arrange(match(age_group, age_groups_vec), sex) 
    
  return(final_next_year_pop_df)
}

# --- Main Projection Loop ---

# Define additional global constants / parameters
first_age_group_name <- age_groups[1]
open_age_group_name <- age_groups[length(age_groups)]
# Define female reproductive ages based on the 'age_groups' vector
# Assuming standard 5-year groups starting "0-4", "5-9", "10-14", "15-19"...
# Typically 15-19 is the 4th group, 45-49 is the 10th group.
female_reproductive_ages_vec <- age_groups[4:10] 

# Load Data
# Define required columns for validation
req_cols_initial_pop <- c("region", "age_group", "sex", "population")
req_cols_asfr <- c("region", "age_group", "asfr_rate")
req_cols_survival <- c("region", "age_group", "sex", "survival_rate")
req_cols_migration <- c("region", "year", "age_group", "sex", "net_migrants")

# Load data using the helper function
initial_pop_df <- load_csv_data("data/initial_population.csv", req_cols_initial_pop)
asfr_df <- load_csv_data("data/asfr.csv", req_cols_asfr)
survival_rates_df <- load_csv_data("data/survival_rates.csv", req_cols_survival)
net_migration_df <- load_csv_data("data/net_migration.csv", req_cols_migration)

# Initialize Results Storage
all_projections_list <- list()

# Store initial population with a 'year' column
# The first element of the list (index 1) will correspond to the base year (year 0)
initial_pop_df_with_year <- initial_pop_df %>% mutate(year = 0)
all_projections_list[[1]] <- initial_pop_df_with_year

# Get unique region names from initial population
regions <- unique(initial_pop_df$region)

# Main Projection Loop
cat("Starting population projection...\n")

for (proj_year in 1:projection_horizon) { # proj_year = 1 means projecting for Year 1 from Year 0 data
  cat(paste0("Projecting year: ", proj_year, "\n"))
  
  # Get population from the PREVIOUS year. 
  # list index is proj_year for previous year's data (e.g. for proj_year 1, use list index 1 which is year 0)
  previous_year_pop_list_df <- all_projections_list[[proj_year]] 
  
  current_year_pop_storage_list <- list() # To store results for all regions for the current proj_year

  for (reg_name in regions) {
    cat(paste0("  Processing region: ", reg_name, "\n"))
    
    current_pop_region_df <- previous_year_pop_list_df %>% 
                              filter(region == reg_name) %>%
                              # Ensure correct columns for functions, without year and region for calculation steps
                              select(age_group, sex, population) 

    # Core projection steps:
    # 1. Calculate survivors
    survivors_df <- calculate_survivors(
      current_pop_df = current_pop_region_df,
      survival_rates_df = survival_rates_df,
      region_name = reg_name,
      current_projection_year = proj_year, 
      age_groups_vec = age_groups,
      open_age_group_name = open_age_group_name
    )
    
    # 2. Calculate births
    births_list <- calculate_births(
      current_pop_df = current_pop_region_df, 
      asfr_df = asfr_df,
      region_name = reg_name,
      current_projection_year = proj_year,
      srb_ratio = SRB,
      female_reproductive_ages_vec = female_reproductive_ages_vec
    )
    
    # 3. Apply migration (to survivors)
    pop_after_migration_df <- apply_net_migration(
      pop_df_survived = survivors_df, 
      net_migration_df = net_migration_df,
      region_name = reg_name,
      current_projection_year = proj_year, 
      age_groups_vec = age_groups
    )
    
    # 4. Age population and add newborns
    next_year_region_pop_df <- age_population_and_add_newborns(
      pop_after_survival_migration_df = pop_after_migration_df,
      births_list_result = births_list,
      survival_rates_df = survival_rates_df,
      region_name = reg_name,
      current_projection_year = proj_year, 
      age_groups_vec = age_groups,
      first_age_group_name = first_age_group_name,
      open_age_group_name = open_age_group_name
    )
    
    # Add region and year columns for storage
    next_year_region_pop_df_final <- next_year_region_pop_df %>% 
                                  mutate(region = reg_name, year = proj_year) %>% # This is the population at the END of proj_year
                                  select(year, region, age_group, sex, population) # Standard order

    current_year_pop_storage_list[[length(current_year_pop_storage_list) + 1]] <- next_year_region_pop_df_final
  } # End region loop
  
  # Combine results for all regions for the current projection year
  all_regions_current_year_df <- bind_rows(current_year_pop_storage_list)
  # Store in list: index proj_year + 1 will store data for year `proj_year`
  # e.g. year 1 results stored at index 2, year 30 results stored at index 31
  all_projections_list[[proj_year + 1]] <- all_regions_current_year_df

} # End year loop

cat("Projection finished.\n")

# Combine All Results
# all_projections_list contains data for year 0 (at index 1) up to year `projection_horizon` (at index `projection_horizon + 1`)
final_projection_df <- bind_rows(all_projections_list)
cat("All projection results combined.\n")

# Optional: print a summary of the final data frame
# print(head(final_projection_df))
# print(tail(final_projection_df))
# print(summary(final_projection_df))

# --- Output Generation ---
output_file_path <- "output/projected_population.csv"

# Use write.csv to save the final_projection_df
# Set row.names = FALSE to avoid writing R row numbers to the CSV
write.csv(final_projection_df, file = output_file_path, row.names = FALSE)

cat(paste0("Projection results saved to: ", output_file_path, "\n"))
