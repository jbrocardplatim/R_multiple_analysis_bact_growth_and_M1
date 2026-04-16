################################################################################
# TEMPORAL ANALYSIS - SIMPLIFIED VERSION - JANUARY 2026
# J. Brocard (PLATIM) for G. Lamiral & M. Faure (CIRI)
################################################################################

# ==============================================================================
# CONFIGURATION 
# ==============================================================================

# Time periods (in hours post infection)
PERIODS <- list(
  c(0.667, 2.000),   # 40-120 min
  c(2.000, 3.333),   # 120-200 min
  c(3.333, 4.000),   # 200-240 min
  c(4.000, 6.000)    # 240-360 min
)
PERIOD_LABELS <- c("40-120", "120-200", "200-240", "240-360")

# Y-axis limits for slope panels
YLIM_SLOPE_M1 <- c(-0.6, 0.6)      
YLIM_SLOPE_BACT <- c(-0.1, 0.2)

# Y-axis limits for panel A (regression plot on 240-360)
YLIM_REGRESSION_BACT <- c(0.6, 1.1)
YLIM_REGRESSION_M1 <- NULL  # automatic scaling

# Grayscale colors
GRAY_COLORS <- c("#A0A0A0", "#707070")

# ==============================================================================
# PACKAGES
# ==============================================================================

required_packages <- c("ggplot2", "dplyr", "tidyr", "gridExtra")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

# ==============================================================================
# FUNCTIONS
# ==============================================================================

read_and_process <- function(file_path) {
  cat("Reading:", file_path, "\n")
  
  df <- read.delim(file_path, sep = "\t", dec = ",", check.names = FALSE)
  colnames(df)[1] <- "frame"
  
  # Rename columns to unique identifiers
  n_cols <- ncol(df) - 1
  if (n_cols > 0) {
    colnames(df)[-1] <- paste0("Cell_", 1:n_cols)
  }
  
  df_long <- df %>%
    pivot_longer(cols = -frame, names_to = "cell", values_to = "value") %>%
    filter(!is.na(value))
  
  # Time conversion
  df_long$time_hours <- 0.667 + df_long$frame * (2/60)
  df_long$time_minutes <- df_long$time_hours * 60
  
  # Extract condition and data type
  filename <- basename(file_path)
  condition <- sub("_M1\\.txt$|_bact\\.txt$", "", filename)
  df_long$condition <- condition
  
  if (grepl("_M1\\.txt$", filename)) {
    df_long$data_type <- "M1"
  } else if (grepl("_bact\\.txt$", filename)) {
    df_long$data_type <- "bact"
  } else {
    stop("File must end with '_M1.txt' or '_bact.txt'")
  }
  
  # Apply transformations
  if (df_long$data_type[1] == "M1") {
    df_long <- df_long %>%
      mutate(
        value_adjusted = case_when(
          value == 0 ~ 0.001,
          value == 1 ~ 0.999,
          TRUE ~ value
        ),
        value_transformed = log(value_adjusted / (1 - value_adjusted))
      )
  } else {
    df_long <- df_long %>%
      mutate(value_transformed = log10(value + 1))
  }
  
  return(df_long)
}

calculate_slopes <- function(df, periods) {
  results <- list()
  
  for (i in seq_along(periods)) {
    period <- periods[[i]]
    
    df_period <- df %>%
      filter(time_hours >= period[1] & time_hours < period[2])
    
    slopes <- df_period %>%
      group_by(condition) %>%
      do({
        model <- lm(value_transformed ~ time_hours, data = .)
        data.frame(
          slope = coef(model)[2],
          intercept = coef(model)[1],
          se = summary(model)$coefficients[2, 2]
        )
      }) %>%
      ungroup() %>%
      mutate(
        period_index = i,
        period_label = PERIOD_LABELS[i],
        period_start = period[1],
        period_end = period[2]
      )
    
    results[[i]] <- slopes
  }
  
  return(bind_rows(results))
}

create_simplified_figure <- function(df, slopes_df, condition1, condition2, data_type) {
  
  # Filter for these two conditions
  df_subset <- df %>% filter(condition %in% c(condition1, condition2))
  slopes_subset <- slopes_df %>% filter(condition %in% c(condition1, condition2))
  
  # Set factor order
  df_subset$condition <- factor(df_subset$condition, 
                                 levels = c(condition1, condition2))
  slopes_subset$condition <- factor(slopes_subset$condition, 
                                     levels = c(condition1, condition2))
  
  colors <- setNames(GRAY_COLORS[1:2], c(condition1, condition2))
  
  # Labels and limits
  if (data_type == "M1") {
    ylabel_transformed <- "logit(M1)"
    ylim_slope <- YLIM_SLOPE_M1
    ylim_regression <- YLIM_REGRESSION_M1
  } else {
    ylabel_transformed <- "log(intensity + 1)"
    ylim_slope <- YLIM_SLOPE_BACT
    ylim_regression <- YLIM_REGRESSION_BACT
  }
  
  # ===== PANEL A - Regression plot =====
  stats_transformed <- df_subset %>%
    group_by(condition, time_minutes) %>%
    summarise(mean_val = mean(value_transformed, na.rm = TRUE),
              .groups = "drop")
  
  panel_a <- ggplot() +
    geom_point(data = stats_transformed, 
               aes(x = time_minutes, y = mean_val, fill = condition),
               alpha = 0.4, size = 2.5, shape = 21, color = "black", stroke = 0.3) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    scale_x_continuous(breaks = seq(240, 360, 30), limits = c(240, 360)) +
    labs(x = "", y = ylabel_transformed) +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 16, hjust = 0),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 10),
      panel.grid.minor = element_blank()
    ) +
    ggtitle("A")
  
  # Apply Y-axis limits if specified
  if (!is.null(ylim_regression)) {
    panel_a <- panel_a + coord_cartesian(ylim = ylim_regression)
  }
  
  # Add regression lines for last period only (240-360)
  last_period_idx <- length(PERIODS)
  last_period <- PERIODS[[last_period_idx]]
  
  for (cond in c(condition1, condition2)) {
    slope_data <- slopes_subset %>%
      filter(condition == cond, period_index == last_period_idx)
    
    if (nrow(slope_data) > 0 && !is.na(slope_data$slope)) {
      time_seq_hours <- seq(last_period[1], last_period[2], length.out = 50)
      time_seq_min <- time_seq_hours * 60
      pred_vals <- slope_data$intercept + slope_data$slope * time_seq_hours
      
      pred_df <- data.frame(
        time_minutes = time_seq_min, 
        pred = pred_vals, 
        condition = cond
      )
      
      panel_a <- panel_a +
        geom_line(data = pred_df, 
                  aes(x = time_minutes, y = pred, color = condition),
                  linewidth = 1, alpha = 1.0)
    }
  }
  
  # ===== PANEL B - Slopes barplot (last period only: 240-360) =====
  slopes_last <- slopes_subset %>% filter(period_label == "240-360")
  
  panel_b <- ggplot(slopes_last, 
                    aes(x = period_label, y = slope, fill = condition)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8),
             color = "black", linewidth = 0.5, width = 0.7) +
    geom_errorbar(aes(ymin = slope - se, ymax = slope + se),
                  position = position_dodge(width = 0.8), width = 0.3,
                  linewidth = 0.8) +
    geom_hline(yintercept = 0, linewidth = 1) +
    scale_fill_manual(values = colors) +
    labs(x = "minutes post-infection", 
         y = paste0("Slope (Δ", ylabel_transformed, ")/hour")) +
    theme_bw() +
    theme(
      legend.position = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.background = element_rect(fill = "white", color = "black", linewidth = 1),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      plot.title = element_text(face = "bold", size = 16, hjust = 0),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 10),
      panel.grid.minor = element_blank()
    ) +
    coord_cartesian(ylim = ylim_slope) +
    ggtitle("B")
  
  # Combine panels with equal widths
  combined <- grid.arrange(panel_a, panel_b, ncol = 2, widths = c(1, 1))
  
  # Save figure with reduced width (50% of previous)
  output_file <- paste0(condition1, "_", condition2, "_", data_type, "_simplified.png")
  ggsave(output_file, combined, width = 4, height = 4, dpi = 300)
  cat("Figure saved:", output_file, "\n\n")
}

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

cat("\n", strrep("=", 80), "\n")
cat("SIMPLIFIED TEMPORAL ANALYSIS (2 panels per figure)\n")
cat(strrep("=", 80), "\n\n")

# Define all files to analyze
all_files <- c(
  "siCTL_bact.txt",
  "siRNF11_bact.txt",
  "siWWP2_bact.txt",
  "siCTL_M1.txt",
  "siRNF11_M1.txt",
  "siWWP2_M1.txt"
)

# Read and process all data
all_data <- lapply(all_files, read_and_process)
df_all <- bind_rows(all_data)

# Separate by data type
df_bact <- df_all %>% filter(data_type == "bact")
df_m1 <- df_all %>% filter(data_type == "M1")

# Calculate slopes
cat("\nCalculating slopes for bacterial data...\n")
slopes_bact <- calculate_slopes(df_bact, PERIODS)

cat("\nCalculating slopes for M1 data...\n")
slopes_m1 <- calculate_slopes(df_m1, PERIODS)

# Generate 4 figures
cat("\n", strrep("=", 80), "\n")
cat("GENERATING FIGURES\n")
cat(strrep("=", 80), "\n\n")

cat("Figure 1: siCTL vs siRNF11 - bacterial data\n")
create_simplified_figure(df_bact, slopes_bact, "siCTL", "siRNF11", "bact")

cat("Figure 2: siCTL vs siRNF11 - M1 data\n")
create_simplified_figure(df_m1, slopes_m1, "siCTL", "siRNF11", "M1")

cat("Figure 3: siCTL vs siWWP2 - bacterial data\n")
create_simplified_figure(df_bact, slopes_bact, "siCTL", "siWWP2", "bact")

cat("Figure 4: siCTL vs siWWP2 - M1 data\n")
create_simplified_figure(df_m1, slopes_m1, "siCTL", "siWWP2", "M1")

cat(strrep("=", 80), "\n")
cat("ANALYSIS COMPLETE\n")
cat(strrep("=", 80), "\n\n")
