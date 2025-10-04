#!/usr/bin/env Rscript
#
# Wim Otte (w.m.otte@umcutrecht.nl)
# Oct 3, 2025
#
# Sensitivity analysis for prior specifications
# concern about heavy dependence on prior specifications
################################################################################
library( 'readr' )
library( 'brms' )
library( 'ggplot2' )
library( 'dplyr' )
library( 'tidyr' )

# set readr progress bar to silent mode
options( readr.show_progress = FALSE )

################################################################################
# CUSTOM FUNCTIONS                                                             #
################################################################################

###
# Get data
##
get_data <- function()
{
    # read first tab sheet
    infile <- 'out.02.simulate.refractory/data__weekly-seizures-baseline-35__placebo-verum-risk-0.1--0.6.tsv.gz'
    df <- as.data.frame( readr::read_tsv( infile, show_col_types = FALSE ) )

    # make null_group (with only 'baseline' and 'post-baseline'
    df$groupnull <- 'postbaseline'
    df[ df$group == 'baseline', 'groupnull' ] <- 'baseline'

    return( df )
}

###
# Get bayes factor with specific priors
##
get_bayes_factor_with_priors <- function( data, prior_intercept_sd, prior_b_sd, prior_label )
{
    # truncate normal priors for the intercept and slope
    # Use eval(parse()) to substitute numeric values into prior specification
    prior_string_intercept <- paste0( "normal(0, ", prior_intercept_sd, ")" )
    prior_string_b <- paste0( "normal(0, ", prior_b_sd, ")" )

    priors <- c(
        brms::set_prior( prior_string_intercept, class = "Intercept", lb = -10, ub = 10 ),
        brms::set_prior( prior_string_b, class = "b", lb = -5, ub = 5 )
    )

    print( paste0( "  Fitting models with ", prior_label, " priors..." ) )
    print( paste0( "    Intercept prior: ", prior_string_intercept ) )
    print( paste0( "    b prior: ", prior_string_b ) )

    # fit Bayesian Poisson model [H1]
    bmodel_H1 <- brms::brm( y ~ group, data = data,
                            iter = 20000, warmup = 2000,
                            family = poisson(),
                            prior = priors,
                            save_pars = save_pars( all = TRUE ), refresh = 0 )

    # update Bayesian Poisson model [H0]
    bmodel_H0 <- update( bmodel_H1, formula = y ~ groupnull, newdata = data, refresh = 0 )

    # Bayes Factor, through bridge sampling
    BF_brms_bridge <- brms::bayes_factor( bmodel_H1, bmodel_H0, log = TRUE )

    # output
    output <- data.frame(
        day = data$day[ nrow( data ) ],
        cycle = data$cycle[ nrow( data ) ],
        log_bf = BF_brms_bridge$bf,
        prior_label = prior_label,
        prior_intercept_sd = prior_intercept_sd,
        prior_b_sd = prior_b_sd
    )

    return( output )
}

################################################################################
# END CUSTOM FUNCTIONS                                                         #
################################################################################

# set seed
set.seed( 999 )

outdir <- 'out.14.sensitivity.priors'
dir.create( outdir, showWarnings = FALSE )

# get data
df <- get_data()

####################################
# Define different prior specifications to test
####################################

# Original: Normal(0, 100) with bounds
# Test more informative and less informative priors

prior_specs <- data.frame(
    label = c(
        "Very informative (SD=10)",
        "Informative (SD=50)",
        "Original (SD=100)",
        "Weakly informative (SD=200)",
        "Vague (SD=500)"
    ),
    intercept_sd = c( 10, 50, 100, 200, 500 ),
    b_sd = c( 10, 50, 100, 200, 500 )
)

####################################
# Run BF analysis at key timepoints for each prior
####################################

# Key timepoints: end of each cycle, plus day 152 (MCID crossing)
timepoints <- c(
    7 * 16 - 1,  # end of cycle 1
    7 * 24 - 1,  # end of cycle 2
    7 * 32 - 1,  # end of cycle 3
    152,         # MCID crossing day
    7 * 40 - 1,  # end of cycle 4
    nrow( df )   # end of cycle 5
)

all_results <- NULL

for( i in 1:nrow( prior_specs ) )
{
    prior_label <- prior_specs[ i, 'label' ]
    intercept_sd <- prior_specs[ i, 'intercept_sd' ]
    b_sd <- prior_specs[ i, 'b_sd' ]

    print( paste0( "*** Testing prior: ", prior_label, " ***" ) )

    for( timepoint in timepoints )
    {
        print( paste0( "  Timepoint: day ", timepoint ) )

        # evaluation data
        data <- df[ 1:timepoint, ]

        # get Bayes Factor
        current_bf <- get_bayes_factor_with_priors( data, intercept_sd, b_sd, prior_label )

        # add to container
        all_results <- rbind( all_results, current_bf )
    }
}

# Save results
write.csv( all_results, file = paste0( outdir, '/sensitivity_priors_results.csv' ), row.names = FALSE, quote = TRUE )

####################################
# Create visualization
####################################

# Prepare data for plotting
plot_data <- all_results %>%
    mutate(
        prior_label = factor( prior_label, levels = prior_specs$label ),
        timepoint_label = case_when(
            day == 7 * 16 - 1 ~ "End Cycle 1",
            day == 7 * 24 - 1 ~ "End Cycle 2",
            day == 7 * 32 - 1 ~ "End Cycle 3",
            day == 152 ~ "Day 152 (MCID)",
            day == 7 * 40 - 1 ~ "End Cycle 4",
            day == nrow( df ) ~ "End Cycle 5",
            TRUE ~ as.character( day )
        )
    )

# Plot: log BF across timepoints for different priors
p1 <- ggplot( plot_data, aes( x = day, y = log_bf, color = prior_label, group = prior_label ) ) +
    geom_line( size = 1.2 ) +
    geom_point( size = 3 ) +
    geom_hline( yintercept = log( 10 ), linetype = 2, color = 'gray50' ) +
    geom_hline( yintercept = log( 100 ), linetype = 2, color = 'gray70' ) +
    annotate( "text", x = max( plot_data$day ) * 0.95, y = log( 10 ),
              label = "BF=10 (strong)", vjust = -0.5, hjust = 1, size = 3, color = 'gray40' ) +
    annotate( "text", x = max( plot_data$day ) * 0.95, y = log( 100 ),
              label = "BF=100 (decisive)", vjust = -0.5, hjust = 1, size = 3, color = 'gray40' ) +
    scale_color_brewer( palette = "Set1" ) +
    labs(
        x = "Time (day)",
        y = "log(Bayes Factor) in favor of H1",
        color = "Prior specification",
        title = "Sensitivity to Prior Specification"
    ) +
    theme_classic() +
    theme(
        legend.position = 'top',
        legend.title = element_text( face = 'bold' )
    )

ggsave( plot = p1, file = paste0( outdir, '/sensitivity_priors_logBF.png' ), dpi = 600, height = 6, width = 10, bg = 'white' )

# Plot 2: Relative difference from original prior
original_results <- all_results %>%
    filter( prior_label == "Original (SD=100)" ) %>%
    select( day, log_bf_original = log_bf )

plot_data2 <- all_results %>%
    left_join( original_results, by = "day" ) %>%
    mutate(
        log_bf_diff = log_bf - log_bf_original,
        prior_label = factor( prior_label, levels = prior_specs$label )
    ) %>%
    filter( prior_label != "Original (SD=100)" )

p2 <- ggplot( plot_data2, aes( x = day, y = log_bf_diff, color = prior_label, group = prior_label ) ) +
    geom_line( size = 1.2 ) +
    geom_point( size = 3 ) +
    geom_hline( yintercept = 0, linetype = 1, color = 'black' ) +
    scale_color_brewer( palette = "Set1" ) +
    labs(
        x = "Time (day)",
        y = "Difference in log(BF) from original prior",
        color = "Prior specification",
        title = "Deviation from Original Prior Results"
    ) +
    theme_classic() +
    theme(
        legend.position = 'top',
        legend.title = element_text( face = 'bold' )
    )

ggsave( plot = p2, file = paste0( outdir, '/sensitivity_priors_deviation.png' ), dpi = 600, height = 6, width = 10, bg = 'white' )

####################################
# Create summary table
####################################

summary_table <- plot_data %>%
    group_by( prior_label, prior_intercept_sd, prior_b_sd ) %>%
    summarise(
        log_bf_cycle1 = log_bf[ day == 7 * 16 - 1 ],
        log_bf_cycle2 = log_bf[ day == 7 * 24 - 1 ],
        log_bf_cycle3 = log_bf[ day == 7 * 32 - 1 ],
        log_bf_day152 = log_bf[ day == 152 ],
        log_bf_cycle5 = log_bf[ day == max( day ) ],
        .groups = 'drop'
    )

write.csv( summary_table, file = paste0( outdir, '/sensitivity_priors_summary.csv' ), row.names = FALSE, quote = TRUE )

print( "Summary of log(BF) by prior specification:" )
print( summary_table )

print( paste0( "Results saved to: ", outdir ) )
