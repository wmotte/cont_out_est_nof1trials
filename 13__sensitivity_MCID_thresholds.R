#!/usr/bin/env Rscript
#
# Wim Otte (w.m.otte@umcutrecht.nl)
# Oct 3, 2025
#
# Sensitivity analysis for MCID thresholds
# Addresses concern about exploring stopping rules with varying MCID
################################################################################
library( 'readr' )
library( 'ggplot2' )
library( 'rstanarm' )
library( 'dplyr' )
library( 'tidyr' )

# set readr progress bar to silent mode
options( readr.show_progress = FALSE )

################################################################################
# CUSTOM FUNCTIONS                                                             #
################################################################################

###
# Helper function plot.
##
number_ticks <- function( n ) { function( limits ) pretty( limits, n ) }

###
# Get data
##
get_data <- function()
{
    # read first tab sheet
    infile <- 'out.02.simulate.refractory/data__weekly-seizures-baseline-35__placebo-verum-risk-0.1--0.6.tsv.gz'
    df <- readr::read_tsv( infile, show_col_types = FALSE )

    return( df )
}

###
# get seizure reduction draws and probability of exceeding MCID
##
get_mcid_probability <- function( bmodel, mcid_threshold )
{
    # posterior post-warmup samples [5000]
    draws <- as.data.frame( bmodel )

    # get mean seizure reduction
    draws$seizure_reduction <- exp( draws$`(Intercept)` + draws$groupverum ) - exp( draws$`(Intercept)` + draws$groupplacebo )

    # calculate probability that reduction exceeds MCID threshold
    prob_exceeds_mcid <- mean( draws$seizure_reduction < mcid_threshold )

    # get quantiles
    values <- quantile( draws$seizure_reduction, c( 0.025, 0.5, 0.975 ) )

    out <- data.frame(
        median = values[ 2 ],
        CI_low = values[ 1 ],
        CI_high = values[ 3 ],
        prob_exceeds_mcid = prob_exceeds_mcid
    )

    return( out )
}

################################################################################
# END CUSTOM FUNCTIONS                                                         #
################################################################################

outdir <- 'out.13.sensitivity.MCID'
dir.create( outdir, showWarnings = FALSE )

# get data
df <- get_data()

# Define MCID thresholds to test (as proportion of baseline seizure rate)
# Baseline = 5 seizures/day = 35 seizures/week
baseline_seizures_per_day <- 5

# Test multiple MCID thresholds: 20%, 30%, 40%, reduction
mcid_proportions <- c( 0.20, 0.30, 0.40 )
mcid_thresholds <- -1 * baseline_seizures_per_day * mcid_proportions

# Stopping criteria: probability threshold for decision
prob_threshold <- 0.95  # 95% posterior probability

####################################
# Run analysis for each MCID threshold
####################################

all_results <- NULL

for( mcid_idx in 1:length( mcid_thresholds ) )
{
    mcid_threshold <- mcid_thresholds[ mcid_idx ]
    mcid_prop <- mcid_proportions[ mcid_idx ]

    print( paste0( "*** Analyzing MCID threshold: ", mcid_prop * 100, "% reduction (", mcid_threshold, " seizures/day) ***" ) )

    bmodel <- NULL
    single_mcid_results <- NULL

    # first row with data from both groups
    start_row <- 85
    end_row <- nrow( df )

    stopping_day <- NA

    for( i in start_row:end_row )
    {
        print( i )
        # set seed
        set.seed( 123 + i )

        # get data upto time == i
        data <- df[ 1:i, ]

        day <- df$day[ i ]
        cycle <- df$cycle[ i ]

        if( i == start_row )
        {
            # fit Bayesian Poisson model (treatment)
            bmodel <- rstanarm::stan_glm( y ~ group, data = data, iter = 5000, refresh = 0, family = poisson() )
        } else
        {
            # update Bayesian Poisson model
            bmodel <- update( bmodel, formula = y ~ group, data = data, iter = 5000, refresh = 0, family = poisson() )
        }

        # get MCID probability
        single <- get_mcid_probability( bmodel, mcid_threshold )
        single$day <- day
        single$cycle <- as.character( cycle )
        single$mcid_threshold <- mcid_threshold
        single$mcid_proportion <- mcid_prop

        # Check stopping criteria
        if( is.na( stopping_day ) && single$prob_exceeds_mcid >= prob_threshold )
        {
            stopping_day <- day
            single$stopped <- TRUE
            print( paste0( "  -> Stopping criterion met at day ", stopping_day, " (prob = ", round( single$prob_exceeds_mcid, 3 ), ")" ) )
        } else
        {
            single$stopped <- FALSE
        }

        # add to container
        single_mcid_results <- rbind( single_mcid_results, single )

        # If stopped, break the loop (but still compute for comparison)
        # Actually continue to see what would happen
    }

    # Add summary for this MCID - replicate stopping_day for all rows
    single_mcid_results$stopping_day <- rep( stopping_day, nrow( single_mcid_results ) )

    # add to overall container
    all_results <- rbind( all_results, single_mcid_results )
}

# Save results
#   mcid_label    mcid_threshold  mcid_proportion stopping_day days_saved_vs_fixed cycles_completed_at_stop
#   <chr>                  <dbl>            <dbl>        <dbl>               <dbl>                    <dbl>
# 1 20% reduction           -1               0.2           85                 251                        0
# 2 30% reduction           -1.5             0.3           88                 248                        0
# 3 40% reduction           -2               0.4          332                   4                        4
write.csv( all_results, file = paste0( outdir, '/sensitivity_mcid_results.csv' ), row.names = FALSE, quote = TRUE )

####################################
# Create visualization
####################################

# Prepare data for plotting
plot_data <- all_results %>%
    mutate(
        mcid_label = paste0( mcid_proportion * 100, "% reduction" ),
        time = day
    )

# Plot 1: Probability of exceeding MCID over time for different thresholds
p1 <- ggplot( plot_data, aes( x = time, y = prob_exceeds_mcid, color = mcid_label, group = mcid_label ) ) +
    geom_line( size = 1.1 ) +
    geom_hline( yintercept = prob_threshold, linetype = 2, color = 'black' ) +
    # Add vertical lines for stopping days
    geom_vline( data = plot_data %>% filter( stopped == TRUE ) %>% distinct( mcid_label, stopping_day ),
                aes( xintercept = stopping_day, color = mcid_label ), linetype = 3, alpha = 0.6 ) +
    # cycle boundaries
    geom_vline( xintercept = ( 7 * 8 ) + 0.5, linetype = 2, colour = 'gray80', alpha = 0.5 ) +
    geom_vline( xintercept = ( 7 * 16 ) + 0.5, linetype = 2, colour = 'gray80', alpha = 0.5 ) +
    geom_vline( xintercept = ( 7 * 24 ) + 0.5, linetype = 2, colour = 'gray80', alpha = 0.5 ) +
    geom_vline( xintercept = ( 7 * 32 ) + 0.5, linetype = 2, colour = 'gray80', alpha = 0.5 ) +
    geom_vline( xintercept = ( 7 * 40 ) + 0.5, linetype = 2, colour = 'gray80', alpha = 0.5 ) +
    scale_color_brewer( palette = "Set1" ) +
    labs(
        x = "Time (day)",
        y = "Posterior probability of exceeding MCID",
        color = "MCID threshold",
        title = "Sensitivity to MCID Threshold Choice"
    ) +
    scale_x_continuous( breaks = seq( from = 0, to = 336, by = 28 ) ) +
    theme_classic() +
    coord_cartesian( expand = FALSE, xlim = c( 80, 338 ), ylim = c( 0, 1 ) ) +
    theme(
        legend.position = 'top',
        legend.title = element_text( face = 'bold' )
    )

ggsave( plot = p1, file = paste0( outdir, '/sensitivity_mcid_probability.png' ), dpi = 600, height = 5, width = 8, bg = 'white' )

  
####################################
# Create summary table
####################################

summary_table <- plot_data %>%
    group_by( mcid_label, mcid_threshold, mcid_proportion ) %>%
    summarise(
        stopping_day = first( stopping_day ),
        .groups = 'drop'
    ) %>%
    mutate(
        days_saved_vs_fixed = 336 - stopping_day,
        cycles_completed_at_stop = floor( ( stopping_day - 56 ) / 56 )
    )

write.csv( summary_table, file = paste0( outdir, '/sensitivity_mcid_summary.csv' ), row.names = FALSE, quote = TRUE )

print( "Summary of stopping days by MCID threshold:" )
print( summary_table )

print( paste0( "Results saved to: ", outdir ) )
