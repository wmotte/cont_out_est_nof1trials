#!/usr/bin/env Rscript
#
# Wim Otte (w.m.otte@umcutrecht.nl)
# Oct 3, 2025
#
# Sensitivity analysis for treatment efficacy assumptions
# Addresses concern about varying treatment efficacy assumptions
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
# get study design
##
get_study_design <- function( ncycles, ndays_baseline, ndays_verum, ndays_placebo )
{
    data <- NULL

    # first add baseline
    single <- data.frame(
        cycle = rep( 'prior', ndays_baseline ),
        group = c( rep( 'baseline', ndays_baseline ) ) )

    # add to container
    data <- rbind( data, single )

    # now add cycles
    for( ncycle in 1:ncycles )
    {
        single <- data.frame(
            cycle = rep( paste0( 'cycle-', ncycle ), ndays_verum + ndays_placebo ),
            group = c( rep( 'placebo', ndays_placebo ), rep( 'verum', ndays_verum ) ) )

        # add to container
        data <- rbind( data, single )
    }

    # add day index for full design
    data$day <- 1:nrow( data )

    return( data )
}

###
# Add lambda (used for Poisson seizure simulation) to each day
##
add_lambda_to_design <- function( design, weekly_seizures_baseline, prob_drop_verum, prob_drop_placebo )
{
    # get daily seizures
    daily_seizures_baseline <- weekly_seizures_baseline / 7

    # get daily seizures after drop in risk
    daily_seizures_placebo <- daily_seizures_baseline - ( daily_seizures_baseline * prob_drop_placebo )
    daily_seizures_verum <- daily_seizures_baseline - ( daily_seizures_baseline * prob_drop_verum )

    # prepare
    design$lambda <- NA

    design[ design$group %in% 'baseline', 'lambda' ] <- daily_seizures_baseline
    design[ design$group %in% 'placebo', 'lambda' ] <- daily_seizures_placebo
    design[ design$group %in% 'verum', 'lambda' ] <- daily_seizures_verum

    return( design )
}

###
# Simulate y
##
simulate_y <- function( design, nsubjects )
{

    data <- NULL

    for( nsubject in 1:nsubjects )
    {
        sdesign <- design
        sdesign$subject <- paste0( 'subject-', nsubject )
        sdesign$y <- NA


        for( i in 1:nrow( sdesign ) )
        {
            sdesign[ i, 'y' ] <- rpois( 1, lambda = sdesign[ i, 'lambda' ] )

        }

        data <- rbind( data, sdesign )
    }

    return( data )
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

    # calculate probability that reduction exceeds MCID threshold (negative values)
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

# specify seed for random number generation
set.seed( 777 )

outdir <- 'out.15.sensitivity.efficacy'
dir.create( outdir, showWarnings = FALSE )

# meta
nsubjects <- 1
ncycles <- 5
ndays_baseline <- 7 * 8
ndays_verum <- 7 * 4
ndays_placebo <- 7 * 4

# get study design
design <- get_study_design( ncycles, ndays_baseline, ndays_verum, ndays_placebo )

# Baseline seizure rate
weekly_seizures_baseline <- 35  # 5 per day
baseline_seizures_per_day <- 5

# Test different treatment efficacies
# Original: 60% reduction for verum, 10% for placebo
efficacy_scenarios <- data.frame(
    label = c(
        "Efficacy 40%",
        "Efficacy 50%",
        "Efficacy 60% (Original)",
        "Efficacy 70%",
        "Efficacy 80%"
    ),
    prob_drop_verum = c( 0.40, 0.50, 0.60, 0.70, 0.80 ),
    prob_drop_placebo = rep( 0.10, 5 )  # Keep placebo constant
)

# MCID threshold (30% reduction)
mcid_threshold <- -1 * baseline_seizures_per_day * 0.30
prob_threshold <- 0.95

####################################
# Run analysis for each efficacy scenario
####################################

all_results <- NULL

for( scenario_idx in 1:nrow( efficacy_scenarios ) )
{
    label <- efficacy_scenarios[ scenario_idx, 'label' ]
    prob_drop_verum <- efficacy_scenarios[ scenario_idx, 'prob_drop_verum' ]
    prob_drop_placebo <- efficacy_scenarios[ scenario_idx, 'prob_drop_placebo' ]

    print( paste0( "*** Simulating scenario: ", label, " ***" ) )

    # Create lambda design
    designl <- add_lambda_to_design( design, weekly_seizures_baseline, prob_drop_verum, prob_drop_placebo )

    # simulate seizures based on design
    df <- simulate_y( designl, nsubjects )

    # Now run sequential analysis
    bmodel <- NULL
    single_scenario_results <- NULL

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
        single$scenario <- label
        single$true_efficacy <- prob_drop_verum

        # Check stopping criteria
        if( is.na( stopping_day ) && single$prob_exceeds_mcid >= prob_threshold )
        {
            stopping_day <- day
            single$stopped <- TRUE
        } else
        {
            single$stopped <- FALSE
        }

        # add to container
        single_scenario_results <- rbind( single_scenario_results, single )
    }

    # Add summary for this scenario - replicate stopping_day for all rows
    single_scenario_results$stopping_day <- rep( stopping_day, nrow( single_scenario_results ) )

    # add to overall container
    all_results <- rbind( all_results, single_scenario_results )

    print( paste0( "  -> Stopping day: ", ifelse( is.na( stopping_day ), "Never", stopping_day ) ) )
}

# Save results
write.csv( all_results, file = paste0( outdir, '/sensitivity_efficacy_results.csv' ), row.names = FALSE, quote = TRUE )

####################################
# Create visualization
####################################

# Prepare data for plotting
plot_data <- all_results %>%
    mutate(
        scenario = factor( scenario, levels = efficacy_scenarios$label ),
        time = day
    )

# Plot 1: Probability of exceeding MCID over time for different efficacies
p1 <- ggplot( plot_data, aes( x = time, y = prob_exceeds_mcid, color = scenario, group = scenario ) ) +
    geom_line( size = 1.1 ) +
    geom_hline( yintercept = prob_threshold, linetype = 2, color = 'black' ) +
    # Add vertical lines for stopping days
    geom_vline( data = plot_data %>% filter( stopped == TRUE ) %>% distinct( scenario, stopping_day ),
                aes( xintercept = stopping_day, color = scenario ), linetype = 3, alpha = 0.6 ) +
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
        color = "Treatment efficacy",
        title = "Sensitivity to Treatment Efficacy Assumptions"
    ) +
    scale_x_continuous( breaks = seq( from = 0, to = 336, by = 28 ) ) +
    theme_classic() +
    coord_cartesian( expand = FALSE, xlim = c( 80, 338 ), ylim = c( 0, 1 ) ) +
    theme(
        legend.position = 'top',
        legend.title = element_text( face = 'bold' )
    )

ggsave( plot = p1, file = paste0( outdir, '/sensitivity_efficacy_probability.png' ), dpi = 600, height = 6, width = 10, bg = 'white' )

# Plot 2: Seizure reduction estimates over time
p2 <- ggplot( plot_data, aes( x = time, y = median, color = scenario, group = scenario ) ) +
    geom_ribbon( aes( ymin = CI_low, ymax = CI_high, fill = scenario ), alpha = 0.2, color = NA ) +
    geom_line( size = 1.1 ) +
    geom_hline( yintercept = mcid_threshold, linetype = 2, color = 'black', alpha = 0.6 ) +
    # cycle boundaries
    geom_vline( xintercept = ( 7 * 8 ) + 0.5, linetype = 2, colour = 'gray80', alpha = 0.5 ) +
    geom_vline( xintercept = ( 7 * 16 ) + 0.5, linetype = 2, colour = 'gray80', alpha = 0.5 ) +
    geom_vline( xintercept = ( 7 * 24 ) + 0.5, linetype = 2, colour = 'gray80', alpha = 0.5 ) +
    geom_vline( xintercept = ( 7 * 32 ) + 0.5, linetype = 2, colour = 'gray80', alpha = 0.5 ) +
    geom_vline( xintercept = ( 7 * 40 ) + 0.5, linetype = 2, colour = 'gray80', alpha = 0.5 ) +
    scale_color_brewer( palette = "Set1" ) +
    scale_fill_brewer( palette = "Set1" ) +
    labs(
        x = "Time (day)",
        y = "Seizure reduction (verum-induced, per day)",
        color = "Treatment efficacy",
        fill = "Treatment efficacy",
        title = "Seizure Reduction Estimates with Different Treatment Efficacies"
    ) +
    scale_x_continuous( breaks = seq( from = 0, to = 336, by = 28 ) ) +
    scale_y_continuous( breaks = number_ticks( 10 ) ) +
    theme_classic() +
    coord_cartesian( expand = FALSE, xlim = c( 80, 338 ) ) +
    theme(
        legend.position = 'top',
        legend.title = element_text( face = 'bold' )
    )

ggsave( plot = p2, file = paste0( outdir, '/sensitivity_efficacy_reduction.png' ), dpi = 600, height = 6, width = 10, bg = 'white' )

####################################
# Create summary table
####################################

summary_table <- plot_data %>%
    group_by( scenario, true_efficacy ) %>%
    summarise(
        stopping_day = first( stopping_day ),
        .groups = 'drop'
    ) %>%
    mutate(
        days_saved_vs_fixed = ifelse( is.na( stopping_day ), 0, 336 - stopping_day ),
        cycles_completed_at_stop = ifelse( is.na( stopping_day ), 5, floor( ( stopping_day - 56 ) / 56 ) )
    )

write.csv( summary_table, file = paste0( outdir, '/sensitivity_efficacy_summary.csv' ), row.names = FALSE, quote = TRUE )

print( "Summary of stopping days by treatment efficacy:" )
print( summary_table )

print( paste0( "Results saved to: ", outdir ) )
