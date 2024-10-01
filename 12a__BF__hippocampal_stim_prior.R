#!/usr/bin/env Rscript
#
# Wim Otte (w.m.otte@umcutrecht.nl)
# Oct 1, 2024
#
################################################################################
library( 'readr' ) # read and write (large) tsv-files
library( 'brms' )
library( 'ggplot2' )
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
    infile <- 'datasets/Tellez-Zenteno.csv'   
    df <- read.csv( infile, stringsAsFactors = TRUE, row.names = 1 )
    df$time <- 1:nrow( df )
    
    rownames( df ) <- NULL
    
    return( df )
}
    

################################################################################
# CUSTOM FUNCTIONS                                                             #
################################################################################

set.seed( 565 )

# outdir
outdir <- 'out.12a.BF.hippocampal.stim.prior'
dir.create( outdir, showWarnings = FALSE )

# get data
head( df <- get_data() )

# select
data <- df[ df$month <= 48, ]

# truncate normal priors for the intercept and slope within reasonable seizure loads 
priors <- c(
    prior( normal( 0, 100 ), class = "Intercept", lb = -10, ub = 10 ),
    prior( normal( 0, 100 ), class = "b", lb = -5, ub = 5 )  
)

# plot normal
curve( dnorm( x, 0, 100 ), from = -500, to = 500, main = "Normal(0, 100)", ylab = "Density" )

# fit Bayesian Poisson model (treatment)
bmodel_H1_prior <- brms::brm( seizures ~ treatment, data = data, 
                        iter = 25000, warmup = 3000,
                        family = poisson(), 
                        sample_prior = "only",  # Sample from priors only
                        prior = priors,
                        save_pars = save_pars( all = TRUE ), refresh = 0 )

 
# posterior post-warm-up samples [5000]
draws <- as.data.frame( bmodel_H1_prior )

median( draws$b_Intercept )
median( draws$b_treatmenton )

min( draws$b_Intercept )
min( draws$b_treatmenton )

max( draws$b_Intercept )
max( draws$b_treatmenton )

# plot

# reshape the data from wide to long format for plotting
draws_long <- draws %>%
    pivot_longer( cols = starts_with( "b_" ), names_to = "Parameter", values_to = "Value" )

# create the plot
p <- ggplot( draws_long, aes(x = Value ) ) +
    geom_density( fill = "blue", alpha = 0.5 ) +
    facet_wrap( ~Parameter, scales = "free" ) +
    labs( x = "Value", y = "Prior distribution (density)" ) +
    theme_minimal()

plot( p )

# save to disk
ggsave( plot = p, file = paste0( outdir, '/prior_distributions.png' ), dpi = 300, height = 4, width = 8, bg = 'white' )

# Assuming 'p' is your original plot
x_limits <- ggplot_build( p )$layout$panel_params[[ 1 ]]$x.range

#################
# Fit with data
#################

# fit Bayesian Poisson model (treatment)
bmodel_H1 <- brms::brm( seizures ~ treatment, data = data, 
                              iter = 25000, warmup = 3000,
                              family = poisson(), 
                             
                              prior = priors,
                              save_pars = save_pars( all = TRUE ), refresh = 0 )

# posterior post-warm-up samples [5000]
draws2 <- as.data.frame( bmodel_H1 )

median( draws2$b_Intercept )
median( draws2$b_treatmenton )

min( draws2$b_Intercept )
min( draws2$b_treatmenton )

max( draws2$b_Intercept )
max( draws2$b_treatmenton )

# plot

# reshape the data from wide to long format for plotting
draws2_long <- draws2 %>%
    pivot_longer( cols = starts_with( "b_" ), names_to = "Parameter", values_to = "Value" )

# create the plot
p2 <- ggplot( draws2_long, aes(x = Value ) ) +
    geom_density( fill = "orange", alpha = 0.5 ) +
    facet_wrap( ~Parameter, scales = "free" ) +
    labs( x = "Value", y = "Posterior distribution (density)" ) +
    theme_minimal() +
    coord_cartesian(xlim = x_limits)  # Apply the same x-axis limits

plot( p2 )

# save to disk
ggsave( plot = p2, file = paste0( outdir, '/posterior_distributions.png' ), dpi = 300, height = 4, width = 8, bg = 'white' )

