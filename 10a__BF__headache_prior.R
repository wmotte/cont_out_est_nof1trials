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
    infile <- 'datasets/SMS.csv'   
    df <- read.csv( infile, stringsAsFactors = TRUE, row.names = 1 )

    df[ df$group == 'injection', 'group' ] <- 'baseline'

    rownames( df ) <- NULL
    
    return( df )
}


###
# Get BF
##
get_bayes_factor <- function( data )
{
    # fit Bayesian model (treatment)
    # refresh = 0 turns off iteration info print
    bmodel_H1 <- brms::brm( delta_severity ~ group, data = data, 
                            prior = prior( normal( 0, 5 ), class = b ),
                            iter = 25000, warmup = 3000,
                            save_pars = save_pars( all = TRUE ), refresh = 0 )
    
    # update Bayesian model [H0], with intercept only
    bmodel_H0 <- update( bmodel_H1, formula = delta_severity ~ 1, refresh = 0 ) 
    
    # Bayes Factor, through bridge sampling
    BF_brms_bridge <- brms::bayes_factor( bmodel_H1, bmodel_H0, log = TRUE )
    
    print( BF_brms_bridge )
    
    # output, last row of data
    output <- data[ nrow( data ), ]
    output$log_bf <- BF_brms_bridge$bf
    
    return( output )
    
}

################################################################################
# CUSTOM FUNCTIONS                                                             #
################################################################################

# set seed
set.seed( 3333 )

outdir <- 'out.10a.BF.headache.prior'
dir.create( outdir, showWarnings = FALSE )

# get data
head( df <- get_data() )

# remove NA's
df <- df[ !is.na( df$delta_severity ), ]

# get until day 4
data <- df[ df$day <= 4, ]

# truncate gamma and normal priors for the intercept and slope within reasonable 
# severity ranges (0-1 for intercept and -1 to +1 for group) 
priors <- c(
    prior( gamma( 2, 3 ), class = "Intercept", lb = 0, ub = 1 ),
    prior( normal( 0, 5 ), class = "b", lb = -1, ub = 1 )  
)

# plot gamma
curve( dgamma( x, shape = 2, rate = 3 ), from = -2, to = 4, main = "Gamma(2, 3)", ylab = "Density" )

# plot normal
curve( dnorm( x, 0, 5 ), from = -25, to = 25, main = "Normal(0, 5)", ylab = "Density" )

# fit Bayesian model (treatment); refresh = 0 turns off iteration info print
bmodel_H1_prior <- brms::brm( delta_severity ~ group, data = data, 
                        iter = 25000, 
                        warmup = 3000,
                        sample_prior = "only",  # Sample from priors only
                        prior = priors,
                        save_pars = save_pars( all = TRUE ), refresh = 0 )


# posterior post-warm-up samples [5000]
draws <- as.data.frame( bmodel_H1_prior )

median( draws$b_Intercept )
median( draws$b_groupSMS )

min( draws$b_Intercept )
min( draws$b_groupSMS )

max( draws$b_Intercept )
max( draws$b_groupSMS )

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

# fit Bayesian model (treatment); refresh = 0 turns off iteration info print
bmodel_H1 <- brms::brm( delta_severity ~ group, data = data, 
                              iter = 25000, 
                              warmup = 3000,
                              
                              prior = priors,
                              save_pars = save_pars( all = TRUE ), refresh = 0 )


# posterior post-warm-up samples [5000]
draws2 <- as.data.frame( bmodel_H1 )

median( draws2$b_Intercept )
median( draws2$b_groupSMS )

min( draws2$b_Intercept )
min( draws2$b_groupSMS )

max( draws2$b_Intercept )
max( draws2$b_groupSMS )

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


