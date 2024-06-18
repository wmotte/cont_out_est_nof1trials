#!/usr/bin/env Rscript
#
# Visualize n-of-1 trial data - headache data
#
# Wim Otte (w.m.otte@umcutrecht.nl)
#
################################################################################
library( 'ggplot2' )        # for plotting

###
# Get data
##
get_data <- function()
{
    # read first tab sheet
    infile <- 'datasets/SMS.csv'   
    df <- read.csv( infile, stringsAsFactors = TRUE, row.names = 1 )
    df$day <- paste0( "day ", df$day )
    
    return( df )
}

### 
# Helper function plot.
##
number_ticks <- function( n ) { function( limits ) pretty( limits, n ) }

################################################################################
# END FUNCTIONS
################################################################################

# output dir
outdir <- 'out.01.plot.headache'
dir.create( outdir, showWarnings = FALSE )

# get data
head( df <- get_data() )

# colors
custom_colours <- c( "gray70", "#E34A33", "#E69F00", "#56B4E9" )

####################################

# plot
p <- ggplot( data = df, aes( x = time, y = severity, group = group, fill = group, colour = group ) ) +
    facet_wrap( ~day, ncol = 1 ) +
    geom_point( size = 3 ) + # for the fat-black strokes
    geom_point( shape = 21, colour = 'gray40', stroke = 0.4, size = 3 ) +
    ylab( "Headache rating (visual analogue scale)" ) +
    xlab( 'Clock time (h)' ) +
    scale_x_continuous( breaks = number_ticks( 16 ) ) +
    scale_y_continuous( breaks = number_ticks( 7 ) ) +
    scale_color_manual( values = custom_colours ) +
    scale_fill_manual( values = custom_colours ) +
    theme_classic() +
    coord_cartesian( expand = FALSE, xlim = c( 8.6, 24.4 ), ylim = c( -0.04, 1.04 ) ) +
    theme( legend.position = 'top', strip.text.x = element_text( face = 'bold' ), 
                                    strip.text.y = element_text( face = 'bold' ) )


# save
ggsave( plot = p, file = paste0( outdir, '/plot_SMS.png' ), dpi = 600, height = 8, width = 4 )

## plot lines

# colors
custom_colours <- c( "gray70", "#E69F00", "#56B4E9" )

df[ df$group == 'injection', 'group' ] <- 'baseline'

# plot
p_black <- ggplot( data = df, aes( x = time, y = severity ) ) +
    facet_wrap( ~day, ncol = 1 ) +
    geom_line( colour = 'gray20' ) +
    geom_point( aes( fill = group ), size = 3 ) + # for the fat-black strokes
    geom_point( aes( fill = group ), shape = 21, colour = 'gray40', stroke = 0.4, size = 3 ) +
    ylab( "Headache rating (visual analogue scale)" ) +
    xlab( 'Clock time (h)' ) +
    scale_x_continuous( breaks = number_ticks( 16 ) ) +
    scale_y_continuous( breaks = number_ticks( 7 ) ) +
    scale_color_manual( values = custom_colours ) +
    scale_fill_manual( values = custom_colours ) +
    labs( fill = 'period' ) +
    theme_classic() +
    coord_cartesian( expand = FALSE, xlim = c( 8.6, 24.4 ), ylim = c( -0.04, 1.04 ) ) +
    theme( legend.position = 'top', strip.text.x = element_text( face = 'bold' ), 
           strip.text.y = element_text( face = 'bold' ) )

p_black
    
# save
ggsave( plot = p_black, file = paste0( outdir, '/plot_SMS_black.png' ), dpi = 600, height = 8, width = 4 )

##############
# delta
##############

# blue + orange
custom_colours <- c( "#E69F00", "#56B4E9" )

# remove 'baseline day'
df_delta <- df[ df$group %in% c( 'placebo', 'SMS' ), ]

# label name 'group' -> 'period'
df_delta$period <- df_delta$group

# plot
p_delta <- ggplot( data = df_delta, aes( x = time, y = delta_severity, fill = period, colour = period ) ) +
    facet_wrap( ~day, ncol = 1 ) +
    geom_hline( yintercept = 0, colour = "gray80", linetype = "dashed" ) +
    geom_point( size = 3 ) + # for the fat-black strokes
    geom_point( shape = 21, colour = 'gray40', stroke = 0.4, size = 3 ) +
    ylab( "∆ headache rating" ) +
    xlab( 'Clock time (h)' ) +
    scale_x_continuous( breaks = number_ticks( 16 ) ) +
    scale_y_continuous( breaks = number_ticks( 7 ) ) +
    scale_color_manual( values = custom_colours ) +
    scale_fill_manual( values = custom_colours ) +
    theme_classic() +
    coord_cartesian( expand = FALSE, xlim = c( 8.6, 24.4 ), ylim = c( -1.1, 1.1 ) ) +
    theme( legend.position = 'top', strip.text.x = element_text( face = 'bold' ), 
           strip.text.y = element_text( face = 'bold' ) )



# save
ggsave( plot = p_delta, file = paste0( outdir, '/plot_SMS_delta.png' ), dpi = 600, height = 7, width = 4 )



#### REPLICATE ####

m <- glm( delta_severity ~ group, data = df )
sum <- summary( m )
emm <- emmeans::emmeans( m, specs = pairwise ~ group, adjust = "none" )

# get separate sets for 1-hour / 3-hour post injection
df_1_h <- df[ df$measurement == 'one', ]
df_3_h <- df[ df$measurement == 'three', ]

# 1 hour
m <- glm( delta_severity ~ group, data = df_1_h )
sum <- summary( m )
emm <- emmeans::emmeans( m, specs = pairwise ~ group, adjust = "none" )

# 3 hours
m <- glm( delta_severity ~ group, data = df_3_h )
sum <- summary( m )
emm <- emmeans::emmeans( m, specs = pairwise ~ group, adjust = "none" )


########## STATS #########


subset <- df[ df$day == 'day 1', ]

m <- glm( severity ~ group, data = subset )
sum <- summary( m )
emm <- emmeans::emmeans( m, specs = pairwise ~ group, adjust = "none" )


subset <- df[ df$day != 'day 0', ]

m <- glm( severity ~ group, data = subset )
sum <- summary( m )
emm <- emmeans::emmeans( m, specs = pairwise ~ group, adjust = "none" )



