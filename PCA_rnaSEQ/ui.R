#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
    
    # Application title
    titlePanel("Linear regression with a spline term"),
    
    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            
            h4('Move the sliders below to optimize the R squared value. Remember the closer it is to 1, the better it is'),
            br(),
            
            # Input: Select the random distribution type ----
            radioButtons("Col_grps", "Color By :",
                         c("Activation State" = "Activation_State",
                           "Treatment" = "Treatment",
                           "Time Period" = "Time_hrs")),

            # br() element to introduce extra vertical spacing ----
            br(),
            
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            
            tabsetPanel(type = "pills",
                        tabPanel("Score Plot", plotOutput("Score_Plot", width = "100%")),
                        tabPanel("Clustering", plotOutput("Heatmap", width = "100%" )),
                        tabPanel("Summary fit 2", plotOutput("Scree_Plot", width = "100%"))
            ), br(),
            div( h4(strong('Guide : ')), br(),
                 p(' - Consider the plot in the `Plot` tab, linear models does not actually capture the whole essence of this relationship. Also look at the value of R squared (0.78); which is the percentage variation explained by the model'), br(),
                 p('Play around with the slider values and try to find an optimal spline term for the best fit. Fit is updated each time. Resulting additional parameters can be viewed in the summary tabs'), br(),br()
            )
            
        )
    )
))
