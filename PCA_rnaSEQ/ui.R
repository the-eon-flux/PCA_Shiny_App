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
    titlePanel("Explore RNA seq data with PCA"),
    
    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            
            h4('Select from the sample grouping & principal component options to change the score plot view'),
            br(),
            
            # Input: Select the color group
            radioButtons("Col_grps", "Color By :",
                         c("Activation State" = "Activation_State",
                           "Treatment" = "Treatment",
                           "Time Period (hrs)" = "Time_hrs")),br(),

            # Select the PC for plot            
            h4('Y Axis Principal component'), br(),
            selectInput("yPC", "Choose a Principal Component :",
                        choices = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18"), selected = "PC2"), br(),
            
            h4('X Axis Principal component'), br(),
            selectInput("xPC", "Choose a Principal Component :",
                        choices = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18"), selected = "PC1"),
            
            
        ),
        
        # Creating tabs for multiple plots
        mainPanel(
            
            tabsetPanel(type = "pills",
                        tabPanel("Score Plot", plotOutput("Score_Plot", width = "50%")),
                        tabPanel("Clustering", plotOutput("Heatmap", width = "60%" )),
                        tabPanel("Scree Plot", plotOutput("Scree_Plot", width = "50%"))
            ), br(),
            
            # Adding some text details 
            div( h4(strong('Details : ')), br(),
                 h4("Data source : "),tags$a(href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144352", "NCBI GEO Link"),
                 
                 h4("Title : "),p("Targeted RNASeq using a custom panel of metabolic genes"), 
                 h4("Experiment Summary:"),p("The purpose of this experiment was to determine the time point following activation at which greatest transcriptional differences in metabolic genes were observed between primary human CD4+ T cells treated with metformin + 2-deoxyglucose combination regimen and cells treated with vehicle."), 
                 h4("Methods:"), p(" Targeted RNA-Seq data was generated from primary human CD4+ T cell mRNA using a Qiagen custom designed panel of metabolic genes and the Illumina TruSeq stranded mRNA prep."), 
                 
                 h4("Overall design :"),p(" Targeted RNASeq transcriptomic profiling of primary human CD4+ T cells activated for 0, 3, 6, 9, and 24 hr with soluble anti-CD3 and anti-CD28 tetrameric complexes (Immunocult) in the presence of metformin and 2-deoxyglucose, of 2 individual donors."), 
                 
                 h4("Citations :"),p(" Tan SY, Kelkar Y, Hadjipanayis A, Shipstone A et al. Metformin and 2-Deoxyglucose Collaboratively Suppress Human CD4<sup>+</sup> T Cell Effector Functions and Activation-Induced Metabolic Reprogramming. J Immunol 2020 Aug 15;205(4). PMID :'32641388'"),
                 br(),br()
            )
            )
        )
    )
)



