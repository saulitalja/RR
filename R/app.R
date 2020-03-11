library(shiny)
library(shinythemes) # Theme
library(magrittr) # IS THIS NEEDED?
library(shinyhelper) # Pop up helps
library(ggplot2) # Plotting
library(mc2d) # Pert- probability distribution
library(ipc) # Asyc computing
library(promises) # Asyc computing
library(future) # Asyc computing
plan(multiprocess) # Asyc computing

source("functions.R")
#########################################################
ui <- 
  
  navbarPage("FinnSURV-Assess PWN",
             theme = shinythemes::shinytheme("united"),               
             
             tabPanel("1. Upload data and define parameter values",
                      withMathJax(),
                      tabsetPanel(
                        tabPanel(h5("1.1 Wood sampling"),
                                 fluidRow(h3("")),
                                 fluidRow(column(3,
                                                 wellPanel(h4(tags$b("The number of inspected sites")),
                                                           tags$hr(),
                                                           help_csv_N_w(),
                                                           fileInput("file_N_w","Upload a csv file", 
                                                                     multiple = FALSE,
                                                                     accept = c("text/csv",
                                                                                "text/comma-separated-values,text/plain",
                                                                                ".csv")),
                                                           helpText("Scroll down to see the uploaded table")),
                                                 
                                                 wellPanel(help_site_w(),
                                                           numericInput(inputId = "site_w",
                                                                        label = "The size of inspection site, km\\(^2\\)",
                                                                        value = 0.35, min = 0, max = 5, step = 0.1)),
                                                 wellPanel(
                                                   help_TSe_w(),
                                                   numericInput(inputId = "TSe_w",
                                                                label = "Test sensitivity",
                                                                value = 1, min = 0, max = 1, step = 0.01))),
                                          
                                          column(3,wellPanel(
                                            fluidRow(column(11,
                                                            help_n_w(),
                                                            h4(tags$b("The number of wood objects sampled per inspected site")),
                                                            tags$hr(),
                                                            radioButtons(inputId = "select_n_w_input_type", 
                                                                         label = "", 
                                                                         choices = c("Upload a csv file" = 1, 
                                                                                     "Estimate as a probability distribution" = 2), 
                                                                         selected = 2),
                                                            tags$hr(),
                                                            uiOutput(outputId = "upload_n_w_data"))))),
                                          
                                          column(3,conditionalPanel(condition = "input.select_n_w_input_type == '2'",
                                                                    wellPanel(h5(tags$i("The estimated probability distribution of the number of wood objects sampled per inspected site")),
                                                                              plotOutput(outputId = "n_wood", height=300), align = "center")))
                                 ),
                                 
                                 h5(textOutput("Table_legend_N_w")),
                                 tableOutput("table_N_w"),
                                 
                                 conditionalPanel(condition = "input.select_n_w_input_type == '1'",
                                                  h5(textOutput("Table_legend_n_w")),
                                                  tableOutput("table_n_w"))
                        ),
                        
                        tabPanel(h5("1.2", tags$i("Monochamus"), "trapping"),
                                 fluidRow(h3("")),
                                 fluidRow(column(3,
                                                 wellPanel(h4(tags$b("The number of inspected sites")),
                                                           tags$hr(),
                                                           help_csv_N_M(),
                                                           fileInput("file_N_M", "Upload a csv file",
                                                                     multiple = FALSE,
                                                                     accept = c("text/csv",
                                                                                "text/comma-separated-values,text/plain",
                                                                                ".csv")),
                                                           helpText("Scroll down to see the uploaded table")),
                                                 wellPanel(help_site_M(),
                                                           numericInput(inputId = "site_M",
                                                                        label = "The size of inspection site, km\\(^2\\)",
                                                                        value = 0.25, min = 0, max = 5, step = 0.1)),
                                                 wellPanel(help_TSe_M(),
                                                           numericInput(inputId = "TSe_M",
                                                                        label = "Test sensitivity",
                                                                        value = 1, min = 0, max = 1, step = 0.01))
                                 ),
                                 
                                 column(3,wellPanel(
                                   fluidRow(column(11,
                                                   help_n_M(),
                                                   h4(tags$b("The number of", tags$i("Monochamus"), "sampled per inspected site")),
                                                   tags$hr(),
                                                   radioButtons(inputId = "select_n_M_input_type", 
                                                                label = "", 
                                                                choices = c("Upload a csv file" = 1, "Estimate as a probability distribution" = 2), 
                                                                selected = 1),
                                                   tags$hr(),
                                                   uiOutput(outputId = "upload_n_M_data"))))),
                                 
                                 column(3,conditionalPanel(condition = "input.select_n_M_input_type == '2'",
                                                           wellPanel(h5(tags$i("The estimated probability distribution of the number of", tags$i("Monochamus"), "sampled per inspected site")),
                                                                     plotOutput(outputId = "n_Monochamus",height=300),align = "center")))
                                 ),
                                 
                                 h5(textOutput("Table_legend_N_M")),
                                 tableOutput("table_N_M"),
                                 
                                 conditionalPanel(condition = "input.select_n_M_input_type == '1'",
                                                  h5(textOutput("Table_legend_n_M")),
                                                  tableOutput("table_n_M"))
                        ),
                        
                        tabPanel(h5("1.3 Target population"),
                                 fluidRow(h3("")),
                                 fluidRow(column(6,h4("At the level of inspection sites"), align = "center"),
                                          column(3,h4("At the level of regions"), align = "center")),
                                 
                                 fluidRow(column(3,
                                                 wellPanel(fluidRow(column (11,help_d_w(),
                                                                            h4(tags$b("The density of wood objects suitable for sampling")),
                                                                            tags$hr(),
                                                                            radioButtons(inputId = "select_d_w_input_type", 
                                                                                         label = "", 
                                                                                         choices = c("Use the estimate from Hannunen and Tuomola (2020)" = 1, "Estimate as a probability distribution" = 2), 
                                                                                         selected = 1),
                                                                            tags$hr(),
                                                                            uiOutput(outputId = "upload_d_w_data"))))),
                                          
                                          column(3,       
                                                 wellPanel(fluidRow(column(11,help_d_M(),
                                                                           h4(tags$b("The density of adult", tags$i("Monochamus"), "beetles")),
                                                                           tags$hr(),
                                                                           radioButtons(inputId = "select_d_M_input_type", 
                                                                                        label = "", 
                                                                                        choices = c("Use the estimate from Hannunen and Tuomola (2020)" = 1, "Estimate as a probability distribution" = 2), 
                                                                                        selected = 1),
                                                                           tags$hr(),
                                                                           uiOutput(outputId = "upload_d_M_data"))))),
                                          
                                          
                                          column(3,     
                                                 wellPanel(help_host_area(),
                                                           h4(tags$b("The area with host plants, km", tags$sup("2"))),
                                                           hr(),
                                                           help_csv_host_area(),
                                                           fileInput("file_host_area", "Upload a csv file",
                                                                     multiple = FALSE,
                                                                     accept = c("text/csv",
                                                                                "text/comma-separated-values,text/plain",
                                                                                ".csv")),
                                                           helpText("Scroll down to see the uploaded table")))
                                 ),
                                 fluidRow(column(3,
                                                 wellPanel(
                                                   h5(tags$i("The estimated probability distribution of the density of wood objects suitable for sampling")),
                                                   plotOutput("d_wood", height=300), align = "center")),
                                          column(3,
                                                 wellPanel(
                                                   h5(tags$i("The estimated probability distribution of the density of Monochamus adults")),
                                                   plotOutput("d_Monochamus", height=300), align = "center"))),
                                 
                                 fluidRow(column(12,
                                                 h5(textOutput("Table_legend_host_area")),
                                                 tableOutput("table_host_area")))
                        ),
                        
                        tabPanel(h5("1.4 Entry sites"),
                                 fluidRow(h3("")),
                                 fluidRow(column(3,
                                                 wellPanel(help_entry_sites(),
                                                           h4(tags$b("The area of entry sites, km",tags$sup("2"))),
                                                           tags$hr(),
                                                           help_csv_entry_sites(),
                                                           fileInput("file_entry_sites", "Upload a csv file",
                                                                     multiple = FALSE,
                                                                     accept = c("text/csv",
                                                                                "text/comma-separated-values,text/plain",
                                                                                ".csv")),
                                                           helpText("Scroll down to see the uploaded table")))),
                                 fluidRow(column(12,
                                                 h5(textOutput("Table_legend_entry_sites")),
                                                 tableOutput("table_entry_sites")))
                        ))),
             
             tabPanel("2. Define the aim of the survey",  
                      
                      fluidRow(
                        column(3,
                               wellPanel(help_survey_type(),
                                         h4(tags$b("The aim of the surveys is")),
                                         tags$hr(),
                                         radioButtons(inputId = "select_survey_type", 
                                                      label = "", 
                                                      choices = c("to provide evidence to justify import requirements related to PWN and 
                                                                  to facilitate export to countries with corresponding requirements" = 1,
                                                                  "to detect possible PWN invasions early enough to enable successful eradication" = 2), 
                                                      selected = ""))),
                        column(3,  
                               uiOutput(outputId = "localDP")),
                        column(3,
                               uiOutput(outputId = "globalDP")))),
             
             tabPanel("3. View results",
                      
                      fluidRow(column(2,
                                      wellPanel(help_PriorPfree(),
                                                numericInput(inputId = "Prior_Pfree",
                                                             label = "Initial prior probability of freedom",
                                                             value = 0.5, min = 0, max = 1, step = 0.05)),
                                      wellPanel(help_Pinv(),
                                                h5(tags$b("Mean time between invasions, years")),
                                                numericInput(inputId = "Finv_min",
                                                             label = "Min",
                                                             value = 2, min = 2, max = 100, step = 1),
                                                numericInput(inputId = "Finv_max",
                                                             label = "Max",
                                                             value = 50, min = 2, max = 10000, step = 10)),
                                      wellPanel(help_n_i(),
                                                numericInput(inputId = "n_i",
                                                             label = "Number of iterations",
                                                             value = 10, min = 100, max = 10000, step = 100))),
                               
                               column(1,
                                      actionButton(inputId = "run", label = "RUN"),
                                      tags$hr(),
                                      actionButton(inputId = "cancel", label = "CANCEL")),
                               
                               column(8,
                                      tabsetPanel(
                                        tabPanel("Sensitivity - Country",
                                                 p(),
                                                 plotOutput(outputId = "SSe_c"),
                                                 tags$hr(),
                                                 helpText("The dots denote the medians and the bars the 95% confidence intervals of the assessment results")),
                                        tabPanel("Probability of freedom - Country",
                                                 p(),
                                                 plotOutput(outputId = "Pfree_c"),
                                                 tags$hr(),
                                                 helpText("The colored area shows the 95% confidence intervals of the assessment results")),
                                        tabPanel("Sensitivity - Regions",
                                                 p(),
                                                 plotOutput(outputId = "SSe_r", height = "1000px"),
                                                 tags$hr(),
                                                 helpText("The dots denote the medians and the bars the 95% confidence intervals of the assessment results")),
                                        tabPanel("Probability of freedom - Regions",
                                                 p(),
                                                 plotOutput(outputId = "Pfree_r", height = "1000px"),
                                                 tags$hr(),
                                                 helpText("The colored area shows the 95% confidence intervals of the assessment results"))
                                      )))),
             
             tabPanel("4. Download results",
                      
                      fluidRow(column(3,                                          
                                      wellPanel(h4(tags$b("The sensitivity of the annual surveys")),
                                                tags$hr(),
                                                selectInput("SSe_table_to_dl","2.5, 50 and 97.5% fractiles as separate csv files",
                                                            choices = c("2.5%",
                                                                        "50%",
                                                                        "97.5%")),
                                                downloadButton("download_SSe_tables", "Download"),
                                                tags$hr(),
                                                h5(tags$b("All the fractiles as a rds file")),
                                                downloadButton("download_SSe_fractiles", "Download"),
                                                tags$hr(),
                                                h5(tags$b("All iterations as a rds file")),
                                                downloadButton("download_SSe_iterations", "Download")
                                      )),
                               column(3,
                                      wellPanel(h4(tags$b("The probability of freedom after the last survey")),
                                                tags$hr(),
                                                selectInput("Pfree_table_to_dl","2.5, 50 and 97.5% fractiles as separate csv files",
                                                            choices = c("2.5%",
                                                                        "50%",
                                                                        "97.5%")),
                                                downloadButton("download_Pfree_tables", "Download"),
                                                tags$hr(),
                                                h5(tags$b("All the fractiles as a rds file")),
                                                downloadButton("download_Pfree_fractiles", "Download"),
                                                tags$hr(),
                                                h5(tags$b("All iterations as a rds file")),
                                                downloadButton("download_Pfree_iterations", "Download")
                                      )),
                               column(3,
                                      wellPanel(h4(tags$b("Download figures")),
                                                tags$hr(),
                                                help_dl_figs(),
                                                selectInput("fig_to_dl", "Choose a figure",
                                                            choices = c("Sensitivity - Country", 
                                                                        "Probability of freedom - Country",
                                                                        "Sensitivity - Regions", 
                                                                        "Probability of freedom - Regions")),
                                                downloadButton("download_fig", "Download"))
                               ))),
             
             tabPanel("About FinnSURV-Assess PWN",
                      fluidRow(column(5,
                                      tags$br(),
                                      wellPanel( 
                                        p(tags$b("FinnSURV-Assess PWN is a tool for assessing the confidence in freedom from pine wood nematode 
                                                 (PWN",",", tags$i("Bursaphelenchus xylophilus"),") gained in official quarantine pest surveys 
                                                 in areas where PWN is not expected to cause symptoms.")),
                                        p(),
                                        p("FinnSURV-Assess PWN can be used to assess both the sensitivity of annual surveys and the probability 
                                          of freedom achieved in multiannual surveys. It assumes that a) the surveys are composed of inspections 
                                          that cover a fixed sized area and b) in the inspections one or more wood or Monochamus samples are collected. 
                                          If all the needed data is provided separately for all regions of a country, the assessment will be done 
                                          separately for each region and for the whole country."),
                                        p(),
                                        p("FinnSURV-Assess PWN was developed in the Risk Assessment Unit of the Finnish Food Authority 
                                          by Salla Hannunen and Juha Tuomola in 2020. The methodology used in it is described in detail in 
                                          Hannunen and Tuomola (2020) and references therein, especially", 
                                          tags$a("Cannon (2002)", href = "https://doi.org/10.1016/S0167-5877(01)00262-8", target = "_blank"),",", 
                                          tags$a("Martin et al. (2007)", href = "https://doi.org/10.1016/j.prevetmed.2006.09.008", target = "_blank"),",", 
                                          tags$a("Efsa (2012)", href = "https://doi.org/10.2903/sp.efsa.2012.EN-366", target = "_blank"),"and",
                                          tags$a("Efsa (2018)", href = "https://doi.org/10.2903/sp.efsa.2018.EN-1399", target = "_blank"),
                                          ".")
                                        )
                                        ),
                               column(1),
                               
                               column(5,
                                      tabsetPanel(
                                        tabPanel("Glossary",
                                                 tags$br(),
                                                 p("Apparent prevalence = The proportion of samples testing positive"),
                                                 p(),
                                                 p("Design prevalence = Roughly, design prevalence determines the minimum prevalence that the 
                                                   survey is aimed to detect. If the pest prevalence is equal to or greater than the design prevalence, 
                                                   at least one infested individual will be detected in the survey, with the probability equal to the 
                                                   sensitivity of the survey."),
                                                 p(),
                                                 p("Early detection survey = A survey that aims to detect possible PWN invasions early enough to 
                                                   enable successful eradication"),
                                                 p(),
                                                 p("Entry site = A site where the probability of PWN introduction is elevated, i.e. harbors, 
                                                   industrial areas and landfills"),
                                                 p(),
                                                 p("Import-export survey = A survey that aims to provide evidence to justify import requirements 
                                                   related to PWN and to facilitate export to countries with corresponding requirements"),
                                                 p(),
                                                 p("Initial prior probability of freedom = The probability that the prevalence of the pest is 
                                                   below the design prevalence before the first survey"),
                                                 p(),
                                                 p("Probability of freedom = The probability that the prevalence of the pest is below the design 
                                                   prevalence if the pest is not detected in the surveys"),
                                                 p(),
                                                 p("Sensitivity = Roughly, sensitivity determines the probability with which a survey is expected 
                                                   to succeed in its aim. If the pest prevalence is equal to or greater than the design prevalence, 
                                                   at least one infested individual will be detected in the survey, with the probability equal to 
                                                   the sensitivity of the survey."),
                                                 p(),
                                                 p("Target population at the level of inspection site = Wood objects suitable for sampling / Monochamus adults"),
                                                 p(),
                                                 p("Target population at the level of regions = The area with PWN host plants in the 
                                                   area for which the results of the survey will be generalized"),
                                                 p(),
                                                 p("Test sensitivity = The probability that the pest is detected in the laboratory analysis, given that it 
                                                   was present in the wood object(s) / Monochamus beetle(s) included in the sample")
                                                 ),
                                        tabPanel("References",
                                                 tags$br(),
                                                 p(tags$a("Cannon RM (2002) Demonstrating disease freedom - Combining confidence levels. 
                                                          Preventive Veterinary Medicine 52: 227-249.", 
                                                          href = "https://doi.org/10.1016/S0167-5877(01)00262-8", target="_blank")),
                                                 p("Corine..."),
                                                 p(tags$a("European Food Safety Authority (2012) 
                                                          A framework to substantiate absence of disease: the risk based estimate 
                                                          of system sensitivity tool (RiBESS) using data collated according to the EFSA Standard Sample Description - 
                                                          An example on", tags$i("Echinococcus multilocularis"),". Supporting Publications 2012 9(12):EN-366: 1-44.",
                                                          href = "https://doi.org/10.2903/sp.efsa.2012.EN-366", target="_blank")),
                                                 p(tags$a("European Food Safety Authority, Ciubotaru RM, Cortinas Abrahantes J, Oyedele J, 
                                                          Parnell S, Schrader G, Zancanaro G, Vos S (2018) Work-plan and methodology for EFSA 
                                                          to develop plant pest survey guidelines for EU Member States. EFSA supporting publication 2018 15(3):EN-1399: 1-36",
                                                          href = "https://doi.org/10.2903/sp.efsa.2018.EN-1399", target="_blank")),
                                                 p(tags$a("European Union (2012) Commission Implementing Decision of 26 September 2012 on 
                                                          emergency measures to prevent the spread within the Union of", 
                                                          tags$i("Bursaphelenchus xylophilus"),"(Steiner et Buhrer) Nickle et al. (the pine wood nematode) 
                                                          (notified under document C(2012) 6543) (2012/535/EU). 
                                                          Official Journal of the European Union L 266 2.10.2012: 42-52.",
                                                          href = "https://eur-lex.europa.eu/legal-content/EN/TXT/?uri=CELEX:02012D0535-20170310", target="_blank")),
                                                 p("Hannunen S and Tuomola J (2020) 
                                                   Assessing the probability of freedom 
                                                   from pine wood nematode based on 19 years of surveys. NeoBiota..."),
                                                 p(tags$a("Martin PAJ, Cameron AR, Greiner M (2007) 
                                                          Demonstrating freedom from disease using multiple complex data sources 
                                                          1: A new methodology based on scenario trees. 
                                                          Preventive Veterinary Medicine 79: 71-97.", 
                                              href = "https://doi.org/10.1016/j.prevetmed.2006.09.008", target="_blank"))
                            ),
                            tabPanel("Copyright",
                                     tags$br(),
                                     p("Please refer to FinnSURV-Assess PWN as:"),
                                     p("Hannunen S and Tuomola J 2020. FinnSURV-Assess PWN - A tool for assessing 
                                       the confidence of freedom from the pine wood nematode gained in official surveys. 
                                       Finnish Food Authority, Helsinki, Finland. Available at..."),
                                     tags$br(),
                                     p("FinnSURV-Assess PWN is published under", tags$a(href="https://creativecommons.org/licenses/by-nc-sa/4.0/","the Creative Commons Attribution-NonCommercial-ShareAlike 4.0
                                                                                        International", target="_blank"), "license.")
                            )
                          )
                          )))
                 
)

#####################################################

server <- function(input,output,session){

  observe_helpers()
  
  # Stop session when the window is closed
  session$onSessionEnded(stopApp)
  
  ##########
  # INTERACTIVE PARTS OF THE UI (RENDER UI)
  
  # Inspection site level design prevalences of the different survey types
  output$localDP <- renderUI({
    req(input$select_survey_type)
    switch(input$select_survey_type,
           "1" = wellPanel(help_localDP_ie(),
                           h4(tags$b("Inspection site level design prevalence")),
                           tags$hr(),
                           numericInput(inputId = "DP_w",
                                        label = "Wood",
                                        value = 0.12, min = 0, max = 0.12, step = 0.01),
                           numericInput(inputId = "DP_M",
                                        label = tags$i("Monochamus"),
                                        value = 0.09, min = 0, max = 0.09, step = 0.01)),
           
           "2" = wellPanel(help_localDP_ed(),
                           h4(tags$b("Inspection site level design prevalence")),
                           tags$hr(),
                           numericInput(inputId = "DP_w",
                                        label = "Wood",
                                        value = 0.06, min = 0, max = 0.12, step = 0.01),
                           numericInput(inputId = "DP_M",
                                        label = tags$i("Monochamus"),
                                        value = 0.045, min = 0, max = 0.09, step = 0.01)))
  })
  
  # Region/country level design prevalences of the different survey types 
  output$globalDP <- renderUI({
    req(input$select_survey_type)
    switch(input$select_survey_type,
           "1" = wellPanel(help_globalDP(),
                           h4(tags$b("Country level design prevalence")),
                           tags$hr(),
                           numericInput(inputId = "DPr",
                                        label = "",
                                        value = 0.01, min = 0, max = 1, step = 0.001)),
           "2" = wellPanel(help_max_inf_size(),
                           h4(tags$b("Maximum acceptable size of PWN infestation at detection, km",tags$sup("2"))),
                           tags$hr(),
                           numericInput(inputId = "max_inf_size",
                                        label = "",
                                        value = 314, min = 1, max = 5000, step = 10)))
  })
  
  # Number of wood objects sampled per inspected site: upload data or define as probability distribution
  output$upload_n_w_data <- renderUI({
    
    switch(input$select_n_w_input_type,
           "1" =  fluidRow(column(12,
                                  help_csv_n_w(),
                                  fileInput("file_n_w", "",
                                            multiple = FALSE,
                                            accept = c("text/csv",
                                                       "text/comma-separated-values,text/plain",
                                                       ".csv")),
                                  helpText("Scroll down to see the uploaded table"))),
           
           "2" =  fluidRow(
             help_pert(),            
             column(6,
                    numericInput(inputId = "n_w_min",
                                 label = "Min",
                                 value = 1, min = 1, max = 1000, step = 1)),
             column(6,
                    numericInput(inputId = "n_w_max",
                                 label = "Max",
                                 value = 5, min = 1, max = 1000, step = 1)),
             column(6,
                    numericInput(inputId = "n_w_mode",
                                 label = "Mode",
                                 value = 2, min = 1, max = 1000, step = 1)),
             column(6,
                    numericInput(inputId = "n_w_lambda",
                                 label = "Lambda",
                                 value = 1, min = 1, max =20, step = 1)),
             column(12,helpText("See the plot on the right for the distribution")))
    )})
  
  # Number of Monochamus sampled per inspected site: upload data or define as probability distribution 
  output$upload_n_M_data <- renderUI({
    
    switch(input$select_n_M_input_type,
           "1"=    fluidRow(column(12,
                                   help_csv_n_M(),                                   
                                   fileInput("file_n_M", "",
                                             multiple = FALSE,
                                             accept = c("text/csv",
                                                        "text/comma-separated-values,text/plain",
                                                        ".csv")),
                                   helpText("Scroll down to see the uploaded table"))),
           
           "2" =  fluidRow(
             help_pert(),
             column(6,
                    numericInput(inputId = "n_M_min",
                                 label = "Min",
                                 value = "", min = 1, max = 1000, step = 1)),
             column(6,
                    numericInput(inputId = "n_M_max",
                                 label = "Max",
                                 value = "", min = 1, max = 1000, step = 1)),
             column(6,
                    numericInput(inputId = "n_M_mode",
                                 label = "Mode",
                                 value = "", min = 1, max = 1000, step = 1)),
             column(6,
                    numericInput(inputId = "n_M_lambda",
                                 label = "Lambda",
                                 value = "", min = 1, max =20, step = 1)),
             column(12,helpText("See the plot on the right for the distribution")))
    )})
  
  # Density of wood objects suitable for sampling: upload data or define as porbability distribution  
  output$upload_d_w_data <- renderUI({
    
    switch(input$select_d_w_input_type,
           "1" = fluidRow(column(12,
                                 helpText("See the plot below for the distribution"))),
           
           "2" = fluidRow(
             column(11,h5(tags$b("Number per km", tags$sup("2")))),
             help_pert(),
             column(6,
                    numericInput(inputId = "d_w_min",
                                 label = "Min",
                                 value = 1, min = 1, max = 10000, step = 1)),
             column(6,
                    numericInput(inputId = "d_w_max",
                                 label = "Max",
                                 value = 60, min = 1, max = 10000, step = 1)),
             column(6,
                    numericInput(inputId = "d_w_mode",
                                 label = "Mode",
                                 value = 30, min = 1, max = 10000, step = 1)),
             column(6,
                    numericInput(inputId = "d_w_lambda",
                                 label = "Lambda",
                                 value = 1, min = 1, max = 10, step = 1)),
             column(12,helpText("See the plot below for the distribution")))
    )})
  
  # Density of Monochamus adults: upload data or define as probability distribution  
  output$upload_d_M_data <- renderUI({
    
    switch(input$select_d_M_input_type,
           "1"= fluidRow(column(12,
                                helpText("See the plot below for the distribution"))),
           
           "2" = fluidRow(column(11,h5(tags$b("Number per km", tags$sup("2")))),
                          help_pert(),
                          column(6,
                                 numericInput(inputId = "d_M_min",
                                              label = "Min",
                                              value = 1, min = 1, max = 10000, step = 1)),
                          column(6,
                                 numericInput(inputId = "d_M_max",
                                              label = "Max",
                                              value = 60, min = 1, max = 10000, step = 1)),
                          column(6,
                                 numericInput(inputId = "d_M_mode",
                                              label = "Mode",
                                              value = 30, min = 1, max = 10000, step = 1)),
                          column(6,
                                 numericInput(inputId = "d_M_lambda",
                                              label = "Lambda",
                                              value = 1, min = 1, max = 10, step = 1)),
                          column(12,helpText("See the plot below for the distribution")))
    )})
  
  ##########  
  # UPLOADING DATA AND PRINTING IT AS TABLES  
  
  # The number of inspected sites in the wood sampling component of the survey
  df_N_w <- reactive({
    req(input$file_N_w)
    A <- read.csv(input$file_N_w$datapath,
                  header = TRUE,
                  sep = ",",
                  quote = "")
    A})
  
  output$Table_legend_N_w <- eventReactive(input$file_N_w,{
    "The number of inspected sites"})
  
  output$table_N_w <- renderTable({
    req(input$file_N_w)
    df_N_w()})
  
  # The number of wood objects sampled per inspected site
  df_n_w <- reactive({
    req(input$file_n_w)
    A <- read.csv(input$file_n_w$datapath,
                  header = TRUE,
                  sep = ",",
                  quote = "")
    A})
  
  output$Table_legend_n_w <- eventReactive(input$file_n_w,{
    "The mean number of wood objects sampled per inspected site"})
  
  output$table_n_w <- renderTable({
    req(input$file_n_w)
    df_n_w()})
  
  # The number of inspected sites in the Monochamsu trapping component of the survey 
  df2 <- reactive({
    req(input$file_N_M)
    A <- read.csv(input$file_N_M$datapath,
                  header = TRUE,
                  sep = ",",
                  quote = "")
    A})
  
  output$Table_legend_N_M <- eventReactive(input$file_N_M,{
    "The number of inspected sites"})  
  
  output$table_N_M <- renderTable({
    req(input$file_N_M)
    df2()})
  
  # The number of Monochamus sampled per inspected site
  df2.1 <- reactive({
    req(input$file_n_M)
    A <- read.csv(input$file_n_M$datapath,
                  header = TRUE,
                  sep = ",",
                  quote = "")
    A})
  
  output$Table_legend_n_M <- eventReactive(input$file_n_M,{
    "The mean number of beetles sampled per inspected site"})
  
  output$table_n_M <- renderTable({
    req(input$file_n_M)
    df2.1()})
  
  # The area with host plants
  df_host_area <- reactive({
    req(input$file_host_area)
    A <- read.csv(input$file_host_area$datapath,
                  header = TRUE,
                  sep = ",",
                  quote = "")
    A})
  
  output$Table_legend_host_area <- eventReactive(input$file_host_area,{
    "The area with host plants, km2"})
  
  output$table_host_area <- renderTable({
    req(input$file_host_area)
    df_host_area()})
  
  # The area of entry sites
  df_entry_sites <- reactive({
    req(input$file_entry_sites)
    A <- read.csv(input$file_entry_sites$datapath,
                  header = TRUE,
                  sep = ",",
                  quote = "")
    A})
  
  output$Table_legend_entry_sites <- eventReactive(input$file_entry_sites,{
    "The area of entry sites, km2"})
  
  output$table_entry_sites <- renderTable({
    req(input$file_entry_sites)
    df_entry_sites()})
  
  ##########  
  # VECTORS NEEDED IN THE CALCULATIONS AND FIGURES
  
  # The considered mean times between invasions
  Finv_FI <- eventReactive(input$run,{
    seq(input$Finv_min,input$Finv_max,1)  
  })
  
  # The years when surveys were conducted
  Y <- eventReactive(input$run,{
    a = data.matrix(df_N_w())
    seq(a[1,1],a[nrow(a),1],1)
  })
  
  # The number of years when surveys were conducted
  n_y <- eventReactive(input$run,{
    length(Y())
  })
  
  # The names of the regions included in the survey
  region_names <- eventReactive(input$run,{
    A <- read.csv(input$file_N_w$datapath,
                  header = TRUE,
                  sep = ",",
                  quote = "'")
    colnames(A[2:ncol(A)])
  })
  
  # The number of regions included in the survey + 1 (to include the whole country)
  n_r <- eventReactive(input$run,{
    ncol(df_N_w())
  })
  
  ##########
  # UPLOADED DATA FOR THE CALCULATIONS
  
  # The number of inspected sites in the wood sampling component of the survey
  N_wood <- reactive({
    a = data.matrix(df_N_w())
    b = a[,2:ncol(a)]
    array(b, dim=c(nrow(b),ncol(b),input$n_i))
  })
  
  # The number of inspected sites in the Monochamus trapping component of the survey
  N_Monochamus <- reactive({
    a = data.matrix(df2())
    b = a[,2:ncol(a)]
    array(b, dim=c(nrow(b),ncol(b),input$n_i))
  })
  
  # The number of wood objects sampled per inspected site
  n_wood <- reactive({
    if(input$select_n_w_input_type == 1){
      a = data.matrix(df_n_w())
      b = a[,2:ncol(a)]
      A <- array(b, dim=c(nrow(b),ncol(b),input$n_i))  
    }else{
      req(input$n_w_min<=input$n_w_mode,input$n_w_min<=input$n_w_max,input$n_w_mode<=input$n_w_max)
      A <- round(rpert(input$n_i,input$n_w_min,input$n_w_mode,input$n_w_max,input$n_w_lambda))
    }
    A
  }) 
  
  # The number of Monochamus adults sampled per inspected site
  n_Monochamus <- reactive({
    if(input$select_n_M_input_type == 1){
      a = data.matrix(df2.1())
      b = a[,2:ncol(a)]
      A <- array(b, dim=c(nrow(b),ncol(b),input$n_i))  
    }else{
      req(input$n_M_min<=input$n_M_mode,input$n_M_min<=input$n_M_max,input$n_M_mode<=input$n_M_max)
      A<- round(rpert(input$n_i,input$n_M_min,input$n_M_mode,input$n_M_max,input$n_M_lambda))
    }
    A
  }) 
  
  # The density of wood objects suitable for sampling
  d_wood <- reactive({
    
    if(input$select_d_w_input_type == 1){
      
      # The number of Monochamus-suitable dead wood objects per km2 
      obj <- round(rpert(input$n_i,166,288,398,1))
      
      # The proportion of Monochamus-suitable dead wood objects that is suitable for sampling 
      psam <- runif(input$n_i,0.05,0.95)
      
      # The density of wood objects suitable for sampling, number/km2 
      round(obj*psam)
      
    }else{
      req(input$d_w_min<=input$d_w_mode,input$d_w_min<=input$d_w_max,input$d_w_mode<=input$d_w_max)
      round(rpert(input$n_i,input$d_w_min,input$d_w_mode,input$d_w_max,input$d_w_lambda))
    }
  }) 
  
  # The number of wood objects suitable for sampling per inspection site
  p_wood <- reactive({
    b <- array(d_wood(),dim=c(n_y(),n_r()-1,input$n_i))
    c = round(b*input$site_w)
    c[c < ceiling(1/input$DP_w)] <- ceiling(1/input$DP_w)
    c
  })
  
  # The density of Monochamus adults
  d_Monochamus <- reactive({
    
    if(input$select_d_M_input_type == 1){
      
      # The number of dead wood objects occupied by Monochamus per km2
      obju <- round(rpert(input$n_i,13.28,28.8,47.76,1))
      
      # The number of Monochamus eggs laid per Monochamus-suitable dead-wood object 
      fobj <- rpert(input$n_i,6,31,88,1)
      
      # The proportion of Monochamus surviving from egg to egg-laying adults
      surv <- rpert(input$n_i,0.1,0.25,0.4,1)
      
      # The density of Monochamus adults, number/km2
      obju*fobj*surv
      
    }else{
      req(input$d_M_min<=input$d_M_mode,input$d_M_min<=input$d_M_max,input$d_M_mode<=input$d_M_max)
      round(rpert(input$n_i,input$d_M_min,input$d_M_mode,input$d_M_max,input$d_M_lambda))
    }
  }) 
  
  # The number of Monochamus adults per inspection site
  p_Monochamus <- reactive({
    b <- array(d_Monochamus(),dim=c(n_y(),n_r()-1,input$n_i))
    c = round(b*input$site_M)
    c[c < ceiling(1/input$DP_M)] <- ceiling(1/input$DP_M)
    c
  })
  
  # The area with host plants
  host_area <- reactive({
    a = data.matrix(df_host_area())
    b = a[,2:ncol(a)]
    array(b, dim=c(nrow(b),ncol(b),input$n_i))
  })
  
  # The area of entry sites
  entry_sites <- reactive({
    a = data.matrix(df_entry_sites())
    b = a[,2:ncol(a)]
    matrix(b, nrow(b), ncol(b))
  })
  
  # The relative probability of invasion to the regions
  RP <- reactive({
    a = entry_sites()/rowSums(entry_sites())
    array(a, dim=c(nrow(a),ncol(a),input$n_i))
  })
  
  # The effective probability of infestation for the import-export survey
  # and the regional-level design prevalence for the early detection survey
  DPr_adj <- reactive({
    
    # The proportion of the target population in the regions
    PropPop = host_area()/rowSums(host_area()[,,1])
    
    # Weighted probability of invasion in the regions
    WR = RP()/rowSums(PropPop[,,1]*RP()[,,1])
    
    # Effective probability of infestation for the import-export survey
    if(input$select_survey_type == 1){
      A = input$DPr*WR
      # Regional-level design prevalence for the early detection survey
    }else{
      A = input$max_inf_size/host_area()
    }
    A
  })
  
  ##########  
  # FIGURES FOR THE TABS "UPLOAD DATA"
  
  # The number of wood objects sampled per inspected site
  output$n_wood <- renderPlot({
    
    validate(
      need(all(input$n_w_min<=input$n_w_mode,
               input$n_w_min<=input$n_w_max,
               input$n_w_mode<=input$n_w_max,
               input$n_w_min >= 0,
               input$n_w_mode >= 0,
               input$n_w_max >= 0,
               input$n_w_lambda >= 0),
           "Make sure that min, mode, max and lambda are logical"),
      need(any(input$n_w_min!=input$n_w_mode,
               input$n_w_min!=input$n_w_max,
               input$n_w_mode!=input$n_w_max),
           "Cannot present a probability distribution if min, mode and max are the same, but analysis can be run if min, max and lambda are all defined")
    )
    
    x <- seq(input$n_w_min,input$n_w_max,0.1)
    dist = dpert(x,input$n_w_min,input$n_w_mode,input$n_w_max,input$n_w_lambda)
    plot(x,dist,type = 'l',
         xlab = "Number per inspection", ylab = "Probability")
  })
  
  # the number of Monochamus sampled per inspected site
  output$n_Monochamus <- renderPlot({
    
    validate(
      need(all(input$n_M_min<=input$n_M_mode,
               input$n_M_min<=input$n_M_max,
               input$n_M_mode<=input$n_M_max,
               input$n_M_min >= 0,
               input$n_M_mode >= 0,
               input$n_M_max >= 0,
               input$n_M_lambda >= 0),
           "Make sure that min, mode, max and lambda are logical"),
      need(any(input$n_M_min!=input$n_M_mode,
               input$n_M_min!=input$n_M_max,
               input$n_M_mode!=input$n_M_max),
           "Cannot present a probability distribution if min, mode and max are the same, but analysis can be run if min, max and lambda are all defined")
    )
    
    x <- seq(input$n_M_min,input$n_M_max,0.1)
    dist = dpert(x,input$n_M_min,input$n_M_mode,input$n_M_max,input$n_M_lambda)
    plot(x,dist,type = 'l', 
         xlab = "Number per inspection", ylab = "Probability")
  })
  
  # The density of wood objects suitable for sampling
  output$d_wood <- renderPlot({
    if(input$select_d_w_input_type == 2){
      
      validate(
        need(all(input$d_w_min<=input$d_w_mode,
                 input$d_w_min<=input$d_w_max,
                 input$d_w_mode<=input$d_w_max,
                 input$d_w_min >= 0,
                 input$d_w_mode >= 0,
                 input$d_w_max >= 0,
                 input$d_w_lambda >= 0),
             "Make sure that min, mode, max and lambda are logical"),
        need(any(input$d_w_min!=input$d_w_mode,
                 input$d_w_min!=input$d_w_max,
                 input$d_w_mode!=input$d_w_max),
             "Cannot present a probability distribution if min, mode and max are the same, but analysis can be run if min, max and lambda are all defined")
      )
      
      x <- seq(input$d_w_min,input$d_w_max,0.1)
      dist = dpert(x,input$d_w_min,input$d_w_mode,input$d_w_max,input$d_w_lambda)
      plot(x,dist,type = 'l',
           xlab = expression(paste("Number per" ~ km^{2})), 
           ylab = "Probability")}
    else{
      dist = round(rpert(600000,166,288,398,1))*runif(600000,0.05,0.95)
      res <- hist(dist, freq = FALSE, breaks = 80, xlim=c(0,400))
      a = res$breaks+res$breaks[2]/2
      x = a[1:length(a)-1]
      plot(x,res$density,'l',
           xlab = expression(paste("Number per" ~ km^{2})),
           ylab = "Probability", main = "")
    }
  })
  
  # The density of Monochamus adults
  output$d_Monochamus <- renderPlot({
    if(input$select_d_M_input_type == 2){
      
      validate(
        need(all(input$d_M_min<=input$d_M_mode,
                 input$d_M_min<=input$d_M_max,
                 input$d_M_mode<=input$d_M_max,
                 input$d_M_min >= 0,
                 input$d_M_mode >= 0,
                 input$d_M_max >= 0,
                 input$d_M_lambda >= 0),
             "Make sure that min, mode, max and lambda are logical"),
        need(any(input$d_M_min!=input$d_M_mode,
                 input$d_M_min!=input$d_M_max,
                 input$d_M_mode!=input$d_M_max),
             "Cannot present a probability distribution if min, mode and max are the same, but analysis can be run if min, max and lambda are all defined")
      )
      
      x <- seq(input$d_M_min,input$d_M_max,0.1)
      dist = dpert(x,input$d_M_min,input$d_M_mode,input$d_M_max,input$d_M_lambda)
      plot(x,dist,type = 'l', 
           xlab = expression(paste("Number per" ~ km^{2})),
           ylab = "Probability")}
    else{
      dist = round(rpert(100000,13.28,28.8,47.76,1))*rpert(100000,6,31,88,1)*rpert(100000,0.1,0.25,0.4,1)
      res <- hist(dist, freq = FALSE, breaks = 50, xlim=c(0,1600))
      a = res$breaks+res$breaks[2]/2
      x = a[1:length(a)-1]
      plot(x,res$density,'l',
           xlab = expression(paste("Number per" ~ km^{2})), 
           ylab = "Probability", main = "")
    }
  })
  
  ##########
  # CALCULATING THE SENSITIVITY OF THE ANNUAL SURVEYS
  
  # All iterations
  SSe_iterations <- eventReactive(input$run, {
    
    Sensitivity(N_wood(), n_wood(), p_wood(), input$TSe_w,
                N_Monochamus(), n_Monochamus(), p_Monochamus(), input$TSe_M,
                host_area(), RP(), DPr_adj(),
                input$select_survey_type, input$DP_w, input$DP_M, input$DPr, input$max_inf_size,
                n_r(), n_y(), input$n_i)
  })
  
  # 2.5, 50,97.5% fractiles
  SSe <- eventReactive(input$run, {
    Sensitivity_fractiles(SSe_iterations(),n_r(),n_y())
  })
  
  #####################
  # CALCULATING THE PROBABILITY OF FREEDOM AFTER THE LAST SURVEY (iterations as an Async operation)
  
  nclicks <- reactiveVal(0)
  Pfree_iterations_val <- reactiveVal()
  inter <- AsyncInterruptor$new()
  
  observeEvent(input$run,{
    
    req(N_wood(), n_wood(), p_wood(), input$TSe_w,
        N_Monochamus(),n_Monochamus(),p_Monochamus(),input$TSe_M,
        host_area(), entry_sites(),
        input$select_survey_type, input$DP_w, input$DP_M, any(input$DPr!="", input$max_inf_size!=""),
        input$Prior_Pfree, Finv_FI(), input$n_i)
    
    # Start progress monitor
    if(nclicks() == 0){
      progress <- AsyncProgress$new(session, min=1, max=100, message = "Computing")
      progress$set(message = "Computing", value = 1)
    }

    # If "Run" is clicked several times, don't do anything as analysis is already running
    if(nclicks() != 0){
      showNotification("Already running analysis")
      return(NULL)
    }
    # Increment clicks and prevent concurrent analyses
    nclicks(nclicks() + 1)
    
    Pfree_iterations_val(NULL)
    
    # Create variables that can be imported to the async operation (since reactive values and input$xxs cannot)
    N_wood <- N_wood()
    n_wood <- n_wood()
    p_wood <- p_wood()
    TSe_w = input$TSe_w
    N_Monochamus <- N_Monochamus()
    n_Monochamus <- n_Monochamus()
    p_Monochamus <- p_Monochamus()
    TSe_M = input$TSe_M
    host_area <- host_area()
    RP <- RP()
    DPr_adj <- DPr_adj()
    n_r <- n_r()
    n_y <- n_y()
    Pinv_FI <- 1/Finv_FI()
    select_survey_type = input$select_survey_type
    DP_w = input$DP_w
    DP_M = input$DP_M
    DPr = input$DPr
    max_inf_size = input$max_inf_size
    Prior_Pfree = input$Prior_Pfree
    n_i = input$n_i
    
    # Compute the iterations as an async operation
    Pfree_iterations <- future({
      
      Probability_of_freedom_iterations(progressMonitor1 = function(h) inter$execInterrupts(),
                                        progressMonitor2 = function(h) progress$set(),
                                        N_wood, n_wood, p_wood,TSe_w,
                                        N_Monochamus, n_Monochamus, p_Monochamus,TSe_M,
                                        host_area, RP, DPr_adj,
                                        select_survey_type, DP_w, DP_M, DPr, max_inf_size,
                                        Pinv_FI, Prior_Pfree,
                                        n_r, n_y, n_i)
      
    }) %...>% Pfree_iterations_val()
    
    # After the promise has been evaluated set nclicks to 0 to allow for anlother run
    # and close the progress monitor
    Pfree_iterations <- finally(Pfree_iterations,
                                function(){
                                  nclicks(0)
                                  progress$close()
                                })
    
    # Return something other than the promise so shiny remains responsive
    NULL
  })
  
  # If cancell is clicked, interrupt computing and show notification
  observeEvent(input$cancel,{
    inter$interrupt("Error: Stop Future")
    showNotification("Analysis cancelled")
  })
  
  # Get the 2.5, 50 and 97.5% fractiles
  Pfree <- reactive({
    req(Pfree_iterations_val())
    Pinv_FI <- 1/Finv_FI()
    Probability_of_freedom_fractiles(Pfree_iterations_val(),Pinv_FI,n_r(),n_y())
  })
  
  #############
  # FIGURES FOR THE TAB "VIEW RESULTS"
  
  # SSe - Country
  output$SSe_c <- renderPlot({
    
    # Validate and if not valid show an error message
    # (need(all(input$file_n_w, any(input$1, input$2,input$3)),"Error message") did not work, so I had to use if...else...)
    if(input$select_n_w_input_type == 1 & input$select_n_M_input_type == 1){
      validate(
        need(input$file_N_w, "Upload data on the number of inspected sites in the wood sampling component"),
        need(input$file_n_w,"Provide information on the number of wood objects sampled per inspection site"),
        need(input$site_w, "Define the size of inspection site for wood sampling"),
        need(input$TSe_w, "Define test sensitivity for wood sampling"),
        need(input$file_N_M, "Upload data on the number of inspected sites in the Monochamus trapping component"),
        need(input$file_n_M,"Provide information on the number of Monochamus sampled per inspection site"),
        need(input$site_M, "Define the size of inspection site for Monochamus trapping"),
        need(input$TSe_M,"Define test sensitivity for Monochamus trapping"),
        need(any(input$select_d_w_input_type==1,all(input$d_w_min,input$d_w_mode,input$d_w_max,input$d_w_lambda)),"Porvide information on the density of wood objects suitable for sampling"),
        need(any(input$select_d_M_input_type==1,all(input$d_M_min,input$d_M_mode,input$d_M_max,input$d_M_lambda)),"Porvide information on the density of Monochamus adults"),
        need(all(input$d_w_min<=input$d_w_mode,input$d_w_min<=input$d_w_max,input$d_w_mode<=input$d_w_max),"Check the information given on the density of wood objects suitable for sampling"),
        need(all(input$d_M_min<=input$d_M_mode,input$d_M_min<=input$d_M_max,input$d_M_mode<=input$d_M_max),"Check the information given on the density of Monochamus adults"),
        need(input$file_host_area,"Upload data on the area with host plants"),
        need(input$file_entry_sites,"Upload data on the area of entry sites"),
        need(input$select_survey_type, "Select survey type"),
        # TAALLA
        #need(all(input$DP_w, input$DP_M, any(input$DPr, input$max_inf_size)),"Define design prevalences"),
        need(input$Prior_Pfree, "Determine the prior porbability of freedom"),
        need(all(input$Finv_min,input$Finv_max), "Determine the range of mean time between invasions to be considered"),
        need(input$n_i,"Determine the number of iterations")
      )
    }else if(input$select_n_w_input_type == 1 & input$select_n_M_input_type == 2){
      validate(
        need(input$file_N_w, "Upload data on the number of inspected sites in the wood sampling component"),
        need(input$file_n_w,"Provide information on the number of wood objects sampled per inspection site"),
        need(input$site_w, "Define the size of inspection site for wood sampling"),
        need(input$TSe_w, "Define test sensitivity for wood sampling"),
        need(input$file_N_M, "Upload data on the number of inspected sites in the Monochamus trapping component"),
        need(all(input$n_M_min,input$n_M_mode,input$n_M_max,input$n_M_lambda), "Provide information on the number of Monochamus sampled per inspection site"),
        need(all(input$n_M_min<=input$n_M_mode,input$n_M_min<=input$n_M_max,input$n_M_mode<=input$n_M_max),"Check the information given on the number of Monochamus sampled per inspection site"),
        need(input$site_M, "Define the size of inspection site for Monochamus trapping"),
        need(input$TSe_M,"Define test sensitivity for Monochamus trapping"),
        need(any(input$select_d_w_input_type==1,all(input$d_w_min,input$d_w_mode,input$d_w_max,input$d_w_lambda)),"Porvide information on the density of wood objects suitable for sampling"),
        need(any(input$select_d_M_input_type==1,all(input$d_M_min,input$d_M_mode,input$d_M_max,input$d_M_lambda)),"Porvide information on the density of Monochamus adults"),
        need(all(input$d_w_min<=input$d_w_mode,input$d_w_min<=input$d_w_max,input$d_w_mode<=input$d_w_max),"Check the information given on the density of wood objects suitable for sampling"),
        need(all(input$d_M_min<=input$d_M_mode,input$d_M_min<=input$d_M_max,input$d_M_mode<=input$d_M_max),"Check the information given on the density of Monochamus adults"),
        need(input$file_host_area,"Upload data on the area with host plants"),
        need(input$file_entry_sites,"Upload data on the area of entry sites"),
        need(input$select_survey_type, "Select survey type"),
        # TAALLA
        #need(all(input$DP_w, input$DP_M, any(input$DPr, input$max_inf_size)),"Define design prevalences"),
        need(input$Prior_Pfree, "Determine the prior porbability of freedom"),
        need(all(input$Finv_min,input$Finv_max), "Determine the range of mean time between invasions to be considered"),
        need(input$n_i,"Determine the number of iterations")
      )
    }else if(input$select_n_w_input_type == 2 & input$select_n_M_input_type == 1){
      validate(
        need(input$file_N_w, "Upload data on the number of inspected sites in the wood sampling component"),
        need(all(input$n_w_min,input$n_w_mode,input$n_w_max,input$n_w_lambda),"Provide information on the number of wood objects sampled per inspection site"),
        need(all(input$n_w_min<=input$n_w_mode,input$n_w_min<=input$n_w_max,input$n_w_mode<=input$n_w_max),"Check the information given on the number of wood objects sampled per inspection site"),
        need(input$site_w, "Define the size of inspection site for wood sampling"),
        need(input$TSe_w, "Define test sensitivity for wood sampling"),
        need(input$file_N_M, "Upload data on the number of inspected sites in the Monochamus trapping component"),
        need(input$file_n_M,"Provide information on the number of Monochamus sampled per inspection site"),
        need(input$site_M, "Define the size of inspection site for Monochamus trapping"),
        need(input$TSe_M,"Define test sensitivity for Monochamus trapping"),
        need(any(input$select_d_w_input_type==1,all(input$d_w_min,input$d_w_mode,input$d_w_max,input$d_w_lambda)),"Porvide information on the density of wood objects suitable for sampling"),
        need(any(input$select_d_M_input_type==1,all(input$d_M_min,input$d_M_mode,input$d_M_max,input$d_M_lambda)),"Porvide information on the density of Monochamus adults"),
        need(all(input$d_w_min<=input$d_w_mode,input$d_w_min<=input$d_w_max,input$d_w_mode<=input$d_w_max),"Check the information given on the density of wood objects suitable for sampling"),
        need(all(input$d_M_min<=input$d_M_mode,input$d_M_min<=input$d_M_max,input$d_M_mode<=input$d_M_max),"Check the information given on the density of Monochamus adults"),
        need(input$file_host_area,"Upload data on the area with host plants"),
        need(input$file_entry_sites,"Upload data on the area of entry sites"),
        need(input$select_survey_type, "Select survey type"),
        # TAALLA
        #need(all(input$DP_w, input$DP_M, any(input$DPr, input$max_inf_size)),"Define design prevalences"),
        need(input$Prior_Pfree, "Determine the prior porbability of freedom"),
        need(all(input$Finv_min,input$Finv_max), "Determine the range of mean time between invasions to be considered"),
        need(input$n_i,"Determine the number of iterations") 
      )
    }else if(input$select_n_w_input_type == 2 & input$select_n_M_input_type == 2){
      validate(
        need(input$file_N_w, "Upload data on the number of inspected sites in the wood sampling component"),
        need(all(input$n_w_min,input$n_w_mode,input$n_w_max,input$n_w_lambda),"Provide information on the number of wood objects sampled per inspection site"),
        need(all(input$n_w_min<=input$n_w_mode,input$n_w_min<=input$n_w_max,input$n_w_mode<=input$n_w_max),"Check the information given on the number of wood objects sampled per inspection site"),
        need(input$site_w, "Define the size of inspection site for wood sampling"),
        need(input$TSe_w, "Define test sensitivity for wood sampling"),
        need(input$file_N_M, "Upload data on the number of inspected sites in the Monochamus trapping component"),
        need(all(input$n_M_min,input$n_M_mode,input$n_M_max,input$n_M_lambda), "Provide information on the number of Monochamus sampled per inspection site"),
        need(all(input$n_M_min<=input$n_M_mode,input$n_M_min<=input$n_M_max,input$n_M_mode<=input$n_M_max),"Check the information given on the number of Monochamus sampled per inspection site"),
        need(input$site_M, "Define the size of inspection site for Monochamus trapping"),
        need(input$TSe_M,"Define test sensitivity for Monochamus trapping"),
        need(any(input$select_d_w_input_type==1,all(input$d_w_min,input$d_w_mode,input$d_w_max,input$d_w_lambda)),"Porvide information on the density of wood objects suitable for sampling"),
        need(any(input$select_d_M_input_type==1,all(input$d_M_min,input$d_M_mode,input$d_M_max,input$d_M_lambda)),"Porvide information on the density of Monochamus adults"),
        need(all(input$d_w_min<=input$d_w_mode,input$d_w_min<=input$d_w_max,input$d_w_mode<=input$d_w_max),"Check the information given on the density of wood objects suitable for sampling"),
        need(all(input$d_M_min<=input$d_M_mode,input$d_M_min<=input$d_M_max,input$d_M_mode<=input$d_M_max),"Check the information given on the density of Monochamus adults"),
        need(input$file_host_area,"Upload data on the area with host plants"),
        need(input$file_entry_sites,"Upload data on the area of entry sites"),
        need(input$select_survey_type, "Select survey type"),
        # TAALLA
        #need(all(input$DP_w, input$DP_M, any(input$DPr, input$max_inf_size)),"Define design prevalences"),
        need(input$Prior_Pfree, "Determine the prior porbability of freedom"),
        need(all(input$Finv_min,input$Finv_max), "Determine the range of mean time between invasions to be considered"),
        need(input$n_i,"Determine the number of iterations")
      )
    }
    
    req(Pfree_iterations_val())
    Plot_SSe_country(Y(),SSe(),n_r())
  })
  
  # Pfree - Country
  output$Pfree_c <- renderPlot({
    
    if(input$select_n_w_input_type == 1 & input$select_n_M_input_type == 1){
      validate(
        need(input$file_N_w, "Upload data on the number of inspected sites in the wood sampling component"),
        need(input$file_n_w,"Provide information on the number of wood objects sampled per inspection site"),
        need(input$site_w, "Define the size of inspection site for wood sampling"),
        need(input$TSe_w, "Define test sensitivity for wood sampling"),
        need(input$file_N_M, "Upload data on the number of inspected sites in the Monochamus trapping component"),
        need(input$file_n_M,"Provide information on the number of Monochamus sampled per inspection site"),
        need(input$site_M, "Define the size of inspection site for Monochamus trapping"),
        need(input$TSe_M,"Define test sensitivity for Monochamus trapping"),
        need(any(input$select_d_w_input_type==1,all(input$d_w_min,input$d_w_mode,input$d_w_max,input$d_w_lambda)),"Porvide information on the density of wood objects suitable for sampling"),
        need(any(input$select_d_M_input_type==1,all(input$d_M_min,input$d_M_mode,input$d_M_max,input$d_M_lambda)),"Porvide information on the density of Monochamus adults"),
        need(all(input$d_w_min<=input$d_w_mode,input$d_w_min<=input$d_w_max,input$d_w_mode<=input$d_w_max),"Check the information given on the density of wood objects suitable for sampling"),
        need(all(input$d_M_min<=input$d_M_mode,input$d_M_min<=input$d_M_max,input$d_M_mode<=input$d_M_max),"Check the information given on the density of Monochamus adults"),
        need(input$file_host_area,"Upload data on the area with host plants"),
        need(input$file_entry_sites,"Upload data on the area of entry sites"),
        need(input$select_survey_type, "Select survey type"),
        # TAALLA
        #need(all(input$DP_w, input$DP_M, any(input$DPr, input$max_inf_size)),"Define design prevalences"),
        need(input$Prior_Pfree, "Determine the prior porbability of freedom"),
        need(all(input$Finv_min,input$Finv_max), "Determine the range of mean time between invasions to be considered"),
        need(input$n_i,"Determine the number of iterations")
      )
    }else if(input$select_n_w_input_type == 1 & input$select_n_M_input_type == 2){
      validate(
        need(input$file_N_w, "Upload data on the number of inspected sites in the wood sampling component"),
        need(input$file_n_w,"Provide information on the number of wood objects sampled per inspection site"),
        need(input$site_w, "Define the size of inspection site for wood sampling"),
        need(input$TSe_w, "Define test sensitivity for wood sampling"),
        need(input$file_N_M, "Upload data on the number of inspected sites in the Monochamus trapping component"),
        need(all(input$n_M_min,input$n_M_mode,input$n_M_max,input$n_M_lambda), "Provide information on the number of Monochamus sampled per inspection site"),
        need(all(input$n_M_min<=input$n_M_mode,input$n_M_min<=input$n_M_max,input$n_M_mode<=input$n_M_max),"Check the information given on the number of Monochamus sampled per inspection site"),
        need(input$site_M, "Define the size of inspection site for Monochamus trapping"),
        need(input$TSe_M,"Define test sensitivity for Monochamus trapping"),
        need(any(input$select_d_w_input_type==1,all(input$d_w_min,input$d_w_mode,input$d_w_max,input$d_w_lambda)),"Porvide information on the density of wood objects suitable for sampling"),
        need(any(input$select_d_M_input_type==1,all(input$d_M_min,input$d_M_mode,input$d_M_max,input$d_M_lambda)),"Porvide information on the density of Monochamus adults"),
        need(all(input$d_w_min<=input$d_w_mode,input$d_w_min<=input$d_w_max,input$d_w_mode<=input$d_w_max),"Check the information given on the density of wood objects suitable for sampling"),
        need(all(input$d_M_min<=input$d_M_mode,input$d_M_min<=input$d_M_max,input$d_M_mode<=input$d_M_max),"Check the information given on the density of Monochamus adults"),
        need(input$file_host_area,"Upload data on the area with host plants"),
        need(input$file_entry_sites,"Upload data on the area of entry sites"),
        need(input$select_survey_type, "Select survey type"),
        # TAALLA
        #need(all(input$DP_w, input$DP_M, any(input$DPr, input$max_inf_size)),"Define design prevalences"),
        need(input$Prior_Pfree, "Determine the prior porbability of freedom"),
        need(all(input$Finv_min,input$Finv_max), "Determine the range of mean time between invasions to be considered"),
        need(input$n_i,"Determine the number of iterations")
      )
    }else if(input$select_n_w_input_type == 2 & input$select_n_M_input_type == 1){
      validate(
        need(input$file_N_w, "Upload data on the number of inspected sites in the wood sampling component"),
        need(all(input$n_w_min,input$n_w_mode,input$n_w_max,input$n_w_lambda),"Provide information on the number of wood objects sampled per inspection site"),
        need(all(input$n_w_min<=input$n_w_mode,input$n_w_min<=input$n_w_max,input$n_w_mode<=input$n_w_max),"Check the information given on the number of wood objects sampled per inspection site"),
        need(input$site_w, "Define the size of inspection site for wood sampling"),
        need(input$TSe_w, "Define test sensitivity for wood sampling"),
        need(input$file_N_M, "Upload data on the number of inspected sites in the Monochamus trapping component"),
        need(input$file_n_M,"Provide information on the number of Monochamus sampled per inspection site"),
        need(input$site_M, "Define the size of inspection site for Monochamus trapping"),
        need(input$TSe_M,"Define test sensitivity for Monochamus trapping"),
        need(any(input$select_d_w_input_type==1,all(input$d_w_min,input$d_w_mode,input$d_w_max,input$d_w_lambda)),"Porvide information on the density of wood objects suitable for sampling"),
        need(any(input$select_d_M_input_type==1,all(input$d_M_min,input$d_M_mode,input$d_M_max,input$d_M_lambda)),"Porvide information on the density of Monochamus adults"),
        need(all(input$d_w_min<=input$d_w_mode,input$d_w_min<=input$d_w_max,input$d_w_mode<=input$d_w_max),"Check the information given on the density of wood objects suitable for sampling"),
        need(all(input$d_M_min<=input$d_M_mode,input$d_M_min<=input$d_M_max,input$d_M_mode<=input$d_M_max),"Check the information given on the density of Monochamus adults"),
        need(input$file_host_area,"Upload data on the area with host plants"),
        need(input$file_entry_sites,"Upload data on the area of entry sites"),
        need(input$select_survey_type, "Select survey type"),
        # TAALLA
        #need(all(input$DP_w, input$DP_M, any(input$DPr, input$max_inf_size)),"Define design prevalences"),
        need(input$Prior_Pfree, "Determine the prior porbability of freedom"),
        need(all(input$Finv_min,input$Finv_max), "Determine the range of mean time between invasions to be considered"),
        need(input$n_i,"Determine the number of iterations") 
      )
    }else if(input$select_n_w_input_type == 2 & input$select_n_M_input_type == 2){
      validate(
        need(input$file_N_w, "Upload data on the number of inspected sites in the wood sampling component"),
        need(all(input$n_w_min,input$n_w_mode,input$n_w_max,input$n_w_lambda),"Provide information on the number of wood objects sampled per inspection site"),
        need(all(input$n_w_min<=input$n_w_mode,input$n_w_min<=input$n_w_max,input$n_w_mode<=input$n_w_max),"Check the information given on the number of wood objects sampled per inspection site"),
        need(input$site_w, "Define the size of inspection site for wood sampling"),
        need(input$TSe_w, "Define test sensitivity for wood sampling"),
        need(input$file_N_M, "Upload data on the number of inspected sites in the Monochamus trapping component"),
        need(all(input$n_M_min,input$n_M_mode,input$n_M_max,input$n_M_lambda), "Provide information on the number of Monochamus sampled per inspection site"),
        need(all(input$n_M_min<=input$n_M_mode,input$n_M_min<=input$n_M_max,input$n_M_mode<=input$n_M_max),"Check the information given on the number of Monochamus sampled per inspection site"),
        need(input$site_M, "Define the size of inspection site for Monochamus trapping"),
        need(input$TSe_M,"Define test sensitivity for Monochamus trapping"),
        need(any(input$select_d_w_input_type==1,all(input$d_w_min,input$d_w_mode,input$d_w_max,input$d_w_lambda)),"Porvide information on the density of wood objects suitable for sampling"),
        need(any(input$select_d_M_input_type==1,all(input$d_M_min,input$d_M_mode,input$d_M_max,input$d_M_lambda)),"Porvide information on the density of Monochamus adults"),
        need(all(input$d_w_min<=input$d_w_mode,input$d_w_min<=input$d_w_max,input$d_w_mode<=input$d_w_max),"Check the information given on the density of wood objects suitable for sampling"),
        need(all(input$d_M_min<=input$d_M_mode,input$d_M_min<=input$d_M_max,input$d_M_mode<=input$d_M_max),"Check the information given on the density of Monochamus adults"),
        need(input$file_host_area,"Upload data on the area with host plants"),
        need(input$file_entry_sites,"Upload data on the area of entry sites"),
        need(input$select_survey_type, "Select survey type"),
        need(all(input$DP_w, input$DP_M, any(input$DPr, input$max_inf_size)),"Define design prevalences"),
        need(input$Prior_Pfree, "Determine the prior porbability of freedom"),
        need(all(input$Finv_min,input$Finv_max), "Determine the range of mean time between invasions to be considered"),
        need(input$n_i,"Determine the number of iterations")
      )
    }
    
    req(Pfree_iterations_val())
    Plot_Pfree_country(Finv_FI(),Pfree(),n_r())
  })
  
  # Pfree - Regions
  output$Pfree_r <- renderPlot({
    
    if(input$select_n_w_input_type == 1 & input$select_n_M_input_type == 1){
      validate(
        need(input$file_N_w, "Upload data on the number of inspected sites in the wood sampling component"),
        need(input$file_n_w,"Provide information on the number of wood objects sampled per inspection site"),
        need(input$site_w, "Define the size of inspection site for wood sampling"),
        need(input$TSe_w, "Define test sensitivity for wood sampling"),
        need(input$file_N_M, "Upload data on the number of inspected sites in the Monochamus trapping component"),
        need(input$file_n_M,"Provide information on the number of Monochamus sampled per inspection site"),
        need(input$site_M, "Define the size of inspection site for Monochamus trapping"),
        need(input$TSe_M,"Define test sensitivity for Monochamus trapping"),
        need(any(input$select_d_w_input_type==1,all(input$d_w_min,input$d_w_mode,input$d_w_max,input$d_w_lambda)),"Porvide information on the density of wood objects suitable for sampling"),
        need(any(input$select_d_M_input_type==1,all(input$d_M_min,input$d_M_mode,input$d_M_max,input$d_M_lambda)),"Porvide information on the density of Monochamus adults"),
        need(all(input$d_w_min<=input$d_w_mode,input$d_w_min<=input$d_w_max,input$d_w_mode<=input$d_w_max),"Check the information given on the density of wood objects suitable for sampling"),
        need(all(input$d_M_min<=input$d_M_mode,input$d_M_min<=input$d_M_max,input$d_M_mode<=input$d_M_max),"Check the information given on the density of Monochamus adults"),
        need(input$file_host_area,"Upload data on the area with host plants"),
        need(input$file_entry_sites,"Upload data on the area of entry sites"),
        need(input$select_survey_type, "Select survey type"),
        need(all(input$DP_w, input$DP_M, any(input$DPr, input$max_inf_size)),"Define design prevalences"),
        need(input$Prior_Pfree, "Determine the prior porbability of freedom"),
        need(all(input$Finv_min,input$Finv_max), "Determine the range of mean time between invasions to be considered"),
        need(input$n_i,"Determine the number of iterations")
      )
    }else if(input$select_n_w_input_type == 1 & input$select_n_M_input_type == 2){
      validate(
        need(input$file_N_w, "Upload data on the number of inspected sites in the wood sampling component"),
        need(input$file_n_w,"Provide information on the number of wood objects sampled per inspection site"),
        need(input$site_w, "Define the size of inspection site for wood sampling"),
        need(input$TSe_w, "Define test sensitivity for wood sampling"),
        need(input$file_N_M, "Upload data on the number of inspected sites in the Monochamus trapping component"),
        need(all(input$n_M_min,input$n_M_mode,input$n_M_max,input$n_M_lambda), "Provide information on the number of Monochamus sampled per inspection site"),
        need(all(input$n_M_min<=input$n_M_mode,input$n_M_min<=input$n_M_max,input$n_M_mode<=input$n_M_max),"Check the information given on the number of Monochamus sampled per inspection site"),
        need(input$site_M, "Define the size of inspection site for Monochamus trapping"),
        need(input$TSe_M,"Define test sensitivity for Monochamus trapping"),
        need(any(input$select_d_w_input_type==1,all(input$d_w_min,input$d_w_mode,input$d_w_max,input$d_w_lambda)),"Porvide information on the density of wood objects suitable for sampling"),
        need(any(input$select_d_M_input_type==1,all(input$d_M_min,input$d_M_mode,input$d_M_max,input$d_M_lambda)),"Porvide information on the density of Monochamus adults"),
        need(all(input$d_w_min<=input$d_w_mode,input$d_w_min<=input$d_w_max,input$d_w_mode<=input$d_w_max),"Check the information given on the density of wood objects suitable for sampling"),
        need(all(input$d_M_min<=input$d_M_mode,input$d_M_min<=input$d_M_max,input$d_M_mode<=input$d_M_max),"Check the information given on the density of Monochamus adults"),
        need(input$file_host_area,"Upload data on the area with host plants"),
        need(input$file_entry_sites,"Upload data on the area of entry sites"),
        need(input$select_survey_type, "Select survey type"),
        need(all(input$DP_w, input$DP_M, any(input$DPr, input$max_inf_size)),"Define design prevalences"),
        need(input$Prior_Pfree, "Determine the prior porbability of freedom"),
        need(all(input$Finv_min,input$Finv_max), "Determine the range of mean time between invasions to be considered"),
        need(input$n_i,"Determine the number of iterations")
      )
    }else if(input$select_n_w_input_type == 2 & input$select_n_M_input_type == 1){
      validate(
        need(input$file_N_w, "Upload data on the number of inspected sites in the wood sampling component"),
        need(all(input$n_w_min,input$n_w_mode,input$n_w_max,input$n_w_lambda),"Provide information on the number of wood objects sampled per inspection site"),
        need(all(input$n_w_min<=input$n_w_mode,input$n_w_min<=input$n_w_max,input$n_w_mode<=input$n_w_max),"Check the information given on the number of wood objects sampled per inspection site"),
        need(input$site_w, "Define the size of inspection site for wood sampling"),
        need(input$TSe_w, "Define test sensitivity for wood sampling"),
        need(input$file_N_M, "Upload data on the number of inspected sites in the Monochamus trapping component"),
        need(input$file_n_M,"Provide information on the number of Monochamus sampled per inspection site"),
        need(input$site_M, "Define the size of inspection site for Monochamus trapping"),
        need(input$TSe_M,"Define test sensitivity for Monochamus trapping"),
        need(any(input$select_d_w_input_type==1,all(input$d_w_min,input$d_w_mode,input$d_w_max,input$d_w_lambda)),"Porvide information on the density of wood objects suitable for sampling"),
        need(any(input$select_d_M_input_type==1,all(input$d_M_min,input$d_M_mode,input$d_M_max,input$d_M_lambda)),"Porvide information on the density of Monochamus adults"),
        need(all(input$d_w_min<=input$d_w_mode,input$d_w_min<=input$d_w_max,input$d_w_mode<=input$d_w_max),"Check the information given on the density of wood objects suitable for sampling"),
        need(all(input$d_M_min<=input$d_M_mode,input$d_M_min<=input$d_M_max,input$d_M_mode<=input$d_M_max),"Check the information given on the density of Monochamus adults"),
        need(input$file_host_area,"Upload data on the area with host plants"),
        need(input$file_entry_sites,"Upload data on the area of entry sites"),
        need(input$select_survey_type, "Select survey type"),
        need(all(input$DP_w, input$DP_M, any(input$DPr, input$max_inf_size)),"Define design prevalences"),
        need(input$Prior_Pfree, "Determine the prior porbability of freedom"),
        need(all(input$Finv_min,input$Finv_max), "Determine the range of mean time between invasions to be considered"),
        need(input$n_i,"Determine the number of iterations") 
      )
    }else if(input$select_n_w_input_type == 2 & input$select_n_M_input_type == 2){
      validate(
        need(input$file_N_w, "Upload data on the number of inspected sites in the wood sampling component"),
        need(all(input$n_w_min,input$n_w_mode,input$n_w_max,input$n_w_lambda),"Provide information on the number of wood objects sampled per inspection site"),
        need(all(input$n_w_min<=input$n_w_mode,input$n_w_min<=input$n_w_max,input$n_w_mode<=input$n_w_max),"Check the information given on the number of wood objects sampled per inspection site"),
        need(input$site_w, "Define the size of inspection site for wood sampling"),
        need(input$TSe_w, "Define test sensitivity for wood sampling"),
        need(input$file_N_M, "Upload data on the number of inspected sites in the Monochamus trapping component"),
        need(all(input$n_M_min,input$n_M_mode,input$n_M_max,input$n_M_lambda), "Provide information on the number of Monochamus sampled per inspection site"),
        need(all(input$n_M_min<=input$n_M_mode,input$n_M_min<=input$n_M_max,input$n_M_mode<=input$n_M_max),"Check the information given on the number of Monochamus sampled per inspection site"),
        need(input$site_M, "Define the size of inspection site for Monochamus trapping"),
        need(input$TSe_M,"Define test sensitivity for Monochamus trapping"),
        need(any(input$select_d_w_input_type==1,all(input$d_w_min,input$d_w_mode,input$d_w_max,input$d_w_lambda)),"Porvide information on the density of wood objects suitable for sampling"),
        need(any(input$select_d_M_input_type==1,all(input$d_M_min,input$d_M_mode,input$d_M_max,input$d_M_lambda)),"Porvide information on the density of Monochamus adults"),
        need(all(input$d_w_min<=input$d_w_mode,input$d_w_min<=input$d_w_max,input$d_w_mode<=input$d_w_max),"Check the information given on the density of wood objects suitable for sampling"),
        need(all(input$d_M_min<=input$d_M_mode,input$d_M_min<=input$d_M_max,input$d_M_mode<=input$d_M_max),"Check the information given on the density of Monochamus adults"),
        need(input$file_host_area,"Upload data on the area with host plants"),
        need(input$file_entry_sites,"Upload data on the area of entry sites"),
        need(input$select_survey_type, "Select survey type"),
        need(all(input$DP_w, input$DP_M, any(input$DPr, input$max_inf_size)),"Define design prevalences"),
        need(input$Prior_Pfree, "Determine the prior porbability of freedom"),
        need(all(input$Finv_min,input$Finv_max), "Determine the range of mean time between invasions to be considered"),
        need(input$n_i,"Determine the number of iterations")
      )
    }
    
    req(Pfree_iterations_val())
    Plot_Pfree_regions(Finv_FI(),Pfree(),region_names())
  })
  
  # SSe - Regions
  output$SSe_r <- renderPlot({
    
    if(input$select_n_w_input_type == 1 & input$select_n_M_input_type == 1){
      validate(
        need(input$file_N_w, "Upload data on the number of inspected sites in the wood sampling component"),
        need(input$file_n_w,"Provide information on the number of wood objects sampled per inspection site"),
        need(input$site_w, "Define the size of inspection site for wood sampling"),
        need(input$TSe_w, "Define test sensitivity for wood sampling"),
        need(input$file_N_M, "Upload data on the number of inspected sites in the Monochamus trapping component"),
        need(input$file_n_M,"Provide information on the number of Monochamus sampled per inspection site"),
        need(input$site_M, "Define the size of inspection site for Monochamus trapping"),
        need(input$TSe_M,"Define test sensitivity for Monochamus trapping"),
        need(any(input$select_d_w_input_type==1,all(input$d_w_min,input$d_w_mode,input$d_w_max,input$d_w_lambda)),"Porvide information on the density of wood objects suitable for sampling"),
        need(any(input$select_d_M_input_type==1,all(input$d_M_min,input$d_M_mode,input$d_M_max,input$d_M_lambda)),"Porvide information on the density of Monochamus adults"),
        need(all(input$d_w_min<=input$d_w_mode,input$d_w_min<=input$d_w_max,input$d_w_mode<=input$d_w_max),"Check the information given on the density of wood objects suitable for sampling"),
        need(all(input$d_M_min<=input$d_M_mode,input$d_M_min<=input$d_M_max,input$d_M_mode<=input$d_M_max),"Check the information given on the density of Monochamus adults"),
        need(input$file_host_area,"Upload data on the area with host plants"),
        need(input$file_entry_sites,"Upload data on the area of entry sites"),
        need(input$select_survey_type, "Select survey type"),
        need(all(input$DP_w, input$DP_M, any(input$DPr, input$max_inf_size)),"Define design prevalences"),
        need(input$Prior_Pfree, "Determine the prior porbability of freedom"),
        need(all(input$Finv_min,input$Finv_max), "Determine the range of mean time between invasions to be considered"),
        need(input$n_i,"Determine the number of iterations")
      )
    }else if(input$select_n_w_input_type == 1 & input$select_n_M_input_type == 2){
      validate(
        need(input$file_N_w, "Upload data on the number of inspected sites in the wood sampling component"),
        need(input$file_n_w,"Provide information on the number of wood objects sampled per inspection site"),
        need(input$site_w, "Define the size of inspection site for wood sampling"),
        need(input$TSe_w, "Define test sensitivity for wood sampling"),
        need(input$file_N_M, "Upload data on the number of inspected sites in the Monochamus trapping component"),
        need(all(input$n_M_min,input$n_M_mode,input$n_M_max,input$n_M_lambda), "Provide information on the number of Monochamus sampled per inspection site"),
        need(all(input$n_M_min<=input$n_M_mode,input$n_M_min<=input$n_M_max,input$n_M_mode<=input$n_M_max),"Check the information given on the number of Monochamus sampled per inspection site"),
        need(input$site_M, "Define the size of inspection site for Monochamus trapping"),
        need(input$TSe_M,"Define test sensitivity for Monochamus trapping"),
        need(any(input$select_d_w_input_type==1,all(input$d_w_min,input$d_w_mode,input$d_w_max,input$d_w_lambda)),"Porvide information on the density of wood objects suitable for sampling"),
        need(any(input$select_d_M_input_type==1,all(input$d_M_min,input$d_M_mode,input$d_M_max,input$d_M_lambda)),"Porvide information on the density of Monochamus adults"),
        need(all(input$d_w_min<=input$d_w_mode,input$d_w_min<=input$d_w_max,input$d_w_mode<=input$d_w_max),"Check the information given on the density of wood objects suitable for sampling"),
        need(all(input$d_M_min<=input$d_M_mode,input$d_M_min<=input$d_M_max,input$d_M_mode<=input$d_M_max),"Check the information given on the density of Monochamus adults"),
        need(input$file_host_area,"Upload data on the area with host plants"),
        need(input$file_entry_sites,"Upload data on the area of entry sites"),
        need(input$select_survey_type, "Select survey type"),
        need(all(input$DP_w, input$DP_M, any(input$DPr, input$max_inf_size)),"Define design prevalences"),
        need(input$Prior_Pfree, "Determine the prior porbability of freedom"),
        need(all(input$Finv_min,input$Finv_max), "Determine the range of mean time between invasions to be considered"),
        need(input$n_i,"Determine the number of iterations")
      )
    }else if(input$select_n_w_input_type == 2 & input$select_n_M_input_type == 1){
      validate(
        need(input$file_N_w, "Upload data on the number of inspected sites in the wood sampling component"),
        need(all(input$n_w_min,input$n_w_mode,input$n_w_max,input$n_w_lambda),"Provide information on the number of wood objects sampled per inspection site"),
        need(all(input$n_w_min<=input$n_w_mode,input$n_w_min<=input$n_w_max,input$n_w_mode<=input$n_w_max),"Check the information given on the number of wood objects sampled per inspection site"),
        need(input$site_w, "Define the size of inspection site for wood sampling"),
        need(input$TSe_w, "Define test sensitivity for wood sampling"),
        need(input$file_N_M, "Upload data on the number of inspected sites in the Monochamus trapping component"),
        need(input$file_n_M,"Provide information on the number of Monochamus sampled per inspection site"),
        need(input$site_M, "Define the size of inspection site for Monochamus trapping"),
        need(input$TSe_M,"Define test sensitivity for Monochamus trapping"),
        need(any(input$select_d_w_input_type==1,all(input$d_w_min,input$d_w_mode,input$d_w_max,input$d_w_lambda)),"Porvide information on the density of wood objects suitable for sampling"),
        need(any(input$select_d_M_input_type==1,all(input$d_M_min,input$d_M_mode,input$d_M_max,input$d_M_lambda)),"Porvide information on the density of Monochamus adults"),
        need(all(input$d_w_min<=input$d_w_mode,input$d_w_min<=input$d_w_max,input$d_w_mode<=input$d_w_max),"Check the information given on the density of wood objects suitable for sampling"),
        need(all(input$d_M_min<=input$d_M_mode,input$d_M_min<=input$d_M_max,input$d_M_mode<=input$d_M_max),"Check the information given on the density of Monochamus adults"),
        need(input$file_host_area,"Upload data on the area with host plants"),
        need(input$file_entry_sites,"Upload data on the area of entry sites"),
        need(input$select_survey_type, "Select survey type"),
        need(all(input$DP_w, input$DP_M, any(input$DPr, input$max_inf_size)),"Define design prevalences"),
        need(input$Prior_Pfree, "Determine the prior porbability of freedom"),
        need(all(input$Finv_min,input$Finv_max), "Determine the range of mean time between invasions to be considered"),
        need(input$n_i,"Determine the number of iterations") 
      )
    }else if(input$select_n_w_input_type == 2 & input$select_n_M_input_type == 2){
      validate(
        need(input$file_N_w, "Upload data on the number of inspected sites in the wood sampling component"),
        need(all(input$n_w_min,input$n_w_mode,input$n_w_max,input$n_w_lambda),"Provide information on the number of wood objects sampled per inspection site"),
        need(all(input$n_w_min<=input$n_w_mode,input$n_w_min<=input$n_w_max,input$n_w_mode<=input$n_w_max),"Check the information given on the number of wood objects sampled per inspection site"),
        need(input$site_w, "Define the size of inspection site for wood sampling"),
        need(input$TSe_w, "Define test sensitivity for wood sampling"),
        need(input$file_N_M, "Upload data on the number of inspected sites in the Monochamus trapping component"),
        need(all(input$n_M_min,input$n_M_mode,input$n_M_max,input$n_M_lambda), "Provide information on the number of Monochamus sampled per inspection site"),
        need(all(input$n_M_min<=input$n_M_mode,input$n_M_min<=input$n_M_max,input$n_M_mode<=input$n_M_max),"Check the information given on the number of Monochamus sampled per inspection site"),
        need(input$site_M, "Define the size of inspection site for Monochamus trapping"),
        need(input$TSe_M,"Define test sensitivity for Monochamus trapping"),
        need(any(input$select_d_w_input_type==1,all(input$d_w_min,input$d_w_mode,input$d_w_max,input$d_w_lambda)),"Porvide information on the density of wood objects suitable for sampling"),
        need(any(input$select_d_M_input_type==1,all(input$d_M_min,input$d_M_mode,input$d_M_max,input$d_M_lambda)),"Porvide information on the density of Monochamus adults"),
        need(all(input$d_w_min<=input$d_w_mode,input$d_w_min<=input$d_w_max,input$d_w_mode<=input$d_w_max),"Check the information given on the density of wood objects suitable for sampling"),
        need(all(input$d_M_min<=input$d_M_mode,input$d_M_min<=input$d_M_max,input$d_M_mode<=input$d_M_max),"Check the information given on the density of Monochamus adults"),
        need(input$file_host_area,"Upload data on the area with host plants"),
        need(input$file_entry_sites,"Upload data on the area of entry sites"),
        need(input$select_survey_type, "Select survey type"),
        need(all(input$DP_w, input$DP_M, any(input$DPr, input$max_inf_size)),"Define design prevalences"),
        need(input$Prior_Pfree, "Determine the prior porbability of freedom"),
        need(all(input$Finv_min,input$Finv_max), "Determine the range of mean time between invasions to be considered"),
        need(input$n_i,"Determine the number of iterations")
      )
    }
    
    req(Pfree_iterations_val())
    Plot_SSe_regions(Y(),SSe(),region_names())
  })
  
  ##########
  # DOWNLOAD RESULTS
  
  # SSe tables
  SSetableInput <- reactive({
    
    a <- as.data.frame(SSe()[,,1])
    b <- as.data.frame(SSe()[,,2])
    c <- as.data.frame(SSe()[,,3])
    
    colnames(a)[seq(1,ncol(a)-1,1)] <- c(region_names())
    colnames(a)[ncol(a)] <- "Country"
    colnames(b) <- colnames(a)
    colnames(c) <- colnames(a)
    
    switch(input$SSe_table_to_dl,
           "2.5%" = a,
           "50%" = b,
           "97.5%" = c
    )
  })
  
  output$download_SSe_tables <- downloadHandler(
    filename = function(){
      paste(input$SSe_table_to_dl, 'csv', sep = ".")},
    content = function(file){
      write.csv(SSetableInput(), file, row.names = c(Y()))
    })
  
  output$download_SSe_fractiles <- downloadHandler(
    filename = function(){
      paste("SSe_fractiles.rds")},
    content = function(file){
      saveRDS(SSe(), file=file)
    })
  
  output$download_SSe_iterations <- downloadHandler(
    filename = function(){
      paste("SSe_iterations.rds")},
    content = function(file){
      saveRDS(SSe_iterations(), file=file)
    })
  
  # Pfree tables
  PfreetableInput <- reactive({
    
    a <- as.data.frame(Pfree()[,,1])
    b <- as.data.frame(Pfree()[,,2])
    c <- as.data.frame(Pfree()[,,3])
    
    colnames(a)[seq(1,ncol(a)-1,1)] <- c(region_names())
    colnames(a)[ncol(a)] <- "Country"
    colnames(b) <- colnames(a)
    colnames(c) <- colnames(a)
    
    switch(input$Pfree_table_to_dl,
           "2.5%" = a,
           "50%" = b,
           "97.5%" = c)
  })
  
  output$download_Pfree_tables <- downloadHandler(
    filename = function(){
      paste(input$Pfree_table_to_dl,"csv", sep = ".")},
    content = function(file){
      write.csv(PfreetableInput(), file, row.names = c(Finv_FI()))
    })
  
  output$download_Pfree_fractiles <- downloadHandler(
    filename = function(){
      paste("Pfree_fractiles.rds")},
    content = function(file){
      saveRDS(Pfree(), file=file)
    })
  
  output$download_Pfree_iterations <- downloadHandler(
    filename = function(){
      paste("Pfree_iterations.rds")},
    content = function(file){
      saveRDS(Pfree_iterations_val(), file=file)
    })
  
  # Download figures
  FigInput <- reactive({
    switch(input$fig_to_dl,
           "Sensitivity - Country" = Plot_SSe_country(Y(),SSe(),n_r()),
           "Probability of freedom - Country" = Plot_Pfree_country(Finv_FI(),Pfree(),n_r()),
           "Sensitivity - Regions" = Plot_SSe_regions(Y(),SSe(),region_names()),
           "Probability of freedom - Regions" = Plot_Pfree_regions(Finv_FI(),Pfree(),region_names())
    )}) 
  
  output$download_fig <- downloadHandler(
    filename = function(){
      paste(input$fig_to_dl,"png",sep = ".")
    },
    content = function(file){
      ggsave(FigInput(), file = file)
    })
} 
############################################################
shinyApp(ui=ui,server=server)