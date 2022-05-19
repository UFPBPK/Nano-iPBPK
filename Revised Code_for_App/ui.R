## Loading required r package
library(shiny)          # R package for shiny App
library(shinydashboard) # R package for shiny App UI structure
library(shinyjs)        # R package for hide and show a shiny element
library(mrgsolve)       # R package for PBPK model
library(ggplot2)        # R package for plot output
library(truncnorm)      # R package for Truncated normal distribution
library(EnvStats)       # R Package for Environmental Statistics, Including US EPA Guidance
library(rmarkdown)      # R package for output report
library(dplyr)          # R package for dataframe manupulation
library(knitr)          # R package for output reportsource("function_2.R")
library(berryFunctions) # R package for the function 'insertRows'

## Loading mrgsolve code
source("MrgCode.R")


## Sidebar: make your main item
sidebar <- dashboardSidebar(
  hr(),
  sidebarMenu(menuItem("Main plot", tabName = "plot", icon = icon("microchip"), selected = TRUE),
              menuItem("Output Table", tabName = 'table', icon = icon("table")),
              menuItem("Model Structure", tabName = 'modelstructure', icon = icon("area-chart")),
              menuItem("Code", tabName = 'code', icon = icon("code")),
              #menuItem("Tutorial", tabName = "tutorial", icon = icon("mortar-board")),
              #menuItem("Output Report", tabName = "report", icon = icon("file-pdf-o")),
              menuItem("About", tabName = 'about', icon = icon("question-circle"))
  ),
  hr()
  )


## Body: the content of menuItem
body <- dashboardBody(
  tabItems( # to match upt the menuItem with tabItems 
    
    ## First item content
    tabItem( 
      tabName = 'plot', # the tabName have to match with the tabName of "menuItem"
      sidebarLayout( # sidebarLayout: layout a sidebar (sidebarPanel) and main area (mainPanel)
        sidebarPanel(width = 2,
                     fluidRow(
                       useShinyjs(),
                       tags$h4("Parameters for Treatment Scenario"),
                       div(
                         id = "Treatment",
                         selectInput(inputId = 'route', label = 'Administration route', 
                                     choices = c('Oral', 'IV','IT','IH'), width = 150, selected = 'Oral'
                         ),
                         selectInput(inputId = 'target', label = 'Target tissue', 
                                     choices = c('Blood', 'Lung','Spleen','Liver','Kidney','GI'), width = 150, selected = 'Blood'
                         ),
                         numericInput(inputId = 'doselevel', label = 'Dose level (mg/kg)', value = 20,
                                      min = 0, step = 1, width = 150
                         ),
                         conditionalPanel(condition = "input.route == 'Oral'",
                                          numericInput(inputId = 'size_oral', label = 'NPs Size (nm)',
                                                       value = 5, min = 1, max = 200, step = 1, width = 150
                                          )),
                         conditionalPanel(condition = "input.route == 'IV'",
                                          numericInput(inputId = 'size_IV', label = 'NPs Size (nm)', 
                                                       value = 5, min = 1, max = 200, step = 1, width = 150
                                          )),
                         conditionalPanel(condition = "input.route == 'IT'",
                                          numericInput(inputId = 'size_IT', label = 'NPs Size (nm)', 
                                                       value = 5, min = 1, max = 200, step = 1, width = 150
                                          )),
                         conditionalPanel(condition = "input.route == 'IH'",
                                          selectInput(inputId = 'size_IH', label = 'NPs Size (nm)', choices = c('23 nm'),
                                                      selected = '23 nm', width = 150
                                          )),
                         
                         conditionalPanel(condition = "input.route != 'IH'",
                         numericInput(inputId = 'Zeta', label = 'Zeta potential (mV)', value = -20.6,
                                      min = -90, step = 1, width = 150
                         ),
                         numericInput(inputId = 'HD', label = 'Hydrodynamic diameter (nm)', value = 12.1,
                                      min = 0, step = 1, width = 150
                         ),
                         numericInput(inputId = 'SA', label = div("Surface area (cm", tags$sup("2"), ")"), value = 27,
                                      min = 0, step = 1, width = 150
                         ),
                         numericInput(inputId = 'nNPs', label = 'Admin. number of NPs per rat (10^)', value = 13,
                                      min = 1, step = 1, width = 150
                         )),
                         
                         tags$p(" "),
                         
                         actionButton(inputId = 'action',label = 'Apply Changes',icon("paper-plane"), style="width:200px"),
                         
                         actionButton(inputId = 'reset', label = 'Default Values',icon("refresh"), style="width:200px") 
                       ))
        ), # end of sidebarPanel
        
        mainPanel(
          fluidRow(
            box(title = 'Time-varying profiles for AuNPs in target organ',
                plotOutput(outputId = 'tarplot'), width = 12, status = "primary", solidHeader = TRUE,
                tags$p(" "),
                downloadButton('downb','Download', width = 110),
                selectInput(inputId = 'unit', label = 'Change Unit', choices = c('%ID', 'Conc.'), 
                                        width = 110, selected = '%ID')
            ), 
            box( width = 6, status = "primary", solidHeader = TRUE,        
                 sliderInput("tinterval", "Dose interval (h):", value = 24, min = 1, max = 50, step = 1),
                 sliderInput("numdose", "Number of doses:", value = 1, min = 1, max = 20, step = 1),
                 sliderInput("simu_time", "Simulation time after last administration (day):", value = 1, min = 0, max = 10, step = 1),
            ),
            box( width = 6, status = "primary", solidHeader = TRUE,         
                 sliderInput("n_sim", "Number of iterations:", value = 100, min = 1, max = 2000, step=1),
                 p(actionButton("recalc","Re-run simulation", icon("random"))) #, style = "color: black; background-color: #35e51d"))
            )
          ), width = 10
        ) # end of mainPanel
      ) # end of sidebarLayout
    ), # end of the tabItem of "plot"
    
    
    ## Second item content
    tabItem(
      tabName = 'table',
      box(width = "100%", height = "100%", status = "primary", solidHeader = TRUE, title = "Concentration in tissues",
          downloadButton('downloadtable','Download'),
          tableOutput('table1')
      )
    ), # end ot the tabItem of "table"
    
    ## Third item content
    tabItem(
      tabName = 'modelstructure',
      fluidRow(
        box(width = "100%",status = "primary",title= "Model Structure", align="center", 
            solidHeader = TRUE, tags$img(src = 'ModStructure.png',width = '40%',height = '40%'),
            tags$h2(strong("
                   Fig. 1. A schematic diagram of the physiologically based pharmacokinetic (PBPK)
                   model for 1.4, 5, 18, 23, 80, and 200 nm gold nanoparticles (AuNPs) in adult rats.
                   "), align = 'center', style = "font-family: 'Times', serif; font-weight: 500; 
                   font-size: 500; line-height: 1; color: #404040;")
        )
      )
    ),
    
    ## Fourth item content
    tabItem(
      tabName = 'code',
      box(width = NULL, status = 'primary', solidHeader = TRUE, title = 'title',
          downloadButton('downloadcode', 'Download'),
          br(),br(),
          pre(includeText('MrgCode.R'))
      )
    ), # end of the tabItem of 'code'
    
    # # ## 5th item content
    # tabItem(
    #   tabName = 'tutorial',
    #   tabPanel("tutorial",
    #            tags$iframe(style="height:800px; width:100%; scrolling=yes",
    #                        src="./tutorial3.pdf"))
    # ),
    
    ## 6th item content
    # tabItem(
    #   tabName = 'report',
    #   box(width = 8, status = 'primary', solidHeader = TRUE, title = 'Output Report',
    #       radioButtons('format', 'Document format', c('PDF', 'HTML', 'Word'),
    #                    inline = TRUE),
    #       downloadButton('downloadreport', "Download Report")
    #   )
    #), # end of the tabItem of 'report'
    
    ## 7th item content
    tabItem(
      tabName = 'about',
      fluidRow(
        box(width = 10, status = 'primary', align="left",
            tags$iframe(src = './about.html', 
                       width = '100%', height = '800px',
                       frameborder = 0, scrolling = 'auto'),solidHeader = TRUE)
      )
    ) # end of the tabItem of 'about'
  ) # end of the tabItems
) # end of the dashboardbody

# The dashboardPage() function expects three components: 
# header, sidebar, and body            
ui<-dashboardPage( 
  skin = 'purple',
  dashboardHeader(title = "Nano-iPBPK Model"),
  sidebar,
  body
)