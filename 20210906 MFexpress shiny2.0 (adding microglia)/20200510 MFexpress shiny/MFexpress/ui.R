library(shiny)

mobileDetect <- function(inputId, value = 0) {
  tagList(
    singleton(tags$head(tags$script(src = "js/mobile.js"))),
    tags$input(id = inputId,
               class = "mobile-element",
               type = "hidden")
  )
}

shinyUI(
  fluidPage(
    
    #unshown function    
    
    mobileDetect('isMobile'),
    
    #App out  
    
    headerPanel("MFexpress : shortcut to macrophage gene expression"),
    
    sidebarPanel(
      
      textInput("Gene", "Please enter gene symbol or Entrez ID", "ENT3"),
      actionButton("enter", "Express !"),
      helpText("If showing unmatched, try the official gene symbol."),
      helpText("Microglia new added, several genes might be missing due to different gene name transformation process. Last updated by TTL 2021.09.06")
      
    ),
    
    mainPanel(plotOutput("myPlot")),
    
  )
)