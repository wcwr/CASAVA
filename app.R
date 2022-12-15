#Change slider to shinywidgets?

#A GUI tool to nominate drug combinations from CRISPR+drug screens and other external resources

#Load libraries
library(shiny)
library(shiny.semantic)
library(semantic.dashboard)
library(plotly)
library(dplyr)
library(stringr)
library(DT)
library(googlesheets4)
library(readr)
library(shinyalert)
library(shinybusy)



gs4_deauth() #Needed so that googlesheets4 skips requiring an authenticated Google user



ui <- dashboard_page(
  
  dashboardHeader(title=actionButton(inputId = "load_button",
                                     label="Load data",
                                     style="background-color: #8bc9f0; border-color: #2e6da4"),
                  show_menu_button = TRUE,
                  titleWidth = "thin"),
  
  dashboardSidebar(
    side = "left",
    size = "wide",
    
    sidebarMenu(
      HTML("<br>"), #Space at the top
      fluidRow(
        column(width=16,
               align="center",
               offset=3,
               shiny::sliderInput("effectSize", #Needs to be direct from shiny, because shiny.semantic is buggy when sliderInput is inside sidebarMenu
                                  HTML("<b>Minimum effect size</b>"),
                                  value = -0.14,
                                  min = -0.4073,
                                  max = 0.382))), #End fluidRow
      
      
      HTML("<br><br>"),
      
      
      splitLayout(
        selectInput(inputId ="stat_metric",
                    label = "Choose statistic option):", #Select metric for cutoff (FDR or p-value)
                    choices = c("p-value", "FDR"),  width = "98%"),
        
        
        numericInput(inputId="stat_alpha", #Select alpha value (for corresponding selected stat metric)
                     label= paste0("Significance Threshold (", "\u03B1", ")"), 
                     value=0.05,
                     min=0, 
                     max=1, 
                     step=0.01, 
                     width = "98%" ),
        
        #Background color and height of splitLayout box
        style = "background-color:#FFFFFF; height:100px;"),
      
      
      HTML("<br>"),
      
      
      #Check box to cap the number of known drugs per gene
      fluidRow(
        column(width=16,
               align="center",
               offset=3,
               checkbox_input(input_id="capdrugs",
                              label="Limit maximum drugs per gene?",
                              is_marked=FALSE,
                              type = ""),
               
               #Display popup with info when icon is clicked
               useShinyalert(), #Required for popup
               actionLink("info_button", "", icon = icon("info")),
        ), #End column
      ), #End fluidRow
      
      #If the 'capdrugs' box is checked, Show the corresponding slider
      conditionalPanel(
        condition = "input.capdrugs == 1", #NOTE for condition this is ==1 and not ==TRUE
        fluidRow(
          column(width=16,
                 align="center",
                 offset=3,
                 shiny::sliderInput(inputId = "max_drugs", #Needs to be direct from shiny, because shiny.semantic is buggy when sliderInput is inside sidebarMenu
                                    HTML("<b>Max drugs per gene</b>"),
                                    value = 10,
                                    min = 1,
                                    max = 63,
                                    step=1)))
      ),#End conditionalPanel
      
      
      HTML("<br><br><br>"),
      
      
      checkbox_input(input_id = "na_checkbox",
                     label="Remove control gRNAs",
                     is_marked = TRUE),
      
      checkbox_input(input_id = "cancer_drugs_checkbox",
                     label="Only show cancer drugs", 
                     is_marked = FALSE)
    )#End sidebarMenu)
  ),#End dashboardSidebar
  
  
  dashboardBody(
    
    add_busy_bar(color = "forestgreen", centered=TRUE, height="10px"), #Add busy bar here
    
    
    
    #Display number of total combinations
    fluidRow(
      column(
        width=16,
        align="center", offset=3,
        htmlOutput("mytext")
      )
    ),
    
    fluidRow(
      box(plotlyOutput("plot1", height = 850), width = 8),
      tabBox(
        tabs = list(
          list(menu="Plots", content=column(width=16,plotlyOutput("plot2", height=400),
                                            plotlyOutput("plot3", height=400))),
          list(menu="Table", content=DT::dataTableOutput("mytable", height = 800), width = 8),
          list(menu="DepMap Data", content=column(width=16,plotlyOutput("plot4", height=700))))
      )
      
      ,
    )
  ),
  
  
  theme = "cyborg"
)

server <- shinyServer(function(input, output, session) {
  
  
  #Create an eventReactive which loads in the google sheet data after the load_button is pressed:
  mysheetreactive <- eventReactive(input$load_button, {
    
    mydf <- as.data.frame(read_sheet("https://docs.google.com/spreadsheets/d/1yGP9HJllTrkLbCEn2YcmdHf0UvjfEx8L90VHFMMjp0A/edit?usp=sharing"))
    
    #Read in the 'cancer drugs'
    #CHEMBL database, antineoplastic agents, phase 3 or higher, small molecules
    #NOTE: There are *definitely* a few discrepencies due to drug synonyms, take with a slight grain of salt
    cancer_drugs <- as.data.frame(read_sheet("https://docs.google.com/spreadsheets/d/12ZvjQMUspSaXzEZfhMPJbCc0kj2SJMGcnzaBRakvt6k/edit?usp=sharing"))
    
    
    depmap_data <- as.data.frame(read_sheet("https://docs.google.com/spreadsheets/d/11NX3jelUwtEfDAodaPpsNxS0QplMEWn9dyLr0vc3WQw/edit?usp=sharing"))
    
    
    return(list(mydf=mydf,
                cancer_drugs=cancer_drugs,
                depmap_data=depmap_data))
    
  })#End eventReactive
  
  
  
  
  #Make a reactive that adjusts data according to thresholds
  myreactive <- reactive({
    
    req(mysheetreactive()$mydf)
    req(mysheetreactive()$cancer_drugs)
    req(mysheetreactive()$depmap_data)
    req(input$stat_metric)
    req(input$stat_alpha)
    req(input$max_drugs)
    
    
    mydf <- mysheetreactive()$mydf
    
    cancer_drugs <- mysheetreactive()$cancer_drugs
    
    depmap_data <- mysheetreactive()$depmap_data
    
    
    #mydf2 will be initial adjusted data
    if(input$stat_metric=="p-value"){
      mydf2 <- mydf[mydf$NB_effectSize <= input$effectSize &
                      mydf$NB_Sensitizing_P <= input$stat_alpha,]
      
    }else if(input$stat_metric=="FDR"){
      mydf2 <- mydf[mydf$NB_effectSize <= input$effectSize &
                      mydf$NB_Sensitizing_FDR <= input$stat_alpha,]
    } 
    
    
    mydf3 <- mydf2 %>%
      group_by(tx, Gene) %>%
      summarize(known_drug=unique(unlist(str_split(known_drug, pattern = "\\|"))),
                count=n(),
                concat=paste0(tx, "_", Gene))
    
    #set max colors
    max_drugs_per_gene <- input$max_drugs
    
    #grab the max drugs for each 'concat' (treatment_Gene)
    group1 = mydf3 %>%
      group_by(concat) %>%
      filter(row_number() <= max_drugs_per_gene) %>%
      ungroup() %>%
      arrange(concat)
    
    
    #Put the last part into 'if' statement making sure there is at least 1 row in mydf3
    #If not, make mydf4 NULL otherwise there will be errors in the app
    if(nrow(group1>=1)){ 
      
      
      #identify those with unique drugs not in group1
      unique_drugs = mydf3 %>%
        filter(!known_drug %in% group1$known_drug)
      
      # bind group1 to unique_colors (unique_colors goes first)
      # and keep the "max_colors" number of rows for each name
      mydf4 = bind_rows(unique_drugs, group1) %>%
        group_by(concat) %>%
        filter(row_number() <= max_drugs_per_gene) %>%
        ungroup() %>%
        arrange(concat)
      
      #Make a second concat column 
      mydf4$concat2 <- paste0(mydf4$tx, "_", mydf4$known_drug)
      
      #Remove duplicate entires of concat2 (tx_knowndrug)
      mydf4 <- mydf4[!duplicated(mydf4$concat2),]
      
      mydf4 <- subset(mydf4, select= -c(count))
      
      
      
    } else {
      mydf4 <- NULL
    }
    
    
    #If 'na_checkbox' is checked, remove entries where Gene contains the string "controlGene" (aka this is a control non-targeting gRNA)
    #NOTE: This is specific to this dataset; This will not work if control names are changed to "nt_gRNA1", etc. 
    if(input$na_checkbox==TRUE){
      mydf4 <- mydf4[!grepl("controlGene", mydf4$Gene),]
    } else {
      mydf4 <- mydf4
    }
    
    #If 'cancer_drugs_checkbox' is checked, 
    #Remove all non-cancer drugs (CAVEATS, see above about chembl data)
    if(input$cancer_drugs_checkbox==TRUE){
      mydf4 <- mydf4[mydf4$known_drug %in% cancer_drugs$known_drug,]
    } else {
      mydf4 <- mydf4
    }
    
    
    depmap_data2 <- merge(mydf4, depmap_data, by="Gene")
    
    
    return(list(mydf2=mydf2,
                mydf3=mydf3,
                mydf4=mydf4,
                depmap_data2=depmap_data2))
    
  }) #End myreactive
  
  
  
  output$plot1 <- renderPlotly({
    
    if(input$stat_metric=="p-value"){
      
      
      p <- ggplot(mysheetreactive()$mydf, aes(x=NB_effectSize, y=-log10(NB_Sensitizing_P), text=Gene)) +
        geom_point(color="gray76") +
        facet_wrap(~tx, ncol=2) +
        xlab("NB Effect Size") +
        ylab("-Log10(p-value)") +
        geom_hline(yintercept=-log10(input$stat_alpha)) +
        geom_vline(xintercept=input$effectSize) +
        theme(panel.border=element_rect(fill=NA, color="black"),
              panel.background=element_rect(fill="gray95")) +
        geom_point(data=myreactive()$mydf2, aes(x=NB_effectSize, y=-log10(NB_Sensitizing_P)), color="red")
      
    } else if(input$stat_metric=="FDR"){
      
      
      p <- ggplot(mysheetreactive()$mydf, aes(x=NB_effectSize, y=-log10(NB_Sensitizing_FDR), text=Gene)) +
        geom_point(color="gray76") +
        facet_wrap(~tx, ncol=2) +
        xlab("NB Effect Size") +
        ylab("-Log10(FDR)") +
        geom_hline(yintercept=-log10(input$stat_alpha)) +
        geom_vline(xintercept=input$effectSize) +
        theme(panel.border=element_rect(fill=NA, color="black"),
              panel.background=element_rect(fill="gray95")) +
        geom_point(data=myreactive()$mydf2, aes(x=NB_effectSize, y=-log10(NB_Sensitizing_FDR)), color="red")
    }
    
    
    ggplotly(p) %>% config(displayModeBar = F)
    
  })
  
  #Count how many total drug combinations there are (rows of mydf4) and unique drugs there are (unique entries of 'known_drug' column)
  output$mytext <- renderText({
    req(myreactive()$mydf4)
    mydf4 <- myreactive()$mydf4
    
    
    myresult <- paste0(HTML("<b><font size='+2'>",nrow(mydf4),"</b></font> total combinations <b>", "<br><br>",
                            length(unique(mydf4$known_drug)), "</b> unique compounds"))
    
    return(myresult)
    
  })
  
  
  #Render data table
  output$mytable = DT::renderDataTable({
    req(myreactive()$mydf4)
    DT::datatable(myreactive()$mydf4, 
                  extensions = c('Responsive', 'Buttons'),
                  options = list(scrollY='700px',
                                 paging=FALSE,
                                 ordering=TRUE,
                                 searching=TRUE,
                                 dom = 'Bfrtip',
                                 buttons = c('csv')))
  })
  
  
  #The popup that appears when the actionlink 'info_button' is clicked
  observeEvent(input$info_button, {
    # Show a modal when the button is pressed
    shinyalert("Drugs per Gene", 
               "The drugs per gene are capped at 10 by default, but this can be changed. When value is >1, drugs are selected to maximize the unique total drug combinations.",
               type = "info")
  })
  
  
  output$plot2 <- renderPlotly({
    req(myreactive()$mydf4)
    
    mydf4 <- myreactive()$mydf4
    
    p <- mydf4 %>% 
      group_by(tx) %>%
      summarize(count=n()) %>%
      ggplot(aes(x=tx, y=count)) + 
      geom_bar(stat="identity", fill="lightblue", color="black")+
      theme(panel.border=element_rect(fill=NA, color="black"),
            axis.text.x = element_text(angle = 45, vjust=0.5),
            plot.title=element_text(hjust=0.5)) +
      xlab("") +
      ggtitle("Total # of drugs per treatment")
    
    ggplotly(p) %>% config(displayModeBar = F)
  })
  
  
  
  
  output$plot3 <- renderPlotly({
    req(myreactive()$mydf4)
    
    mydf4 <- myreactive()$mydf4
    
    p <- mydf4 %>% 
      group_by(Gene) %>%
      summarize(count=n()) %>%
      ggplot(aes(x=Gene, y=count)) + 
      geom_bar(stat="identity", fill="darkolivegreen2", color="black")+
      theme(panel.border=element_rect(fill=NA, color="black"),
            axis.text.x = element_text(angle = 45, vjust=0.5),
            plot.title=element_text(hjust=0.5)) +
      xlab("") +
      ggtitle("Total # of drugs per gene (all treatments)")
    
    ggplotly(p) %>% config(displayModeBar = F)
  })
  
  
  
  
  
  
  
  output$plot4 <- renderPlotly({
    req(myreactive()$depmap_data2)
    
    depmap_data2 <- myreactive()$depmap_data2
    
    p <-ggplot(depmap_data2, aes(x=reorder(Gene, -log10(`P-Value`)), y=-log10(`P-Value`))) +
      geom_point() +
      xlab("Gene") +
      ylab("-Log10(p-value)") +
      theme(panel.border=element_rect(fill=NA, color="black"),
            plot.title = element_text(hjust=0.5))+
      ggtitle("Hits referenced against dependencies in Neuroblastoma (DepMap)")
    
    ggplotly(p) %>% config(displayModeBar = F)
  })
  
  
  
  
  
})

shinyApp(ui, server)