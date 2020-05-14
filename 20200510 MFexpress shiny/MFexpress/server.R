
#####################################
##                                 ##
##   This shiny app is for gene    ##
## expression quick search of MFs  ##
##                                 ##
#####################################

#
# This is a online shiny app version. For detail information, please see in lab R project version.
# Dataset includes tsm_gene (default, tissue resident macropahge, ImmGen)
#                  m1m2_gene (M1, M2 data from Abhishek K. Jha, et al. Immunity(2015))
#                  lact_gene (LPS treated macrophage from Di Zhang, et al. Nature(2019)).
#
# last updated : 2020.05.10 Developed by Tsai,Tsung-Lin (ttsunglin@gmail.com)
#                2020.05.13 Device detection added (mobile detection method from https://github.com/g3rv4/mobileDetect)
#

#----library----
suppressMessages(library(readr))
suppressMessages(library(rvest))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))

#----input----

tsm_gene <- read.csv("tsm_gene.csv")
m1m2_gene <- read.csv("m1m2_gene.csv")
lact_gene <- read.csv("lact_gene.csv")

colnames(tsm_gene)[2] <- "Entrez Gene ID"
colnames(m1m2_gene)[1] <- "Entrez Gene ID"
colnames(lact_gene)[1] <- "Entrez Gene ID"

#----grouping-----

tsm_group <- c(rep("THY",4),
               rep("PC",4),
               rep("KC",3),
               rep("ALV",3),
               rep("SPL",2))

m1m2_group <- c(rep("M0",2),
                rep("M2",3),
                rep("M1",3))

lact_group <- c(rep("0h",4),
                rep("4h",4),
                rep("8h",4),
                rep("16h",4),
                rep("24h",4))

tsm_order <- c("THY",
               "PC",
               "KC",
               "ALV",
               "SPL")

m1m2_order <- c("M0",
                "M1",
                "M2")


lact_order <- c("0h",
                "4h",
                "8h",
                "16h",
                "24h")


#----entrez id search-----

ENTREZCONVERT <- function(x){
  
  if(class(x) == "character"){
    
    url <- paste("https://www.ncbi.nlm.nih.gov/gene/?term=",
                 x,
                 "+and+mus+musculus",
                 sep = "")
    doc <- read_html(url)           #node of ID in interest, found with SelectorGaget in Chrome
    
    node <- '//*[contains(concat( " ", @class, " " ), concat( " ", "ncbi-doc-id", " " ))]//li'
    html_node <- html_nodes(doc, xpath = node)      
    html_catch <- html_text(html_node)
    
    node_n <- '//*[contains(concat( " ", @class, " " ), concat( " ", "highlight", " " ))]'
    html_node_n <- html_nodes(doc, xpath = node_n)      
    html_catch_n <- html_text(html_node_n)
    
    if(identical(html_catch, character(0)) == TRUE){
      
      node <- '//*[contains(concat( " ", @class, " " ), concat( " ", "ncbi-doc-id", " " ))]//li'
      html_node <- html_nodes(doc, xpath = node)      
      html_catch <- html_text(html_node)
      
      node_n <- '//*[contains(concat( " ", @class, " " ), concat( " ", "rprt", " " )) and (((count(preceding-sibling::*) + 1) = 1) and parent::*)]//a'
      html_node_n <- html_nodes(doc, xpath = node_n)      
      html_catch_n <- html_text(html_node_n)
    }
    
    if(identical(html_catch, character(0)) == TRUE){
      
      node <- '//*[contains(concat( " ", @class, " " ), concat( " ", "geneid", " " ))]'
      html_node <- html_nodes(doc, xpath = node)      
      html_catch <- html_text(html_node)
      
      node_n <- '//*[contains(concat( " ", @class, " " ), concat( " ", "gn", " " ))]'
      html_node_n <- html_nodes(doc, xpath = node_n)      
      html_catch_n <- html_text(html_node_n)}
    
    if(is.na(html_catch_n[1]) == FALSE) {
      if(as.character(html_catch_n[1]) == "Mus musculus"){
        
        node_n <- '//*[contains(concat( " ", @class, " " ), concat( " ", "rprt", " " )) and (((count(preceding-sibling::*) + 1) = 1) and parent::*)]//a'
        html_node_n <- html_nodes(doc, xpath = node_n)      
        html_catch_n <- html_text(html_node_n)
        
      }
    }
    
    if(identical(html_catch, character(0)) == TRUE){
      
      out_id <- NULL
      out_name <- NULL
      message(paste("input error", "input gene symbol unmatched", sep=": "))
      
    }else{
      
      id_start <- gregexpr(":",html_catch)[[1]][1]+2
      id_end <- gregexpr(",",html_catch)[[1]][1]-1
      
      
      if(id_end < 0){
        
        id_end <- 100
        
      } 
      
      out_id <- substring(html_catch, id_start, id_end)
      out_name <- html_catch_n[1]   
      
    }
    
  }else if(class(x) == "numeric"){
    
    url <- paste("https://www.ncbi.nlm.nih.gov/gene?Db=gene&Cmd=DetailsSearch&Term=",
                 x,
                 sep = "")
    doc <- read_html(url)           #node of ID in interest, found with SelectorGaget in Chrome
    
    node_n <- '//*[contains(concat( " ", @class, " " ), concat( " ", "gn", " " ))]'
    html_node_n <- html_nodes(doc, xpath = node_n)      #output name 
    html_catch_n <- html_text(html_node_n)
    
    if(identical(html_catch_n, character(0)) == TRUE){
      
      out_id <- NULL
      out_name <- NULL
      message(paste("input error", "input Entrez ID unmatched", sep=": "))
      
    }else{
      
      out_name <- html_catch_n
      out_id <- x
      
    }
    
  }else{warning(paste("input error", "input type undefined", sep=": "))}
  
  out <- data.frame(gene = as.character(out_name), id = as.numeric(out_id))
  out
  
}

#----drawing plot----

DOTPLOT <- function(input, dataset = tsm_gene){
  
  if(dataset[3,1] == tsm_gene[3,1]){
    if(class(input) == "character"){
      
      input <- paste(toupper(substring(input,1,1)),substring(input,2),sep = "")
      
    }else if(class(input) == "numeric"){
      
      input <- ENTREZCONVERT(input)$gene %>% as.character()
      
    }
    
    mx_entrez <- dataset[which(dataset$`gene`== input), -c(1:2)] %>% as.matrix()
    mx_entrez_t_df <- data.frame(t(mx_entrez))
    
  }else{
    
    mx_entrez <- dataset[which(dataset$`Entrez Gene ID`== input), -c(1:2)] %>% as.matrix()
    mx_entrez_t_df <- data.frame(t(mx_entrez))
    
  }
  
  if(dataset[3,1] == tsm_gene[3,1]){
    
    name <- "Tissue resident MF"
    level_order <- tsm_order
    group <- tsm_group
    unit <- "RPKM"
    
  }else if(dataset[3,1] == m1m2_gene[3,1]){
    
    name <- "M0, M1, M2"
    level_order <- m1m2_order
    group <- m1m2_group
    unit <- "Log2_QNorm"
    
  }else if(dataset[3,1] == lact_gene[3,1]){
    
    name <- "after LPS treatment"
    level_order <- lact_order
    group <- lact_group
    unit <- "DESeq2"
    
  }else{ warning("dataset error") }
  
  if(ncol(mx_entrez_t_df) < 1 ){
    
    NULL
    
  }else{
    mx_entrez_t_df$type <- group
    colnames(mx_entrez_t_df) <- c("value","type")
    
    p <- ggplot(data = mx_entrez_t_df, aes(x = factor(type, level = level_order), y = value)) +
      geom_boxplot() +
      theme_classic() +
      theme(text = element_text(size=15),plot.title = element_text(size=15)) +
      expand_limits(y = 0) +
      geom_dotplot(binaxis='y', stackdir='center')+
      labs(title = name,
           y = unit,
           x = NULL)+
      guides(size = FALSE)
    p
  }
}

#----shiny server----

library(shiny)

shinyServer(
  function(input, output, session) {
    
    gene_input <- eventReactive(input$enter, {
      
      input$Gene
      
    }
    )
    
    output$myPlot <- renderPlot({
      
      message("------------------------------------------")
      
      suppressWarnings(
        if(is.na(as.numeric(gene_input())) == FALSE){
          
          gene_input <- as.numeric(gene_input())
          
        }else{
          
          gene_input <- gene_input()
          
        }
      )
      withProgress(message = 'Plotting', value = 0, {
        
        time = proc.time()
        
        #plot start
        
        incProgress(1/3, detail = paste("Converting"))
        
        message(paste("input", gene_input(), sep = ": "))
        
        if(class(gene_input) == "numeric"){
          
          conv <- ENTREZCONVERT(gene_input)
          name <- conv$gene %>% as.character()
          id <- conv$id
          message(paste("input type", "Entrez ID", sep = ": "))
          
          if(identical(id, numeric(0)) == TRUE){
            
            name <- paste("input error", "undefined Entrez ID", sep=": \n")
          }
          
        }else{
          
          conv <- ENTREZCONVERT(gene_input)
          name <- conv$gene %>% as.character()
          id <- conv$id
          message(paste("input type", "gene symbol", sep = ": "))
          
          if(identical(id, numeric(0)) == TRUE){
            
            name <- paste("input error", "undefined gene symbol", sep=": \n")
            
          }
          
        } 
        
        incProgress(2/3, detail = paste("Plotting"))
        
        p1 <- DOTPLOT(name, tsm_gene)
        p2 <- DOTPLOT(id, m1m2_gene)
        p3 <- DOTPLOT(id, lact_gene)
        gene_name <- name
        
        if(input$isMobile == TRUE){
          
          ncol <- 1
          nrow <- 3
          message(paste("device", "mobile", sep = ": "))
          
        }else{
          
          ncol <- 3
          nrow <- 1 
          message(paste("device", "PC", sep = ": "))
          
        }
        
        p <- suppressMessages(ggarrange(p1, p2, p3,
                                        labels = c("A", "B", "C"),
                                        ncol = ncol, nrow = nrow))
        
        plot <- annotate_figure(p, 
                                top = text_grob(paste(id,
                                                      gene_name,
                                                      sep = "  "),
                                                color = "red", 
                                                face = "bold", 
                                                size = 25)) 
        #plot end 
        
        print(proc.time()-time)
        
        incProgress(3/3, detail = paste("Express !"))
        
        message("------------------------------------------")
        
        plot
      })
    }
    )
  }
)