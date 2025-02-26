  # Load your package

# UI: User Interface
# Load necessary libraries
library(shiny)
library(NucleoFeature)  # Load your custom packageF

# UI: User Interface
ui <- fluidPage(
  titlePanel("Nucleotide Sequence Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      textAreaInput("sequence", "Enter DNA Sequence:", 
                    value = "ACGTACGTACGT", rows = 3),
      numericInput("kmer", "Select k-mer size:", value = 2, min = 1, max = 6),
      actionButton("analyze", "Analyze"),
      hr(),
      fileInput("file", "Upload Sequence File (.txt)", accept = ".txt")
    ),
    
    mainPanel(
      h4("Nucleotide Composition"),
      tableOutput("composition"),
      
      h4("GC Content"),
      verbatimTextOutput("gcContent"),
      
      h4("Shannon Entropy"),
      verbatimTextOutput("entropy"),
      
      h4("K-mer Frequency"),
      tableOutput("kmerTable"),
      
      plotOutput("compositionPlot"),
      plotOutput("kmerPlot")
    )
  )
)

# Server: Backend Processing
server <- function(input, output) {
  
  analyze_sequence <- reactive({
    req(input$analyze)  # Ensure button is clicked
    seq_input <- input$sequence  # Get user input
    
    if (!is.null(input$file)) {  
      # Read file if uploaded
      seq_input <- readLines(input$file$datapath)
      seq_input <- paste(seq_input, collapse = "")
    }
    
    # Compute features using your package functions
    comp_result <- ACTG_composition(seq_input)  # Nucleotide composition
    gc_result <- gc_content(seq_input)  # GC content
    entropy_result <- shannon_entropy(seq_input)  # Sequence entropy
    kmer_result <- frequency(seq_input, input$kmer)  # k-mer frequencies
    
    return(list(
      composition = comp_result,
      gc_content = gc_result,
      entropy = entropy_result,
      kmer_freq = kmer_result
    ))
  })
  
  # Output: Nucleotide Composition Table
  output$composition <- renderTable({
    analyze_sequence()$composition
  })
  
  # Output: GC Content
  output$gcContent <- renderPrint({
    analyze_sequence()$gc_content
  })
  
  # Output: Shannon Entropy
  output$entropy <- renderPrint({
    analyze_sequence()$entropy
  })
  
  # Output: K-mer Frequency Table
  output$kmerTable <- renderTable({
    kmer_data <- analyze_sequence()$kmer_freq
    
    if (is.null(kmer_data)) {
      return(data.frame(Kmer = NA, Frequency = NA))  # Avoid NULL output
    }
    
    # Convert named vector to a data frame
    kmer_df <- data.frame(Kmer = names(kmer_data), Frequency = as.numeric(kmer_data))
    
    return(kmer_df)  # Return as a data frame for proper table display
  })
  
  
  # Plot: Nucleotide Composition
  output$compositionPlot <- renderPlot({
    barplot(as.numeric(analyze_sequence()$composition), 
            names.arg = names(analyze_sequence()$composition), 
            col = "skyblue", main = "Nucleotide Composition")
  })
  
  output$kmerPlot <- renderPlot({
    kmer_data <- analyze_sequence()$kmer_freq
    
    # Check if kmer_data is valid
    if (is.null(kmer_data) || length(kmer_data) == 0) {
      plot.new()
      text(0.5, 0.5, "No valid k-mers found!", cex = 1.5, col = "red")
      return()
    }
    
    # Ensure k-mer data is structured properly
    kmer_df <- data.frame(Kmer = names(kmer_data), Frequency = as.numeric(kmer_data))
    
    # Handle missing or zero values
    kmer_df$Frequency[is.na(kmer_df$Frequency)] <- 0
    
    # Prevent empty plot issue
    if (all(kmer_df$Frequency == 0)) {
      plot.new()
      text(0.5, 0.5, "All k-mer frequencies are zero!", cex = 1.5, col = "red")
      return()
    }
    
    # Plot the k-mer frequencies
    barplot(kmer_df$Frequency, 
            names.arg = kmer_df$Kmer, 
            col = "lightgreen", 
            las = 2,  # Rotate labels for better visibility
            ylim = c(0, max(kmer_df$Frequency, na.rm = TRUE) + 5),  # Adjust ylim dynamically
            main = paste(input$kmer, "-mer Frequency"),
            ylab = "Frequency (%)",
            xlab = "K-mers")
  })
  
}

# Run the application
shinyApp(ui = ui, server = server)
#shiny::runApp("app.R")


