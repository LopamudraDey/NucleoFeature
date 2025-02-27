library(shiny)
library(NucleoFeature)
# UI: User Interface
ui <- fluidPage(
  titlePanel("Nucleotide Sequence Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      textAreaInput("sequence", "Enter DNA Sequence:", 
                    value = "ACGTACGTACGT", rows = 12),
      numericInput("kmer", "Select k-mer size:", value = 2, min = 1, max = 6),
      actionButton("analyze", "Analyze"),
      hr(),
      fileInput("file", "Upload Sequence File (.txt)", accept = ".txt")
    ),
    
    mainPanel(
      h2("Sequence-Based Features"),
      
      h4("Nucleotide Composition"),
      tableOutput("composition"),
      
      h4("GC Content"),
      verbatimTextOutput("gcContent"),
      
      h4("Shannon Entropy"),
      verbatimTextOutput("entropy"),
      
      h4("K-mer Frequency"),
      tableOutput("kmerTable"),
      h2("Structural Features"),
      
      h4("DNA Stability Index"),
      verbatimTextOutput("stabilityIndex"),
      
      h4("Z-DNA Probability"),
      verbatimTextOutput("zDnaProbability"),
      h4("Melting temperature"),
      verbatimTextOutput("Meltingtemperature"),
      
      h2("Physicochemical Features"),
      
      
      h4("Hydrophobicity, stacking energy, bendability, flexibility"),
      verbatimTextOutput("Hydrophobicity"),
      
      h4("Nucleosome Positioning Signals"),
      verbatimTextOutput("Nucleosome"),
      
      h4("DNA Bending Stiffness"),
      verbatimTextOutput("DNA Bending Stiffness"),
      
      h2("Motif-Based Features"),
      
      h4("Transcription Factor Binding Sites"),
      tableOutput("Transcription Factor Binding Sites"),
      
      h4("CpG Islands"),
      verbatimTextOutput("CpG Islands"),
       
      h4("DNA Conserved Elements"),
      verbatimTextOutput("DNA Conserved Elements"),
      plotOutput("compositionPlot"),
      plotOutput("kmerPlot")
      
    
        
          )
  )
)

# Empty server function
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
    stabilityIndex<-compute_stability(seq_input)
    zDnaProbability<-compute_z_dna_probability(seq_input)
    Meltingtemperature<-compute_tm(seq_input)
    Hydrophobicity<-compute_physicochemical_properties(seq_input)
    Nucleosome<-compute_nucleosome_positioning(seq_input)
    return(list(
      composition = comp_result,
      gc_content = gc_result,
      entropy = entropy_result,
      kmer_freq = kmer_result,
      stabilityIndex= stabilityIndex,
      zDnaProbability=zDnaProbability,
      Meltingtemperature= Meltingtemperature,
      Hydrophobicity=Hydrophobicity,
      Nucleosome=Nucleosome
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
  output$stabilityIndex <- renderPrint({
    analyze_sequence()$stabilityIndex
  })
  
  output$zDnaProbability <- renderPrint({
    analyze_sequence()$zDnaProbability
  })
  
  output$Meltingtemperature <- renderPrint({
    analyze_sequence()$Meltingtemperature
  }) 
  output$Hydrophobicity <- renderPrint({
    analyze_sequence()$Hydrophobicity
  }) 
  output$Nucleosome <- renderPrint({
    analyze_sequence()$Nucleosome
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
