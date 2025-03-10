library(shiny)
#install.packages("shinythemes")
library(shinythemes)
library(NucleoFeature)

# UI: User Interface
ui <- fluidPage(
  tags$style(HTML("
    .my-custom-button {
      background-color: #4499ff !important; /* Custom orange color */
      color: white !important; /* White text */
      border-radius: 10px; /* Rounded corners */
      font-size: 16px;
      font-weight: bold;
      padding: 10px 20px;
      border: none;
    }
    .my-custom-button:hover {
      background-color: #44aa19 !important; /* Darker color on hover */
    }
  ")),
  #theme = shinythemes::shinytheme("superhero"),
  titlePanel("Nucleotide Sequence Analysis"),

  sidebarLayout(
    sidebarPanel(
      textAreaInput("sequence", "Enter DNA Sequence:",
                    value = "ACGTACGTACGT", rows = 12),
      numericInput("kmer", "Select k-mer size:", value = 2, min = 1, max = 6),
      actionButton("analyze", "Analyze",class = "my-custom-button"),
      hr(),
      fileInput("file", "Upload Sequence File (.txt)", accept = ".txt"),
      downloadButton("downloadData", "Download Results as CSV",class = "my-custom-button")  # <-- Download Button in UI
    ),


    mainPanel(
      h2("Sequence-Based Compositional Features"),

      h4("Nucleotide Composition"),
      tableOutput("composition"),

      h4("GC Content"),
      verbatimTextOutput("gcContent"),

      h4("Shannon Entropy"),
      verbatimTextOutput("entropy"),

      h4("K-mer Frequency"),
      tableOutput("kmerTable"),
      h2("Pseudo K-tuple Nucleotide Composition (PseKNC)"),
      tableOutput("psekncTable"),
      h2("Structural Features"),

      h4("DNA Stability Index"),
      verbatimTextOutput("stabilityIndex"),

      h4("Z-DNA Probability"),
      verbatimTextOutput("zDnaProbability"),
      h4("DNA Curvature"),
      verbatimTextOutput("curvature"),
      h4("Helical Twist"),
      verbatimTextOutput("Helicaltwist"),
      h4("Correlation Features"),
      verbatimTextOutput("correlation"),
      h4("Melting temperature"),
      verbatimTextOutput("Meltingtemperature"),

      h2("Physicochemical Features"),


      h4("Hydrophobicity, stacking energy, bendability, flexibility"),
      verbatimTextOutput("Structure"),

      h4("Nucleosome Positioning Signals"),
      verbatimTextOutput("Nucleosome"),

      h4("DNA Bending Stiffness"),
      verbatimTextOutput("BendingStiffness"),

      h2("Motif-Based Features"),

      h4("Transcription Factor Binding Sites"),
      tableOutput("BindingSites"),


      h4("CpG Islands"),
      verbatimTextOutput("CpGIslands"),

     # h4("DNA Conserved Elements"),
     # verbatimTextOutput("DNAConserved"),
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
    kmerTable <- kmer_frequency(seq_input, input$kmer)
    psekncTable<-pseknc_composition(seq_input, input$kmer)# k-mer frequencies
    stabilityIndex<-compute_stability(seq_input)
    zDnaProbability<-compute_z_dna_probability(seq_input)
    curvature<-curvature_value(seq_input,input$kmer)
    Helicaltwist<-helical_twist(seq_input)
    correlation<-compute_autocorrelation(seq_input)
    Meltingtemperature<-compute_tm(seq_input)
    Structure<-compute_physicochemical_properties(seq_input)
    Nucleosome<-compute_nucleosome_positioning(seq_input)
    BendingStiffness<-compute_bending_stiffness(seq_input)
    BindingSites<-detect_tfbs(seq_input)
    CpGIslands<-detect_cpg_island(seq_input)
    #DNAConserved<-detect_conserved_elements(seq_input)
    return(list(
      composition = comp_result,
      gc_content = gc_result,
      entropy = entropy_result,
      kmerTable = kmerTable,
      psekncTable=psekncTable,
      stabilityIndex= stabilityIndex,
      zDnaProbability=zDnaProbability,
      curvature=curvature,
      Helicaltwist=Helicaltwist,
      correlation=correlation,
      Meltingtemperature= Meltingtemperature,
      Structure=Structure,
      Nucleosome=Nucleosome,
      BendingStiffness=BendingStiffness,
      BindingSites=BindingSites,
      CpGIslands=CpGIslands
      #DNAConserved=DNAConserved
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
  output$psekncTable <- renderTable({
    analyze_sequence()$ psekncTable
  })
  output$stabilityIndex <- renderPrint({
    analyze_sequence()$stabilityIndex
  })

  output$zDnaProbability <- renderPrint({
    analyze_sequence()$zDnaProbability
  })
  output$curvature <- renderPrint({
    analyze_sequence()$curvature
  })
  output$Helicaltwist <- renderPrint({
    analyze_sequence()$Helicaltwist
  })
  output$correlation <- renderPrint({
    analyze_sequence()$correlation
  })
  output$Meltingtemperature <- renderPrint({
    analyze_sequence()$Meltingtemperature
  })
  output$Structure <- renderPrint({
    analyze_sequence()$Structure
  })
  output$Nucleosome <- renderPrint({
    analyze_sequence()$Nucleosome
  })

  output$BendingStiffness <- renderPrint({
    analyze_sequence()$BendingStiffness
  })


  output$BindingSites <- renderTable({
    binding <- analyze_sequence()$BindingSites  # Extract binding sites data

    # Debugging: Print structure
    print(str(binding))

    # Check if it's a named numeric vector
    if (is.null(binding) || length(binding) == 0) {
      return(data.frame(BindingSite = "No data", Frequency = 0))
    }

    if (is.numeric(binding) && !is.data.frame(binding)) {
      binding <- data.frame(BindingSite = names(binding), Frequency = as.numeric(binding))
    }

    # Ensure there are no NA values
    binding[is.na(binding)] <- 0

    return(binding)  # Return properly formatted data frame
  })




    output$CpGIslands<- renderPrint({
      analyze_sequence()$CpGIslands
  })


  # Output: K-mer Frequency Table
    output$kmerTable <- renderTable({
      kmer_data <- analyze_sequence()$kmerTable  # Extract k-mer frequencies

      # Debugging: Check structure
      print(str(analyze_sequence()))
      print(head(analyze_sequence()$kmerTable))


      # Handle NULL or incorrectly formatted data
      if (is.null(kmer_data) || nrow(kmer_data) == 0) {
        return(data.frame(Kmer = "No data", Frequency = 0))
      }

      # Ensure the data frame has proper column names
      if (!("Kmer" %in% colnames(kmer_data)) || !("Frequency" %in% colnames(kmer_data))) {
        return(data.frame(Kmer = "Error", Frequency = "NA"))
      }

      return(kmer_data)  # Directly return data frame
    })




  # Plot: Nucleotide Composition
  output$compositionPlot <- renderPlot({
    barplot(as.numeric(analyze_sequence()$composition),
            names.arg = names(analyze_sequence()$composition),
            col = "skyblue", main = "Nucleotide Composition")
  })

 # output$kmerPlot <- renderPlot({
 #   kmer_data <- analyze_sequence()$kmer_freq


    # Ensure k-mer data is structured properly
 #   kmer_df <- data.frame(Kmer = names(kmer_data), Frequency = as.numeric(kmer_data))

    # Handle missing or zero values
 #   kmer_df$Frequency[is.na(kmer_df$Frequency)] <- 0

    # Prevent empty plot issue
    #if (all(kmer_df$Frequency == 0)) {
  #    plot.new()
  #    text(0.5, 0.5, "All k-mer frequencies are zero!", cex = 1.5, col = "red")
  #    return()
  #  }

    # Plot the k-mer frequencies
 #   barplot(kmer_df$Frequency,
      #      names.arg = kmer_df$Kmer,
      #      col = "lightgreen",
      #      las = 2,  # Rotate labels for better visibility
      #      ylim = c(0, max(kmer_df$Frequency, na.rm = TRUE) + 5),  # Adjust ylim dynamically
      #      main = paste(input$kmer, "-mer Frequency"),
        #    ylab = "Frequency (%)",
       #     xlab = "K-mers")
 # })

  output$downloadData <- downloadHandler(
    filename = function() {
      paste("sequence_analysis_results", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      # Retrieve the analysis results
      results_list <- analyze_sequence()

      # Convert composition (which may be a named vector or data frame) to a data frame
      composition_df <- as.data.frame(t(results_list$composition))  # Transpose for CSV format

      # Create a data frame with other results
      results <- data.frame(
        GC_Content = results_list$gc_content,
        Shannon_Entropy = results_list$entropy,
        Stability_Index = results_list$stabilityIndex,
        Z_DNA_Probability = results_list$zDnaProbability,
        Curvature = results_list$curvature,
        Helical_Twist = results_list$Helicaltwist,
        Melting_Temperature = results_list$Meltingtemperature,
        Hydrophobicity = results_list$Structure,
        Nucleosome_Positioning = results_list$Nucleosome,
        BendingStiffness= results_list$BendingStiffness,
        BindingSites=results_list$BindingSites,
        CpGIslands=results_list$CpGIslands


      )

      # Merge both data frames
      final_results <- cbind(composition_df, results)

      # Write the complete data frame to CSV
      write.csv(final_results, file, row.names = FALSE)
    }
  )
}
# Run the application
shinyApp(ui = ui, server = server)
