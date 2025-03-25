library(shiny)
#install.packages("shinythemes")
library(shinythemes)
library(NucleoFeature)

# UI: User Interface
ui <- fluidPage(
  tags$style(HTML("
    .my-custom-button {
      background-color: #4499ff;
      color: white;
      border-radius: 10px;
      font-size: 16px;
      font-weight: bold;
      padding: 10px 20px;
      border: none;
    }
    .my-custom-button:hover {
      background-color: #44aa19 ;
    }
  ")),
 # theme = shinythemes::shinytheme("superhero"),
  titlePanel("NucleoFeature: Shiny Interface for Nucleotide Sequence Analysis"),

  sidebarLayout(
    sidebarPanel(
      textAreaInput("sequence", "Enter DNA Sequence:",
                    value = "TACGTACGTACGT", rows = 12),
      numericInput("kmer", "Select k-mer size:", value = 2, min = 1, max = 6),
      actionButton("analyze", "Analyze",class = "my-custom-button"),
      hr(),
      fileInput("file", "Upload Sequence File (.txt)", accept = ".txt"),
      downloadButton("downloadData", "Download Results as CSV",class = "my-custom-button")    ),


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
      h2("Reverse Complement and Kmer"),

      verbatimTextOutput("reverse_complement_output"),

      tableOutput("kmer_compositionreverse"),
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


      tableOutput("BindingSites"),


      h4("CpG Islands"),
      verbatimTextOutput("CpGIslands"),

     # h4("DNA Conserved Elements"),
     # verbatimTextOutput("DNAConserved"),
      plotOutput("compositionPlot"),
      plotOutput("kmerPlot")



          )

  ),
 HTML("<hr><p style='text-align:center; font-weight:bold;'>
       Developed by Lopamudra Dey | © 2025</p>")
)

server <- function(input, output) {

  analyze_sequence <- reactive({
    req(input$analyze)  # Ensure button is clicked
    seq_input <- input$sequence  # Get user input

    if (!is.null(input$file)) {
      seq_input <- readLines(input$file$datapath)
      seq_input <- paste(seq_input, collapse = "")
    }

    comp_result <- ACTG_composition(seq_input)  
    gc_result <- gc_content(seq_input)  
    entropy_result <- shannon_entropy(seq_input)  
    kmerTable <- kmer_frequency(seq_input, input$kmer)
    reverse_complement_output<-reverse_seqcomplement(seq_input)
    kmer_compositionreverse<-kmer_compositionreverse(seq_input, input$kmer)
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
   # BindingSites<-detect_tfbs(seq_input)
    CpGIslands<-detect_cpg_island(seq_input)
    #DNAConserved<-detect_conserved_elements(seq_input)
    return(list(
      composition = comp_result,
      gc_content = gc_result,
      entropy = entropy_result,
      kmerTable = kmerTable,
      psekncTable=psekncTable,
      reverse_complement_output=reverse_complement_output,
      kmer_compositionreverse=kmer_compositionreverse,
      stabilityIndex= stabilityIndex,
      zDnaProbability=zDnaProbability,
      curvature=curvature,
      Helicaltwist=Helicaltwist,
      correlation=correlation,
      Meltingtemperature= Meltingtemperature,
      Structure=Structure,
      Nucleosome=Nucleosome,
      BendingStiffness=BendingStiffness,
      #BindingSites=BindingSites,
      CpGIslands=CpGIslands
      #DNAConserved=DNAConserved
    ))
  })

  output$composition <- renderTable({
    analyze_sequence()$composition
  })

  output$gcContent <- renderPrint({
    analyze_sequence()$gc_content
  })

  output$entropy <- renderPrint({
    analyze_sequence()$entropy
  })
  output$psekncTable <- renderTable({
    analyze_sequence()$ psekncTable
  })
  output$reverse_complement_output <- renderText({
    analyze_sequence()$reverse_complement_output
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

  output$kmer_compositionreverse <- renderTable({
    analyze_sequence()$kmer_compositionreverse
  })

    output$CpGIslands<- renderPrint({
      analyze_sequence()$CpGIslands
  })


    output$kmerTable <- renderTable({
      kmer_data <- analyze_sequence()$kmerTable  

      print(str(analyze_sequence()))
      print(head(analyze_sequence()$kmerTable))


      if (is.null(kmer_data) || nrow(kmer_data) == 0) {
        return(data.frame(Kmer = "No data", Frequency = 0))
      }

      if (!("Kmer" %in% colnames(kmer_data)) || !("Frequency" %in% colnames(kmer_data))) {
        return(data.frame(Kmer = "Error", Frequency = "NA"))
      }

      return(kmer_data)
    })




  # Plot: Nucleotide Composition
  output$compositionPlot <- renderPlot({
    barplot(as.numeric(analyze_sequence()$composition),
            names.arg = names(analyze_sequence()$composition),
            col = "skyblue", main = "Nucleotide Composition")
  })


 # })

  output$downloadData <- downloadHandler(
    filename = function() {
      paste("sequence_analysis_results", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {

      results_list <- analyze_sequence()
      composition_df <- as.data.frame(t(results_list$composition))

     
      reverse_complement_str <- ifelse(is.null(results_list$reverse_complement_output),
                                       NA,
                                       paste(results_list$reverse_complement_output, collapse = "; "))
      #df2<- as.data.frame(t(results_list$kmer_compositionreverse))
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
        CpGIslands=results_list$CpGIslands,
       reverse_complement_output=reverse_complement_str

      )

      final_results <- cbind(composition_df, results)
      str(composition_df)

      write.csv(final_results, file, row.names = FALSE)
    }
  )
}

shinyApp(ui = ui, server = server)
