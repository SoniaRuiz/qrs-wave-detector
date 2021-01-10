pageWithSidebar(
  headerPanel('Detection of QRS complex and PT waves'),
  sidebarPanel(
   
    #numericInput('numSamples', 'Num. muestras', 500,
     #            min = 500, max = 5000)
    fileInput("fichero", "Open CSV file", multiple = FALSE, 
        accept = c(
          "text/csv",
          "text/comma-separated-values,text/plain",
          ".csv")#,
        #buttonLabel = "Abrir...",
        #placeholder = "Por favor, seleccione un fichero .csv",
          
        ),
    wellPanel(
      sliderInput("samples", "Number of samples", 2000, 5000,
                  value = 2000, step = 500)
    )
  ),
  mainPanel(
    plotOutput('plot1')
  )
)