library(bs4Dash)


ui <- dashboardPage(dark = NULL,help = NULL,
	header = dashboardHeader(title = "NGS"),
	sidebar = dashboardSidebar(
		sidebarMenu(
			sidebarHeader("Menu"),
			menuItem(
				"FASTQ download",
				tabName = "iGE3",
				icon = icon("download")
			),
			menuItem(
				"Read quantification",
				tabName = "rnaseq",
				icon = icon("calculator")
			)
		)
	),
	body = dashboardBody(
		tabItems(
			tabItem(
				tabName = "iGE3",
				fluidRow(
					box(title="Parameters",width=12,collapsible = FALSE,
							textInput("ige3_link","iGE3 link",placeholder = "https://data.ige3.genomics.unige.ch/dataset/download/xxxxxxxxxx"),
							textInput("ige3_fastq_outdir","Destination subfolder (in data/fastq/)",value = "test"),
							actionButton("download","Download",icon = icon("download")),
							verbatimTextOutput("ige3_terminal_stdout")
					)
				)
			),
			tabItem(
				tabName = "rnaseq",
				fluidRow(
					box(title="Run",width=12,collapsible = FALSE,
							selectInput("rnaseq_genome","GENOME",choices = c("Dd+Mm","Dd","Mm","GRCh38-r45","GRCm39-M34")),
							textInput("rnaseq_fastq_dir","FASTQ directory",value = "data/fastq/"),
							actionButton("run_rnaseq","Run",icon=icon("play")),
							verbatimTextOutput("rnseq_terminal_stdout")
					)
				)				
			)
		)
	)
)

server <- function(input, output, session) {
	
	# Update rnaseq_fastq_dir when changing ige3 download dir
	observe({
		updateTextInput(session,"rnaseq_fastq_dir",value = sprintf("data/fastq/%s",input$ige3_fastq_outdir))
	})
	
	
	#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
	# IGE3 panel logic
	#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
	output$ige3_terminal_stdout <- eventReactive(input$download,{
		cmd <- sprintf(
			"mkdir -p '%s' && wget -P '%s' --content-disposition --trust-server-names -i '%s' 2>&1",
			file.path("data/fastq",input$ige3_fastq_outdir),
			file.path("data/fastq",input$ige3_fastq_outdir),
			sub("(\\.txt)?$",".txt",input$ige3_link)
		)
		print(cmd)
		withProgress(message="Downloading the data",{
			pipe(cmd) |> readLines() |> paste(collapse = "\n")
		})
	})
	
	
	#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
	# RNASEQ panel pipeline logic
	#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
	output$rnseq_terminal_stdout <- eventReactive(input$run_rnaseq,{
		cmd <- sprintf("make GENOME='%s' '%s/RNASEQ.ALL' 2>&1",input$rnaseq_genome,input$rnaseq_fastq_dir)
		print(cmd)
		withProgress(message="Running the pipeline",{
			pipe(cmd) |> readLines() |> paste(collapse = "\n")
		})
	})
}

shinyApp(ui, server)