library(bs4Dash)


ui <- dashboardPage(dark = NULL,help = NULL,
	header = dashboardHeader(title = "NGS"),
	sidebar = dashboardSidebar(
		sidebarMenu(
			sidebarHeader("Menu"),
			menuItem(
				"FASTQ Download from iGE3",
				tabName = "iGE3",
				icon = icon("download")
			),
			menuItem(
				"RNAseq Quantification",
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
					box(title="Parameters",
							textInput("ige3_link","iGE3 link"),
							textInput("ige3_fastq_outdir","Destination subfolder (in data/fastq/)",value = "project1")
					),
					box(title="Run",width=12,collapsible = FALSE,
							p(textOutput("ige3_run_command",container = code)),
							actionButton("download","Download",icon = icon("download"))
					)
				)
			),
			tabItem(
				tabName = "rnaseq",
				fluidRow(
					box(title = "Parameters",
							selectInput("rnaseq_genome","GENOME",choices = c("Dd+Mm","Dd","Mm","GRCh38-r45","GRCm39-M34")),
							textInput("rnaseq_fastq_dir","FASTQ directory",value = "data/fastq/")
					),
					box(title="Run",width=12,collapsible = FALSE,
							p(textOutput("rnaseq_run_command",container=code)),
							actionButton("run_rnaseq","Run",icon=icon("play"))
					)
				)				
			)
		)
	)
)

server <- function(input, output, session) {
	observe({
		updateTextInput(session,"rnaseq_fastq_dir",value = sprintf("data/fastq/%s",input$ige3_fastq_outdir))
	})
	
	output$ige3_run_command <- renderText({
		sprintf(
			"mkdir -p '%s' && wget -P '%s' --content-disposition --trust-server-names -i '%s'",
			file.path("data/fastq",input$ige3_fastq_outdir),
			file.path("data/fastq",input$ige3_fastq_outdir),
			input$ige3_link
		)
	})
	output$rnaseq_run_command <- renderText({
		sprintf(
			"make GENOME='%s' '%s/all'",
			input$rnaseq_genome,
			input$rnaseq_fastq_dir
		)
	})	
}

shinyApp(ui, server)