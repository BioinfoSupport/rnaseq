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
					box(title="Parameters",width=12,
							textInput("ige3_link","iGE3 link",placeholder = "https://data.ige3.genomics.unige.ch/dataset/download/xxxxxxxxxx"),
							textInput("ige3_fastq_outdir","Destination subfolder (in data/fastq/)",value = "test")
					),
					box(title="Run",width=12,collapsible = FALSE,
							p(textOutput("ige3_cmd_line",container = code)),
							actionButton("download","Download",icon = icon("download"))
					),
					box(title="Log",icon = icon("terminal"),width=12,
							textOutput("ige3_terminal_stdout",container = pre)
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
							p(textOutput("rnaseq_cmd_line",container=code)),
							actionButton("run_rnaseq","Run",icon=icon("play"))
					),
					box(title="Log",icon = icon("terminal"),width=12,
							textOutput("rnseq_terminal_stdout",container = pre)
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
	ige3_cmd_stdout <- reactiveVal(NULL)
	ige3_cmd <- reactive({
		sprintf(
			"mkdir -p '%s' && wget -P '%s' --content-disposition --trust-server-names -i '%s' 2>&1",
			file.path("data/fastq",input$ige3_fastq_outdir),
			file.path("data/fastq",input$ige3_fastq_outdir),
			sub("(\\.txt)?$",".txt",input$ige3_link)
		)
	})
	output$ige3_cmd_line <- renderText({ige3_cmd()})
	
	# On Action Button Click
	observeEvent(input$download,{ige3_cmd_stdout(pipe(ige3_cmd()))})
	output$ige3_terminal_stdout <- renderText({
		validate(need(!is.null(ige3_cmd_stdout()),"Not run yet"))
		print("bbb")
		paste(readLines(ige3_cmd_stdout()),collapse = "\n")
	})
	
	
	#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
	# RNASEQ panel pipeline logic
	#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
	rnaseq_cmd_stdout <- reactiveVal(NULL)
	rnaseq_cmd <- reactive({sprintf("make GENOME='%s' '%s/all' 2>&1",input$rnaseq_genome,input$rnaseq_fastq_dir)})
	output$rnaseq_cmd_line <- renderText({rnaseq_cmd()})
	
	# On Action Button Click
	observeEvent(input$run_rnaseq,{rnaseq_cmd_stdout(pipe(rnaseq_cmd()))})
	#timer <- reactiveTimer(1000)
	output$rnseq_terminal_stdout <- renderText({
		validate(need(!is.null(rnaseq_cmd_stdout()),"Not run yet"))
		#timer()
		print("aaa")
		paste(readLines(rnaseq_cmd_stdout()),collapse = "\n")
	})
	
}

shinyApp(ui, server)