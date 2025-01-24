# A web app to process the MacDonald lab growth assay XLSX files
# Alastair Droop & Sandy MacDonald, 2023-10-25

version_id <- "1.0.3"

# Get the necessary libraries:
library(shiny)
library(shinyjs)
library(growthcurver)
library(tidyverse)
library(readxl)
library(DT)
library(zip)

# Define the Application UI
ui <- fluidPage(
  shinyjs::useShinyjs(),
  # Application title:
  titlePanel("Macdonald Lab Growth Rate Data"),
  p("Code available at: ",
    a("https://github.com/uoy-research/cm_growthcurves", href = "https://github.com/uoy-research/cm_growthcurves", target = "_blank")
  ),
  # Sidebar with options:
  sidebarLayout(
    fluid = TRUE,
    sidebarPanel(
      # The input file to process:
      fileInput(
        "input_xlsx",
        "Input XLSX growth data file",
        multiple = FALSE,
        accept = c("application/vnd.openxmlformatsofficedocument.spreadsheetml.sheet")
      ),
      tags$hr(),
      # The output dataset prefix:
      textInput(
        "output_prefix",
        "Output data prefix",
        "",
        placeholder = "output_prefix"
      ),
      # Include full data in plots:
      checkboxInput(
        "show_full_sample_info",
        "Show full sample info in well plot",
        value = FALSE
      ),
      # The reference strain to analyse against prefix:
      textInput(
        "reference_strain",
        "Reference strain",
        "WT",
        placeholder = "WT"
      ),
      # The number of lines to skip:
      numericInput("skip", "Number of header rows to skip", 10, min=0),
      tags$hr(),
      # The run analysis button:
      downloadButton(
        "generate_download",
        "Download Analysis Dataset",
        icon = icon("download", class="fa-solid")
      ),
      # The version number:
      tags$hr(),
      tags$p(sprintf("version %s", version_id))
    ),
    mainPanel(
      uiOutput(
        "startup_message",
        container = tags$div,
        fill = TRUE,
        style = "text-align:center; font-weight: bold"
      ),
      plotOutput("well_plot"),
    )
  )
)

# A function to read times in the form "X h Y min" to minutes:
to_minutes <- function(x) {
  n_hours <- as.numeric(
    gsub("(\\d+) h(?: (\\d+) min)?", "\\1", x, perl = TRUE)
  )
  n_minutes <- as.numeric(
    gsub("(\\d+) h(?: (\\d+) min)?", "\\2", x, perl = TRUE)
  )
  n_minutes[is.na(n_minutes)] <- 0
  return((n_hours * 60) + n_minutes)
}

# Define a logistic curve:
logistic_curve <- function(t, n0, k, r) {
  return((n0 * k) / (n0 + (k - n0) * exp(-r * t)))
}

# Load the growth data from XLSX:
read_growth_xlsx <- function(filename, n_skip, show_full_sample_info) {
  # Read the times directly from XLSX:
  # NB: This skips the 4 empty columns at the start of the time row.
  suppressMessages({
    time_raw <- read_excel(
      path = filename,
      sheet = "All Cycles",
      col_names = FALSE,
      skip = n_skip,
      n_max = 1
    )
    times <- as.character(time_raw[1, ])
    if (times[1] == "Time") {
      times <- times[-1]
    }
    times <- to_minutes(times)
  })
  # Read the raw data directly from XLSX:
  # NB: This is both the metadata & OD data
  suppressMessages({
    data_raw <- read_excel(
      path = filename,
      sheet = "All Cycles",
      col_names = FALSE,
      skip = n_skip + 1,
      na = "blank"
    )
  })
  # Make sure we have 96 wells:
  stopifnot("incorrect well count" = nrow(data_raw) == 96)
  # Split metadata away and tidy it up:
  metadata <- as.data.frame(data_raw[, 1:4])
  colnames(metadata) <- c("row", "col", "treatment", "strain")
  metadata$row <- factor(trimws(metadata$row), ordered=TRUE, levels=LETTERS[1:8])
  metadata$col <- factor(trimws(metadata$col), ordered=TRUE, levels=1:12)
  metadata$treatment <- factor(trimws(metadata$treatment))
  metadata$strain <- factor(trimws(metadata$strain))
  # Make sure the rows are in order:
  row_order <- order(metadata$row, metadata$col)
  metadata <- metadata[row_order, ]
  # Build the sample IDs:
  rownames(metadata) <- sprintf("%s%02d", metadata$row, as.numeric(as.character(metadata$col)))
  # Pull out the OD data:
  od <- data.frame(
    times,
    t(data_raw[, 5:ncol(data_raw)])
  )
  colnames(od) <- c("time", rownames(metadata))
  rownames(od) <- NULL
  # Normalise the OD data:
  min_od <- apply(od, 2, min)
  norm_data <- sweep(od, 2, min_od)
  # Generate the growth curves:
  growth_data <- SummarizeGrowthByPlate(od)
  rownames(growth_data) <- growth_data$sample
  growth_data$fit_ok <- factor(
    c("FAILED", "OK")[as.numeric(growth_data$note == "") + 1],
    levels = c("OK", "FAILED")
  )
  # Build a comprehensive ggplot dataset:
  plot_data <- pivot_longer(
    norm_data,
    cols = !time,
    names_to = "sample",
    values_to = "OD"
  )
  plot_data$prediction <- logistic_curve(
    t = plot_data$time,
    n0 = growth_data[plot_data$sample, "n0"],
    k = growth_data[plot_data$sample, "k"],
    r = growth_data[plot_data$sample, "r"]
  )
  plot_data$row <- metadata[plot_data$sample, "row"]
  plot_data$col <- metadata[plot_data$sample, "col"]
  plot_data$strain <- metadata[plot_data$sample, "strain"]
  plot_data$treatment <- metadata[plot_data$sample, "treatment"]
  
  # Get the maximum OD, rounded to the greatest 1.5:
  max_y <- ceiling(max(plot_data$OD) / 1.5) * 1.5
  
  # Build the label data:
  label_data <- data.frame(
    sample = rownames(metadata),
    row = metadata$row,
    col = metadata$col,
    strain = metadata$strain,
    treatment = metadata$treatment,
    x = 0,
    y = max_y,
    status = growth_data[rownames(metadata), "fit_ok"]
  )
  if (identical(show_full_sample_info, TRUE)) {
    label_data$label <- sprintf("%s (%s/%s)", label_data$sample, label_data$strain, label_data$treatment)
  } else {
    label_data$label <- label_data$sample
  }
  label_data[is.na(label_data$strain), "label"] <- "empty"
  # Plot the 96-well growth data:
  metadata <- metadata[complete.cases(metadata),]
  # suppressWarnings({
  g <- ggplot(plot_data, aes(x = time, y = OD))
  g <- g + geom_point(
    size = 0.25,
    colour = "grey50"
  )
  g <- g + geom_line(
    aes(x = time, y = prediction),
    colour = "red"
  )
  g <- g + ylim(0, max_y)
  g <- g + geom_text(
    aes(x = x, y = y, label = label),
    data = label_data,
    vjust = 1,
    hjust = 0,
    size = 3,
    colour = "navy"
  )
  g <- g + geom_text(
    aes(y = y, label = status),
    x = max(plot_data$time),
    data = label_data,
    vjust = 1,
    hjust = 1,
    size = 3
  )
  g <- g + facet_grid(
    row ~ col,
    switch = "y"
  )
  g <- g + theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 16, face = "bold"),
    strip.background = element_rect(fill = "grey75", colour = "black"),
    strip.text.y.left = element_text(angle = 0),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_rect(fill = "white", colour = "black")
  )
  
  # Return the raw data:
  return(list(
    "metadata" = metadata,
    "OD" = od,
    "norm_data" = norm_data,
    "plot_data" = plot_data,
    "growth_data" = growth_data,
    "growth_data_plot" = g
  ))
}

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  shinyjs::disable("generate_download")
  shinyjs::disable("well_plot")
  shinyjs::disable("output_prefix")
  output$startup_message <- renderText("Load an XLSX file to start analysis")
  # Observe when the input file is selected:
  datasetInput <- reactive({
    req(input$input_xlsx) # Watch for a change to the input file picker
    req(input$skip) # Watch for a change to the skip value
    output$startup_message <- renderText("")
    withProgress(
      message = "loading data from XLSX", value = 0,
      {
        res <- read_growth_xlsx(input$input_xlsx$datapath, input$skip, input$show_full_sample_info)
      }
    )
    
    shinyjs::enable("generate_download")
    shinyjs::enable("output_prefix")
    return(res)
  })
  # Update the prefix when the filename changes:
  observeEvent(input$input_xlsx, {
    updateTextInput(session, "output_prefix", value = gsub(".xlsx", "", basename(input$input_xlsx$name), fixed=TRUE))
  })
  # Build the download from the loaded data:
  output$generate_download <- downloadHandler(
    filename = function(){
      sprintf("%s.zip", input$output_prefix)
    },
    content = function(file){
      shinyjs::disable("generate_download")
      d <- datasetInput()
      # Count the number of things to do:
      n_actions <- 7 + (3 * nlevels(d$metadata$strain)) + (nlevels(d$metadata$treatment) * (nlevels(d$metadata$strain) - 1))
      # Create a Progress object
      progress <- shiny::Progress$new()
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())
          progress$inc(1/n_actions, detail = "creating output directory structure")
          temp_directory <- file.path(tempdir(), as.integer(Sys.time()))
          dir.create(temp_directory)
          dir.create(file.path(temp_directory, "strains"))
          dir.create(file.path(temp_directory, "aggregated"))
          dir.create(file.path(temp_directory, "averaged"))
          dir.create(file.path(temp_directory, "parameters"))
          message(sprintf("temporary directory: %s", temp_directory))
          # Save the input data:
          progress$inc(1/n_actions, detail = "writing metadata")
          write.csv(d$metadata, file=file.path(temp_directory, sprintf("%s-metadata.csv", input$output_prefix)))
          progress$inc(1/n_actions, detail = "writing OD data")
          write.csv(d$OD, file=file.path(temp_directory, sprintf("%s-od_data.csv", input$output_prefix)), row.names=FALSE)
          # Save the fitted growth data:
          progress$inc(1/n_actions, detail = "writing growth data")
          write.csv(d$growth_data, file=file.path(temp_directory, sprintf("%s-growth_data.csv", input$output_prefix)), row.names=FALSE)
          # Save the growth_data plot:
          progress$inc(1/n_actions, detail = "writing well plot")
          ggsave(file.path(temp_directory, sprintf("%s-growth_wells.pdf", input$output_prefix)), plot=d$growth_data_plot, width=18, height=12)
          # Build the per-sample plots:
          for (sel_strain in levels(d$metadata$strain)) {
            progress$inc(1/n_actions, detail = sprintf("generating sample plot for %s", sel_strain))
            sel_label <- gsub("[/ ]", "_", sel_strain, perl = TRUE)
            plot_filename <- file.path(temp_directory, "strains", sprintf("%s-sample_%s.pdf", input$output_prefix, sel_label))
            # Select the per-sample data to plot:
            sel_metadata <- d$metadata[d$metadata$strain == sel_strain, , drop = FALSE]
            sel_plot_data <- d$plot_data[d$plot_data$sample %in% rownames(sel_metadata), ]
            sel_plot_data$sample <- factor(sel_plot_data$sample)
            sel_plot_data$treatment <- droplevels(sel_plot_data$treatment)
            # Plot the data:
            suppressWarnings({
              g <- ggplot(
                sel_plot_data,
                aes(x = time, y = OD, colour = treatment, group = sample)
              )
              g <- g + geom_point(
                size = 0.5
              )
              g <- g + geom_line()
              g <- g + xlab("Time (minutes)")
              g <- g + ylab(expression(OD[600]))
              g <- g + ggtitle(sel_strain)
              g <- g + theme(
                panel.grid = element_blank(),
                panel.background = element_rect(fill = "white", colour = "black")
              )
              ggsave(
                file = plot_filename,
                width = 12,
                height = 4
              )
            })
          }
          
          # Plot the aggregated sample plots:
          for(sel_strain in levels(d$metadata$strain)) {
            progress$inc(1/n_actions, detail = sprintf("generating aggregated plot for %s", sel_strain))
            sel_label <- gsub("[/ ]", "_", sel_strain, perl = TRUE)
            plot_filename <- file.path(temp_directory, "aggregated", sprintf("%s-aggregated_%s.pdf", input$output_prefix, sel_label))
            # Select the per-sample data to plot:
            sel_metadata <-d$metadata[d$metadata$strain == sel_strain, , drop = FALSE]
            sel_metadata$treatment <- droplevels(sel_metadata$treatment)
            # Calculate the mean of each group:
            sel_summary_data <- do.call(
              rbind.data.frame, 
              lapply(levels(sel_metadata$treatment), function(treatment){
                sel_samples <- rownames(sel_metadata[sel_metadata$treatment == treatment, ])
                sel_data <- d$norm_data[, sel_samples, drop = FALSE]
                sel_summary <- data.frame(
                  treatment = treatment,
                  time = d$norm_data$time,
                  n = ncol(sel_data),
                  mean = apply(sel_data, 1, mean),
                  sd = apply(sel_data, 1, sd)
                )
                sel_summary[is.na(sel_summary$sd), "sd"] <- 0
                return(sel_summary)
              })
            )
            sel_summary_data$treatment <- factor(sel_summary_data$treatment)
            sel_summary_data$ymin <- sel_summary_data$mean - sel_summary_data$sd
            sel_summary_data$ymax <- sel_summary_data$mean + sel_summary_data$sd
            # Plot the data:
            suppressWarnings({
              g <- ggplot(
                sel_summary_data,
                aes(x = time, y = mean, colour = treatment)
              )
              g <- g + geom_ribbon(
                aes(ymin = ymin, ymax = ymax, fill = treatment),
                alpha = 0.10,
                linewidth = 0.1
              )
              g <- g + geom_point(
                size = 0.5
              )
              g <- g + geom_line()
              g <- g + xlab("Time (minutes)")
              g <- g + ylab(expression(OD[600]))
              g <- g + ggtitle(sel_strain)
              g <- g + theme(
                panel.grid = element_blank(),
                panel.background = element_rect(fill = "white", colour = "black")
              )
              ggsave(
                file = plot_filename,
                width = 12,
                height = 4
              )
            })
          }
          # Plot the average plots against wildtype:
          # Loop over each treatment, plotting strain vs WT for each one:
          for(sel_treatment in levels(d$metadata$treatment)) {
            for(sel_strain in setdiff(levels(d$metadata$strain), input$reference_strain)) {
              progress$inc(1/n_actions, detail = sprintf("generating treatment plot for %s %s", sel_treatment, sel_strain))
              sel_strain_label <- gsub("[/ ]", "_", sel_strain, perl = TRUE)
              sel_treatment_label <- gsub("[/ ]", "_", sel_treatment, perl = TRUE)
              plot_filename <- file.path(temp_directory, "averaged", sprintf("%s-averaged_%s_%s.pdf", input$output_prefix, sel_strain_label, sel_treatment_label))
              # Pull out the data & metadata for this (strain, treatment) pair and (WT, treatment):
              sel_metadata <- d$metadata[d$metadata$strain %in% c(input$reference_strain, sel_strain) & d$metadata$treatment == sel_treatment,]
              if(nrow(sel_metadata) > 0) {
                sel_metadata$treatment <- droplevels(sel_metadata$treatment)
                sel_metadata$strain <- droplevels(sel_metadata$strain)
                sel_summary_data <- do.call(
                  rbind.data.frame,
                  lapply(levels(sel_metadata$strain), function(strain){
                    sel_samples <- rownames(sel_metadata[sel_metadata$strain == strain, ])
                    sel_data <- d$norm_data[, sel_samples, drop = FALSE]
                    sel_summary <- data.frame(
                      strain = strain,
                      time = d$norm_data$time,
                      n = ncol(sel_data),
                      mean = apply(sel_data, 1, mean),
                      sd = apply(sel_data, 1, sd)
                    )
                    sel_summary[is.na(sel_summary$sd), "sd"] <- 0
                    return(sel_summary)
                  })
                )
                sel_summary_data$strain <- factor(sel_summary_data$strain)
                sel_summary_data$ymin <- sel_summary_data$mean - sel_summary_data$sd
                sel_summary_data$ymax <- sel_summary_data$mean + sel_summary_data$sd
                # Plot the data:
                palette <- c("grey", "blue")
                names(palette) <- c("WT", sel_strain)
                suppressWarnings({
                  g <- ggplot(
                    sel_summary_data,
                    aes(x = time, y = mean, colour = strain)
                  )
                  g <- g + geom_ribbon(
                    aes(ymin = ymin, ymax = ymax, fill = strain),
                    alpha = 0.10,
                    linewidth = 0.1
                  )
                  g <- g + geom_point(
                    size = 0.5
                  )
                  g <- g + geom_line()
                  g <- g + scale_fill_manual(values=palette)
                  g <- g + scale_colour_manual(values=palette)
                  g <- g + xlab("Time (minutes)")
                  g <- g + ylab(expression(OD[600]))
                  g <- g + ggtitle(sprintf("strain %s ( treatment: %s)", sel_strain, sel_treatment))
                  g <- g + theme(
                    panel.grid = element_blank(),
                    panel.background = element_rect(fill = "white", colour = "black")
                  )
                  ggsave(
                    file = plot_filename,
                    width = 12,
                    height = 4
                  )
                })
              }
            }
          }
          
          # Test the treatment parameters:
          strain_tests <- do.call(rbind.data.frame, lapply(levels(d$metadata$strain), function(sel_strain){
            progress$inc(1/n_actions, detail = sprintf("generating treatment parameter plot for %s", sel_strain))
            
            sel_label <- gsub("[/ ]", "_", sel_strain, perl = TRUE)
            # Select the per-sample data:R
            sel_metadata <- d$metadata[d$metadata$strain == sel_strain, , drop = FALSE]
            # Pull out the growth rate parameters:
            sel_params <- data.frame(
              well_row = sel_metadata$row,
              well_col = sel_metadata$col,
              treatment = droplevels(sel_metadata$treatment),
              strain = sel_metadata$strain,
              k = d$growth_data[rownames(sel_metadata), "k"],
              n0 = d$growth_data[rownames(sel_metadata), "n0"],
              r = d$growth_data[rownames(sel_metadata), "r"],
              row.names = rownames(sel_metadata)
            )
            if(nlevels(sel_params$treatment) != 2) {
              return(data.frame(
                strain = sel_strain,
                n_none = length(which(sel_params$treatment == "none")),
                n_K28 = length(which(sel_params$treatment == "K28")),
                k_stat = NA,
                k_pvalue = NA,
                n0_stat = NA,
                n0_pvalue = NA,
                r_stat = NA,
                r_pvalue = NA
              ))
            }
            # Plot the parameter boxplots:
            for(parameter in c("r", "n0", "k")) {
              g <- ggplot(
                sel_params,
                aes(x = treatment, y = .data[[parameter]], fill = treatment)
              )
              g <- g + geom_boxplot()
              g <- g + ylab(parameter)
              g <- g + ggtitle(sprintf("%s - parameter %s", sel_strain, parameter))
              g <- g + theme(
                panel.grid = element_blank(),
                panel.background = element_rect(fill = "white", colour = "black")
              )
              ggsave(
                file = file.path(temp_directory, "parameters", sprintf("%s-parameter_%s_%s.pdf", input$output_prefix, parameter, sel_label)),
                width = 4,
                height = 6
              )
            }
            # Perform the tests:
            k_test <- t.test(k ~ treatment, data = sel_params)
            n0_test <- t.test(n0 ~ treatment, data = sel_params)
            r_test <- t.test(r ~ treatment, data = sel_params)
            sel_test_res <- data.frame(
              strain = sel_strain,
              n_none = length(which(sel_params$treatment == "none")),
              n_K28 = length(which(sel_params$treatment == "K28")),
              k_stat = k_test$statistic,
              k_pvalue = k_test$p.value,
              n0_stat = n0_test$statistic,
              n0_pvalue = n0_test$p.value,
              r_stat = r_test$statistic,
              r_pvalue = r_test$p.value
            )
            return(sel_test_res)
          }))
          rownames(strain_tests) <- NULL
          progress$inc(1/n_actions, detail = "writing parameter test data")
          write.csv(
            strain_tests,
            file = file.path(temp_directory, sprintf("%s-parameter_tests.csv", input$output_prefix)),
            row.names = FALSE
          )
          progress$inc(1/n_actions, detail = "building output ZIP")
          zip::zip(
            zipfile = file,
            files = dir(temp_directory),
            root = temp_directory
          )
      shinyjs::enable("generate_download")
      shinyjs::enable("output_prefix")
      
    },
    contentType = "application/zip"
  )
  
  output$well_plot <- renderPlot({
    datasetInput()$growth_data_plot
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
