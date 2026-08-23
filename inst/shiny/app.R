# ============================================================================
# EpiTracer web app -- upload copy-number + structural-variant calls, plot them
# with plot_sv_linear(), and (optionally) call amplicon-formation mechanisms.
# Input reading/reformatting is done by the EpiTracer package (read_cnv/read_sv/
# read_purple_sv_vcf/gene_locus/prepare_amplicon_inputs); this file is UI + glue.
#
# Run locally:  shiny::runApp("inst/shiny")
# Deploy:       source("inst/shiny/deploy.R")   # -> shinyapps.io
# ============================================================================

suppressPackageStartupMessages({ library(shiny); library(bslib); library(DT); library(data.table) })
source("helpers.R", local = TRUE)   # onc_gr, seeds_for_sample, call_mechanisms
options(shiny.maxRequestSize = 200 * 1024^2)
nz <- function(x) if (is.null(x) || !nzchar(trimws(x))) NULL else trimws(x)

ui <- page_sidebar(
  title = "EpiTracer — plot & call amplicon mechanisms",
  theme = bs_theme(version = 5, primary = "#1E7A6F",
                   base_font = font_google("IBM Plex Sans"),
                   heading_font = font_google("IBM Plex Serif"),
                   code_font = font_google("IBM Plex Mono")),
  sidebar = sidebar(
    width = 340,
    fileInput("cn", "Copy-number segments", accept = c(".tsv",".csv",".txt",".rds",".gz")),
    fileInput("sv", "Structural variants (BEDPE or VCF)", accept = c(".tsv",".csv",".txt",".bedpe",".vcf",".rds",".gz")),
    helpText("CN: generic TSV/CSV, PURPLE-style .rds, or a raw PURPLE",
             tags$code(".purple.cnv.somatic.tsv"), "(sample & ploidy inferred).",
             "SV: a bedpe with VF/JCN, or a PURPLE/GRIDSS", tags$code(".vcf(.gz)"),
             "(VF & PURPLE_JCN preserved)."),
    div(class = "small text-success", textOutput("cninfo")),
    textInput("samplename", "Sample name (optional)", placeholder = "one name for both files"),
    numericInput("ploidy", "Ploidy (0 = auto-estimate)", value = 0, min = 0, step = 0.1),
    selectInput("genome", "Genome", c("hg38","hg19","mm10")),
    selectInput("sample", "Sample", choices = NULL),
    radioButtons("target", "Target", c("None (untargeted)" = "none", "Amplicons" = "amp", "Homozygous deletions" = "homdel"), selected = "amp"),
    radioButtons("region", "Region",
      c("Auto: all target loci" = "auto", "By gene(s)" = "gene",
        "Pick a detected amplicon" = "amp", "Whole chromosome" = "chr", "Custom locus" = "custom"),
      selected = "auto"),
    conditionalPanel("input.region == 'gene'", textInput("genes", "Gene symbol(s)", value = "EGFR", placeholder = "EGFR, CDKN2A")),
    conditionalPanel("input.region == 'amp'", selectInput("ampchoice", "Amplicon", choices = NULL)),
    conditionalPanel("input.region == 'chr'", selectInput("chrchoice", "Chromosome", choices = NULL)),
    conditionalPanel("input.region == 'custom'", textInput("locus", "chr:start-end", placeholder = "chr7:54000000-56000000")),
    checkboxInput("showonc", "Show oncogene labels", FALSE),
    checkboxInput("mech", "Call amplicon-formation mechanisms", TRUE),
    actionButton("go", "Plot", class = "btn-primary w-100"),
    downloadButton("dl", "Download PDF", class = "btn-outline-secondary w-100 mt-2")
  ),
  navset_card_tab(
    nav_panel("Plot",
      conditionalPanel("input.go == 0",
        div(class = "text-muted", style = "padding:2rem;",
            "Upload copy-number + SV files, name the sample if needed, choose target and region, then press Plot.")),
      conditionalPanel("input.go > 0",
        div(style = "padding-top:3rem;",
          div(class = "d-flex align-items-center gap-2 mb-2 flex-wrap",
            actionButton("panL", "<< Pan left",  class = "btn-sm btn-outline-secondary"),
            actionButton("zoomOut", "- Zoom out", class = "btn-sm btn-outline-secondary"),
            actionButton("zoomIn",  "+ Zoom in",  class = "btn-sm btn-outline-secondary"),
            actionButton("panR", "Pan right >>", class = "btn-sm btn-outline-secondary"),
            actionButton("resetview", "Reset view", class = "btn-sm btn-outline-secondary"),
            span(class = "small text-muted ms-2 font-monospace", textOutput("winlab", inline = TRUE)),
            div(class = "ms-auto d-flex align-items-center gap-2",
              tags$span(class = "small text-muted", "width"),
              sliderInput("pw", NULL, min = 900, max = 4500, value = 1800, step = 100, ticks = FALSE, width = "180px"))),
          div(id = "dragbox",
              style = "overflow-x:auto; cursor:grab; user-select:none; border:1px solid var(--bs-border-color); border-radius:8px;",
              uiOutput("plotbox")),
          div(class = "small text-muted mt-1", "Tip: click and drag the plot to pan, double-click to zoom in, and use Reset view to return to the starting window."),
          tags$script(HTML(
"(function(){function init(){var box=document.getElementById('dragbox');if(!box){setTimeout(init,300);return;}",
"var down=false,x0=0,dx=0,img=null;",
"box.addEventListener('mousedown',function(e){img=box.querySelector('img');if(!img)return;img.draggable=false;down=true;x0=e.clientX;dx=0;box.style.cursor='grabbing';e.preventDefault();});",
"document.addEventListener('mousemove',function(e){if(!down||!img)return;dx=e.clientX-x0;img.style.transform='translateX('+dx+'px)';});",
"document.addEventListener('mouseup',function(){if(!down)return;down=false;box.style.cursor='grab';if(img)img.style.transform='';",
"if(img&&Math.abs(dx)>3){Shiny.setInputValue('drag_dx',{dx:dx,w:img.clientWidth,n:Date.now()},{priority:'event'});}});}",
"if(document.readyState!=='loading')init();else document.addEventListener('DOMContentLoaded',init);})();"))))),
    nav_panel("Mechanisms", DTOutput("mtab")),
    nav_panel("Detected amplicons", DTOutput("stab")),
    nav_panel("Data preview", h6("Copy number"), DTOutput("cnprev"), h6("Structural variants"), DTOutput("svprev"))
  )
)

server <- function(input, output, session) {

  cn <- reactive({ req(input$cn)
    read_cnv(input$cn$datapath, sample = nz(input$samplename),
             ploidy = if (isTRUE(input$ploidy > 0)) input$ploidy else NULL, name = input$cn$name) })
  sv <- reactive({ req(input$sv)
    read_sv(input$sv$datapath, sample = nz(input$samplename), name = input$sv$name) })

  output$cninfo <- renderText({ req(input$cn); c <- cn()
    sprintf("CN: sample '%s' · ploidy %s%s · %d segments", unique(c$sample)[1], unique(c$ploidy)[1],
            if (isTRUE(input$ploidy > 0)) " (set)" else " (auto)", nrow(c)) })

  observeEvent(list(cn(), sv()), {
    updateSelectInput(session, "sample", choices = unique(cn()$sample), selected = unique(cn()$sample)[1])
    updateSelectInput(session, "chrchoice", choices = sort(unique(cn()$seqnames)))
  })

  ## SV rows for the selected sample, relabelled to it when the file is single-sample
  sv_for <- reactive({ req(input$sample); x <- copy(sv())
    if (uniqueN(x$sample) == 1) x[, sample := input$sample]
    y <- x[sample == input$sample]; if (nrow(y)) y else x })

  inp <- reactive({ req(input$sample)
    prepare_amplicon_inputs(cn(), sv_for(), sample = input$sample) })

  seeds <- reactive({ req(input$sample); seeds_for_sample(inp()) })

  observeEvent(seeds(), {
    sd <- seeds()
    if (length(sd)) updateSelectInput(session, "ampchoice",
      choices = stats::setNames(seq_along(sd), sprintf("%s  %s:%s-%s", sd$ID, as.character(seqnames(sd)),
        format(start(sd), big.mark=",", trim=TRUE), format(end(sd), big.mark=",", trim=TRUE))))
    else updateSelectInput(session, "ampchoice", choices = c("(none detected)" = ""))
  })

  plot_region <- reactive({
    switch(input$region,
      auto = {
        if (input$target == "none") {
          # Untargeted auto view: show every chromosome that carries or is linked by an SV.
          svd <- sv_for(); chrs <- unique(c(.pfx(as.character(svd$chrom1)), .pfx(as.character(svd$chrom2))))
          chrs <- chrs[!is.na(chrs) & nzchar(chrs)]
          cl <- tryCatch(names(load_chrom_lengths(input$genome)), error = function(e) NULL)
          chrs <- if (!is.null(cl)) chrs[order(match(chrs, cl))] else sort(chrs)
          validate(need(length(chrs) > 0, "No SVs found to display. Pick a target or region."))
          list(loci = NULL, chromosome = chrs, genes = NULL)
        } else list(loci = NULL, chromosome = NULL, genes = NULL)
      },
      gene = { gs <- strsplit(input$genes %||% "", "[,;\\s]+")[[1]]; gs <- gs[nzchar(gs)]
        validate(need(length(gs) > 0, "Enter one or more gene symbols, e.g. EGFR, CDKN2A"))
        loc <- tryCatch(gene_locus(gs, genome = input$genome, cnv = cn(), sample = input$sample,
                                   target = if (input$target == "none") "any" else input$target),
                        error = function(e) validate(need(FALSE, conditionMessage(e))))
        list(loci = loc[, c("chr","start","end")], chromosome = NULL, genes = loc$gene) },
      amp = { sd <- seeds(); i <- as.integer(input$ampchoice); req(length(sd), !is.na(i)); s <- sd[i]
        list(loci = data.frame(chr = as.character(seqnames(s)), start = pmax(1, start(s) - 5e5),
             end = end(s) + 5e5), chromosome = NULL, genes = NULL) },
      chr = list(loci = NULL, chromosome = input$chrchoice, genes = NULL),
      custom = { m <- regmatches(input$locus, regexec("^\\s*(chr[^:]+|[^:]+):([0-9,]+)-([0-9,]+)", input$locus))[[1]]
        validate(need(length(m) == 4, "Enter a locus like chr7:54000000-56000000"))
        list(loci = data.frame(chr = .pfx(m[2]), start = as.numeric(gsub(",","",m[3])),
             end = as.numeric(gsub(",","",m[4]))), chromosome = NULL, genes = NULL) })
  })

  ## The plotted view is a single genomic window we can pan/zoom, captured on Plot.
  base_args <- reactiveVal(NULL); view <- reactiveVal(NULL); home_view <- reactiveVal(NULL)
  chrom_len <- function(chr) {
    cl <- tryCatch(load_chrom_lengths(input$genome), error = function(e) NULL)
    if (!is.null(cl) && chr %in% names(cl)) cl[[chr]] else suppressWarnings(max(cn()[seqnames == chr]$end, na.rm = TRUE)) }

  observeEvent(input$go, {
    req(input$sample); reg <- plot_region()
    svd <- copy(sv_for()); svd[, sample := input$sample]
    ev <- if (input$target == "none") c("amp","gain","loh","homdel") else input$target
    base_args(list(sample = input$sample, cnv_data = as.data.frame(cn()[sample == input$sample]),
                   sv_data = as.data.frame(svd), genome = input$genome, events = ev,
                   genes = reg$genes, chromosome = reg$chromosome))
    v0 <- if (!is.null(reg$loci) && nrow(reg$loci) >= 1)
      list(chr = reg$loci$chr[1], start = min(reg$loci$start), end = max(reg$loci$end))
    else if (!is.null(reg$chromosome) && length(reg$chromosome) == 1)
      list(chr = reg$chromosome, start = 1, end = chrom_len(reg$chromosome))
    else NULL     # multi-chromosome / whole-genome: nothing single to pan
    view(v0); home_view(v0)
  })

  shift_win <- function(frac) { v <- view(); if (is.null(v)) return()
    sp <- v$end - v$start; view(list(chr = v$chr, start = max(1, v$start + frac*sp), end = max(2, v$end + frac*sp))) }
  zoom_win  <- function(k)    { v <- view(); if (is.null(v)) return()
    c0 <- (v$start + v$end)/2; sp <- (v$end - v$start) * k
    view(list(chr = v$chr, start = max(1, c0 - sp/2), end = c0 + sp/2)) }
  # Zoom keeping a chosen genomic point centred (used by the double-click handler).
  zoom_at <- function(cx, k) { v <- view(); if (is.null(v)) return()
    cx <- max(v$start, min(v$end, cx)); sp <- (v$end - v$start) * k
    view(list(chr = v$chr, start = max(1, cx - sp/2), end = cx + sp/2)) }
  observeEvent(input$panL,    shift_win(-0.1)); observeEvent(input$panR, shift_win(0.1))
  observeEvent(input$zoomIn,  zoom_win(1/1.7)); observeEvent(input$zoomOut, zoom_win(1.7))
  observeEvent(input$resetview, view(home_view()))
  # Double-click to zoom in, centred on the clicked position. The plot x-axis is a
  # concatenated coord where, for a single-window view, gx = pos - view_start.
  observeEvent(input$plot_dblclick, { v <- view(); dc <- input$plot_dblclick
    if (!is.null(v) && !is.null(dc$x)) zoom_at(v$start + dc$x, 1/1.7) })
  # Click-and-drag panning: dx = pixels dragged, w = rendered image width in px.
  # Dragging content right (dx > 0) pulls earlier coordinates into view (shift left).
  observeEvent(input$drag_dx, {
    d <- input$drag_dx; w <- as.numeric(d$w); dx <- as.numeric(d$dx)
    if (!is.null(view()) && is.finite(w) && w > 0) shift_win(-(dx / w))
  }, ignoreInit = TRUE)

  output$winlab <- renderText({ v <- view()
    if (is.null(v)) "whole-genome view (pan/zoom off)"
    else sprintf("%s:%s-%s", v$chr, format(round(v$start), big.mark=","), format(round(v$end), big.mark=",")) })

  current_plot <- reactive({ ba <- base_args(); req(ba); v <- view()
    args <- ba[c("sample","cnv_data","sv_data","genome","events")]
    if (!is.null(v)) args$loci <- data.frame(chr = v$chr, start = round(max(1, v$start)), end = round(v$end))
    else if (!is.null(ba$chromosome)) args$chromosome <- ba$chromosome
    if (!is.null(ba$genes)) args$genes_to_highlight <- ba$genes
    # Oncogene panel off by default: feed an empty gene table unless the user ticks
    # "Show oncogene labels" (gene-region mode keeps its own highlighted genes).
    if (is.null(ba$genes) && !isTRUE(input$showonc))
      args$gene_coord <- data.frame(chr = character(0), start = numeric(0),
                                    end = numeric(0), strand = character(0), gene = character(0))
    tryCatch(do.call(plot_sv_linear, args),
             error = function(e) structure(list(), class = "epitracer_error", msg = conditionMessage(e))) })

  output$plotbox <- renderUI({ req(input$go > 0)
    plotOutput("plot", width = paste0(input$pw, "px"), height = "440px", dblclick = "plot_dblclick") })
  output$plot <- renderPlot({ p <- current_plot()
    validate(need(!inherits(p, "epitracer_error"), paste("Could not plot:", attr(p, "msg"))))
    print(p) }, res = 96)

  mech_tab <- eventReactive(input$go, { req(input$mech)
    withProgress(message = "Calling mechanisms…", value = 0.5, { call_mechanisms(inp(), input$genome) }) })
  output$mtab <- renderDT({
    if (!isTRUE(input$mech)) return(datatable(data.frame(Note = "Mechanism calling is switched off.")))
    m <- mech_tab(); validate(need(nrow(m) > 0, "No focal amplicons detected in this sample."))
    datatable(m, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE)) })

  output$stab <- renderDT({ sd <- seeds()
    validate(need(length(sd) > 0, "No focal amplicon seeds detected (need copyNumber > 3 x ploidy)."))
    datatable(data.frame(ID = sd$ID, chr = as.character(seqnames(sd)), start = start(sd), end = end(sd),
              width_Mb = round(width(sd)/1e6, 2)), rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE)) })

  output$cnprev <- renderDT(datatable(head(cn(), 200), rownames = FALSE, options = list(scrollX = TRUE)))
  output$svprev <- renderDT(datatable(head(sv(), 200), rownames = FALSE, options = list(scrollX = TRUE)))

  output$dl <- downloadHandler(
    filename = function() sprintf("epitracer_%s.pdf", input$sample %||% "plot"),
    content = function(file) { p <- current_plot(); validate(need(!inherits(p, "epitracer_error"), "Nothing to download."))
      ggplot2::ggsave(file, p, width = 16, height = 6, device = grDevices::cairo_pdf, bg = "white") })
}

shinyApp(ui, server)
