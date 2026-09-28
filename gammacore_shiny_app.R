## =============================================================================
## gammaCore (nVNS) for chronic cluster headache -- Shiny app
## Interactive version of the base-R cost-effectiveness (cost-minimisation)
## model replicating the electroCore sponsor submission to NICE, MT323.
##
## Run with:  shiny::runApp("app.R")
##
## Structure deliberately mirrors the companion Excel workbook:
##   - "Inputs" tab  <-> Excel "Inputs" sheet: every hardcoded parameter is
##      editable here (unit costs, probabilities, the dose table for all 5
##      responder definitions, and the uncertainty ranges used by OWSA).
##   - Sidebar scenario selectors <-> Excel's yellow "scenario selector" cells.
##   - "Base case" / "Scenarios" / "OWSA" tabs <-> the equivalent Excel sheets,
##      all driven live off whatever is currently set on the "Inputs" tab.
## =============================================================================

library(shiny)
library(shinythemes)
library(ggplot2)
library(DT)
library(scales)

## -----------------------------------------------------------------------------
## STATIC STRUCTURE (not user-editable): which cell each parameter maps to.
## Point values and uncertainty ranges themselves live in reactive tables
## below and ARE user-editable, exactly like the blue cells in the Excel
## Inputs sheet.
## -----------------------------------------------------------------------------

days_per_month <- 365.25 / 12

default_unit_costs <- data.frame(
  Parameter = c("Zolmitriptan nasal - cost per unit (\u00a3)",
                "Sumatriptan s.c. - cost per unit (\u00a3)",
                "Sumatriptan nasal - cost per unit (\u00a3)",
                "Oxygen static - cost per treatment (\u00a3)",
                "Oxygen portable - cost per treatment (\u00a3)",
                "% of sumatriptan treatments that are s.c.",
                "% of oxygen treatments that are portable"),
  key = c("cost_zolm", "cost_suma_sc", "cost_suma_nasal", "cost_o2_static", "cost_o2_portable", "frac_sc", "frac_o2p"),
  Value = c(6.08, 19.75, 7.08, 0.56, 0.79, 0.867, 0.50),
  stringsAsFactors = FALSE
)

default_gcore <- data.frame(
  Parameter = c("gammaCore cost per 3-month refill (\u00a3)", "Number of free initial months"),
  key = c("gc_price", "gc_free"),
  Value = c(625, 3),
  stringsAsFactors = FALSE
)

default_probs <- data.frame(
  Parameter = c("Probability of response (\u226550%) - SoC, month 1", "Probability of discontinued response per month - initial gCore responders"),
  key = c("p_resp_soc", "p_loss"),
  Value = c(0.083, 0.31),
  stringsAsFactors = FALSE
)

## Dose table: point estimates per responder definition (Excel Inputs Sec. 5)
default_dose_grid <- data.frame(
  Definition = c("25", "40", "50", "65", "50m"),
  `P(response) gCore` = c(0.60, 0.47, 0.40, 0.24, 0.40),
  `Responder zolm` = c(0.80, 1.00, 0.60, 0.00, 1.60),
  `Responder suma` = c(2.50, 3.00, 2.50, 2.00, 2.80),
  `Responder oxy` = c(3.50, 2.80, 2.20, 0.70, 6.50),
  `NonRespOnTx zolm` = c(3.80, 2.50, 2.50, 2.20, 1.30),
  `NonRespOnTx suma` = c(5.80, 3.70, 4.10, 3.80, 7.50),
  `NonRespOnTx oxy` = c(16.20, 12.20, 11.20, 9.20, 10.80),
  `NonRespBaseline zolm` = c(4.80, 3.30, 3.80, 3.80, 1.30),
  `NonRespBaseline suma` = c(3.80, 2.90, 4.50, 5.00, 7.50),
  `NonRespBaseline oxy` = c(16.30, 19.30, 18.60, 31.00, 10.80),
  check.names = FALSE, stringsAsFactors = FALSE
)

default_socnr <- data.frame(zolm = 1.30, suma = 7.50, oxy = 10.80)

## Uncertainty ranges for OWSA (Excel OWSA sheet / Table C13). These are the
## 10 parameters (all keyed to the 50% responder definition) with a documented
## one-way sensitivity range in the sponsor submission.
default_ranges <- data.frame(
  key    = c("p_resp_gc", "p_loss",
             "resp_zolm", "resp_suma",
             "onTx_zolm", "onTx_suma",
             "socnr_zolm", "socnr_suma", "socnr_oxy",
             "frac_sc"),
  Parameter = c("P(response) - gammaCore", "P(discontinued response)/month",
                "Responder dose/14d - zolmitriptan", "Responder dose/14d - sumatriptan",
                "Non-resp (on Tx) dose/14d - zolmitriptan", "Non-resp (on Tx) dose/14d - sumatriptan",
                "SoC non-responder dose/14d - zolmitriptan", "SoC non-responder dose/14d - sumatriptan", "SoC non-responder dose/14d - oxygen",
                "% sumatriptan treatments that are s.c."),
  Mean  = c(0.40, 0.31, 0.60, 2.50, 2.50, 4.10, 1.30, 7.50, 10.80, 0.867),
  Lower = c(0.263, 0.162, 0.10, 1.04, 0.32, 1.00, 0.45, 4.88, 6.68, 0.661),
  Upper = c(0.545, 0.541, 1.52, 4.59, 6.89, 9.33, 2.59, 10.67, 15.90, 0.982),
  stringsAsFactors = FALSE
)

def_labels <- c("25" = "\u226525% responders", "40" = "\u226540% responders",
                "50" = "\u226550% responders (base case)", "65" = "\u226565% responders",
                "50m" = "\u226550% (whole-arm means)")
loss_labels <- c("none" = "Single loss after month 1 (base case)",
                 "constant" = "Constant monthly loss",
                 "diminishing" = "Diminishing monthly loss")
nonresp_labels <- c("soc" = "Revert to SoC use (base case)", "baseline" = "Revert to PREVA baseline use")

## selectInput() shows the NAMES of `choices` as labels and returns the VALUES.
## We want the long descriptive text shown, and the short code returned -- so
## build the choices vector the opposite way round from the lookup vectors above.
def_choices      <- setNames(names(def_labels), def_labels)
loss_choices     <- setNames(names(loss_labels), loss_labels)
nonresp_choices  <- setNames(names(nonresp_labels), nonresp_labels)

## -----------------------------------------------------------------------------
## MODEL ENGINE -- now takes the full dose table / SoC dose / unit costs as
## explicit arguments (defaulting to the published Table C10 values), so the
## Shiny server can pass through whatever the user has edited on the Inputs tab.
## -----------------------------------------------------------------------------

med_cost <- function(dose, frac_sc, frac_o2_portable, cost_zolm, cost_suma_sc, cost_suma_nasal, cost_o2_static, cost_o2_portable) {
  f <- days_per_month / 14
  suma_unit <- frac_sc * cost_suma_sc + (1 - frac_sc) * cost_suma_nasal
  o2_unit   <- frac_o2_portable * cost_o2_portable + (1 - frac_o2_portable) * cost_o2_static
  unname(dose["zolm"] * f * cost_zolm + dose["suma"] * f * suma_unit + dose["oxy"] * f * o2_unit)
}

#' @param dose_table_arg named list (by definition code) of list(resp=,onTx=,base=,p=) -- from the Inputs tab
#' @param soc_nr_arg named numeric vector c(zolm=,suma=,oxy=) -- from the Inputs tab
run_model <- function(responder_def = "50",
                      loss_type = c("none", "constant", "diminishing"),
                      nonresp_use = c("soc", "baseline"),
                      soc_no_response = FALSE,
                      gcore_first_n_free, gcore_cost_3m,
                      p_resp_soc_use, p_loss_use,
                      dose_table_arg, soc_nr_arg,
                      frac_sc, frac_o2_portable,
                      cost_zolm, cost_suma_sc, cost_suma_nasal, cost_o2_static, cost_o2_portable,
                      p_resp_gc_override = NULL, dose_override = NULL, ...) {
  
  loss_type   <- match.arg(loss_type)
  nonresp_use <- match.arg(nonresp_use)
  n_cycles <- 12
  
  dd <- dose_table_arg[[responder_def]]
  d_resp <- dd$resp; d_onTx <- dd$onTx; d_base <- dd$base
  p_resp_gc <- if (!is.null(p_resp_gc_override)) p_resp_gc_override else dd$p
  d_soc_nr <- soc_nr_arg
  
  if (!is.null(dose_override)) {
    if (!is.null(dose_override$resp))   d_resp   <- dose_override$resp
    if (!is.null(dose_override$onTx))   d_onTx   <- dose_override$onTx
    if (!is.null(dose_override$base))   d_base   <- dose_override$base
    if (!is.null(dose_override$soc_nr)) d_soc_nr <- dose_override$soc_nr
  }
  
  d_nonresp_post <- if (nonresp_use == "soc") d_soc_nr else d_base
  gcore_month_cost <- gcore_cost_3m / 3
  mc <- function(d) med_cost(d, frac_sc, frac_o2_portable, cost_zolm, cost_suma_sc, cost_suma_nasal, cost_o2_static, cost_o2_portable)
  
  cyc <- data.frame(cycle = 1:n_cycles, r_gc = NA_real_, cost_gc = NA_real_, r_soc = NA_real_, cost_soc = NA_real_)
  r_prev <- p_resp_gc
  for (t in 1:n_cycles) {
    if (t == 1) {
      r_gc <- p_resp_gc
    } else if (loss_type == "none") {
      r_gc <- if (t == 2) p_resp_gc * (1 - p_loss_use) else r_prev
    } else if (loss_type == "constant") {
      r_gc <- r_prev * (1 - p_loss_use)
    } else {
      p_t <- p_loss_use * (0.9)^(t - 2)
      r_gc <- r_prev * (1 - p_t)
    }
    r_prev <- r_gc
    in_free_trial <- t <= gcore_first_n_free
    if (in_free_trial) {
      cost_resp <- mc(d_resp); cost_nonresp <- mc(d_onTx)
      cost_gc_month <- r_gc * cost_resp + (1 - r_gc) * cost_nonresp
    } else {
      cost_resp <- mc(d_resp) + gcore_month_cost; cost_nonresp <- mc(d_nonresp_post)
      cost_gc_month <- r_gc * cost_resp + (1 - r_gc) * cost_nonresp
    }
    if (t == 1 && !soc_no_response) {
      r_soc <- p_resp_soc_use
      cost_soc_month <- r_soc * mc(d_resp) + (1 - r_soc) * mc(d_soc_nr)
    } else {
      r_soc <- 0
      cost_soc_month <- mc(d_soc_nr)
    }
    cyc$r_gc[t] <- r_gc; cyc$cost_gc[t] <- cost_gc_month
    cyc$r_soc[t] <- r_soc; cyc$cost_soc[t] <- cost_soc_month
  }
  list(cycles = cyc, total_gc = sum(cyc$cost_gc), total_soc = sum(cyc$cost_soc),
       diff = sum(cyc$cost_gc) - sum(cyc$cost_soc))
}

disaggregate <- function(model_result, responder_def, nonresp_use, dose_table_arg, soc_nr_arg,
                         frac_sc, frac_o2_portable, cost_zolm, cost_suma_sc, cost_suma_nasal,
                         cost_o2_static, cost_o2_portable, gcore_cost_3m, gcore_first_n_free, ...) {
  dd <- dose_table_arg[[responder_def]]
  d_resp <- dd$resp; d_onTx <- dd$onTx; d_base <- dd$base; d_soc_nr <- soc_nr_arg
  d_nonresp_post <- if (nonresp_use == "soc") d_soc_nr else d_base
  f <- days_per_month / 14
  suma_unit <- frac_sc * cost_suma_sc + (1 - frac_sc) * cost_suma_nasal
  o2_unit   <- frac_o2_portable * cost_o2_portable + (1 - frac_o2_portable) * cost_o2_static
  gcore_month_cost <- gcore_cost_3m / 3
  
  gc_zolm <- 0; gc_suma <- 0; gc_oxy <- 0; gc_device <- 0
  soc_zolm <- 0; soc_suma <- 0; soc_oxy <- 0
  for (t in 1:12) {
    r_gc <- model_result$cycles$r_gc[t]; r_soc <- model_result$cycles$r_soc[t]
    in_free <- t <= gcore_first_n_free
    d_nr_gc <- if (in_free) d_onTx else d_nonresp_post
    gc_zolm <- gc_zolm + r_gc * d_resp["zolm"] * f * cost_zolm + (1 - r_gc) * d_nr_gc["zolm"] * f * cost_zolm
    gc_suma <- gc_suma + r_gc * d_resp["suma"] * f * suma_unit + (1 - r_gc) * d_nr_gc["suma"] * f * suma_unit
    gc_oxy  <- gc_oxy  + r_gc * d_resp["oxy"]  * f * o2_unit   + (1 - r_gc) * d_nr_gc["oxy"]  * f * o2_unit
    if (!in_free) gc_device <- gc_device + r_gc * gcore_month_cost
    if (t == 1) {
      soc_zolm <- soc_zolm + r_soc * d_resp["zolm"] * f * cost_zolm + (1 - r_soc) * d_soc_nr["zolm"] * f * cost_zolm
      soc_suma <- soc_suma + r_soc * d_resp["suma"] * f * suma_unit + (1 - r_soc) * d_soc_nr["suma"] * f * suma_unit
      soc_oxy  <- soc_oxy  + r_soc * d_resp["oxy"]  * f * o2_unit   + (1 - r_soc) * d_soc_nr["oxy"]  * f * o2_unit
    } else {
      soc_zolm <- soc_zolm + d_soc_nr["zolm"] * f * cost_zolm
      soc_suma <- soc_suma + d_soc_nr["suma"] * f * suma_unit
      soc_oxy  <- soc_oxy  + d_soc_nr["oxy"]  * f * o2_unit
    }
  }
  data.frame(Item = c("gammaCore device", "Sumatriptan", "Zolmitriptan", "Oxygen", "Total"),
             gammaCore_SoC = round(c(gc_device, gc_suma, gc_zolm, gc_oxy, gc_device + gc_suma + gc_zolm + gc_oxy), 2),
             SoC_alone = round(c(0, soc_suma, soc_zolm, soc_oxy, soc_suma + soc_zolm + soc_oxy), 2))
}

dose_override_helper <- function(dose_table_arg, soc_nr_arg, which_dose, drug, value) {
  base_doses <- list(resp = dose_table_arg[["50"]]$resp, onTx = dose_table_arg[["50"]]$onTx,
                     base = dose_table_arg[["50"]]$base, soc_nr = soc_nr_arg)
  d <- base_doses[[which_dose]]; d[drug] <- value
  ov <- list(); ov[[which_dose]] <- d; ov
}

## Maps each of the 16 uncertainty-range keys to which dose-table cell (or
## scalar model argument) it perturbs. Static structure, not user-editable.
param_map <- list(
  p_resp_gc  = list(type = "scalar", target = "p_resp_gc"),
  p_resp_soc = list(type = "scalar", target = "p_resp_soc"),
  p_loss     = list(type = "scalar", target = "p_loss"),
  resp_zolm  = list(type = "dose", which = "resp", drug = "zolm"),
  resp_suma  = list(type = "dose", which = "resp", drug = "suma"),
  resp_oxy   = list(type = "dose", which = "resp", drug = "oxy"),
  onTx_zolm  = list(type = "dose", which = "onTx", drug = "zolm"),
  onTx_suma  = list(type = "dose", which = "onTx", drug = "suma"),
  onTx_oxy   = list(type = "dose", which = "onTx", drug = "oxy"),
  socnr_zolm = list(type = "dose", which = "soc_nr", drug = "zolm"),
  socnr_suma = list(type = "dose", which = "soc_nr", drug = "suma"),
  socnr_oxy  = list(type = "dose", which = "soc_nr", drug = "oxy"),
  frac_sc    = list(type = "scalar", target = "frac_sc"),
  frac_o2p   = list(type = "scalar", target = "frac_o2p"),
  o2_static  = list(type = "scalar", target = "o2_static"),
  o2_portable= list(type = "scalar", target = "o2_portable")
)

## -----------------------------------------------------------------------------
## UI
## -----------------------------------------------------------------------------

ui <- fluidPage(
  theme = shinytheme("flatly"),
  tags$head(tags$style(HTML("
    .well { background-color: #f4f6f8; }
    .value-box { background: #1F4E78; color: white; border-radius: 6px; padding: 14px; text-align: center; margin-bottom: 10px;}
    .value-box .big { font-size: 26px; font-weight: 700; }
    .value-box .lbl { font-size: 12px; opacity: 0.85; }
    .save-box { background: #1F9E89; }
    .caveat { background-color: #FFF2CC; border-left: 4px solid #E9A23B; padding: 10px 14px; margin-bottom: 14px; font-size: 13px; }
    .section-hdr { background-color: #1F4E78; color: white; padding: 6px 10px; font-weight: 600; margin-top: 14px; margin-bottom: 8px; border-radius: 4px; }
    .well .help-block { font-size: 13px; color: #6c757d; margin-top: -6px; margin-bottom: 12px; line-height: 1.35; }
    /* Make DT's inline cell editor clearly visible: explicit dark text on a white
       background with a visible border, regardless of the surrounding row's
       banding color or the Bootstrap theme's default (unstyled) input look. */
    table.dataTable td input, table.dataTable td textarea {
      color: #1a1a1a !important;
      background-color: #ffffff !important;
      border: 2px solid #1F4E78 !important;
      border-radius: 3px !important;
      font-size: 14px !important;
      padding: 4px 6px !important;
      box-shadow: 0 0 5px rgba(31, 78, 120, 0.5) !important;
    }
  "))),
  titlePanel("gammaCore for Chronic Cluster Headache \u2014 Interactive Cost Model (MT323)"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Scenario"),
      selectInput("responder_def", "Responder definition", choices = def_choices, selected = "50"),
      helpText("Threshold used to classify a patient as a treatment responder (e.g. 50% = \u226550% reduction in weekly attack frequency from baseline). \u201850% (whole-arm means)\u2019 uses the same threshold but a different dose-averaging method. OWSA always use the 50% definition, regardless of this setting."),
      selectInput("loss_type", "Response-loss assumption", choices = loss_choices, selected = "none"),
      helpText("How the gammaCore responder fraction changes over the 12 months: a single one-off drop after month 1 (base case), a constant loss applied every month, or a loss rate that shrinks by 10% each month."),
      selectInput("nonresp_use", "Non-responder use after discontinuation", choices = nonresp_choices, selected = "soc"),
      helpText("From month 4 onward, what medication-use level non-responders revert to after stopping gammaCore: the standard-of-care-arm average (base case), or their own pre-trial (PREVA baseline) level."),
      checkboxInput("soc_no_response", "SoC arm: no initial response (scenario)", value = FALSE),
      helpText("When checked, no patients in the SoC-alone arm are assumed to respond in month 1 (overrides the 8.3% base-case SoC response rate). Explores how sensitive the result is to that assumption."),
      hr(),
      helpText("All other parameters \u2014 unit costs, probabilities, dose tables, and uncertainty ranges \u2014 ",
               "are editable on the ", strong("Inputs"), " tab, exactly like the blue cells in the companion Excel workbook."),
      actionButton("reset_all", "Reset all inputs to published values", icon = icon("rotate-left"), class = "btn-warning btn-sm")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        
        tabPanel("Inputs",
                 br(),
                 
                 div(class = "section-hdr", "1. Unit costs (NHS Drug Tariff, March 2019 / East of England PAC, 2017)"),
                 DTOutput("unitCostTable"),
                 
                 div(class = "section-hdr", "2. gammaCore device pricing"),
                 DTOutput("gcoreTable"),
                 
                 div(class = "section-hdr", "3. Response probabilities (point estimates used in the base case)"),
                 DTOutput("probsTable"),
                 
                 div(class = "section-hdr", "4. Dose table by responder definition (average dose per patient per 14 days)"),
                 DTOutput("doseGridTable"),
                 
                 div(class = "section-hdr", "5. Standard-of-care non-responder dose (average dose per patient per 14 days)"),
                 DTOutput("socnrTable"),
                 
                 div(class = "section-hdr", "6. Uncertainty ranges for OWSA (mean, 95% CI, Table C13)"),
                 DTOutput("rangesTable")
        ),
        
        tabPanel("Base case",
                 br(),
                 fluidRow(
                   column(4, div(class = "value-box", div(class = "lbl", "gammaCore + SoC"), div(class = "big", textOutput("bc_gc", inline = TRUE)))),
                   column(4, div(class = "value-box", div(class = "lbl", "SoC alone"), div(class = "big", textOutput("bc_soc", inline = TRUE)))),
                   column(4, div(class = "value-box save-box", div(class = "lbl", "Difference"), div(class = "big", textOutput("bc_diff", inline = TRUE))))
                 ),
                 fluidRow(
                   column(6, h4("Cost by category"), plotOutput("catPlot", height = 320)),
                   column(6, h4("Responder % by month"), plotOutput("tracePlot", height = 320))
                 ),
                 fluidRow(
                   column(6, DTOutput("catTable")),
                   column(6, DTOutput("cycleTable"))
                 )
        ),
        
        tabPanel("Scenarios (all 30 combinations)",
                 br(),
                 p("All 5 responder-definition \u00d7 3 loss-type \u00d7 2 non-responder-use combinations, computed with the current Inputs-tab settings. The row matching your current sidebar selection is highlighted."),
                 plotOutput("scenarioPlot", height = 380),
                 DTOutput("scenarioTable")
        ),
        
        tabPanel("OWSA (tornado)",
                 br(),
                 p("One-way sensitivity analysis on the ten parameters listed in the Inputs tab's uncertainty-range table. ",
                   "Always uses the 50% responder definition (the only definition with published ranges); loss-type, non-responder-use, and pricing follow your sidebar/Inputs selections."),
                 plotOutput("tornadoPlot", height = 420),
                 DTOutput("owsaTable")
        ),
        
        tabPanel("About & validation",
                 br(),
                 htmlOutput("about_html")
        )
      )
    )
  )
)

## -----------------------------------------------------------------------------
## SERVER
## -----------------------------------------------------------------------------

server <- function(input, output, session) {
  
  ## ---- reactive editable data stores (the Shiny equivalent of Excel's blue cells) ----
  unitCosts <- reactiveVal(default_unit_costs)
  gcoreVals <- reactiveVal(default_gcore)
  probVals  <- reactiveVal(default_probs)
  doseGrid  <- reactiveVal(default_dose_grid)
  socnrVals <- reactiveVal(default_socnr)
  rangesVals<- reactiveVal(default_ranges)
  
  ## ---- editable DT renderers + cell-edit observers ----
  render_editable <- function(df, editable_cols) {
    datatable(df, rownames = FALSE, options = list(dom = "t", paging = FALSE, ordering = FALSE),
              editable = list(target = "cell", disable = list(columns = setdiff(seq_len(ncol(df)) - 1, editable_cols))))
  }
  
  output$unitCostTable <- renderDT(render_editable(unitCosts()[, c("Parameter", "Value")], 1))
  observeEvent(input$unitCostTable_cell_edit, {
    d <- unitCosts(); e <- input$unitCostTable_cell_edit
    d[e$row, "Value"] <- as.numeric(e$value); unitCosts(d)
  })
  
  output$gcoreTable <- renderDT(render_editable(gcoreVals()[, c("Parameter", "Value")], 1))
  observeEvent(input$gcoreTable_cell_edit, {
    d <- gcoreVals(); e <- input$gcoreTable_cell_edit
    d[e$row, "Value"] <- as.numeric(e$value); gcoreVals(d)
  })
  
  output$probsTable <- renderDT(render_editable(probVals()[, c("Parameter", "Value")], 1))
  observeEvent(input$probsTable_cell_edit, {
    d <- probVals(); e <- input$probsTable_cell_edit
    d[e$row, "Value"] <- as.numeric(e$value); probVals(d)
  })
  
  output$doseGridTable <- renderDT(render_editable(doseGrid(), 1:10))
  observeEvent(input$doseGridTable_cell_edit, {
    d <- doseGrid(); e <- input$doseGridTable_cell_edit
    d[e$row, e$col + 1] <- as.numeric(e$value); doseGrid(d)
  })
  
  output$socnrTable <- renderDT(render_editable(socnrVals(), 0:2))
  observeEvent(input$socnrTable_cell_edit, {
    d <- socnrVals(); e <- input$socnrTable_cell_edit
    d[e$row, e$col + 1] <- as.numeric(e$value); socnrVals(d)
  })
  
  output$rangesTable <- renderDT(render_editable(rangesVals()[, c("Parameter", "Mean", "Lower", "Upper")], 1:3))
  observeEvent(input$rangesTable_cell_edit, {
    d <- rangesVals(); e <- input$rangesTable_cell_edit
    col_name <- c("Parameter", "Mean", "Lower", "Upper")[e$col + 1]
    d[e$row, col_name] <- as.numeric(e$value); rangesVals(d)
  })
  
  observeEvent(input$reset_all, {
    unitCosts(default_unit_costs); gcoreVals(default_gcore); probVals(default_probs)
    doseGrid(default_dose_grid); socnrVals(default_socnr); rangesVals(default_ranges)
  })
  
  ## ---- assemble the current parameter state into the shapes run_model() expects ----
  currentDoseTable <- reactive({
    g <- doseGrid()
    setNames(lapply(seq_len(nrow(g)), function(i) {
      list(resp = c(zolm = g[i, "Responder zolm"], suma = g[i, "Responder suma"], oxy = g[i, "Responder oxy"]),
           onTx = c(zolm = g[i, "NonRespOnTx zolm"], suma = g[i, "NonRespOnTx suma"], oxy = g[i, "NonRespOnTx oxy"]),
           base = c(zolm = g[i, "NonRespBaseline zolm"], suma = g[i, "NonRespBaseline suma"], oxy = g[i, "NonRespBaseline oxy"]),
           p = g[i, "P(response) gCore"])
    }), g$Definition)
  })
  currentSocNr <- reactive({ v <- socnrVals(); c(zolm = v$zolm[1], suma = v$suma[1], oxy = v$oxy[1]) })
  uc <- function(key) unitCosts()$Value[unitCosts()$key == key]
  gc_ <- function(key) gcoreVals()$Value[gcoreVals()$key == key]
  pv <- function(key) probVals()$Value[probVals()$key == key]
  
  ## common argument list shared by every run_model() call
  commonArgs <- reactive({
    list(gcore_first_n_free = gc_("gc_free"), gcore_cost_3m = gc_("gc_price"),
         p_resp_soc_use = pv("p_resp_soc"), p_loss_use = pv("p_loss"),
         dose_table_arg = currentDoseTable(), soc_nr_arg = currentSocNr(),
         frac_sc = uc("frac_sc"), frac_o2_portable = uc("frac_o2p"),
         cost_zolm = uc("cost_zolm"), cost_suma_sc = uc("cost_suma_sc"), cost_suma_nasal = uc("cost_suma_nasal"),
         cost_o2_static = uc("cost_o2_static"), cost_o2_portable = uc("cost_o2_portable"))
  })
  
  base <- reactive({
    do.call(run_model, c(list(responder_def = input$responder_def, loss_type = input$loss_type,
                              nonresp_use = input$nonresp_use, soc_no_response = input$soc_no_response),
                         commonArgs()))
  })
  
  output$bc_gc   <- renderText(sprintf("\u00a3%s", format(round(base()$total_gc, 2), big.mark = ",", nsmall = 2)))
  output$bc_soc  <- renderText(sprintf("\u00a3%s", format(round(base()$total_soc, 2), big.mark = ",", nsmall = 2)))
  output$bc_diff <- renderText(sprintf("\u00a3%s", format(round(base()$diff, 2), big.mark = ",", nsmall = 2)))
  
  catData <- reactive({
    do.call(disaggregate, c(list(model_result = base(), responder_def = input$responder_def, nonresp_use = input$nonresp_use),
                            commonArgs()))
  })
  
  output$catPlot <- renderPlot({
    d <- catData(); d <- d[d$Item != "Total", ]
    dl <- data.frame(Item = rep(d$Item, 2), Arm = rep(c("gammaCore + SoC", "SoC alone"), each = nrow(d)), Cost = c(d$gammaCore_SoC, d$SoC_alone))
    dl$Item <- factor(dl$Item, levels = d$Item)
    ggplot(dl, aes(x = Item, y = Cost, fill = Arm)) +
      geom_col(position = position_dodge(width = 0.7), width = 0.65) +
      geom_text(aes(label = sprintf("\u00a3%.0f", Cost)), position = position_dodge(width = 0.7), vjust = -0.4, size = 3.4) +
      scale_fill_manual(values = c("gammaCore + SoC" = "#1F4E78", "SoC alone" = "#3D6E8C")) +
      labs(x = NULL, y = "Annual cost per patient (\u00a3)", fill = NULL) +
      theme_minimal(base_size = 12) + theme(legend.position = "bottom")
  })
  
  output$tracePlot <- renderPlot({
    cyc <- base()$cycles
    dl <- data.frame(cycle = rep(cyc$cycle, 2), Arm = rep(c("gammaCore + SoC", "SoC alone"), each = 12), Responder = c(cyc$r_gc, cyc$r_soc))
    ggplot(dl, aes(x = cycle, y = Responder, color = Arm)) +
      geom_line(linewidth = 1) + geom_point(size = 2) +
      scale_color_manual(values = c("gammaCore + SoC" = "#1F4E78", "SoC alone" = "#3D6E8C")) +
      scale_x_continuous(breaks = 1:12) + scale_y_continuous(labels = percent) +
      labs(x = "Month", y = "Responder %", color = NULL) +
      theme_minimal(base_size = 12) + theme(legend.position = "bottom")
  })
  
  output$catTable <- renderDT({
    datatable(catData(), rownames = FALSE, options = list(dom = "t", paging = FALSE),
              colnames = c("Item", "gammaCore + SoC (\u00a3)", "SoC alone (\u00a3)")) |>
      formatCurrency(c("gammaCore_SoC", "SoC_alone"), currency = "\u00a3")
  })
  
  output$cycleTable <- renderDT({
    cyc <- base()$cycles
    show <- data.frame(Cycle = cyc$cycle, `gCore resp %` = sprintf("%.1f%%", cyc$r_gc * 100),
                       `gCore cost` = sprintf("\u00a3%.2f", cyc$cost_gc), `SoC resp %` = sprintf("%.1f%%", cyc$r_soc * 100),
                       `SoC cost` = sprintf("\u00a3%.2f", cyc$cost_soc), check.names = FALSE)
    datatable(show, rownames = FALSE, options = list(dom = "t", paging = FALSE))
  })
  
  scenarioData <- reactive({
    grid <- expand.grid(responder_def = names(def_labels), loss_type = names(loss_labels), nonresp_use = names(nonresp_labels), stringsAsFactors = FALSE)
    res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
      g <- grid[i, ]
      r <- do.call(run_model, c(list(responder_def = g$responder_def, loss_type = g$loss_type, nonresp_use = g$nonresp_use), commonArgs()))
      data.frame(responder_def = g$responder_def, loss_type = g$loss_type, nonresp_use = g$nonresp_use,
                 gammaCore_SoC = round(r$total_gc, 0), SoC = round(r$total_soc, 0), Difference = round(r$diff, 0))
    }))
    res$is_current <- res$responder_def == input$responder_def & res$loss_type == input$loss_type & res$nonresp_use == input$nonresp_use
    res
  })
  
  output$scenarioPlot <- renderPlot({
    d <- scenarioData()
    d$label <- paste(def_labels[d$responder_def], "|", loss_labels[d$loss_type])
    d$label <- factor(d$label, levels = unique(d$label))
    ggplot(d, aes(x = label, y = Difference, fill = nonresp_use, alpha = is_current)) +
      geom_col(position = position_dodge(width = 0.75), width = 0.65) +
      scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.45), guide = "none") +
      scale_fill_manual(values = c("soc" = "#1F4E78", "baseline" = "#C1443C"),
                        labels = c("soc" = "Revert to SoC use", "baseline" = "Revert to baseline use"), name = NULL) +
      coord_flip() +
      labs(x = NULL, y = "Incremental cost, gammaCore + SoC vs. SoC (\u00a3)") +
      theme_minimal(base_size = 11) + theme(legend.position = "bottom")
  })
  
  output$scenarioTable <- renderDT({
    d <- scenarioData()
    d$responder_def <- def_labels[d$responder_def]; d$loss_type <- loss_labels[d$loss_type]; d$nonresp_use <- nonresp_labels[d$nonresp_use]
    d$Current <- ifelse(d$is_current, "\u2605 current", "")
    show <- d[, c("Current", "responder_def", "loss_type", "nonresp_use", "gammaCore_SoC", "SoC", "Difference")]
    names(show) <- c("Current?", "Responder definition", "Loss assumption", "Non-responder use", "gammaCore + SoC (\u00a3)", "SoC (\u00a3)", "Difference (\u00a3)")
    datatable(show, rownames = FALSE, options = list(pageLength = 10))
  }, server = FALSE)
  
  run_owsa_point <- function(key, value) {
    pm <- param_map[[key]]
    args <- c(list(responder_def = "50", loss_type = input$loss_type, nonresp_use = input$nonresp_use), commonArgs())
    if (pm$type == "scalar") {
      if (pm$target == "p_resp_gc") args$p_resp_gc_override <- value
      else if (pm$target == "p_resp_soc") args$p_resp_soc_use <- value
      else if (pm$target == "p_loss") args$p_loss_use <- value
      else if (pm$target == "frac_sc") args$frac_sc <- value
      else if (pm$target == "frac_o2p") args$frac_o2_portable <- value
      else if (pm$target == "o2_static") args$cost_o2_static <- value
      else if (pm$target == "o2_portable") args$cost_o2_portable <- value
    } else {
      args$dose_override <- dose_override_helper(currentDoseTable(), currentSocNr(), pm$which, pm$drug, value)
    }
    do.call(run_model, args)$diff
  }
  
  owsaData <- reactive({
    rv <- rangesVals()
    res <- do.call(rbind, lapply(seq_len(nrow(rv)), function(i) {
      lo <- run_owsa_point(rv$key[i], rv$Lower[i]); hi <- run_owsa_point(rv$key[i], rv$Upper[i])
      data.frame(Parameter = rv$Parameter[i], Base = rv$Mean[i], Cost_at_lower = round(lo, 0), Cost_at_upper = round(hi, 0))
    }))
    res$Range <- abs(res$Cost_at_upper - res$Cost_at_lower)
    res[order(res$Range), ]
  })
  
  output$tornadoPlot <- renderPlot({
    d <- owsaData()
    d$Parameter <- factor(d$Parameter, levels = d$Parameter)
    base_diff <- base()$diff
    ggplot(d) +
      geom_segment(aes(x = Parameter, xend = Parameter, y = Cost_at_lower, yend = Cost_at_upper), linewidth = 6, color = "#3D6E8C") +
      geom_hline(yintercept = base_diff, linetype = "dashed", color = "#E9A23B") +
      coord_flip() +
      labs(x = NULL, y = "Incremental cost, gammaCore + SoC vs. SoC (\u00a3)",
           caption = "Dashed line = current base-case difference (50% definition, current loss/non-responder-use/pricing settings).") +
      theme_minimal(base_size = 11)
  })
  
  output$owsaTable <- renderDT({
    d <- owsaData()[order(-owsaData()$Range), ]
    names(d) <- c("Parameter", "Base value", "Cost diff. \u2014 lower bound (\u00a3)", "Cost diff. \u2014 upper bound (\u00a3)", "Range (\u00a3)")
    datatable(d, rownames = FALSE, options = list(dom = "t", paging = FALSE))
  })
  
  output$about_html <- renderUI({
    HTML(paste0(
      "<h4>Model overview</h4>",
      "<p>This is a 2-state Markov <b>cost-minimisation</b> model (no QALYs \u2014 the sponsor did not model health utilities, only NHS resource costs). ",
      "States are Responder and Non-responder, defined as achieving a threshold reduction in weekly cluster-headache attack frequency vs. baseline. ",
      "Two arms are compared: gammaCore (nVNS) added to standard-of-care abortive medication, vs. standard of care (SoC) alone. ",
      "Cycle length is 1 month, time horizon is 12 months, perspective is NHS, and no discounting is applied.</p>",
      "<h4>Editable inputs</h4>",
      "<p>Every parameter used by the model is editable on the <b>Inputs</b> tab \u2014 unit costs, gammaCore pricing, response probabilities, ",
      "the dose table for all five responder definitions, the SoC non-responder dose, and the uncertainty ranges used by OWSA. This mirrors ",
      "the companion Excel workbook, where the same values sit in blue (point estimate) cells on the Inputs sheet. ",
      "Use the 'Reset all inputs to published values' button in the sidebar to restore the Table C10/C13 defaults at any time.</p>",
      "<h4>Validation against the published results</h4>",
      "<p>Base case (50% definition, single-loss, revert-to-SoC, all inputs at published defaults): this model returns <b>\u00a33,458.96</b> for ",
      "gammaCore + SoC and <b>\u00a33,913.04</b> for SoC alone (difference <b>\u2212\u00a3454.08</b>), versus the company-reported \u00a33,448.45 / \u00a33,898.86 ",
      "(difference \u2212\u00a3450.42) \u2014 within about 1%. Every alternative scenario and every OWSA parameter reproduces the same sign, ",
      "ranking, and order of magnitude as the published results (Tables C19\u2013C20), generally within 1\u20135%.</p>",
      "<p>The residual gap is not a structural error \u2014 it comes from the sponsor's own tables, which round resource-use inputs to ",
      "2\u20133 significant figures (e.g. \u201c2.50 doses per 14 days\u201d). The true unrounded PREVA post-hoc data is not published in the ",
      "submission, so a small rounding-driven gap of this size is the practical ceiling on replication accuracy achievable from the public submission alone.</p>",
      "<h4>Why there's no PSA / CE plane / CEAC</h4>",
      "<p>The sponsor's submission is a <b>cost-minimisation model with no QALYs or other formal health-outcome measure</b> \u2014 gammaCore was ",
      "submitted to NICE on the basis that it is cost-saving, not cost-effective per QALY gained. A genuine cost-effectiveness plane or CEAC ",
      "would require a fabricated effect measure not present in the original evidence, which risks being mistaken for a NICE-endorsed outcome. ",
      "This app therefore deliberately omits probabilistic sensitivity analysis and sticks to the one-way (OWSA) sensitivity analysis that the ",
      "sponsor's own submission actually reports (Table C13).</p>",
      "<h4>Source</h4>",
      "<p>electroCore. Supporting documentation \u2013 Committee papers: gammaCore for cluster headache [MT323]. Sponsor submission of evidence, ",
      "Section 9 (pp. 132\u2013177); NICE External Assessment Centre report, Section 4.6 (pp. 62\u201374).</p>"
    ))
  })
}

shinyApp(ui, server)