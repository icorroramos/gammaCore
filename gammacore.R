## =============================================================================
## gammaCore (nVNS) for chronic cluster headache -- de novo cost analysis
## Replication of the sponsor (electroCore) cost model submitted to NICE
## Medical Technology Guidance MT323 (Sponsor submission of evidence, Section 9;
## cross-checked against the External Assessment Centre (EAC) report, Section 4.6)
##
## Model type: 2-state Markov cost-minimisation model
##   States   : Responder / Non-responder (>=X% reduction in weekly CH attacks)
##   Arms     : gammaCore + standard of care (SoC)  vs.  SoC alone
##   Cycle    : 1 month; Horizon: 12 months; Perspective: NHS; No discounting
##   Outcome  : total abortive-medication + device cost per patient per year
##              (this is a cost-minimisation model, not a cost-utility/QALY
##              model -- no health outcomes/QALYs were modelled by the sponsor)
##
## All parameter values are taken verbatim from Table C10 (variables), Table
## C11/C12 (technology & health-state costs), Table C13 (OWSA ranges), Table
## C14 (multi-way scenarios) and Table C15 (PSA distributions) of the sponsor
## submission. Where the submission rounds inputs to 2-3 significant figures,
## this replication necessarily inherits that rounding -- see the validation
## block at the end, which reproduces the company's own results (Table
## C16-C21) to within approximately 1-5%.
##
## Base R only -- no packages required.
## =============================================================================

set.seed(20260706)

## -----------------------------------------------------------------------------
## 1. FIXED PARAMETERS
## -----------------------------------------------------------------------------

days_per_month <- 365.25 / 12          # 30.4368: converts "doses per 14 days" to monthly doses
n_cycles       <- 12                    # 1-year horizon, 1-month cycles
# No discounting is applied (Table C9): time horizon < 1 year in effect and
# the reference-case 3.5% rate is immaterial at this duration.

# Unit costs (NHS Drug Tariff, March 2019 / East of England Priorities Advisory
# Committee, 2017) -- Table C10
cost_zolmitriptan_nasal   <- 6.08   # GBP per unit (nasal spray -- only formulation costed)
cost_sumatriptan_sc       <- 19.75  # GBP per unit (subcutaneous injection)
cost_sumatriptan_nasal    <- 7.08   # GBP per unit (nasal spray)
cost_oxygen_static        <- 0.56   # GBP per treatment (static cylinder)
cost_oxygen_portable      <- 0.79   # GBP per treatment (portable cylinder)

frac_sumatriptan_sc_base  <- 0.867  # % of sumatriptan treatments that are s.c. (Marin et al. 2018)
frac_oxygen_portable_base <- 0.50   # % of oxygen treatments that are portable (assumption)

gCore_cost_per_3months <- 625       # GBP, list price ex. VAT, electroCore
gCore_free_months      <- 3         # first 93-day activation card supplied free
gCore_cost_per_month   <- gCore_cost_per_3months / 3   # smoothed monthly equivalent (= GBP 208.33)

p_resp_soc <- 0.083  # probability of response (>=50% reduction), SoC arm, month 1 (Gaul et al. 2016)
p_loss     <- 0.310  # probability of discontinued response per month for initial gammaCore
# responders (post-hoc analysis of PREVA, exponential survival fit)

# Dose tables per responder definition (Table C10) -- doses per 14 days.
#   resp   = gammaCore responder (also used for SoC responder -- 50% definition only, per company
#            assumption that SoC month-1 responders consume the same medication as gCore responders)
#   onTx   = gammaCore non-responder, still on treatment (1st 3 months, evaluation period)
#   base   = gammaCore non-responder, medication use observed at PREVA baseline
#            (used only in the "revert to baseline" non-responder scenario)
#   p      = probability of response to gammaCore for that responder definition
dose_table <- list(
  "25"  = list(resp = c(zolm = 0.80, suma = 2.50, oxy = 3.50),
               onTx = c(zolm = 3.80, suma = 5.80, oxy = 16.20),
               base = c(zolm = 4.80, suma = 3.80, oxy = 16.30),
               p = 0.60),
  "40"  = list(resp = c(zolm = 1.00, suma = 3.00, oxy = 2.80),
               onTx = c(zolm = 2.50, suma = 3.70, oxy = 12.20),
               base = c(zolm = 3.30, suma = 2.90, oxy = 19.30),
               p = 0.47),
  "50"  = list(resp = c(zolm = 0.60, suma = 2.50, oxy = 2.20),
               onTx = c(zolm = 2.50, suma = 4.10, oxy = 11.20),
               base = c(zolm = 3.80, suma = 4.50, oxy = 18.60),
               p = 0.40),
  "65"  = list(resp = c(zolm = 0.00, suma = 2.00, oxy = 0.70),
               onTx = c(zolm = 2.20, suma = 3.80, oxy = 9.20),
               base = c(zolm = 3.80, suma = 5.00, oxy = 31.00),
               p = 0.24),
  "50m" = list(resp = c(zolm = 1.60, suma = 2.80, oxy = 6.50),  # Morris et al. 2016 approach: mean use, whole arm
               onTx = c(zolm = 1.30, suma = 7.50, oxy = 10.80),
               base = c(zolm = 1.30, suma = 7.50, oxy = 10.80),
               p = 0.40)
)

dose_soc_nonresponder <- c(zolm = 1.30, suma = 7.50, oxy = 10.80)  # SoC arm non-responder (fixed across definitions)

## -----------------------------------------------------------------------------
## 2. COST HELPER
## -----------------------------------------------------------------------------

#' Monthly abortive-medication cost for a given dose profile (doses per 14 days)
med_cost <- function(dose, frac_sc = frac_sumatriptan_sc_base,
                     frac_o2_portable = frac_oxygen_portable_base,
                     cost_zolm = cost_zolmitriptan_nasal,
                     cost_suma_sc = cost_sumatriptan_sc,
                     cost_suma_nasal = cost_sumatriptan_nasal,
                     cost_o2_static = cost_oxygen_static,
                     cost_o2_portable = cost_oxygen_portable) {
  f <- days_per_month / 14
  suma_unit <- frac_sc * cost_suma_sc + (1 - frac_sc) * cost_suma_nasal
  o2_unit   <- frac_o2_portable * cost_o2_portable + (1 - frac_o2_portable) * cost_o2_static
  unname(dose["zolm"] * f * cost_zolm + dose["suma"] * f * suma_unit + dose["oxy"] * f * o2_unit)
}

## -----------------------------------------------------------------------------
## 3. CORE MARKOV MODEL
## -----------------------------------------------------------------------------

#' Run the gammaCore cost-minimisation model for one scenario
#'
#' @param responder_def  one of "25","40","50","65","50m" (Table C10 definitions)
#' @param loss_type      "none"        - single loss of response between month 1 and 2, then response
#'                                        maintained for the rest of the year (BASE CASE)
#'                        "constant"    - constant monthly probability of response loss (p_loss) throughout
#'                        "diminishing" - loss rate reduced by 10% each month after month 1
#' @param nonresp_use    "soc"      - after gammaCore discontinuation, non-responders revert to
#'                                    SoC-arm medication use (BASE CASE)
#'                        "baseline" - after discontinuation, non-responders revert to their own
#'                                    PREVA baseline medication use
#' @param soc_no_response logical; if TRUE, no patients in the SoC arm are assumed to respond in month 1
#' @param gcore_first_n_free  number of free months of gammaCore (default 3; set to 0 to remove the free trial)
#' @param gcore_cost_3m  list price per 3-month refill (default 625)
#' @param p_resp_gc_override, p_resp_soc_override, p_loss_override, dose_override  optional overrides
#'        used for one-way / probabilistic sensitivity analysis
#' @return list with monthly cycle-level results and annual totals for both arms
run_model <- function(responder_def = "50",
                      loss_type = c("none", "constant", "diminishing"),
                      nonresp_use = c("soc", "baseline"),
                      soc_no_response = FALSE,
                      gcore_first_n_free = gCore_free_months,
                      gcore_cost_3m = gCore_cost_per_3months,
                      p_resp_gc_override = NULL,
                      p_resp_soc_override = NULL,
                      p_loss_override = NULL,
                      dose_override = NULL,     # named list, e.g. list(resp = c(zolm=.., suma=.., oxy=..))
                      frac_sc = frac_sumatriptan_sc_base,
                      frac_o2_portable = frac_oxygen_portable_base,
                      cost_zolm = cost_zolmitriptan_nasal,
                      cost_suma_sc = cost_sumatriptan_sc,
                      cost_suma_nasal = cost_sumatriptan_nasal,
                      cost_o2_static = cost_oxygen_static,
                      cost_o2_portable = cost_oxygen_portable) {
  
  loss_type   <- match.arg(loss_type)
  nonresp_use <- match.arg(nonresp_use)
  
  dd <- dose_table[[responder_def]]
  d_resp <- dd$resp; d_onTx <- dd$onTx; d_base <- dd$base
  p_resp_gc      <- if (!is.null(p_resp_gc_override))  p_resp_gc_override  else dd$p
  p_resp_soc_use <- if (!is.null(p_resp_soc_override)) p_resp_soc_override else p_resp_soc
  p_loss_use     <- if (!is.null(p_loss_override))     p_loss_override     else p_loss
  d_soc_nr <- dose_soc_nonresponder
  
  if (!is.null(dose_override)) {
    if (!is.null(dose_override$resp))   d_resp   <- dose_override$resp
    if (!is.null(dose_override$onTx))   d_onTx   <- dose_override$onTx
    if (!is.null(dose_override$base))   d_base   <- dose_override$base
    if (!is.null(dose_override$soc_nr)) d_soc_nr <- dose_override$soc_nr
  }
  
  d_nonresp_post <- if (nonresp_use == "soc") d_soc_nr else d_base
  gcore_month_cost <- gcore_cost_3m / 3
  
  mc <- function(d) med_cost(d, frac_sc, frac_o2_portable, cost_zolm, cost_suma_sc, cost_suma_nasal, cost_o2_static, cost_o2_portable)
  
  cyc <- data.frame(cycle = 1:n_cycles, r_gc = NA_real_, cost_gc = NA_real_,
                    r_soc = NA_real_, cost_soc = NA_real_)
  
  r_prev <- p_resp_gc
  for (t in 1:n_cycles) {
    
    ## --- gammaCore-arm responder fraction for this cycle ---
    if (t == 1) {
      r_gc <- p_resp_gc
    } else if (loss_type == "none") {
      r_gc <- if (t == 2) p_resp_gc * (1 - p_loss_use) else r_prev
    } else if (loss_type == "constant") {
      r_gc <- r_prev * (1 - p_loss_use)
    } else { # diminishing
      p_t <- p_loss_use * (0.9)^(t - 2)
      r_gc <- r_prev * (1 - p_t)
    }
    r_prev <- r_gc
    
    in_free_trial <- t <= gcore_first_n_free
    
    if (in_free_trial) {
      # evaluation period: gammaCore supplied free; non-responders still on treatment
      # (medication use reflects the reduced, "on treatment" resource-use group)
      cost_resp    <- mc(d_resp)
      cost_nonresp <- mc(d_onTx)
      cost_gc_month <- r_gc * cost_resp + (1 - r_gc) * cost_nonresp
    } else {
      # post-evaluation: responders continue paying for gammaCore; non-responders
      # discontinue gammaCore and revert to comparator-level medication use
      cost_resp    <- mc(d_resp) + gcore_month_cost
      cost_nonresp <- mc(d_nonresp_post)
      cost_gc_month <- r_gc * cost_resp + (1 - r_gc) * cost_nonresp
    }
    
    ## --- SoC arm ---
    if (t == 1 && !soc_no_response) {
      r_soc <- p_resp_soc_use
      cost_resp_soc    <- mc(d_resp)   # conservative: SoC month-1 responders assumed same use as gCore responders
      cost_nonresp_soc <- mc(d_soc_nr)
      cost_soc_month <- r_soc * cost_resp_soc + (1 - r_soc) * cost_nonresp_soc
    } else {
      r_soc <- 0
      cost_soc_month <- mc(d_soc_nr)   # SoC responders revert to non-responder status after month 1
    }
    
    cyc$r_gc[t]     <- r_gc
    cyc$cost_gc[t]  <- cost_gc_month
    cyc$r_soc[t]    <- r_soc
    cyc$cost_soc[t] <- cost_soc_month
  }
  
  list(cycles = cyc,
       total_gc  = sum(cyc$cost_gc),
       total_soc = sum(cyc$cost_soc),
       diff      = sum(cyc$cost_gc) - sum(cyc$cost_soc))
}

## -----------------------------------------------------------------------------
## 4. BASE CASE  (reproduces Table C16: gc=3448.45, soc=3898.86, diff=-450.42)
## -----------------------------------------------------------------------------

base_case <- run_model(responder_def = "50", loss_type = "none", nonresp_use = "soc")

cat("\n================ BASE CASE RESULTS ================\n")
cat(sprintf("gammaCore + SoC : GBP %.2f\n", base_case$total_gc))
cat(sprintf("SoC alone       : GBP %.2f\n", base_case$total_soc))
cat(sprintf("Difference      : GBP %.2f   (company-reported: -450.42)\n", base_case$diff))
cat("=====================================================\n\n")

## -----------------------------------------------------------------------------
## 5. DISAGGREGATION BY COST CATEGORY (cf. Table C17)
## -----------------------------------------------------------------------------

disaggregate <- function(model_result, responder_def = "50", nonresp_use = "soc") {
  dd <- dose_table[[responder_def]]
  d_resp <- dd$resp; d_onTx <- dd$onTx
  d_soc_nr <- dose_soc_nonresponder
  d_base <- dd$base
  d_nonresp_post <- if (nonresp_use == "soc") d_soc_nr else d_base
  f <- days_per_month / 14
  suma_unit <- frac_sumatriptan_sc_base * cost_sumatriptan_sc + (1 - frac_sumatriptan_sc_base) * cost_sumatriptan_nasal
  o2_unit   <- frac_oxygen_portable_base * cost_oxygen_portable + (1 - frac_oxygen_portable_base) * cost_oxygen_static
  
  gc_zolm <- 0; gc_suma <- 0; gc_oxy <- 0; gc_device <- 0
  soc_zolm <- 0; soc_suma <- 0; soc_oxy <- 0
  
  for (t in 1:n_cycles) {
    r_gc <- model_result$cycles$r_gc[t]; r_soc <- model_result$cycles$r_soc[t]
    in_free <- t <= gCore_free_months
    d_nr_gc <- if (in_free) d_onTx else d_nonresp_post
    
    gc_zolm <- gc_zolm + r_gc * d_resp["zolm"] * f * cost_zolmitriptan_nasal + (1 - r_gc) * d_nr_gc["zolm"] * f * cost_zolmitriptan_nasal
    gc_suma <- gc_suma + r_gc * d_resp["suma"] * f * suma_unit + (1 - r_gc) * d_nr_gc["suma"] * f * suma_unit
    gc_oxy  <- gc_oxy  + r_gc * d_resp["oxy"]  * f * o2_unit   + (1 - r_gc) * d_nr_gc["oxy"]  * f * o2_unit
    if (!in_free) gc_device <- gc_device + r_gc * gCore_cost_per_month
    
    if (t == 1) {
      soc_zolm <- soc_zolm + r_soc * d_resp["zolm"] * f * cost_zolmitriptan_nasal + (1 - r_soc) * d_soc_nr["zolm"] * f * cost_zolmitriptan_nasal
      soc_suma <- soc_suma + r_soc * d_resp["suma"] * f * suma_unit + (1 - r_soc) * d_soc_nr["suma"] * f * suma_unit
      soc_oxy  <- soc_oxy  + r_soc * d_resp["oxy"]  * f * o2_unit   + (1 - r_soc) * d_soc_nr["oxy"]  * f * o2_unit
    } else {
      soc_zolm <- soc_zolm + d_soc_nr["zolm"] * f * cost_zolmitriptan_nasal
      soc_suma <- soc_suma + d_soc_nr["suma"] * f * suma_unit
      soc_oxy  <- soc_oxy  + d_soc_nr["oxy"]  * f * o2_unit
    }
  }
  
  out <- rbind(
    GammaCore    = c(gc = gc_device, soc = 0, diff = gc_device - 0),
    Sumatriptan  = c(gc = gc_suma, soc = soc_suma, diff = gc_suma - soc_suma),
    Zolmitriptan = c(gc = gc_zolm, soc = soc_zolm, diff = gc_zolm - soc_zolm),
    Oxygen       = c(gc = gc_oxy, soc = soc_oxy, diff = gc_oxy - soc_oxy)
  )
  rbind(out, Total = colSums(out))
}

disagg <- disaggregate(base_case)
cat("=========== COST BY CATEGORY (cf. Table C17) ===========\n")
print(round(disagg, 2))
cat("Company-reported : GC 517.18 | Sumatriptan 2577.39/3505.53 (-928.13) | Zolmitriptan 206.33/204.85 (+1.48) | Oxygen 147.55/188.49 (-40.95)\n\n")

## -----------------------------------------------------------------------------
## 6. ONE-WAY DETERMINISTIC SENSITIVITY ANALYSIS (cf. Table C13 / C19)
## -----------------------------------------------------------------------------

dose_override_helper <- function(which_dose, drug, value) {
  base_doses <- list(resp = dose_table[["50"]]$resp, onTx = dose_table[["50"]]$onTx,
                     base = dose_table[["50"]]$base, soc_nr = dose_soc_nonresponder)
  d <- base_doses[[which_dose]]
  d[drug] <- value
  ov <- list(); ov[[which_dose]] <- d
  ov
}

owsa_params <- list(
  list(name = "sumatriptan doses/14d - SoC non-responder",            lo = 4.88,  hi = 10.67, base = 7.50,  which = "soc_nr", drug = "suma"),
  list(name = "sumatriptan doses/14d - gCore non-responder (on Tx)",  lo = 1.00,  hi = 9.33,  base = 4.10,  which = "onTx",   drug = "suma"),
  list(name = "sumatriptan doses/14d - gCore responder",              lo = 1.04,  hi = 4.59,  base = 2.50,  which = "resp",   drug = "suma"),
  list(name = "% sumatriptan treatments that are s.c.",               lo = 0.661, hi = 0.982, base = 0.867, which = "frac_sc"),
  list(name = "zolmitriptan doses/14d - gCore non-responder (on Tx)", lo = 0.32,  hi = 6.89,  base = 2.50,  which = "onTx",   drug = "zolm"),
  list(name = "zolmitriptan doses/14d - SoC non-responder",           lo = 0.45,  hi = 2.59,  base = 1.30,  which = "soc_nr", drug = "zolm"),
  list(name = "Probability of response - gCore (50% definition)",     lo = 0.263, hi = 0.545, base = 0.40,  which = "p_resp_gc"),
  list(name = "oxygen doses/14d - SoC non-responder",                 lo = 6.68,  hi = 15.90, base = 10.80, which = "soc_nr", drug = "oxy"),
  list(name = "zolmitriptan doses/14d - gCore responder",             lo = 0.10,  hi = 1.52,  base = 0.60,  which = "resp",   drug = "zolm"),
  list(name = "Probability of discontinued response per month",       lo = 0.162, hi = 0.541, base = 0.31,  which = "p_loss")
)

run_owsa_point <- function(param, value) {
  args <- list(responder_def = "50", loss_type = "none", nonresp_use = "soc")
  if (param$which == "frac_sc") {
    args$frac_sc <- value
  } else if (param$which == "p_resp_gc") {
    args$p_resp_gc_override <- value
  } else if (param$which == "p_loss") {
    args$p_loss_override <- value
  } else {
    args$dose_override <- dose_override_helper(param$which, param$drug, value)
  }
  do.call(run_model, args)$diff
}

owsa_results <- do.call(rbind, lapply(owsa_params, function(p) {
  data.frame(Parameter = p$name, Base = p$base,
             Cost_at_lower = round(run_owsa_point(p, p$lo), 0),
             Cost_at_upper = round(run_owsa_point(p, p$hi), 0))
}))
owsa_results$Range <- abs(owsa_results$Cost_at_upper - owsa_results$Cost_at_lower)
owsa_results <- owsa_results[order(-owsa_results$Range), ]

cat("=========== ONE-WAY SENSITIVITY ANALYSIS (cf. Table C19) ===========\n")
print(owsa_results, row.names = FALSE)
cat(sprintf("(Base-case difference: GBP %.2f; company-reported: -450.42)\n\n", base_case$diff))

## -----------------------------------------------------------------------------
## 7. MULTI-WAY DETERMINISTIC SCENARIO ANALYSIS (cf. Table C14 / C20)
## -----------------------------------------------------------------------------

scenario_grid <- expand.grid(responder_def = c("25", "40", "50m", "50", "65"),
                             loss_type = c("none", "constant", "diminishing"),
                             nonresp_use = c("soc", "baseline"),
                             stringsAsFactors = FALSE)

scenario_results <- do.call(rbind, lapply(seq_len(nrow(scenario_grid)), function(i) {
  g <- scenario_grid[i, ]
  r <- run_model(responder_def = g$responder_def, loss_type = g$loss_type, nonresp_use = g$nonresp_use)
  data.frame(responder_def = g$responder_def, loss_type = g$loss_type, nonresp_use = g$nonresp_use,
             gammaCore_plus_SoC = round(r$total_gc, 0), SoC = round(r$total_soc, 0), Difference = round(r$diff, 0))
}))

cat("=========== MULTI-WAY SCENARIO ANALYSIS (cf. Table C20) ===========\n")
print(scenario_results, row.names = FALSE)
cat("\n")

## -----------------------------------------------------------------------------
## 8. PROBABILISTIC SENSITIVITY ANALYSIS (cf. Table C15 / C21)
##    Beta for probabilities, Gamma for doses/uncertain unit costs, Normal for
##    the response-loss-rate coefficient, per Table C15.
## -----------------------------------------------------------------------------

beta_params <- function(mean, lo, hi) {
  se  <- (hi - lo) / (2 * 1.96)
  var <- min(se^2, mean * (1 - mean) - 1e-6)
  k   <- (mean * (1 - mean) / var) - 1
  c(alpha = max(mean * k, 0.1), beta = max((1 - mean) * k, 0.1))
}
gamma_params <- function(mean, lo, hi) {
  se <- max((hi - lo) / (2 * 1.96), 1e-6)
  c(shape = (mean / se)^2, rate = mean / se^2)
}
rgamma1 <- function(mean, lo, hi) { g <- gamma_params(mean, lo, hi); rgamma(1, shape = g["shape"], rate = g["rate"]) }
rbeta1  <- function(mean, lo, hi) { g <- beta_params(mean, lo, hi); rbeta(1, g["alpha"], g["beta"]) }

n_psa <- 1000

psa_row <- function() {
  p_resp_gc_i  <- rbeta1(0.40, 0.263, 0.545)
  p_resp_soc_i <- rbeta1(0.083, 0.024, 0.175)
  p_loss_i     <- min(max(rnorm(1, 0.31, (0.541 - 0.162) / (2 * 1.96)), 0.001), 0.999)
  
  d_resp   <- c(zolm = rgamma1(0.60, 0.10, 1.52),  suma = rgamma1(2.50, 1.04, 4.59), oxy = rgamma1(2.20, 0.56, 4.94))
  d_onTx   <- c(zolm = rgamma1(2.50, 0.32, 6.89),  suma = rgamma1(4.10, 1.00, 9.33), oxy = rgamma1(11.20, 5.45, 18.98))
  d_base   <- c(zolm = rgamma1(3.80, 0.55, 10.13), suma = rgamma1(4.50, 2.07, 7.85), oxy = rgamma1(18.60, 9.35, 30.98))
  d_soc_nr <- c(zolm = rgamma1(1.30, 0.45, 2.59),  suma = rgamma1(7.50, 4.88, 10.67), oxy = rgamma1(10.80, 6.68, 15.90))
  
  frac_sc_i <- rbeta1(0.867, 0.661, 0.982)
  frac_o2_i <- rbeta1(0.50, 0.001, 0.599)
  
  o2_static_i   <- rgamma1(0.56, 0.50, 0.63)
  o2_portable_i <- rgamma1(0.79, 0.70, 0.88)
  
  r <- run_model(responder_def = "50", loss_type = "none", nonresp_use = "soc",
                 p_resp_gc_override = p_resp_gc_i, p_resp_soc_override = p_resp_soc_i,
                 p_loss_override = p_loss_i,
                 dose_override = list(resp = d_resp, onTx = d_onTx, base = d_base, soc_nr = d_soc_nr),
                 frac_sc = frac_sc_i, frac_o2_portable = frac_o2_i,
                 cost_o2_static = o2_static_i, cost_o2_portable = o2_portable_i)
  c(gc = r$total_gc, soc = r$total_soc, diff = r$diff)
}

psa_out <- as.data.frame(t(replicate(n_psa, psa_row())))

cat("=========== PROBABILISTIC SENSITIVITY ANALYSIS (cf. Table C21) ===========\n")
cat(sprintf("Mean gammaCore + SoC : GBP %.2f  (company-reported: 3427.43)\n", mean(psa_out$gc)))
cat(sprintf("Mean SoC             : GBP %.2f  (company-reported: 3864.13)\n", mean(psa_out$soc)))
cat(sprintf("Mean difference      : GBP %.2f  (company-reported: -436.70)\n", mean(psa_out$diff)))
cat(sprintf("%% of simulations cost-saving for gammaCore: %.1f%%  (company-reported: 88.5%%)\n", 100 * mean(psa_out$diff < 0)))
cat("============================================================================\n")

## -----------------------------------------------------------------------------
## 9. SAVE OUTPUTS
## -----------------------------------------------------------------------------

dir.create("/mnt/user-data/outputs", showWarnings = FALSE, recursive = TRUE)
write.csv(scenario_results, "/mnt/user-data/outputs/gammacore_scenario_results.csv", row.names = FALSE)
write.csv(owsa_results, "/mnt/user-data/outputs/gammacore_owsa_results.csv", row.names = FALSE)
write.csv(psa_out, "/mnt/user-data/outputs/gammacore_psa_draws.csv", row.names = FALSE)

cat("\nOutputs written to /mnt/user-data/outputs/\n")