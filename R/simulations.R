# Source R/dp_richness.R first. Dependencies: Rcpp, posterior; C++11 compiler.
# Gamma distributions use shape/rate. All random-number states are restored.
# These reusable experiments use the manuscript sampler; their default small
# settings are demonstrations, not the manuscript's saved production results.

simulate_species_survey <- function(Lambda, truth = "gamma", shape = 1,
                                    rate = 0.5, seed = 1L) {
  .dp_scalar(Lambda, "Lambda")
  .dp_scalar(shape, "shape"); .dp_scalar(rate, "rate")
  .dp_scalar(seed, "seed", integer = TRUE, inclusive = TRUE)
  truth <- match.arg(truth, c("gamma", "two_point", "rare_group"))
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection"); set.seed(seed)
  N <- rpois(1, Lambda)
  rates <- switch(truth,
    gamma = rgamma(N, shape = shape, rate = rate),
    two_point = sample(c(0.5, 3.5), N, replace = TRUE),
    rare_group = {
      rare <- runif(N) < 0.2
      x <- rgamma(N, shape = 2, rate = 0.8); x[rare] <- 0.05; x
    })
  counts <- rpois(N, rates)
  positive <- counts[counts > 0]
  a_true <- switch(truth,
    gamma = (rate / (rate + 1))^shape,
    two_point = mean(exp(-c(0.5, 3.5))),
    rare_group = 0.2 * exp(-0.05) + 0.8 * (0.8 / 1.8)^2)
  list(N = N, K = length(positive), unseen = sum(counts == 0),
       data = if (length(positive)) frequency_data(positive) else NULL,
       a_true = a_true, truth = truth, Lambda = Lambda, shape = shape,
       rate = rate, seed = seed)
}

# Use the same simulated population for each prior, with distinct chain seeds.
# Failed fits and empty surveys remain explicit rows. No surviving-only pooling.
# Optionally retain complete fit objects using fit_directory.
simulate_richness_study <- function(
    design = data.frame(truth = "gamma", Lambda = 100,
                        shape = c(0.1, 1), rate = c(0.05, 0.5)),
    priors = data.frame(label = c("Ga(1,1)", "Ga(0.1,0.1)"),
                        shape = c(1, 0.1), rate = c(1, 0.1), concentration = 1),
    replicates = 2L, chains = 2L, H = 20L,
    sweeps = 4000L, burn = 1000L, thin = 5L,
    seed = 20261005L, lower_rate = 0, fit_directory = NULL) {
  stopifnot(is.data.frame(design), nrow(design) > 0,
            all(c("truth", "Lambda", "shape", "rate") %in% names(design)),
            is.data.frame(priors), nrow(priors) > 0,
            all(c("label", "shape", "rate", "concentration") %in% names(priors)),
            !anyNA(priors$label), !anyDuplicated(priors$label))
  .dp_scalar(replicates, "replicates", integer = TRUE)
  .dp_scalar(chains, "chains", integer = TRUE)
  .dp_scalar(seed, "seed", integer = TRUE, inclusive = TRUE)
  .dp_scalar(H, "H", lower = 2, integer = TRUE, inclusive = TRUE)
  .dp_scalar(sweeps, "sweeps", integer = TRUE)
  .dp_scalar(burn, "burn", integer = TRUE, inclusive = TRUE)
  .dp_scalar(thin, "thin", integer = TRUE)
  .dp_scalar(lower_rate, "lower_rate", inclusive = TRUE)
  if (burn >= sweeps || (sweeps - burn) %% thin || (sweeps - burn) / thin < 4)
    stop("Invalid sweep, burn-in or thinning settings.")
  if (anyNA(design) || any(!design$truth %in% c("gamma", "two_point", "rare_group")))
    stop("Invalid simulation design.")
  for (nm in c("Lambda", "shape", "rate"))
    for (x in design[[nm]]) .dp_scalar(x, paste("design", nm))
  for (nm in c("shape", "rate", "concentration"))
    for (x in priors[[nm]]) .dp_scalar(x, paste("prior", nm))
  stride <- 1 + nrow(priors) * chains
  if (seed + nrow(design) * replicates * stride > .Machine$integer.max)
    stop("Seed sequence would exceed R's integer range.")
  if (!requireNamespace("posterior", quietly = TRUE)) stop("Install the posterior package.")
  .dp_backend()  # Check the compiler and dependency before starting the study.
  if (!is.null(fit_directory)) {
    if (dir.exists(fit_directory) && length(list.files(fit_directory, all.files = TRUE,
                                                     no.. = TRUE)))
      stop("fit_directory must be new or empty to preserve earlier runs.")
    dir.create(fit_directory, recursive = TRUE, showWarnings = FALSE)
  }
  rows <- list(); job <- 0L
  for (i in seq_len(nrow(design))) for (rep in seq_len(replicates)) {
    job <- job + 1L; survey_seed <- seed + (job - 1) * stride
    survey <- simulate_species_survey(design$Lambda[i], design$truth[i],
                                       design$shape[i], design$rate[i], survey_seed)
    for (p in seq_len(nrow(priors))) {
      seeds <- survey_seed + (p - 1) * chains + seq_len(chains)
      row <- data.frame(design = i, truth = design$truth[i],
        truth_shape = design$shape[i], truth_rate = design$rate[i],
        Lambda = design$Lambda[i], replicate = rep, prior = priors$label[p],
        prior_shape = priors$shape[p], prior_rate = priors$rate[p],
        concentration = priors$concentration[p], H = H, lower_rate = lower_rate,
        sweeps = sweeps, burn = burn, thin = thin, chains = chains,
        survey_seed = survey_seed, chain_seeds = paste(seeds, collapse = ";"),
        N_true = survey$N, observed = survey$K, a_true = survey$a_true,
        status = if (survey$K) "pending" else "no_observations",
        N_lo = NA_real_, N_median = NA_real_, N_hi = NA_real_,
        covered = NA, rhat = NA_real_, ess_bulk = NA_real_, ess_tail = NA_real_,
        mcse_median = NA_real_, mcse_upper = NA_real_, error = "")
      fit <- NULL
      if (survey$K) {
        fit <- fit_dp_richness(survey$data, shape = priors$shape[p],
          rate = priors$rate[p], concentration = priors$concentration[p],
          H = H, sweeps = sweeps, burn = burn, thin = thin,
          seeds = seeds, lower_rate = lower_rate)
        failed <- vapply(fit$chains, function(ch) ch$status != "complete", logical(1))
        if (any(failed)) {
          row$status <- "failed"
          row$error <- paste(vapply(fit$chains[failed], function(ch)
            paste0("seed ", ch$seed, ": ", ch$error), character(1)), collapse = " | ")
        } else {
          s <- summary(fit); s <- s[s$quantity == "richness", ]
          row$status <- "complete"
          row$N_lo <- s$q025; row$N_median <- s$median; row$N_hi <- s$q975
          row$covered <- s$q025 <= survey$N && survey$N <= s$q975
          row$rhat <- s$rhat; row$ess_bulk <- s$ess_bulk; row$ess_tail <- s$ess_tail
          row$mcse_median <- s$mcse_median; row$mcse_upper <- s$mcse_q975
        }
      }
      if (!is.null(fit_directory))
        saveRDS(list(survey = survey, fit = fit, summary = row),
                file.path(fit_directory, sprintf("design_%03d_rep_%03d_prior_%02d.rds", i, rep, p)))
      rows[[length(rows) + 1L]] <- row
    }
  }
  do.call(rbind, rows)
}
