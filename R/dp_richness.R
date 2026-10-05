# Public R functions for the positive-count species-richness model.
# Source this file; no analysis runs until fit_dp_richness() is called.
# The unchanged, validated Rcpp kernel lives in src/samplers.cpp.

.dp_source_file <- normalizePath(sys.frame(1)$ofile, mustWork = TRUE)
.dp_backend <- local({
  kernel <- file.path(dirname(dirname(.dp_source_file)),
                      "src", "samplers.cpp")
  env <- new.env(parent = globalenv())
  function(compile = TRUE) {
    if (!file.exists(kernel)) stop("Cannot find the bundled sampler: ", kernel)
    if (compile && !exists("dp_mcmc_cpp", envir = env, inherits = FALSE)) {
      if (!requireNamespace("Rcpp", quietly = TRUE))
        stop("Install Rcpp before fitting: install.packages('Rcpp')")
      Rcpp::sourceCpp(kernel, env = env,
                     cacheDir = file.path(tempdir(), "dp-richness-rcpp"),
                     showOutput = FALSE)
    }
    list(env = env, file = kernel)
  }
})
rm(.dp_source_file)

.dp_scalar <- function(x, name, lower = 0, integer = FALSE, inclusive = FALSE) {
  valid <- is.numeric(x) && length(x) == 1L && is.finite(x)
  if (valid) valid <- if (inclusive) x >= lower else x > lower
  if (valid && integer) valid <- x == floor(x) && x <= .Machine$integer.max
  if (!valid) stop(name, " must be a finite ", if (integer) "integer " else "number ",
                   if (inclusive) ">= " else "> ", lower, call. = FALSE)
  invisible(x)
}

# Convert positive counts, or a frequency table, to the sampler's data format.
# pooled_frequency counts species observed at least pooled_from times.
frequency_data <- function(count, frequency = NULL, pooled_from = NULL,
                           pooled_frequency = 0) {
  if (!is.numeric(count) || any(!is.finite(count)) ||
      any(count < 1 | count != floor(count) | count > .Machine$integer.max))
    stop("count must contain positive integers within R's integer range.")
  if (is.null(frequency)) frequency <- rep(1, length(count))
  if (!is.numeric(frequency) || length(frequency) != length(count) ||
      any(!is.finite(frequency)) || any(frequency < 0 | frequency != floor(frequency)))
    stop("frequency must contain nonnegative integers, one per count.")
  .dp_scalar(pooled_frequency, "pooled_frequency", inclusive = TRUE)
  if (pooled_frequency != floor(pooled_frequency))
    stop("pooled_frequency must be an integer.")
  if (is.null(pooled_from)) {
    if (pooled_frequency > 0) stop("Specify pooled_from for a pooled category.")
    pooled_from <- 0L
  } else {
    .dp_scalar(pooled_from, "pooled_from", integer = TRUE)
    if (any(count[frequency > 0] >= pooled_from))
      stop("Exact counts must be below pooled_from; avoid counting a species twice.")
  }
  keep <- frequency > 0
  count <- count[keep]; frequency <- frequency[keep]
  j <- sort(unique(count))
  f <- vapply(j, function(value) sum(frequency[count == value]), numeric(1))
  K <- sum(f) + pooled_frequency
  if (!is.finite(K) || K < 1 || K > 2^53 - 1)
    stop("The total observed species count must be between 1 and 2^53 - 1.")
  structure(list(j = as.integer(j), f = f, cv = as.integer(pooled_from),
                 cf = pooled_frequency, K = K), class = "dp_frequency_data")
}

# Exact ideal-DP result for a Gamma base with fixed parameters and prior 1/Lambda.
# It does not describe the moment threshold of a fixed H-atom approximation.
dp_moment_criterion <- function(concentration, shape, orders = c(1, 2)) {
  .dp_scalar(concentration, "concentration")
  .dp_scalar(shape, "shape")
  if (!is.numeric(orders) || !length(orders) ||
      any(!is.finite(orders) | orders <= 0)) stop("orders must be positive and finite.")
  data.frame(order = orders, threshold = concentration + shape,
             finite = orders < concentration + shape)
}

# Fit the conditional positive-count likelihood under the reciprocal 1/Lambda
# prior. Gamma uses shape/rate, and H is the finite stick-breaking atom count.
# A positive lower_rate truncates the Gamma base below that rate.
# Each seed starts one chain. Errors are retained; failed chains are never dropped.
fit_dp_richness <- function(data, shape = 1, rate = 1, concentration = 1,
                            H = 60L, sweeps = 100000L, burn = 10000L,
                            thin = 10L, seeds = 101:104, lower_rate = 0,
                            mh = TRUE) {
  if (!is.list(data) || !all(c("j", "f", "cv", "cf") %in% names(data)))
    stop("Create data with frequency_data().")
  # Revalidate even if the caller has modified a previously constructed object.
  data <- frequency_data(data$j, data$f,
                         pooled_from = if (identical(data$cv, 0L) ||
                                           identical(data$cv, 0)) NULL else data$cv,
                         pooled_frequency = data$cf)
  .dp_scalar(shape, "shape"); .dp_scalar(rate, "rate")
  .dp_scalar(concentration, "concentration")
  .dp_scalar(lower_rate, "lower_rate", inclusive = TRUE)
  .dp_scalar(H, "H", lower = 2, integer = TRUE, inclusive = TRUE)
  .dp_scalar(sweeps, "sweeps", integer = TRUE)
  .dp_scalar(burn, "burn", integer = TRUE, inclusive = TRUE)
  .dp_scalar(thin, "thin", integer = TRUE)
  if (burn >= sweeps || (sweeps - burn) %% thin != 0 || (sweeps - burn) / thin < 4)
    stop("Require burn < sweeps and at least four retained draws; (sweeps-burn) must divide by thin.")
  if (!is.numeric(seeds) || !length(seeds) || any(!is.finite(seeds)) ||
      any(seeds < 0 | seeds != floor(seeds) | seeds > .Machine$integer.max) ||
      anyDuplicated(seeds)) stop("seeds must be distinct, nonnegative R integers.")
  if (!is.logical(mh) || length(mh) != 1L || is.na(mh)) stop("mh must be TRUE or FALSE.")
  backend <- .dp_backend()
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  chains <- lapply(as.integer(seeds), function(seed) {
    RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection"); set.seed(seed)
    tryCatch({
      draws <- backend$env$dp_mcmc_cpp(
        j = data$j, f = data$f, shape = shape, rate = rate,
        alpha = concentration, H = as.integer(H), sweeps = as.integer(sweeps),
        burn = as.integer(burn), thin = as.integer(thin), truncated = TRUE,
        cv = data$cv, cf = data$cf, lower = lower_rate, mh = mh,
        mh_scans = 1L, swaps = TRUE,
        thresholds = c(1 / data$K, 1 / sqrt(data$K)))
      if (any(!is.finite(draws$N)) || any(draws$N < data$K | draws$N > 2^53 - 1) ||
          any(!is.finite(draws$a)) || any(draws$a < 0 | draws$a > 1))
        stop("Draws outside the supported numerical range.")
      list(status = "complete", seed = seed, draws = draws, error = NULL)
    }, error = function(e) list(status = "failed", seed = seed,
                                draws = NULL, error = conditionMessage(e)))
  })
  settings <- list(shape = shape, rate = rate, concentration = concentration,
                   H = H, sweeps = sweeps, burn = burn, thin = thin,
                   seeds = seeds, lower_rate = lower_rate, mh = mh,
                   population_size_prior = "1/Lambda",
                   representation = "finite stick-breaking; final weight is residual")
  structure(list(data = data, settings = settings, chains = chains,
                 kernel_md5 = unname(tools::md5sum(backend$file)),
                 rng = "L'Ecuyer-CMRG / Inversion / Rejection"),
            class = "dp_richness_fit")
}

print.dp_richness_fit <- function(x, ...) {
  n <- sum(vapply(x$chains, function(ch) ch$status == "complete", logical(1)))
  cat("Finite-DP richness fit:", x$data$K, "observed species; H =", x$settings$H,
      "\n", n, "of", length(x$chains), "chains completed.\n")
  if (n < length(x$chains)) cat("Inspect $chains for errors; pooled summaries are unavailable.\n")
  invisible(x)
}

# Quantiles, rank-normalized diagnostics and quantile Monte Carlo errors.
# No pooled output is returned when any requested chain has failed.
summary.dp_richness_fit <- function(object, ...) {
  if (any(vapply(object$chains, function(ch) ch$status != "complete", logical(1))))
    stop("At least one chain failed. Inspect fit$chains; do not pool only surviving chains.")
  if (!requireNamespace("posterior", quietly = TRUE))
    stop("Install posterior for diagnostics: install.packages('posterior')")
  ans <- lapply(c("a", "N"), function(variable) {
    m <- do.call(cbind, lapply(object$chains, function(ch) ch$draws[[variable]]))
    q <- unname(quantile(as.vector(m), c(.025, .5, .975), type = 7))
    mc <- posterior::mcse_quantile(m, probs = c(.025, .5, .975))
    data.frame(quantity = if (variable == "a") "missing_probability" else "richness",
               q025 = q[1], median = q[2], q975 = q[3],
               rhat = if (ncol(m) > 1) posterior::rhat(m) else NA_real_,
               ess_bulk = posterior::ess_bulk(m), ess_tail = posterior::ess_tail(m),
               mcse_q025 = unname(mc[1]), mcse_median = unname(mc[2]),
               mcse_q975 = unname(mc[3]), chains = ncol(m),
               draws_per_chain = nrow(m))
  })
  do.call(rbind, ans)
}
