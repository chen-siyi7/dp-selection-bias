
source("R/dp_richness.R")
source("R/simulations.R")
args <- commandArgs(TRUE)
mode <- if (length(args)) args[1] else "demo"
stopifnot(mode %in% c("demo", "full"))
output <- if (length(args) > 1) args[2] else file.path("results", mode)
if (dir.exists(output) && length(list.files(output, all.files = TRUE, no.. = TRUE)))
  stop("Choose a new output directory; existing results will not be overwritten.")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
if (mode == "demo") {
  result <- simulate_richness_study(fit_directory = file.path(output, "fits"))
} else {
  truths <- data.frame(truth = c("gamma", "gamma", "two_point", "rare_group"),
                       shape = c(0.1, 1, 1, 1), rate = c(0.05, 0.5, 1, 1))
  design <- do.call(rbind, lapply(c(500, 2000, 10000), function(Lambda)
    transform(truths, Lambda = Lambda)))
  # shape/rate in the last two design rows are placeholders; those laws are fixed
  # in simulate_species_survey(). Both Gamma truths have mean species rate two.
  result <- simulate_richness_study(design = design, replicates = 20,
    chains = 4, H = 60, sweeps = 100000, burn = 10000, thin = 10,
    fit_directory = file.path(output, "fits"))
}
write.csv(result, file.path(output, "simulation_results.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output, "session_info.txt"))
print(table(result$status))
cat("Saved results, all fit statuses and simulation settings to", output, "\n")
