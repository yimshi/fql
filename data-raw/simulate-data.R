# Script to generate the sim_microbiome example dataset
# Run this script to regenerate data/sim_microbiome.rda

set.seed(2023)

n <- 200  # samples
n_taxa <- 10  # taxa

# Covariates
group <- factor(rep(c("Control", "Treatment"), each = n / 2))
age <- rnorm(n, mean = 50, sd = 10)

# Simulate OTU counts from negative binomial with varying dispersion
# to demonstrate heteroscedasticity
otu_counts <- matrix(0, nrow = n, ncol = n_taxa)
colnames(otu_counts) <- paste0("OTU", seq_len(n_taxa))

# True effects: first 3 OTUs are differentially abundant
beta0_vals <- c(3, 4, 2.5, 3.5, 4.5, 3, 2, 4, 3.5, 3)
beta_group <- c(0.5, -0.4, 0.3, 0, 0, 0, 0, 0, 0, 0)
beta_age <- c(0.01, -0.005, 0, 0.01, 0, 0, 0, 0, 0, 0)
size_vals <- c(2, 3, 1, 5, 4, 2, 3, 6, 4, 2)

x_group <- as.numeric(group == "Treatment")

for (j in seq_len(n_taxa)) {
  log_mu <- beta0_vals[j] + beta_group[j] * x_group + beta_age[j] * age
  mu <- exp(log_mu)
  otu_counts[, j] <- MASS::rnegbin(n, mu = mu, theta = size_vals[j])
}

# Add some zeros to mimic zero-inflation in rare taxa
for (j in 8:10) {
  zero_idx <- sample(n, size = floor(n * 0.3))
  otu_counts[zero_idx, j] <- 0
}

# Create the dataset
sim_microbiome <- list(
  otu_table = otu_counts,
  sample_data = data.frame(
    group = group,
    age = age,
    total_reads = rowSums(otu_counts)
  ),
  taxa_info = data.frame(
    taxon = colnames(otu_counts),
    differential = c(rep(TRUE, 3), rep(FALSE, 7)),
    true_effect = beta_group
  )
)

usethis::use_data(sim_microbiome, overwrite = TRUE)
