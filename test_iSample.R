# Quick R script to 
# Input: h, d
# Output: Exact p-value, write bounds, write correlation matrix ('newbounds/newsig.txt')
# Correlation matrix is random ~0.3

# Call it with Rscript test_iSample.R H D


library(mvtnorm)

# Inputs
args <- commandArgs(trailingOnly=TRUE)
h <- as.numeric(args[1])
d <- as.numeric(args[2])

# Random correlation matrix
start_sig <- matrix(data=0.3, nrow=d, ncol=d)
diag(start_sig) <- 1
temp_samp <- rmvnorm(n=2*d, sigma=start_sig)
random_sig <- cor(temp_samp)

# Explicit inverse of HC to find the p-value bounds
i_vec <- 1:d
HC_p_bounds <- ((2*i_vec+h^2)/d - sqrt((2*i_vec/d+h^2/d)^2 - 4*i_vec^2/d^2 - 4*i_vec^2*h^2/d^3)) /				(2*(1+h^2/d))
HC_z_bounds <- qnorm(1-HC_p_bounds/2)
HC_z_bounds <- sort(HC_z_bounds, decreasing=F)

# qnorm can't handle more precision than 10^-16
HC_z_bounds[which(HC_z_bounds > 8.2)]= 8.2

# Write 
write.table(HC_z_bounds, 'newbounds.txt', append=F, quote=F, row.names=F, col.names=F)
write.table(random_sig[upper.tri(random_sig)], 'newsig.txt', append=F, quote=F, row.names=F, col.names=F)

# Exact p-value
system2(command="./GOF_exact_pvalue", args=c(d, 'newbounds.txt', 
							'newsig.txt', 0))


