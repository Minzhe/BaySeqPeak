###########################################################################
###                           peak.calling.R                            ###
###########################################################################

# The following script is used to fit MeRIP-Seq data to the Bayesian negative binomial 
# hidden Markov model proposed in the submitted manuscript titled "A Bayesian 
# Hierarchical Model for Analyzing Methylated RNA Immunoprecipitation Sequencing Data."

# Before running the following code, please first load MeRIP-Seq data. The necessary inputs should be 
# (1) a n-by-W count matrix Y, where n is the number of samples and W is the number of bins
# (2) a n-dimensional logical vector c, indiating the allocation of IP samples.
# (3) size factor obtained from function size.factors.estimator()


# Note that the notations in the code and data follow the notations in the manuscript.



############################  call peaks  ##############################
peak.calling <- function(Y, c, s, plot = FALSE) {
      
      # ------------------------ Check for data format and size ------------------------ #
      if (ncol(Y) != length(c)) stop("Wrong condition vector c size!")
      if (ncol(Y) != length(s)) stop("Wrong size factor size!")
      
      # ------------------------ Check for zero count ------------------------ #
      est <- data.frame(fc = rep(NA, nrow(Y)), ppi = rep(NA, nrow(Y)))
      if (sum(rowSums(Y) >= 10) == 0) {
            cat("Too many zero count\n")
            return(est)
      }
      
      # ------------------------ Load data ------------------------ #
      Y <- as.matrix(t(Y))
      n <- dim(Y)[1]
      W <- dim(Y)[2]
      
      # ------------------------ Load algorithm setting ------------------------ #
      iter <- 10 * W
      upp_eta_1 <- 1          # Load hyperparameters, please adjust the other hyperparameters in "zinb.hmm.cpp"
      z_start <- rbinom(W, 1, 0.5)  # Load initial configuration
      
      # ------------------------ Implement MCMC algorithm ------------------------ #
      g <- colSums(Y)
      start_time <- proc.time()
      res <- tryCatch({
            M <- zinb_hmm(Y, c, s, g, upp_eta_1, z_start, iter)
            est <- data.frame(fc = M$delta, ppi = M$z[2,])
            # The MCMC outputs are stored in M
            # $d_0:        the mean of MCMC samples of d_0
            # $delta:         the mean of MCMC samples of delta
            # $phi:           the mean of MCMC samples of phi
            # $pi:            the mean of MCMC samples of pi
            # $H:             the marginal posterior probability of inclusion (PPI) of H
            # $accept_d_0:    the acceptance rate for updating d_0
            # $accept_delta:  the acceptance rate for updating delta
            # $accept_phi:    the acceptance rate for updating phi
            # $A:             the mean of MCMC samples of A
            # $mu:            the mean of MCMC samples of mu
            # $sigma2:        the mean of MCMC samples of sigma2
            # $z:             the marginal posterior probability of inclusion (PPI) of z
            # $z_sum:         the total number of selected methylated bins
      }, error = function(e) {
            cat("Error!\n")
            return(est)
      })
      end_time <- proc.time()
      time <- end_time - start_time
      cat(paste0("Runtime = ", round(time[3], 1), "s\n"))
      
      if (plot) {
            plot(M$z_sum, type = "l", xlab = "Iterations", ylab = "Number of selected methylated bins")
            plot(M$z[2,], type = "h", xlab = "Bins", ylab = "Posterior probabilities of inclusion") 
      }
      
      return(est)
}




############################  function  #################################

# The function of estimating size factors from sequencing data
# Y: a n-by-W count matrix, where n is the number of samples and W is the number of features
# s.mode: "total", "median", "quantile"

size.factors.estimator <- function(Y, s.mode) {
      Y <- as.matrix(t(Y))
      if (s.mode == "total") {
            s <- rowSums(Y)/sum(Y)
      } else if (s.mode == "median") {
            geoMeans = apply(t(Y), 1, gm_mean)
            s <- estimateSizeFactorsForMatrix(t(Y), geoMeans = geoMeans)
            s <- s/sum(s)
      } else if (s.mode == "quantile") {
            s <- apply(Y, 1, quantile, 0.75)
            s <- s/sum(s)
      }
      g <- colSums(Y);
      if (sum(g == 0) > 0) {
            g[which(g == 0)] <- 1
      }
      return(s)
}
