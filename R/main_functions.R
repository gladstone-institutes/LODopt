#' simulate_cellCounts_fromTissue
#'
#' Simulate single cell assay counts of cell-types derived from a given tissue for each of nsamp samples
#' starting by simulating tissue level counts and then downsampling to obtain single cell assay counts
#'
#'
#' @param props A vector of proportions of cell-types in the tissue. Values should be positive and add up to 1.
#' @param nsamp An integer representing the number of samples
#' @param size A numeric negative binomial distribution parameter representing the biological coefficient of variation for cell counts in the tissue
#' @param depth Expected number of cells in tissue (default = 1e9)
#' @param change_mean a positive numeric vector of the same size as proportion representing the change in the proportion of individual cell-types. A value of 1 implies no change, a value less than 1 implies depletion while values above 1 represent increase in the number of the associated cell-type
#'
#' @return A list with the following elements:
#'   \itemize{
#'     \item counts: count matrix with rows corresponding to the clusters and nsamp columns
#'     \item observed_log_odds_ratios: observed log odds ratios based on simulated tissue level counts
#'     \item theoretical_log_odds_ratios: theoretical log odds ratios based on provided parameters
#'     \item theoretical_log_absolute_abundance_ratios: theoretical log of absolute abundance ratios based on provided parameters
#'     \item observed_log_absolute_abundance_ratios: observed log of absolute abundance ratios based on simulated tissue level counts
#'   }
#'
#' @examples
#' # Create example data
#' set.seed(123)
#' # No of clusters in single cell dataset
#' K=25
#' nsamp = 30
#' alpha = 10^runif(K, min=log10(0.5), max = log10(10))
#' p <- dirmult::rdirichlet(alpha = alpha) |> sort()
#' p <- p[p > 0.001]
#' p <- p/sum(p)
#' size <- rep(10, length(p))
#' change_mean = rep(1, length(p))
#' change_mean[c(1,3,8,15)] = c(0.2, 2, 0.2, 2)
#' depth = 1e9
#' # Simulate counts
#' counts_res <- simulate_cellCounts_fromTissue(props=p,nsamp=nsamp,depth=depth, size = size, change_mean = change_mean)
#' print(counts_res$counts[1:3,1:3])
#'
#' @import stats
#' @importFrom dirmult rdirichlet
#' @export
#'
simulate_cellCounts_fromTissue <- function(props,
                                         nsamp,
                                         size=NULL,
                                         depth = 1e9,
                                         change_mean){

  numcells.tissue <- matrix(NA, length(props),nsamp)
  # Generate total cell counts for each cluster
  for(c in 1:length(props)) {
    mu0 = round(depth*props[c])
    mu1 = round((change_mean[c])*depth*props[c])
    for(s in 1:(nsamp/2)) {
      #mu <- rnorm(1, mean = mu0, sd = 0.2*mu0)
      numcells.tissue[c,s] <- rnbinom(1,size=size[c],mu=mu0)
    }
    for(s in ((nsamp/2) + 1):nsamp) {
      #mu <- rnorm(1, mean = mu1, sd = 0.2*mu1)
      numcells.tissue[c,s] <- rnbinom(1,size=size[c],mu=mu1)
    }
  }

  grp <- rep(c(0,1), each=nsamp/2)

  total.numcells.tissue <- colSums(numcells.tissue)
  logodds.numcells.tissue <- log(numcells.tissue/(total.numcells.tissue - numcells.tissue))
  log.abundance.numcells.tissue <- log(numcells.tissue)

  theoretical_log_odds_ratios = log(change_mean*(1-props)/(1 - (change_mean*props)))
  observed_log_odds_ratios <- rowMeans(logodds.numcells.tissue[,((nsamp/2) + 1):nsamp]) -
    rowMeans(logodds.numcells.tissue[,1:nsamp/2])


  theoretical_log_absolute_abundance_ratios <- log(change_mean)
  observed_log_absolute_abundance_ratios <- rowMeans(log.abundance.numcells.tissue[,((nsamp/2) + 1):nsamp]) -
    rowMeans(log.abundance.numcells.tissue[,1:nsamp/2])

  ###sampling fraction
  sampling_fraction <- 1/(1e5*runif(nsamp, min = 0.5, max = 1.5))
  numcells.10x <- round(sampling_fraction*total.numcells.tissue)

  #numcells.10x <- pmax(rnorm(nsamp, mean = mean_10x_cells, sd = sd_10x_cells), 1000) %>% round(.)
  counts <- matrix(NA,ncol=nsamp, nrow=length(props))
  rownames(counts) <- paste("c",0:(length(props)-1), sep="")
  colnames(counts) <- paste0("S", 1:nsamp)

  for(s in 1:nsamp){
    #print(c(numcells.10x[s], numcells.tissue[,s]))
    counts[,s] <- rmultinom(1, size=numcells.10x[s], prob=numcells.tissue[,s])

  }

  #print(theoretical_log_odds_ratios)
  res <- list(counts = counts,
              observed_log_odds_ratios = observed_log_odds_ratios,
              theoretical_log_odds_ratios = theoretical_log_odds_ratios,
              theoretical_log_absolute_abundance_ratios = theoretical_log_absolute_abundance_ratios,
              observed_log_absolute_abundance_ratios = observed_log_absolute_abundance_ratios
  )
  return(res)
}


#' estimate_variance_w_args
#'
#' Total (across a specified set of clusters) of the sample variance of the log odds of clusters membership across all samples within each of those clusters
#'
#'
#' @param pars A numeric vector of length equal to no. of samples, representing their corresponding normalization factors
#' @param counts A cell count matrix with number of columns representing the number of samples and number of rows representing the number of cell-types or clusters
#' @param total_cells An integer vector representing the total number of cells per sample
#' @param optim_clusters An integer vector representing the indices of the clusters over which the total variance is to be calculated
#'
#' @return A numeric value representing the sum of variances of log odds of cluster membership across all clusters specified in optim_clusters
#'
#' @examples
#' @export
#'
estimate_variance_w_args <- function(pars, counts, total_cells, optim_clusters) {

  total_var.norm <- vector(mode = "numeric")
  nsamp <- ncol(counts)
  nclust <- nrow(counts)
  log_odds <- matrix(NA, nclust, nsamp)
  pars <- pars/exp(mean(log(pars)))
  for(s in 1:ncol(counts)) {
    if(sum(pars[s]*total_cells[s] - counts[,s] < 0) > 0) {
      print("non-feasible solution")
      print((pars[s]*total_cells[s] - counts[,s]))
    }
    log_odds[,s] <- log((counts[,s]+1)/(pars[s]*total_cells[s] - counts[,s]))
  }
  for(c in optim_clusters) {
    temp_var_estimate <- (1/(nsamp - 1))*t(log_odds[c,])%*%(diag(nsamp) - (1/nsamp)*(rep(1, nsamp))%*%t(rep(1, nsamp)))%*%((log_odds[c,]))
    total_var.norm[c] <- temp_var_estimate[1,1]
  }
  return(sum(total_var.norm[optim_clusters]))
}

#' ineqfun_data
#'
#' Normalized number of total cells per sample
#'
#'
#' @param pars A numeric vector of length equal to no. of samples, representing their corresponding normalization factors
#' @param total_cells An integer vector representing the total number of cells per sample
#'
#' @return A numeric vector representing the normalized number of total cells per sample
#'
#' @examples
#' @export
#'
ineqfun_data <- function(pars, total_cells) {
  norm_pars <- pars/exp(mean(log(pars)))
  ineqres <- (norm_pars*total_cells)
  return(ineqres)
}

#' logodds_optimized_normFactors
#'
#' The primary LODopt method
#'
#'
#' @param cellcomp_se a SummarizedExperiment object that includes
#' \itemize{
#' \item counts: an assay named counts with the number of cells from each sample in each of the clusters. The counts is a data frame with rownames set to the cluster names (without spaces or special characters) and column names representing the sample names
#' \item data for the association analyses assigned to colData as a data frame with row names the same as the column names of the counts data frame
#' \item Five parameters assigned to the metadata slot
#' \itemize{
#' \item modelFormula representing the model to be tested using the variable in the colData
#' \item coef_of_interest_index representing the index of the coefficient of interest to be tested in the model
#' \item reference_levels_of_variables is a list of vectors.Each vector will have two elements - the name of a categorical variable (a column name in colData) and the reference level to set this variable to
#' \item random_seed integer random seed
#' \item unchanged_cluster_indices cluster indices of the clusters known not to change, set to NULL if not known
#' }
#' }
#' @return A list with two elements
#' \itemize{
#' \item results: data frame of results of differential analyses with 4 columns - cluster_id, comparison, estimates (log odds ratio) and estimate_significance (pvalue)
#' \item optim_factors: a numeric vector with the optimal normalization factors for each sample in the data
#' }
#'
#' @examples
#' @export
#' @importFrom magrittr %<>%
#' @importFrom dplyr mutate filter arrange slice all_of starts_with select across group_by summarise
#' @importFrom lme4 glmer
#' @importFrom SummarizedExperiment SummarizedExperiment
logodds_optimized_normFactors <- function(cellcomp_se) {

  set.seed(cellcomp_se@metadata$random_seed)

  clusters <- rownames(cellcomp_se)
  model_formula <- cellcomp_se@metadata$modelFormula
  coef_of_interest <- cellcomp_se@metadata$coef_of_interest_index

  if(is.null(cellcomp_se@metadata$unchanged_cluster_indices))
    optim_clusters <- 1:length(clusters)
  else
    optim_clusters <- cellcomp_se@metadata$unchanged_cluster_indices

  nsamp <- ncol(cellcomp_se)

  pheno_data <- colData(cellcomp_se) %>% as.data.frame() %>% tibble::rownames_to_column("sampleid")
  if(!is.null(cellcomp_se@metadata$reference_levels_of_variables)) {
    reference_levels_of_variables <- cellcomp_se@metadata$reference_levels_of_variables
    for(i in 1:length(reference_levels_of_variables)) {
      pheno_data[[reference_levels_of_variables[[i]][1]]] <- as.factor(pheno_data[[reference_levels_of_variables[[i]][1]]])
      pheno_data[[reference_levels_of_variables[[i]][1]]] <- relevel(pheno_data[[reference_levels_of_variables[[i]][1]]],
                                                                     ref = reference_levels_of_variables[[i]][2])
    }
  }
  counts <- assays(cellcomp_se)$counts

  total_cells = data.frame(total_cells = colSums(counts),
                           colnames(counts))
  colnames(total_cells)[2] <- "sampleid"


  counts %<>% as.data.frame() %>% tibble::rownames_to_column(var = "clusterid")

  counts_long <- counts %>% tidyr::pivot_longer( cols = -dplyr::starts_with("clusterid"), # Specify the columns to pivot
                                                 names_to = "sampleid",     # New column to hold variable names
                                                 values_to = "count"                # New column to hold values
  )

  Total_Cells <- total_cells$total_cells

  counts %<>% tibble::column_to_rownames("clusterid")


  ineqfun <- function(pars) {
    ineqres <- ineqfun_data(pars, Total_Cells)
    return(ineqres)
  }


  found_stable_solution <- FALSE
  estimates <- vector(mode = "numeric")
  estimates_significance <- vector(mode = "numeric")

  old_outlier_clusters <- 1
  count_iter <- 1

  while(!found_stable_solution & count_iter < 10) {
    print(paste0("Running iteration no. ", count_iter))
    count_iter <- count_iter + 1


    counts_long_norm <- merge(counts_long, total_cells)
    counts_long_norm %<>% merge(., pheno_data)

    counts_long_norm[["sampleid"]] <- as.factor(counts_long_norm[["sampleid"]])



    estimate_variance <- function(pars) {
      estimate_variance_w_args(pars, counts, Total_Cells, optim_clusters)
    }



    initial_guess <- runif(nsamp, min = 0.7, max = 1.4)
    infeasible_solution <- TRUE
    while(infeasible_solution) {
      #print(infeasible_solution)
      norm_pars <- initial_guess/exp(mean(log(initial_guess)))
      if(sum((norm_pars*counts_long_norm$total_cells - counts_long_norm$count) < 0) > 0) {
        #print("New initial guess")
        initial_guess <- runif(nsamp, min = 0.7, max = 1.4)
      }else{
        infeasible_solution <- FALSE
      }
    }
    #print(infeasible_solution)
    temp_res <- Rsolnp::solnp(pars = initial_guess,
                              estimate_variance,
                              LB = rep(0.1, nsamp),
                              UB = rep(50, nsamp),
                              ineqfun = ineqfun,
                              ineqLB = matrixStats::colMaxs(counts %>% as.matrix(.)),
                              ineqUB = rep(2*max(Total_Cells), nsamp),
                              control = list(tol = 1e-12,delta = 1e-9, trace = 0),
    )

    optim.pars <- temp_res$pars/exp(mean(log(temp_res$pars)))
    if(temp_res$convergence != 0)
      print(paste0("Warning: Rsolnp optimization did not converge"))

    optim_factor <- data.frame(colnames(counts),
                               optim.norm.factor = optim.pars)
    colnames(optim_factor)[1] <- "sampleid"
    counts_long_norm %<>% merge(., optim_factor)

    total_var.norm <- vector(mode = "numeric")
    temp_count <- 1
    temp_count_cluster <- 1
    comparison <- vector(mode = "character")
    cluster_id <- vector(mode = "character")
    cluster_significance <- vector(mode = "numeric")
    for(clust in clusters) {
      formula1 <- paste0("cbind(count+1, (round(total_cells*optim.norm.factor) - count)) ~ (1|sampleid) + ", model_formula) %>% as.formula()
      formula0 <- paste0("cbind(count+1, (round(total_cells*optim.norm.factor) - count)) ~ (1|sampleid)") %>% as.formula()

      glmerFit <- lme4::glmer(formula = formula1,
                        data = counts_long_norm %>% dplyr::filter(clusterid == clust) ,
                        family = "binomial")

      sglmerFit <- summary(glmerFit)
      estimates[temp_count] <- sglmerFit$coefficients[coef_of_interest,1]

      glmerFit0 <- lme4::glmer(formula = formula0,
                         data = counts_long_norm %>% dplyr::filter(clusterid == clust) ,
                         family = "binomial")

      anova_res <- anova(glmerFit, glmerFit0)

      cluster_significance[temp_count_cluster] <- anova_res$`Pr(>Chisq)`[2]
      estimates_significance[temp_count] <- sglmerFit$coefficients[coef_of_interest,4]
      comparison[temp_count] <- row.names(sglmerFit$coefficients)[-1]
      cluster_id[temp_count] <- rep(clust, 1)
      temp_count <- temp_count + 1
      temp_count_cluster <- temp_count_cluster + 1
    }

    results <- data.frame(cluster_id, comparison, estimates, estimates_significance)

    temp_res <- data.frame(cluster = 1:nrow(counts), cluster_significance) %>% dplyr::arrange(cluster_significance)
    outlier_clusters <- temp_res %>%
      dplyr::filter(cluster_significance < 0.05) %>%
      dplyr::slice(1:min(c(nrow(.), floor(nrow(counts)/2)))) %>%
      .$cluster %>%
      sort()

    if(identical(old_outlier_clusters, outlier_clusters) | !is.null(cellcomp_se@metadata$unchanged_cluster_indices)) {
      found_stable_solution <- TRUE
      print("Found stable solution")
    }else{
      old_outlier_clusters <- outlier_clusters
      optim_clusters <- setdiff(1:nrow(counts), outlier_clusters)
    }

  }
  return(list(results=results, optim_factor = optim_factor))

}
