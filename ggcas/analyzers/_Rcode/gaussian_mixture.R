library(mclust)

GaussianMixture <- function(
  data,
  n_clusters = 2,
  model_name = "VII",
  n_init = 10,
  max_iter = 1000,
  tol = 1e-6,
  verbose = FALSE
) {
  # Fit the model
  model <- Mclust(data,
                  G = n_clusters,
                  modelNames = model_name,
                  nstart = n_init,
                  maxiter = max_iter,
                  tol = tol,
                  verbose = verbose)
  # Return the model
  return(model)
}
