# Construct the ECDM model object
new_ecdm_model = function(model_mcmc, details = list()) {

  structure(
    list(
      # Iterates over each parameter and obtains summary information
      estimates = lapply(model_mcmc, summarize_model),
      # Store the MCMC chain results for diagnostics
      chain = model_mcmc,
      # Provide overarching details on model construction
      details = details
      ),
    class = c("exploratory_model")
  )
}
