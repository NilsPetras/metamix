methods::setClass("metamix",
         representation(
           data = "list", # raw data underlying the model
           model = "list", # model definition; might change to call instead
           estimates = "list", # model parameters including SEs, CIs, chisqs, and p-values
           model_fit_test = "list", # ks-test results
           theoretical_distribution = "list" # nrep x n1 sized matrix of bootstrapped t-values and publication status
         ))