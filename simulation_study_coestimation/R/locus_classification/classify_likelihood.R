
    # Prepare a jaatha_model and jaatha_data with only one locus
    options <- jaatha_data$get_options()
    model_sel_1 <- model_sel
    model_sel_1$loci <- list()
    model_sel_1 <- model_sel_1 + locus_single(get_locus_length(model_sel))

    jaatha_model_1 <- create_jaatha_model(model_sel_1, test = FALSE)
    stats_1 <- list(jsfs = test_data$data[[run]]$jsfs_per_locus[[1]],
                    mcmf = test_data$data[[run]]$mcmf[1])
    jaatha_data_1 <- create_jaatha_data(stats_1, jaatha_model_1)
    suppressWarnings(jaatha_data_1$set_options(options))

    # Prepare parameters
    est <- jaatha_model_1$get_par_ranges()$normalize(result$param)
    pars_sel <- matrix(est, n_sim, length(est), byrow = TRUE)
    colnames(pars_sel) <- names(result$param)
    pars_sel[, "s"] <- 1
    pars_neutr <- pars_sel
    pars_neutr[, "s"] <- 0

    # Simulate n_sim loci each
    sim_data_sel <- jaatha_model_1$simulate(pars_sel, jaatha_data_1, cores = 1)
    sim_data_neutr <- jaatha_model_1$simulate(pars_neutr, jaatha_data_1, cores = 1)

    # Calculate the likelihood-ratio for each locus
    llh <- sapply(1:get_locus_number(model_sel), function(locus) {
      stats_1 <- list(jsfs = test_data$data[[run]]$jsfs_per_locus[[locus]],
                      mcmf = test_data$data[[run]]$mcmf[locus])
      jaatha_data_1 <- create_jaatha_data(stats_1, jaatha_model_1)

      # Calculate the log-likelihoods
      llh_sel <- estimate_llh(jaatha_model_1, jaatha_data_1, est,
                              sim = n_sim, sim_data = sim_data_sel)
      llh_neutr <- estimate_llh(jaatha_model_1, jaatha_data_1, est,
                                sim = n_sim, sim_data = sim_data_neutr)
      llh_sel$value - llh_neutr$value
    })
