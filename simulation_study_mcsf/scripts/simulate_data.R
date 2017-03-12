library(parallel)
library(coala)

simulate_data <- function(name, models, reps, cores, 
                          pars = numeric(),
                          with_false_demography = FALSE) {
  
  folder <- file.path("data", name)
  dir.create(folder, showWarnings = FALSE, recursive = TRUE)
  seeds <- sample.int(2 ^ 20, reps)
  
  null <- mclapply(seeds, function(seed) {
    set.seed(seed)
    
    if (with_false_demography) {
      pars_neutr <- pars * rnorm(length(pars), 1, 0.2)
    } else {
      pars_neutr <- pars
    }
    
    # Simulate
    sim_data <- lapply(models, simulate, pars = c(pars, s = 0.05))
    sim_data_neutr <- lapply(models, simulate, pars = c(pars_neutr, s = 0))
    
    # Get loci simulated with a sweep
    is_sel_true <- lapply(sim_data, function(x) {
      if ("trees" %in% names(x$cmds)) {
        x$cmds <- x$cmds$trees
      }
      sel_grp <- grepl("-SA 1000", x$cmds)
      loci_in_grp <- sapply(strsplit(x$cmds, " "), function(x) as.numeric(x[3]))
      stopifnot(sum(loci_in_grp) == n_loci)
      do.call(c, lapply(seq(along = sel_grp), function(i) {
        rep(sel_grp[i], loci_in_grp[i])
      }))
    })
    
    # Convert to data.frames
    sim_data <- do.call(rbind, lapply(names(sim_data), function(model) {
      rbind(data.frame(model = model, seed = seed, sel = is_sel_true[[model]], stat = "mcmf", 
                       value_min = NA, value_max = sim_data[[model]]$mcmf),
            data.frame(model = model, seed = seed, sel = is_sel_true[[model]], stat = "ihs", 
                       value_min = sim_data[[model]]$ihh[ , 1], value_max = sim_data[[model]]$ihh[ , 2]),
            data.frame(model = model, seed = seed, sel = is_sel_true[[model]], stat = "omega", 
                       value_min = NA, value_max = sim_data[[model]]$omega)
      )
    }))
    
    sim_data_neutr <- do.call(rbind, lapply(names(sim_data_neutr), function(model) {
      rbind(data.frame(model = model, seed = seed, stat = "mcmf", 
                       value_min = NA, value_max = sim_data_neutr[[model]]$mcmf),
            data.frame(model = model, seed = seed, stat = "omega",
                       value_min = NA, value_max = sim_data_neutr[[model]]$omega))
    }))
    
    save(sim_data, sim_data_neutr,
         file = paste0(folder, "/", seed, ".Rda"))
    NULL
  }, mc.cores = cores, mc.preschedule = FALSE)
}
