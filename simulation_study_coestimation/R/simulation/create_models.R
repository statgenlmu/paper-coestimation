#!/usr/bin/Rscript --no-site-file --no-save

library(coala)

model_iwm_neutr <- coal_model(c(20, 25)) +
  feat_migration(par_range('m', .2, 5), symmetric = TRUE) +
  feat_pop_merge(par_range('tau', .15, 3.75), 2, 1) +
  feat_size_change(2, population = 1, time = "tau") +
  sumstat_jsfs() +
  sumstat_four_gamete(population = 1) +
  par_range("rho", 0.2, 5) +
  feat_recombination(par_expr(rho * locus_length / 1000)) +
  par_range('theta', 0.2, 5) +
  feat_mutation(par_expr(theta * locus_length / 1000)) +
  locus_averaged(200, 150000)

model_iwm_sel <- model_iwm_neutr +
  par_range("s", 0, 1) +
  feat_selection(strength_A = par_zero_inflation(1000, par_expr(1 - s), FALSE),
                 population = 2, time = 0.03) +
  sumstat_mcmf(population = 2)


model_growth_neutr <- model_iwm_neutr +
  feat_size_change(par_range("q", 2, 50), 2, time = 0) +
  feat_growth(par_expr(log(q)/tau), population = 2, time = 0)

model_growth_sel <- model_iwm_sel +
  feat_size_change(par_range("q", 2, 50), 2, time = 0) +
  feat_growth(par_expr(log(q)/tau), population = 2, time = 0)


model_bn_neutr <- model_iwm_neutr +
  feat_size_change(par_range("q", 0.02, 0.5), 2, time = par_expr(tau - 0.1))

model_bn_sel <- model_iwm_sel +
  feat_size_change(par_range("q", 0.02, 0.5), 2, time = par_expr(tau - 0.1))

save(model_bn_sel, model_bn_neutr,
     model_iwm_sel, model_iwm_neutr,
     model_growth_sel, model_growth_neutr,
     file = "cache/models.Rda")
