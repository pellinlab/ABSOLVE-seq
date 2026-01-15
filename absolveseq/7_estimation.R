option_list = list(
  optparse::make_option("--maindir", type = "character", default = NULL, help = "Main working directory (default: current directory)"),
  optparse::make_option("--data_folder", type = "character", default = "target_umi_barcode_table_full/dedudS_filtered_v3"),
  optparse::make_option("--barcode_file", type = "character", default = "barcodeList.csv"),
  optparse::make_option("--baselevel_treat", type = "character", default = "NoEP"),
  optparse::make_option("--number_cores", type = "integer", default = 1L),
  optparse::make_option("--include_nontested_level", type = "logical", default = TRUE),
  optparse::make_option("--include_replicate", type = "logical", default = FALSE),
  optparse::make_option("--number_permutations", type = "integer", default = 0L)
)
opt = optparse::parse_args(optparse::OptionParser(option_list = option_list))

if(is.null(opt$maindir)) maindir = getwd() else maindir = opt$maindir
setwd(maindir)
cat("current working directory:", getwd(), "\n")

barcode_file = opt$barcode_file
data_folder = opt$data_folder
baselevel_treat = opt$baselevel_treat
number_cores = opt$number_cores
include_nontested_level = opt$include_nontested_level
number_permutations = opt$number_permutations
include_replicate_in_analysis = opt$include_replicate

bar_list = data.table::fread(barcode_file, header = T)
cat("imported bar list (head):\n")
print(head(bar_list))
cat("number of sequences:", nrow(bar_list),"\n")

resu_folder = file.path(data_folder, "results")
if(!dir.exists(resu_folder)) {
  dir.create(resu_folder)
  cat(paste0("\nresult stored in new directory '", resu_folder, "'\n"))
} else {
  cat(paste0("\nresult stored in existing directory '", resu_folder, "'\n"))
}

glm_formula = "cbind(cInd, cNet) ~ treat + donor"
if(include_replicate_in_analysis) glm_formula = paste(glm_formula, "+ replicate")
cat("binomial regression estimated with formula: ", glm_formula, "\n")
if(include_replicate_in_analysis) cat("files are saved with 'adjustedForReplicate' tag\n")
glm_formula = as.formula(glm_formula)

glm_permutation = function(data, formula, start = NULL, include_nontested = TRUE) {
  treat_levels = levels(data$treat)
  ref_level = head(treat_levels, 1)
  alt_level = tail(treat_levels,-1)
  
  treat_coef = vector("double", length(alt_level))
  for(i in seq_along(alt_level)) { #i=1
    di = data.table::copy(data)
    da = di[!(treat %in% c(ref_level, alt_level[i])),]
    di = di[treat %in% c(ref_level, alt_level[i]),]
    di = di[, treat := sample(treat), by = donor] # permute treat between reference and target alternative level
    
    if(include_nontested) {
      di = rbind(di, da)
      si = start
    } else {
      si = start[-which(names(start) %in% paste0("treat", setdiff(alt_level, alt_level[i])))]
    }
    
    mi = glm(formula, data = di, family = binomial(), start = si, method = brglm2::"brglmFit")
    treat_coef[i] = mi$coefficients[paste0("treat", alt_level[i])] # coefficient of target alternative level
  }
  names(treat_coef) = paste0("treat", alt_level)
  return(treat_coef)
}





cat("\n\nstart loop\n\n")
for(qq in 1:nrow(bar_list)) { #qq = 1L
  barcode = bar_list$barcode[qq]
  cat("- processing barcode", barcode, paste0("(", qq, "/", nrow(bar_list), ")"),"\n")
  
  files = list.files(path = file.path(maindir, data_folder), pattern = paste0("^", barcode), full.names = TRUE)
  
  all_data = data.table::data.table()
  for(i in seq_along(files)) all_data = rbind(all_data, data.table::fread(files[i], header=T))
  
  names(all_data)[names(all_data) == "#Reads"] = "reads"
  all_data$edited = as.numeric(!all_data$Unedited)
  
  result_tot = all_data[, .(umiTotReads = sum(reads)), by = .(UMI, replicate, donor, group)]
  result_edt  = all_data[edited==1, .(editingReads = sum(reads)), by = .(UMI, replicate, donor, group)]
  result_all = merge(result_tot, result_edt, by = c("UMI", "replicate", "donor", "group"), all.x = T)
  result_all = result_all[is.na(editingReads), editingReads := 0]
  result_all$uneditingReads = result_all$umiTotReads - result_all$editingReads
  
  levels_treat = sort(unique(result_all$group))
  if(baselevel_treat %in% levels_treat) {
    levels_treat = unique(c(baselevel_treat, levels_treat))
  } else {
    warning("baselevel", baselevel_treat, "not in data, levels data (first is default base level):", paste0(levels_treat, collapse = ", "))
  }
  
  result_all$group = factor(result_all$group, levels = levels_treat)
  names(result_all) = c("UMI","replicate", "donor", "treat" ,"cTot", "cInd", "cNet")
  result_all[, rep := paste0(replicate, donor, treat)]
  result_all$donor = factor(result_all$donor)
  
  mod = glm(glm_formula, data = result_all, family = binomial(), method = brglm2::"brglmFit")
  eff_joint = ggeffects::ggpredict(mod, c("donor", "treat"))
  eff_marg  = ggeffects::ggeffect (mod, c("treat"))
  
  estimates = mod$coefficients
  treat_est = estimates[grepl("^treat", names(estimates))]
  test_stat = coef(summary(mod))[paste0("treat", setdiff(levels_treat, baselevel_treat)), 3]
  asym_pvalue = 1 - pnorm(test_stat)
  
  if(number_permutations > 0L) {
    # start_parameter = estimates # non stable, as permutation estimates are very different
    start_parameter = rep(0, length(estimates)) # better
    
    cat("    start permutation tests\n")
    permuted_treat_est = parallel::mclapply(
      1:number_permutations, 
      function(x) glm_permutation(result_all, glm_formula, start_parameter, include_nontested_level), 
      mc.cores = number_cores
    )
    permuted_treat_est = do.call("rbind", permuted_treat_est)
    cat("    end permutation tests\n")
    
    unstable_permutations = which(apply(permuted_treat_est, 1, function(x) any(abs(x) > 1e4)))
    if(length(unstable_permutations)) {
      if(length(unstable_permutations) / nrow(permuted_treat_est) > .01) warning(length(unstable_permutations), " permutations removed\n")
      permuted_treat_est = permuted_treat_est[-unstable_permutations,]
    } 
    
    perm_pvalue = colMeans(permuted_treat_est > do.call("rbind", replicate(nrow(permuted_treat_est), treat_est, simplify = F)))
  } else {
    perm_pvalue = rep(NA_real_, length(asym_pvalue))
  } 
  
  idej = paste0(eff_joint$x, "-", eff_joint$group)
  idem = as.character(eff_marg$x)
  
  output = data.table::data.table(t(c(
    estimates, test_stat, 
    asym_pvalue, perm_pvalue, 
    eff_joint$predicted, eff_joint$std.error, eff_joint$conf.high, eff_joint$conf.low,
    eff_marg$predicted,  eff_marg$std.error,  eff_marg$conf.high,  eff_marg$conf.low
  )))
  names(output) = c(
    paste0("estim-", names(estimates)), paste0("zScore-", names(test_stat)), 
    paste0("asymPval-", names(asym_pvalue)), paste0("permPval-", names(perm_pvalue)),
    paste0("predicted-", idej), paste0("stdError-", idej), paste0("confHigh-", idej), paste0("confLow-", idej),
    paste0("predicted-", idem), paste0("stdError-", idem), paste0("confHigh-", idem), paste0("confLow-", idem)
  )
  output = cbind(barcode = barcode, output)
  outnam = file.path(resu_folder, paste0(barcode, "_output_", c("", "adjustedForReplicate")[1+include_replicate_in_analysis],".csv")) 
  
  data.table::fwrite(x = output, file = outnam, sep = "\t", col.names = T, row.names = F, quote = F)
  cat("    stored file:    ", outnam, "\n")
}
cat("end loop\n")
