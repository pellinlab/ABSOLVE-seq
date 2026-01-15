option_list = list(
  optparse::make_option("--maindir", type = "character", default = NULL, help = "Main working directory (default: current directory)"),
  optparse::make_option("--data_folder", type = "character", default = "target_umi_barcode_table_full/dedudS_filtered_v3"),
  optparse::make_option("--baselevel_treat", type = "character", default = "NoEP"),
  optparse::make_option("--number_cores", type = "integer", default = 1L),
  optparse::make_option("--include_replicate", type = "logical", default = FALSE),
  optparse::make_option("--number_simulations", type = "integer", default = 1000L),
  optparse::make_option("--significance_level", type = "double", default = 0.05),
  optparse::make_option("--power_sizes", type = "character", default = "0.01,0.1,1")
)
opt = optparse::parse_args(optparse::OptionParser(option_list = option_list))

if(is.null(opt$maindir)) maindir = getwd() else maindir = opt$maindir
setwd(maindir)
cat("current working directory:", getwd(), "\n")

data_folder = opt$data_folder
baselevel_treat = opt$baselevel_treat
number_cores = opt$number_cores
number_simulations = opt$number_simulations
include_replicate_in_analysis = opt$include_replicate
power_sizes = as.double(strsplit(opt$power_sizes, ",")[[1]])
significance_level = opt$significance_level

files = list.files(
  file.path(data_folder, "results"), 
  pattern = paste0("output_", c(".csv", "adjustedForReplicate.csv")[1+include_replicate_in_analysis]), 
  full.names = F
)
barcodes = stringr::str_split(files, "_output_")
barcodes = purrr::map_chr(barcodes, ~ .x[[1]][1])

stopifnot(length(barcodes) == length(unique(barcodes)))
if(include_replicate_in_analysis) warning("not yet implemented: 'replicate' is not uses ad covariate model")

analysis = purrr::map(file.path(data_folder, "results", files), data.table::fread)
analysis = do.call("rbind", analysis)

analysis$adjustedForReplicate = include_replicate_in_analysis
cat("imported data with", nrow(analysis), "barcodes\n")
cat("output will be store in files:  ", file.path(data_folder, "results", paste0("powerAnalysis_", "*", "_.csv")), "\n")

simulate_model = function(tab, alpha = significance_level) {
  treat_levels = levels(tab$treat)
  ref_level = head(treat_levels, 1)
  alt_level = tail(treat_levels,-1)
  
  # simulate data from null
  tab$cInd = purrr::map2_int(tab$cTot, tab$powerProb, ~ sum(rbinom(n = .x, size = 1, prob = .y)))
  
  tab$cNet = tab$cTot - tab$cInd
  mod = glm(cbind(cInd,cNet) ~ treat + donor, data = tab, family = binomial, method = brglm2::"brglmFit")
  zval = coef(summary(mod))[paste0("treat", alt_level), "z value"]
  pval = 1 - pnorm(zval)
  
  return(pval < alpha)
}





cat("\n\nstart power analysis")
for(ps in seq_along(power_sizes)){ # ps=1
  cat("\n\npower size", paste0(power_sizes[ps], "%"), "\n\n")
  
  effect_size = power_sizes[ps] / 100
  povtab = data.table::data.table()
  
  for(bc in seq_along(barcodes)) { #bc=1
    cat("- processing barcode", bc, "/", length(barcodes), "\n")
    analysis_bc = analysis[barcode == barcodes[bc]]
    files = list.files(path = file.path(maindir, data_folder), pattern = paste0("^", barcodes[bc]), full.names = TRUE)
    
    tab = data.table::data.table(
      donor = rep(x = "NA", length(files)),
      treat = rep(x = "NA", length(files)),
      cInd  = rep(x = 0L, length(files)),
      cTot  = rep(x = 0L, length(files))
    )
    
    for(fi in seq_along(files)) { #fi=1
      data = data.table::fread(files[fi], header=T)
      tab[fi,1] = unique(data$donor)
      tab[fi,2] = unique(data$group)
      tab[fi,3] = sum(data$`#Reads`[data$Unedited == "FALSE"])
      tab[fi,4] = sum(data$`#Reads`)
    }
    
    levels_treat = sort(unique(tab$treat))
    if(baselevel_treat %in% levels_treat) {
      levels_treat = unique(c(baselevel_treat, levels_treat))
    } else {
      warning("baselevel", baselevel_treat, "not in data, levels data (first is default base level):", paste0(levels_treat, collapse = ", "))
    }
    tab$treat = factor(tab$treat, levels = levels_treat)
    tab$donor = factor(tab$donor)
    tab$cNet  = tab$cTot - tab$cInd
    
    ndonor=length(unique(tab$donor))
    donors=unique(tab$donor)
    
    whpred = purrr::map2_int(tab$donor, tab$treat, ~ which(paste("predicted", .x, .y, sep = "-") == names(analysis_bc)))
    tab$predicted = purrr::map_dbl(whpred, ~ analysis_bc[[.x]])
    tab$powerProb = 0.0
    for(dn in 1:ndonor) { #dn=1
      tab$powerProb[tab$donor == donors[dn]] = tab$predicted[tab$donor==donors[dn] & tab$treat==baselevel_treat][1]
    }
    for(tr in setdiff(levels_treat, baselevel_treat)) {
      tab$powerProb[tab$treat == tr] = tab$powerProb[tab$treat == tr] + effect_size
    }
    
    powerV = parallel::mclapply(1:number_simulations, function(x) simulate_model(tab), mc.cores = number_cores)
    powerV = colMeans(do.call("rbind", powerV))
    
    povtab = rbind(povtab, cbind(data.table::data.table(barcode = barcodes[bc]), t(powerV)))
  }
  
  outnam = file.path(data_folder, "results", paste0("powerAnalysis_", power_sizes[ps], "_.csv"))
  data.table::fwrite(povtab, file = outnam, sep = "\t", row.names = F)
  cat(" stored file:    ", outnam, "\n")
}
cat("end loop\n")
