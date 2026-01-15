option_list = list(
  optparse::make_option("--maindir", type = "character", default = NULL, help = "Main working directory (default: current directory)"),
  optparse::make_option("--data_folder", type = "character", default = "target_umi_barcode_table_full"),
  optparse::make_option("--experiments", type = "character", default = "dedudS_filtered_v3,dedudS_filtered_v3_dedup,raw,raw_dedup"),
  optparse::make_option("--barcode_file", type = "character", default = "barcodeList.csv"),
  optparse::make_option("--baselevel_treat", type = "character", default = "NoEP"),
  optparse::make_option("--format_plot", type = "character", default = "svg"),
  optparse::make_option("--height_plot", type = "double", default = 5.0),
  optparse::make_option("--width_plot", type = "double", default = 5.0),
  optparse::make_option("--y_min", type = "double", default = 1e-3),
  optparse::make_option("--shift_predicted_reference", type = "logical", default = T),
  optparse::make_option("--shift_x_treat", type = "double", default = 0.2),
  optparse::make_option("--colors_plot", type = "character", default = "grey20,green4,deeppink3"),
  optparse::make_option("--include_replicate", type = "logical", default = FALSE)
)
opt = optparse::parse_args(optparse::OptionParser(option_list = option_list))

if(is.null(opt$maindir)) maindir = getwd() else maindir = opt$maindir
setwd(maindir)
cat("current working directory:", getwd(), "\n")

barcode_file = opt$barcode_file
data_folder = opt$data_folder
experiments = strsplit(opt$experiments, ",")[[1]] # experiments = c("dedudS_filtered_v3", "dedudS_filtered_v3_dedup", "raw", "raw_dedup")
baselevel_treat = opt$baselevel_treat
include_replicate_in_analysis = opt$include_replicate

format_plot = opt$format_plot
height_plot = opt$height_plot
width_plot  = opt$width_plot
shift_predicted_reference = opt$shift_predicted_reference
shift_x_treat = opt$shift_x_treat
colors_plot = strsplit(opt$colors_plot, ",")[[1]] 
y_min = opt$y_min

plot_folder = file.path(maindir, "plots")
individual_plot_folders = file.path(plot_folder, paste0("individual_", c("joint", "marginal")))
if(!dir.exists(plot_folder)) {
  dir.create(plot_folder)
  cat(paste0("\nplots stored in new directory '", plot_folder, "'\n"))
} else {
  cat(paste0("\nplots stored in existing directory '", plot_folder, "'\n"))
}
for(i in seq_along(individual_plot_folders)) {
  if(!dir.exists(individual_plot_folders[i])) {
    dir.create(individual_plot_folders[i])
    cat(paste0(c("joint", "marginal")[i], " individual plots stored in new directory '", individual_plot_folders[i], "'\n"))
  } else {
    cat(paste0(c("joint", "marginal")[i], " individual plots stored in existing directory '", individual_plot_folders[i], "'\n"))
  }
}

files = purrr::map(experiments, ~ list.files(
  file.path(data_folder, .x, "results"), 
  pattern = paste0("output_", c(".csv", "adjustedForReplicate.csv")[1 + include_replicate_in_analysis]), 
  full.names = F
))
barcodes = purrr::map(files, stringr::str_split, pattern = "_output_")
barcodes = purrr::map(barcodes, ~ purrr::map_chr(.x, ~ .x[[1]][1]))
nbarcodes = length(barcodes[[1]])

for(i in 1:length(barcodes)) stopifnot(length(barcodes[[i]]) == length(unique(barcodes[[i]])))

analysis = purrr::map2(experiments, files, ~ purrr::map(file.path(data_folder, .x, "results", .y), data.table::fread))
analysis = purrr::map(analysis, ~ do.call("rbind", .x))

names(barcodes) = names(analysis) = experiments

analysis = purrr::map(analysis, ~ data.table::setkey(.x, barcode))

name_cols_analysis = names(analysis[[1]])
treat_levels = c(baselevel_treat, sort(unique(unlist(purrr::map(strsplit(name_cols_analysis, "treat"), ~ .x[2])))))

if(length(colors_plot) != length(treat_levels)) {
  warning("length 'colors_plot' ", paste0("(", length(colors_plot), ")"), 
          " plot is not equal to the number of treatments ",  paste0("(", length(treat_levels), ")"),
          "\nnew colors generated with 'rainbow(length(treat_levels))'")
  colors_plot = rainbow(length(treat_levels))
}

x0coord = shift_x_treat + (1 - 2*shift_x_treat) / (length(treat_levels)-1) * (0:(length(treat_levels)-1))
names(x0coord) = treat_levels

log10_minor_break = function (...){
  function(x) {
    minx         = floor(min(log10(x), na.rm=T))-1;
    maxx         = ceiling(max(log10(x), na.rm=T))+1;
    n_major      = maxx-minx+1;
    major_breaks = seq(minx, maxx, by=1)
    minor_breaks = 
      rep(log10(seq(1, 9, by=1)), times = n_major)+
      rep(major_breaks, each = 9)
    return(10^(minor_breaks))
  }
}




for(ex in seq_along(experiments)) { # ex=1
  cat("\nprocessing data experiment", experiments[ex], "\n")
  
  anex = analysis[[ex]]
  barcodes = anex$barcode
  
  for(bc in seq_along(barcodes)) { # bc=1
    cat(" - processing barcode", barcodes[bc], "\n")
    files = list.files(path = file.path(maindir, data_folder, experiments[ex]), pattern = paste0("^", barcodes[bc]), full.names = TRUE)
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
    
    data_plot = result_all[, .(sum_cInd = sum(cInd), sum_cTot = sum(cTot)), by = .(replicate, donor, treat)]
    data_plot$prop = data_plot$sum_cInd / data_plot$sum_cTot
    data_plot$x = as.numeric(data_plot$donor) + shift_x_treat * c(-1:1)[as.integer(data_plot$treat)]
    
    select_columns_analysis = c(
      grep("^predicted-", names(anex), value = T), 
      grep("^confLow-",   names(anex), value = T), 
      grep("^confHigh-",  names(anex), value = T)
    )  
    anbc = anex[barcode == barcodes[bc], ..select_columns_analysis]
    split_names = strsplit(names(anbc), "-")
    
    whjo = which(sapply(split_names, function(x) length(x) == 3))
    anjo = anbc[, .SD, .SDcols = whjo]
    wpjo = grep("predicted-", names(anjo), value = T)
    wljo = grep("confLow-",   names(anjo), value = T)
    whjo = grep("confHigh-",  names(anjo), value = T)
    anjo = data.table::data.table(
      treat = purrr::map_chr(strsplit(wpjo, "-"), ~ .x[3]), 
      donor = purrr::map_chr(strsplit(wpjo, "-"), ~ .x[2]), 
      predicted = as.double(anjo[,..wpjo]),
      confLow   = as.double(anjo[,..wljo]),
      confHigh  = as.double(anjo[,..whjo])
    )
    anjo$treat = factor(anjo$treat, levels_treat)
    anjo$donor = factor(anjo$donor)
    
    whma = which(sapply(split_names, function(x) length(x) == 2))
    anma = anbc[, .SD, .SDcols = whma]
    wpma = grep("predicted-", names(anma), value = T)
    wlma = grep("confLow-",   names(anma), value = T)
    whma = grep("confHigh-",  names(anma), value = T)
    anma = data.table::data.table(
      treat = purrr::map_chr(strsplit(wpma, "-"), ~ .x[2]), 
      predicted = as.double(anma[,..wpma]),
      confLow   = as.double(anma[,..wlma]),
      confHigh  = as.double(anma[,..whma])
    )
    anma$treat = factor(anma$treat, levels_treat)
    
    anjo$x = as.numeric(anjo$donor) + shift_x_treat * c(-1:1)[as.integer(anjo$treat)]
    
    ggplot2::ggplot() +
      ggplot2::geom_hline(ggplot2::aes(yintercept = predicted * 100, color = treat), data = anma, alpha = .6) + 
      ggplot2::geom_hline(ggplot2::aes(yintercept = predicted * 100 + .01, color = treat), data = anma[treat == baselevel_treat], alpha = .6, linetype = "dotdash") + 
      ggplot2::geom_hline(ggplot2::aes(yintercept = predicted * 100 + 0.1, color = treat), data = anma[treat == baselevel_treat], alpha = .6, linetype = "dashed") + 
      ggplot2::geom_point(ggplot2::aes(x = x, y = predicted * 100, color = treat), data = anjo, alpha = .6) +
      ggplot2::geom_point(ggplot2::aes(x = x, y = prop * 100, color = treat), data = data_plot, pch = 3) +
      ggplot2::scale_color_discrete(type = colors_plot) +
      ggplot2::scale_y_continuous(limits = c(y_min, 100), trans = "log10", breaks = 10^(round(log10(y_min)):2), minor_breaks = log10_minor_break()) +
      ggplot2::scale_x_continuous(breaks = seq_len(length(levels(anjo$donor))), minor_breaks = NULL, labels = levels(anjo$donor)) +
      ggplot2::theme_minimal() + ggplot2::ylab("predicted (percentage, log scale)") + ggplot2::xlab("donor") +
      ggplot2::theme(
        #legend.position = "none", 
        plot.background=ggplot2::element_rect(fill="white")
      ) -> plt
    
    ggplot2::ggsave(
      file.path(individual_plot_folders[1], paste0("dot-", experiments[ex], "-", barcodes[bc], "-.", format_plot)), 
      plt, height = height_plot, width = width_plot
    )
    
    ggplot2::ggplot() +
      ggplot2::geom_hline(ggplot2::aes(yintercept = predicted * 100, color = treat), data = anma) +
      ggplot2::geom_segment(ggplot2::aes(x = treat, y = pmax(confLow * 100, y_min), yend = pmax(pmin(confHigh * 100, 100), y_min), color = treat), data = anma) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = predicted * 100 + 0.0, color = treat), data = anma[treat == baselevel_treat], alpha = .6, linetype = "solid") + 
      ggplot2::geom_hline(ggplot2::aes(yintercept = predicted * 100 + .01, color = treat), data = anma[treat == baselevel_treat], alpha = .6, linetype = "dotdash") + 
      ggplot2::geom_hline(ggplot2::aes(yintercept = predicted * 100 + 0.1, color = treat), data = anma[treat == baselevel_treat], alpha = .6, linetype = "dashed") + 
      ggplot2::geom_point(ggplot2::aes(x = treat, y = predicted * 100, color = treat, shape = donor), data = anjo) +
      ggplot2::scale_color_discrete(type = colors_plot) +
      ggplot2::scale_shape_manual(values = c(2,5,6)) + 
      ggplot2::scale_y_continuous(limits = c(y_min, 100), trans = "log10", breaks = 10^(round(log10(y_min)):2), minor_breaks = log10_minor_break()) +
      #scale_y_continuous(trans = "log10", minor_breaks = log10_minor_break()) +
      #scale_x_continuous(breaks = 1:3, minor_breaks = NULL) +
      ggplot2::theme_minimal() + ggplot2::ylab("predicted (percentage, log scale)") + ggplot2::xlab("donor") +
      ggplot2::theme(
        #legend.position = "none", 
        plot.background=ggplot2::element_rect(fill="white")
      ) -> plt
    
    ggplot2::ggsave(
      file.path(individual_plot_folders[2], paste0("dot-", experiments[ex], "-", barcodes[bc], "-.", format_plot)), 
      plt, height = height_plot, width = width_plot
    )
  }
}  