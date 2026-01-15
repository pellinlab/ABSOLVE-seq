option_list = list(
  optparse::make_option("--maindir", type = "character", default = NULL, help = "Main working directory (default: current directory)"),
  optparse::make_option("--data_folder", type = "character", default = "target_umi_barcode_table_full"),
  optparse::make_option("--experiments", type = "character", default = "dedudS_filtered_v3,dedudS_filtered_v3_dedup,raw,raw_dedup"),
  optparse::make_option("--baselevel_treat", type = "character", default = "NoEP"),
  optparse::make_option("--power_sizes", type = "character", default = "0.01,0.1,1"),
  optparse::make_option("--significance_level", type = "double", default = 0.05),
  optparse::make_option("--format_plot", type = "character", default = "svg"),
  optparse::make_option("--height_plot", type = "double", default = 8.0),
  optparse::make_option("--width_plot", type = "double", default = 17.0),
  optparse::make_option("--filled_color_heat", type = "character", default = "deeppink3,green4"),
  optparse::make_option("--size_barcodes", type = "double", default = 3.0),
  optparse::make_option("--heigh_top_proportion", type = "double", default = 0.8),
  optparse::make_option("--significance_powers", type = "character", default = "0.0001,0.001"),
  optparse::make_option("--shift_power", type = "double", default = 0.6),
  optparse::make_option("--include_replicate", type = "logical", default = FALSE)
)
opt = optparse::parse_args(optparse::OptionParser(option_list = option_list))

if(is.null(opt$maindir)) maindir = getwd() else maindir = opt$maindir
setwd(maindir)
cat("current working directory:", getwd(), "\n")

data_folder = opt$data_folder
experiments = strsplit(opt$experiments, ",")[[1]] # experiments = c("dedudS_filtered_v3", "dedudS_filtered_v3_dedup", "raw", "raw_dedup")
baselevel_treat = opt$baselevel_treat
include_replicate_in_analysis = opt$include_replicate
power_sizes = as.double(strsplit(opt$power_sizes, ",")[[1]])
significance_level = opt$significance_level
format_plot = opt$format_plot
height_plot = opt$height_plot
width_plot  = opt$width_plot
size_barcodes = opt$size_barcodes
filled_color_heat = strsplit(opt$filled_color_heat, ",")[[1]]
heigh_top_proportion = opt$heigh_top_proportion
significance_powers = as.double(strsplit(opt$significance_powers, ",")[[1]])
shift_power = opt$shift_power

#shift_folders = (0:3) * 7 + 1 # +1 for power 0 / not used 

plot_folder = file.path(maindir, "plots")
if(!dir.exists(plot_folder)) {
  dir.create(plot_folder)
  cat(paste0("\nplots stored in new directory '", plot_folder, "'\n"))
} else {
  cat(paste0("\nplots stored in existing directory '", plot_folder, "'\n"))
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

powers = purrr::map(experiments, function(e) purrr::map(
  power_sizes, ~ data.table::fread(paste0(file.path(data_folder, e, "results", "powerAnalysis"), "_", .x, "_.csv"))
))

names(barcodes) = names(analysis) = names(powers) = experiments

for(j in seq_along(experiments)) for(i in seq_along(power_sizes)) { #i=j=1
  names(powers[[j]][[i]])[-1] = paste0("power-", power_sizes[i], "-", names(powers[[j]][[i]])[-1])
}

analysis = purrr::map(analysis, ~ data.table::setkey(.x, barcode))
powers = purrr::map(powers, ~ purrr::map(.x, ~ data.table::setkey(.x, barcode)))

for(i in seq_along(powers)) if(i > 1) stopifnot(all(powers[[i]]$barcode == powers[[1]]$barcode)) # possible merge with cbind
powers = purrr::map(powers, ~ do.call("cbind", .x))
for(i in seq_along(powers)) {
  stopifnot(all(powers[[i]][[1]] == analysis[[i]]$barcode)) # possible merge with cbind
  analysis[[i]] = cbind(analysis[[i]], powers[[i]][, .SD, .SDcols = grep("^power", colnames(powers[[i]]), value = TRUE)])
} 

rm(powers, barcodes)

name_cols_analysis = names(analysis[[1]])
treat_levels = c(baselevel_treat, sort(unique(unlist(purrr::map(strsplit(name_cols_analysis, "treat"), ~ .x[2])))))
names(filled_color_heat) = treat_levels[-1]

if(length(treat_levels) != 3L) stop("heat plot requires 3 treatiment levels")

cat("\n compute heatplot: \n")

select_pvalues = c("barcode", grep("^asymPval", name_cols_analysis, value = TRUE))
pvalues = purrr::map(analysis, ~ .x[,..select_pvalues])
pvalues = purrr::map2(pvalues, names(pvalues), ~ .x[, experiment := .y])
if(length(treat_levels) > 1L) for(i in 2:length(treat_levels)) {
  pvalues = purrr::map(pvalues, ~ .x[, VNEW := p.adjust(.x[[i]], method = "fdr")])
  pvalues = purrr::map(pvalues, ~ data.table::setnames(.x, "VNEW", paste0("adjusted-treat", treat_levels[i])))
}
pvalues = do.call("rbind", pvalues)
if(length(treat_levels) > 1L) for(i in 2:length(treat_levels)) { #i=2
  pvalues[, VNEW := pvalues[[paste0("adjusted-treat", treat_levels[i])]] < significance_level]
  data.table::setnames(pvalues, "VNEW", paste0("significant-treat", treat_levels[i]))
}

barcodes = sort(unique(pvalues$barcode))
barcodes_ids = 1:length(barcodes)
names(barcodes_ids) = barcodes
exper_ids = 1:length(experiments)
names(exper_ids) = experiments

pvalues$barcode_id = barcodes_ids[pvalues$barcode]
pvalues$exper_id = exper_ids[pvalues$experiment]

pvalues_long = data.table::melt(
  pvalues,
  id.vars = c("barcode", "experiment", "barcode_id", "exper_id"),
  measure = patterns(
    asymPval   = "^asymPval-treat",
    adjusted   = "^adjusted-treat",
    significant = "^significant-treat"
  ),
  variable.name = "treat"
)
pvalues_long[, treat_id := as.integer(treat)]
pvalues_long[, treat := factor(treat_levels[1+treat_id], levels = treat_levels)]

pvalues_long$x1 = as.double(pvalues_long$barcode_id)
pvalues_long$x2 = pvalues_long$x1
pvalues_long$x3 = pvalues_long$x1 + 1.0
pvalues_long$y2 = pvalues_long$y1 = 0; pvalues_long$y3 = 1

treat_upper = treat_levels[2]
treat_lower = treat_levels[3]
pvalues_long$y2[pvalues_long$treat == treat_upper] = pvalues_long$y2[pvalues_long$treat == treat_upper] + 1
pvalues_long$x2[pvalues_long$treat == treat_lower] = pvalues_long$x2[pvalues_long$treat == treat_lower] + 1

pvalues_long = data.table::melt(
  pvalues_long,
  id.vars = setdiff(names(pvalues_long), c("x1","x2","x3","y1","y2","y3")),
  measure = patterns(
    x = "^x",
    y = "^y"
  ),
  variable.name = "coord_id"
)
pvalues_long$coord_id = NULL
data.table::setorder(pvalues_long, exper_id, barcode_id, treat_id)
pvalues_long[, triangles_id := rep(1:(nrow(pvalues_long)/3), each = 3)]

select_powers = c(
  "barcode", 
  tail(grep("^predicted-", name_cols_analysis, value = T), length(treat_levels)), 
  grep("^power", name_cols_analysis, value = TRUE)
)
powers = purrr::map(analysis, ~ .x[,..select_powers])
powers = purrr::map2(powers, names(powers), ~ .x[, experiment := .y])
powers = do.call("rbind", powers)
if(length(treat_levels) > 1L) for(i in 2:length(treat_levels)) { #i=2
  powers[, VNEW := powers[[paste0("predicted-", treat_levels[i])]] - powers[[paste0("predicted-", treat_levels[1])]]]
  data.table::setnames(powers, "VNEW", paste0("delta-", treat_levels[i]))
}
if(length(treat_levels) > 1L) for(i in 2:length(treat_levels)) for(j in seq_along(significance_powers)) { #i=2
  powers[, VNEW := (powers[[paste0("delta-", treat_levels[i])]] > significance_powers[j])]
  data.table::setnames(powers, "VNEW", paste0("significant-", format(significance_powers[j], scientific = F), "-", treat_levels[i]))
}  
powers$barcode_id = barcodes_ids[powers$barcode]
powers$exper_id = exper_ids[powers$experiment]

select_powers_heat = c(
  'barcode', 'barcode_id', 'experiment','exper_id',
  'power-0.01-treatHiFi', 'power-0.01-treatWTCas9', 'power-0.1-treatHiFi', 'power-0.1-treatWTCas9', 
  'significant-0.001-HiFi', 'significant-0.001-WTCas9', 'significant-0.0001-HiFi', 'significant-0.0001-WTCas9'
)
powers_heat = powers[,..select_powers_heat]
powers_heat_long = data.table::melt(
  powers_heat,
  id.vars = c("barcode", "experiment", "barcode_id", "exper_id"),
  measure = patterns(
    pow0001 = "^power-0.01-treat",
    pow001  = "^power-0.1-treat",
    sig0001 = "^significant-0.0001-",
    sig001  = "^significant-0.001-"
  ),
  variable.name = "treat"
)
powers_heat_long[, treat_id := as.integer(treat)]
powers_heat_long[, treat := factor(treat_levels[1+treat_id], levels = treat_levels)]
powers_heat_long = data.table::melt(
  powers_heat_long,
  id.vars = c("barcode", "experiment", "treat", "barcode_id", "exper_id", "treat_id"),
  measure = patterns(
    power  = "^pow00",
    signif = "^sig00"
  ),
  variable.name = "size"
)
powers_heat_long[, size := factor(c("s0001", "s001")[size], levels = c("s0001", "s001"))]

powers_heat_long$x1 = as.double(powers_heat_long$barcode_id)
powers_heat_long$x2 = powers_heat_long$x1
powers_heat_long$x3 = powers_heat_long$x1 + 1.0
powers_heat_long$y2 = powers_heat_long$y1 = 0; powers_heat_long$y3 = 1

treat_upper = treat_levels[2]
treat_lower = treat_levels[3]
powers_heat_long$y2[powers_heat_long$treat == treat_upper] = powers_heat_long$y2[powers_heat_long$treat == treat_upper] + 1
powers_heat_long$x2[powers_heat_long$treat == treat_lower] = powers_heat_long$x2[powers_heat_long$treat == treat_lower] + 1

powers_heat_long = data.table::melt(
  powers_heat_long,
  id.vars = setdiff(names(powers_heat_long), c("x1","x2","x3","y1","y2","y3")),
  measure = patterns(
    x = "^x",
    y = "^y"
  ),
  variable.name = "coord_id"
)
powers_heat_long$coord_id = NULL
data.table::setorder(powers_heat_long, exper_id, barcode_id, treat_id, size)
powers_heat_long[1:13]
powers_heat_long[, triangles_id := rep(1:(nrow(powers_heat_long)/3), each = 3)]

# shift y
#pvalues_long$y = pvalues_long$y + shift_folders[pvalues_long$exper_id]
#powers_heat_long$y = powers_heat_long$y + shift_folders[powers_heat_long$exper_id]
powers_heat_long$y = powers_heat_long$y + c(s0001 = 1L, s001 = 2L)[powers_heat_long$size] + shift_power

background_shade = data.table::data.table()
background_shade$xmin = 1:nbarcodes
background_shade$xmax = 1:nbarcodes + 1L
background_shade$ymin = min(pvalues_long$y) - shift_power
background_shade$ymax = max(powers_heat_long$y) + 2 + shift_power + shift_power
background_shade$brid = 1:nbarcodes
background_shade$col  = rep(paste0("c", 1:(length(treat_levels)-1)), times = ceiling(nbarcodes/(length(treat_levels)-1)))[1:nbarcodes]

ggplot2::ggplot() +
  ggplot2::geom_rect(ggplot2::aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax, fill = col), data = background_shade) +
  ggplot2::scale_fill_manual(values = c("c1" = "grey95", "c2" = "white")) +
  ggnewscale::new_scale_fill() +
  ggplot2::geom_polygon(ggplot2::aes(x = x, y = y, group = triangles_id, alpha = as.numeric(significant), fill = treat), data = pvalues_long, color = "black", linewidth = .2) +
  ggplot2::facet_wrap(ggplot2::vars(experiment), ncol = 1, strip.position = "left") + 
  ggplot2::geom_polygon(ggplot2::aes(x = x, y = y, group = triangles_id, alpha = as.numeric(signif), fill = treat), data = powers_heat_long, color = "black", linewidth = .2) +
  ggplot2::facet_wrap(ggplot2::vars(experiment), ncol = 1, strip.position = "left") + 
  ggplot2::scale_fill_manual(values = filled_color_heat) +
  ggplot2::scale_alpha_continuous(range = c(0,1)) +
  ggnewscale::new_scale_fill() +
  ggplot2::geom_polygon(ggplot2::aes(x = x, y = y + 2 + shift_power, fill = power, group = triangles_id), data = powers_heat_long, color = "black", linewidth = .2) +
  ggplot2::facet_wrap(ggplot2::vars(experiment), ncol = 1, strip.position = "left") +
  ggplot2::scale_fill_gradient2(low = "blue3", mid = "grey70", high = "red3", midpoint = 0.5) +
  ggplot2::scale_x_continuous(limits = c(1, nbarcodes+1), breaks = NULL, minor_breaks = NULL, expand = c(0,0)) +
  #ggplot2::scale_y_continuous(breaks = NULL, minor_breaks = NULL) +
  ggplot2::scale_y_continuous(
    breaks = c(shift_power, 1+2*shift_power, 2+2*shift_power, 4+shift_power, 5+shift_power), 
    labels = c(1:5),
    minor_breaks = NULL) +
  ggplot2::ylab(NULL) + ggplot2::xlab(NULL) +
  #scale_fill_gradient(high = "red3", low = "grey50") +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    strip.placement = "outside"
    #legend.position = "none"
  ) -> heatplt
#heatplt

background_shade_text = background_shade
background_shade_text$ymin = -1
background_shade_text$ymax = 1
ggplot2::ggplot() +
  ggplot2::geom_rect(ggplot2::aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax, fill = col), data = background_shade_text) +
  ggplot2::scale_fill_manual(values = c("c1" = "grey95", "c2" = "white")) +
  ggplot2::geom_text(ggplot2::aes(x = barcode_id + .5, y = 0, label = barcode), 
                     data = unique(pvalues[,.(barcode, barcode_id)]), angle = 90, fontface = "bold", size = size_barcodes) +
  ggplot2::scale_x_continuous(limits = c(1,nbarcodes+1), breaks = NULL, minor_breaks = NULL, expand = c(0,0)) +
  ggplot2::scale_y_continuous(breaks = 0, minor_breaks = NULL, labels = "barcode") +
  ggplot2::xlab(NULL) + ggplot2::ylab(NULL) +
  ggplot2::theme_minimal() + 
  ggplot2::theme(
    legend.position = "none",
    axis.text.y = ggplot2::element_text(angle = 90, hjust = 0.5),
    axis.ticks.y = ggplot2::element_blank()
  ) -> textplt
#textplt

plt = heatplt + textplt + patchwork :: plot_layout(ncol = 1, height = c(heigh_top_proportion, 1-heigh_top_proportion), axes = "collect_x")
#plt
ggplot2::ggsave(file.path(plot_folder, paste0("heat.", format_plot)), plt, height = height_plot, width = width_plot)
cat("stored plot:   ", file.path(plot_folder, paste0("heat.", format_plot)), "\n")


