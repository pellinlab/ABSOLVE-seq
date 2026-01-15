option_list = list(
  optparse::make_option("--maindir", type = "character", default = NULL, help = "Main working directory (default: current directory)"),
  optparse::make_option("--data_folder", type = "character", default = "target_umi_barcode_table_full"),
  optparse::make_option("--experiments", type = "character", default = "dedudS_filtered_v3,dedudS_filtered_v3_dedup,raw,raw_dedup"),
  optparse::make_option("--baselevel_treat", type = "character", default = "NoEP"),
  optparse::make_option("--format_plot", type = "character", default = "svg"),
  optparse::make_option("--height_plot", type = "double", default = 8.0),
  optparse::make_option("--width_plot", type = "double", default = 17.0),
  optparse::make_option("--shift_predicted_reference", type = "logical", default = TRUE),
  optparse::make_option("--shift_x_treat", type = "double", default = 0.2),
  optparse::make_option("--colors_plot", type = "character", default = "grey20,green4,deeppink3"),
  optparse::make_option("--breaks_plot", type = "character", default = "0.04,0.4"),
  optparse::make_option("--breaks_y_by", type = "character", default = "0.02,0.1,20"),
  optparse::make_option("--power_hlines", type = "character", default = "0,0.01,0.1"),
  optparse::make_option("--size_barcodes", type = "double", default = 3.0),
  optparse::make_option("--prop_barcodes", type = "double", default = 0.2),
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

format_plot = opt$format_plot
height_plot = opt$height_plot
width_plot  = opt$width_plot
shift_predicted_reference = opt$shift_predicted_reference
shift_x_treat = opt$shift_x_treat
colors_plot = strsplit(opt$colors_plot, ",")[[1]] 
breaks_plot = as.double(strsplit(opt$breaks_plot, ",")[[1]])
breaks_y_by = as.double(strsplit(opt$breaks_y_by, ",")[[1]])
power_hlines = as.double(strsplit(opt$power_hlines, ",")[[1]])
size_barcodes = opt$size_barcodes
prop_barcodes = opt$prop_barcodes

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

names(barcodes) = names(analysis)

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

select_columns_analysis = c("barcode", paste0("predicted-", treat_levels), paste0("confLow-", treat_levels), paste0("confHigh-", treat_levels))

background_shade = data.table::data.table()
background_shade$xmin = 1:nbarcodes
background_shade$xmax = 1:nbarcodes + 1L
background_shade$brid = 1:nbarcodes
background_shade$col  = rep(paste0("c", 1:(length(treat_levels)-1)), times = ceiling(nbarcodes/(length(treat_levels)-1)))[1:nbarcodes]

for(ex in seq_along(experiments)) { # ex=1
  cat("\n processing data experiment", experiments[ex], "\n")
  
  anex = analysis[[ex]][,..select_columns_analysis]
  if(shift_predicted_reference) {
    stopifnot(names(anex)[1] == "barcode" & names(anex)[2] == paste0("predicted-", baselevel_treat))
    anex = cbind(anex[,1], anex[,-1] - anex[[2]]) # shift predicted reference to 0
  }
  anex$barcode_id = 1:nrow(anex)
  
  anex_full = data.table::melt(
    anex,
    id.vars = c("barcode", "barcode_id"),
    measure = patterns(
      predicted = "^predicted-",
      confLow   = "^confLow-",
      confHigh  = "^confHigh-"
    ),
    variable.name = "treat"
  )
  anex_full = anex_full[,.(
    barcode, barcode_id, 
    treat = factor(treat_levels[as.integer(treat)], levels = treat_levels), treat_id = as.integer(treat),
    predicted = predicted * 100, confLow = confLow * 100, confHigh = confHigh * 100
  )]
  
  anex_full$x = anex_full$barcode_id + x0coord[anex_full$treat]
  data.table::setorder(anex_full, x)
  
  plt = vector("list", (1+length(breaks_plot)))
  for(pl in 1:(1+length(breaks_plot))) { #pl=1
    if(pl > 1L & pl <= length(breaks_plot)) { 
      y_limits = breaks_plot[pl-c(1,0)]
      y_breaks = unique(c(y_limits[1], rev(seq(from = y_limits[2], to = y_limits[1], by = -breaks_y_by[pl]))))
    } else if(pl == 1L) { 
      y_limits = c(min(anex_full$confLow), breaks_plot[1]) 
      y_breaks = rev(seq(from = y_limits[2], to = y_limits[1], by = -breaks_y_by[pl]))
    } else { 
      y_limits = c(breaks_plot[length(breaks_plot)], 100) 
      y_breaks = unique(c(y_limits[1], rev(seq(from = y_limits[2], to = y_limits[1], by = -breaks_y_by[pl]))))
    }
    background_shade$ymin = y_limits[1]
    background_shade$ymax = y_limits[2]
    line_power = power_hlines[power_hlines >= y_limits[1] & power_hlines <= y_limits[2]]
    
    anex_range = anex_full[confHigh > y_limits[1] & confLow < y_limits[2]]
    anex_range$confHigh = pmin(anex_range$confHigh, y_limits[2])
    anex_range$confLow  = pmax(anex_range$confLow,  y_limits[1])
    #plot(anex_range$confLow)
    #plot(anex_range$confHigh)
    anex_range$predicted[anex_range$predicted < anex_range$confLow | anex_range$predicted > anex_range$confHigh] = NA
    
    ggplot() +
      geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = col), data = background_shade) +
      scale_fill_manual(values = c("c1" = "grey95", "c2" = "white")) +
      geom_hline(yintercept = line_power, color = "grey60") +
      geom_point(aes(x = x, y = predicted, colour = treat), size = 1.5, data = anex_range) +
      geom_segment(aes(x = x, y = confLow, yend = confHigh, colour = treat), data = anex_range) + 
      scale_color_discrete(type = colors_plot, drop = F) +
      scale_x_continuous(limits = c(1,nbarcodes+1), breaks = NULL, minor_breaks = NULL, expand = c(0,0)) +
      scale_y_continuous(limits = y_limits, breaks = y_breaks, minor_breaks = NULL, labels = y_breaks) +
      xlab(NULL) + ylab(NULL) +
      theme_minimal() + theme(
        #legend.position = "none", 
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line(colour = "grey95")) -> plt[[pl]]
  }
  
  background_shade_text = background_shade
  background_shade_text$ymin = -1
  background_shade_text$ymax = 1
  ggplot2::ggplot() +
    ggplot2::geom_rect(ggplot2::aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax, fill = col), data = background_shade_text) +
    ggplot2::scale_fill_manual(values = c("c1" = "grey95", "c2" = "white")) +
    ggplot2::geom_text(ggplot2::aes(x = barcode_id + .5, y = 0, label = barcode), 
                       data = unique(anex_full[,.(barcode, barcode_id)]), angle = 90, fontface = "bold", size = size_barcodes) +
    ggplot2::scale_x_continuous(limits = c(1,nbarcodes+1), breaks = NULL, minor_breaks = NULL, expand = c(0,0)) +
    ggplot2::scale_y_continuous(breaks = 0, minor_breaks = NULL, labels = "barcode") +
    ggplot2::xlab(NULL) + ggplot2::ylab(NULL) +
    ggplot2::theme_minimal() + 
    ggplot2::theme(
      legend.position = "none",
      axis.text.y = ggplot2::element_text(angle = 90, hjust = 0.5),
      axis.ticks.y = ggplot2::element_blank()
    ) -> plt_text
  
  plt[[3]] + plt[[2]] + plt[[1]] + plt_text + patchwork :: plot_layout(
    ncol = 1, height = c((1-prop_barcodes) * rep(1/3,3), prop_barcodes), axes = "collect_x"
  ) -> plt
  
  ggplot2::ggsave(file.path(plot_folder, paste0("bar_", experiments[ex], ".", format_plot)), plt, height = height_plot, width = width_plot)
  cat("stored plot:   ", file.path(plot_folder, paste0("bar_", experiments[ex], ".", format_plot)), "\n")
}

