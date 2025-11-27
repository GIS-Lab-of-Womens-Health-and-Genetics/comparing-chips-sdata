library(tidyverse)

# define constants
arr_prefix <- c("Axiom", "GDA", "GSA", "OncoArray")
arr_size <- c(581424, 1904599, 654027, 499170)
names(arr_size) <- arr_prefix

# helper function to read stats outputs
get_arr_df <- function(arr, ext_type, colname){
  filename <- paste0(arr, ext_type)
  df <- read.table(filename, header = TRUE)
  df$arr <- arr
  return(df[,c("arr", colname)])
}

# plot histograms, faceted by array
grid_histogram <- function(arr_prefix, ext_type, colname, xlab, bins, xticks, inv = FALSE, scales = "fixed", threshold = NULL, annot_left = TRUE, log_trf = FALSE, filename = NULL){
  df <- bind_rows(lapply(arr_prefix, get_arr_df, ext_type, colname))
  if (inv) df[,colname] <- 1 - df[,colname]
  df$arr <- factor(df$arr, levels = arr_prefix) # retain order of arrays
  # print(summary(df[[colname]]))
  arr_stats <- df %>% group_by(arr) %>% summarize(
    min = min(get(colname), na.rm = TRUE),
    Q1 = quantile(get(colname), 0.25, na.rm = TRUE),
    med = median(get(colname), na.rm = TRUE),
    Q3 = quantile(get(colname), 0.75, na.rm = TRUE),
    max = max(get(colname), na.rm = TRUE),
    mean = mean(get(colname), na.rm = TRUE),
    sd = sd(get(colname), na.rm = TRUE)
  )

  p <- ggplot(df, aes(x = get(colname))) + theme_bw() + xlab(xlab) + ylab("")
  if (!is.null(threshold)) {
    df$thres_mask <- (df[[colname]] > threshold) # 0: false, 1: true
    p <- p + geom_histogram(aes(fill = as.factor(df$thres_mask)), breaks=bins) +
      geom_vline(xintercept = threshold, colour = "red2", linetype="solid", linewidth = 0.7) +
      scale_fill_manual(
        values = c("tomato", "darkslategrey"),
        labels = c("Excluded", "Retained after QC"),
        name = NULL
      ) + theme(legend.position = "bottom")

    arr_stats <- merge(
      arr_stats,
      df %>% group_by(arr) %>% summarize(
        num_ret = sum(thres_mask, na.rm = TRUE),
        num_total = n()
      ), by = "arr"
    )
    arr_stats$annot_label <- with(arr_stats, paste0(
      format(num_ret, big.mark=",", scientific=FALSE),
      " above ", threshold,
      " (", round(num_ret/num_total*100, 1), "%", " of ",
      format(num_total, big.mark=",", scientific=FALSE), ")"
    ))
  } else{
    p <- p + geom_histogram(breaks=bins, fill = "darkslategrey")

    arr_stats$annot_label <- with(arr_stats, paste0(
      "Median ", format(med, digits = 3, big.mark=",", scientific=FALSE),
      " (IQR ", format(Q1, digits = 3, big.mark=",", scientific=FALSE), " - ",
      format(Q3, digits = 3, big.mark=",", scientific=FALSE), "), ",
      "Mean ", format(mean, digits = 3, big.mark=",", scientific=FALSE),
      " (SD ", format(sd, digits = 3, big.mark=",", scientific=FALSE), ")"
    ))
  }

  print(arr_stats)
  p <- p + facet_grid(rows = vars(arr), scales = scales) +
    geom_label(
      data=arr_stats,
      aes(
        x = if_else(annot_left, -Inf, Inf),
        y = Inf,
        label= annot_label
      ),
      vjust = "top",
      hjust = if_else(annot_left, "left", "right"),
      # hjust = -0.05, vjust = 1.2,
      alpha = 0.3, size = 3, label.r = unit(0, "lines"),
      inherit.aes=FALSE
    )

  if (log_trf) p <- p + scale_y_log10() #+ guides(y = "axis_logticks")
  p <- p + scale_x_continuous(breaks=xticks)

  if (!is.null(filename)) {
    png(file = filename, width = 7, height = 5, units = "in", res = 600)
    dev.off()
  } else print(p)
}


# sample-level call rate
grid_histogram(arr_prefix, ".imiss", "F_MISS",
               xlab = "Sample-level call rate",
               bins = seq(0.97, 1, by = 0.001),
               xticks = seq(0.97, 1, by = 0.005),
               inv = TRUE, # 1 - missing rate
               filename = "figures/fig3-sample_level_call_rate.png"
)


# variant-level call rate
grid_histogram(arr_prefix, ".lmiss", "F_MISS",
               xlab = "Variant-level call rate",
               bins = seq(0, 1, by = 0.01),
               xticks = seq(0, 1, by = 0.1),
               # scales = "free_y",
               threshold = 0.95, inv = TRUE,
               log_trf = TRUE,
               filename = "figures/fig4-variant_level_call_rate.png"
)


# minor allele frequency (note: some variants MAF = NA, ignored in plot)
grid_histogram(arr_prefix, ".frq", "MAF",
               xlab = "Minor allele frequency",
               bins = seq(0, 0.5, by = 0.005),
               xticks = seq(0, 0.5, by = 0.05),
               scales = "free_y",
               threshold = 0.01,
               annot_left = FALSE,
               filename = "figures/fig5-minor_allele_freq.png"
)
