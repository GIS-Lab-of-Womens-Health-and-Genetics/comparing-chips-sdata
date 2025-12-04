# admixture analysis using radmixture package

#### 1. General set up
# check R version, must be < 4.3.0 for radmixture compatibility
# otherwise, install from CRAN archives and change version in RStudio: Tools > Global options > General > R version
getRversion()

# check renv
# may need to run renv::restore() on first time use
renv::status()

library(tidyverse)
library(radmixture)
library(RColorBrewer)
library(cowplot)
library(ggnewscale)

# initialise data constants
geno_data_dir <- "data/genotyping_arrays"
arr_prefix <- c("PreCIPV1", "GDA", "GSA", "OncoArray")
arr_size <- c(581424, 1904599, 654027, 499170)
names(arr_size) <- arr_prefix

file_suffix <- "_2025-02-17_clean"

pt_withdrawn <- c("WHS408", "WHS526")
ctrl_cell_lines <- c("WHS427", "WHS428", "WHS452", "WHS453", "WHS502", "WHS503", "WHS527", "WHS528")


#### 2. Set up radmixture
# reformat plink binary data into radmixture compatible format
map_table <- read.table(
  file.path(geno_data_dir, "GDA_QC_recode.map"),
  col.names = c(
    "chr", # Chromosome code. PLINK 1.9 also permits contig names here, but most older programs do not.
    "varID", # Variant identifier
    "pos_cm", # Position in morgans or centimorgans (optional; also safe to use dummy value of '0')
    "bp_coords" # Base-pair coordinate
  )
)
head(map_table)
summary(map_table)

# radmixture reference datasets
radmixture_github_data_url <- "https://github.com/wegene-llc/radmixture/raw/master/data"
local_radmixture_dir <- "data/radmixture_datasets"
radmixture_config <- rbind(
  list(
    dataset = "K7b",
    allele_map = "K7b.alleles",
    allele_freq = "K7b.7.F",
    num_pop = 7
  )
  # add other datasets here if desired (e.g. K12b, world9, globe13)
)
rownames(radmixture_config) <- radmixture_config[,"dataset"]


#### 3. Define functions
# radmixture functions
init_radmixture_file <- function(filename) {
  local_filepath <- file.path(local_radmixture_dir, filename)
  if (!file.exists(local_filepath)) {
    download.file(
      url = paste0(radmixture_github_data_url, "/", filename),
      destfile = local_filepath
    )
  }
  load(local_filepath)
}

init_radmixture_dataset <- function(dataset) {
  dataset_config <- radmixture_config[dataset,]
  init_radmixture_file(paste0(dataset_config$allele_map, ".RData"))
  init_radmixture_file(paste0(dataset_config$allele_freq, ".RData"))
  return(dataset_config)
}

process_line <- function(line, map_table) {
  line_data <- unlist(strsplit(line, "\\s", fixed = FALSE)) # split line by whitespace

  pt_data <- as.list(line_data[1:6])
  names(pt_data) <- c(
    "FID", # Family ID ('FID')
    "IID", # Within-family ID ('IID'; cannot be '0')
    "fatherIID", # Within-family ID of father ('0' if father isn't in dataset)
    "motherIID", # Within-family ID of mother ('0' if mother isn't in dataset)
    "sex", # Sex code ('1' = male, '2' = female, '0' = unknown)
    "phenotype" # Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)
  )
  pt_data$sex <- factor(pt_data$sex, levels = 1:2, labels = c("male", "female"))

  # remove control cell lines / withdrawn consent
  if (pt_data[["FID"]] %in% c(ctrl_cell_lines, pt_withdrawn)) gene_df <- NULL

  # 2 x allele calls for each variant in the .map file ('0' = no call);
  gene_data <- sapply(1:nrow(map_table), function(i) paste(line_data[6 + (2*i-1):(2*i)], collapse = ''))
  names(gene_data) <- map_table$varID

  # radmixture input table: var, chr, pos, genotype
  gene_df <- map_table[,c("varID", "chr", "bp_coords")]
  row.names(gene_df) <- gene_df$varID
  gene_df$geno <- gene_data
  row.names(gene_df) <- NULL

  return(list(pt = pt_data, gene = gene_df))
}

get_radmixture_anc <- function(genotype, dataset) {
  dataset_config <- init_radmixture_dataset(dataset)
  res <- tfrdpub(genotype, dataset_config$num_pop, get(dataset_config$allele_map), get(dataset_config$allele_freq))
  ances <- fFixQN(res$g, res$q, res$f, tol = 1e-4, method = "BR", pubdata = dataset)
  return(ances$q)
}

# main function - compare radmixture results to ethnic group
get_anc_res <- function(ped_filepath, demog_filepath, dataset, verbose = FALSE) {
  ped_con <- file(ped_filepath, "r")
  anc_res <- NULL
  while ( TRUE ) {
    line <- readLines(ped_con, n = 1)
    if (length(line) == 0) break # end of file

    data <- process_line(line, map_table)
    if (!is.null(data$gene)) {
      if (verbose) cat("Calculating ancestry for pt", data$pt[["IID"]], "...\n")
      anc <- get_radmixture_anc(data$gene, dataset)
      if (is.null(anc_res)) {
        anc_res <- data.frame(c(data$pt, anc))
      } else {
        anc_res[nrow(anc_res)+1,] <- c(data$pt, anc)
      }
    } else {
      cat("Skipping pt", data$pt[["IID"]], "...\n")
    }
  }
  close(ped_con)

  demog_df <- read.csv(demog_filepath)
  anc_res <- merge(
    within(demog_df, rm("sex")),
    anc_res,
    by.x = "sample_ID", by.y = "FID",
    all.y = TRUE # keep ancestry of control cell lines, with no demographic data
  )
  anc_res$control <- anc_res$sample_ID %in% ctrl_cell_lines

  res_filepath <- file.path("data/radmixture_results", paste0(dataset, "_radmixture.csv"))
  write.csv(anc_res, res_filepath, row.names = FALSE)

  return(res_filepath)
}

# plot functions
sort_anc_cat <- function(anc_cat) {
  # Asian categories first, then other categories alphabetically
  return(anc_cat[order(!grepl("Asian", anc_cat), anc_cat)])
}

plot_anc_res <- function(anc_res, dataset) {
  # sort by self-ident ethnic group then radmixture values (asian cols first)
  anc_col_seq <- sort_anc_cat(colnames(anc_res)[9:(ncol(anc_res)-1)])
  anc_res <- anc_res %>%
    dplyr::filter(!control) %>%
    mutate(ethnic_group = factor(
      ethnic_group,
      levels = c("Chinese", "Malay", "Indian")
    )) %>%
    select(1:8, anc_col_seq) %>%
    group_by(ethnic_group) %>%
    arrange(across(anc_col_seq, desc), .by_group = TRUE)
  anc_res$plot_seq <- as.numeric(row.names(anc_res))

  # prepare data for plotting
  plot_data <- pivot_longer(anc_res, anc_col_seq, names_to = "anc_cat", values_to = "admixture")
  plot_data$anc_cat <- gsub("_", " ", plot_data$anc_cat)

  # prepare annotations for self-reported ethnic group
  ethnic_group_annot <- pivot_longer(
    anc_res %>% group_by(ethnic_group) %>% summarise_at(
      8:(ncol(anc_res)-2),
      mean
    ),
    anc_col_seq,
    names_to = "anc_cat",
    values_to = "mean_prop"
  ) %>%
    mutate(anc_cat = gsub("_", " ", anc_cat)) %>%
    arrange(desc(mean_prop)) %>%
    group_by(ethnic_group) %>%
    slice(1:3) %>%
    summarize(annot = paste(anc_cat, round(mean_prop, 3), collapse = "; \n")) %>%
    mutate(
      start_x = anc_res[!duplicated(anc_res$ethnic_group),]$plot_seq, # 1st entry of each self-reported ethnic group
      end_x = c(start_x[-1] - 1, nrow(anc_res)),
      mid_x = floor((end_x - start_x)/2) + start_x
    )

  # plot
  p <- ggplot(
    plot_data,
    aes(x = plot_seq, y = admixture, fill = factor(anc_cat, levels = gsub("_", " ", anc_col_seq)))
  ) + geom_bar(position = "stack", stat = "identity", just = 0, width = 0.9) +
    geom_vline(xintercept = ethnic_group_annot$start_x[-1], colour= "grey30", linetype = "solid") +
    annotate(
      geom = "text",
      x = ethnic_group_annot$mid_x, hjust = 0.4,
      y = -0.03, vjust = 0.5,
      label = ethnic_group_annot$ethnic_group,
      size = 3.5,
      lineheight = 1
    ) +
    labs(
      x = "Samples, grouped by self-reported ethnic group",
      y = "Ancestry contribution",
      fill = "Ancestry population"
    ) +
    scale_x_discrete(breaks = NULL, expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(-0.06, 1), labels = scales::label_percent(), expand = c(0, 0), oob = scales::squish) +
    scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Set2"))(length(anc_col_seq))) +
    theme_bw() + theme(axis.title.x = element_text(margin = margin(t = 5)))

  png(file = file.path("figures/admixture_", paste0(dataset, ".png")), width = 9, height = 5, units = "in", res = 600)
  print(p)
  dev.off()

  print(ethnic_group_annot)
}


#### 4. Call functions for admixture analysis + plot resullts
for (ds in rownames(radmixture_config)) {
  res_filepath <- get_anc_res(
    ped_filepath = file.path(geno_data_dir, "GDA_QC_recode.ped"),
    demog_filepath = "data/demographics.csv",
    dataset = ds
  )
  anc_res <- read.csv(res_filepath)

  plot_anc_res(anc_res, dataset = ds)
}
