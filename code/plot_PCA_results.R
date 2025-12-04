# read PCA results
vec <- read.table("data/plink_pca_1000g/1000GxGDA_pca.eigenvec", header=T, comment.char="")
colnames(vec) <- sub("^X.", "", colnames(vec))

# read 1000G sample data (super-population phenotype)
sample <- read.table("data/plink_pca_1000g/all_phase3.psam", header = T, comment.char="")
colnames(sample) <- sub("^X.", "", colnames(sample))


# create dataframe for plotting
pca_df <- merge(sample, vec, all.y = T) # merge sample characteristics and PCs by IID
head(pca_df)

# subset just samples' GDA PCA results
pca_df$sample <- grepl("^WHS", pca_df$IID) # N = 90
sample_pca <- pca_df[pca_df$sample,]
dim(sample_pca)

# subset just 1000G reference PCA results
ref_pca <- pca_df[!pca_df$sample,]
dim(ref_pca)
colSums(is.na(ref_pca)) # check for missing population data


# read sample ethnic group (demographics)
demog_df <- read.csv("data/demographics.csv")
sample_pca <- merge(demog_df, sample_pca[,-6:-2], by.x = "sample_ID", by.y = "IID")
# sample_pca$ethnic_group <- factor(sample_pca$ethnic_group, levels = c("Chinese", "Malay", "Indian"))


# convert super-population from abbreviation to full name
superpop_abbr_map <- list(
  "AFR" = "AFR: African",
  "AMR" = "AMR: Ad Mixed American",
  "EAS" = "EAS: East Asian",
  "EUR" = "EUR: European",
  "SAS" = "SAS: South Asian"
)
pop_abbr_map <- list(
  "CHB" = "CHB: Han Chinese in Beijing, China",
  "JPT" = "JPT: Japanese in Tokyo, Japan",
  "CHS" = "CHS: Southern Han Chinese",
  "CDX" = "CDX: Chinese Dai in Xishuangbanna, China",

  "KHV" = "KHV: Kinh in Ho Chi Minh City, Vietnam",
  "CEU" = "CEU: Utah European Residents (CEPH)",
  "TSI" = "TSI: Toscani in Italia",
  "FIN" = "FIN: Finnish in Finland",

  "GBR" = "GBR: British in England and Scotland",
  "IBS" = "IBS: Iberian Population in Spain",
  "YRI" = "YRI: Yoruba in Ibadan, Nigeria",
  "LWK" = "FIN: Luhya in Webuye, Kenya",

  "GWD" = "GWD: Gambian in Western Divisions in the Gambia",
  "MSL" = "MSL: Mende in Sierra Leone",
  "ESN" = "ESN: Esan in Nigeria",
  "ASW" = "ASW: Americans of African Ancestry in SW USA",

  "ACB" = "ACB: African Caribbeans in Barbados",
  "MXL" = "MXL: Mexican Ancestry from Los Angeles USA",
  "PUR" = "PUR: Puerto Ricans from Puerto Rico",
  "CLM" = "CLM: Colombians from Medellin, Colombia",

  "PEL" = "PEL: Peruvians from Lima, Peru",
  "GIH" = "GIH: Gujarati Indian from Houston, Texas",
  "PJL" = "PJL: Punjabi from Lahore, Pakistan",
  "BEB" = "BEB: Bengali from Bangladesh",

  "STU" = "STU: Sri Lankan Tamil from the UK",
  "ITU" = "ITU: Indian Telugu from the UK"
)

table(ref_pca$SuperPop)
table(ref_pca$Population)

ref_pca$SuperPop <- as.factor(unlist(superpop_abbr_map[ref_pca$SuperPop], use.names = F))
ref_pca$Population <- as.factor(unlist(pop_abbr_map[ref_pca$Population], use.names = F))


# PCA plot
plot_pca <- function(ref_pca, sample_pca, pc = 1:2, xlims = NULL, ylims = NULL, legend_pos = "right") {
  pc <- sapply(pc, function(x) paste0("PC", x))
  p <- ggplot() +
    geom_point(data = ref_pca, mapping = aes(x = get(pc[1]), y = get(pc[2]), colour = SuperPop), shape = 18, size = 2, alpha = 0.5) +
    # scale_colour_brewer(palette = "Set2", name = "1000G super-population") +
    scale_colour_manual(values = brewer.pal(name="Set2", n=8)[c(1, 8, 3, 4, 2)], name = "1000G super-population") +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
    new_scale_colour() +
    # geom_point(data = sample_pca, mapping = aes(x = get(pc[1]), y = get(pc[2]), shape = ethnic_group), colour = "grey20", size = 2.5, alpha = 1, stroke = 1.05) +
    geom_point(data = sample_pca, mapping = aes(x = get(pc[1]), y = get(pc[2]), shape = ethnic_group, colour = ethnic_group), size = 2.2, alpha = 0.9, stroke = 1) +
    scale_shape_manual(values = 2:0, name = "Sample ethnic group", breaks=c("Chinese", "Malay", "Indian")) +
    scale_colour_brewer(palette = "Set1", name = "Sample ethnic group", breaks=c("Chinese", "Malay", "Indian")) +
    labs(x = pc[1], y = pc[2]) +
    theme_bw() + theme(aspect.ratio = 1, legend.position = legend_pos)

  if (!is.null(xlims)) p <- p + xlim(xlims[1], xlims[2])
  if (!is.null(ylims)) p <- p + ylim(ylims[1], ylims[2])

  return(p)
}

# PC1 vs PC2
p1 <- plot_pca(ref_pca, sample_pca, 1:2)
p1_zoom <- plot_pca(ref_pca, sample_pca, 1:2, c(NA, -0.005), c(NA, 0.02))

# PC1 vs PC3
p2 <- plot_pca(ref_pca, sample_pca, c(1, 3))
p2_zoom <- plot_pca(ref_pca, sample_pca, c(1, 3), c(NA, -0.005), c(NA, 0.01))


# arrange subplots
NUM_SUB <- 2
panels <- plot_grid(
  p1 + theme(legend.position = "none"), NULL,
  p2 + theme(legend.position = "none"),
  p1_zoom + theme(legend.position = "none"), NULL,
  p2_zoom + theme(legend.position = "none"),
  nrow = 2,
  rel_widths = rep(rep(c(0.025, 1), NUM_SUB)[-1], 2),
  labels = unlist(lapply(1:2, function(r) sapply(letters[((r-1)*NUM_SUB+1):(r*NUM_SUB)], function(l) c("", paste0("(", l, ")")))[-1])),
  label_size = 12, hjust = -0.05, vjust = 2.3 # shift label left
)

legend <- get_legend(
  p1 + theme(
    legend.box.margin = margin(24, 0, 0, 0), # top spacer
    legend.justification = "top"
  )
)


# save plot
SIZE_UNIT <- 4
rw <- c(0.05, NUM_SUB, 0.5, 0.05)
png(file = sprintf("figures/fig1-pca_ancestry_1000g_superpop.png", NUM_SUB), width = sum(rw) * SIZE_UNIT, height = 2 * SIZE_UNIT, units = "in", res = 600)
plot_grid(
  NULL, panels, legend, NULL,
  align = "v", nrow = 1,
  rel_widths = rw
)
dev.off()
