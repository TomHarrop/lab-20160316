library(data.table)
library(DESeq2)

tpm <- readRDS("data/tpm_with_calls.Rds")

### CROP PHOTOS

png.files <- c("data/N1r.png", "data/N2-1.png", "data/N3-1.png", "data/N4-1.png")
CropAndRaster <- function(png.file) {
  tmp <- tempfile(fileext = ".png")
  system(paste("convert", png.file,
               "-resize 600x600^ -gravity center -background transparent",
               "-crop  600x600+0+0 +repage", tmp))
  png::readPNG(tmp)
}
pngs <- lapply(png.files, CropAndRaster)

### PANICLE PHENOTYPES

# mung data
raw.data <- data.table(read.csv('data/Phenotype_Panicle_corrected.csv',
                                stringsAsFactors = FALSE))
raw.data <- raw.data[!grepl("^Echelle", file_name) ]
raw.data[, Accession := gsub("^([[:alnum:]]+).*", "\\1", file_name)]
raw.data[Accession == 'tog5681', Accession := "Tog5681"]
raw.data[Accession == 'Nip', Accession := "Nipponbare"]
raw.data[, Accession := factor(
  Accession, levels = c('IR64', 'Nipponbare', 'W1654','Tog5681', 'B88'))]
raw.data[, Plant.no := gsub(
  "^[[:alnum:]]+_([[:digit:]]+)_.*", "\\1", file_name)]
raw.data[, Panicle.no := gsub(
  "^[[:alnum:]]+_[[:digit:]]+_([[:digit:]]+)_.*", "\\1", file_name)]
raw.data[, Species := plyr::mapvalues(
  Accession,
  from = c("Nipponbare", "IR64", "B88", "Tog5681", "W1654"),
  to = c("O. sativa", "O. sativa", "O. barthii", "O. glaberrima",
         "O. rufipogon"))]
setkey(raw.data, "file_name")

# set up plot
plot.data.wide <- unique(raw.data[, .(
  name = file_name,
  Accession, Plant.no, Panicle.no, Species,
  Spikelets = Sp_nb, `Secondary branches` = TA_nb)])
plot.data.phenotypes <- reshape2::melt(plot.data.wide,
                                       id.vars = c("name", "Accession", "Plant.no",
                                                   "Panicle.no", "Species"))

## MAPPING STATS

star.logs <- readRDS('data/starLogs.Rds')
setkey(star.logs, "Library")
tpm.with.calls <- readRDS('data/tpm_with_calls.Rds')

genes.per.sample <- tpm.with.calls[, .(
  `Genes detected` = sum(call)),
  by = sample]
setkey(genes.per.sample, "sample")

plot.data.stats.wide <- genes.per.sample[star.logs,.(
  Sample = sample, 
  `Input reads (M)` = `Number of input reads` / 1e6,
  `Uniquely mapped reads (M)` = `Uniquely mapped reads number` / 1e6,
  `Reads in genes (M)` = `Number of reads in genes` / 1e6,
  `Genes detected / 1000` = `Genes detected`/1e3
)]

plot.data.stats.wide[, Accession := factor(plyr::mapvalues(
  substr(Sample, 1, 1), from = c("I", "J", "R", "G", "B"),
  to = c("IR64", "Nipponbare", "W1654", "Tog5681", "B88")),
  levels = c("IR64", "Nipponbare", "W1654", "Tog5681", "B88"))]

plot.data.stats <- reshape2::melt(plot.data.stats.wide,
                                  id.vars = c("Accession", "Sample",
                                              "Input reads (M)"))

### WALD SPECIES

dds.species <- readRDS("data/dds.species.Rds")
species.contrast.results <- readRDS("data/species_contrast_results.Rds")

# summary stats
species.contrast.results[log2FoldChange > 0, direction := "up"]
species.contrast.results[log2FoldChange < 0, direction := "down"]
species.contrast.results[, s1 := gsub("\\..*", "", contrast)]
species.contrast.results[, s2 := gsub(".*\\.", "", contrast)]
species.contrast.results[, s1 := plyr::mapvalues(
  s1, from = c("indica", "japonica", "rufipogon", "glaberrima", "barthii"),
  to = c("IR64", "Nipponbare", "W1654", "Tog5681", "B88"))]
species.contrast.results[, s2 := plyr::mapvalues(
  s2, from = c("indica", "japonica", "rufipogon", "glaberrima", "barthii"),
  to = c("IR64", "Nipponbare", "W1654", "Tog5681", "B88"))]
species.contrast.results[, contrast := paste(s1, s2, sep = "–")]

species.contrast.results[, .(
  `padj < 0.05` = sum(padj < 0.05, na.rm = TRUE),
  `L2FC > 0` = sum(padj < 0.05 & log2FoldChange > 0, na.rm = TRUE),
  `L2FC < 0` = sum(padj < 0.05 & log2FoldChange < 0, na.rm = TRUE)
),
by = contrast]

# illustrate some genes
IllustratePlot <- function(x, contrast){
  plot.data <- tpm[x,]
  plot.data[, gene.name := oryzr::LocToGeneName(gene)$symbols, by = gene]
  plot.data[is.na(gene.name), gene.name := gene]
  g <- ggplot(plot.data, aes(x = stage, y = tpm, colour = species,
                             group = species)) +
    theme_grey(base_size = 6, base_family = "Lato") +
    theme(plot.background = element_rect(fill = "transparent", colour = NA),
          legend.background = element_rect(fill = "transparent", colour = NA)) +
    scale_colour_brewer(palette = "Set1", guide=FALSE) +
    ggtitle(plot.data[, paste(gene.name, contrast, sep = " · ")]) +
    ylab(NULL) + xlab(NULL) +
    geom_smooth(method = "lm", se = FALSE, size = 0.5, alpha = 0.5) +
    geom_point(position = position_jitter(width = 0.2))
  ggplotGrob(g)
}

illustrative.genes <- species.contrast.results[!is.na(direction),
                                               .SD[which.min(padj)],
                                               by = .(contrast, direction)]

# list of plots
plots <- lapply(1:dim(illustrative.genes)[1], function(i)
  illustrative.genes[i, IllustratePlot(x = gene, contrast = contrast)])

### WALD DOMESTICATION ###
dds.domestication <- readRDS('data/dds_domestication.Rds')
results.table.domestication <- readRDS(
  'data/wald_domestication_results_table.Rds')

sig.domestication <- results.table.domestication[padj < 0.05]
setkey(sig.domestication, "gene")

sig.domestication.tpm <- tpm[sig.domestication, .(
  gene, sample, call, tpm, species, stage,  padj
)]

sig.domestication.tpm[, Accession := factor(plyr::mapvalues(
  species,
  from = c("I", "J", "R", "G", "B"),
  to = c("IR64", "Nipponbare", "W1654", "Tog5681", "B88")),
  levels = c("IR64", "Nipponbare", "W1654", "Tog5681", "B88"))]
sig.domestication.tpm[, gene.name :=
                        oryzr::LocToGeneName(gene)$symbols, by = gene]
sig.domestication.tpm[is.na(gene.name), gene.name := gene]
sig.domestication.tpm[, gene.name := paste0(gene.name, ", p = ", round(padj, 4))]
sig.domestication.tpm[, gene.name :=
                        factor(gene.name,
                               levels = unique(gene.name[order(padj)]))]

### DOMESTICATION BY CONTINENT
dds.domestication.by.continent <- readRDS(
  'data/dds_domestication_by_continent.Rds')
results.table.dbc <- readRDS(
  "data/wald_domestication_by_continent_results_table.Rds")
setkey(results.table.dbc, "gene") 
dbc.tpm <- tpm[results.table.dbc, .(
  gene, sample, call, tpm, species, stage, domestication, padj
)]
dbc.tpm[, Accession := factor(plyr::mapvalues(
  species,
  from = c("I", "J", "R", "G", "B"),
  to = c("IR64", "Nipponbare", "W1654", "Tog5681", "B88")),
  levels = c("IR64", "Nipponbare", "W1654", "Tog5681", "B88"))]

dbc.asia <- dbc.tpm[domestication == "asia" & padj < 0.05 &
                      species %in% c("R", "I", "J")]
dbc.asia[, gene.name := oryzr::LocToGeneName(gene)$symbols, by = gene]
dbc.asia[is.na(gene.name), gene.name := gene]
dbc.asia[, gene.name := paste0(gene.name, ", p = ", round(padj, 4))]
dbc.asia[, gene.name :=
           factor(gene.name, levels = unique(gene.name[order(padj)]))]

dbc.africa <- dbc.tpm[domestication == "africa" & padj < 0.05 &
                      species %in% c("B", "G")]
dbc.africa[, gene.name := oryzr::LocToGeneName(gene)$symbols, by = gene]
dbc.africa[is.na(gene.name), gene.name := gene]
dbc.africa[, gene.name := paste0(gene.name, ", p = ", round(padj, 4))]
dbc.africa[, gene.name :=
           factor(gene.name, levels = unique(gene.name[order(padj)]))]
 
### ERF GENES

sig.erfs <- c('LOC_Os04g52090', 'LOC_Os05g41760', 'LOC_Os05g41780',
              'LOC_Os08g31580')
erfs.tpm <- tpm[sig.erfs]
erfs.tpm[, Accession := factor(plyr::mapvalues(
  species,
  from = c("I", "J", "R", "G", "B"),
  to = c("IR64", "Nipponbare", "W1654", "Tog5681", "B88")),
  levels = c("IR64", "Nipponbare", "W1654", "Tog5681", "B88"))]
erfs.tpm[, gene.name := oryzr::LocToGeneName(gene)$symbols, by = gene]
erfs.tpm[is.na(gene.name), gene.name := gene]