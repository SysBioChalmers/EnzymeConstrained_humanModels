# Script to generate ecGEM plots used in Fig. 5

library(ggplot2)
library(reshape2)

fpath <- '../Results/11_cellLines_NCI60/'
fpath_fva <-  '../Results/FVA/'

# ======================== Comparison of fluxes and growth rates ========================

# load experimental and predicted exchange fluxes at lowest constraint level (defining only media composition)
fluxes <- read.table(paste0(fpath, 'ecModels_const_0_exchangeFluxesComp.txt'), stringsAsFactors=F, header=T)
fluxes <- subset(fluxes, select = -exchangeIDs)  # remove this unneeded column

samples <- colnames(fluxes)
sample.types <- factor(startsWith(samples,'exp_'), labels=c('Predicted','Measured'))

# rearrange data
fluxdata <- NULL
for (i in 1:(ncol(fluxes)/2)) {
  add_data <- fluxes[, c(1, (2*i):(2*i+1))]
  add_data$cell.line <- colnames(fluxes)[2*i+1]
  colnames(add_data) <- c('Metabolite','Measured','Predicted','Cell.Line')
  fluxdata <- rbind(fluxdata, add_data)
}
fluxdata$Cell.Line <- as.factor(fluxdata$Cell.Line)
fluxdata$Metabolite <- as.factor(fluxdata$Metabolite)

# log-transform fluxes
log_fluxdata <- fluxdata
log_fluxdata[, 2:3] <- log10(abs(fluxdata[, 2:3]))
log_fluxdata[, 2:3] <- do.call(data.frame,lapply(log_fluxdata[, 2:3], function(x) replace(x, is.infinite(x), -8)))


# ---------------- Plotting ----------------

cpalette <- c('#f44336','#c2185b','#ab47bc','#4527a0','#3f51b5','#1e88e5',
              '#4db6ac','#2e7d32','#c0ca33','#f9a825','#795548')

# plot all log-transformed fluxes
pdf(paste0(fpath,'exchFlux_compare_ecGEM.pdf'), width=4, height=4.2)
ggplot(log_fluxdata, aes(x=Measured, y=Predicted, color=Cell.Line)) +
  geom_abline(intercept=0, slope=1, alpha=0.5) + 
  geom_point(alpha=0.6, size=2) +
  scale_color_manual(values=cpalette) +
  theme_classic() + 
  theme(text=element_text(size=12),
        legend.position='bottom',
        legend.title=element_blank(),
        axis.text=element_text(size=12, color='black'))
dev.off()


# ======================== Comparison of errors ========================

# load data
ec_error <- read.table(paste0(fpath,'ecModels_error_gRates.txt'), header=T, row.names=1)
reg_error <- read.table(paste0(fpath,'models_error_gRates.txt'), header=T, row.names=1)

# re-arrange data
ec_error <- melt(as.matrix(ec_error), id.vars=NULL)
ec_error$model.type <- 'ecGEM'
reg_error <- melt(as.matrix(reg_error), id.vars=NULL)
reg_error$model.type <- 'GEM'

# merge data
error_data <- rbind(ec_error, reg_error)
colnames(error_data) <- c('cell.lines','constraint.level','relative.error','model.type')

# box-plot of growth rate errors
cpalette <- c('#A4C8D6','#296EA0')  # light blue and blue
pdf(paste0(fpath,'relative_error_boxplots.pdf'), width=3.5, height=2.5)
ggplot(error_data, aes(x=constraint.level, y=relative.error, fill=model.type)) +
  geom_boxplot() +
  scale_fill_manual(values=cpalette) +
  theme_classic() + 
  theme(text=element_text(size=12),
        axis.text=element_text(color='black')) +
  coord_cartesian(ylim=c(0,2), xlim=c(0.5,4.5), expand=F)
dev.off()


# ======================== FVA ========================

# load data
# cell lines: HOP62, HOP92, HS_578T, HT29, MALME_3M, MDMAMB_231, NCI_H226, O_786, RPMI_8226, SR, UO_31
cell_line <- 'HOP62'
filename <- paste0('FVA_comp_',cell_line)
fva <- read.table(paste0(fpath_fva,filename,'.txt'), header=T, sep='\t')
fva <- subset(fva, select=-c(formulas, subSystems, rxns))

# sort FVA ranges and add "fraction" column
fva$model_ranges <- sort(fva$model_ranges)
fva$ecModel_ranges <- sort(fva$ecModel_ranges)
fva$fraction <- seq(0, 1, length.out=nrow(fva))

# melt dataframe for plotting
fva_data <- melt(fva, id.vars='fraction')
colnames(fva_data) <- c('fraction','type','range')
fva_data$range[fva_data$range == 0] <- min(fva_data$range[fva_data$range > 0])
fva_data$range <- log10(fva_data$range)

# plot
cpalette <- c('#105284','#8EB9D0')  # light blue and blue
pdf(paste0(fpath_fva, filename, '.pdf'), width=3.1, height=2.8)
# pdf(paste0(fpath_fva, filename, '.pdf'), width=2.8, height=2.5)
ggplot(fva_data, aes(x=range, y=fraction, group=type, color=type)) +
  geom_path(size=1.2, lineend='round') +
  scale_color_manual(values=cpalette) +
  theme_classic() +
  theme(text=element_text(size=12),
        axis.text=element_text(color='black'),
        legend.position='none') +
  coord_cartesian(ylim=c(0,1.01), xlim=c(-12.5,3.5), expand=F) +
  scale_x_discrete(limits=c(-12, -9, -6, -3, 0, 3))
dev.off()


