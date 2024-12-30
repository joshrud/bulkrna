library(Rtsne)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(gplots)
library(edgeR)
library(tidyverse)
library(BiocParallel)
library(reshape2)
library(ggrepel)
library(fgsea)

scatter_corr <- function(v1, 
                         v2, 
                         n1=NULL,
                         n2=NULL) {
  
  cor.result <- cor.test(v1,v2)
  cur.df <- data.frame(v1,v2)
  if (!is.null(n1)) {
    colnames(cur.df)[1] <- n1
  }
  if (!is.null(n2)) {
    colnames(cur.df)[2] <- n2
  }
  
  # get the names of the cols, they could be V1/V2 or n1/n2
  col1.n <- colnames(cur.df)[1]
  col2.n <- colnames(cur.df)[2]
  
  # plot 
  ggplot(data=cur.df, aes_string(x=col1.n,
                                 y=col2.n)) + 
    geom_point() + 
    geom_smooth(formula = y~x, se=F, method = "lm") + 
    # geom_abline(slope=1, intercept = 0) + 
    labs(title = paste0("Pearson  R: ", round(cor.result$estimate, 3))) + 
    theme_bw()
  
}




# copied from here becuase I'm dumbb: https://stackoverflow.com/questions/35194048/using-r-how-to-calculate-the-distance-from-one-point-to-a-line
# uses three points: a, b, c of which each are concatenated arrays, ex: a=c(1,2), b=c(2,3), c=c(3,4)
# caclculates the distance of the point a to the line formed by points b and c
dist2d <- function(a,b,c) {
  v1 <- b - c
  v2 <- a - b
  m <- cbind(v1,v2)
  d <- abs(det(m))/sqrt(sum(v1*v1))
} 


#' personalis_boxplot
#' @description For simply making a ggplot summary plot  \
#' only with median instead of boxplot
#' 
#' @param data the sampleform usually  
#' @param groupName column in sampleform that will be used for grouping
#' @param valueName column in sampleform that will be used for grabbing y axis
#' @param colorList a named list for the grouping in data$groupName \
#'  if we want a specific order
#' @return the ggplot object 
personalis_boxplot <- function(data, 
                               groupName, 
                               valueName, 
                               colorList=NULL) {
  require(dplyr)
  require(ggplot2)
  
  indorder <- data %>% 
    group_by(groupname) %>% 
    summarise(stat=median(valuename)) %>%
    arrange(desc(stat))
  sampleform$groupname <- 
    factor(sampleform$groupname, 
           levels = indorder$groupname)
  
  p1 <- ggplot(data=sampleform, aes(x=groupName, 
                                    y=valueName, 
                                    color=groupName)) +
    geom_jitter(position=position_jitterdodge(jitter.width = 0.15,
                                              jitter.height = 0, 
                                              dodge.width = 0.75)) +
    stat_summary(
      mapping = aes(x = groupName, 
                    y = valueName, 
                    color = groupName),
      fun.min = function(a) { median(a) },
      fun.max = function(a) { median(a) }, # if we remove then we get warnings...
      fun = median,
      geom = "pointrange",
      position = position_dodge(width = 0.75), 
      shape = 95, size=0.5, fatten = 15) + 
    theme_bw()+
    theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))
  
  # add the colors in the way we want if it was supplied to the function  
  if (!is.null(colorList)) {
    p1 <- p1 + scale_color_manual(values=colorList)
  }
  
  return(p1)
}

plotDists = function (vsd.obj) {
  sampleDists <- dist(t(assay(vsd.obj)))

  sampleDistMatrix <- as.matrix( sampleDists )
  rownames(sampleDistMatrix) <- paste( vsd.obj$celltype )

  # colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  colors <- rev(viridis_pal(option = "D")(255))
  print(pheatmap(sampleDistMatrix,
                 clustering_distance_rows = sampleDists,
                 clustering_distance_cols = sampleDists,
                 col = colors,
                 border_color = NA))
}

# joshR's modified version which uses prcomp instead of the DESeq2 implemention of pca
name.plotPCA = function (obj, topn = NULL, coi="condition", label=T) {
  if (length(coi) != 1 && length(coi) != 2) {
    print("coi needs to be either a column name in vsd.obj coldata, or 2 of them")
    return(NULL)
  }
  targets <- as.data.frame(obj@colData)
  if (!is.null(topn)) {
    a <- assay(obj)
    a <- head(a[order(rowVars(a), decreasing = T),], n=topn) #use topn genes
  } else {
    a <- assay(obj)
    pca <- prcomp(t(a)) #use all genes
  }
  pca <- prcomp(t(a))
  pca.df <- as.data.frame(pca$x[,c(1,2)]) #first 2 columns are PC1 and PC2
  if (length(coi) > 1) {
    if (!is.factor(targets[,coi[1]])) {
      targets[,coi[1]] <- as.factor(targets[,coi[1]])
    }
    if (!is.factor(targets[,coi[2]])) {
      targets[,coi[2]] <- as.factor(targets[,coi[2]])
    }
    pca.df[,coi[1]] <- targets[,coi[1]]
    pca.df[,coi[2]] <- targets[,coi[2]]
  } else {
    pca.df[,coi] <- targets[,coi]
  }
  pca.df$samples <- rownames(pca.df)
  labels <- as.vector(pca$sdev^2/sum(pca$sdev^2))
  pc1.label <- labels[1]
  pc2.label <- labels[2]
  if (length(coi) > 1) {
    p <- ggplot(pca.df, aes(x = PC1, y = PC2)) +
      geom_point(aes_string(color=coi[1], shape=coi[2]), size=3) +
      xlab(paste0("PC1: ", as.character(round(pc1.label, 3)*100), "% variance")) +
      ylab(paste0("PC2: ", as.character(round(pc2.label, 3)*100), "% variance")) +
      theme_bw()
  } else {
    p <- ggplot(pca.df, aes(x = PC1, y = PC2)) +
      geom_point(aes_string(color=coi), size=3) +
      xlab(paste0("PC1: ", as.character(round(pc1.label, 3)*100), "% variance")) +
      ylab(paste0("PC2: ", as.character(round(pc2.label, 3)*100), "% variance")) +
      theme_bw()
  }
  if (label) {
    p <- p + geom_text_repel(aes_string(label = "samples"),
                            color = "black")
  }
  return(p)
}

write.pcatable = function(obj, topn=3000) {
  a <- assay(obj)
  a <- head(a[order(rowVars(a), decreasing = T),], n=topn) #use topn genes
  pca <- prcomp(t(a))
  pca.rot <- pca$rotation
  pca.rot <- merge(annotation, pca.rot[,1:5], by.x="Gene.stable.ID", by.y = "row.names")
  rownames(pca.rot) <- pca.rot$Gene.stable.ID
  pca.rot <- pca.rot[,-1]
  pca.rot <- pca.rot[order(abs(pca.rot[,grep("PC1", colnames(pca.rot))]), decreasing = T),]
  write.table(pca.rot, file="pca_rotations.tsv", quote = F, sep = "\t", col.names = NA)
}

# for plotting tsne (using input perplexity)
plotTSNE = function (obj, perp, coi) {
  targets <- as.data.frame(obj@colData)
  set.seed(777)
  tsne <- Rtsne(t(assay(obj)), perplexity = perp)
  tsne.df <- as.data.frame(tsne$Y)
  colnames(tsne.df) <- c("tSNE_1", "tSNE_2")
  tsne.df$samples <- rownames(targets)
  if (length(coi) > 1) {
    tsne.df[,coi[1]] <- targets[,coi[1]]
    tsne.df[,coi[2]] <- targets[,coi[2]]
  } else {
    tsne.df[,coi] <- targets[,coi]
  }
  if (length(coi) > 1) {
    p <- ggplot(tsne.df, aes_string(x="tSNE_1", y="tSNE_2", color=coi[1], shape=coi[2])) + geom_point() +
      geom_text_repel(aes_string(label = "samples"),
                      color = "black") +
      labs(title = "tSNE Colored by Condition") +
      theme(legend.title = element_text(size = 15),
            legend.text = element_text(size = 12),
            axis.text.x = element_text(size=12),
            axis.text.y = element_text(size = 12))
  } else {
    p <- ggplot(tsne.df, aes_string(x="tSNE_1", y="tSNE_2", color=coi)) + geom_point() +
      geom_text_repel(aes_string(label = "samples"),
                      color = "black") +
      labs(title = "tSNE Colored by Condition") +
      theme(legend.title = element_text(size = 15),
            legend.text = element_text(size = 12),
            axis.text.x = element_text(size=12),
            axis.text.y = element_text(size = 12))
  }

  print(p)
}


#' plotHeatmap
#' @description Plots a QC heatmap using top n variable genes
#' @param vsd.obj the stabilized object from deseq2
#' @param with_labels whether to show labels
#' @param TOPN number of variable genes to use in heatmap
#' @param quantile_smoothing whether to set value in the heatmap higher than the 99th percentile or lower \
#' than the 1st percentile to the 99th or 1st, respectively
#' @param write_heatmap_txt whether to print a tsv of the heatmap info
#'
#' @return
plotHeatmap <- function (vsd.obj, 
                         with_labels=F,
                         TOPN=1000, 
                         quantile_smoothing=T, 
                         write_heatmap_txt=F,
                         targets=NULL) {
  brewer_palette <- "RdBu"
  ramp <- colorRampPalette(brewer.pal(11, brewer_palette))# Ramp the color in order to get the scale.
  mr <- ramp(256)[256:1]
  stabilized <- assay(vsd.obj)
  rv <- rowVars(stabilized)# calculate the variances by row(gene) to find out which genes are the most variable across the samples.
  select <- order(rv, decreasing=T)[seq_len(TOPN)]# order matrix from top to bottom by the rowVars vairable, take the top 1000 most variable
  M <- stabilized[select,]
  M <- M - rowMeans(M, na.rm=T)# subtract out the means from each row, leaving the variances for each gene.

  allquantile <- quantile(as.matrix(M), 0:100/100)
  quantile_99th <- allquantile[grep("^99%$", names(allquantile))]
  quantile_1st <- allquantile[grep("^1%$", names(allquantile))]
  M[M > quantile_99th] <- quantile_99th
  M[M < quantile_1st] <- quantile_1st

  # annotate the matrix and remove miscellaneous columns
  annot_M <- left_join(cbind(rownames(M),as_tibble(M)), annotation, by = c("rownames(M)" = "Gene.name"))
  #added the next 7 lines of code to avoid duplicate rownames errors
  colnames(annot_M)[1] <- "Gene.name"
  annot_M <- annot_M[!duplicated(annot_M$Gene.name),]
  rownames(annot_M) <- annot_M$Gene.name
  annot_M <- annot_M[, !(colnames(annot_M) %in% c("rownames(M)", "Gene.name", "Gene.description", "Gene.type", "Gene.stable.ID"))]
  # # re-color each of the conditions / treatments
  # ann_colors <-  brewer.pal(length(unique(vsd.obj$condition)), "Paired")
  # names(ann_colors) <- unique(vsd.obj$condition)
  
  if (is.null(targets)) {
    targets <- as.data.frame(vsd.obj@colData)
    targets <- targets[,grep("replaceable|sizeFactor|replicate|sample", colnames(targets), invert = T)]
    rownames(targets) <- rownames(vsd.obj@colData)
    targets <- as.data.frame(targets)
  }

  # call the pheatmap function, providing the raw data (annot_M), annotation info (targets), and color info (ann_colors)
  if (with_labels) {
    my_heatmap <- pheatmap(annot_M, color = mr, annotation_col = targets,
                           fontsize_col = 12, fontsize_row = 6, fontsize = 5 )
  } else {
    my_heatmap <- pheatmap(annot_M, color = mr, annotation_col = targets, show_rownames = F,
                           fontsize_col = 12, cellheight = 1)
  }
  if (write_heatmap_txt) {
    write.table(annot_M, paste0("QCheatmap_",TOPN,"genes.txt"), quote = F, sep = "\t", col.names = NA)
  }
  print(my_heatmap)
}

apply_annotation <- function (res, 
                              my_annotation, 
                              cpms, 
                              joinCPMcol1,
                              joinCPMcol2) {
  res1 <- merge(as.data.frame(res), annotation, 
               by.x="row.names",
               by.y="Gene.name")
  res1 <- res1[,c("Row.names", "Gene.stable.ID", "Gene.description", "Gene.type", 
                "baseMean", "log2FoldChange", "pvalue", "padj")]
  res2 <- merge(res1, cpms, 
                by.x=joinCPMcol1,
                by.y=joinCPMcol2)
  colnames(res2)[grep("Row.names", colnames(res2))] <- "Gene.name"
  return(res2)
}


generate_DEseq <- function (my_data, groups) {
  data_subset1 <- my_data[,grep(str_c("^", groups[1]), colnames(my_data))]
  data_subset2 <- my_data[,grep(str_c("^", groups[2]), colnames(my_data))]
  data_subset <- cbind(data_subset1, data_subset2)
  condition <- c(rep(groups[1],ncol(data_subset1)), rep(groups[2],ncol(data_subset2)))
  targets <- as.data.frame(condition)
  rownames(targets) <- colnames(data_subset)

  print(targets)
  dds <- DESeqDataSetFromMatrix(countData = data_subset,
                                colData = targets,
                                design = ~ condition)
  dds <- DESeq(dds)
  return(dds)
}

#' run_DE **untested**
#' 
#' needs to be fixed so that only samples from the current comparison are 
#' ... outputted to tables
#' 
#' 
#' @description Prints tables of DEGs
#' @param dds a dds object where DESeq2() was already run
#' @param comp a vector of c('condition', 'condition1', 'condition2') in resultsNames(dds)
#' @param annot path to a gene annotation file that will get merged with results
#' @param padj_cutoff float, a cutoff for adjusted p.value
#' @param logfc_cutoff float, a cutoff for log2fc
#' @param nc_cutoff float, a cutoff for mean normalized counts
#'
#' @return NULL
run_DE <- function(dds, 
                   comp,
                   annot="~/Data/Genomes/Human/biomart_dls/human_ensID_genename_240627.txt",
                   padj_cutoff = 0.05,
                   logfc_cutoff = 1,
                   nc_cutoff = 10) {
  if (length(comp) != 3) {
    stop("comp must be length 3")
  } else {
    res <- results(dds, name = paste0(comp[1], "_",
                                        comp[2], "_vs_",
                                        comp[3])) %>%
      data.frame(.) %>% 
      rownames_to_column('Gene stable ID') %>% 
      as_tibble()
  }
   
  # add the annotated gene names and descriptions 
  transl.tab <- read_delim(annot,delim='\t')%>% 
    dplyr::arrange("Gene stable ID")  
  
  # create the normalized counts table
  nc <- as.data.frame(DESeq2::counts(dds, normalized=T))
  nc$mean_norm_counts <- rowMeans(nc) 
  nc <- nc %>% 
    rownames_to_column('Gene stable ID')
  
  # join everything 
  res.gn <- transl.tab %>%
    left_join(res, by='Gene stable ID') %>%
    filter(!is.na(log2FoldChange)) %>% 
    select(c("Gene stable ID", 
             "Gene name", 
             "Gene description", 
             'Gene type',
             'log2FoldChange',
             'pvalue',
             'padj')) %>% 
    left_join(nc, by='Gene stable ID')
  
  # extract DE genes
  de_genes_fdr <- res.gn %>%
    filter(padj < padj_cutoff) %>%
  de_genes_log2f <- de_genes_fdr %>%
    filter(log2FoldChange > logfc_cutoff)
  de_genes_nc <- de_genes_log2f %>%
    filter(mean_norm_counts > nc_cutoff)
  rnk_all <- res.gn %>%
    filter(`Gene type` == 'protein_coding') %>%
    select(c("Gene name", "log2FoldChange")) %>%
    arrange(desc(log2FoldChange))
  rnk_sig <- de_genes_fdr %>%
    filter(`Gene type` == 'protein_coding') %>%
    select(c("Gene name", "log2FoldChange")) %>%
    arrange(desc(log2FoldChange))
  
  # write output to files
  compdir <- paste0("DEG_", comp[1], '_', comp[2], "_vs_", comp[3])
  dir.create(compdir)
  write_delim(de_genes_fdr, file = file.path(compdir, paste0(comp[1], '_', comp[2], 
                                          "_vs_", comp[3], "_fdrgenes.csv")))
  write_delim(de_genes_log2f, file = file.path(compdir, paste0(comp[1], '_', comp[2], 
                                            "_vs_", comp[3], "_log2f.csv")))
  write_delim(de_genes_nc, file = file.path(paste0(comp[1], '_', comp[2], 
                                         "_vs_", comp[3], "_normcounts_cutoff.csv")))
  write_delim(res.gn, file = file.path(paste0(comp[1], '_', comp[2], 
                                    "_vs_", comp[3], "_allgenes.csv")))
  write.table(rnk_all, file = file.path(paste0(comp[1], '_', comp[2], 
                                     "_vs_", comp[3], "_rank_all.rnk")), 
              sep = "\t", row.names = F, quote = F)
  write.table(rnk_sig, file = file.path(paste0(comp[1], '_', comp[2], 
                                     "_vs_", comp[3], "_rank_sig.rnk")), 
              sep = "\t", row.names = F, quote = F)
  
  return(NULL)
}


#' gsea_analysis
#'
#'..... finish documenting this ...
#'
#' @param rnk_path 
#' @param sigs 
#' @param outdir 
#' @param SEED 
#'
#' @return
#' @export
#'
#' @examples
gsea_analysis <- function(allgenes_path,
                          sigs,
                          outdir,
                          SEED=415) {
  analysisdir = file.path(outdir, paste0("GSEA_Results_",format(Sys.time(), "%Y%m%d_%H%M%S")))
  dir.create(analysisdir)
  for (x in c('ALL', 'UP', 'DOWN')) {
    dir.create(file.path(analysisdir, x))
    dir.create(file.path(analysisdir, x, "Figures"))    
  }
  
  rnk <- read_delim(allgenes_path,delim='\t') %>%
    filter(!is.na(`Gene name`))  %>%
    filter(!is.na(log2FoldChange)) %>%
    filter(`Gene type` == 'protein_coding') %>%
    arrange(desc(log2FoldChange)) 
  rnk <- rnk[!duplicated(rnk$`Gene name`),]
  rnk.up <- rnk %>%
    filter(log2FoldChange > 0) %>%
    mutate(rank.up = rank(log2FoldChange, ties.method = "random"))
  rnk.dn <- rnk %>%
    filter(log2FoldChange < 0) %>%
    mutate(rank.dn = -(rank(desc(log2FoldChange), ties.method = 'random')))
  
  rnk.up.v <- setNames(rnk.up$rank.up, 
                       nm=rnk.up$`Gene name`) 
  rnk.dn.v <- setNames(-(rnk.dn$rank.dn), 
                       nm=rnk.dn$`Gene name`) 
  rnk.v <- setNames(rnk$log2FoldChange, 
                    nm=rnk$`Gene name`)
  
  # output rank files 
  write.table(rnk[,c('Gene name', "log2FoldChange")], file.path(analysisdir,'ALL','input_ALL.rnk.txt'),
              quote=F,sep='\t',col.names = NA)
  write.table(rnk.up[,c('Gene name', "rank.up")], file.path(analysisdir,'UP','input_UP.rnk.txt'),
              quote=F,sep='\t',col.names = NA)
  write.table(rnk.dn[,c('Gene name', "rank.dn")], file.path(analysisdir,'DOWN','input_DOWN.rnk.txt'),
              quote=F,sep='\t',col.names = NA)
  
  # Create Tables First ...
  set.seed(SEED)
  
  # BOTH (log2foldchange)
  res.gsea <- fgsea(pathways = sigs, 
                    stats = rnk.v,
                    scoreType = 'std')
  res.gsea$leadingEdge <- sapply(res.gsea$leadingEdge, paste, collapse=', ')
  write.table(res.gsea, file.path(analysisdir, 'ALL', 'gsea_results_ALL.tsv'),
              sep='\t', quote=F, col.names = NA)
  
  # UP (log2foldchange)
  res.gsea <- fgsea(pathways = sigs, 
                    stats = rnk.up.v,
                    scoreType = 'pos')
  res.gsea$leadingEdge <- sapply(res.gsea$leadingEdge, paste, collapse=', ')
  write.table(res.gsea, file.path(analysisdir, 'UP', 'gsea_results_UP.tsv'), 
              sep='\t', quote=F, col.names = NA)
  
  # DOWN (log2foldchange)
  res.gsea <- fgsea(pathways = sigs, 
                    stats = rnk.dn.v,
                    scoreType = 'neg')
  res.gsea$leadingEdge <- sapply(res.gsea$leadingEdge, paste, collapse=', ')
  write.table(res.gsea, file.path(analysisdir, 'DOWN', 'gsea_results_DN.tsv'), 
              sep='\t', quote=F, col.names = NA)
  
  # Create Figures ...
  for (i in seq_along(sigs)) {
    
    # ALL
    pdf(file.path(analysisdir, "ALL", 'Figures', 
                  paste0(names(sigs)[i], "_ALL.pdf")))
    print(plotEnrichment(pathway=sigs[[i]], 
                         stats=rnk.v))
    dev.off()
    
    # UP
    pdf(file.path(analysisdir, "UP", 'Figures', 
                  paste0(names(sigs)[i], "_UP.pdf")))
    print(plotEnrichment(pathway=sigs[[i]], 
                         stats=rnk.up.v))
    dev.off()
    
    # DOWN
    pdf(file.path(analysisdir, "DOWN", 'Figures', 
                  paste0(names(sigs)[i], "_DOWN.pdf")))
    print(plotEnrichment(pathway=sigs[[i]], 
                         stats=rnk.dn.v))
    dev.off()
  }
  return(NULL)
}



# old run_DE function
run_DE_old <- function (dds, 
                    comparisons, 
                    coi) {
  i <- comparisons

  raw_counts <- counts(dds, normalized = F)
  cpms <- as_tibble(rowMeans(cpm(raw_counts)))
  cpms <- cbind(rownames(raw_counts), cpms)
  colnames(cpms) <- c("ensembl", "cpm")

  res <- results(dds, contrast = c(coi, i[1], i[2]))[,-c(3,4)]
  res <- apply_annotation(res, annotation, cpms,
                          joinCPMcol1 = "Row.names",
                          joinCPMcol2 = "row.names")

  # extract DE genes
  de_genes_fdr <- sort_genes_fdr(res, fdrcutoff)
  de_genes_rawp <- sort_genes_rawp(res, rawpcutoff)
  de_genes_log2f <- sort_genes_log2f(res, fdrcutoff, log2cutoff)
  de_genes_cpm <- sort_genes_cpm(res, fdrcutoff, log2cutoff, cpmcutoff)

  # combine normalized counts with entire DE list
  normalized_counts <- round(counts(dds, normalized = TRUE),3)
  pattern <- str_c(i[1], "|", i[2])
  combined_data <- as_tibble(cbind(res, normalized_counts[,grep(pattern, colnames(normalized_counts))] ))
  combined_data <- combined_data[order(combined_data$log2FoldChange, decreasing = T),]

  # make ordered rank file for GSEA, selecting only protein coding genes
  res_prot <- res[which(res$Gene.type == "protein_coding"),]
  resOrdered <- res_prot[order(res_prot$log2FoldChange, decreasing = T),c("Gene.name", "log2FoldChange", "padj")]
  resOrdered <- na.omit(resOrdered)
  resOrdered$Gene.name <- str_to_upper(resOrdered$Gene.name)
  resOrdered <- resOrdered[,-3]

  # write output to files
  write.csv (de_genes_fdr, file = paste0(i[1], "_vs_", i[2], "_fdrgenes.csv"))
  write.csv (de_genes_rawp, file = paste0(i[1], "_vs_", i[2], "_rawpgenes.csv"))
  write.csv (de_genes_log2f, file = paste0(i[1], "_vs_", i[2], "_log2f.csv"))
  write.csv (de_genes_cpm, file = paste0(i[1], "_vs_", i[2], "_cpm_cutoff.csv"))
  write.csv (combined_data, file = paste0(i[1], "_vs_", i[2], "_allgenes.csv"))
  write.table (resOrdered, file = paste0(i[1], "_vs_", i[2], "_rank.rnk"), sep = "\t", row.names = F, quote = F)
  return ( c(nrow(de_genes_rawp), nrow(de_genes_fdr), nrow(de_genes_log2f), nrow(de_genes_cpm)))
}


#' plot_volcano
#' @description Prints a ggplot2 volcano plot for allgenes data
#' @param file The allgenes.csv file generated from deseq2 de results being used to make the figure
#' @param withlabels Whether to print labels on the figure
#' @param labels The concatenated list of labels to label on the figure
#' @param p_opts plot options, a ggplot2 layer which will be added to the plot before printing at the end
#' @return
plot_volcano <- function(file, withlabels=T, labels=NULL, noOthers=F, p_opts=NULL) {
  if (grepl("allgenes.csv",file)) {
    name <- gsub(".*/|_allgenes.csv", "", file)
    # name_num <- strsplit(name, "_vs_")[[1]][1]
    # name_den <- strsplit(name, "_vs_")[[1]][2]
    comp <- read.csv(file, header=T)
  } else if (grepl("*.tsv", file)) {
    name <- gsub(".*/|\\.*", "",file)
    comp <- read_delim(file, delim = "\t")
  }
  
  comp <- comp %>%
    select(c( "Gene name", "log2FoldChange", "padj"))
  comp$log10fdr <- -log10(comp$padj)
  comp$fdrpass <- ifelse(comp$padj < fdrcutoff, "passing FDR", "not passing")
  comp$l2fcdir <- ifelse(comp$log2FoldChange > 0, "l2fc+", "l2fc-")
  comp$fdr_l2fc_color <- paste0(comp$fdrpass, " and ", comp$l2fcdir)
  comp$distfromorigin <- NA
  select <- which(!is.na(comp$log2FoldChange) & !is.na(comp$log10fdr))
  comp$distfromorigin[select] <- sqrt((comp$log10fdr[select])^2 + (comp$log2FoldChange[select])^2 )
  dist_cutoff_high <- tail(head(sort(comp$distfromorigin[!is.na(comp$distfromorigin)], decreasing = T), n=10), n=1)
  if (withlabels && !is.null(labels)) {
    comp$label <- NA
    comp$label[which(comp$`Gene name` %in% labels)] <- comp$Gene.name[which(comp$`Gene name` %in% labels)]
  } else if (withlabels && is.null(labels)) {
    comp$label <- NA
    comp$label[which(comp$distfromorigin >= dist_cutoff_high)] <-
      comp$`Gene name`[which(comp$distfromorigin >= dist_cutoff_high)]
  }
  decolors <- c(rep("#000000", 5), "#ff595e", "#1982c4")
  names(decolors) <- c("NA and l2fc+", "NA and l2fc-",
                       "NA and NA", "not passing and l2fc+",
                       "not passing and l2fc-", "passing FDR and l2fc+",
                       "passing FDR and l2fc-")
  comp <- comp[!is.na(comp$padj),]
  if (noOthers) {
    comp <- comp[which(!is.na(comp$label) | comp$fdrpass == "passing FDR"),]
  }
  p <- ggplot(data=comp, aes(x=log2FoldChange,y=log10fdr,color=fdr_l2fc_color)) +
    geom_hline(yintercept = 0, color="gray") +
    geom_vline(xintercept = 0, color="gray") +
    geom_point() +
    scale_color_manual(values=decolors) +
    labs(title=paste0(name)) + 
    theme(axis.line.x = element_line(color="black"),
          axis.line.y = element_line(color="black"),
          panel.background = element_rect(color="black", fill = "white"))
  if (withlabels) {
    p <- p + geom_text_repel(aes(label=label),size= 2.5, color="black",
                             segment.size = 0.5,point.padding = 0.2,
                             min.segment.length = 0, seed=415)
  }
  if (!is.null(p_opts)) {
    p <- p + p_opts
  }
  print(p)
  return(comp)
}


#' plot_MM
#' @description Plots an MM plot from two allgenes.txt DE files
#' @param file1 log2fc for X axis in MM plot
#' @param file2 log2fc for Y axis in MM plot
#' @param nonDEremoval TRUE if we want to remove the genes that don't pass FDR threshold for either comparison
#' @param max_axis a maximum value for x and y axes (for manually setting into a square)
#' @param labels a concatenated list of genes used for labeling the figure
#' @return An MMplot using ggplot scatter plot
plot_MM <- function(file1, file2, nonDEremoval=T, max_axis=NULL, labels=NULL) {
  name1 <- gsub(".*/|_allgenes.csv", "", file1)
  name2 <- gsub(".*/|_allgenes.csv", "", file2)
  comp1 <- read.csv(file1, header=T)
  comp1 <- comp1[,grep("log2FoldChange|padj|Gene\\.name", colnames(comp1))]
  comp2 <- read.csv(file2, header=T)
  comp2 <- comp2[,grep("log2FoldChange|padj|Gene\\.name", colnames(comp2))]
  lowercutoff <- floor(min(min(comp1$log2FoldChange), min(comp2$log2FoldChange)))
  uppercutoff <- ceiling(max(max(comp1$log2FoldChange), max(comp2$log2FoldChange)))
  padj.col1 <- paste0("padj_", name1)
  padj.col2 <- paste0("padj_", name2)
  l2fc.col1 <- paste0("log2FoldChange_", name1)
  l2fc.col2 <- paste0("log2FoldChange_", name2)
  comps <- merge(comp1, comp2, by = "Gene.name", suffixes=c(paste0("_",name1),paste0("_", name2)))
  comps <- comps[which(!is.na(comps[,l2fc.col1]) & !is.na(comps[,l2fc.col2])),]
  # comps$dist <- 1
  # comps$dist[which(!is.na(comps[,padj.col1]) & !is.na(comps[,padj.col2]))] <-
  #   apply(comps[which(!is.na(comps[,padj.col1]) & !is.na(comps[,padj.col2])),
  #               grep(paste0(padj.col1, "|", padj.col2), colnames(comps))], 1, max)
  comps$passing_fdr <- "not_DE"
  comps$passing_fdr[which(comps[,padj.col1] < fdrcutoff)] <- paste0(name1, "_DE")
  comps$passing_fdr[which(comps[,padj.col2] < fdrcutoff)] <- paste0(name2, "_DE")
  comps$passing_fdr[which(comps[,padj.col1] < fdrcutoff & comps[,padj.col2] < fdrcutoff)] <- "both_DE"
  decolors <- c("#000000", "#ff595e", "#1982c4", "#8ac926")
  names(decolors) <- c("not_DE", paste0(name1, "_DE"), paste0(name2, "_DE"), "both_DE")
  comps$labels <- NA
  if (is.null(labels)) {
    if (length(comps$Gene.name[which(comps$passing_FDR=="both_DE")]) > 8) {
      comps_filtered <- comps[which(comps$passing_FDR == "both_DE"),]
      genes2use <- head(comps_filtered[order(comps_filtered[,l2fc.col1]),'Gene.name'], n=8)
      comps$labels[which(comps$Gene.name %in% genes2use)] <- comps$Gene.name[which(comps$Gene.name %in% genes2use)]
    } else {
      comps$labels[which(comps$passing_fdr == "both_DE")] <- comps$Gene.name[which(comps$passing_fdr == "both_DE")]
    }
  } else {
    comps$labels[which(comps$Gene.name %in% labels)] <- comps$Gene.name[which(comps$Gene.name %in% labels)]
  }

  if (nonDEremoval) {
    comps <- comps[which(comps$passing_fdr != "not_DE"),]
  }
  if (!is.null(max_axis)) {
    uppercutoff = max_axis
    lowercutoff = -max_axis
  }
  p <- ggplot(data=comps, aes_string(x=paste0("log2FoldChange_", name1),
                                y=paste0("log2FoldChange_", name2),
                                color = "passing_fdr",
                                # size="dist",
                                alpha = 0.8)) +
    geom_point() +
    scale_size(breaks=c(0.05),
      range=c(0.3,3), trans = "reverse") +
    scale_color_manual(values=decolors) +
    geom_text_repel(aes(label=labels),size= 2.5, color="black",
                    segment.size = 0.5,point.padding = 0.2,
                    min.segment.length = 0, seed=415, alpha=1) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_abline(slope=1, intercept=0) +
    # geom_abline(slope=-1, intercept=0) +
    scale_x_continuous(limits = c(lowercutoff, uppercutoff)) +
    scale_y_continuous(limits = c(lowercutoff, uppercutoff)) +
    theme(axis.line.x = element_line(color="black"),
          axis.line.y = element_line(color="black"),
          panel.background = element_rect(color="black", fill = "white"),
          axis.text = element_text(size = 12))
  print(p)
  return(comps)
}


sort_genes_fdr = function(df, fdrcutoff) {
  df = df[which(df$padj < fdrcutoff),]
  #df = df[order(abs(df$log2FoldChange), decreasing=TRUE),]
  df = df[order(df$log2FoldChange, decreasing=TRUE ),]
  return(df)
}

sort_genes_rawp = function(df, rawpcutoff) {
  df = df[which(df$pvalue < rawpcutoff),]
  #df = df[order(abs(df$log2FoldChange), decreasing=TRUE),]
  df = df[order(df$log2FoldChange, decreasing=TRUE ),]
  return(df)
}

sort_genes_log2f = function(df, fdrcutoff, log2cutoff) {
  df = df[which(df$padj < fdrcutoff),]
  df = df[which(abs(df$log2FoldChange) > log2cutoff),]
  #df = df[order(abs(df$log2FoldChange), decreasing=TRUE),]
  df = df[order(df$log2FoldChange, decreasing=TRUE ),]
  return(df)
}

sort_genes_cpm = function(df, fdrcutoff, log2cutoff, cpmcutoff) {
  df = df[which(df$padj < fdrcutoff),]
  df = df[which(df$cpm > cpmcutoff),]
  df = df[which(abs(df$log2FoldChange) > log2cutoff),]
  #df = df[order(abs(df$log2FoldChange), decreasing=TRUE),]
  df = df[order(df$log2FoldChange, decreasing=TRUE ),]
  return(df)
}

#' @title zip
#' @description Equivalent to Python's zip() definition
zip_paste <- function(x,y,...,sep=" ") {
  z = data.frame(x,y,...)
  return(apply(z, 1, paste, collapse=sep))
}



#' boxplot_panel_personalis
#' @description same as boxplot_panel but we're assuming the gene name is already in the counts matrix
#'
#' @param counts The gene expression matrix used for values in the boxplot
#' @param gene_input The genes to be made a panel of (something should be standard)
#' @param coi The condition of interest
#' @param targets The condition model including a column that matches coi
#' @param cond_order The factor order of the coi for sorting
#' @param output_folder The folder used for saving all boxplot pdfs
#' @return NULL (just prints boxplots to folder)
boxplot_panel_personalis <- function(counts,
                          gene_input,
                          coi,
                          gene_display=F,
                          targets=NULL,
                          cond_order=NULL,
                          output_folder=NULL) {
  if (is.null(targets)) {
    print("include targets..")
    return(NULL)
  }
  if (!is.null(output_folder)) {
    output_folder <- file.path(output_folder, "boxplot_panel")
    dir.create(output_folder)
  }
  if (is.null(cond_order)) {
    cond_order <- sort(unique(targets[,coi]))
  }
  if (gene_display) {
    if (length(gene_input) == 1) {
      stop("need more than 1 gene for gene_display=T")
    }
  }
  gene_input.notincounts <- gene_input[which(!(gene_input %in% rownames(counts)))]
  message(paste0(paste(gene_input.notincounts, collapse=", "), " aren't in counts, removing..."))
  gene_input <- gene_input[which(gene_input %in% rownames(counts))]
  
  if (length(gene_input) == 1) {
    print("only 1 gene...")
    if (!all(rownames(targets) == colnames(counts))) {
      print("rownames of targets and colnames of counts aren't the same...")
      return(NULL)
    }
    if (length(unique(targets[,coi])) == 1 || length(unique(targets[,coi])) > 2) {
      result.p = NA
    } else if (length(unique(targets[,coi])) == 2) {
      result <- wilcox.test(formula = counts[gene_input,] ~ targets[,coi])
      result.p <- result$p.value
    } 
    cur_counts=counts[gene_input,]
    cur_targets <- cbind(targets, counts=t(cur_counts))
    cur_targets[,coi] <- factor(cur_targets[,coi], levels=cond_order) # set the order for neatness
    gene_name <- gene_input
    p <- ggplot(data=cur_targets,aes_string(x=coi,y="counts",
                                            group=coi)) +
      geom_boxplot() +
      geom_jitter(width=0.1,height=0) +
      labs(title=paste0("Expression of ",gene_name, " across ", coi),
           subtitle = paste0("Mann-whitney p = ", round(result.p,4)),
           x=coi,y=gene_name)
    
    pdf(paste0(output_folder, gene_name, "_by_", coi,".pdf"))
    print(p)
    dev.off()
    return(NULL)
  }
  
  if (gene_display) {
    cur_counts=counts[gene_input,]
    cur_targets <- cbind(targets, t(cur_counts))
    cur_targets.melt <- reshape2::melt(cur_targets,
                                       measure.vars = gene_input)
    for (i in seq_along(cond_order)) {
      if (!all(rownames(targets) == colnames(counts))) {
        print("rownames of targets and colnames of counts aren't the same...")
        return(NULL)
      }
      cur <- cur_targets.melt[which(cur_targets.melt[[coi]] == cond_order[i]),]
      p <- ggplot(data=cur,aes_string(x="variable",y="value",
                                      group="variable")) +
        geom_boxplot() +
        geom_jitter(width=0.1,height=0) +
        labs(title=paste0("Expression of genes across indication: ", cond_order[i]),
             x="Genes",y="TPM") +
        theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))
      
      rng <- range(cur$value)
      if ((rng[2] - rng[1]) > 1e3) {
        p <- p + scale_y_log10()
      }
      
      pdf(file.path(output_folder, paste0(cond_order[i], ".pdf")))
      print(p)
      dev.off()
    }
  } else {
    for (i in 1:length(gene_input)) {
      # append counts[genes[i]] to targets
      # (optional) append DE stats info to separate dataframe
      # plot & save
      
      if (!all(rownames(targets) == colnames(counts))) {
        print("rownames of targets and colnames of counts aren't the same...")
        return(NULL)
      }
      if (length(unique(targets[,coi])) == 1 || length(unique(targets[,coi])) > 2) {
        result.p = NA
      } else if (length(unique(targets[,coi])) == 2) {
        result <- wilcox.test(formula = counts[gene_input,] ~ targets[,coi])
        result.p <- result$p.value
      } 
      cur_counts=counts[gene_input[i],]
      cur_targets <- cbind(targets, counts=cur_counts)
      cur_targets[,coi] <- factor(cur_targets[,coi], levels=cond_order) # set the order for neatness
      p <- ggplot(data=cur_targets,aes_string(x=coi,y="counts",
                                              group=coi)) +
        geom_boxplot() +
        geom_jitter(width=0.1,height=0) +
        labs(title=paste0("Expression of ",gene_input[i], " across ", coi),
             subtitle = paste0("Mann-whitney p = ", round(result.p,4)),
             x=coi,y=gene_input[i])
      
      pdf(file.path(output_folder, paste0(gene_input[i], "_by_", coi,".pdf")))
      print(p)
      dev.off()
    }
  }
  return(NULL)
}

#' boxplot_panel
#'
#' @param counts The gene expression matrix used for values in the boxplot
#' @param gene_input The genes to be made a panel of (something should be standard)
#' @param coi The condition of interest
#' @param annot The annotation matrix which connects gene names with ensembl ids existing in the counts table
#' @param targets The condition model including a column that matches coi
#' @param cond_order The factor order of the coi for sorting
#' @param output_folder The folder used for saving all boxplot pdfs
#' @return NULL (just prints boxplots to folder)
boxplot_panel <- function(counts,
                          gene_input,
                          coi,
                          annot=annotation,
                          targets=NULL,
                          cond_order=NULL,
                          output_folder=NULL) {
  if (is.null(targets)) {
    print("include targets..")
    return(NULL)
  }
  if (!is.null(output_folder)) {
    output_folder <- "boxplot_panel/"
    dir.create(output_folder)
  }
  check_genes <- gene_input %in% annot$Gene.name
  names(check_genes) <- gene_input
  print(check_genes)
  gene_annot <- annot[annot$Gene.name %in% gene_input,]
  gene_annot <- gene_annot[!duplicated(gene_annot$Gene.name),]
  gene_list <- gene_annot$Gene.stable.ID
  names(gene_list) <- gene_annot$Gene.name
  if (is.null(cond_order)) {
    cond_order <- sort(unique(targets[,coi]))
  }
  if (length(gene_list) == 1) {
    print("only 1 gene...")
    if (!all(rownames(targets) == colnames(counts))) {
      print("rownames of targets and colnames of counts aren't the same...")
      return(NULL)
    }
    cur_counts=counts[gene_list,]
    cur_targets <- cbind(targets, counts=cur_counts)
    cur_targets[,coi] <- factor(cur_targets[,coi], levels=cond_order) # set the order for neatness
    gene_name <- gene_input
    p <- ggplot(data=cur_targets,aes_string(x=coi,y="counts")) +
      geom_boxplot() +
      labs(title=paste0("Expression of ",gene_name, " across ", coi),
           x=coi,y=gene_name)

    pdf(paste0(output_folder, gene_name, "_by_", coi,".pdf"))
    print(p)
    dev.off()
    return(NULL)
  }

  for (i in 1:length(gene_list)) {
    # append counts[genes[i]] to targets
    # (optional) append DE stats info to separate dataframe
    # plot & save

    if (!all(rownames(targets) == colnames(counts))) {
      print("rownames of targets and colnames of counts aren't the same...")
      return(NULL)
    }
    cur_counts=counts[gene_list[i],]
    cur_targets <- cbind(targets, counts=cur_counts)
    cur_targets[,coi] <- factor(cur_targets[,coi], levels=cond_order) # set the order for neatness
    gene_name <- names(gene_list)[i]
    p <- ggplot(data=cur_targets,aes_string(x=coi,y="counts")) +
      geom_boxplot() +
      labs(title=paste0("Expression of ",gene_name, " across ", coi),
           x=coi,y=gene_name)

    pdf(paste0(output_folder, gene_name, "_by_", coi,".pdf"))
    print(p)
    dev.off()
  }
  return(NULL)
}

#' mediannormsig
#' @description Creates a gene signature to add to a model matrix (like targets) as a covariate
#' @param gene_input The gene names in annotation to use for signature
#' @param annot The annotation object to use with Gene.name as a column
#' @param cpm_temp The gene expression data, usually cpm is used
#' @return A vector of covariates, with length = # input samples
mediannormsig <- function(gene_input,annot=annotation,cpm_temp=cpm) {
  gene_annot <- annot[annot$Gene.name %in% gene_input,]
  gene_annot <- gene_annot[!duplicated(gene_annot$Gene.name),]
  gene_list <- gene_annot$Gene.stable.ID
  names(gene_list) <- gene_annot$Gene.name
  cpm_temp <- as.data.frame(cpm_temp)
  # should add something here that states which genes didn't give complete cases
  cpm_temp <- cpm_temp[complete.cases(cpm_temp),]
  mydf <- t(cpm_temp[gene_list,])
  cmed <- apply(mydf, 2, median, na.rm = TRUE)
  cmed <- cmed / mean(cmed)
  data2 <- t(t(mydf)/cmed)
  s <- as.data.frame(matrix(nrow=nrow(data2), ncol=ncol(data2)))
  rownames(s) <- rownames(mydf)
  colnames(s) <- colnames(mydf)
  for (i in c(1:ncol(s))) {
    f <- (data2[,i] - median(data2[,i], na.rm=T)) / sd(data2[,i], na.rm=T)
    s[,i] <- f
  }
  mysig <- rowMeans(s, na.rm=T)
  return(mysig)
}


#' mediannormsig_simple
#' @description Creates a gene signature to add to a model matrix (like targets) as a covariate
#' @param gene_input The gene names in annotation to use for signature
#' @param cpm_temp The gene expression data that gene_input will look up in 
#' @return A vector of covariates, with length = # input samples
mediannormsig_simple <- function(gene_input, cpm_temp) {
  in_cpm_temp <- gene_input[gene_input %in% rownames(cpm_temp)]
  notin_cpm_temp <- gene_input[which(!(gene_input %in% rownames(cpm_temp)))]
  if (length(in_cpm_temp) == 0) {
    stop("no genes from gene_input are present in rownames of cpm_temp")
  }
  if (length(in_cpm_temp) < length(gene_input)) {
    message(paste0("these were not found in rownames of cpm_temp: \n", 
                   paste(notin_cpm_temp, collapse = ", "),
            " \ncontinuing with leftovers: \n",
            paste(in_cpm_temp, collapse = ", ")))
    gene_input <- in_cpm_temp
  }
  mydf <- t(cpm_temp[gene_input,])
  cmed <- apply(mydf, 2, median, na.rm = TRUE)
  cmed <- cmed / mean(cmed)
  data2 <- t(t(mydf)/cmed)
  s <- as.data.frame(matrix(nrow=nrow(data2), ncol=ncol(data2)))
  rownames(s) <- rownames(mydf)
  colnames(s) <- colnames(mydf)
  for (i in c(1:ncol(s))) {
    f <- (data2[,i] - median(data2[,i], na.rm=T)) / sd(data2[,i], na.rm=T)
    s[,i] <- f
  }
  mysig <- rowMeans(s, na.rm=T)
  return(mysig)
}


#' norm_ratio_sig
#' @description Creates a gene signature to add to a model matrix (like targets) as a covariate
#' @param gene_input The gene names in annotation to use for signature
#' @param cpm_temp The gene expression data that gene_input will look up in 
#' @return A vector of covariates, with length = # input samples
norm_ratio_sig <- function(gene_list1, gene_list2, cpm_temp) {
  if (!all(gene_list1 %in% rownames(cpm_temp))) {
    stop("not all of gene_list1 are present in rownames of cpm_temp")
  }
  if (!all(gene_list2 %in% rownames(cpm_temp))) {
    stop("not all of gene_list2 are present in rownames of cpm_temp")
  }
  mydf <- t(cpm_temp[c(gene_list1, gene_list2),])
  cmed <- apply(mydf, 2, median, na.rm = TRUE)
  cmed <- cmed / mean(cmed)
  data2 <- t(t(mydf)/cmed)
  s <- as.data.frame(matrix(nrow=nrow(data2), ncol=ncol(data2)))
  rownames(s) <- rownames(mydf)
  colnames(s) <- colnames(mydf)
  for (i in c(1:ncol(s))) {
    f <- (data2[,i] - median(data2[,i], na.rm=T)) / sd(data2[,i], na.rm=T)
    s[,i] <- f
  }
  mysig <- rowMeans(s, na.rm=T)
  return(mysig)
}

# GSEA ----
processed_ranks <- function(rank) {
  if (class(rank$log2FoldChange) != "numeric") { #if the log2FoldChange is no numeric- makes it numeric
    rank$log2FoldChange <- as.numeric(rank$log2FoldChange)

  }
  if(dim(rank[duplicated(rank$Gene.name),]) [1] > 0) { #checks if there are duplicates
    dupe = rank[,"Gene.name"]
    dupes <- rank[duplicated(dupe) | duplicated(dupe, fromLast=TRUE),] #finds the gene name
    print("Data has duplicates, average of duplicate names was taken.")
    show(dupes) #shows the duplicates
    rank <- aggregate(.~Gene.name, FUN = mean, data = rank) #takes the mean of the duplicated rows
  }

  if (any(is.na(rank)) == T) { #if there are any NA's in the ranks this will remove them
    rank <- rank[complete.cases(rank), ]
  }
  rank <- tibble::deframe(rank) #makes ranks into a vector
  rank
}

# for inputting plotenrichment easier
plot.enrichment <- function (geneset, pathway, ranked_list) {
  plotEnrichment(geneset[[pathway]], ranked_list)+labs (title = pathway)
}

#ranks <- tibble::deframe(ranks)
fgsea_object <- function(mypathway, rank, write_table = F, minSize=15) { #creates the fgsea object
  fgsea_obj <- fgsea(pathways = mypathway,
                     stats = rank,
                     minSize = minSize,
                     maxSize = 500,
                     nperm= 1000)

}

signif_barplot <- function(fgseas_obj, mytitle,p_val=0.05) {
  fgseas_obj$pathway <- gsub("^[^_]*_*", "\\1", fgseas_obj$pathway) #removes up to and including the first underscore of each pathway
  fgseaRes_tidy <- fgseas_obj %>% as_tibble() %>% arrange(desc(NES))
  ggplot(fgseaRes_tidy, aes(reorder(pathway, NES), NES))+
    geom_col(aes(fill = padj < p_val)) + coord_flip() +
    labs(x = "Pathway", y = "Normalized Enrichment Score", title = mytitle, fill = "") +
    theme(axis.text.y = element_text(size = 5)) +
    scale_fill_discrete(name = paste0("padj <", p_val))
}

fgsea_barplot_hallmark <- function(rankfile_loc,
                                   hallmark.pathway=hallmark.pathway,
                                   GO.pathway=GO.pathway,
                                   fromfile=T
) {
  if (!fromfile) {
    ranks <- rankfile_loc
  } else {
    ranks <- read.table(rankfile_loc, header = T,
                        colClasses = c("character", "numeric"), stringsAsFactors = F)
  }
  rank <- processed_ranks(ranks)
  hallmark_object <- fgsea_object(hallmark.pathway, ranks, write_table = T)
  signif_barplot(hallmark_object, "hallmark genes")
}

# it just does it all.
entire_fgsea_analysis <- function(folder, rawpcutoff=0.05, pathways = "both",
                                  pattern=NULL, alt_pathways=NULL,output=NULL) {
  if (is.null(output)) {
    output <- paste0(folder, "/fgsea_analysis")
    dir.create(paste0(folder, "/fgsea_analysis"))
  } else {
    dir.create(output)
  }
  if (pathways == "hallmark") {
    hallmark.pathway <- gmtPathways("~/Data/h.all.v7.0.symbols.gmt.txt")
    if (!is.null(pattern)) {
      cur.pathways <- hallmark.pathway[grep(pattern, names(hallmark.pathway))]
      message("pattern: ", pattern, " found ",
              length(names(cur.pathways)), " pathways")
      if (length(cur.pathways) == 0) {
        stop(paste0("pattern: ", pattern , " found no matches"))
      }
    } else {
      cur.pathways <- hallmark.pathway
    }
  } else if (pathways == "go") {
    GO.pathway <- gmtPathways("~/Data/c5.all.v7.0.symbols.gmt.txt")
    if (!is.null(pattern)) {
      cur.pathways <- GO.pathway[grep(pattern, names(GO.pathway))]
      message("pattern: ", pattern, " found ",
              length(names(cur.pathways)), " pathways")
      if (length(cur.pathways) == 0) {
        stop(paste0("pattern: ", pattern , " found no matches"))
      }
    } else {
      cur.pathways <- GO.pathway
    }
  } else {
    if (!is.null(alt_pathways)) {
      if (file.exists(alt_pathways)) {
        message(paste0("not hallmark or go, using: ", alt_pathways, " as pathway file"))
        cur.pathways <- gmtPathways(alt_pathways)
        if (!is.null(pattern)) {
          cur.pathways <- cur.pathways[grep(pattern, names(GO.pathway))]
          message("pattern: ", pattern, " found ",
                  length(names(cur.pathways)), " pathways")
        }
      } else {
        stop("pathway file not found...")
      }
    }
  }
  if (length(cur.pathways) > 100) {
    input=readline(paste0("Processing over 100 pathways (", length(cur.pathways) ,
                          "), are you sure you want to continue?(Y/N):"))
    if (tolower(input) != "y") {
      stop("stopping...")
    }
  }
  allfiles=list.files(folder)
  rankfiles=allfiles[grep("*\\.rnk", allfiles)]
  if (length(rankfiles) == 0) {
    print("can't find any rankfiles in this dir, please specify another.")
    return(NULL)
  }
  for (i in seq_along(rankfiles)) {
    compname <- gsub("_rank\\.rnk", "", rankfiles[i])
    ranks <- read.table(paste0(folder, "/", rankfiles[i]), header=T,
                        colClasses = c("character", "numeric"), stringsAsFactors = F)
    ranks_current <- processed_ranks(ranks)
    ranks_current_loaded <- fgsea_object(cur.pathways, ranks_current, write_table = T)
    ranks_current_loaded$leadingEdge <- do.call(rbind, lapply(as.list(ranks_current_loaded$leadingEdge), paste0, collapse=", "))
    ranks_current_loaded <- ranks_current_loaded %>% arrange (desc(NES))
    write.csv(ranks_current_loaded, file = paste0(output, "/", compname,"_ALLpathways_gsea_table.csv"))
    passing <- ranks_current_loaded[which(ranks_current_loaded$pval < rawpcutoff),]
    if (nrow(passing) == 0) {
      message(paste0("no significant pathways found for ", rankfiles[i]))
      next #don't make plots for ranks that aren't even significant
    }
    dir.create(paste0(output, "/enrichment"))
    for (j in seq_along(passing$pathway)) {
      message(paste0("printing enrichment plot: ", passing$pathway[j], " for: ", rankfiles[i]))
      dir.create(paste0(output, "/enrichment/", compname))
      pdf(paste0(output, "/enrichment/", compname, "/", passing$pathway[j], ".pdf"),width=5,height=4)
      print(plot.enrichment(cur.pathways, passing$pathway[j] , ranks_current))
      dev.off()
    }
    message(paste0("printing barplot for: ", rankfiles[i]))
    # print the barplot
    pdf(paste0(output, "/", compname, "signif_barplot.pdf"))
    print(signif_barplot(ranks_current_loaded, "pathways", p_val <- rawpcutoff))
    dev.off()

    passing$leadingEdge <- do.call(rbind, lapply(as.list(passing$leadingEdge), paste0, collapse=", "))
    passing <- passing %>% arrange (desc(NES))
    write.csv(passing, file = paste0(output, "/", compname,"_pathway_gsea_table.csv"))
  }
}


find_in_annot <- function(g, annotation) {

}

#' closest_interaction
#' @description naive implementation of a shortest path algorithm for finding \
#' shortest path between two proteins in the stringdb graph
#' @note TMW you write BFS and it's implemented in igraph already...
#' @param gene1 root gene in protein search
#' @param gene2 end gene in protein search
#' @param db the string_db name (named string_db on the vignette)
#' @param debug whether we want to print messages
#'
#' @return NULL if not found, a concatenated string of protein ids if found
closest_interaction <- function(gene1, gene2, db, debug=F) {
  requires(STRINGdb)
  # create an array of all nodes here with labels of "explored" or "unexplored"
  #  so we can know whether to include or not include in the queue
  if (is.null(db)) {
    string_db <- STRINGdb$new( version="11",
                               species=9606,
                               score_threshold=200,
                               input_directory="")
  } else string_db <- db
  all_proteins <- string_db$get_proteins()
  all_proteins.visit <- new.env()
  message("creating ht of visited nodes...")
  for (i in 1:nrow(all_proteins)) {
    assign(all_proteins$protein_external_id[i], new.env(), all_proteins.visit) # initialize visited graph
    all_proteins.visit[[all_proteins$protein_external_id[i]]]$preferred_name <-
      all_proteins$preferred_name[i]
    all_proteins.visit[[all_proteins$protein_external_id[i]]]$visited <- F
  }
  message("getting gene1's protein...")
  gene1.p <- string_db$mp(gene1)
  message("getting gene2's protein...")
  gene2.p <- string_db$mp(gene2)
  gene1.p.env <- rlang::env(protein=gene1.p,path=gene1.p)
  nq <- list(gene1.p.env)
  all_proteins.visit[[gene1.p]]$visited <- T # note that the first protein has been visited
  curpath <- c()
  message("starting graph search...")
  while (length(nq) != 0) {
    n <- nq[[1]]$protein #use as current neighbor to search for gene2
    curpath <- nq[[1]]$path
    if (debug) {
      # message("length of nq: ", length(nq))
      # print(paste0("current protein: ", n, ", ", all_proteins.visit[[n]]$preferred_name))
      # message(paste0("current path:  "))
      # message(curpath)
    }
    if (n == gene2.p) {
      message("found path: ")
      return(paste(curpath, collapse=", "))
    } else {
      for (nn in string_db$get_neighbors(n)) {
        if (!all_proteins.visit[[nn]]$visited) {
          nq <- c(nq, rlang::env(protein=nn, path=c(curpath, nn)))
          all_proteins.visit[[nn]]$visited <- T # note that this protein has been visited
        }
      }
      nq <- nq[-1] #remove the front of queue
    }
  }
  print("no path...")
  v.p <- 0
  for (p in all_proteins$protein_external_id) {
    if (all_proteins.visit[[p]]$visited) {
      v.p <- v.p + 1
    }
  }
  message(paste0("visited ",
                 length(v.p),
                 " proteins out of ",
                 length(all_proteins$protein_external_id),
                 "..."))
  return(NULL)
}


# 
output_de_tables <- function(res, annotation,cpms, outdir,name) {
  
  res <- apply_annotation(res, annotation, cpms,
                          joinCPMcol1 = "Row.names",
                          joinCPMcol2 = "row.names")
  
  de_genes_fdr <- sort_genes_fdr(res, fdrcutoff)
  de_genes_rawp <- sort_genes_rawp(res, rawpcutoff)
  de_genes_log2f <- sort_genes_log2f(res, fdrcutoff, log2cutoff)
  de_genes_cpm <- sort_genes_cpm(res, fdrcutoff, log2cutoff, cpmcutoff)
  
  # make ordered rank file for GSEA, selecting only protein coding genes
  res_prot <- res[which(res$Gene.type == "protein_coding"),]
  resOrdered <- res_prot[order(res_prot$log2FoldChange, decreasing = T),c("Gene.name", "log2FoldChange", "padj")]
  resOrdered <- na.omit(resOrdered)
  resOrdered$Gene.name <- str_to_upper(resOrdered$Gene.name)
  resOrdered <- resOrdered[,-3]
  
  # write output to files
  write.csv (de_genes_fdr, file = file.path(outdir, paste0(name, "_fdrgenes.csv")))
  write.csv (de_genes_rawp, file = file.path(outdir, paste0(name, "_rawpgenes.csv")))
  write.csv (de_genes_log2f, file = file.path(outdir, paste0(name, "_log2f.csv")))
  write.csv (de_genes_cpm, file = file.path(outdir, paste0(name, "_cpm_cutoff.csv")))
  write.table (resOrdered, file = file.path(outdir, paste0(name, "_rank.rnk")), sep = "\t", row.names = F, quote = F)
  return ( c(nrow(de_genes_rawp), nrow(de_genes_fdr), nrow(de_genes_log2f), nrow(de_genes_cpm)))
}
