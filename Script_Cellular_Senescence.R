# 细胞衰老的相关基因
# cellular senescence-related gene (CSRG)

read_path = "/home/alan_turing/Projects/PD/NM_PD/data/"
save_path = "/home/alan_turing/Projects/PD/ND_CL_EE/result/"

# CellAge
CSRG_CellAge = readr::read_tsv("cellage3.tsv")
CSRG_CellAge = as.data.frame(CSRG_CellAge)
CSRG_CellAge = CSRG_CellAge$`Gene symbol`

# Cell Senescence Database
CSRG_CSGene = read.table("csgene_human.txt", header = TRUE, sep = "\t", dec = ".")
CSRG_CSGene = CSRG_CSGene$GeneSymb

# GeneCards
CSRG_GeneCard = read.csv("GeneCards-SearchResults.csv")
# 选取Relevance Score＞2.5者
CSRG_GeneCard = CSRG_GeneCard[CSRG_GeneCard$Relevance.score >= 2.5,]
CSRG_GeneCard = CSRG_GeneCard$Gene.Symbol

CSRG_total = Reduce(union, c(CSRG_CellAge,CSRG_CSGene,CSRG_GeneCard))

# 加载数据集
library(dplyr)


# 差异基因分析
cohort_loading = function(cohort,type){
  if (type == "cohort"){
    colnames(cohort)[1] = "gene"
    rownames(cohort) = cohort$gene
    cohort$gene = NULL
    return(cohort)
  }
  else if (type == "group"){
    cohort$V1 = NULL
    rownames(cohort) = cohort$sample
    return(cohort)
  }
  else {
    print("Either parameter or option is illegal.")
  }
}

edgeR_DEG_automatic_screen = function(cohort,cohort_group,p_cutoff,foldChange){
  library(edgeR)
  library(dplyr)
  group = c(rep('control', nrow(cohort_group[cohort_group$group == "control",])), rep('PD', nrow(cohort_group[cohort_group$group == "PD",])))
  group = factor(group, levels = c("PD", "control"))
  #创建DGEList对象以存储基因表达数据和组信息
  d <- DGEList(counts = cohort, group = group)
  #根据每个基因的CPM值去除低表达基因
  keep <- rowSums(cpm(d) > 1) >= 2
  #从DGEList对象中筛选出符合条件的基因
  d <- d[keep, , keep.lib.sizes = FALSE]
  #更新样本的库大小信息
  d$samples$lib.size <- colSums(d$counts)
  #归一化（TMM法）并把结果赋值给dge变量
  d <- calcNormFactors(d);dge = d
  #创建设计矩阵以指定差异分析模型
  design <- model.matrix(~0 + factor(group))
  design
  rownames(design) <- colnames(dge)
  colnames(design) <- levels(factor(group))
  # 估计数据的离散度 —— common离散度、trended离散度、tagwise离散度
  dge <- estimateGLMCommonDisp(dge, design)
  dge <- estimateGLMTrendedDisp(dge, design)
  dge <- estimateGLMTagwiseDisp(dge, design)
  # 在估计的模型基础上进行 广义线性模型 (GLM) 拟合
  fit <- glmQLFit(dge, design)
  ## bulk RNA-seq 选择 quasi-likelihood(QL) F-test test
  ## scRNA-seq 或是没有重复样品的数据选用 likelihood ratio test。
  # Design的分组顺序为PD control（设置前须先查看design中的各组顺序，再设置）
  lrt <- glmQLFTest(fit, contrast = c(1, -1))
  ## 如此设置contrast的效果为PD-control
  # 从 LRT 计算结果中获取前 nrow(dge) 个顶部差异表达基因
  nrDEG <- topTags(lrt, n = nrow(dge))
  # 将差异表达基因结果转换为数据框形式
  DEG_edgeR <- as.data.frame(nrDEG)
  # 根据阈值筛选出差异基因
  k1 <- (DEG_edgeR$FDR < p_cutoff) & (DEG_edgeR$logFC < -foldChange)
  k2 <- (DEG_edgeR$FDR < p_cutoff) & (DEG_edgeR$logFC > foldChange)
  DEG_edgeR <- mutate(DEG_edgeR, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))
  table(DEG_edgeR$change)
  #输出结果
  All_diffSig = DEG_edgeR[DEG_edgeR$change != "stable",]
  All_diffSig$gene = rownames(All_diffSig)
  return(All_diffSig)
}

limma_DEG_automatic_screen = function(array_matrix, array_group, foldChange, p_cutoff){
  library(limma)
  library(dplyr)
  group <- model.matrix(~factor(array_group$group)+0)
  colnames(group) <- c("control", "PD")
  limma_fit <- lmFit(array_matrix, group)  ## 数据与list进行匹配
  # 计算logFC、p value等各项指标
  limma_matrix <- makeContrasts(PD - control , levels = group)
  fit <- contrasts.fit(limma_fit, limma_matrix)
  fit <- eBayes(fit)
  tempOutput <- topTable(fit,n = Inf, adjust = "fdr")
  #根据logFC和adjusted p value筛选差异基因
  diffsig = tempOutput
  All_diffSig <- diffsig[(diffsig$P.Value < p_cutoff & (diffsig$logFC>foldChange | diffsig$logFC < (-foldChange))),]
  All_diffSig$gene = rownames(All_diffSig)
  k1 <- (All_diffSig$P.Value < p_cutoff) & (All_diffSig$logFC < -foldChange)
  k2 <- (All_diffSig$P.Value < p_cutoff) & (All_diffSig$logFC > foldChange)
  All_diffSig <- mutate(All_diffSig, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))
  table(All_diffSig$change)
  return(All_diffSig)
}

DESeq2_DEG_automatic_screen = function(cohort,cohort_group,p_cutoff,foldChange){
  # DESeq2法
  library(dplyr)
  group = c(rep('control', nrow(cohort_group[cohort_group$group == "control",])),rep('PD', nrow(cohort_group[cohort_group$group == "PD",])))
  group = factor(group, levels = c("PD", "control"))
  
  library(DESeq2)
  colData <- data.frame(row.names = colnames(cohort),group = group)
  colData$group <- factor(colData$group, levels = c("control", "PD"))
  # 构建DESeqDataSet对象，也就是dds矩阵，将基因计数数据、样本分组信息和设计矩阵关联起来
  dds <- DESeqDataSetFromMatrix(countData = cohort, # 表达矩阵
                                colData = colData,        # 表达矩阵列名和分组信息的对应关系
                                design = ~ group)         # group为colData中的group，也就是分组信息
  # 进行差异表达分析
  dds <- DESeq(dds); resultsNames(dds)
  # 提取差异表达结果（contrast参数指定对比组别，必须写成下面三个元素的向量格式且顺序一致）
  res <- results(dds, contrast = c("group", "PD", "control"))
  # 按照padj（调整后的p值）的大小对差异结果进行排序
  resOrdered <- res[order(res$padj), ]
  # 将差异表达结果转换为数据框
  DEG <- as.data.frame(resOrdered)
  #去除缺失值
  DEG_deseq2 <- na.omit(DEG)
  
  k1 <- (DEG_deseq2$padj < p_cutoff) & (DEG_deseq2$log2FoldChange < -foldChange)
  k2 <- (DEG_deseq2$padj < p_cutoff) & (DEG_deseq2$log2FoldChange > foldChange)
  
  DEG_deseq2 <- mutate(DEG_deseq2, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))
  table(DEG_deseq2$change)
  
  All_diffSig = DEG_deseq2[DEG_deseq2$change != "stable",]
  All_diffSig$gene = rownames(All_diffSig)
  return(All_diffSig)
}

# 加载数据
# 加载数据
GSE8397 = as.data.frame(data.table::fread(paste0(read_path,"GSE8397_GPL96_annotated.txt")))
GSE8397 = cohort_loading(GSE8397,"cohort")

GSE8397_group = as.data.frame(data.table::fread(paste0(read_path,"GSE8397_group_info.txt")))
GSE8397_group = cohort_loading(GSE8397_group,"group")

GSE20168 = as.data.frame(data.table::fread(paste0(read_path,"GSE20168_annotated.txt")))
GSE20168 = cohort_loading(GSE20168,"cohort")

GSE20168_group = as.data.frame(data.table::fread(paste0(read_path,"GSE20168_group_info.txt")))
GSE20168_group = cohort_loading(GSE20168_group,"group")

GSE20291 = as.data.frame(data.table::fread(paste0(read_path,"GSE20291_annotated.txt")))
GSE20291 = cohort_loading(GSE20291,"cohort")

GSE20291_group = as.data.frame(data.table::fread(paste0(read_path,"GSE20291_group_info.txt")))
GSE20291_group = cohort_loading(GSE20291_group,"group")

GSE20292 = as.data.frame(data.table::fread(paste0(read_path,"GSE20292_annotated.txt")))
GSE20292 = cohort_loading(GSE20292,"cohort")

GSE20292_group = as.data.frame(data.table::fread(paste0(read_path,"GSE20292_group_info.txt")))
GSE20292_group = cohort_loading(GSE20292_group,"group")

GSE20163 = as.data.frame(data.table::fread(paste0(read_path,"GSE20163_annotated.txt")))
GSE20163 = cohort_loading(GSE20163,"cohort")

GSE20163_group = as.data.frame(data.table::fread(paste0(read_path,"GSE20163_group_info.txt")))
GSE20163_group = cohort_loading(GSE20163_group,"group")

GSE20164 = as.data.frame(data.table::fread(paste0(read_path,"GSE20164_annotated.txt")))
GSE20164 = cohort_loading(GSE20164,"cohort")

GSE20164_group = as.data.frame(data.table::fread(paste0(read_path,"GSE20164_group_info.txt")))
GSE20164_group = cohort_loading(GSE20164_group,"group")

GSE6613 = as.data.frame(data.table::fread(paste0(read_path,"GSE6613_annotated.txt")))
GSE6613 = cohort_loading(GSE6613,"cohort")

GSE6613_group = as.data.frame(data.table::fread(paste0(read_path,"GSE6613_group_info.txt")))
GSE6613_group = cohort_loading(GSE6613_group,"group")


GSE57475 = as.data.frame(data.table::fread(paste0(read_path,"GSE57475_annotated.txt")))
GSE57475 = cohort_loading(GSE57475,"cohort")

GSE57475_group = as.data.frame(data.table::fread(paste0(read_path,"GSE57475_group_info.txt")))
GSE57475_group = cohort_loading(GSE57475_group,"group")


GSE22491 = as.data.frame(data.table::fread(paste0(read_path,"GSE22491_annotated.txt")))
GSE22491 = cohort_loading(GSE22491,"cohort")

GSE22491_group = as.data.frame(data.table::fread(paste0(read_path,"GSE22491_group_info.txt")))
GSE22491_group = cohort_loading(GSE22491_group,"group")


GSE72267 = as.data.frame(data.table::fread(paste0(read_path,"GSE72267_annotated.txt")))
GSE72267 = cohort_loading(GSE72267,"cohort")

GSE72267_group = as.data.frame(data.table::fread(paste0(read_path,"GSE72267_group_info.txt")))
GSE72267_group = cohort_loading(GSE72267_group,"group")

GSE7621 = as.data.frame(data.table::fread(paste0(read_path,"GSE7621_annotated.txt")))
GSE7621 = cohort_loading(GSE7621,"cohort")

GSE7621_group = as.data.frame(data.table::fread(paste0(read_path,"GSE7621_group_info.txt")))
GSE7621_group = cohort_loading(GSE7621_group,"group")

GSE20141 = as.data.frame(data.table::fread(paste0(read_path,"GSE20141_annotated.txt")))
GSE20141 = cohort_loading(GSE20141,"cohort")

GSE20141_group = as.data.frame(data.table::fread(paste0(read_path,"GSE20141_group_info.txt")))
GSE20141_group = cohort_loading(GSE20141_group,"group")

GSE99039 = as.data.frame(data.table::fread(paste0(read_path,"GSE99039_annotated.txt")))
GSE99039 = cohort_loading(GSE99039,"cohort")

GSE99039_group = as.data.frame(data.table::fread(paste0(read_path,"GSE99039_group_info.txt")))
GSE99039_group = cohort_loading(GSE99039_group,"group")


GSE20146 = as.data.frame(data.table::fread(paste0(read_path,"GSE20146_annotated.txt")))
GSE20146 = cohort_loading(GSE20146,"cohort")

GSE20146_group = as.data.frame(data.table::fread(paste0(read_path,"GSE20146_group_info.txt")))
GSE20146_group = cohort_loading(GSE20146_group,"group")


GSE18838 = as.data.frame(data.table::fread(paste0(read_path,"GSE18838_annotated.txt")))
GSE18838 = cohort_loading(GSE18838,"cohort")

GSE18838_group = as.data.frame(data.table::fread(paste0(read_path,"GSE18838_group_info.txt")))
GSE18838_group = cohort_loading(GSE18838_group,"group")



GSE133347 = as.data.frame(data.table::fread(paste0(read_path,"GSE133347_annotated.txt")))
GSE133347 = cohort_loading(GSE133347,"cohort")

GSE133347_group = as.data.frame(data.table::fread(paste0(read_path,"GSE133347_group_info.txt")))
GSE133347_group = cohort_loading(GSE133347_group,"group")


GSE20333 = as.data.frame(data.table::fread(paste0(read_path,"GSE20333_annotated.txt")))
GSE20333 = cohort_loading(GSE20333,"cohort")

GSE20333_group = as.data.frame(data.table::fread(paste0(read_path,"GSE20333_group_info.txt")))
GSE20333_group = cohort_loading(GSE20333_group,"group")


GSE150696 = as.data.frame(data.table::fread(paste0(read_path,"GSE150696_annotated.txt")))
GSE150696 = cohort_loading(GSE150696,"cohort")

GSE150696_group = as.data.frame(data.table::fread(paste0(read_path,"GSE150696_group_info.txt")))
GSE150696_group = cohort_loading(GSE150696_group,"group")


GSE24378 = as.data.frame(data.table::fread(paste0(read_path,"GSE24378_annotated.txt")))
GSE24378 = cohort_loading(GSE24378,"cohort")

GSE24378_group = as.data.frame(data.table::fread(paste0(read_path,"GSE24378_group_info.txt")))
GSE24378_group = cohort_loading(GSE24378_group,"group")

GSE169755_Count = as.data.frame(data.table::fread(paste0(read_path,"GSE169755_Count_annotated.txt")))
GSE169755_Count = cohort_loading(GSE169755_Count,"cohort")

GSE169755_TPM = as.data.frame(data.table::fread(paste0(read_path,"GSE169755_TPM_annotated.txt")))
GSE169755_TPM = cohort_loading(GSE169755_TPM,"cohort")

GSE169755_group = as.data.frame(data.table::fread(paste0(read_path,"GSE169755_group_info.txt")))
GSE169755_group = cohort_loading(GSE169755_group,"group")

GSE135036_Count = as.data.frame(data.table::fread(paste0(read_path,"GSE135036_Count_annotated.txt")))
GSE135036_Count = cohort_loading(GSE135036_Count,"cohort")

GSE135036_TPM = as.data.frame(data.table::fread(paste0(read_path,"GSE135036_TPM_annotated.txt")))
GSE135036_TPM = cohort_loading(GSE135036_TPM,"cohort")

GSE135036_group = as.data.frame(data.table::fread(paste0(read_path,"GSE135036_group_info.txt")))
GSE135036_group = cohort_loading(GSE135036_group,"group")


GSE182622_Count = as.data.frame(data.table::fread(paste0(read_path,"GSE182622_Count_annotated.txt")))
GSE182622_Count = cohort_loading(GSE182622_Count,"cohort")

GSE182622_TPM = as.data.frame(data.table::fread(paste0(read_path,"GSE182622_TPM_annotated.txt")))
GSE182622_TPM = cohort_loading(GSE182622_TPM,"cohort")

GSE182622_group = as.data.frame(data.table::fread(paste0(read_path,"GSE182622_group_info.txt")))
GSE182622_group = cohort_loading(GSE182622_group,"group")


GSE110716_Count = as.data.frame(data.table::fread(paste0(read_path,"GSE110716_Count_annotated.txt")))
GSE110716_Count = cohort_loading(GSE110716_Count,"cohort")

GSE110716_TPM = as.data.frame(data.table::fread(paste0(read_path,"GSE110716_TPM_annotated.txt")))
GSE110716_TPM = cohort_loading(GSE110716_TPM,"cohort")

GSE110716_group = as.data.frame(data.table::fread(paste0(read_path,"GSE110716_group_info.txt")))
GSE110716_group = cohort_loading(GSE110716_group,"group")


GSE136666_Count = as.data.frame(data.table::fread(paste0(read_path,"GSE136666_Count_annotated.txt")))
GSE136666_Count = cohort_loading(GSE136666_Count,"cohort")

GSE136666_TPM = as.data.frame(data.table::fread(paste0(read_path,"GSE136666_TPM_annotated.txt")))
GSE136666_TPM = cohort_loading(GSE136666_TPM,"cohort")

GSE136666_group = as.data.frame(data.table::fread(paste0(read_path,"GSE136666_group_info.txt")))
GSE136666_group = cohort_loading(GSE136666_group,"group")


GSE68719_Count = as.data.frame(data.table::fread(paste0(read_path,"GSE68719_Count_annotated.txt")))
GSE68719_Count = cohort_loading(GSE68719_Count,"cohort")

GSE68719_TPM = as.data.frame(data.table::fread(paste0(read_path,"GSE68719_TPM_annotated.txt")))
GSE68719_TPM = cohort_loading(GSE68719_TPM,"cohort")

GSE68719_group = as.data.frame(data.table::fread(paste0(read_path,"GSE68719_group_info.txt")))
GSE68719_group = cohort_loading(GSE68719_group,"group")

# 芯片DEG
foldChange = 0.5
p_cutoff = 0.05

GSE8397_DEG = limma_DEG_automatic_screen(GSE8397, GSE8397_group, foldChange = foldChange, p_cutoff = p_cutoff)
GSE20168_DEG = limma_DEG_automatic_screen(GSE20168, GSE20168_group, foldChange = foldChange, p_cutoff = p_cutoff)
GSE20291_DEG = limma_DEG_automatic_screen(GSE20291, GSE20291_group, foldChange = foldChange, p_cutoff = p_cutoff)
GSE20292_DEG = limma_DEG_automatic_screen(GSE20292, GSE20292_group, foldChange = foldChange, p_cutoff = p_cutoff)
GSE20163_DEG = limma_DEG_automatic_screen(GSE20163, GSE20163_group, foldChange = foldChange, p_cutoff = p_cutoff)
GSE20164_DEG = limma_DEG_automatic_screen(GSE20164, GSE20164_group, foldChange = foldChange, p_cutoff = p_cutoff)
GSE6613_DEG = limma_DEG_automatic_screen(GSE6613, GSE6613_group, foldChange = foldChange, p_cutoff = p_cutoff)
GSE57475_DEG = limma_DEG_automatic_screen(GSE57475, GSE57475_group, foldChange = foldChange, p_cutoff = p_cutoff)
GSE22491_DEG = limma_DEG_automatic_screen(GSE22491, GSE22491_group, foldChange = foldChange, p_cutoff = p_cutoff)
GSE72267_DEG = limma_DEG_automatic_screen(GSE72267, GSE72267_group, foldChange = foldChange, p_cutoff = p_cutoff)
GSE7621_DEG = limma_DEG_automatic_screen(GSE7621, GSE7621_group, foldChange = foldChange, p_cutoff = p_cutoff)
GSE20141_DEG = limma_DEG_automatic_screen(GSE20141, GSE20141_group, foldChange = foldChange, p_cutoff = p_cutoff)
GSE99039_DEG = limma_DEG_automatic_screen(GSE99039, GSE99039_group, foldChange = foldChange, p_cutoff = p_cutoff)
GSE20146_DEG = limma_DEG_automatic_screen(GSE20146, GSE20146_group, foldChange = foldChange, p_cutoff = p_cutoff)
GSE18838_DEG = limma_DEG_automatic_screen(GSE18838, GSE18838_group, foldChange = foldChange, p_cutoff = p_cutoff)
GSE133347_DEG = limma_DEG_automatic_screen(GSE133347, GSE133347_group, foldChange = foldChange, p_cutoff = p_cutoff)
GSE20333_DEG = limma_DEG_automatic_screen(GSE20333, GSE20333_group, foldChange = foldChange, p_cutoff = p_cutoff)
GSE150696_DEG = limma_DEG_automatic_screen(GSE150696, GSE150696_group, foldChange = foldChange, p_cutoff = p_cutoff)
GSE24378_DEG = limma_DEG_automatic_screen(GSE24378, GSE24378_group, foldChange = foldChange, p_cutoff = p_cutoff)

GSE169755_DEG = DESeq2_DEG_automatic_screen(GSE169755_Count,GSE169755_group,p_cutoff = p_cutoff, foldChange = foldChange)

GSE182622_DEG = edgeR_DEG_automatic_screen(GSE182622_Count,GSE182622_group,p_cutoff = p_cutoff, foldChange = foldChange)

GSE136666_DEG = DESeq2_DEG_automatic_screen(GSE136666_Count,GSE136666_group,p_cutoff = p_cutoff, foldChange = foldChange)

# GSE8397 + GSE20163
CSDEG = Reduce(intersect,list(CSRG_total,GSE8397_DEG$gene,GSE20163_DEG$gene))

# 火山图（标记基因）
# GSE8397
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tibble)

limma_DEG_automatic_screen_volcano_data = function(array_matrix, array_group, foldChange, padj){
  library(limma)
  group <- model.matrix(~factor(array_group)+0)
  colnames(group) <- c("control", "sepsis")
  limma_fit <- lmFit(array_matrix, group)  ## 数据与list进行匹配
  # 计算logFC、p value等各项指标
  limma_matrix <- makeContrasts(sepsis - control , levels = group)
  fit <- contrasts.fit(limma_fit, limma_matrix)
  fit <- eBayes(fit)
  tempOutput <- topTable(fit,n = Inf, adjust = "fdr")
  #根据logFC和adjusted p value筛选差异基因
  diffsig = tempOutput
  return(diffsig)
}

# GSE8397
GSE8397_volcano_data = limma_DEG_automatic_screen_volcano_data(GSE8397,GSE8397_group$group,foldChange= foldChange, padj = p_cutoff)
GSE8397_k1 <- (GSE8397_volcano_data$P.Value < p_cutoff) & (GSE8397_volcano_data$logFC < -foldChange)
GSE8397_k2 <- (GSE8397_volcano_data$P.Value < p_cutoff) & (GSE8397_volcano_data$logFC > foldChange)
GSE8397_volcano_data <- mutate(GSE8397_volcano_data, change = ifelse(GSE8397_k1, "down", ifelse(GSE8397_k2, "up", "stable")))

GSE8397_volcano_data$gene = rownames(GSE8397_volcano_data)


# 设定火山图纵轴的最高点
GSE8397_volcano_data$yvalue = -log(GSE8397_volcano_data$P.Value)


GSE8397_volcano_plot <- ggplot(data = GSE8397_volcano_data, 
                               aes(x = logFC, 
                                   y = -log10(P.Value))) +
  geom_point(alpha = 0.4, size = 3.5, 
             aes(color = change)) +
  ylab("-log10(P-value)")+
  xlab("log2(FoldChange)")+
  labs(title="GSE8397") + 
  theme(axis.title = element_text(size = 25), 
        axis.ticks = element_text(size = 25),
        legend.key.size = unit(5, "in")) + 
  scale_color_manual(values = c("#339DB5", "grey", "#C9342B"))+
  geom_vline(xintercept = c(-0.5, 0.5), lty = 4, col = "black", lwd = 0.8) +
  geom_hline(yintercept = -log10(0.05), lty = 4, col = "black", lwd = 0.8) +
  theme_bw() + 
  xlim(c(-5, 5)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_label_repel(data = filter(GSE8397_volcano_data, ((logFC > 1.7 & logFC < 2.5) | (logFC < -2.1 & logFC > - 3.5)) & -log10(P.Value) > 2.3), max.overlaps = getOption("ggrepel.max.overlaps", default = 90),  aes(label = gene, color = change), size = 5)
GSE8397_volcano_plot

pdf(paste0(save_path,"GSE8397_volcano_plot.pdf"), width=8, height=6) # 开启图形设备
GSE8397_volcano_plot
dev.off() #关闭图形设备

# GSE20163
GSE20163_volcano_data = limma_DEG_automatic_screen_volcano_data(GSE20163,GSE20163_group$group,foldChange= foldChange, padj = p_cutoff)
GSE20163_k1 <- (GSE20163_volcano_data$P.Value < p_cutoff) & (GSE20163_volcano_data$logFC < -foldChange)
GSE20163_k2 <- (GSE20163_volcano_data$P.Value < p_cutoff) & (GSE20163_volcano_data$logFC > foldChange)
GSE20163_volcano_data <- mutate(GSE20163_volcano_data, change = ifelse(GSE20163_k1, "down", ifelse(GSE20163_k2, "up", "stable")))

GSE20163_volcano_data$gene = rownames(GSE20163_volcano_data)


# 设定火山图纵轴的最高点
GSE20163_volcano_data$yvalue = -log(GSE20163_volcano_data$P.Value)


GSE20163_volcano_plot <- ggplot(data = GSE20163_volcano_data, 
                                aes(x = logFC, 
                                    y = -log10(P.Value))) +
  geom_point(alpha = 0.4, size = 3.5, 
             aes(color = change)) +
  ylab("-log10(P-value)")+
  xlab("log2(FoldChange)")+
  labs(title="GSE20163") + 
  theme(axis.title = element_text(size = 25), 
        axis.ticks = element_text(size = 25),
        legend.key.size = unit(5, "in")) + 
  scale_color_manual(values = c("#339DB5", "grey", "#C9342B"))+
  geom_vline(xintercept = c(-0.5, 0.5), lty = 4, col = "black", lwd = 0.8) +
  geom_hline(yintercept = -log10(0.05), lty = 4, col = "black", lwd = 0.8) +
  theme_bw() + 
  xlim(c(-4, 4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_label_repel(data = filter(GSE20163_volcano_data, ((logFC > 1.2 & logFC < 2) | (logFC < -1.9 & logFC > - 2.9)) & -log10(P.Value) > 2.2), max.overlaps = getOption("ggrepel.max.overlaps", default = 90),  aes(label = gene, color = change), size = 5)
GSE20163_volcano_plot

pdf(paste0(save_path,"GSE20163_volcano_plot.pdf"), width=8, height=6) # 开启图形设备
GSE20163_volcano_plot
dev.off() #关闭图形设备

# 热图（前10个差异基因）
# GSE8397
## 准备数据
GSE8397_heatmap_data_log <- log(GSE8397)+1
GSE8397_volcano_data$gene = rownames(GSE8397_volcano_data)
GSE8397_DEGs_selected_down = GSE8397_volcano_data[GSE8397_volcano_data$change == "down", ]
GSE8397_DEGs_selected_down = GSE8397_DEGs_selected_down[order(GSE8397_DEGs_selected_down$logFC),]

GSE8397_DEGs_selected_up = GSE8397_volcano_data[GSE8397_volcano_data$change == "up", ]
GSE8397_DEGs_selected_up = GSE8397_DEGs_selected_up[order(GSE8397_DEGs_selected_up$logFC,decreasing = T),]

GSE8397_DEGs_selected = rbind(GSE8397_DEGs_selected_up[1:10,],GSE8397_DEGs_selected_down[1:10,])

GSE8397_heatmap_data_log_selected <- GSE8397_heatmap_data_log[GSE8397_DEGs_selected$gene,]

## 准备分组信息
GSE8397_annotation_col = GSE8397_group
rownames(GSE8397_annotation_col) = GSE8397_annotation_col$sample
GSE8397_annotation_col$sample = NULL

library(pheatmap)
GSE8397_heatmap_color = colorRampPalette(c("navy", "white", "firebrick3"))(100)
GSE8397_ann_colors = list(group = c(control = "#339DB5", PD = "#C9342B"))
GSE8397_heatmap_DEGs <- pheatmap::pheatmap(GSE8397_heatmap_data_log_selected,
                                           show_colnames = F,
                                           show_rownames = T,
                                           scale = "row",
                                           cluster_cols = F,
                                           annotation_col = GSE8397_annotation_col,
                                           breaks = seq(-3, 3, length.out = 100),
                                           annotation_legend = TRUE,
                                           border_color = NA,
                                           annotation_colors = GSE8397_ann_colors,
                                           fontsize_row = 6.5,
                                           main = "GSE8397")
GSE8397_heatmap_DEGs
#保存结果
pdf(paste0(save_path,"GSE8397_heatmap.pdf"), width=6, height=6)
GSE8397_heatmap_DEGs
dev.off()


# GSE20163

GSE20163_heatmap_data_log <- log(GSE20163)+1
GSE20163_volcano_data$gene = rownames(GSE20163_volcano_data)
GSE20163_DEGs_selected_down = GSE20163_volcano_data[GSE20163_volcano_data$change == "down", ]
GSE20163_DEGs_selected_down = GSE20163_DEGs_selected_down[order(GSE20163_DEGs_selected_down$logFC),]

GSE20163_DEGs_selected_up = GSE20163_volcano_data[GSE20163_volcano_data$change == "up", ]
GSE20163_DEGs_selected_up = GSE20163_DEGs_selected_up[order(GSE20163_DEGs_selected_up$logFC,decreasing = T),]

GSE20163_DEGs_selected = rbind(GSE20163_DEGs_selected_up[1:10,],GSE20163_DEGs_selected_down[1:10,])

GSE20163_heatmap_data_log_selected <- GSE20163_heatmap_data_log[GSE20163_DEGs_selected$gene,]

## 准备分组信息
GSE20163_annotation_col = GSE20163_group
rownames(GSE20163_annotation_col) = GSE20163_annotation_col$sample
GSE20163_annotation_col$sample = NULL

library(pheatmap)
GSE20163_heatmap_color = colorRampPalette(c("navy", "white", "firebrick3"))(100)
GSE20163_ann_colors = list(group = c(control = "#339DB5", PD = "#C9342B"))
GSE20163_heatmap_DEGs <- pheatmap::pheatmap(GSE20163_heatmap_data_log_selected,
                                            show_colnames = F,
                                            show_rownames = T,
                                            scale = "row",
                                            cluster_cols = F,
                                            annotation_col = GSE20163_annotation_col,
                                            breaks = seq(-3, 3, length.out = 100),
                                            annotation_legend = TRUE,
                                            border_color = NA,
                                            annotation_colors = GSE20163_ann_colors,
                                            fontsize_row = 6.5,
                                            main = "GSE20163")
GSE20163_heatmap_DEGs
#保存结果
pdf(paste0(save_path,"GSE20163_heatmap.pdf"), width=6, height=6)
GSE20163_heatmap_DEGs
dev.off()

# 富集分析
# GSE8397
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(msigdbr)
library(GSEABase)
library(GseaVis)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(enrichplot)

GSE8397_All_diffSig = GSE8397_volcano_data[GSE8397_volcano_data$change != "stable",]

GSE8397_All_diffSig$gene = rownames(GSE8397_All_diffSig)
GSE8397_IDs_symbol_and_ENTREZID <- bitr(geneID = GSE8397_All_diffSig$gene, 
                                        fromType = "SYMBOL", 
                                        toType = "ENTREZID", 
                                        OrgDb = org.Hs.eg.db)

GSE8397_All_diffSig$EZTREZID <- GSE8397_IDs_symbol_and_ENTREZID[match(GSE8397_All_diffSig$gene,GSE8397_IDs_symbol_and_ENTREZID$SYMBOL),2]
GSE8397_All_diffSig <- na.omit(GSE8397_All_diffSig)

# GO富集分析
GSE8397_GO_enrich = enrichGO(gene = GSE8397_All_diffSig$EZTREZID, #待富集的基因列表
                             OrgDb = org.Hs.eg.db, 
                             keyType = "ENTREZID", #输入数据的类型
                             ont = "ALL", #可以输入CC/MF/BP/ALL
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.05,
                             readable = T)

GSE8397_GO_enrich_result  = data.frame(GSE8397_GO_enrich)


# dot plot
library("enrichplot")

pdf(paste0(save_path,"GSE8397_GO_enrichment_dotplot.pdf"),width = 7,height = 11) 
dotplot(GSE8397_GO_enrich, x = "GeneRatio", split="ONTOLOGY",showCategory = 5) + facet_grid(ONTOLOGY~., scale="free")
dev.off()

# KEGG富集分析
# GSE8397
GSE8397_KEGG_enrich = enrichKEGG(gene = GSE8397_All_diffSig$EZTREZID, #待富集的基因列表
                                 keyType = "kegg",
                                 organism= "human", #可根据你自己要研究的物种更改（参阅https://www.kegg.jp/brite/br08611）
                                 qvalueCutoff = 0.05, #指定p值阈值
                                 pvalueCutoff = 0.05) #指定q值阈值

# dot plot
library("enrichplot")
pdf(paste0(save_path,"GSE8397_KEGG_enrichment_dotplot.pdf"),width = 7,height = 10) 
dotplot(GSE8397_KEGG_enrich, x = "GeneRatio", color = "p.adjust", #默认参数
        showCategory = 10) #只显示前5
dev.off()

# chord diagram
library(GOplot)
library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)

source("Helper_v2.R")
## 准备数据

#将富集分析结果中的ENTREZID再转换成基因名字
GSE8397_KEGG_enrich_result=setReadable(GSE8397_KEGG_enrich,OrgDb=org.Hs.eg.db,keyType  ="ENTREZID")

#从富集结果中挑选"ONTOLOGY"，"ID"，"Description"，"geneID"，"p.adjust" 
GSE8397_kegg=data.frame(GSE8397_KEGG_enrich_result)[,c("ID","Description","geneID","p.adjust" )]

#KEGG富集分析结果中没有category这里一列，这里我们自己造一个
GSE8397_kegg$Category=rep("kegg",nrow(GSE8397_kegg))

#重新给列命名，满足circle_dat函数的要求
names(GSE8397_kegg)=c("ID","Term","Genes", "adj_pval","Category")
#将/隔开的基因名字转换成,隔开
GSE8397_kegg$Genes=gsub("/",",",GSE8397_kegg$Genes)
#将差异表达分析结果的第一列symbol改成ID，满足circle_dat函数的要求
GSE8397_All_diffSig$ID = rownames(GSE8397_All_diffSig)
GSE8397_All_diffSig = GSE8397_All_diffSig %>% dplyr::select("ID",everything())

#创建circ对象，第一个参数为KEGG富集分析结果，第二个参数为差异表达分析结果
GSE8397_circ <- circle_dat(GSE8397_kegg, GSE8397_All_diffSig)

#根据差异表达基因的logFC排序
GSE8397_genes=GSE8397_All_diffSig$logFC
names(GSE8397_genes)=GSE8397_All_diffSig$ID
GSE8397_genes=sort(GSE8397_genes)
#挑选logFC最大和最小的各100各基因，100这个数字可以自己修改
#也可以自己挑选基因，需要带上logFC
GSE8397_chord_genes=c(head(GSE8397_genes,250),tail(GSE8397_genes,250))
#构建数据框，第一列为基因名字，第二列为对应的logFC
GSE8397_chord_genes_df=data.frame(ID=names(GSE8397_chord_genes),logFC=GSE8397_chord_genes)


#取前10个显著附件的KEGG通路，这个10可以自己修改
GSE8397_process=head(GSE8397_kegg$Term,10)

#构建chord对象，第一个参数为circ对象，第二个参数为带有logFC的基因名字，
#第三个是想展示的KEGG通路的名字
GSE8397_chord <- chord_dat(GSE8397_circ, GSE8397_chord_genes_df, GSE8397_process)
GSE8397_color = c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","#DC0000FF","#7E6148FF","#B09C85FF")
#保存图片
pdf(paste0(save_path,"GSE8397_KEGG_enrichment_chord_diagram.pdf"),width = 20,height = 16)
KEGGChord(GSE8397_chord,   #chord对象
          limit=c(1,3), #第一个数每个基因至少需要根几个term相连，第二个数每个term至少需要根几个基因相连
          space = 0.02,  #右侧色块之间的间距
          gene.order = 'logFC',   #基因展示顺序根据logFC来
          gene.space = 0.25,  #基因名字和色块之间的距离
          gene.size = 5,
          ribbon.col = GSE8397_color, #颜色调整（根据你的通路个数-1确定）
          process.label = 15)
dev.off()

# GSEA
## 1.1 在data中添加一列：ENTREZID，并使数据框按照log2FoldChange序列
GSE8397_All_diffSig$gene = rownames(GSE8397_All_diffSig)
GSE8397_IDs_symbol_and_ENTREZID <- bitr(geneID = GSE8397_All_diffSig$gene, 
                                        fromType = "SYMBOL", 
                                        toType = "ENTREZID", 
                                        OrgDb = org.Hs.eg.db)

GSE8397_All_diffSig$EZTREZID <- GSE8397_IDs_symbol_and_ENTREZID[match(GSE8397_All_diffSig$gene,GSE8397_IDs_symbol_and_ENTREZID$SYMBOL),2]
GSE8397_All_diffSig <- na.omit(GSE8397_All_diffSig)
GSE8397_All_diffSig <- GSE8397_All_diffSig[order(GSE8397_All_diffSig$logFC,decreasing = T),]
## 1.2 构建包含ENTREZID名称的log2FoldChange向量
GSE8397_gene_list <- GSE8397_All_diffSig$logFC 
names(GSE8397_gene_list) <- GSE8397_All_diffSig$EZTREZID
# 2.GSEA分析
## 2.1 GSEA分析
GSE8397_gsea_start.time <- Sys.time()
GSE8397_gsea_start.time

GSE8397_GSEAresult <- gseGO(geneList = GSE8397_gene_list, 
                            OrgDb = org.Hs.eg.db, #相应物种的数据库
                            ont = "ALL", #可以选择ALL, BP, CC, MF
                            nPerm = 1000, 
                            pvalueCutoff = 0.05,
                            verbose = TRUE, 
                            by = "DOSE")
# GSEAresult <- gseKEGG(geneList = gene_list,  organism = "hsa", keyType = "kegg",  verbose = TRUE, pvalueCutoff = 0.05, by = "fgsea")

GSE8397_gsea_end.time <- Sys.time()
GSE8397_gsea_end.time


GSE8397_gsea_time.taken <- GSE8397_gsea_end.time - GSE8397_gsea_start.time
GSE8397_gsea_time.taken

## 2.2 输出并提取GSEA分析结果
GSE8397_GSEA_result = GSE8397_GSEAresult@result

## 2.3 可视化并保存结果（分别取NES＞0和NES＜0时的Top 5予以可视化）
pdf(paste0(save_path,"GSE8397_GSEA_positive_NES_top_5.pdf"),width = 6, height= 4)
gseaplot2(x = GSE8397_GSEAresult,
          title = "GSEA in PD (GSE8397)",
          geneSetID = c("GO:0061077","GO:0006457","GO:0061687","GO:0097501","GO:0071214"),
          pvalue_table = F)
dev.off()

pdf(paste0(save_path,"GSE8397_GSEA_negative_NES_top_5.pdf"),width = 6, height= 4)
gseaplot2(x = GSE8397_GSEAresult,
          title = "GSEA in Control (GSE8397)",
          geneSetID = c("GO:0098793","GO:0150034","GO:0007268","GO:0098916","GO:0044306"),
          pvalue_table = F)
dev.off()

# GSVA
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(msigdbr)
library(GSEABase)
library(GSVA)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

# 2.加载数据
# 2.1 原始的表达矩阵
GSE8397_GSVA_data = GSE8397
## 2.2 构建meta data,即sample与group的对应关系
#GSVA_meta <- data.frame(sample = colnames(GSVA_data), group = group_gsva, function(x) x[1])
GSE8397_GSVA_meta = GSE8397_group

## 3.1 基因集数据下载（一般会选择GO和【或】KEGG做GSVA分析）
#后续选择GO进行GSVA分析
GO_df_all <- msigdbr(species = "Homo sapiens", 
                     # Homo sapiens or Mus musculus
                     category = "C5")
GO_df <- dplyr::select(GO_df_all, gs_name, gene_symbol, gs_exact_source, gs_subcat)
GO_df <- GO_df[GO_df$gs_subcat!="HPO",]
#按照gs_name给gene_symbol分
go_list <- split(GO_df$gene_symbol, GO_df$gs_name)

##3.2 基因表达矩阵
GSE8397_GSVA_data <- as.matrix(GSE8397_GSVA_data)

##3.3 GSVA分析
GSE8397_gsva_mat <- gsvaParam(exprData=GSE8397_GSVA_data, geneSets=go_list, kcdf="Gaussian", minSize = 2)
#"Gaussian" for microarray fluorescent units in logarithmic scale, RNA-seq log-CPMs, log-RPKMs or log-TPMs.
#"Poisson" for interger counts.
GSE8397_gsva_mat = gsva(GSE8397_gsva_mat,verbose = TRUE)

library(limma)
# 设置或导入分组
rownames(GSE8397_GSVA_meta) = GSE8397_GSVA_meta$sample
GSE8397_GSVA_meta$sample = NULL
GSE8397_GSVA_meta = factor(GSE8397_GSVA_meta$group)

GSE8397_group_gsva = GSE8397_GSVA_meta
GSE8397_design_gsva <- model.matrix(~0+GSE8397_group_gsva)
colnames(GSE8397_design_gsva) = levels(factor(GSE8397_group_gsva))
rownames(GSE8397_design_gsva) = colnames(GSE8397_gsva_mat)
GSE8397_design_gsva
# PD VS control
GSE8397_compare_gsva <- makeContrasts(PD - control, levels=GSE8397_design_gsva)
GSE8397_fit_gsva <- lmFit(GSE8397_gsva_mat, GSE8397_design_gsva)
GSE8397_fit2_gsva <- contrasts.fit(GSE8397_fit_gsva, GSE8397_compare_gsva)
GSE8397_fit3_gsva <- eBayes(GSE8397_fit2_gsva)
GSE8397_Diff_gsva <- topTable(GSE8397_fit3_gsva, coef=1, number=200)
head(GSE8397_Diff_gsva)

# 可视化（柱形偏差图）
library(ggprism)
library(ggplot2)
library(tidyverse)

#数据准备
GSE8397_gsva_limma_visual= GSE8397_Diff_gsva
#去掉"GOBP_"等前缀
library(stringr)
GSE8397_gsva_limma_visual$id = rownames(GSE8397_gsva_limma_visual)
rownames(GSE8397_gsva_limma_visual) = NULL
GSE8397_gsva_limma_visual$id <- str_replace(GSE8397_gsva_limma_visual$id , "GOBP_","")
GSE8397_gsva_limma_visual$id <- str_replace(GSE8397_gsva_limma_visual$id , "GOCC_","")
GSE8397_gsva_limma_visual$id <- str_replace(GSE8397_gsva_limma_visual$id , "GOMF_","")
GSE8397_gsva_limma_visual = dplyr::select(GSE8397_gsva_limma_visual,"id",everything())
# 新增一列（根据t阈值分类，以2为准）
GSE8397_gsva_limma_visual$threshold = factor(ifelse(GSE8397_gsva_limma_visual$t  >-2, ifelse(GSE8397_gsva_limma_visual$t >= 2 ,'Up','NoSignifi'),'Down'),levels=c('Up','Down','NoSignifi'))
# 排序
GSE8397_gsva_limma_visual <- GSE8397_gsva_limma_visual %>% arrange(t)
# 变成因子类型
GSE8397_gsva_limma_visual$id <- factor(GSE8397_gsva_limma_visual$id,levels = GSE8397_gsva_limma_visual$id)
backup_GSE8397_gsva_limma_visual = GSE8397_gsva_limma_visual

# 由于GSVA结果太多，只保留特定范围的t值的以作可视化（后续代码复用时注意自行修改）
# GSE8397_gsva_limma_visual = backup_GSE8397_gsva_limma_visual #恢复数据以作调整
GSE8397_gsva_limma_visual_negative = GSE8397_gsva_limma_visual[1:10,]
GSE8397_gsva_limma_visual_positive = GSE8397_gsva_limma_visual[191:200,]
GSE8397_gsva_limma_visual = rbind(GSE8397_gsva_limma_visual_negative,GSE8397_gsva_limma_visual_positive)
rm(GSE8397_gsva_limma_visual_negative,GSE8397_gsva_limma_visual_positive)

# 绘图
GSE8397_gsva_limma_picture <- ggplot(data = GSE8397_gsva_limma_visual,aes(x = id,y = t,fill = threshold)) +
  geom_col() +
  xlab('') + ylab('') +
  # 主题
  theme_prism(border = T) +
  # 填充颜色
  scale_fill_manual(values = c('Up'= '#C9342B','NoSignifi'='grey','Down'='#339DB5')) +
  #scale_fill_gradient2(low = '#FFC074',mid = '#B6C867',high = '#01937C') +
  # 竖线
  geom_hline(yintercept = c(-1,1),color = 'white',linewidth = 1,lty='dashed') +
  # 翻转坐标轴
  coord_flip() + ylim(-10,10) +
  # 去除坐标轴标签
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        # 移动图例位置
        legend.position.inside = c(0.1,0.9)) +
  ylab("t value of GSVA score, PD verus control in GSE8397")

#针对The legend.text.align theme element is not defined in the element hierarchy错误的解决方法
GSE8397_gsva_limma_picture$theme[c("legend.text.align", "legend.title.align")] <- NULL

# 添加标签在图里

# 小于-2的数量
low1 <- GSE8397_gsva_limma_visual %>% filter(t < -2) %>% nrow()
# 小于0总数量
low0 <- GSE8397_gsva_limma_visual %>% filter(t < 0) %>% nrow()
# 小于2总数量
high0 <- GSE8397_gsva_limma_visual %>% filter(t < 2) %>% nrow()
# 总的柱子数量
high1 <- nrow(GSE8397_gsva_limma_visual)

# 依次从下到上添加标签
## 此处删除了灰色标签（原始数据中无NoSignifi分组者，而且不删就放不下），后续代码复用时注意自行修改
GSE8397_gsva_limma_picture = GSE8397_gsva_limma_picture + geom_text(data = GSE8397_gsva_limma_visual[1:low1,],aes(x = id,y = 0.1,label = id),
                                                                    hjust = 0,color = 'black') + # 小于-1的为黑色标签
  #geom_text(data = GSE8397_gsva_limma_visual[(low1 +1):low0,],aes(x = id,y = 0.1,label = id), hjust = 0,color = 'grey') + # 灰色标签
  # geom_text(data = GSE8397_gsva_limma_visual[(low0 + 1):high0,],aes(x = id,y = -0.1,label = id), hjust = 1,color = 'grey') + # 灰色标签
  geom_text(data = GSE8397_gsva_limma_visual[(high0 +1):high1,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'black') # 大于1的为黑色标签
pdf(paste0(save_path,"GSE8397_gsva_barplot.pdf"), width = 13.5, height = 6)
GSE8397_gsva_limma_picture
dev.off()

# GSE20163
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(msigdbr)
library(GSEABase)
library(GseaVis)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(enrichplot)

GSE20163_All_diffSig = GSE20163_volcano_data[GSE20163_volcano_data$change != "stable",]

GSE20163_All_diffSig$gene = rownames(GSE20163_All_diffSig)
GSE20163_IDs_symbol_and_ENTREZID <- bitr(geneID = GSE20163_All_diffSig$gene, 
                                        fromType = "SYMBOL", 
                                        toType = "ENTREZID", 
                                        OrgDb = org.Hs.eg.db)

GSE20163_All_diffSig$EZTREZID <- GSE20163_IDs_symbol_and_ENTREZID[match(GSE20163_All_diffSig$gene,GSE20163_IDs_symbol_and_ENTREZID$SYMBOL),2]
GSE20163_All_diffSig <- na.omit(GSE20163_All_diffSig)

# GO富集分析
GSE20163_GO_enrich = enrichGO(gene = GSE20163_All_diffSig$EZTREZID, #待富集的基因列表
                             OrgDb = org.Hs.eg.db, 
                             keyType = "ENTREZID", #输入数据的类型
                             ont = "ALL", #可以输入CC/MF/BP/ALL
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.05,
                             readable = T)

GSE20163_GO_enrich_result  = data.frame(GSE20163_GO_enrich)


# dot plot
library("enrichplot")

pdf(paste0(save_path,"GSE20163_GO_enrichment_dotplot.pdf"),width = 7,height = 11) 
dotplot(GSE20163_GO_enrich, x = "GeneRatio", split="ONTOLOGY",showCategory = 5) + facet_grid(ONTOLOGY~., scale="free")
dev.off()

# KEGG富集分析
# GSE20163
GSE20163_KEGG_enrich = enrichKEGG(gene = GSE20163_All_diffSig$EZTREZID, #待富集的基因列表
                                 keyType = "kegg",
                                 organism= "human", #可根据你自己要研究的物种更改（参阅https://www.kegg.jp/brite/br08611）
                                 qvalueCutoff = 0.05, #指定p值阈值
                                 pvalueCutoff = 0.05) #指定q值阈值

# chord diagram
library(GOplot)
library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)

source("Helper_v2.R")
## 准备数据

#将富集分析结果中的ENTREZID再转换成基因名字
GSE20163_KEGG_enrich_result=setReadable(GSE20163_KEGG_enrich,OrgDb=org.Hs.eg.db,keyType  ="ENTREZID")

#从富集结果中挑选"ONTOLOGY"，"ID"，"Description"，"geneID"，"p.adjust" 
GSE20163_kegg=data.frame(GSE20163_KEGG_enrich_result)[,c("ID","Description","geneID","p.adjust" )]

#KEGG富集分析结果中没有category这里一列，这里我们自己造一个
GSE20163_kegg$Category=rep("kegg",nrow(GSE20163_kegg))

#重新给列命名，满足circle_dat函数的要求
names(GSE20163_kegg)=c("ID","Term","Genes", "adj_pval","Category")
#将/隔开的基因名字转换成,隔开
GSE20163_kegg$Genes=gsub("/",",",GSE20163_kegg$Genes)
#将差异表达分析结果的第一列symbol改成ID，满足circle_dat函数的要求
GSE20163_All_diffSig$ID = rownames(GSE20163_All_diffSig)
GSE20163_All_diffSig = GSE20163_All_diffSig %>% dplyr::select("ID",everything())

#创建circ对象，第一个参数为KEGG富集分析结果，第二个参数为差异表达分析结果
GSE20163_circ <- circle_dat(GSE20163_kegg, GSE20163_All_diffSig)

#根据差异表达基因的logFC排序
GSE20163_genes=GSE20163_All_diffSig$logFC
names(GSE20163_genes)=GSE20163_All_diffSig$ID
GSE20163_genes=sort(GSE20163_genes)
#挑选logFC最大和最小的各100各基因，100这个数字可以自己修改
#也可以自己挑选基因，需要带上logFC
GSE20163_chord_genes=c(head(GSE20163_genes,250),tail(GSE20163_genes,250))
#构建数据框，第一列为基因名字，第二列为对应的logFC
GSE20163_chord_genes_df=data.frame(ID=names(GSE20163_chord_genes),logFC=GSE20163_chord_genes)


#取前10个显著附件的KEGG通路，这个10可以自己修改
GSE20163_process=head(GSE20163_kegg$Term,10)

#构建chord对象，第一个参数为circ对象，第二个参数为带有logFC的基因名字，
#第三个是想展示的KEGG通路的名字
source("Helper_v2.R")
GSE20163_chord <- chord_dat(GSE20163_circ, GSE20163_chord_genes_df, GSE20163_process)
GSE20163_color = c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","#DC0000FF","#7E6148FF")

#保存图片
pdf(paste0(save_path,"GSE20163_KEGG_enrichment_chord_diagram.pdf"),width = 20,height = 16)
KEGGChord(GSE20163_chord,   #chord对象
          limit=c(1,3), #第一个数每个基因至少需要根几个term相连，第二个数每个term至少需要根几个基因相连
          space = 0.02,  #右侧色块之间的间距
          gene.order = 'logFC',   #基因展示顺序根据logFC来
          gene.space = 0.25,  #基因名字和色块之间的距离
          gene.size = 5,
          ribbon.col = GSE20163_color, #颜色调整（根据你的通路个数-1确定）
          process.label = 15)
dev.off()

# GSEA
## 1.1 在data中添加一列：ENTREZID，并使数据框按照log2FoldChange序列
GSE20163_All_diffSig$gene = rownames(GSE20163_All_diffSig)
GSE20163_IDs_symbol_and_ENTREZID <- bitr(geneID = GSE20163_All_diffSig$gene, 
                                        fromType = "SYMBOL", 
                                        toType = "ENTREZID", 
                                        OrgDb = org.Hs.eg.db)

GSE20163_All_diffSig$EZTREZID <- GSE20163_IDs_symbol_and_ENTREZID[match(GSE20163_All_diffSig$gene,GSE20163_IDs_symbol_and_ENTREZID$SYMBOL),2]
GSE20163_All_diffSig <- na.omit(GSE20163_All_diffSig)
GSE20163_All_diffSig <- GSE20163_All_diffSig[order(GSE20163_All_diffSig$logFC,decreasing = T),]
## 1.2 构建包含ENTREZID名称的log2FoldChange向量
GSE20163_gene_list <- GSE20163_All_diffSig$logFC 
names(GSE20163_gene_list) <- GSE20163_All_diffSig$EZTREZID
# 2.GSEA分析
## 2.1 GSEA分析
GSE20163_gsea_start.time <- Sys.time()
GSE20163_gsea_start.time

GSE20163_GSEAresult <- gseGO(geneList = GSE20163_gene_list, 
                            OrgDb = org.Hs.eg.db, #相应物种的数据库
                            ont = "ALL", #可以选择ALL, BP, CC, MF
                            nPerm = 1000, 
                            pvalueCutoff = 0.05,
                            verbose = TRUE, 
                            by = "DOSE")
# GSEAresult <- gseKEGG(geneList = gene_list,  organism = "hsa", keyType = "kegg",  verbose = TRUE, pvalueCutoff = 0.05, by = "fgsea")

GSE20163_gsea_end.time <- Sys.time()
GSE20163_gsea_end.time


GSE20163_gsea_time.taken <- GSE20163_gsea_end.time - GSE20163_gsea_start.time
GSE20163_gsea_time.taken

## 2.2 输出并提取GSEA分析结果
GSE20163_GSEA_result = GSE20163_GSEAresult@result

# ## 2.3 可视化并保存结果（分别取NES＞0和NES＜0时的Top 5予以可视化）
# pdf(paste0(save_path,"GSE20163_GSEA_positive_NES_top_5.pdf"),width = 6, height= 4)
# gseaplot2(x = GSE20163_GSEAresult,
#           title = "GSEA in PD (GSE20163)",
#           geneSetID = c("GO:0061077","GO:0006457","GO:0061687","GO:0097501","GO:0071214"),
#           pvalue_table = F)
# dev.off()
# 
# pdf(paste0(save_path,"GSE20163_GSEA_negative_NES_top_5.pdf"),width = 6, height= 4)
# gseaplot2(x = GSE20163_GSEAresult,
#           title = "GSEA in Control (GSE20163)",
#           geneSetID = c("GO:0098793","GO:0150034","GO:0007268","GO:0098916","GO:0044306"),
#           pvalue_table = F)
# dev.off()

# GSVA
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(msigdbr)
library(GSEABase)
library(GSVA)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

# 2.加载数据
# 2.1 原始的表达矩阵
GSE20163_GSVA_data = GSE20163
## 2.2 构建meta data,即sample与group的对应关系
#GSVA_meta <- data.frame(sample = colnames(GSVA_data), group = group_gsva, function(x) x[1])
GSE20163_GSVA_meta = GSE20163_group

## 3.1 基因集数据下载（一般会选择GO和【或】KEGG做GSVA分析）
#后续选择GO进行GSVA分析
GO_df_all <- msigdbr(species = "Homo sapiens", 
                     # Homo sapiens or Mus musculus
                     category = "C5")
GO_df <- dplyr::select(GO_df_all, gs_name, gene_symbol, gs_exact_source, gs_subcat)
GO_df <- GO_df[GO_df$gs_subcat!="HPO",]
#按照gs_name给gene_symbol分
go_list <- split(GO_df$gene_symbol, GO_df$gs_name)

##3.2 基因表达矩阵
GSE20163_GSVA_data <- as.matrix(GSE20163_GSVA_data)

##3.3 GSVA分析
GSE20163_gsva_mat <- gsvaParam(exprData=GSE20163_GSVA_data, geneSets=go_list, kcdf="Gaussian", minSize = 2)
#"Gaussian" for microarray fluorescent units in logarithmic scale, RNA-seq log-CPMs, log-RPKMs or log-TPMs.
#"Poisson" for interger counts.
GSE20163_gsva_mat = gsva(GSE20163_gsva_mat,verbose = TRUE)

library(limma)
# 设置或导入分组
rownames(GSE20163_GSVA_meta) = GSE20163_GSVA_meta$sample
GSE20163_GSVA_meta$sample = NULL
GSE20163_GSVA_meta = factor(GSE20163_GSVA_meta$group)

GSE20163_group_gsva = GSE20163_GSVA_meta
GSE20163_design_gsva <- model.matrix(~0+GSE20163_group_gsva)
colnames(GSE20163_design_gsva) = levels(factor(GSE20163_group_gsva))
rownames(GSE20163_design_gsva) = colnames(GSE20163_gsva_mat)
GSE20163_design_gsva
# PD VS control
GSE20163_compare_gsva <- makeContrasts(PD - control, levels=GSE20163_design_gsva)
GSE20163_fit_gsva <- lmFit(GSE20163_gsva_mat, GSE20163_design_gsva)
GSE20163_fit2_gsva <- contrasts.fit(GSE20163_fit_gsva, GSE20163_compare_gsva)
GSE20163_fit3_gsva <- eBayes(GSE20163_fit2_gsva)
GSE20163_Diff_gsva <- topTable(GSE20163_fit3_gsva, coef=1, number=200)
head(GSE20163_Diff_gsva)

# 可视化（柱形偏差图）
library(ggprism)
library(ggplot2)
library(tidyverse)

#数据准备
GSE20163_gsva_limma_visual= GSE20163_Diff_gsva
#去掉"GOBP_"等前缀
library(stringr)
GSE20163_gsva_limma_visual$id = rownames(GSE20163_gsva_limma_visual)
rownames(GSE20163_gsva_limma_visual) = NULL
GSE20163_gsva_limma_visual$id <- str_replace(GSE20163_gsva_limma_visual$id , "GOBP_","")
GSE20163_gsva_limma_visual$id <- str_replace(GSE20163_gsva_limma_visual$id , "GOCC_","")
GSE20163_gsva_limma_visual$id <- str_replace(GSE20163_gsva_limma_visual$id , "GOMF_","")
GSE20163_gsva_limma_visual = dplyr::select(GSE20163_gsva_limma_visual,"id",everything())
# 新增一列（根据t阈值分类，以2为准）
GSE20163_gsva_limma_visual$threshold = factor(ifelse(GSE20163_gsva_limma_visual$t  >-2, ifelse(GSE20163_gsva_limma_visual$t >= 2 ,'Up','NoSignifi'),'Down'),levels=c('Up','Down','NoSignifi'))
# 排序
GSE20163_gsva_limma_visual <- GSE20163_gsva_limma_visual %>% arrange(t)
# 变成因子类型
GSE20163_gsva_limma_visual$id <- factor(GSE20163_gsva_limma_visual$id,levels = GSE20163_gsva_limma_visual$id)
backup_GSE20163_gsva_limma_visual = GSE20163_gsva_limma_visual

# 由于GSVA结果太多，只保留特定范围的t值的以作可视化（后续代码复用时注意自行修改）
# GSE20163_gsva_limma_visual = backup_GSE20163_gsva_limma_visual #恢复数据以作调整
GSE20163_gsva_limma_visual_negative = GSE20163_gsva_limma_visual[1:10,]
GSE20163_gsva_limma_visual_positive = GSE20163_gsva_limma_visual[191:200,]
GSE20163_gsva_limma_visual = rbind(GSE20163_gsva_limma_visual_negative,GSE20163_gsva_limma_visual_positive)
rm(GSE20163_gsva_limma_visual_negative,GSE20163_gsva_limma_visual_positive)

# 绘图
GSE20163_gsva_limma_picture <- ggplot(data = GSE20163_gsva_limma_visual,aes(x = id,y = t,fill = threshold)) +
  geom_col() +
  xlab('') + ylab('') +
  # 主题
  theme_prism(border = T) +
  # 填充颜色
  scale_fill_manual(values = c('Up'= '#C9342B','NoSignifi'='grey','Down'='#339DB5')) +
  #scale_fill_gradient2(low = '#FFC074',mid = '#B6C867',high = '#01937C') +
  # 竖线
  geom_hline(yintercept = c(-1,1),color = 'white',linewidth = 1,lty='dashed') +
  # 翻转坐标轴
  coord_flip() + ylim(-10,10) +
  # 去除坐标轴标签
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        # 移动图例位置
        legend.position.inside = c(0.1,0.9)) +
  ylab("t value of GSVA score, PD verus control in GSE20163")

#针对The legend.text.align theme element is not defined in the element hierarchy错误的解决方法
GSE20163_gsva_limma_picture$theme[c("legend.text.align", "legend.title.align")] <- NULL

# 添加标签在图里

# 小于-2的数量
low1 <- GSE20163_gsva_limma_visual %>% filter(t < -2) %>% nrow()
# 小于0总数量
low0 <- GSE20163_gsva_limma_visual %>% filter(t < 0) %>% nrow()
# 小于2总数量
high0 <- GSE20163_gsva_limma_visual %>% filter(t < 2) %>% nrow()
# 总的柱子数量
high1 <- nrow(GSE20163_gsva_limma_visual)

# 依次从下到上添加标签
## 此处删除了灰色标签（原始数据中无NoSignifi分组者，而且不删就放不下），后续代码复用时注意自行修改
GSE20163_gsva_limma_picture = GSE20163_gsva_limma_picture + geom_text(data = GSE20163_gsva_limma_visual[1:low1,],aes(x = id,y = 0.1,label = id),
                                                                    hjust = 0,color = 'black') + # 小于-1的为黑色标签
  #geom_text(data = GSE20163_gsva_limma_visual[(low1 +1):low0,],aes(x = id,y = 0.1,label = id), hjust = 0,color = 'grey') + # 灰色标签
  # geom_text(data = GSE20163_gsva_limma_visual[(low0 + 1):high0,],aes(x = id,y = -0.1,label = id), hjust = 1,color = 'grey') + # 灰色标签
  geom_text(data = GSE20163_gsva_limma_visual[(high0 +1):high1,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'black') # 大于1的为黑色标签
pdf(paste0(save_path,"GSE20163_gsva_barplot.pdf"), width = 13.5, height = 6)
GSE20163_gsva_limma_picture
dev.off()

# Venn图
library(VennDiagram) 
library(venn)

# cellular senescence-related gene (CSG)
library(ggvenn)
library(tidyverse)
library(ggtext)
Candidate_CSG = list("CellAge" = CSRG_CellAge, "Cell Senescence Database" = CSRG_CSGene, "GeneCards" = CSRG_GeneCard)
CSG_All = ggvenn(data = Candidate_CSG, 
                 show_percentage = FALSE, 
                 show_elements = FALSE, 
                 label_sep = ",", 
                 digits = 1,
                 stroke_color = "white",
                 text_size = 7,
                 fill_color = c("#C9342B","#FAAC90", "#339DB5"), set_name_color = c("#C9342B","#FAAC90", "#339DB5"))
CSG_All

pdf(paste0(save_path,"cellular_senescence_related_genes.pdf"), width = 8, height = 6)
CSG_All
dev.off()

# cellular senescence-related differentially expressed gene (CSDEG)
Candidate_CSDEG = list("GSE8397" = GSE8397_DEG$gene, "GSE20163" = GSE20163_DEG$gene, "CSRG" = CSRG_total)
CSDEG_All = ggvenn(data = Candidate_CSDEG, 
                   show_percentage = FALSE, 
                   show_elements = FALSE, 
                   label_sep = ",", 
                   digits = 1,
                   stroke_color = "white",
                   text_size = 7,
                   fill_color = c("#C9342B","#FAAC90", "#339DB5"), set_name_color = c("#C9342B","#FAAC90", "#339DB5"))
CSDEG_All

pdf(paste0(save_path,"cellular_senescence_related_CSDEG.pdf"), width = 8, height = 6)
CSDEG_All
dev.off()

# 前10个CSDEG（上调/下调）的热图
gene_set_intersection = function(DEG,GeneSet,type){
  if (type == "up"){
    target_DEG = DEG[DEG$change == "up", ]
    target_DEG = rownames(target_DEG)
    up_inter_DEG = intersect(target_DEG,GeneSet)
    return(up_inter_DEG)
  }
  if (type == "down"){
    target_DEG = DEG[DEG$change == "down", ]
    target_DEG = rownames(target_DEG)
    down_inter_DEG = intersect(target_DEG,GeneSet)
    return(down_inter_DEG)
  }
  else {
    print("Either parameter or option is illegal.")
  }
}

# GSE8397
up_CSDEG_GSE8397 = gene_set_intersection(GSE8397_DEG,CSRG_total,"up")
down_CSDEG_GSE8397 = gene_set_intersection(GSE8397_DEG,CSRG_total,"down")
CSDEG_GSE8397 = c(up_CSDEG_GSE8397,down_CSDEG_GSE8397)

up_CSDEG_GSE8397_ordered = GSE8397_volcano_data[GSE8397_volcano_data$change == "up", ]
up_CSDEG_GSE8397_ordered = up_CSDEG_GSE8397_ordered[order(up_CSDEG_GSE8397_ordered$logFC,decreasing = T),]

down_CSDEG_GSE8397_ordered = GSE8397_volcano_data[GSE8397_volcano_data$change == "down", ]
down_CSDEG_GSE8397_ordered = down_CSDEG_GSE8397_ordered[order(down_CSDEG_GSE8397_ordered$logFC),]

GSE8397_CSDEG_top10 = rbind(up_CSDEG_GSE8397_ordered[1:10,],down_CSDEG_GSE8397_ordered[1:10,])

GSE8397_heatmap_CSDEG_top10 <- GSE8397_heatmap_data_log[GSE8397_CSDEG_top10$gene,]

annotation_col = GSE8397_group
rownames(annotation_col) = annotation_col$sample
annotation_col$sample = NULL

library(pheatmap)
heatmap_color = colorRampPalette(c("navy", "white", "firebrick3"))(100)
ann_colors = list(group = c(control = "#339DB5", PD = "#C9342B"))
GSE8397_CSDEG_top10_heatmap <- pheatmap::pheatmap(GSE8397_heatmap_CSDEG_top10,
                                            show_colnames = F,
                                            show_rownames = T,
                                            scale = "row",
                                            cluster_cols = F,
                                            annotation_col = annotation_col,
                                            breaks = seq(-3, 3, length.out = 100),
                                            annotation_legend = TRUE,
                                            border_color = NA,
                                            annotation_colors = ann_colors,
                                            fontsize_row = 6.5,
                                            main = "GSE8397 & CSRG")
GSE8397_CSDEG_top10_heatmap
#保存结果
pdf(paste0(save_path,"GSE8397_CSDEG_top_10_heatmap.pdf"), width=6, height=8)
GSE8397_CSDEG_top10_heatmap
dev.off()

# GSE20163
up_CSDEG_GSE20163 = gene_set_intersection(GSE20163_DEG,CSRG_total,"up")
down_CSDEG_GSE20163 = gene_set_intersection(GSE20163_DEG,CSRG_total,"down")

CSDEG_GSE20163 = c(up_CSDEG_GSE20163,down_CSDEG_GSE20163)

up_CSDEG_GSE20163_ordered = GSE20163_volcano_data[GSE20163_volcano_data$change == "up", ]
up_CSDEG_GSE20163_ordered = up_CSDEG_GSE20163_ordered[order(up_CSDEG_GSE20163_ordered$logFC,decreasing = T),]

down_CSDEG_GSE20163_ordered = GSE20163_volcano_data[GSE20163_volcano_data$change == "down", ]
down_CSDEG_GSE20163_ordered = down_CSDEG_GSE20163_ordered[order(down_CSDEG_GSE20163_ordered$logFC),]

GSE20163_CSDEG_top10 = rbind(up_CSDEG_GSE20163_ordered[1:10,],down_CSDEG_GSE20163_ordered[1:10,])

GSE20163_heatmap_CSDEG_top10 <- GSE20163_heatmap_data_log[GSE20163_CSDEG_top10$gene,]

annotation_col = GSE20163_group
rownames(annotation_col) = annotation_col$sample
annotation_col$sample = NULL

library(pheatmap)
heatmap_color = colorRampPalette(c("navy", "white", "firebrick3"))(100)
ann_colors = list(group = c(control = "#339DB5", PD = "#C9342B"))
GSE20163_CSDEG_top10_heatmap <- pheatmap::pheatmap(GSE20163_heatmap_CSDEG_top10,
                                                   show_colnames = F,
                                                   show_rownames = T,
                                                   scale = "row",
                                                   cluster_cols = F,
                                                   annotation_col = annotation_col,
                                                   breaks = seq(-3, 3, length.out = 100),
                                                   annotation_legend = TRUE,
                                                   border_color = NA,
                                                   annotation_colors = ann_colors,
                                                   fontsize_row = 6.5,
                                                   main = "GSE20163 & CSRG")
GSE20163_CSDEG_top10_heatmap
#保存结果
pdf(paste0(save_path,"GSE20163_CSDEG_top_10_heatmap.pdf"), width=6, height=8)
GSE20163_CSDEG_top10_heatmap
dev.off()

# write.csv(CSDEG,paste0(save_path,"CSDEG.csv"))


# 机器学习（LASSO + SVM）
##设置随机数种子以求结果可重复
set.seed(76523)

# GSE8397
# 准备数据
GSE8397_select_ml_expr = GSE8397[CSDEG,]
GSE8397_select_ml_expr = t(GSE8397_select_ml_expr)

GSE8397_select_ml_group = GSE8397_group
rownames(GSE8397_select_ml_group) = GSE8397_select_ml_group$sample
GSE8397_select_ml_group$sample = NULL
GSE8397_select_ml_group = as.matrix(GSE8397_select_ml_group)

# A. LASSO法
library(glmnet)
#通过glmnet函数中的设置family参数定义采用的算法模型（binomial为二分类logistics回归，cox包含生存分析）
GSE8397_mod <- glmnet(x = GSE8397_select_ml_expr,y = GSE8397_select_ml_group,family = "binomial") 
GSE8397_mod
##Lasso回归最重要的就是选择合适的λ值，可以通过cv.glmnet函数实现
#交叉验证
GSE8397_cvmod <- cv.glmnet(x = GSE8397_select_ml_expr,y = GSE8397_select_ml_group,family = "binomial")
GSE8397_cvmod
#保存结果

pdf(paste0(save_path,"GSE8397_lasso_mod.pdf"), width = 7, height = 5)
plot(GSE8397_mod,label = T,lwd=2)
dev.off()

pdf(paste0(save_path,"GSE8397_lasso_lambda.pdf"), width = 7, height = 5)
plot(GSE8397_mod,xvar = "lambda",label = T,lwd=2)
dev.off()

pdf(paste0(save_path,"GSE8397_lasso_cvmod.pdf"), width = 7, height = 7)
plot(GSE8397_cvmod)
dev.off()
#基于该图选择最佳的λ，一般可以通过cvmod$lambda.min和cvmod$lambda.1se实现
GSE8397_cvmod$lambda.min
#基因筛选，采用coef函数，有相应参数的gene则被保留，采用λ使用的是lambda.min
GSE8397_coef.min <- coef(GSE8397_cvmod,s="lambda.min")
GSE8397_coef.min #第二列有数值是非点号的则代表被选择的基因

#汇总结果
GSE8397_gene_lasso = c("HMGCR","HSPB1","ITPR1","SDC1","CDKN2C")

# SVM
library(e1071)
library(caret)
source('msvmRFE.R')

#SVM-REF算法输入数据
GSE8397_select_ml_expr = as.data.frame(GSE8397_select_ml_expr)
GSE8397_select_ml_expr$sample = rownames(GSE8397_select_ml_expr)

GSE8397_select_ml_group = as.data.frame(GSE8397_select_ml_group)
GSE8397_select_ml_group$sample = rownames(GSE8397_select_ml_group)

GSE8397_select_svm_ml_ci = merge(GSE8397_select_ml_expr,GSE8397_select_ml_group,by.x = "sample", by.y = "sample")
rownames(GSE8397_select_svm_ml_ci) = GSE8397_select_svm_ml_ci$sample
GSE8397_select_svm_ml_ci$sample = NULL

GSE8397_select_svm_ml_ci = dplyr::select(GSE8397_select_svm_ml_ci,"group",everything())

GSE8397_svm_input <- GSE8397_select_svm_ml_ci
GSE8397_svm_input$group = as.factor(GSE8397_svm_input$group)

#采用十折交叉验证 (k-fold crossValidation）
svmRFE(GSE8397_svm_input, k = 10, halve.above = 100) #分割数据，分配随机数
nfold = 10
nrows = nrow(GSE8397_svm_input)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))
GSE8397_results = lapply(folds, svmRFE.wrap, GSE8397_svm_input, k=10, halve.above=100) #特征选择
top.features = WriteFeatures(GSE8397_results, GSE8397_svm_input, save=F) #查看主要变量
head(top.features)

library(parallel)

# Calculate the number of cores
no_cores <- detectCores() - 1

# Initiate cluster
cl <- makeCluster(no_cores)

# Estimate generalization error using a varying number of top features
clusterExport(cl, list("top.features","GSE8397_results", "tune","tune.control","svm"))

start.time <- Sys.time()
start.time

GSE8397_featsweep = parLapply(cl,1:length(CSDEG), FeatSweep.wrap, GSE8397_results, GSE8397_svm_input)

stopCluster(cl)
end.time <- Sys.time()
GSE8397_time.taken <- end.time - start.time
GSE8397_time.taken

# 画图
GSE8397_no.info = min(prop.table(table(GSE8397_svm_input[,1])))
GSE8397_errors = sapply(GSE8397_featsweep, function(x) ifelse(is.null(x), NA, x$error))

#绘制基于SVM-REF算法的错误率曲线图
PlotErrors(GSE8397_errors, no.info=GSE8397_no.info) #查看错误率
pdf(paste0(save_path,"GSE8397_svm-error.pdf"),width = 5,height = 5)
PlotErrors(GSE8397_errors, no.info=GSE8397_no.info) #查看错误率
dev.off()

#Plotaccuracy函数（SVM-RFE有关的代码没有，此处特地补充）
PlotAccuracy <- function(errors, errors2=NULL, no.info=0.5, ylim=range(c(errors, errors2), na.rm=T), xlab='Number of Features',  ylab='5x CV Accuracy') {
  # Makes a plot of average generalization error vs. number of top features
  AddLine <- function(x, col='black') {
    lines(which(!is.na(errors)), na.omit(x), col=col)
    points(which.max(x), max(x, na.rm=T), col='red')
    text(which.max(x), max(x, na.rm=T), paste(which.max(x), '-', format(max(x, na.rm=T), dig=3)), pos=4, col='red', cex=0.75)
  }
  plot(errors, type='n', ylim=ylim, xlab=xlab, ylab=ylab)
  AddLine(errors)
  if(!is.null(errors2)) AddLine(errors2, 'gray30')
  abline(h=no.info, lty=3)
}

#绘制基于SVM-REF算法的正确率曲线图
PlotAccuracy(1-GSE8397_errors,no.info=GSE8397_no.info) #查看准确率
pdf(paste0(save_path,"GSE8397_svm-accuracy.pdf"),width = 5,height = 5)
PlotAccuracy(1-GSE8397_errors,no.info=GSE8397_no.info) #查看准确率
dev.off()

#根据图示裁定阈值为n（取前n个基因）
GSE8397_svm_n = 14
GSE8397_gene_svm = top.features[1:GSE8397_svm_n,1]

# 随机森林法
GSE8397_rf_input = as.data.frame(t(GSE8397))
GSE8397_rf_input = dplyr::select(GSE8397_rf_input,all_of(CSDEG))

identical(rownames(GSE8397_rf_input),GSE8397_group$sample)
GSE8397_rf_input = as.matrix(GSE8397_rf_input)

GSE8397_rf_group = factor(GSE8397_group$group)

library(randomForest)

# 直接使用默认参数
GSE8397_rf <- randomForest(GSE8397_rf_input, GSE8397_rf_group)
#查看随机森林的模型
#查看变量的重要性
GSE8397_importance <- importance(x=GSE8397_rf)
#可视化特征的重要性并保存
plot(GSE8397_rf,main = "Random Forest")
pdf(paste0(save_path,"GSE8397_rf_trees.pdf"), width = 6, height = 6)
plot(GSE8397_rf,main = "Random Forest")
dev.off()

varImpPlot(GSE8397_rf,main = "Random Forest")
pdf(paste0(save_path,"GSE8397_rf_MeanDecreaseGini.pdf"), width = 6, height = 6)
varImpPlot(GSE8397_rf,main = "Random Forest")
dev.off()

#汇总结果（MeanDecreaseGini>2）
GSE8397_importance_result = as.data.frame(GSE8397_importance)
GSE8397_gene_rf = names(GSE8397_importance[GSE8397_importance_result$MeanDecreaseGini >= 2,])

# GSE20163
# 准备数据
GSE20163_select_ml_expr = GSE20163[CSDEG,]
GSE20163_select_ml_expr = t(GSE20163_select_ml_expr)

GSE20163_select_ml_group = GSE20163_group
rownames(GSE20163_select_ml_group) = GSE20163_select_ml_group$sample
GSE20163_select_ml_group$sample = NULL
GSE20163_select_ml_group = as.matrix(GSE20163_select_ml_group)

# A. LASSO法
library(glmnet)
#通过glmnet函数中的设置family参数定义采用的算法模型（binomial为二分类logistics回归，cox包含生存分析）
GSE20163_mod <- glmnet(x = GSE20163_select_ml_expr,y = GSE20163_select_ml_group,family = "binomial") 
GSE20163_mod
##Lasso回归最重要的就是选择合适的λ值，可以通过cv.glmnet函数实现
#交叉验证
GSE20163_cvmod <- cv.glmnet(x = GSE20163_select_ml_expr,y = GSE20163_select_ml_group,family = "binomial")
GSE20163_cvmod
#保存结果

pdf(paste0(save_path,"GSE20163_lasso_mod.pdf"), width = 7, height = 5)
plot(GSE20163_mod,label = T,lwd=2)
dev.off()

pdf(paste0(save_path,"GSE20163_lasso_lambda.pdf"), width = 7, height = 5)
plot(GSE20163_mod,xvar = "lambda",label = T,lwd=2)
dev.off()

pdf(paste0(save_path,"GSE20163_lasso_cvmod.pdf"), width = 7, height = 7)
plot(GSE20163_cvmod)
dev.off()
#基于该图选择最佳的λ，一般可以通过cvmod$lambda.min和cvmod$lambda.1se实现
GSE20163_cvmod$lambda.min
#基因筛选，采用coef函数，有相应参数的gene则被保留，采用λ使用的是lambda.min
GSE20163_coef.min <- coef(GSE20163_cvmod,s="lambda.min")
GSE20163_coef.min #第二列有数值是非点号的则代表被选择的基因

#汇总结果
GSE20163_gene_lasso = c("CDK5R1","HMGCR","HSPB1","YWHAZ")

# SVM
library(e1071)
library(caret)
source('msvmRFE.R')

#SVM-REF算法输入数据
GSE20163_select_ml_expr = as.data.frame(GSE20163_select_ml_expr)
GSE20163_select_ml_expr$sample = rownames(GSE20163_select_ml_expr)

GSE20163_select_ml_group = as.data.frame(GSE20163_select_ml_group)
GSE20163_select_ml_group$sample = rownames(GSE20163_select_ml_group)

GSE20163_select_svm_ml_ci = merge(GSE20163_select_ml_expr,GSE20163_select_ml_group,by.x = "sample", by.y = "sample")
rownames(GSE20163_select_svm_ml_ci) = GSE20163_select_svm_ml_ci$sample
GSE20163_select_svm_ml_ci$sample = NULL

GSE20163_select_svm_ml_ci = dplyr::select(GSE20163_select_svm_ml_ci,"group",everything())

GSE20163_svm_input <- GSE20163_select_svm_ml_ci
GSE20163_svm_input$group = as.factor(GSE20163_svm_input$group)

#采用十折交叉验证 (k-fold crossValidation）
svmRFE(GSE20163_svm_input, k = 10, halve.above = 100) #分割数据，分配随机数
nfold = 10
nrows = nrow(GSE20163_svm_input)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))
GSE20163_results = lapply(folds, svmRFE.wrap, GSE20163_svm_input, k=10, halve.above=100) #特征选择
top.features = WriteFeatures(GSE20163_results, GSE20163_svm_input, save=F) #查看主要变量
head(top.features)

library(parallel)

# Calculate the number of cores
no_cores <- detectCores() - 1

# Initiate cluster
cl <- makeCluster(no_cores)

# Estimate generalization error using a varying number of top features
clusterExport(cl, list("top.features","GSE20163_results", "tune","tune.control","svm"))

start.time <- Sys.time()
start.time

GSE20163_featsweep = parLapply(cl,1:length(CSDEG), FeatSweep.wrap, GSE20163_results, GSE20163_svm_input)

stopCluster(cl)
end.time <- Sys.time()
GSE20163_time.taken <- end.time - start.time
GSE20163_time.taken

# 画图
GSE20163_no.info = min(prop.table(table(GSE20163_svm_input[,1])))
GSE20163_errors = sapply(GSE20163_featsweep, function(x) ifelse(is.null(x), NA, x$error))

#绘制基于SVM-REF算法的错误率曲线图
PlotErrors(GSE20163_errors, no.info=GSE20163_no.info) #查看错误率
pdf(paste0(save_path,"GSE20163_svm-error.pdf"),width = 5,height = 5)
PlotErrors(GSE20163_errors, no.info=GSE20163_no.info) #查看错误率
dev.off()

#Plotaccuracy函数（SVM-RFE有关的代码没有，此处特地补充）
PlotAccuracy <- function(errors, errors2=NULL, no.info=0.5, ylim=range(c(errors, errors2), na.rm=T), xlab='Number of Features',  ylab='5x CV Accuracy') {
  # Makes a plot of average generalization error vs. number of top features
  AddLine <- function(x, col='black') {
    lines(which(!is.na(errors)), na.omit(x), col=col)
    points(which.max(x), max(x, na.rm=T), col='red')
    text(which.max(x), max(x, na.rm=T), paste(which.max(x), '-', format(max(x, na.rm=T), dig=3)), pos=4, col='red', cex=0.75)
  }
  plot(errors, type='n', ylim=ylim, xlab=xlab, ylab=ylab)
  AddLine(errors)
  if(!is.null(errors2)) AddLine(errors2, 'gray30')
  abline(h=no.info, lty=3)
}

#绘制基于SVM-REF算法的正确率曲线图
PlotAccuracy(1-GSE20163_errors,no.info=GSE20163_no.info) #查看准确率
pdf(paste0(save_path,"GSE20163_svm-accuracy.pdf"),width = 5,height = 5)
PlotAccuracy(1-GSE20163_errors,no.info=GSE20163_no.info) #查看准确率
dev.off()

#根据图示裁定阈值为n（取前n个基因）
GSE20163_svm_n = 9
GSE20163_gene_svm = top.features[1:GSE20163_svm_n,1]

# 随机森林法
GSE20163_rf_input = as.data.frame(t(GSE20163))
GSE20163_rf_input = dplyr::select(GSE20163_rf_input,all_of(CSDEG))

identical(rownames(GSE20163_rf_input),GSE20163_group$sample)
GSE20163_rf_input = as.matrix(GSE20163_rf_input)

GSE20163_rf_group = factor(GSE20163_group$group)

library(randomForest)

# 直接使用默认参数
GSE20163_rf <- randomForest(GSE20163_rf_input, GSE20163_rf_group)
#查看随机森林的模型
#查看变量的重要性
GSE20163_importance <- importance(x=GSE20163_rf)
#可视化特征的重要性并保存
plot(GSE20163_rf,main = "Random Forest")
pdf(paste0(save_path,"GSE20163_rf_trees.pdf"), width = 6, height = 6)
plot(GSE20163_rf,main = "Random Forest")
dev.off()

varImpPlot(GSE20163_rf,main = "Random Forest")
pdf(paste0(save_path,"GSE20163_rf_MeanDecreaseGini.pdf"), width = 6, height = 6)
varImpPlot(GSE20163_rf,main = "Random Forest")
dev.off()

#汇总结果（MeanDecreaseGini>0.6）
GSE20163_importance_result = as.data.frame(GSE20163_importance)
GSE20163_gene_rf = names(GSE20163_importance[GSE20163_importance_result$MeanDecreaseGini >= 0.6,])

# 获取hubCSDEG
GSE8397_hubCSDEG = Reduce(intersect,list(GSE8397_gene_lasso,GSE8397_gene_svm,GSE8397_gene_rf))

GSE20163_hubCSDEG = Reduce(intersect,list(GSE20163_gene_lasso,GSE20163_gene_svm,GSE20163_gene_rf))

hubCSDEG = union(GSE8397_hubCSDEG,GSE20163_hubCSDEG)

# 可视化
# GSE8397
library(ggvenn)
library(tidyverse)
library(ggtext)
GSE8397_CSDEG = list("LASSO" = GSE8397_gene_lasso, "SVM-RFE" = GSE8397_gene_svm, "Random Forest" = GSE8397_gene_rf)
GSE8397_CSDEG_venn = ggvenn(data = GSE8397_CSDEG, 
                            show_percentage = FALSE, 
                            show_elements = FALSE, 
                            label_sep = ",", 
                            digits = 1,
                            stroke_color = "white",
                            text_size = 7,
                            fill_color = c("#C9342B","#FAAC90", "#339DB5"), set_name_color = c("#C9342B","#FAAC90", "#339DB5"))
GSE8397_CSDEG_venn

pdf(paste0(save_path,"GSE8397_CSDEG_ml_venn.pdf"), width = 8, height = 6)
GSE8397_CSDEG_venn
dev.off()

# GSE20163
library(ggvenn)
library(tidyverse)
library(ggtext)
GSE20163_CSDEG = list("LASSO" = GSE20163_gene_lasso, "SVM-RFE" = GSE20163_gene_svm, "Random Forest" = GSE20163_gene_rf)
GSE20163_CSDEG_venn = ggvenn(data = GSE20163_CSDEG, 
                             show_percentage = FALSE, 
                             show_elements = FALSE, 
                             label_sep = ",", 
                             digits = 1,
                             stroke_color = "white",
                             text_size = 7,
                             fill_color = c("#C9342B","#FAAC90", "#339DB5"), set_name_color = c("#C9342B","#FAAC90", "#339DB5"))
GSE20163_CSDEG_venn

pdf(paste0(save_path,"GSE20163_CSDEG_ml_venn.pdf"), width = 8, height = 6)
GSE20163_CSDEG_venn
dev.off()

# hubCSDEG
hub_CSDEG_list = list("hub-CSDEG in GSE8397" = GSE8397_hubCSDEG, "hub-CSDEG in GSE20163" = GSE20163_hubCSDEG)
hub_CSDEG_venn = ggvenn(data = hub_CSDEG_list, 
                        show_percentage = FALSE, 
                        show_elements = F, 
                        label_sep = ",", 
                        digits = 1,
                        stroke_color = "white",
                        text_size = 7,
                        fill_color = c("#C9342B","#339DB5"), set_name_color = c("#C9342B","#339DB5"))
hub_CSDEG_venn

pdf(paste0(save_path,"hub_CSDEG_ml_venn.pdf"), width = 8, height = 6)
hub_CSDEG_venn
dev.off()

# 免疫浸润分析（合并执行）
# 对GSE8397和GSE20163合并去批次
library(FactoMineR)
library(factoextra)
library(sva)
library(dplyr)

sva_GSE8397 = GSE8397
sva_GSE8397$gene = rownames(sva_GSE8397)

sva_GSE20163 = GSE20163
sva_GSE20163$gene = rownames(sva_GSE20163)


Exp_entire = merge(sva_GSE8397,sva_GSE20163,by.x = "gene", by.y = "gene", all = TRUE)
#删除含NA值的行
Exp_entire = na.omit(Exp_entire)
#基因名作行名
rownames(Exp_entire) = Exp_entire$gene
Exp_entire$gene = NULL
#准备分组信息
sva_GSE8397_group = GSE8397_group
sva_GSE8397_group$batch = "GPL5175"

sva_GSE20163_group = GSE20163_group
sva_GSE20163_group$batch = "GPL1352"

Exp_entire_group = rbind(sva_GSE8397_group,sva_GSE20163_group)
Exp_entire_group_control = Exp_entire_group[Exp_entire_group$group == "control",]
Exp_entire_group_PD = Exp_entire_group[Exp_entire_group$group == "PD",]
Exp_entire_group = rbind(Exp_entire_group_control,Exp_entire_group_PD)
rm(Exp_entire_group_control,Exp_entire_group_PD)

identical(Exp_entire_group$sample,colnames(Exp_entire))


Exp_entire = select(Exp_entire,all_of(Exp_entire_group$sample))

identical(Exp_entire_group$sample,colnames(Exp_entire))

model <- model.matrix(~as.factor(Exp_entire_group$group))


combat_Expr <- ComBat(dat = Exp_entire,batch = Exp_entire_group$batch, mod = model)
combat_Expr = as.data.frame(combat_Expr)
identical(Exp_entire_group$sample,colnames(Exp_entire))
Merged_Matrix = as.data.frame(t(combat_Expr))
Merged_Matrix$sample = rownames(Merged_Matrix)
Merged_Matrix = dplyr::select(Merged_Matrix,"sample",everything())
Merged_Matrix = merge(Merged_Matrix,Exp_entire_group,by.x = "sample",by.y = "sample")
Merged_Matrix$batch = NULL
Merged_Matrix = dplyr::select(Merged_Matrix,c("sample","group"),everything())
rownames(Merged_Matrix) = Merged_Matrix$sample
Merged_Matrix$sample = NULL
Merged_Matrix = as.data.frame(t(Merged_Matrix))

# write.table(Merged_Matrix,file = paste0(save_path,"Merged_Matrix.txt"), quote = FALSE, sep = "  ",)


# ssGSEA分析
library(ggplot2)
library(tinyarray)
library(GSVA)
library(dplyr)
library(Hmisc)
library(pheatmap)
library(ggpubr)

# 加载28种免疫细胞的cell marker
temp = read.csv("cell_marker.csv")
# 提取细胞种类
cell_type = unique(temp[,2])
cell_marker = c()
# 转换各细胞种类的全部marker为列表
for (cell in cell_type) {
  cell_name = cell
  all_markers = temp[temp$Cell.type == cell_name,]
  all_markers = lapply(all_markers, as.list)
  all_markers[2] = NULL
  names(all_markers) = cell_name
  cell_marker = append(cell_marker, all_markers)
  rm(cell_name,all_markers)
}
rm(temp,cell_type,cell)

library(genefilter)
library(GSVA)
library(Biobase)
library(stringr)


# 基于ssGSEA计算免疫浸润分数
im_ssgsea = ssgseaParam(exprData=as.matrix(combat_Expr), geneSets=cell_marker)
im_ssgsea = gsva(im_ssgsea,verbose = TRUE, BPPARAM = MulticoreParam(workers = 30, progressbar = TRUE))
gc()
im_ssgsea_backup = im_ssgsea
# 0-1标准化（可视化的通例）
for (coln in colnames(im_ssgsea)) {
  im_ssgsea[,coln] <- (im_ssgsea[,coln] -min(im_ssgsea[,coln]))/(max(im_ssgsea[,coln] )-min(im_ssgsea[,coln] ))
}
im_ssgsea_backup = im_ssgsea
# 针对全体ssGSEA结果可视化
all_im_ssgsea = im_ssgsea
all_im_ssgsea = as.data.frame(t(all_im_ssgsea))
all_im_ssgsea$sample = rownames(all_im_ssgsea)
all_im_ssgsea = merge(all_im_ssgsea,Exp_entire_group,by.x = "sample", by.y = "sample")
all_im_ssgsea = all_im_ssgsea %>% tidyr::pivot_longer(cols = c(2:29),names_to = "cell_type",values_to = "fraction")
colnames(all_im_ssgsea)[2] = "Group"

ssGSEA_output = as.data.frame(im_ssgsea)
ssGSEA_output$immunocyte = rownames(ssGSEA_output)
ssGSEA_output = dplyr::select(ssGSEA_output,"immunocyte",everything())
# openxlsx::write.xlsx(ssGSEA_output,paste0(save_path,"Table S20.xlsx"))

# a.箱线图
library(forcats)
library(gghalves)
all_im_ssgsea$cell_type <- factor(all_im_ssgsea$cell_type,levels = c("Activated B cell","Mast cell","Natural killer cell","Monocyte", "Natural killer T cell","Activated CD4 T cell","Activated CD8 T cell","Activated dendritic cell","CD56bright natural killer cell","CD56dim natural killer cell","Central memory CD4 T cell","Central memory CD8 T cell","Effector memeory CD4 T cell","Effector memeory CD8 T cell","Eosinophil","Gamma delta T cell","Immature  B cell", "Immature dendritic cell", "Macrophage", "MDSC", "Memory B cell", "Neutrophil", "Plasmacytoid dendritic cell", "Regulatory T cell", "T follicular helper cell", "Type 1 T helper cell", "Type 17 T helper cell", "Type 2 T helper cell"))
ssGSEA_boxplot_grouped <- all_im_ssgsea %>%
  ggplot(aes(cell_type,fraction)) +
  geom_half_boxplot(aes(fill=Group),color="black",
                    side="l",errorbar.draw = T, 
                    outlier.shape = NA, width=0.8) +            
  geom_half_point(aes(fill=Group),
                  side = "r", size = 0.1)+
  scale_fill_manual(values = c("#339DB5","#C9342B")) +
  theme_bw() + 
  labs(x = NULL, y = "Estimated Expression") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=60,hjust = 1, size = 12),
        axis.text = element_text(color = "black",size = 12),
        axis.text.y = element_text(size = 12))+
  ggpubr::stat_compare_means(aes(group = Group,label = ..p.signif..), method = "t.test",label.y = 1.005, show.legend = TRUE) +
  theme(panel.grid=element_blank())+
  coord_cartesian()
ssGSEA_boxplot_grouped
pdf(paste0(save_path,"ssGSEA_boxplot_grouped.pdf"), width = 9, height = 6)
ssGSEA_boxplot_grouped
dev.off()

# b.相关性热图
library(corrplot)

cor_im_ssGSEA_input = as.data.frame(t(im_ssgsea))


#使用ggcorrplot包的cor_pmat()函数计算p值：
cor_im_ssGSEA <- round(cor(cor_im_ssGSEA_input), 3)
pvalue_cor_im_ssGSEA <- ggcorrplot::cor_pmat(cor_im_ssGSEA_input)


mycol <- colorRampPalette(c("#339DB5", "white", "#C9342B"), alpha = TRUE)

pdf(paste0(save_path,"ssGSEA_cor_heatmap.pdf"), width = 15.5, height = 10)
corrplot(cor_im_ssGSEA, method = c('pie'), 
         type = c('upper'), 
         col = mycol(100),
         outline = 'grey', 
         order = c('AOE'), 
         diag = TRUE,
         tl.cex = 1.2,
         cl.cex = 1,
         tl.col = 'black',
         tl.pos = 'tp',
         tl.srt = 90,
         p.mat = pvalue_cor_im_ssGSEA,
         sig.level = c(.001, .01, .05),
         insig = "label_sig", #显著性标注样式："pch", "p-value", "blank", "n", "label_sig"
         pch.cex = 1.2, #显著性标记大小
         pch.col = 'black' #显著性标记颜色
)
#下三角图添加不显著叉号：
corrplot(cor_im_ssGSEA, add = TRUE,
         method = c('number'), 
         type = c('lower'),
         col = mycol(100),
         order = c('AOE'), 
         diag = FALSE, 
         number.cex = 0.7,
         tl.pos = 'n', 
         cl.pos = 'n',
         p.mat = pvalue_cor_im_ssGSEA,
         insig = "pch")
dev.off()

# c.基因-免疫细胞相关性热图
re = im_ssgsea

exp = combat_Expr

library(Hmisc)
identical(colnames(re),colnames(exp))

nc = t(rbind(re,exp[hubCSDEG,]))

m = rcorr(nc)$r[1:nrow(re),(ncol(nc)-length(hubCSDEG)+1):ncol(nc)]

p = rcorr(nc)$P[1:nrow(re),(ncol(nc)-length(hubCSDEG)+1):ncol(nc)]

library(dplyr)

ssGSEA_m_output = as.data.frame(m)
ssGSEA_m_output$immunocyte = rownames(ssGSEA_m_output)
ssGSEA_m_output = dplyr::select(ssGSEA_m_output,"immunocyte",everything())
# openxlsx::write.xlsx(ssGSEA_m_output,paste0(save_path,"Table S22A.xlsx"))


ssGSEA_p_output = as.data.frame(p)
ssGSEA_p_output$immunocyte = rownames(ssGSEA_p_output)
ssGSEA_p_output = dplyr::select(ssGSEA_p_output,"immunocyte",everything())
# openxlsx::write.xlsx(ssGSEA_p_output,paste0(save_path,"Table S22B.xlsx"))

tmp=ifelse(p>0.1,"",ifelse(p>0.05,"",ifelse(p>0.01,"*",ifelse(p>0.001,"**","***"))))

source("modified_pheatmap.R")

pdf(paste0(save_path,"ssGSEA_correlation_hubCSDEG_immuocytes.pdf"), width = 8, height = 6)
pheatmap(m,
         display_numbers =tmp,
         angle_col =45,
         color = colorRampPalette(c("#339DB5", "white", "#C9342B"))(100),
         border_color = "white",
         treeheight_col = 0,
         treeheight_row = 0)
dev.off()

# d. 堆叠柱状图
stackplot_data = as.data.frame(t(im_ssgsea))
identical(rownames(stackplot_data),Exp_entire_group$sample)
stackplot_data$sample = rownames(stackplot_data)
stackplot_data = merge(stackplot_data,Exp_entire_group, by.x = "sample", by.y = "sample")
stackplot_data$batch = NULL
stackplot_data = dplyr::select(stackplot_data,c("sample","group"),everything())
stackplot_data_part_1 = stackplot_data[stackplot_data$group == "control", ]
stackplot_data_part_2 = stackplot_data[stackplot_data$group == "PD", ]
stackplot_data = rbind(stackplot_data_part_1, stackplot_data_part_2)
rm(stackplot_data_part_1, stackplot_data_part_2)
library(dplyr)
library(tidyr)
stackplot_data = stackplot_data %>% 
  pivot_longer(c(3:30),names_to = "cell_type",values_to = "fraction")
## 把ssGSEA的数值转换为比例
library(tidyverse)
stackplot_input = stackplot_data %>% 
  group_by(sample) %>% 
  mutate(ratio=fraction/sum(fraction))
## 绘制堆叠柱状图
library(ggplot2)
library(wesanderson)

stack_color = wes_palette("Zissou1", 28, type = "continuous")
imssgsea_stackplot <- stackplot_input %>% 
  ggplot(aes(sample,ratio))+
  geom_bar(stat = "identity",position = "stack",aes(fill=cell_type))+
  labs(x=NULL)+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = stack_color,name=NULL)+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom"
  )
#最后保存结果
pdf(paste0(save_path,"imssgsea_stackplot.pdf"), width =10, height = 8)
imssgsea_stackplot
dev.off()

# 单细胞部分
rm(list = ls());gc() #删除所有环境变量整理内存空间

hubCSDEG = c("HMGCR", "HSPB1", "ITPR1", "SDC1")
save_path = "/home/alan_turing/Projects/PD/ND_CL_EE/result/"

# brain
# GSE178265

# 1. R包安装与加载
if(!require(multtest))BiocManager::install("multtest")
if(!require(Seurat))install.packages("Seurat")
if(!require(dplyr))install.packages("dplyr")
if(!require(patchwork))install.packages("patchwork")
if(!require(R.utils))install.packages("R.utils")
# devtools::install_github('satijalab/seurat-data')
library(SeuratData)
# 
library(SeuratDisk)
# 注意：Seurat v5不支持future包！不要使用future包进行平行计算了！ 

# 2.加载数据

PD_brain.data <- Read10X(data.dir = "/home/alan_turing/Projects/PD/NM_PD/data/GSE178265/")
PD_brain <- CreateSeuratObject(counts = PD_brain.data, project = "PD_brain_snRNA")
objectSize(PD_brain)
rm(PD_brain.data)
gc()

# 加载metadata并删掉其中的无用内容
PD_brain_metadata <- read.table('/home/alan_turing/Projects/PD/NM_PD/data/GSE178265_METADATA_PD.tsv',
                                row.names = 1,
                                sep = '\t',
                                header = T)
objectSize(PD_brain_metadata)
PD_brain_metadata = PD_brain_metadata[-1,]
PD_brain_metadata = PD_brain_metadata[,c("libname","biosample_id","Status","organ__ontology_label")]
objectSize(PD_brain_metadata)

gc()

# 把metadata添加到Seurat对象之中（同时移除分组为LBD的细胞）
PD_brain <- AddMetaData(PD_brain,metadata = PD_brain_metadata)
# 去掉LBD样本与caudate nucleus的样本（即仅保留源于substantia nigra pars compacta的样本）
PD_targetsamples <- rownames(PD_brain@meta.data)[PD_brain@meta.data[["Status"]] != "LBD" & PD_brain@meta.data[["organ__ontology_label"]] != "caudate nucleus"]
PD_brain = subset(PD_brain, cells = PD_targetsamples)
rm(PD_targetsamples,PD_brain_metadata)
gc()

## 2.2 初步查看Seurat对象
PD_brain
ncol(PD_brain) # 评估细胞数目
# ncol(as.data.frame(PD_brain[["RNA"]]@counts)) # 评估细胞数目
# ncol(PD_brain.data)
# PD_brain@meta.data # 基本信息或质控数据

# 3.质控数据及可视化
# 计算线粒体RNA的含量
PD_brain[["percent.mt"]] <- PercentageFeatureSet(PD_brain, pattern = "^MT-")
# 计算线粒体RNA的比例
PD_brain$mitoRatio <- PD_brain[["percent.mt"]] / 100
gc()
# 人源的数据为MT，鼠源的需要换成mt

## 3.1 查看质控情况
head(PD_brain@meta.data,5)

## 3.2 筛选细胞（过滤掉不符合阈值的细胞）
# 原文筛掉了线粒体基因比例＜10%且单个细胞UMI计数（即nCount_RNA）＞650的细胞
PD_brain <- subset(PD_brain, subset = nCount_RNA > 650 & mitoRatio < 0.1)
gc()


# 4.细胞分群
## 4.1 表达数据标准化（转换定量的Count表达矩阵为基因表达量）

PD_brain <- NormalizeData(PD_brain,
                          normalization.method = "LogNormalize", 
                          scale.factor = 10000)
gc()
# scale.factor = 缩放因子
# 标准化值存放于PD_brain[["RNA"]]@data之中


## 4.2 寻找高变基因（为后续的PCA与UMAP图等分析做准备）
PD_brain <- FindVariableFeatures(PD_brain,
                                 selection.method = "vst",
                                 nfeatures = 2000)
gc()
# nfeatures = 高变基因的数量

## 4.4 scale（针对每一种基因在所有细胞中的标准化）
# 它可以减少由于某基因的表达量过于离群所导致的误差
# 基于高变基因进行scale（默认）
PD_brain <- ScaleData(PD_brain)
gc()

## 4.5 PCA 
PD_brain <- RunPCA(PD_brain)
gc()

# 由于本部分数据已经在NM_PD中使用，故不再另行计算维度
# ## 4.6 PCA结果可视化（确定后续选择的PCA维度）
# ### 展现出前30个维度的PCA结果（舍弃掉矩阵分离不显著的维度）
# pdf("DimHeatmap.pdf", width = 10, height =8)
# DimHeatmap(PD_brain, dims = 1:30, cells = 500, balanced = TRUE)
# dev.off()
# 
# ### 在标准误（Standard Deviation）下降至完全水平的拐点，选择相应的PCA维度
# pdf("ElbowPlot.pdf", width = 10, height =8)
# ElbowPlot(PD_brain)
# dev.off()

## 4.7 根据选取的PCA维度进行细胞分群
# dim = 10
PD_brain <- FindNeighbors(PD_brain, dims = 1:10)
# install.packages("clustree")
gc()

# 确定resolution的值
# 低分辨率（res小）：大群（如T细胞、B细胞等粗分类）
# 高分辨率（res大）： 亚群（如CD4+ T细胞细分亚型）
PD_brain <- FindClusters(PD_brain,
                         resolution = 0.05)
gc()

# 5.分群后的可视化（绘制降维图，降维并观察细胞分群的情况）
## 5.1 降维分群
### A. UMAP（适合大型数据，速度快）
PD_brain <- RunUMAP(PD_brain, dims = 1:10, n.neighbors = 60, min.dist = 0.7, spread = 0.6)

Idents(PD_brain)="RNA_snn_res.0.05"
PD_brain$RNA_snn_res.0.05 <-Idents(PD_brain)

# saveRDS(PD_brain,'ND_CL_EE_Annotated_GSE178265.rds')

library(ggsci)

mycol <- ggsci::pal_npg()(9)
mycol

pdf(paste0(save_path,"ND_CL_EE_umap_before_cell_annotation.pdf"), width = 10, height =8)
DimPlot(PD_brain, reduction = "umap", label = T, cols = mycol, raster = F)
dev.off()
gc()

## 5.2 寻找marker基因（此处细胞亚群 = cluster）
# 寻找所有分群的所有高表达marker基因
# install.packages('devtools')
# devtools::install_github('immunogenomics/presto')
# 提取安装好presto包可加速Seurat中的Wilcoxon秩和检验
# 如果看到推荐安装presto的提示，可以直接安装
PD_brain.markers <- FindAllMarkers(PD_brain,
                                   only.pos = TRUE, #只会去找特定上调的基因
                                   min.pct = 0.1, #基因至少在该cluster多少比例细胞中表达 
                                   logfc.threshold = 0.25 #倍数的对数的阈值，默认是0.25
)
gc()
## 获得的结果中，cluster为目标细胞分群
## 本质上是通过对比目标细胞分群（ident.1）和其他所有细胞分群（ident.2），寻找marker

### 提取细胞亚群的特异性marker
if(!require(dplyr))install.packages("dplyr")
# N = 50
top_N = PD_brain.markers %>% group_by(cluster) %>% top_n(n = 75, wt = avg_log2FC) # 提取各细胞亚群的前2个特异性marker

## 2-1. 可视化（热图，前10个标记基因）
# pdf("PD_cluster_top10_markers.pdf", width = 10, height =8)
# DoHeatmap(PD_brain, features = top_N$gene) + NoLegend()
# dev.off()
## 2-2. 可视化（小提琴图，前20个标记基因）
# VlnPlot(pbmc, features = top10$gene[1:20],pt.size=0)
## 可视化（小提琴图，前20个标记基因，显示marker表达量）
# VlnPlot(pbmc, features = top10$gene[1:20])

write.csv(top_N,paste0(save_path,"ND_CL_EE_top_N_markers.csv"))
write.csv(PD_brain.markers,paste0(save_path,"ND_CL_EE_PD_brain.markers.csv"))

# 手动重命名各细胞亚群（向量法）
GSE178265_celltype = c("Olig", "GABAn","Astro","MG","DAn","GLUn", "Endo","OPC","Astro")

names(GSE178265_celltype) <- levels(PD_brain)
PD_brain <- RenameIdents(PD_brain, GSE178265_celltype)
PD_brain$cell_type <- PD_brain@active.ident

Idents(PD_brain)="cell_type"

library(ggsci)
mycol <- ggsci::pal_aaas()(9)
mycol

pdf(paste0(save_path,"ND_CL_EE_umap_after_cell_annotation.pdf"), width = 10, height = 8)
DimPlot(PD_brain, reduction = "umap", label = T, cols = mycol, raster = F)
dev.off()

# saveRDS(PD_brain,'ND_CL_EE_Annotated_GSE178265.rds')
# PD_brain = readRDS('ND_CL_EE_Annotated_GSE178265.rds')

PD_brain@meta.data$Status = gsub("Ctrl","control",PD_brain@meta.data$Status)

pdf(paste0(save_path,"GSE178265_HMGCR_Feature_Plot.pdf"),width = 8, height = 4)
FeaturePlot(object = PD_brain, features =  "HMGCR", pt.size = 0.5,order = T,split.by = 'Status', cols = c("white", "#3C5488FF"), raster = F) & ggplot2::theme(legend.position = "right")
dev.off()

pdf(paste0(save_path,"GSE178265_HSPB1_Feature_Plot.pdf"),width = 8, height = 4)
FeaturePlot(object = PD_brain, features =  "HSPB1", pt.size = 0.5,order = T,split.by = 'Status', cols = c("white", "#3C5488FF"), raster = F) & ggplot2::theme(legend.position = "right")
dev.off()

pdf(paste0(save_path,"GSE178265_ITPR1_Feature_Plot.pdf"),width = 8, height = 4)
FeaturePlot(object = PD_brain, features =  "ITPR1", pt.size = 0.5,order = T,split.by = 'Status', cols = c("white", "#3C5488FF"), raster = F) & ggplot2::theme(legend.position = "right")
dev.off()

pdf(paste0(save_path,"GSE178265_SDC1_Feature_Plot.pdf"),width = 8, height = 4)
FeaturePlot(object = PD_brain, features =  "SDC1", pt.size = 0.5,order = T,split.by = 'Status', cols = c("white", "#3C5488FF"), raster = F) & ggplot2::theme(legend.position = "right")
dev.off()

library(ggpubr)
for(i in hubCSDEG){
  PD_brain <-AddModuleScore(PD_brain,features = i)
  data_long_m <-PD_brain@meta.data
  p.group=ggviolin(data_long_m, x = "Status", y = "Cluster1",
                   color = "Status",add = 'mean_sd',fill = 'Status',
                   add.params = list(color = "black")) + 
    stat_compare_means(comparisons = list(c("control" ,"PD")),label = "p.signif",method = 't.test') +
    scale_fill_manual(values = c("#339DB5","#C9342B")) +
    facet_wrap(~cell_type,nrow = 1) + 
    NoLegend() + labs(x = '', y = i)  + theme(axis.text.x.bottom = element_text(angle = 60,size = 12,hjust = 1))
  ggsave(filename=paste0(save_path,"GSE178265_boxplot_", i, ".pdf"),plot = p.group,width = 8,height =6)
}


# GSE141578
# CSF
rm(list = ls());gc() #删除所有环境变量整理内存空间

hubCSDEG = c("HMGCR", "HSPB1", "ITPR1", "SDC1")
if(!require(multtest))BiocManager::install("multtest")
if(!require(Seurat))install.packages("Seurat")
if(!require(dplyr))install.packages("dplyr")
if(!require(patchwork))install.packages("patchwork")
if(!require(R.utils))install.packages("R.utils")
# devtools::install_github('satijalab/seurat-data')
library(SeuratData)
library(SeuratDisk)

save_path = "/home/alan_turing/Projects/PD/ND_CL_EE/result/"

# 加载数据
library(Seurat)
setwd("/home/alan_turing/Projects/PD/ND_CL_EE/data/")
dir <- list.dirs("GSE141578")[-1]
names(dir) <- list.files("GSE141578",recursive = F)
sc_list <- list()
for(i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i])
  sc_list[[i]] <- CreateSeuratObject(counts,
                                     project = dir[i],
                                     min.cells=10,
                                     min.features=200)
}
scRNAlist <- merge(x=sc_list[[1]],y=sc_list[-1])
scRNAlist <- JoinLayers(scRNAlist)

# 初步查看Seurat对象
scRNAlist
ncol(scRNAlist)

# 数据质控（物种为鼠源）
scRNAlist[["percent.mt"]] <- PercentageFeatureSet(scRNAlist, pattern = "^MT-")
# 计算线粒体RNA的比例
scRNAlist$mitoRatio <- scRNAlist[["percent.mt"]] / 100
gc()

## 3.1 查看质控情况
head(scRNAlist@meta.data,5)

## 3.2 筛选细胞（过滤掉不符合阈值的细胞）
# 原文删除细胞的标准（满足其一即可）
# 线粒体基因比例> 10%
# 细胞的基因数（即nFeature_RNA）＞2500
# 单个细胞UMI计数（即nCount_RNA）＞10000
scRNAlist <- subset(scRNAlist, subset = nCount_RNA < 10000 & nFeature_RNA < 2500 & mitoRatio < 0.1)
gc()

ncol(scRNAlist)

# 4.细胞分群
## 4.1 表达数据标准化（转换定量的Count表达矩阵为基因表达量）

scRNAlist <- NormalizeData(scRNAlist,
                          normalization.method = "LogNormalize", 
                          scale.factor = 10000)
gc()

## 4.2 寻找高变基因（为后续的PCA与UMAP图等分析做准备）
scRNAlist <- FindVariableFeatures(scRNAlist,
                                 selection.method = "vst",
                                 nfeatures = 2000)
gc()

## 4.4 scale（针对每一种基因在所有细胞中的标准化）
# 它可以减少由于某基因的表达量过于离群所导致的误差
# 基于高变基因进行scale（默认）
scRNAlist <- ScaleData(scRNAlist)
gc()

## 4.5 PCA 
scRNAlist <- RunPCA(scRNAlist)
gc()

## 4.6 PCA结果可视化（确定后续选择的PCA维度）
### 展现出前30个维度的PCA结果（舍弃掉矩阵分离不显著的维度）
pdf(paste0(save_path, "DimHeatmap.pdf"), width = 10, height =8)
DimHeatmap(scRNAlist, dims = 1:30, cells = 500, balanced = TRUE)
dev.off()


pdf(paste0(save_path, "ElbowPlot.pdf"), width = 10, height =8)
ElbowPlot(scRNAlist)
dev.off()

## 4.7 根据选取的PCA维度进行细胞分群
# dim = 10
scRNAlist <- FindNeighbors(scRNAlist, dims = 1:10)

gc()

## 4.8 分群与细胞注释
# install.packages("clustree")
library(clustree)

# 预计有9个细胞类群（PMID：36503256）
resolutions <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1)
for (res in resolutions) {
  scRNAlist <- FindClusters(scRNAlist, resolution = res)
}

# Clustree动态分群评估
pdf(paste0(save_path, "clustree_plot.pdf"), width = 12, height = 12)
clustree_Plot <- clustree(scRNAlist@meta.data, prefix = "RNA_snn_res.")
clustree_Plot
dev.off()
clustree_Plot


# resolution = 0.2
# 5.分群后的可视化（绘制降维图，降维并观察细胞分群的情况）
## 5.1 降维分群
### A. UMAP（适合大型数据，速度快）
scRNAlist <- RunUMAP(scRNAlist, dims = 1:10)

Idents(scRNAlist)="RNA_snn_res.0.2"
scRNAlist$RNA_snn_res.0.2 <-Idents(scRNAlist)
library(ggsci)

# mycol <- ggsci::pal_npg()(9)
mycol = c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","#DC0000FF","#7E6148FF","#B09C85FF","#BA5B59","#7976A3")

pdf(paste0(save_path, "GSE141578_umap_before_cell_annotation.pdf"), width = 10, height =8)
DimPlot(scRNAlist, reduction = "umap", label = F, cols = mycol, raster = F)
dev.off()
gc()


## 5.2 寻找marker基因（此处细胞亚群 = cluster）
# 寻找所有分群的所有高表达marker基因
# install.packages('devtools')
# devtools::install_github('immunogenomics/presto')
# 提取安装好presto包可加速Seurat中的Wilcoxon秩和检验
scRNAlist.markers <- FindAllMarkers(scRNAlist,
                                    only.pos = TRUE, #只会去找特定上调的基因
                                    min.pct = 0.1, #基因至少在该cluster多少比例细胞中表达 
                                    logfc.threshold = 0.25 #倍数的对数的阈值，默认是0.25
)
gc()
## 获得的结果中，cluster为目标细胞分群
## 本质上是通过对比目标细胞分群（ident.1）和其他所有细胞分群（ident.2），寻找marker

### 提取细胞亚群的特异性marker
if(!require(dplyr))install.packages("dplyr")
# N = 75
top_N = scRNAlist.markers %>% group_by(cluster) %>% top_n(n = 75, wt = avg_log2FC) # 提取各细胞亚群的前2个特异性marker

write.csv(top_N,paste0(save_path,"GSE141578_top_N_markers.csv"))
write.csv(scRNAlist.markers,paste0(save_path,"GSE141578_scRNAlist.markers.csv"))

# 人工细胞注释（结合文献与CellMarker数据库）
# 参考文献（PMID）：34648304（原文献），31937773，36516855

# cluster 0: IL7R (CD4+ T cell)
# cluster 1: LTB (CD4+ T cell)
# cluster 2: CD8A + CD8B(CD8+ T cell)
# cluster 3: LTB (CD4+ T cell),
# cluster 4: CD8B + CCL5 (CD8+ T cell)
# cluster 5: PTPRCAP, PPA1 (dendritic cell)
# cluster 6: CD1C, FCER1A, CLEC10A (Dendritic cells)
# cluster 7: S100A8, S100A9 (Dendritic cells)
# cluster 8: NKG7, GNLY (NK cell)
# cluster 9: NKG7, GNLY (NK cell)
# clutser 10: JCHAIN (Plasma cell)
# cluster 11: PDE3B (CD8+ T cell)

# 手动重命名各细胞亚群（采用向量法而非表格法）
# 表格法（即NM_PD针对GSE178265的所用方法）容易出现cluster下细胞为NA问题（尤其是cluster数目和细胞种类数目相差较大时）
# 下面的向量法可对齐各cluster（只要位置一一对应即可）
new.cluster.ids <- c("CD4+ T cell",
                     "CD4+ T cell",
                     "CD8+ T cell",
                     "CD4+ T cell",
                     "CD8+ T cell",
                     "Dendritic cell",
                     "Dendritic cell",
                     "Dendritic cell",
                     "NK cell",
                     "NK cell",
                     "Plasma cell",
                     "CD8+ T cell")

names(new.cluster.ids) <- levels(scRNAlist)
scRNAlist <- RenameIdents(scRNAlist, new.cluster.ids)
scRNAlist$celltype <-Idents(scRNAlist)

library(ggsci)
mycol = ggsci::pal_npg()(5)
mycol 

# Idents(scRNAlist)="celltype"
pdf(paste0(save_path, "GSE141578_umap_after_cell_annotation.pdf"), width = 10, height =8)
DimPlot(scRNAlist, reduction = "umap", label = F, cols = mycol, raster = F)
dev.off()
gc()

# marker基因的dotplot绘制及美化
final_annotation_markers = c("IL7R","LTB","CD8A","CD8B","CCL5","PDE3B","PTPRCAP","PPA1","CD1C","FCER1A","CLEC10A","S100A8","S100A9","NKG7","GNLY","JCHAIN")

pdf(paste0(save_path,"GSE141578_cell_marker_dot_plot.pdf"), width = 10, height = 8)
DotPlot(scRNAlist, features = final_annotation_markers, cols = c("white", "#3C5488FF")) + ggplot2::coord_flip() + ggplot2::xlab(label = "") + ggplot2::ylab(label = "")
dev.off()

# 添加组别信息
scRNAlist$group <- scRNAlist$orig.ident
Idents(scRNAlist) <- "group"
scRNAlist <- RenameIdents(scRNAlist,
                          "GSM4208766" = "PD",
                          "GSM4208768" = "PD",
                          "GSM4208769" = "PD",
                          "GSM4208770" = "PD",
                          "GSM4208771" = "PD",
                          "GSM4208772" = "control",
                          "GSM4208773" = "control",
                          "GSM4208774" = "control",
                          "GSM4208775" = "control",
                          "GSM4208776" = "control",
                          "GSM4208777" = "control",
                          "GSM4208778" = "control",
                          "GSM4208779" = "control",
                          "GSM4208780" = "control",
                          "GSM4404053" = "PD",
                          "GSM4404054" = "PD",
                          "GSM4404055" = "control",
                          "GSM4404056" = "PD",
                          "GSM4404057" = "PD",
                          "GSM4404059" = "control")
scRNAlist$group <- scRNAlist@active.ident
Idents(scRNAlist) = "orig.ident"
Idents(scRNAlist) = "celltype"

pdf("GSE141578_HMGCR_Feature_Plot.pdf",width = 8, height = 4)
FeaturePlot(object = scRNAlist, features =  "HMGCR", pt.size = 0.5,order = T,split.by = 'group', cols = c("white", "#3C5488FF"), raster = F) & ggplot2::theme(legend.position = "right")
dev.off()

pdf("GSE141578_HSPB1_Feature_Plot.pdf",width = 8, height = 4)
FeaturePlot(object = scRNAlist, features =  "HSPB1", pt.size = 0.5,order = T,split.by = 'group', cols = c("white", "#3C5488FF"), raster = F) & ggplot2::theme(legend.position = "right")
dev.off()

pdf("GSE141578_ITPR1_Feature_Plot.pdf",width = 8, height = 4)
FeaturePlot(object = scRNAlist, features =  "ITPR1", pt.size = 0.5,order = T,split.by = 'group', cols = c("white", "#3C5488FF"), raster = F) & ggplot2::theme(legend.position = "right")
dev.off()

pdf("GSE141578_SDC1_Feature_Plot.pdf",width = 8, height = 4)
FeaturePlot(object = scRNAlist, features =  "BCL2", pt.size = 0.5,order = T,split.by = 'group', cols = c("white", "#3C5488FF"), raster = F) & ggplot2::theme(legend.position = "right")
dev.off()

# Zenodo 3993994
# CSF + PBMC
rm(list = ls());gc() #删除所有环境变量整理内存空间

hubCSDEG = c("HMGCR", "HSPB1", "ITPR1", "SDC1")
if(!require(multtest))BiocManager::install("multtest")
if(!require(Seurat))install.packages("Seurat")
if(!require(dplyr))install.packages("dplyr")
if(!require(patchwork))install.packages("patchwork")
if(!require(R.utils))install.packages("R.utils")
# devtools::install_github('satijalab/seurat-data')
library(SeuratData)
library(SeuratDisk)

read_path = "/home/alan_turing/Projects/PD/ND_CL_EE/data/"
save_path = "/home/alan_turing/Projects/PD/ND_CL_EE/result/"

# 加载数据
library(Seurat)

Zenodo_3993994_data <- readRDS(paste0(read_path,"Zenodo_3993994_counts.rds"))
Zenodo_3993994_metadata <- read.table(paste0(read_path,"Zenodo_3993994_metadata.txt"),row.names = 1,sep = '\t',header = T)

Zenodo_3993994_scRNA <- CreateSeuratObject(counts = Zenodo_3993994_data, project = "Zenodo_3993994")



# PRJNA1145007
rm(list = ls());gc() #删除所有环境变量整理内存空间

hubCSDEG = c("HMGCR", "HSPB1", "ITPR1", "SDC1")
if(!require(multtest))BiocManager::install("multtest")
if(!require(Seurat))install.packages("Seurat")
if(!require(dplyr))install.packages("dplyr")
if(!require(patchwork))install.packages("patchwork")
if(!require(R.utils))install.packages("R.utils")
# devtools::install_github('satijalab/seurat-data')
library(SeuratData)
library(SeuratDisk)

read_path = "/home/alan_turing/Projects/PD/NM_PD/data/"
save_path = "/home/alan_turing/Projects/PD/ND_CL_EE/result/"

# 加载数据
library(Seurat)

PD_PBMC.data <- readRDS(paste0(read_path,"scPD_AllCells_Annotated.rds"))
ncol(PD_PBMC.data)
gc()


# 3. 去掉MSA样本与PSP的样本（即仅保留Ctrl与PD的样本）
PD_PBMC_targetsamples <- rownames(PD_PBMC.data@meta.data)[PD_PBMC.data@meta.data[["Group"]] == "CTRL" | PD_PBMC.data@meta.data[["Group"]] == "PD"]
PD_PBMC.data = subset(PD_PBMC.data, cells = PD_PBMC_targetsamples)
gc()

ncol(PD_PBMC.data)
PD_PBMC.data@meta.data$Group = gsub("CTRL","control",PD_PBMC.data@meta.data$Group)

# 4. 分群后的可视化（绘制降维图，降维并观察细胞分群的情况）
# UMAP（适合大型数据，速度快）

Idents(PD_PBMC.data)="final_broad"
PD_PBMC.data$final_broad <-Idents(PD_PBMC.data)

# 重新命名类群
# 合并单核细胞：CD14 mono/CD16 mono = monocytes, cDC/pDC = dendritic cells
# 合并B细胞：B cells/Plasmablast = B cells
# 合并血小板生成细胞和淋巴细胞：Platelets Prolif/lymphocytes = Other myeloïd
# 合并gdT、MAIT、NKT cells：gdT/MAIT/NKT cells = unconventional T cells (UC t)

new.cluster.ids <- c("B",
                     "Basophil",
                     "Monocyte",
                     "Monocyte",
                     "CD4 T",
                     "CD8 T",
                     "DC",
                     "UC T",
                     "HSC",
                     "UC T",
                     "NK",
                     "UC T",
                     "DC",
                     "B",
                     "Other myeloïd",
                     "Other myeloïd",
                     "Treg")

names(new.cluster.ids) <- levels(PD_PBMC.data)
PD_PBMC.data <- RenameIdents(PD_PBMC.data, new.cluster.ids)
PD_PBMC.data$new_cell_type <-Idents(PD_PBMC.data)


library(ggsci)
mycol_PBMC <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","#DC0000FF","#7E6148FF","#B09C85FF","#329845")

pdf(paste0(save_path,"ND_CL_EE_PRJNA1145007_umap_with_cell_annotation.pdf"), width = 10, height =8)
DimPlot(PD_PBMC.data, reduction = "umap", label = F, cols = mycol_PBMC, raster = F)
dev.off()
gc()

pdf(paste0(save_path,"ND_CL_EE_PRJNA1145007_HMGCR_Feature_Plot.pdf"),width = 8, height = 4)
FeaturePlot(object = PD_PBMC.data, features =  "HMGCR", reduction = "umap", pt.size = 0.5,order = T,split.by = 'Group', cols = c("white", "#3C5488FF"), raster = F) & ggplot2::theme(legend.position = "right")
dev.off()

pdf(paste0(save_path,"ND_CL_EE_PRJNA1145007_HSPB1_Feature_Plot.pdf"),width = 8, height = 4)
FeaturePlot(object = PD_PBMC.data, features =  "HSPB1", reduction = "umap", pt.size = 0.5,order = T,split.by = 'Group', cols = c("white", "#3C5488FF"), raster = F) & ggplot2::theme(legend.position = "right")
dev.off()

pdf(paste0(save_path,"ND_CL_EE_PRJNA1145007_ITPR1_Feature_Plot.pdf"),width = 8, height = 4)
FeaturePlot(object = PD_PBMC.data, features =  "ITPR1", reduction = "umap",  pt.size = 0.5,order = T,split.by = 'Group', cols = c("white", "#3C5488FF"), raster = F) & ggplot2::theme(legend.position = "right")
dev.off()

pdf(paste0(save_path,"ND_CL_EE_PRJNA1145007_SDC1_Feature_Plot.pdf"),width = 8, height = 4)
FeaturePlot(object = PD_PBMC.data, features =  "SDC1", reduction = "umap", pt.size = 0.5,order = T,split.by = 'Group', cols = c("white", "#3C5488FF"), raster = F) & ggplot2::theme(legend.position = "right")
dev.off()


library(ggpubr)
for(i in hubCSDEG){
  PD_PBMC.data <-AddModuleScore(PD_PBMC.data,features = i)
  data_long_m <-PD_PBMC.data@meta.data
  p.group=ggviolin(data_long_m, x = "Group", y = "Cluster1",
                   color = "Group",add = 'mean_sd',fill = 'Group',
                   add.params = list(color = "black")) + 
    stat_compare_means(comparisons = list(c("control" ,"PD")),label = "p.signif",method = 't.test') +
    scale_fill_manual(values = c("#339DB5","#C9342B")) +
    facet_wrap(~new_cell_type,nrow = 1) + 
    NoLegend() + labs(x = '', y = i)  + theme(axis.text.x.bottom = element_text(angle = 60,size = 12,hjust = 1))  + scale_x_discrete(labels=c('control', 'PD'))
  ggsave(filename=paste0(save_path,"PRJNA1145007_", i, "_boxplot",".pdf"),plot = p.group,width = 12, height = 6)
}

saveRDS(PD_PBMC.data,'ND_CL_EE_Annotated_CTRL_PD_PRJNA1145007.rds')
rm(PD_PBMC.data);gc()


