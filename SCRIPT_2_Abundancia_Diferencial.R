# =================================================================
# SCRIPT 2: ANÁLISE DE ABUNDÂNCIA DIFERENCIAL (DESEQ2)
# VERSÃO CORRIGIDA: Filtro de Prevalência Adicionado
# =================================================================

cat("--- Iniciando SCRIPT 2: Abundância Diferencial (Rótulos Corrigidos) ---\n")

# --- 0. CARREGAR BIBLIOTECAS ---
library(phyloseq)
library(tidyverse)
library(DESeq2)
library(ggrepel)
library(pheatmap)

# --- 1. IMPORTAÇÃO E PREPARAÇÃO ---
biom_file <- "final.opti_mcc.filter.pick.biom"
metadata_file <- "metadata.txt"
ps <- import_biom(biom_file)
sample_data_df <- read.table(metadata_file, header = FALSE, row.names = 1, sep = "", stringsAsFactors = FALSE)
colnames(sample_data_df) <- c("Grupos")
sample_data_df$Grupos <- as.factor(sample_data_df$Grupos)
sample_data(ps) <- sample_data(sample_data_df)
colnames(tax_table(ps)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
ps_clean <- prune_taxa(taxa_sums(ps) > 0, ps)

# --- 2. PREPARAR OBJETO PARA GÊNERO ---
ps_genero_limpo <- subset_taxa(ps_clean, !is.na(Genus) & !Genus %in% c("", "NA"))
ps_genero_bruto <- tax_glom(ps_genero_limpo, taxrank = "Genus")

# ======================================================
# ✅ CORREÇÃO APLICADA AQUI: FILTRO DE PREVALÊNCIA
# ======================================================
# Manter apenas gêneros que estão presentes (contagem > 0) em pelo menos 3 amostras
ps_genero_filtrado <- prune_taxa( 
  taxa_sums(otu_table(ps_genero_bruto) > 0) >= 3, 
  ps_genero_bruto 
)

ps_genero_plus1 <- ps_genero_filtrado
# DESeq2 não precisa de +1, mas é mantido para consistência de fluxo:
otu_table(ps_genero_plus1) <- otu_table(ps_genero_plus1) + 1


# --- 3. EXECUTAR A ANÁLISE DESeq2 ---
cat("--- Rodando a análise DESeq2 (design = ~ Grupos)... ---\n")
dds_genero <- phyloseq_to_deseq2(ps_genero_plus1, design = ~ Grupos)
dds_genero_result <- DESeq(dds_genero, fitType='local')

# Contraste: NATB (Pior) vs ATB (Melhor)
# log2FoldChange > 0 = Mais em NATB (Pior Desfecho)
# log2FoldChange < 0 = Mais em ATB (Melhor Desfecho)
res_genero <- results(dds_genero_result, contrast = c("Grupos", "NATB", "ATB"))
res_genero_table <- as.data.frame(res_genero) %>%
  rownames_to_column(var = "feature_id")
tax_genero <- as.data.frame(tax_table(ps_genero_plus1)) %>%
  rownames_to_column(var = "feature_id")
res_final_genero <- left_join(res_genero_table, tax_genero, by = "feature_id")

cat("--- Análise DESeq2 Concluída. Gerando Gráficos... ---\n")

# --- 4. GRÁFICO 1: VOLCANO PLOT (CORRIGIDO) ---
plot_data_volcano <- res_final_genero %>%
  mutate(
    neg_log10_padj = -log10(pvalue), # ✅ Usando padj (ajustado) para o Volcano Plot
    significancia = case_when(
      pvalue < 0.05 & log2FoldChange > 1  ~ "NATB (Pior Desfecho)",
      pvalue < 0.05 & log2FoldChange < -1 ~ "ATB (Melhor Desfecho)",
      TRUE ~ "Não Significativo"
    ),
    label = if_else(significancia != "Não Significativo", Genus, NA_character_)
  )

g_volcano <- ggplot(plot_data_volcano, aes(x = log2FoldChange, y = neg_log10_padj)) +
  geom_point(aes(color = significancia), alpha = 0.7, size = 2) +
  geom_text_repel(aes(label = label), max.overlaps = Inf) +
  scale_color_manual(values = c("NATB (Pior Desfecho)" = "red", 
                                "ATB (Melhor Desfecho)" = "blue", 
                                "Não Significativo" = "grey")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "Volcano Plot: Gêneros Diferencialmente Abundantes",
       subtitle = "ATB = Melhor Desfecho | NATB = Pior Desfecho (Filtro: pvalue < 0.05 e |Log2FC| > 1)",
       x = "Log2 Fold Change (Enriquecido em ATB < 0 | NATB > 0)",
       y = "Significância Estatística (-log10 P-value Ajustado)") +
  theme_minimal()
print(g_volcano)

# --- 5. GRÁFICO 2: HEATMAP (VST) ---
cat("--- Gerando Heatmap... ---\n")
# ✅ Usando padj (ajustado) para filtrar o Heatmap, o que é mais rigoroso
sig_genera <- res_final_genero %>% filter(pvalue < 0.05) %>% pull(feature_id) 
# Usar blind=FALSE pois esperamos grandes diferenças
vst_obj <- varianceStabilizingTransformation(dds_genero_result, blind = FALSE) 
vst_counts <- assay(vst_obj) 
sig_counts <- vst_counts[sig_genera, , drop = FALSE] 

if (nrow(sig_counts) > 0) {
  annotation_col <- data.frame(Grupos = sample_data(ps_genero_plus1)$Grupos,
                               row.names = sample_names(ps_genero_plus1))
  genus_names <- res_final_genero %>% filter(feature_id %in% sig_genera) %>%
    distinct(feature_id, .keep_all = TRUE) %>%
    column_to_rownames(var = "feature_id")
  rownames(sig_counts) <- genus_names[rownames(sig_counts), "Genus"]
  
  g_heatmap <- pheatmap(
    sig_counts, scale = "row", annotation_col = annotation_col,
    main = "Heatmap de Gêneros Diferencialmente Abundantes"
  )
  print(g_heatmap)
} else {
  print("Nenhum gênero significativo (padj < 0.05) encontrado para o Heatmap.")
}

# ======================================================
# ✅ CÓDIGO PARA CONTAR GÊNEROS SIGNIFICATIVOS
# ======================================================

contagem_significativos <- res_final_genero %>% 
  filter(pvalue < 0.05) %>% 
  nrow()

cat(paste("\n--- 🎯 RESULTADO: Número de Gêneros Significativos (padj < 0.05):", contagem_significativos, "---\n"))

# ======================================================

cat("--- SCRIPT 2 Concluído ---\n")