# =================================================================
# SCRIPT: ANÁLISE DE REDE COMPLEXA (GRÁFICOS INDIVIDUAIS)
# =================================================================
# Objetivo: 1. Identificar Promotores, Antagonistas, Mediadores e Guarda-Costas.
#           2. Gerar Box Plots INDIVIDUAIS (por Antagonista)
#              da abundância (VST) dos Mediadores e Guarda-Costas.

cat("--- Iniciando Análise de Rede Complexa (Gráficos Individuais) ---\n")

# --- 0. CARREGAR BIBLIOTECAS ---
library(phyloseq)
library(tidyverse)
library(compositions)
library(Hmisc) # Para rcorr (Rede de Correlação)
library(ggplot2)
library(DESeq2) # Para VST

# --- 1. SETUP INICIAL E PREPARAÇÃO DOS DADOS ---

# 1.1. Importação Básica
biom_file <- "final.opti_mcc.filter.pick.biom"
metadata_file <- "metadata.txt"
ps <- import_biom(biom_file) 
# Correção para o erro de 'row.names = 1'
sample_data_df <- read.table(metadata_file, header = FALSE, sep = "", stringsAsFactors = FALSE)
rownames(sample_data_df) <- sample_data_df$V1
sample_data_df$V1 <- NULL
colnames(sample_data_df) <- c("Grupos")
sample_data_df$Grupos <- as.factor(sample_data_df$Grupos)
sample_data(ps) <- sample_data(sample_data_df)

# 1.2. Integração e Filtros
desenvolvimento_df <- data.frame(
  NAME = c("JC_1", "JC_10", "JC_11", "JC_12", "JC_13", "JC_14", "JC_2", "JC_3", "JC_4", "JC_5", "JC_6", "JC_7", "JC_8", "JC_9"),
  premorula = c(0, 0, 4, 0, 0, 0, 0, 0, 7, 0, 14, 0, 0, 1),
  morula = c(0, 0, 0, 15, 1, 4, 15, 11, 5, 0, 0, 0, 0, 1),
  blastocisto = c(12, 11, 5, 1, 5, 1, 0, 2, 0, 0, 0, 9, 7, 0),
  row.names = 1
)
desenvolvimento_df <- desenvolvimento_df %>% mutate(Total_Embrioes = premorula + morula + blastocisto)
novos_metadados <- sample_data(ps)
novos_metadados <- cbind(novos_metadados, desenvolvimento_df[rownames(novos_metadados), c("premorula", "morula", "blastocisto", "Total_Embrioes")])
sample_data(ps) <- sample_data(novos_metadados)
amostras_para_remover <- rownames(sample_data(ps))[sample_data(ps)$Total_Embrioes == 0]
ps_filtrado_desenv <- prune_samples(!sample_names(ps) %in% amostras_para_remover, ps)
ps_filtrado_taxa <- prune_taxa(taxa_sums(ps_filtrado_desenv) >= 1, ps_filtrado_desenv)
ps_genero_bruto <- tax_glom(ps_filtrado_taxa, taxrank = "Rank6", NArm = FALSE)
ps_genero_filtrado <- prune_taxa(taxa_sums(otu_table(ps_genero_bruto) > 0) >= 3, ps_genero_bruto)
otu_bruto_genero_pseudo <- as.data.frame(otu_table(ps_genero_filtrado)) + 1
clr_data_matrix <- t(clr(t(otu_bruto_genero_pseudo)))
taxa_abund_clr <- as.data.frame(clr_data_matrix)
metadados_filtrados <- as(sample_data(ps_genero_filtrado), "data.frame") %>% mutate(Taxa_Blastocisto = blastocisto / Total_Embrioes)
amostras_atb <- rownames(metadados_filtrados[metadados_filtrados$Grupos == "ATB", ])
amostras_natb <- rownames(metadados_filtrados[metadados_filtrados$Grupos == "NATB", ])
taxa_abund_atb <- taxa_abund_clr[, amostras_atb]
taxa_abund_natb <- taxa_abund_clr[, amostras_natb]

# Funções essenciais (Corrigido)
get_tax_name <- function(otu_id, ps_obj, rank_level = "Rank6") {
  tax_info <- as(tax_table(ps_obj)[otu_id, ], "matrix"); name <- tax_info[1, rank_level]
  if (is.na(name) || name == "" || grepl("_unclassified", name) || grepl("_group", name)) {
    name_rank5 <- try(as.character(tax_info[1, "Rank5"]), silent = TRUE)
    if (!inherits(name_rank5, "try-error") && !is.na(name_rank5) && name_rank5 != "") { name <- paste0(name_rank5, " (Família)") } else { name <- otu_id }
  }
  return(as.character(name))
}
run_spearman_analysis <- function(count_vector, taxa_abund_df) {
  correlacoes <- apply(taxa_abund_df, 1, function(x) {
    if (sd(x, na.rm = TRUE) == 0) return(list(estimate = NA, p.value = NA)); suppressWarnings(cor.test(x, count_vector, method = "spearman"))
  }); corr_df <- as.data.frame(t(sapply(correlacoes, function(x) {
    c(Rho = x$estimate, P_valor_Bruto = x$p.value)
  }))); return(corr_df %>% na.omit())
}

# --- 2. IDENTIFICAR PROMOTORES E ANTAGONISTAS ---
cat("--- 2. Identificando Promotores (ATB) e Antagonistas (NATB)... ---\n")
TopN_select <- 6 
blastocisto_taxa_atb <- metadados_filtrados[amostras_atb, "Taxa_Blastocisto"]
corr_atb_selecao <- run_spearman_analysis(blastocisto_taxa_atb, taxa_abund_atb)
colnames(corr_atb_selecao)<-c("Rho","P_valor_Bruto")
top_promotores_ids <- corr_atb_selecao %>% filter(Rho > 0) %>% arrange(desc(Rho)) %>% head(TopN_select) %>% rownames()

taxa_abund_natb_t <- t(taxa_abund_natb)
cor_matrix_natb <- Hmisc::rcorr(as.matrix(taxa_abund_natb_t), type = "spearman")
antagonistas_encontrados <- data.frame()
for (promotor_id in top_promotores_ids) {
  rho_vector <- cor_matrix_natb$r[promotor_id, ]; p_value_vector <- cor_matrix_natb$P[promotor_id, ]
  df_parceiros <- data.frame(Parceiro_ID = names(rho_vector), Rho_NATB = as.numeric(rho_vector), P_valor_Bruto_NATB = as.numeric(p_value_vector)) %>%
    filter(Parceiro_ID != promotor_id) %>% na.omit()
  antagonistas_do_promotor <- df_parceiros %>%
    filter(Rho_NATB < 0, P_valor_Bruto_NATB < 0.05) %>% mutate(Promotor_ATB_ID = promotor_id)
  antagonistas_encontrados <- rbind(antagonistas_encontrados, antagonistas_do_promotor)
}
if (nrow(antagonistas_encontrados) == 0) { stop("Nenhum Antagonista (Rho < 0, P < 0.05) foi encontrado. Abortando.") }
antagonistas_encontrados$Nome_Antagonista <- sapply(antagonistas_encontrados$Parceiro_ID, get_tax_name, ps_obj = ps_genero_filtrado)
antagonistas_encontrados$Nome_Promotor <- sapply(antagonistas_encontrados$Promotor_ATB_ID, get_tax_name, ps_obj = ps_genero_filtrado)


# --- 3. BUSCAR "MEDIADORES" (Parceiros dos Antagonistas na Rede NATB) ---
cat("--- 3. Buscando Mediadores (Parceiros de Antagonistas no NATB)... ---\n")
mediadores_encontrados <- data.frame()
antagonistas_ids_unicos_natb <- unique(antagonistas_encontrados$Parceiro_ID)

for (antagonista_id in antagonistas_ids_unicos_natb) {
  rho_vector <- cor_matrix_natb$r[antagonista_id, ]; p_vector <- cor_matrix_natb$P[antagonista_id, ]
  df_parceiros_med <- data.frame(Mediador_ID = names(rho_vector), Rho_NATB = as.numeric(rho_vector), P_Valor = as.numeric(p_vector)) %>%
    filter(!Mediador_ID %in% top_promotores_ids, Mediador_ID != antagonista_id) %>% na.omit()
  mediadores_do_antagonista <- df_parceiros_med %>%
    filter(Rho_NATB > 0, P_Valor < 0.05) %>% 
    mutate(Antagonista_Pai_ID = antagonista_id)
  mediadores_encontrados <- rbind(mediadores_encontrados, mediadores_do_antagonista)
}
if (nrow(mediadores_encontrados) > 0) {
  mediadores_encontrados$Nome_Mediador <- sapply(mediadores_encontrados$Mediador_ID, get_tax_name, ps_obj = ps_genero_filtrado)
  mediadores_encontrados$Nome_Antagonista_Pai <- sapply(mediadores_encontrados$Antagonista_Pai_ID, get_tax_name, ps_obj = ps_genero_filtrado)
}

# --- 4. BUSCAR "GUARDA-COSTAS" (Inibidores dos Antagonistas na Rede ATB) ---
cat("--- 4. Buscando Guarda-Costas (Inibidores de Antagonistas no ATB)... ---\n")
taxa_abund_atb_t <- t(taxa_abund_atb)
cor_matrix_atb <- Hmisc::rcorr(as.matrix(taxa_abund_atb_t), type = "spearman")
guardas_encontrados <- data.frame()

for (antagonista_id in antagonistas_ids_unicos_natb) { # Usamos os mesmos antagonistas
  if(antagonista_id %in% rownames(cor_matrix_atb$r)) {
    rho_vector <- cor_matrix_atb$r[antagonista_id, ]; p_vector <- cor_matrix_atb$P[antagonista_id, ]
    df_parceiros_gua <- data.frame(Guarda_ID = names(rho_vector), Rho_ATB = as.numeric(rho_vector), P_Valor = as.numeric(p_vector)) %>%
      filter(!Guarda_ID %in% top_promotores_ids, Guarda_ID != antagonista_id) %>% na.omit()
    guardas_do_antagonista <- df_parceiros_gua %>%
      filter(Rho_ATB < 0, P_Valor < 0.05) %>% 
      mutate(Antagonista_Alvo_ID = antagonista_id)
    guardas_encontrados <- rbind(guardas_encontrados, guardas_do_antagonista)
  }
}
if (nrow(guardas_encontrados) > 0) {
  guardas_encontrados$Nome_Guarda <- sapply(guardas_encontrados$Guarda_ID, get_tax_name, ps_obj = ps_genero_filtrado)
  guardas_encontrados$Nome_Antagonista_Alvo <- sapply(guardas_encontrados$Antagonista_Alvo_ID, get_tax_name, ps_obj = ps_genero_filtrado)
}

# --- 5. RODAR DESEQ2/VST (Necessário para os Box Plots) ---
cat("--- 5. Rodando DESeq2/VST no objeto filtrado... ---\n")
ps_para_deseq <- ps_genero_filtrado
otu_table(ps_para_deseq) <- otu_table(ps_para_deseq) + 1 
dds_obj_corrigido <- phyloseq_to_deseq2(ps_para_deseq, design = ~ Grupos)
dds_result_corrigido <- DESeq(dds_obj_corrigido, fitType='local')
vst_obj <- varianceStabilizingTransformation(dds_result_corrigido, blind = FALSE) 
vst_counts_matrix <- assay(vst_obj)
metadados_grupos <- data.frame(sample_data(ps_genero_filtrado)) %>% rownames_to_column("Sample") %>% select(Sample, Grupos)

# --- 6. (NOVO) LOOP PARA GRÁFICOS (MEDIADORES) ---
if (nrow(mediadores_encontrados) > 0) {
  cat("--- 6. Gerando Box Plots Individuais de Mediadores... ---\n")
  mediadores_ids_unicos <- unique(mediadores_encontrados$Mediador_ID)
  vst_mediadores <- vst_counts_matrix[mediadores_ids_unicos, , drop = FALSE]
  
  plot_data_vst_med <- vst_mediadores %>%
    t() %>% as.data.frame() %>% rownames_to_column("Sample") %>%
    left_join(metadados_grupos, by = "Sample") %>%
    pivot_longer(cols = all_of(mediadores_ids_unicos), names_to = "Mediador_ID", values_to = "Abundancia_VST")
  
  # Loop pelos Antagonistas que TÊM mediadores
  lista_antagonistas_pais <- unique(mediadores_encontrados$Antagonista_Pai_ID)
  
  for (antagonista_pai_id in lista_antagonistas_pais) {
    nome_antagonista_pai <- get_tax_name(antagonista_pai_id, ps_genero_filtrado)
    
    # Identificar os mediadores deste antagonista
    pares_mediadores_filtrado <- mediadores_encontrados %>%
      filter(Antagonista_Pai_ID == antagonista_pai_id)
    
    df_plot_med_individual <- plot_data_vst_med %>%
      filter(Mediador_ID %in% pares_mediadores_filtrado$Mediador_ID) %>%
      left_join(pares_mediadores_filtrado %>% select(Mediador_ID, Nome_Mediador), by = "Mediador_ID", relationship = "many-to-many") %>%
      distinct()
    
    g_boxplot_med <- ggplot(df_plot_med_individual, aes(x = Nome_Mediador, y = Abundancia_VST, fill = Grupos)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      geom_jitter(position = position_jitterdodge(jitter.width = 0.2), size = 1.5, alpha = 0.7) +
      labs(title = paste("Abundância (VST) dos 'Mediadores' do Antagonista:", nome_antagonista_pai),
           subtitle = "Hipótese: Mediadores (Eixo X) são MAIORES no NATB (Vermelho)",
           x = "Gênero Mediador (Parceiro do Antagonista)", y = "Abundância Normalizada (VST)") +
      scale_fill_manual(values = c("NATB" = "red", "ATB" = "blue")) +
      theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    nome_arquivo_plot <- paste0("grafico_abundancia_mediadores_DE_", gsub(" |\\(|\\)|/", "_", nome_antagonista_pai), ".png")
    ggsave(nome_arquivo_plot, g_boxplot_med, width = 8, height = 6, dpi = 300)
    cat(paste("--- Gráfico salvo:", nome_arquivo_plot, "---\n"))
  }
} else {
  cat("--- 6. Nenhum Mediador encontrado. ---\n")
}

# --- 7. (NOVO) LOOP PARA GRÁFICOS (GUARDA-COSTAS) ---
if (nrow(guardas_encontrados) > 0) {
  cat("--- 7. Gerando Box Plots Individuais de Guarda-Costas... ---\n")
  guardas_ids_unicos <- unique(guardas_encontrados$Guarda_ID)
  vst_guardas <- vst_counts_matrix[guardas_ids_unicos, , drop = FALSE]
  
  plot_data_vst_gua <- vst_guardas %>%
    t() %>% as.data.frame() %>% rownames_to_column("Sample") %>%
    left_join(metadados_grupos, by = "Sample") %>%
    pivot_longer(cols = all_of(guardas_ids_unicos), names_to = "Guarda_ID", values_to = "Abundancia_VST")
  
  # Loop pelos Antagonistas que TÊM guarda-costas
  lista_antagonistas_alvo <- unique(guardas_encontrados$Antagonista_Alvo_ID)
  
  for (antagonista_alvo_id in lista_antagonistas_alvo) {
    nome_antagonista_alvo <- get_tax_name(antagonista_alvo_id, ps_genero_filtrado)
    
    pares_guardas_filtrado <- guardas_encontrados %>%
      filter(Antagonista_Alvo_ID == antagonista_alvo_id)
    
    df_plot_gua_individual <- plot_data_vst_gua %>%
      filter(Guarda_ID %in% pares_guardas_filtrado$Guarda_ID) %>%
      left_join(pares_guardas_filtrado %>% select(Guarda_ID, Nome_Guarda), by = "Guarda_ID", relationship = "many-to-many") %>%
      distinct()
    
    g_boxplot_gua <- ggplot(df_plot_gua_individual, aes(x = Nome_Guarda, y = Abundancia_VST, fill = Grupos)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      geom_jitter(position = position_jitterdodge(jitter.width = 0.2), size = 1.5, alpha = 0.7) +
      labs(title = paste("Abundância (VST) dos 'Guarda-Costas' do Antagonista:", nome_antagonista_alvo),
           subtitle = "Hipótese: Guarda-Costas (Eixo X) são MAIORES no ATB (Azul)",
           x = "Gênero Guarda-Costas (Inibidor do Antagonista)", y = "Abundância Normalizada (VST)") +
      scale_fill_manual(values = c("NATB" = "red", "ATB" = "blue")) +
      theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    nome_arquivo_plot <- paste0("grafico_abundancia_guardacostas_DE_", gsub(" |\\(|\\)|/", "_", nome_antagonista_alvo), ".png")
    ggsave(nome_arquivo_plot, g_boxplot_gua, width = 8, height = 6, dpi = 300)
    cat(paste("--- Gráfico salvo:", nome_arquivo_plot, "---\n"))
  }
} else {
  cat("--- 7. Nenhum Guarda-Costas encontrado. ---\n")
}

cat("\n--- ANÁLISE COMPLETA DE REDE COMPLEXA (INDIVIDUAL) CONCLUÍDA ---\n")