# =================================================================
# SCRIPT: ANÁLISE COMPLETA DE ANTAGONISTAS (INIBIDORES NATB)
# (Versão Corrigida do Erro '$p.value')
# =================================================================

cat("--- Iniciando Análise Completa de Antagonistas Inibidores ---\n")

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
sample_data_df <- read.table(metadata_file, header = FALSE, row.names = 1, sep = "", stringsAsFactors = FALSE)
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

# Funções essenciais
get_tax_name <- function(otu_id, ps_obj, rank_level = "Rank6") {
  tax_info <- as(tax_table(ps_obj)[otu_id, ], "matrix"); name <- tax_info[1, rank_level]
  if (is.na(name) || name == "" || grepl("_unclassified", name) || grepl("_group", name)) {
    name_rank5 <- try(as.character(tax_info[1, "Rank5"]), silent = TRUE)
    if (!inherits(name_rank5, "try-error") && !is.na(name_rank5) && name_rank5 != "") { name <- paste0(name_rank5, " (Família)") } else { name <- otu_id }
  }
  return(as.character(name))
}

# ======================================================
# --- ✅ CORREÇÃO DO ERRO '$' APLICADA AQUI ---
# ======================================================
run_spearman_analysis <- function(count_vector, taxa_abund_df) {
  correlacoes <- apply(taxa_abund_df, 1, function(x) {
    if (sd(x, na.rm = TRUE) == 0) return(list(estimate = NA, p.value = NA)); suppressWarnings(cor.test(x, count_vector, method = "spearman"))
  }); corr_df <- as.data.frame(t(sapply(correlacoes, function(x) {
    # CORRIGIDO: Era $p.value, agora é x$p.value
    c(Rho = x$estimate, P_valor_Bruto = x$p.value) 
  }))); return(corr_df %>% na.omit())
}
# ======================================================

# --- 2. IDENTIFICAR PROMOTORES ATB E SEUS ANTAGONISTAS ---
cat("--- 2. Identificando Promotores ATB e seus Antagonistas (Rede NATB)... ---\n")
TopN_select <- 6 

blastocisto_taxa_atb <- metadados_filtrados[amostras_atb, "Taxa_Blastocisto"]
corr_atb_selecao <- run_spearman_analysis(blastocisto_taxa_atb, taxa_abund_atb)
colnames(corr_atb_selecao)<-c("Rho","P_valor_Bruto")
top_promotores_ids <- corr_atb_selecao %>%
  filter(Rho > 0) %>% 
  arrange(desc(Rho)) %>%
  head(TopN_select) %>%
  rownames()

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

if (nrow(antagonistas_encontrados) == 0) {
  stop("Nenhum Antagonista Sinergista (Rho < 0, P < 0.05) foi encontrado na rede NATB. Abortando.")
}
antagonistas_encontrados$Nome_Antagonista <- sapply(antagonistas_encontrados$Parceiro_ID, get_tax_name, ps_obj = ps_genero_filtrado)
antagonistas_encontrados$Nome_Promotor <- sapply(antagonistas_encontrados$Promotor_ATB_ID, get_tax_name, ps_obj = ps_genero_filtrado)


# --- 3. GRÁFICO 1: LOOP DE CORRELAÇÃO (SCATTER PLOT POR PROMOTOR) ---
cat("--- 3. Gerando Gráficos de Correlação (Promotor vs. Antagonista)... ---\n")

# Dataframe de abundância CLR (todas as amostras) e Grupos (COM CORREÇÃO)
df_abund_grupos <- t(taxa_abund_clr) %>% 
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  left_join(
    metadados_filtrados %>% 
      rownames_to_column("Sample"), # Usar o rowname se 'NAME' falhou
    by = "Sample"
  )

for (promotor_id in top_promotores_ids) {
  
  nome_promotor_atual <- get_tax_name(promotor_id, ps_genero_filtrado)
  
  antagonistas_deste_promotor_df <- antagonistas_encontrados %>%
    filter(Promotor_ATB_ID == promotor_id)
  
  if (nrow(antagonistas_deste_promotor_df) == 0) {
    cat(paste("--- Aviso: Promotor", nome_promotor_atual, "não possui antagonistas na rede NATB. Pulando gráfico de correlação. ---\n"))
    next
  }
  
  antagonistas_ids_para_plotar <- antagonistas_deste_promotor_df$Parceiro_ID
  
  df_plot_scatter <- df_abund_grupos %>%
    select(Sample, Grupos, 
           Abund_Promotor = all_of(promotor_id), 
           all_of(antagonistas_ids_para_plotar)) %>%
    pivot_longer(
      cols = all_of(antagonistas_ids_para_plotar),
      names_to = "OTU_ID_Antagonista",
      values_to = "Abund_Antagonista"
    ) %>%
    left_join(antagonistas_encontrados %>% select(OTU_ID_Antagonista = Parceiro_ID, Nome_Antagonista), 
              by = "OTU_ID_Antagonista")
  
  g_scatter_antagonista <- ggplot(df_plot_scatter, 
                                  aes(x = Abund_Antagonista, y = Abund_Promotor, color = Grupos)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, fullrange = TRUE) +
    facet_wrap(~ Nome_Antagonista, scales = "free") +
    labs(
      title = paste("Antagonismo (Rede NATB):", nome_promotor_atual, "vs. Seus Inibidores"),
      subtitle = "Curva Vermelha (NATB) deve ser negativa; Curva Azul (ATB) deve ser plana/positiva (Quebra de Inibição)",
      x = "Abundância (CLR) do Gênero Antagonista",
      y = paste("Abundância (CLR) de", nome_promotor_atual),
      color = "Grupo"
    ) +
    scale_color_manual(values = c("NATB" = "red", "ATB" = "blue")) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "italic"))
  
  nome_arquivo_scatter <- paste0("grafico_correlacao_antagonistas_de_", gsub(" |\\(|\\)|/", "_", nome_promotor_atual), ".png")
  ggsave(nome_arquivo_scatter, g_scatter_antagonista, width = 10, height = 8, dpi = 300)
  cat(paste("--- Gráfico salvo:", nome_arquivo_scatter, "---\n"))
}


# =================================================================
# SCRIPT: BOX PLOT DE ANTAGONISTAS (ANINHADO POR PROMOTOR)
# =================================================================
# Objetivo: 1. Identificar Promotores ATB e seus Antagonistas (Rede NATB).
#           2. Plotar a abundância (Box Plot) dos Antagonistas,
#              separados (faceted) por Promotor.

cat("--- Iniciando Análise de Antagonistas (Aninhado por Promotor) ---\n")

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
sample_data_df <- read.table(metadata_file, header = FALSE, row.names = 1, sep = "", stringsAsFactors = FALSE)
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

# --- 2. IDENTIFICAR PROMOTORES ATB E SEUS ANTAGONISTAS ---
cat("--- 2. Identificando Promotores ATB e seus Antagonistas (Rede NATB)... ---\n")
TopN_select <- 6 

blastocisto_taxa_atb <- metadados_filtrados[amostras_atb, "Taxa_Blastocisto"]
corr_atb_selecao <- run_spearman_analysis(blastocisto_taxa_atb, taxa_abund_atb)
colnames(corr_atb_selecao)<-c("Rho","P_valor_Bruto")
top_promotores_ids <- corr_atb_selecao %>%
  filter(Rho > 0) %>% 
  arrange(desc(Rho)) %>%
  head(TopN_select) %>%
  rownames()

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

if (nrow(antagonistas_encontrados) == 0) {
  stop("Nenhum Antagonista (Rho < 0, P < 0.05) foi encontrado na rede NATB. Abortando.")
}
antagonistas_encontrados$Nome_Antagonista <- sapply(antagonistas_encontrados$Parceiro_ID, get_tax_name, ps_obj = ps_genero_filtrado)
antagonistas_encontrados$Nome_Promotor <- sapply(antagonistas_encontrados$Promotor_ATB_ID, get_tax_name, ps_obj = ps_genero_filtrado)


# --- 3. RODAR DESEQ2/VST NO OBJETO FILTRADO ---
cat("--- 3. Rodando DESeq2/VST no objeto filtrado (ps_genero_filtrado)... ---\n")

ps_para_deseq <- ps_genero_filtrado
otu_table(ps_para_deseq) <- otu_table(ps_para_deseq) + 1 
dds_obj_corrigido <- phyloseq_to_deseq2(ps_para_deseq, design = ~ Grupos)
dds_result_corrigido <- DESeq(dds_obj_corrigido, fitType='local')
vst_obj <- varianceStabilizingTransformation(dds_result_corrigido, blind = FALSE) 
vst_counts_matrix <- assay(vst_obj)

# --- 4. PREPARAR DATAFRAME PARA PLOTAGEM (ANINHADO) ---
cat("--- 4. Preparando dados aninhados para ggplot... ---\n")

# A. Identificar IDs únicos dos Antagonistas
antagonistas_ids_unicos <- unique(antagonistas_encontrados$Parceiro_ID)

# B. Filtrar a matriz VST para conter apenas os Antagonistas
vst_antagonistas <- vst_counts_matrix[antagonistas_ids_unicos, , drop = FALSE]

# C. Obter Metadados (Grupos)
metadados_grupos <- data.frame(sample_data(ps_genero_filtrado)) %>%
  rownames_to_column("Sample") %>%
  select(Sample, Grupos)

# D. Criar o Dataframe "Melted" de VST (Abundância dos Antagonistas)
plot_data_vst_antag <- vst_antagonistas %>%
  t() %>% as.data.frame() %>%
  rownames_to_column("Sample") %>%
  left_join(metadados_grupos, by = "Sample") %>%
  pivot_longer(
    cols = all_of(antagonistas_ids_unicos),
    names_to = "OTU_ID_Antagonista",
    values_to = "Abundancia_VST"
  )

# E. Criar a Tabela de Pares (O que aninha o quê)
pares_antagonismo <- antagonistas_encontrados %>%
  select(OTU_ID_Antagonista = Parceiro_ID, Nome_Promotor, Nome_Antagonista) %>%
  distinct() # Garantir que cada par (Antagonista-Promotor) seja único

# F. Juntar a Abundância VST com a Tabela de Pares
df_plot_aninhado <- plot_data_vst_antag %>%
  left_join(pares_antagonismo, by = "OTU_ID_Antagonista", relationship = "many-to-many") %>%
  # Se um Antagonista não tiver Promotor (não deveria acontecer, mas por segurança)
  filter(!is.na(Nome_Promotor)) 

# --- 5. GERAR O GRÁFICO BOX PLOT (ANINHADO) ---
cat("--- 5. Gerando Box Plots Aninhados por Promotor... ---\n")

g_boxplot_aninhado <- ggplot(df_plot_aninhado, 
                             # Eixo X = Antagonista, Cor = Grupo
                             aes(x = Nome_Antagonista, y = Abundancia_VST, fill = Grupos)) +
  
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  # Jitter separado por cor
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), size = 1.5, alpha = 0.7) +
  
  # ✅ A CHAVE: Facet (Painel) por Nome_Promotor
  facet_wrap(~ Nome_Promotor, scales = "free_x") + 
  
  labs(
    title = "Abundância (VST) dos Antagonistas, Agrupados por Promotor (Alvo)",
    subtitle = "Hipótese: Antagonistas (Eixo X) são mais abundantes no NATB (Vermelho), inibindo o Promotor (Painel)",
    x = "Gênero Antagonista (Inibidor)",
    y = "Abundância Normalizada (VST)"
  ) +
  scale_fill_manual(values = c("NATB" = "red", "ATB" = "blue")) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "italic"),
        # Rotacionar nomes no eixo X se forem muitos
        axis.text.x = element_text(angle = 45, hjust = 1)) 

print(g_boxplot_aninhado)
ggsave("grafico_abundancia_antagonistas_ANINHADO_POR_PROMOTOR.png", g_boxplot_aninhado, width=12, height=10, dpi=300)

cat("\n--- ANÁLISE COMPLETA DE ANTAGONISTAS (ANINHADO) CONCLUÍDA ---\n")