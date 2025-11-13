# =================================================================
# SCRIPT: CURVA DE CINÉTICA POR RHO (COMPARATIVO ATB/NATB)
# =================================================================
# Objetivo: Plotar a curva de Correlação (Rho) pelos 3 estágios,
# comparando ATB vs. NATB, para os Top Promotores de ATB.

cat("--- Iniciando Geração da Curva Cinética de Rho (Comparativa)... ---\n")

# --- 0. CARREGAR BIBLIOTECAS ---
library(phyloseq)
library(tidyverse)
library(compositions)
library(ggplot2)
library(dplyr)
# --- 1. SETUP INICIAL E PREPARAÇÃO DOS DADOS (Necessário para criar os objetos) ---

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
metadados_filtrados <- as(sample_data(ps_genero_filtrado), "data.frame") %>%
  mutate(
    Taxa_Premorula = premorula / Total_Embrioes,
    Taxa_Morula = morula / Total_Embrioes,
    Taxa_Blastocisto = blastocisto / Total_Embrioes
  )

# Funções essenciais
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

# 1.3. Separar por Grupos
amostras_atb <- rownames(metadados_filtrados[metadados_filtrados$Grupos == "ATB", ])
amostras_natb <- rownames(metadados_filtrados[metadados_filtrados$Grupos == "NATB", ])
taxa_abund_atb <- taxa_abund_clr[, amostras_atb]
taxa_abund_natb <- taxa_abund_clr[, amostras_natb]


# --- 2. SELEÇÃO DE GÊNEROS (PROMOTORES ATB) ---
cat("--- 2. Selecionando Promotores ATB (Top Positivos vs. Taxa de Blastocisto)... ---\n")
TopN_select <- 12

blastocisto_taxa_atb <- metadados_filtrados[amostras_atb, "Taxa_Blastocisto"]
corr_atb_selecao <- run_spearman_analysis(blastocisto_taxa_atb, taxa_abund_atb)

colnames(corr_atb_selecao)<-c("Rho","P_valor_Bruto")

top_promotores_ids <- corr_atb_selecao %>%
  filter(Rho > 0) %>% 
  arrange(desc(Rho)) %>%
  head(TopN_select) %>%
  rownames()

cat(paste("---", length(top_promotores_ids), "Gêneros Promotores ATB Selecionados. ---\n"))

# --- 3. CÁLCULO DAS 6 CORRELAÇÕES CINÉTICAS ---
cat("--- 3. Calculando as 6 curvas de correlação (ATB/NATB vs Pré/Mórula/Blasto)... ---\n")

# A. Vetores de Taxa ATB
taxa_atb_pre <- metadados_filtrados[amostras_atb, "Taxa_Premorula"]
taxa_atb_mor <- metadados_filtrados[amostras_atb, "Taxa_Morula"]
taxa_atb_bla <- metadados_filtrados[amostras_atb, "Taxa_Blastocisto"]
# B. Vetores de Taxa NATB
taxa_natb_pre <- metadados_filtrados[amostras_natb, "Taxa_Premorula"]
taxa_natb_mor <- metadados_filtrados[amostras_natb, "Taxa_Morula"]
taxa_natb_bla <- metadados_filtrados[amostras_natb, "Taxa_Blastocisto"]

# Rodar as 6 análises
corr_atb_pre <- run_spearman_analysis(taxa_atb_pre, taxa_abund_atb)
corr_atb_mor <- run_spearman_analysis(taxa_atb_mor, taxa_abund_atb)
corr_atb_bla <- run_spearman_analysis(taxa_atb_bla, taxa_abund_atb)
corr_natb_pre <- run_spearman_analysis(taxa_natb_pre, taxa_abund_natb)
corr_natb_mor <- run_spearman_analysis(taxa_natb_mor, taxa_abund_natb)
corr_natb_bla <- run_spearman_analysis(taxa_natb_bla, taxa_abund_natb)

# --- 4. (CORRIGIDO) CRIAR O DATAFRAME "LONGO" PARA PLOTAGEM ---
cat("--- 4. (Corrigido) Criando Dataframe Longo para Plotagem... ---\n")

# A. Juntar todas as 6 tabelas de Rho
df_rho_total <- rbind(
  corr_atb_pre %>% rownames_to_column("OTU_ID") %>% mutate(Estagio = "Pré-Mórula", Grupo = "ATB"),
  corr_atb_mor %>% rownames_to_column("OTU_ID") %>% mutate(Estagio = "Mórula", Grupo = "ATB"),
  corr_atb_bla %>% rownames_to_column("OTU_ID") %>% mutate(Estagio = "Blastocisto", Grupo = "ATB"),
  
  corr_natb_pre %>% rownames_to_column("OTU_ID") %>% mutate(Estagio = "Pré-Mórula", Grupo = "NATB"),
  corr_natb_mor %>% rownames_to_column("OTU_ID") %>% mutate(Estagio = "Mórula", Grupo = "NATB"),
  corr_natb_bla %>% rownames_to_column("OTU_ID") %>% mutate(Estagio = "Blastocisto", Grupo = "NATB")
)
df_rho_total<- df_rho_total %>%
  rename(Rho.rho ="Rho")


# B. Criar a tabela de nomes (apenas para os IDs selecionados)
nomes_tax_df <- data.frame(
  OTU_ID = top_promotores_ids,
  Nome_Taxonomico = sapply(top_promotores_ids, get_tax_name, ps_obj = ps_genero_filtrado)
)

# C. Filtrar o Rho total APENAS para os promotores e adicionar os nomes
df_curva_rho_comparativa <- df_rho_total %>%
  filter(OTU_ID %in% top_promotores_ids) %>%
  left_join(nomes_tax_df, by = "OTU_ID") %>%
  # Remove Gêneros que não têm nome (caso get_tax_name falhe)
  na.omit() 

# D. Reordenar os estágios
df_curva_rho_comparativa$Estagio <- factor(df_curva_rho_comparativa$Estagio, 
                                           levels = c("Pré-Mórula", "Mórula", "Blastocisto"))

cat("--- Dataframe Longo criado com sucesso. ---\n")


# --- 5. GERAR O GRÁFICO DE CURVA CINÉTICA COMPARATIVA ---
# (Esta etapa permanece idêntica à do script anterior)
cat("--- 5. Gerando Gráfico de Discrepância Cinética (Rho)... ---\n")

g_cinetica_rho_comparativa <- ggplot(df_curva_rho_comparativa, 
                                     aes(x = Estagio, y = Rho, 
                                         group = Grupo, color = Grupo)) + 
  
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  
  facet_wrap(~ Nome_Taxonomico) + 
  
  labs(
    title = paste("Discrepância Cinética (Top", length(top_promotores_ids), "Promotores ATB)"),
    subtitle = "Comparação da Curva de Correlação (Rho) pelos Estágios (Taxa de Sucesso)",
    x = "Estágio de Desenvolvimento (Baseado na Taxa)",
    y = "Coeficiente de Correlação (Spearman's Rho)",
    color = "Grupo"
  ) +
  scale_color_manual(values = c("NATB" = "red", "ATB" = "blue")) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "italic")) 

print(g_cinetica_rho_comparativa)
ggsave("grafico_cinetica_rho_comparativa_atb_natb.png", g_cinetica_rho_comparativa, width=12, height=10, dpi=300)
