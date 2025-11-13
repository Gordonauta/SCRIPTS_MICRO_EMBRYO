# =================================================================
# SCRIPT 1 (MODIFICADO): ANÁLISE CRUZADA (Interação vs. Abundância)
#
# OBJETIVO: Encontrar gêneros onde o IMPACTO é diferente (Interação)
# E a ABUNDÂNCIA é diferente (P-VALOR BRUTO < 0.05).
# =================================================================

cat("--- Iniciando Script de Análise Cruzada (MODIFICADO) ---\n")

# --- 0. CARREGAR BIBLIOTECAS ---
library(phyloseq)
library(tidyverse)
library(compositions) # Para CLR
library(DESeq2)       # Para Abundância Diferencial

# --- 1. IMPORTAÇÃO E PREPARAÇÃO DE DADOS (Baseado em 'analise-taxa-xlr.txt') ---

# 1.1. Importação Básica
biom_file <- "final.opti_mcc.filter.pick.biom"
metadata_file <- "metadata.txt"
ps <- import_biom(biom_file)
sample_data_df <- read.table(metadata_file, header = FALSE, row.names = 1, sep = "", stringsAsFactors = FALSE)
colnames(sample_data_df) <- c("Grupos")
sample_data_df$Grupos <- as.factor(sample_data_df$Grupos)
sample_data(ps) <- sample_data(sample_data_df)

# 1.2. Integrar Dados de Embriões
desenvolvimento_df <- data.frame(
  NAME = c("JC_1", "JC_10", "JC_11", "JC_12", "JC_13", "JC_14", "JC_2", "JC_3", "JC_4", "JC_5", "JC_6", "JC_7", "JC_8", "JC_9"),
  premorula = c(0, 0, 4, 0, 0, 0, 0, 0, 7, 0, 14, 0, 0, 1),
  morula = c(0, 0, 0, 15, 1, 4, 15, 11, 5, 0, 0, 0, 0, 1),
  blastocisto = c(12, 11, 5, 1, 5, 1, 0, 2, 0, 0, 0, 9, 7, 0),
  row.names = 1
)
desenvolvimento_df <- desenvolvimento_df %>%
  mutate(Total_Embrioes = premorula + morula + blastocisto)
novos_metadados <- sample_data(ps)
novos_metadados <- cbind(novos_metadados, desenvolvimento_df[rownames(novos_metadados), c("premorula", "morula", "blastocisto", "Total_Embrioes")])
sample_data(ps) <- sample_data(novos_metadados)

# 1.3. Filtros e Agrupamento (GLOM)
amostras_para_remover <- rownames(sample_data(ps))[sample_data(ps)$Total_Embrioes == 0]
ps_filtrado_desenv <- prune_samples(!sample_names(ps) %in% amostras_para_remover, ps)
ps_filtrado_taxa <- prune_taxa(taxa_sums(ps_filtrado_desenv) >= 1, ps_filtrado_desenv)
ps_genero_bruto <- tax_glom(ps_filtrado_taxa, taxrank = "Rank6", NArm = FALSE)

# IMPORTANTE: Usamos o mesmo objeto filtrado para AMBAS as análises
ps_genero_filtrado <- prune_taxa(
  taxa_sums(otu_table(ps_genero_bruto) > 0) >= 3,
  ps_genero_bruto
)

cat(paste("--- Objeto 'ps_genero_filtrado' pronto com", ntaxa(ps_genero_filtrado), "gêneros ---\n"))

# --- 2. ANÁLISE PARTE A: TESTE DE INTERAÇÃO (Impacto Diferencial) ---
# (Lógica de 'analise-diversidade-alfa-e-grupos-tratamento.txt') [cite_start][cite: 75-81]

cat("--- Iniciando Parte A: Teste de Interação (LM)... ---\n")

# 2.1. [cite_start]Preparar dados CLR (para o LM) [cite: 22-27]
otu_bruto_genero_pseudo <- as.data.frame(otu_table(ps_genero_filtrado)) + 1
clr_data_matrix <- t(clr(t(otu_bruto_genero_pseudo)))
taxa_abund_clr <- as.data.frame(clr_data_matrix)
metadados_filtrados <- as(sample_data(ps_genero_filtrado), "data.frame") %>%
  mutate(Taxa_Blastocisto = blastocisto / Total_Embrioes)

# 2.2. [cite_start]Loop do Modelo Linear de Interação [cite: 77-81]
generos_para_testar <- rownames(taxa_abund_clr)
resultados_lm_interacao <- data.frame(OTU_ID = character(), P_valor_Interacao = numeric())

for (genero_id in generos_para_testar) {
  df_modelo <- data.frame(
    Taxa_Blastocisto = metadados_filtrados$Taxa_Blastocisto,
    Grupos = metadados_filtrados$Grupos,
    Abundancia_CLR = as.numeric(taxa_abund_clr[genero_id, ])
  )
  
  modelo_lm <- try(lm(Taxa_Blastocisto ~ Abundancia_CLR * Grupos, data = df_modelo), silent = TRUE)
  
  if (!inherits(modelo_lm, "try-error")) {
    sumario_modelo <- summary(modelo_lm)
    if ("Abundancia_CLR:GruposNATB" %in% rownames(sumario_modelo$coefficients)) {
      p_valor_interacao <- sumario_modelo$coefficients["Abundancia_CLR:GruposNATB", "Pr(>|t|)"]
      resultados_lm_interacao <- rbind(resultados_lm_interacao, data.frame(
        OTU_ID = genero_id,
        P_valor_Interacao = p_valor_interacao
      ))
    }
  }
}

# 2.3. Filtrar resultados da Interação
df_interacao_sig <- resultados_lm_interacao %>%
  na.omit() %>%
  filter(P_valor_Interacao < 0.05) %>%
  arrange(P_valor_Interacao)

cat(paste("--- Parte A Concluída:", nrow(df_interacao_sig), "gêneros com Interação Significativa (P<0.05) ---\n"))
print("Gêneros-Chave (Interação):")
print(df_interacao_sig)


# --- 3. ANÁLISE PARTE B: TESTE DE ABUNDÂNCIA DIFERENCIAL (DESeq2) ---
# (Lógica de 'graficos.txt', mas adaptada para 'ps_genero_filtrado')

cat("\n--- Iniciando Parte B: Teste de Abundância (DESeq2)... ---\n")

# 3.1. Adicionar pseudocontagem para DESeq2
ps_genero_filtrado_plus1 <- ps_genero_filtrado
otu_table(ps_genero_filtrado_plus1) <- otu_table(ps_genero_filtrado_plus1) + 1

# 3.2. Rodar DESeq2
dds_obj <- phyloseq_to_deseq2(ps_genero_filtrado_plus1, design = ~ Grupos)
dds_result <- DESeq(dds_obj)
res_deseq <- results(dds_result, contrast = c("Grupos", "NATB", "ATB"))

# 3.3. ***** MUDANÇA ESTÁ AQUI *****
# Filtrar resultados do DESeq2 usando P-VALOR BRUTO
df_deseq_sig_bruto <- as.data.frame(res_deseq) %>%
  rownames_to_column(var = "OTU_ID") %>%
  filter(pvalue < 0.05) %>% # <-- MUDANÇA: de 'padj' para 'pvalue'
  select(OTU_ID, log2FoldChange, pvalue, padj) %>% # Selecionamos ambos os p-valores
  arrange(pvalue)

cat(paste("--- Parte B Concluída:", nrow(df_deseq_sig_bruto), "gêneros com Abundância Diferencial (P-Bruto < 0.05) ---\n"))
print("Gêneros Diferencialmente Abundantes (P-Bruto):")
print(df_deseq_sig_bruto)


# --- 4. ANÁLISE CRUZADA (O "JOIN") ---
cat("\n--- Iniciando Parte C: Análise Cruzada (O 'Cruzamento')... ---\n")

df_crossover <- inner_join(
  df_interacao_sig,     # Gêneros da Parte A (Interação)
  df_deseq_sig_bruto,   # Gêneros da Parte B (Abundância P-Bruto)
  by = "OTU_ID"
)

# 4.1. Adicionar Taxonomia para leitura
if (nrow(df_crossover) > 0) {
  df_crossover_taxa <- cbind(
    df_crossover,
    as(tax_table(ps_genero_filtrado)[df_crossover$OTU_ID, ], "matrix")
  )
  
  cat("\n--- 🎯 SUCESSO! Gêneros 'Duplamente Importantes' Encontrados (com P-Bruto): ---\n")
  cat("Estes gêneros MUDAM DE ABUNDÂNCIA (P-Bruto) e também MUDAM SEU IMPACTO no desfecho:\n")
  print(df_crossover_taxa)
  
} else {
  cat("\n--- Nenhum gênero 'Duplamente Importante' encontrado. ---\n")
  cat("Não houve sobreposição entre os gêneros com interação significativa e os com abundância diferencial (p-bruto).\n")
}

cat("\n--- Script de Análise Cruzada (MODIFICADO) Concluído ---\n")