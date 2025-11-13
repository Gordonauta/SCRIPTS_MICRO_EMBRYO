# =================================================================
# SCRIPT: CRUZAMENTO ML (VIP) vs. IMPACTO (L2FC e Rho) - V2 CORRIGIDO
# =================================================================
# (Correção aplicada na Etapa 3 para adicionar a coluna 'Rank_ML'
#  e corrigir o 'select()' )

cat("--- Iniciando Cruzamento: Biomarcadores (ML) vs. Impacto (DESeq2/Rho) ---\n")

# --- 0. CARREGAR BIBLIOTECAS ---
library(tidyverse)

# --- 1. VERIFICAR PRÉ-REQUISITOS (As tabelas devem existir no R) ---
if (!exists("tabela_biomarcadores_pls")) {
  stop("ERRO: 'tabela_biomarcadores_pls' não encontrada. Rode o SCRIPT_MACHINE_LEARNING (V5) primeiro.")
}
if (!exists("res_final_genero")) {
  stop("ERRO: 'res_final_genero' não encontrada. Rode o SCRIPT_2_Abundancia_Diferencial.R primeiro.")
}
if (!exists("corr_geral_selecao")) {
  stop("ERRO: 'corr_geral_selecao' não encontrada. Rode o SCRIPT_CURVA_CINÉTICA_RHO_IMPACTO_GERAL.R primeiro.")
}

cat("--- 1. Pré-requisitos encontrados. Preparando tabelas... ---\n")

# --- 2. PREPARAR TABELAS PARA O JOIN ---

# Tabela de Abundância Diferencial (DESeq2)
tabela_deseq2 <- res_final_genero %>%
  select(OTU_ID = feature_id, log2FoldChange, padj) %>%
  mutate(
    Impacto_DESeq2 = case_when(
      log2FoldChange > 0 ~ "Maléfico (NATB)",
      log2FoldChange < 0 ~ "Benéfico (ATB)",
      TRUE ~ "Neutro"
    )
  )

# Tabela de Correlação (Spearman vs. Taxa de Blastocisto Geral)
tabela_correlacao <- corr_geral_selecao %>%
  rownames_to_column("OTU_ID") %>%
  select(OTU_ID, Rho_Geral = Rho, P_Valor_Rho = P_valor_Bruto) %>%
  mutate(
    Impacto_Correlacao = case_when(
      Rho_Geral > 0 ~ "Benéfico (Positivo)",
      Rho_Geral < 0 ~ "Maléfico (Negativo)",
      TRUE ~ "Neutro"
    )
  )

# ======================================================
# --- 3. (CORRIGIDO) EXECUTAR O CRUZAMENTO (JOIN) ---
# ======================================================
cat("--- 2. Cruzando ML (VIP) com DESeq2 (L2FC) e Correlação (Rho)... ---\n")

tabela_impacto_final <- tabela_biomarcadores_pls %>%
  # ✅ CORREÇÃO 1: Adiciona o Rank (já que a tabela está ordenada por VIP)
  mutate(Rank_ML = row_number()) %>% 
  
  # Juntar com os resultados do DESeq2
  left_join(tabela_deseq2, by = "OTU_ID") %>%
  # Juntar com os resultados da Correlação
  left_join(tabela_correlacao, by = "OTU_ID") %>%
  
  # ✅ CORREÇÃO 2: Seleciona a nova coluna 'Rank_ML'
  select(
    Rank_ML, # (Corrigido de 'Rank_ML = Rank')
    Nome_Biomarcador,
    VIP_Score,
    Impacto_DESeq2,
    log2FoldChange,
    Impacto_Correlacao,
    Rho_Geral,
    OTU_ID
  ) %>%
  arrange(Rank_ML)

# --- 4. MOSTRAR A TABELA FINAL ---
cat("\n\n--- 🎯 TABELA DE IMPACTO DOS BIOMARCADORES DE MACHINE LEARNING ---\n")
print(tabela_impacto_final)

cat("\n--- SCRIPT DE CRUZAMENTO CONCLUÍDO ---\n")