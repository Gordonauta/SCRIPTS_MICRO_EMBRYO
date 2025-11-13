# =================================================================
# SCRIPT 2 (CORRIGIDO): ANÁLISE DE REDE DE CORRELAÇÃO
# =================================================================
# OBJETIVO: Comparar as redes de correlação Gênero-vs-Gênero
# nos grupos ATB (Melhor Desfecho) e NATB (Pior Desfecho).
# =================================================================

cat("--- Iniciando Script de Análise de Rede (Ideia 2 - Corrigido) ---\n")

# --- 0. CARREGAR BIBLIOTECAS ---
library(phyloseq)
library(tidyverse)
library(compositions) # Para CLR
library(Hmisc)        # Para rcorr (cálculo de matriz de correlação)

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
ps_genero_filtrado <- prune_taxa(
  taxa_sums(otu_table(ps_genero_bruto) > 0) >= 3,
  ps_genero_bruto
)

# 1.4. [cite_start]Preparar dados CLR (para a correlação) [cite: 47-51]
otu_bruto_genero_pseudo <- as.data.frame(otu_table(ps_genero_filtrado)) + 1
clr_data_matrix <- t(clr(t(otu_bruto_genero_pseudo)))
taxa_abund_clr <- as.data.frame(clr_data_matrix)
metadados_filtrados <- as(sample_data(ps_genero_filtrado), "data.frame")

cat(paste("--- Dados CLR (para", ntaxa(ps_genero_filtrado), "gêneros) prontos ---\n"))

# 1.5. [cite_start]Função para buscar nomes de Gênero (de 'analise-taxa-xlr.txt') [cite: 55-57]
get_tax_name <- function(otu_id, ps_obj, rank_level = "Rank6") {
  tax_info <- as(tax_table(ps_obj)[otu_id, ], "matrix")
  name <- tax_info[1, rank_level]
  if (is.na(name) || name == "" || grepl("_unclassified", name) || grepl("_group", name)) {
    name_rank5 <- as.character(tax_info[1, "Rank5"])
    if (!is.na(name_rank5) && name_rank5 != "") {
      name <- paste0(name_rank5, " (Família)")
    } else {
      name <- otu_id
    }
  }
  return(as.character(name))
}


# --- 2. SEPARAR OS DADOS POR GRUPO ---
cat("--- Separando dados em ATB (Melhor Desfecho) e NATB (Pior Desfecho) ---\n")

amostras_atb <- rownames(metadados_filtrados[metadados_filtrados$Grupos == "ATB", ])
amostras_natb <- rownames(metadados_filtrados[metadados_filtrados$Grupos == "NATB", ])

# Matrizes de abundância CLR separadas
taxa_abund_atb <- taxa_abund_clr[, amostras_atb]
taxa_abund_natb <- taxa_abund_clr[, amostras_natb]

# --- 3. FUNÇÃO PARA RODAR E "ACHATAR" A MATRIZ DE CORRELAÇÃO ---
run_and_flatten_network <- function(taxa_abund_df, ps_obj) {
  
  taxa_abund_t <- t(taxa_abund_df)
  cor_matrix <- Hmisc::rcorr(as.matrix(taxa_abund_t), type = "spearman")
  
  df_long <- data.frame()
  generos <- rownames(cor_matrix$r)
  
  for (i in 1:(length(generos) - 1)) {
    for (j in (i + 1):length(generos)) {
      otu_1 <- generos[i]
      otu_2 <- generos[j]
      
      df_temp <- data.frame(
        OTU_1 = otu_1,
        OTU_2 = otu_2,
        Rho = cor_matrix$r[otu_1, otu_2],
        P_valor_Bruto = cor_matrix$P[otu_1, otu_2]
      )
      df_long <- rbind(df_long, df_temp)
    }
  }
  
  df_long <- df_long %>%
    na.omit() %>%
    mutate(P_valor_BH = p.adjust(P_valor_Bruto, method = "BH"))
  
  df_long$Genero_1 <- sapply(df_long$OTU_1, get_tax_name, ps_obj = ps_obj)
  df_long$Genero_2 <- sapply(df_long$OTU_2, get_tax_name, ps_obj = ps_obj)
  
  return(df_long)
}

# --- 4. EXECUTAR A ANÁLISE PARA AMBOS OS GRUPOS ---

cat("--- Calculando Rede de Correlação NATB (Pior Desfecho)... ---\n")
rede_natb_df <- run_and_flatten_network(taxa_abund_natb, ps_genero_filtrado)

cat("--- Calculando Rede de Correlação ATB (Melhor Desfecho)... ---\n")
rede_atb_df <- run_and_flatten_network(taxa_abund_atb, ps_genero_filtrado)

# --- 5. COMPARAR AS REDES (A "DESCOBERTA") ---
cat("--- Comparando as redes para encontrar correlações 'quebradas'... ---\n")

rede_natb_sig <- rede_natb_df %>% filter(P_valor_BH < 0.05)
rede_atb_sig <- rede_atb_df %>% filter(P_valor_BH < 0.05)

rede_comparativa <- full_join(
  rede_natb_sig,
  rede_atb_sig,
  by = c("Genero_1", "Genero_2"), # O par de gêneros
  suffix = c("_NATB", "_ATB")
)

rede_comparativa <- rede_comparativa %>%
  mutate(
    Status = case_when(
      !is.na(Rho_NATB) & is.na(Rho_ATB)    ~ "Apenas no NATB (Quebrada no ATB)",
      is.na(Rho_NATB)  & !is.na(Rho_ATB)   ~ "Apenas no ATB (Nova no ATB)",
      !is.na(Rho_NATB) & !is.na(Rho_ATB)   ~ "Presente em Ambos"
    )
  )

# --- 6. MOSTRAR RESULTADOS (TÍTULOS CORRIGIDOS) ---

# 6.1. Correlações do Pior Desfecho (NATB) que "sumiram" no ATB
# ***** CORREÇÃO APLICADA AQUI *****
rede_quebrada_no_atb <- rede_comparativa %>%
  filter(Status == "Apenas no NATB (Quebrada no ATB)") %>%
  arrange(Rho_NATB) %>%
  select(Genero_1, Genero_2, Rho_NATB, P_valor_BH_NATB, Status)

if (nrow(rede_quebrada_no_atb) > 0) {
  cat("\n\n--- 🎯 RESULTADO 1: Correlações 'Patogênicas' (NATB) que 'Quebraram' no ATB ---\n")
  cat("Estas relações existiam no Pior Desfecho (NATB) e sumiram no Melhor Desfecho (ATB):\n")
  print(head(rede_quebrada_no_atb, 20)) # Mostrar as 20 primeiras
} else {
  cat("\n\n--- Nenhuma correlação do NATB foi perdida no ATB ---\n")
}

# 6.2. Correlações "novas" que só existem no Melhor Desfecho (ATB)
# ***** CORREÇÃO APLICADA AQUI *****
rede_nova_no_atb <- rede_comparativa %>%
  filter(Status == "Apenas no ATB (Nova no ATB)") %>%
  arrange(Rho_ATB) %>%
  select(Genero_1, Genero_2, Rho_ATB, P_valor_BH_ATB, Status)

if (nrow(rede_nova_no_atb) > 0) {
  cat("\n\n--- 🎯 RESULTADO 2: Correlações 'Benéficas' (Novas) que Apareceram Apenas no ATB ---\n")
  cat("Estas relações só existem no Melhor Desfecho (ATB):\n")
  print(head(rede_nova_no_atb, 20)) # Mostrar as 20 primeiras
} else {
  cat("\n\n--- Nenhuma correlação nova apareceu no ATB ---\n")
}

cat("\n--- Script de Análise de Rede (Ideia 2 - Corrigido) Concluído ---\n")