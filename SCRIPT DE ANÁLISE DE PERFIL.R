# =================================================================
# SCRIPT DE ANÁLISE DE PERFIL (D5):
# "PERFIL DE SUCESSO" vs. "PERFIL DE ATRASO"
# =================================================================
# OBJETIVO:
# 1. Identificar o gênero mais associado ao ATRASO (Pré-Mórula).
# 2. Identificar os gêneros que formam o "Perfil de Atraso" (amigos do Vilão).
# 3. Identificar os gêneros que formam o "Perfil de Sucesso" (inimigos do Vilão).
# =================================================================

# --- 0. CARREGAR BIBLIOTECAS ---
library(phyloseq)
library(tidyverse)
library(compositions)
library(Hmisc) # Para rcorr (Rede de Correlação)
library(ggplot2)

cat("--- Bibliotecas Carregadas ---\n")

# --- 1. SETUP E DADOS (COM LÓGICA D5) ---
biom_file <- "final.opti_mcc.filter.pick.biom"
metadata_file <- "metadata.txt"

ps <- import_biom(biom_file) 
sample_data_df <- read.table(metadata_file, header = FALSE, row.names = 1, sep = "", stringsAsFactors = FALSE)
colnames(sample_data_df) <- c("Grupos")
sample_data_df$Grupos <- as.factor(sample_data_df$Grupos)
sample_data(ps) <- sample_data(sample_data_df)

# Dados de Embriões
desenvolvimento_df <- data.frame(
  NAME = c("JC_1", "JC_10", "JC_11", "JC_12", "JC_13", "JC_14", "JC_2", "JC_3", "JC_4", "JC_5", "JC_6", "JC_7", "JC_8", "JC_9"),
  premorula = c(0, 0, 4, 0, 0, 0, 0, 0, 7, 0, 14, 0, 0, 1),
  morula = c(12, 0, 0, 15, 1, 4, 15, 11, 5, 0, 0, 0, 0, 1),
  blastocisto = c(0, 11, 5, 1, 5, 1, 0, 2, 0, 0, 0, 9, 7, 0),
  row.names = 1
)

# Calcular Taxas D5
desenvolvimento_df <- desenvolvimento_df %>%
  mutate(
    Total_Embrioes = premorula + morula + blastocisto,
    Taxa_Atraso = premorula / Total_Embrioes, 
    Taxa_Sucesso = (morula + blastocisto) / Total_Embrioes
  )

# Integrar e Filtrar
novos_metadados <- sample_data(ps)
novos_metadados <- cbind(novos_metadados, desenvolvimento_df[rownames(novos_metadados), ])
sample_data(ps) <- sample_data(novos_metadados)
ps <- prune_samples(sample_data(ps)$Total_Embrioes > 0, ps)

# Filtros e CLR
ps_filtrado_taxa <- prune_taxa(taxa_sums(ps) >= 1, ps)
ps_genero <- tax_glom(ps_filtrado_taxa, taxrank = "Rank6", NArm = FALSE)
ps_genero_filtrado <- prune_taxa(taxa_sums(otu_table(ps_genero) > 0) >= 3, ps_genero)

otu_table_df <- as.data.frame(otu_table(ps_genero_filtrado))
taxa_abund_clr <- as.data.frame(t(clr(t(otu_table_df + 1))))
meta_final <- as(sample_data(ps_genero_filtrado), "data.frame")

# Função de Nomes
get_tax_name <- function(id, ps_obj) {
  tax <- as(tax_table(ps_obj)[id, ], "matrix")
  name <- tax[1, "Rank6"]
  if(is.na(name) | name == "" | grepl("unclassified", name)) name <- paste(tax[1, "Rank5"], "(Fam)")
  return(name)
}

cat("--- Dados D5 Prontos ---\n")


# =================================================================
# PASSO 1: IDENTIFICAR O GÊNERO "VILÃO" (Associado ao Atraso)
# =================================================================
cat("\n--- Passo 1: Identificando o Gênero 'Vilão' (Pior Impacto no Atraso) ---\n")

# Calcular correlação de Spearman de TODOS os gêneros vs. Atraso
cor_atraso <- apply(taxa_abund_clr, 1, function(x) {
  if (sd(x) == 0) return(NA)
  cor(x, meta_final$Taxa_Atraso, method = "spearman")
})

# O "Vilão" é quem tem o Rho POSITIVO mais alto com Atraso
vilao_id <- names(which.max(cor_atraso))
vilao_nome <- get_tax_name(vilao_id, ps_genero_filtrado)
vilao_rho <- max(cor_atraso, na.rm = TRUE)

cat(paste("--- Gênero 'Vilão' (Pior Impacto) Identificado:", vilao_nome, "(ID:", vilao_id, ") ---\n"))
cat(paste("--- Correlação com Atraso (Rho):", round(vilao_rho, 3), "---\n"))


# =================================================================
# PASSO 2: CONSTRUIR A REDE DE CORRELAÇÃO GERAL
# =================================================================
cat("\n--- Passo 2: Construindo Rede de Correlação Geral (Hmisc) ---\n")

# Usar Hmisc para calcular a matriz de correlação (Rho e P-valor)
cor_matrix_geral <- Hmisc::rcorr(as.matrix(t(taxa_abund_clr)), type = "spearman")

# "Aplanar" a matriz (pegar os pares)
flatten_cor_matrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    OTU_1 = rownames(cormat)[row(cormat)[ut]],
    OTU_2 = rownames(cormat)[col(cormat)[ut]],
    Rho = cormat[ut],
    P_valor = pmat[ut]
  )
}

rede_completa <- flatten_cor_matrix(cor_matrix_geral$r, cor_matrix_geral$P)
cat("--- Rede de Correlação Geral Calculada ---\n")


# =================================================================
# PASSO 3: IDENTIFICAR OS PERFIS DE SUCESSO E ATRASO
# =================================================================
cat(paste("\n--- Passo 3: Encontrando 'Amigos' e 'Inimigos' do Vilão:", vilao_nome, " ---\n"))

# Filtrar a rede para focar apenas no nosso Vilão
rede_vilao <- rede_completa %>%
  filter(OTU_1 == vilao_id | OTU_2 == vilao_id) %>%
  # Garantir que o Vilão esteja sempre na Coluna 1
  mutate(
    Parceiro_ID = ifelse(OTU_1 == vilao_id, OTU_2, OTU_1),
    Parceiro_Nome = sapply(Parceiro_ID, get_tax_name, ps_obj=ps_genero_filtrado)
  ) %>%
  dplyr::select(Parceiro_ID, Parceiro_Nome, Rho, P_valor) %>%
  filter(P_valor < 0.1) # Usamos p-valor mais leniente para perfis

# --- Perfil de ATRASO (Amigos do Vilão) ---
# Gêneros que co-ocorrem (Rho > 0.5) com o Vilão
perfil_atraso <- rede_vilao %>%
  filter(Rho > 0.5) %>%
  arrange(desc(Rho))

# --- Perfil de SUCESSO (Inimigos do Vilão) ---
# Gêneros que inibem (Rho < -0.5) o Vilão
perfil_sucesso <- rede_vilao %>%
  filter(Rho < -0.5) %>%
  arrange(Rho)


# =================================================================
# PASSO 4: VISUALIZAÇÃO DOS PERFIS
# =================================================================
cat("\n--- Passo 4: Gerando Gráfico dos Perfis ---\n")

# Juntar os dois perfis para o gráfico
df_plot_perfis <- rbind(
  perfil_atraso %>% mutate(Perfil = "Perfil de Atraso (Amigos do Vilão)"),
  perfil_sucesso %>% mutate(Perfil = "Perfil de Sucesso (Inimigos do Vilão)")
)

if(nrow(df_plot_perfis) > 0) {
  p_perfis <- ggplot(df_plot_perfis, aes(x = reorder(Parceiro_Nome, Rho), y = Rho, fill = Perfil)) +
    geom_col() +
    coord_flip() + # Barras horizontais
    geom_hline(yintercept = 0, color = "black") +
    scale_fill_manual(values = c("Perfil de Atraso (Amigos do Vilão)" = "#E41A1C", 
                                 "Perfil de Sucesso (Inimigos do Vilão)" = "#377EB8")) +
    labs(
      title = paste("Perfis de Sucesso e Atraso (Baseado no Vilão:", vilao_nome, ")"),
      subtitle = "Correlação (Rho) de outros gêneros com o principal indicador de atraso (Pré-Mórula)",
      x = "Gênero",
      y = "Força da Correlação (Spearman Rho)"
    ) +
    theme_bw(base_size = 14) +
    theme(legend.position = "bottom")
  
  print(p_perfis)
  
} else {
  cat("\nNenhum perfil forte (Rho > 0.5 ou < -0.5) foi encontrado para o gênero vilão.\n")
}

cat("\n--- Análise de Perfis D5 Concluída ---\n")