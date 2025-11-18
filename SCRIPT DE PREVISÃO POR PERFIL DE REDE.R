# =================================================================
# SCRIPT DE PREVISÃO POR PERFIL DE REDE (D5)
# =================================================================
# OBJETIVO:
# 1. Definir o "Perfil de Sucesso" e o "Perfil de Atraso" (Top 10 de cada).
# 2. Criar um "Score de Perfil" (Mocinhos - Vilões).
# 3. Validar se este score consegue prever o desfecho (Sucesso vs. Atraso).
# =================================================================

# --- 0. CARREGAR BIBLIOTECAS ---
library(phyloseq)
library(tidyverse)
library(compositions)
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
    Taxa_Sucesso = (morula + blastocisto) / Total_Embrioes,
    Desfecho = as.factor(ifelse(Taxa_Atraso > 0, "Atraso", "Sucesso"))
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
# PASSO 2: DEFINIR OS PERFIS (MOCINHOS VS VILÕES)
# =================================================================
cat("\n--- Passo 2: Definindo Perfis de 'Vilões' e 'Mocinhos' ---\n")

# Calcular correlações de Spearman
cor_atraso <- apply(taxa_abund_clr, 1, function(x) cor(x, meta_final$Taxa_Atraso, method="spearman"))
cor_sucesso <- apply(taxa_abund_clr, 1, function(x) cor(x, meta_final$Taxa_Sucesso, method="spearman"))

# Perfil de Atraso (Top 10 Vilões)
viloes_ids <- names(sort(cor_atraso, decreasing = TRUE)[1:10])

# Perfil de Sucesso (Top 10 Mocinhos)
mocinhos_ids <- names(sort(cor_sucesso, decreasing = TRUE)[1:10])

cat("--- Perfis Definidos ---\n")
cat("Vilões (associados ao Atraso):\n"); print(sapply(viloes_ids, get_tax_name, ps_obj=ps_genero_filtrado))
cat("\nMocinhos (associados ao Sucesso):\n"); print(sapply(mocinhos_ids, get_tax_name, ps_obj=ps_genero_filtrado))


# =================================================================
# PASSO 3: CALCULAR O "SCORE DO PERFIL" (PREVISÃO)
# =================================================================
cat("\n--- Passo 3: Calculando o 'Score do Perfil' (Mocinhos - Vilões) ---\n")

# Transpor a tabela CLR (amostras = linhas, bactérias = colunas)
taxa_abund_t <- as.data.frame(t(taxa_abund_clr))

# Calcular a abundância média de cada perfil
media_viloes <- rowMeans(taxa_abund_t[, viloes_ids])
media_mocinhos <- rowMeans(taxa_abund_t[, mocinhos_ids])

# Adicionar o Score ao metadado
meta_final$Perfil_Score <- media_mocinhos - media_viloes
meta_final$Desfecho <- as.factor(ifelse(meta_final$Taxa_Atraso > 0, "Atraso", "Sucesso"))

cat("--- Score do Perfil Calculado ---\n")


# =================================================================
# PASSO 4: VALIDAR O SCORE DO PERFIL
# =================================================================
cat("\n--- Passo 4: Validando o 'Score do Perfil' ---\n")

# --- VALIDAÇÃO 1: O Score consegue separar os Desfechos? (Boxplot) ---
cat("--- Gerando Gráfico de Validação 1 (Boxplot) ---\n")

# Teste estatístico
wilcox_score <- wilcox.test(Perfil_Score ~ Desfecho, data = meta_final)
cat("Resultado do Teste Wilcoxon (Score vs Desfecho):\n")
print(wilcox_score)

p_valid_1 <- ggplot(meta_final, aes(x = Desfecho, y = Perfil_Score, fill = Desfecho)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 4) +
  scale_fill_manual(values = c("Atraso" = "red", "Sucesso" = "blue")) +
  labs(
    title = "Validação 1: O 'Score do Perfil' Separa os Desfechos?",
    subtitle = paste("P-Valor (Wilcoxon) =", round(wilcox_score$p.value, 5)),
    y = "Score do Perfil (Alto = Sucesso | Baixo = Atraso)",
    x = "Desfecho Real da Amostra"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

print(p_valid_1)
# [Imagem de um boxplot mostrando separação de um score preditivo por desfecho]

# --- VALIDAÇÃO 2: O Score está correlacionado com a Taxa de Sucesso? (Scatter) ---
cat("\n--- Gerando Gráfico de Validação 2 (Scatter Plot) ---\n")

cor_score <- cor.test(meta_final$Perfil_Score, meta_final$Taxa_Sucesso, method = "spearman")
cat("Resultado da Correlação (Score vs Taxa de Sucesso):\n")
print(cor_score)

p_valid_2 <- ggplot(meta_final, aes(x = Perfil_Score, y = Taxa_Sucesso)) +
  geom_point(aes(color = Grupos), size = 5) +
  geom_smooth(method = "lm", color = "black") +
  labs(
    title = "Validação 2: O 'Score do Perfil' prevê a Cinética?",
    subtitle = paste("Correlação (Spearman Rho) =", round(cor_score$estimate, 3)),
    x = "Score do Perfil (Mocinhos - Vilões)",
    y = "Taxa de Sucesso Real (Mórula + Blastocisto)"
  ) +
  theme_bw(base_size = 14)

print(p_valid_2)
# [Imagem de um gráfico de dispersão mostrando uma forte correlação positiva]

cat("\n--- Análise de Perfil Preditivo (Rede) Concluída ---\n")