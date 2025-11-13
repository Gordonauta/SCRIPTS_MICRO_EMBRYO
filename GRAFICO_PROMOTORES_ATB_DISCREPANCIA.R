# =================================================================
# SCRIPT: PROMOTORES ATB - ANÁLISE DE DISCREPÂNCIA CONDICIONAL (TAXA NO EIXO Y)
# =================================================================
# Objetivo: Selecionar os TOP 12 Promotores ATB (Rho > 0 vs. Taxa)
# e plotar a Discrepância de suas curvas (slopes) usando a TAXA no Eixo Y.

cat("--- Iniciando Geração do Gráfico de Discrepância (Taxa no Y) ---\n")

# --- 0. CARREGAR BIBLIOTERARYS ---
library(phyloseq)
library(tidyverse)
library(compositions)
library(ggplot2)

# --- 1. SETUP INICIAL E PREPARAÇÃO DOS DADOS (Necessário para criar os objetos) ---

# 1.1. Importação Básica
biom_file <- "final.opti_mcc.filter.pick.biom"
metadata_file <- "metadata.txt"
ps <- import_biom(biom_file) 
# Tentando carregar os metadados (Este bloco é o que tem falhado na execução)
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
desenvolvimento_df <- desenvolvimento_df %>%
  mutate(Total_Embrioes = premorula + morula + blastocisto)

novos_metadados <- sample_data(ps)
novos_metadados <- cbind(novos_metadados, desenvolvimento_df[rownames(novos_metadados), c("premorula", "morula", "blastocisto", "Total_Embrioes")])
sample_data(ps) <- sample_data(novos_metadados)

amostras_para_remover <- rownames(sample_data(ps))[sample_data(ps)$Total_Embrioes == 0]
ps_filtrado_desenv <- prune_samples(!sample_names(ps) %in% amostras_para_remover, ps)
ps_filtrado_taxa <- prune_taxa(taxa_sums(ps_filtrado_desenv) >= 1, ps_filtrado_desenv)
ps_genero_bruto <- tax_glom(ps_filtrado_taxa, taxrank = "Rank6", NArm = FALSE)
ps_genero_filtrado <- prune_taxa(
  taxa_sums(otu_table(ps_genero_bruto) > 0) >= 3,
  ps_genero_bruto
)

# 1.3. Preparar CLR e Metadados
otu_bruto_genero_pseudo <- as.data.frame(otu_table(ps_genero_filtrado)) + 1
clr_data_matrix <- t(clr(t(otu_bruto_genero_pseudo)))
taxa_abund_clr <- as.data.frame(clr_data_matrix)
metadados_filtrados <- as(sample_data(ps_genero_filtrado), "data.frame") %>%
  mutate(Taxa_Blastocisto = blastocisto / Total_Embrioes)

# 1.4. Separar por Grupos
amostras_atb <- rownames(metadados_filtrados[metadados_filtrados$Grupos == "ATB", ])
taxa_abund_atb <- taxa_abund_clr[, amostras_atb]


# --- 2. FUNÇÕES ESSENCIAIS ---

# Função para rodar Correlação Spearman
run_spearman_analysis <- function(count_vector, taxa_abund_df) {
  correlacoes <- apply(taxa_abund_df, 1, function(x) {
    if (sd(x, na.rm = TRUE) == 0) return(list(estimate = NA, p.value = NA))
    suppressWarnings(cor.test(x, count_vector, method = "spearman"))
  })
  corr_df <- as.data.frame(t(sapply(correlacoes, function(x) {
    c(Rho = x$estimate, P_valor_Bruto = x$p.value)
  })))
  return(corr_df %>% na.omit())
}

# Função para buscar nomes de Gênero
get_tax_name <- function(otu_id, ps_obj, rank_level = "Rank6") {
  tax_info <- as(tax_table(ps_obj)[otu_id, ], "matrix")
  name <- tax_info[1, rank_level]
  if (is.na(name) || name == "" || grepl("_unclassified", name) || grepl("_group", name)) {
    name_rank5 <- try(as.character(tax_info[1, "Rank5"]), silent = TRUE)
    if (!inherits(name_rank5, "try-error") && !is.na(name_rank5) && name_rank5 != "") {
      name <- paste0(name_rank5, " (Família)")
    } else {
      name <- otu_id
    }
  }
  return(as.character(name))
}


# --- 3. SELEÇÃO DE GÊNEROS (PROMOTORES ATB) ---
cat("--- 3. Selecionando Promotores ATB (Top Positivos vs. Taxa de Blastocisto)... ---\n")
TopN_select <- 12

# A. Correlacionar APENAS o grupo ATB (Taxa de Blastocisto)
blastocisto_taxa_atb <- metadados_filtrados[amostras_atb, "Taxa_Blastocisto"]
corr_atb <- run_spearman_analysis(blastocisto_taxa_atb, taxa_abund_atb)
colnames(corr_atb)<-c("Rho","P_valor_Bruto")
# Gêneros Promotores ATB (Top Positivos)
top_plot_ids <- corr_atb %>%
  filter(Rho > 0) %>% # Apenas Promotores
  arrange(desc(Rho)) %>%
  head(TopN_select) %>%
  rownames()

cat(paste("---", length(top_plot_ids), "Gêneros Promotores ATB Selecionados. ---\n"))


# --- 4. CRIAR O DATAFRAME "LONGO" PARA PLOTAGEM ---
df_plot_total <- data.frame()
nomes_taxonomicos <- sapply(top_plot_ids, get_tax_name, ps_obj = ps_genero_filtrado)

for (otu_id in top_plot_ids) {
  df_temp <- data.frame(
    Nome_Taxonomico = nomes_taxonomicos[otu_id],
    # ✅ Taxa de Blastocisto para o Eixo Y
    Taxa_Blastocisto_Plot = metadados_filtrados$Taxa_Blastocisto, 
    Grupos = metadados_filtrados$Grupos,
    Abundancia_CLR = as.numeric(taxa_abund_clr[otu_id, ])
  )
  df_plot_total <- rbind(df_plot_total, df_temp)
}


# --- 5. GERAR O GRÁFICO DE DISCREPÂNCIA ---
cat("--- 5. Gerando Gráfico de Discrepância Condicional... ---\n")

plot_discrepancia <- ggplot(df_plot_total, 
                            # ✅ Usando a Taxa para o eixo Y
                            aes(x = Abundancia_CLR, y = Taxa_Blastocisto_Plot, 
                                color = Grupos)) + 
  
  geom_point(size = 3, alpha = 0.8) +
  
  # A prova da discrepância: as linhas têm inclinações muito diferentes (ou invertidas)
  geom_smooth(method = "lm", se = FALSE, fullrange = TRUE) +
  
  facet_wrap(~ Nome_Taxonomico, scales = "free_x") +
  
  labs(
    title = paste("Discrepância Condicional (Top", length(top_plot_ids), "Promotores ATB)"),
    subtitle = "A curva vermelha (NATB) deve ser mais plana/negativa, mostrando a falha de eficiência do Promotor ATB fora do seu ambiente.",
    x = "Abundância do Gênero (CLR)",
    # ✅ Eixo Y Corrigido para Taxa
    y = "Taxa de Blastocisto (Proporção de Sucesso)",
    color = "Grupo"
  ) +
  scale_color_manual(values = c("NATB" = "red", "ATB" = "blue")) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "italic")) 

print(plot_discrepancia)
# ggsave("grafico_promotores_atb_taxa_discrepancia.png", plot_discrepancia, width=12, height=10, dpi=300)