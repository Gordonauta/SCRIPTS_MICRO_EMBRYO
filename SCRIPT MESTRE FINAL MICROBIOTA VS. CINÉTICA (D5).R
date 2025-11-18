# =================================================================
# SCRIPT MESTRE FINAL: MICROBIOTA VS. CINÉTICA (D5)
# Versão: Visual Avançado (Gráficos Aninhados/Facetados)
# =================================================================
# OBJETIVO: Gerar gráficos complexos de Rede e Cinética com o estilo
# visual dos seus exemplos (Scatter plots e Boxplots aninhados).
#
# LÓGICA D5:
#   - VILÃO (Antagonista) = Associado à Taxa de Atraso (Pré-Mórula)
#   - MOCINHO (Promotor)  = Associado à Taxa de Sucesso (Mórula+Blast)
# =================================================================

# --- 0. CARREGAR BIBLIOTECAS ---
library(phyloseq)
library(tidyverse)
library(vegan)
library(compositions)
library(ggplot2)
library(Hmisc)
library(DESeq2)
library(ggrepel)

cat("--- Bibliotecas Carregadas ---\n")

# =================================================================
# 1. SETUP E DADOS (LÓGICA D5)
# =================================================================
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

# Calcular Variáveis D5
desenvolvimento_df <- desenvolvimento_df %>%
  mutate(
    Total_Embrioes = premorula + morula + blastocisto,
    Taxa_Premorula = premorula / Total_Embrioes, 
    Taxa_Morula = morula / Total_Embrioes,
    Taxa_Blastocisto = blastocisto / Total_Embrioes,
    Taxa_Atraso = Taxa_Premorula, 
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

cat("--- Dados Prontos ---\n")


# =================================================================
# 2. IDENTIFICAR O "VILÃO Nº 1" (ANTAGONISTA DO SUCESSO)
# =================================================================
# Vilão = Maior correlação positiva com Taxa_Atraso (Pré-Mórula)

cor_atraso <- apply(taxa_abund_clr, 1, function(x) cor(x, meta_final$Taxa_Atraso, method="spearman"))
vilao_id <- names(which.max(cor_atraso)) # O maior Rho positivo
vilao_nome <- get_tax_name(vilao_id, ps_genero_filtrado)

cat(paste("--- Vilão Identificado:", vilao_nome, "---\n"))


# =================================================================
# 3. CONSTRUIR A REDE DO VILÃO (AMIGOS E INIMIGOS)
# =================================================================
# Amigos (Mediadores) = Rho Positivo com Vilão
# Inimigos (Inibidores) = Rho Negativo com Vilão

cormat <- rcorr(as.matrix(t(taxa_abund_clr)), type="spearman")
rede_flat <- data.frame(
  OTU1 = rownames(cormat$r)[row(cormat$r)[upper.tri(cormat$r)]],
  OTU2 = colnames(cormat$r)[col(cormat$r)[upper.tri(cormat$r)]],
  Rho = cormat$r[upper.tri(cormat$r)],
  P = cormat$P[upper.tri(cormat$r)]
)

rede_vilao <- rede_flat %>% 
  filter(OTU1 == vilao_id | OTU2 == vilao_id) %>%
  mutate(Parceiro_ID = ifelse(OTU1 == vilao_id, OTU2, OTU1)) %>%
  filter(P < 0.05 & abs(Rho) > 0.4) %>% # Filtro de força
  mutate(
    Tipo = ifelse(Rho > 0, "Amigo (Mediador do Atraso)", "Inimigo (Inibidor do Atraso)"),
    Nome_Parceiro = sapply(Parceiro_ID, get_tax_name, ps_obj=ps_genero_filtrado)
  )

cat(paste("--- Rede do Vilão construída com", nrow(rede_vilao), "parceiros ---\n"))


# =================================================================
# GRÁFICO TIPO 1: SCATTER PLOTS ANINHADOS (FACET WRAP)
# (Estilo: grafico_correlacao_antagonistas_de_Azospirillum.jpg)
# =================================================================
cat("\n--- Gerando Gráfico 1: Correlação Aninhada (Vilão vs. Parceiros) ---\n")

if(nrow(rede_vilao) > 0) {
  
  # Preparar dados longos para ggplot
  # Queremos: X = Abundância Parceiro, Y = Abundância Vilão
  
  df_abund_parceiros <- taxa_abund_clr[rede_vilao$Parceiro_ID, ] %>%
    as.data.frame() %>% rownames_to_column("Parceiro_ID") %>%
    pivot_longer(-Parceiro_ID, names_to="Sample", values_to="Abund_Parceiro")
  
  df_abund_vilao <- data.frame(Sample = colnames(taxa_abund_clr), Abund_Vilao = as.numeric(taxa_abund_clr[vilao_id, ]))
  
  df_plot_scatter <- df_abund_parceiros %>%
    left_join(df_abund_vilao, by="Sample") %>%
    left_join(meta_final %>% rownames_to_column("Sample") %>% select(Sample, Grupos), by="Sample") %>%
    left_join(rede_vilao %>% select(Parceiro_ID, Nome_Parceiro, Tipo), by="Parceiro_ID")
  
  # Plotar
  g_scatter_rede <- ggplot(df_plot_scatter, aes(x=Abund_Parceiro, y=Abund_Vilao, color=Grupos)) +
    geom_point(size=3, alpha=0.7) +
    geom_smooth(method="lm", se=TRUE, fullrange=TRUE) +
    facet_wrap(~Nome_Parceiro, scales="free") +
    scale_color_manual(values=c("NATB"="red", "ATB"="blue")) +
    labs(
      title = paste("Rede do Vilão:", vilao_nome),
      subtitle = "Correlação com seus Parceiros (Amigos e Inimigos)",
      x = "Abundância (CLR) do Parceiro",
      y = paste("Abundância (CLR) de", vilao_nome)
    ) +
    theme_bw(base_size = 12) +
    theme(strip.text = element_text(face="bold.italic"))
  
  print(g_scatter_rede)
  ggsave(paste0("GRAFICO_1_SCATTER_REDE_", vilao_nome, ".png"), g_scatter_rede, width=12, height=10)
}


# =================================================================
# GRÁFICO TIPO 2: BOX PLOTS ANINHADOS (VST / ABUNDÂNCIA REAL)
# (Estilo: grafico_abundancia_mediadores_DE_Ruminococcus.jpg)
# =================================================================
cat("\n--- Gerando Gráfico 2: Boxplots Aninhados (Abundância Real) ---\n")

# Calcular VST (Abundância Normalizada para visualização limpa)
ps_plus1 <- ps_genero_filtrado; otu_table(ps_plus1) <- otu_table(ps_plus1) + 1 
dds <- phyloseq_to_deseq2(ps_plus1, ~ Grupos)
dds <- DESeq(dds, fitType='local')
vst_counts <- assay(varianceStabilizingTransformation(dds, blind=FALSE))

if(nrow(rede_vilao) > 0) {
  
  # Preparar dados longos de VST
  df_vst_parceiros <- vst_counts[rede_vilao$Parceiro_ID, ] %>%
    as.data.frame() %>% rownames_to_column("Parceiro_ID") %>%
    pivot_longer(-Parceiro_ID, names_to="Sample", values_to="Abundancia_VST") %>%
    left_join(meta_final %>% rownames_to_column("Sample") %>% select(Sample, Grupos), by="Sample") %>%
    left_join(rede_vilao %>% select(Parceiro_ID, Nome_Parceiro, Tipo), by="Parceiro_ID")
  
  # Plotar
  g_boxplot_rede <- ggplot(df_vst_parceiros, aes(x=Nome_Parceiro, y=Abundancia_VST, fill=Grupos)) +
    geom_boxplot(alpha=0.7, outlier.shape=NA) +
    geom_jitter(position=position_jitterdodge(jitter.width=0.2), size=1.5, alpha=0.6) +
    facet_wrap(~Tipo, scales="free_x") + # Separa Amigos de Inimigos
    scale_fill_manual(values=c("NATB"="red", "ATB"="blue")) +
    labs(
      title = paste("Abundância dos Parceiros de", vilao_nome),
      subtitle = "Comparação de abundância (VST) entre grupos",
      x = "Gênero Parceiro",
      y = "Abundância Normalizada (VST)"
    ) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle=45, hjust=1), strip.text = element_text(face="bold"))
  
  print(g_boxplot_rede)
  ggsave(paste0("GRAFICO_2_BOXPLOT_REDE_", vilao_nome, ".png"), g_boxplot_rede, width=12, height=8)
}


# =================================================================
# GRÁFICO TIPO 3: DISCREPÂNCIA CONDICIONAL (TAXA NO EIXO Y)
# (Estilo: PROVA_2_discrepancia_promotores_atb.jpg)
# =================================================================
cat("\n--- Gerando Gráfico 3: Discrepância Condicional (Vilão vs. Taxa) ---\n")

# Dados do Vilão vs Taxa de Atraso
df_disc_vilao <- data.frame(
  Abundancia = as.numeric(taxa_abund_clr[vilao_id, ]),
  Taxa_Atraso = meta_final$Taxa_Atraso,
  Grupos = meta_final$Grupos
)

g_discrepancia <- ggplot(df_disc_vilao, aes(x=Abundancia, y=Taxa_Atraso, color=Grupos)) +
  geom_point(size=4, alpha=0.8) + 
  geom_smooth(method="lm", se=FALSE, size=1.5) +
  scale_color_manual(values=c("NATB"="red", "ATB"="blue")) +
  labs(
    title = paste("Discrepância Condicional:", vilao_nome),
    subtitle = "O impacto do vilão no atraso muda com o antibiótico?",
    x = paste("Abundância (CLR) de", vilao_nome),
    y = "Taxa de Atraso (Pré-Mórula)"
  ) +
  theme_bw(base_size = 14)

print(g_discrepancia)
ggsave(paste0("GRAFICO_3_DISCREPANCIA_", vilao_nome, ".png"), g_discrepancia, width=8, height=6)

cat("\n--- Script Visual Avançado Concluído! ---\n")