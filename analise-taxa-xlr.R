# =================================================================
# SCRIPT COMPLETO E UNIFICADO: MICROBIOTA VS. DESENVOLVIMENTO (D5)
# Gera TODOS os gráficos, incluindo Heatmap específico "Sucesso vs Atraso".
# =================================================================

# -----------------------------------------------------------------
# 0. CARREGAR BIBLIOTECAS
# -----------------------------------------------------------------
library(phyloseq)
library(tidyverse)
library(vegan)
library(compositions)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)

cat("--- Bibliotecas Carregadas ---\n")

# -----------------------------------------------------------------
# 1. SETUP, IMPORTAÇÃO E TRATAMENTO DE DADOS
# -----------------------------------------------------------------
biom_file <- "final.opti_mcc.filter.pick.biom"
metadata_file <- "metadata.txt"

# Importar
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

# Calcular TAXAS (Lógica D5)
desenvolvimento_df <- desenvolvimento_df %>%
  mutate(
    Total_Embrioes = premorula + morula + blastocisto,
    Taxa_Premorula = premorula / Total_Embrioes, 
    Taxa_Morula = morula / Total_Embrioes,       
    Taxa_Blastocisto = blastocisto / Total_Embrioes, 
    
    # VARIÁVEIS AGREGADAS (Para o Heatmap Binário)
    Taxa_Atraso = Taxa_Premorula,
    Taxa_Sucesso = (morula + blastocisto) / Total_Embrioes
  )

# Integrar e Filtrar
novos_metadados <- sample_data(ps)
novos_metadados <- cbind(novos_metadados, desenvolvimento_df[rownames(novos_metadados), ])
sample_data(ps) <- sample_data(novos_metadados)
ps <- prune_samples(sample_data(ps)$Total_Embrioes > 0, ps)

# Filtros de Taxonomia e Prevalência
ps_filtrado_taxa <- prune_taxa(taxa_sums(ps) >= 1, ps)
ps_genero <- tax_glom(ps_filtrado_taxa, taxrank = "Rank6", NArm = FALSE)
ps_genero_filtrado <- prune_taxa(taxa_sums(otu_table(ps_genero) > 0) >= 3, ps_genero)



# Transformação CLR
otu_table_df <- as.data.frame(otu_table(ps_genero_filtrado))
taxa_abund_clr <- as.data.frame(t(clr(t(otu_table_df + 1))))
meta_final <- as(sample_data(ps_genero_filtrado), "data.frame")

# Função auxiliar para nomes
get_tax_name <- function(id, ps_obj) {
  tax <- as(tax_table(ps_obj)[id, ], "matrix")
  name <- tax[1, "Rank6"]
  if(is.na(name) | name == "" | grepl("unclassified", name)) name <- paste(tax[1, "Rank5"], "(Fam)")
  return(name)
}

cat("--- Dados Prontos. Iniciando Gráficos ---\n")


# =================================================================
# PARTE A: GRÁFICOS DE CORRELAÇÃO (RHO)
# =================================================================

# Calcular Spearman
cor_pre <- apply(taxa_abund_clr, 1, function(x) cor(x, meta_final$Taxa_Premorula, method="spearman"))
cor_mor <- apply(taxa_abund_clr, 1, function(x) cor(x, meta_final$Taxa_Morula, method="spearman"))
cor_bla <- apply(taxa_abund_clr, 1, function(x) cor(x, meta_final$Taxa_Blastocisto, method="spearman"))

# Selecionar Top OTUs
top_mocinhos <- names(sort(cor_bla, decreasing = TRUE)[1:4])
top_viloes <- names(sort(cor_pre, decreasing = TRUE)[1:4])
otus_rho <- unique(c(top_mocinhos, top_viloes))

# --- GRÁFICO 1: Slope Chart (RHO) ---
cat("\n--- Gerando Gráfico 1: Slope Chart (Rho) ---\n")
df_slope_rho <- data.frame()
estagios <- c("1. Pré-Mórula", "2. Mórula", "3. Blastocisto")

for(otu in otus_rho) {
  nome <- get_tax_name(otu, ps_genero_filtrado)
  df_slope_rho <- rbind(df_slope_rho, data.frame(OTU=nome, Estagio=estagios[1], Rho=cor_pre[otu]))
  df_slope_rho <- rbind(df_slope_rho, data.frame(OTU=nome, Estagio=estagios[2], Rho=cor_mor[otu]))
  df_slope_rho <- rbind(df_slope_rho, data.frame(OTU=nome, Estagio=estagios[3], Rho=cor_bla[otu]))
}

p1 <- ggplot(df_slope_rho, aes(x = Estagio, y = Rho, group = OTU, color = OTU)) +
  geom_line(linewidth = 1.2, alpha = 0.8) +
  geom_point(size = 4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  labs(title = "1. Dinâmica de Correlação (Rho)",
       subtitle = "Subindo = Associação com Sucesso | Descendo = Associação com Atraso",
       y = "Correlação de Spearman (Rho)") +
  theme_bw(base_size = 14) + scale_color_brewer(palette = "Dark2")
print(p1)

# --- GRÁFICO 2: Heatmap Evolutivo (RHO) ---
cat("\n--- Gerando Gráfico 2: Heatmap Evolutivo (Rho) ---\n")
rho_matrix <- data.frame(Pre_Morula = cor_pre[otus_rho], Morula = cor_mor[otus_rho], Blastocisto = cor_bla[otus_rho])
rownames(rho_matrix) <- sapply(rownames(rho_matrix), get_tax_name, ps_obj=ps_genero_filtrado)

pheatmap(rho_matrix, main = "2. Heatmap de Correlação (Evolução)", cluster_cols = FALSE, display_numbers = TRUE,
         color = colorRampPalette(c("red", "white", "blue"))(100))

# --- GRÁFICO 3: Scatter Plot do Vilão (RHO) ---
cat("\n--- Gerando Gráfico 3: Scatter Plot Validação ---\n")
vilao_id <- names(which.max(cor_pre))
vilao_nome <- get_tax_name(vilao_id, ps_genero_filtrado)
df_val <- data.frame(Abundancia = as.numeric(taxa_abund_clr[vilao_id, ]), Atraso = meta_final$Taxa_Atraso, Grupos = meta_final$Grupos)

p3 <- ggplot(df_val, aes(x = Abundancia, y = Atraso)) +
  geom_point(aes(color = Grupos), size = 5) +
  geom_smooth(method = "lm", color = "darkred", fill = "pink") +
  labs(title = paste("3. O Principal Vilão:", vilao_nome), subtitle = "Correlação Rho com Atraso", x = "Abundância (CLR)", y = "Taxa de Atraso") +
  theme_bw(base_size = 14)
print(p3)


# =================================================================
# PARTE B: GRÁFICOS DE REGRESSÃO (BETA) - IMPACTO
# =================================================================

get_beta <- function(abundance_vector, outcome_vector) {
  if(sd(abundance_vector) == 0) return(0)
  return(coef(lm(outcome_vector ~ abundance_vector))[2])
}

# Calcular Betas para os Estágios
beta_pre <- apply(taxa_abund_clr, 1, get_beta, outcome_vector = meta_final$Taxa_Premorula)
beta_mor <- apply(taxa_abund_clr, 1, get_beta, outcome_vector = meta_final$Taxa_Morula)
beta_bla <- apply(taxa_abund_clr, 1, get_beta, outcome_vector = meta_final$Taxa_Blastocisto)

# --- NOVO CÁLCULO: Betas para SUCESSO vs ATRASO (Agregado) ---
beta_sucesso_agregado <- apply(taxa_abund_clr, 1, get_beta, outcome_vector = meta_final$Taxa_Sucesso)
beta_atraso_agregado  <- apply(taxa_abund_clr, 1, get_beta, outcome_vector = meta_final$Taxa_Atraso)

# Top OTUs (Impacto)
otus_beta <- unique(c(names(sort(abs(beta_pre), decreasing=T)[1:4]), names(sort(abs(beta_bla), decreasing=T)[1:4])))

# --- GRÁFICO 4: Slope Chart (BETA) ---
cat("\n--- Gerando Gráfico 4: Slope Chart (Beta) ---\n")
df_slope_beta <- data.frame()
for(otu in otus_beta) {
  nome <- get_tax_name(otu, ps_genero_filtrado)
  df_slope_beta <- rbind(df_slope_beta, data.frame(OTU=nome, Estagio=estagios[1], Beta=beta_pre[otu]))
  df_slope_beta <- rbind(df_slope_beta, data.frame(OTU=nome, Estagio=estagios[2], Beta=beta_mor[otu]))
  df_slope_beta <- rbind(df_slope_beta, data.frame(OTU=nome, Estagio=estagios[3], Beta=beta_bla[otu]))
}

p4 <- ggplot(df_slope_beta, aes(x = Estagio, y = Beta, group = OTU, color = OTU)) +
  geom_line(linewidth = 1.2, alpha = 0.8) + geom_point(size = 4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  labs(title = "4. Dinâmica de Impacto (Beta / Regressão)", subtitle = "Eixo Y = Coeficiente de Regressão (Impacto)", y = "Beta (Slope)") +
  theme_bw(base_size = 14) + scale_color_brewer(palette = "Paired")
print(p4)

# --- GRÁFICO 5: Heatmap Evolutivo (BETA) ---
cat("\n--- Gerando Gráfico 5: Heatmap Evolutivo (Beta) ---\n")
beta_matrix <- data.frame(Pre_Morula = beta_pre[otus_beta], Morula = beta_mor[otus_beta], Blastocisto = beta_bla[otus_beta])
rownames(beta_matrix) <- sapply(rownames(beta_matrix), get_tax_name, ps_obj=ps_genero_filtrado)

pheatmap(beta_matrix, main = "5. Heatmap de Impacto (Evolução)", cluster_cols = FALSE, display_numbers = TRUE, number_format = "%.3f",
         color = colorRampPalette(c("blue", "white", "red"))(100))

# --- GRÁFICO 5.1 (NOVO): Heatmap CONTRASTE (Sucesso vs Atraso) ---
cat("\n--- Gerando Gráfico 5.1: Heatmap Contraste (Sucesso vs Atraso) ---\n")

# Matriz Binária
beta_contrast_matrix <- data.frame(
  Atraso = beta_atraso_agregado[otus_beta],
  Sucesso = beta_sucesso_agregado[otus_beta]
)
rownames(beta_contrast_matrix) <- sapply(rownames(beta_contrast_matrix), get_tax_name, ps_obj=ps_genero_filtrado)

pheatmap(beta_contrast_matrix, 
         main = "5.1 Heatmap de Contraste: Quem ajuda vs Quem atrapalha?", 
         cluster_cols = FALSE, 
         display_numbers = TRUE, 
         number_format = "%.3f",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         fontsize_number = 12)


# =================================================================
# PARTE C: ECOLOGIA E COMPARAÇÃO
# =================================================================

# --- GRÁFICO 6: Interação (O Gráfico em X) ---
cat("\n--- Gerando Gráfico 6: Interação ---\n")
meta_final$Chao1 <- estimate_richness(ps_filtrado_taxa, measures = "Chao1")$Chao1
# Para contar GÊNEROS especificamente:
meta_final$Riqueza_Generos <- estimate_richness(ps_genero, measures = "Observed")$Observed


# Estatísticas no Console
cat("\n>>> Teste de Interação (Chao1 vs Sucesso) <<<\n")
atb_d <- meta_final %>% filter(Grupos=="ATB"); natb_d <- meta_final %>% filter(Grupos=="NATB")
if(sd(atb_d$Taxa_Sucesso)>0) print(cor.test(atb_d$Chao1, atb_d$Taxa_Sucesso, method="spearman"))
if(sd(natb_d$Taxa_Sucesso)>0) print(cor.test(natb_d$Chao1, natb_d$Taxa_Sucesso, method="spearman"))

p6 <- ggplot(meta_final, aes(x = Chao1, y = Taxa_Sucesso, color = Grupos)) +
  geom_point(size = 5, alpha = 0.8) + geom_smooth(method = "lm", se = FALSE, linewidth = 2) + 
  labs(title = "6. Interação: Riqueza vs. Sucesso", subtitle = "O antibiótico inverteu o papel da riqueza?", x = "Riqueza (Chao1)", y = "Taxa de Sucesso") +
  theme_bw(base_size = 14) + theme(legend.position = "top")
print(p6)

# --- GRÁFICOS 7 e 8: PCoA ---
cat("\n--- Gerando Gráficos 7 e 8: PCoA ---\n")
dist_bray <- phyloseq::distance(ps_filtrado_taxa, method = "bray")
ord_pcoa <- ordinate(ps_filtrado_taxa, "PCoA", "bray")
adonis_res <- adonis2(dist_bray ~ Grupos, data = meta_final)
print(adonis_res) 

p7 <- plot_ordination(ps_filtrado_taxa, ord_pcoa, color = "Grupos") + geom_point(size = 5) + stat_ellipse() +
  labs(title = "7. PCoA: Grupos (ATB vs NATB)") + theme_bw()
print(p7)

sample_data(ps_filtrado_taxa)$Taxa_Premorula <- meta_final$Taxa_Premorula
p8 <- plot_ordination(ps_filtrado_taxa, ord_pcoa, color = "Taxa_Premorula") + geom_point(size = 6) +
  scale_color_viridis_c(option = "magma", direction = -1, name = "Atraso") + 
  labs(title = "8. PCoA: Gradiente de Atraso", subtitle = "Cor Forte = Mais Pré-Mórulas") + theme_bw() + theme(panel.background = element_rect(fill = "gray90"))
print(p8)

# --- GRÁFICO 9: Boxplots Comparativos ---
cat("\n--- Gerando Gráfico 9: Painel de Boxplots ---\n")
# Estatísticas no Console
cat("\n>>> Testes de Wilcoxon (Comparação de Grupos) <<<\n")
print(wilcox.test(Chao1 ~ Grupos, data = meta_final))
print(wilcox.test(Taxa_Premorula ~ Grupos, data = meta_final))
print(wilcox.test(Taxa_Sucesso ~ Grupos, data = meta_final))

df_long <- meta_final %>% select(Grupos, Taxa_Premorula, Taxa_Sucesso, Chao1) %>%
  pivot_longer(cols = c(Taxa_Premorula, Taxa_Sucesso, Chao1), names_to = "Variavel", values_to = "Valor")

p9 <- ggplot(df_long, aes(x = Grupos, y = Valor, fill = Grupos)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) + geom_jitter(width = 0.2, size = 2) +
  facet_wrap(~Variavel, scales = "free_y", labeller = as_labeller(c("Chao1"="Riqueza (Chao1)", "Taxa_Premorula"="Atraso (Pré-Mórula)", "Taxa_Sucesso"="Sucesso (Mórula+Blast)"))) +
  labs(title = "9. Comparação Geral: Grupos ATB vs NATB", y = "Valor") +
  theme_bw(base_size = 14) + theme(legend.position = "none", strip.text = element_text(face="bold"))
print(p9)

cat("\n--- TODOS OS GRÁFICOS GERADOS COM SUCESSO ---\n")