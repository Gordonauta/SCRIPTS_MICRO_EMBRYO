# =================================================================
# SCRIPT COMPLETO: MICROBIOTA UTERINA VS. DESENVOLVIMENTO EMBRIONÁRIO
# Abordagem: Taxas de Desenvolvimento + Transformação CLR
# =================================================================

# -----------------------------------------------------------------
# 0. CARREGAR BIBLIOTECAS
# -----------------------------------------------------------------
# Se faltar alguma, instale com: install.packages("nome_do_pacote")
library(phyloseq)
library(tidyverse)
library(dplyr)
library(vegan)       # Para estimate_richness (Alfa) e adonis (Beta)
library(compositions)  # Para a transformação CLR
library(pheatmap)      # Para o heatmap
library(ggplot2)     # Para os gráficos

cat("--- Bibliotecas Carregadas ---\n")

# =================================================================
# 1. SETUP INICIAL E IMPORTAÇÃO
# =================================================================
# Defina os nomes dos seus arquivos
biom_file <- "final.opti_mcc.filter.pick.biom"
metadata_file <- "metadata.txt"
ps <- import_biom(biom_file) 

# Metadados
sample_data_df <- read.table(metadata_file, header = FALSE, row.names = 1, sep = "", stringsAsFactors = FALSE)
colnames(sample_data_df) <- c("Grupos")
sample_data_df$Grupos <- as.factor(sample_data_df$Grupos)
sample_data(ps) <- sample_data(sample_data_df)

cat("--- (Etapa 1) Objeto Phyloseq 'ps' criado ---\n")

# =================================================================
# 2. INTEGRAÇÃO (EMBRIÕES), FILTRAGEM E AGRUPAMENTO (GLOM)
# =================================================================
# 2.1. Dados de Desenvolvimento Embrionário
desenvolvimento_df <- data.frame(
  NAME = c("JC_1", "JC_10", "JC_11", "JC_12", "JC_13", "JC_14", "JC_2", "JC_3", "JC_4", "JC_5", "JC_6", "JC_7", "JC_8", "JC_9"),
  premorula = c(0, 0, 4, 0, 0, 0, 0, 0, 7, 0, 14, 0, 0, 1),
  morula = c(0, 0, 0, 15, 1, 4, 15, 11, 5, 0, 0, 0, 0, 1),
  blastocisto = c(12, 11, 5, 1, 5, 1, 0, 2, 0, 0, 0, 9, 7, 0),
  row.names = 1
)
desenvolvimento_df <- desenvolvimento_df %>%
  mutate(Total_Embrioes = premorula + morula + blastocisto)

# 2.2. Combinar metadados com dados de embriões
novos_metadados <- sample_data(ps)
novos_metadados <- cbind(novos_metadados, desenvolvimento_df[rownames(novos_metadados), c("premorula", "morula", "blastocisto", "Total_Embrioes")])
sample_data(ps) <- sample_data(novos_metadados)

# 2.3. Remover Amostras (JC_5, com 0 embriões)
amostras_para_remover <- rownames(sample_data(ps))[sample_data(ps)$Total_Embrioes == 0]
ps_filtrado_desenv <- prune_samples(
  !sample_names(ps) %in% amostras_para_remover,
  ps
)

# 2.4. Filtro suave de OTUs (manter todos que aparecem pelo menos 1 vez)
ps_filtrado_taxa <- prune_taxa(taxa_sums(ps_filtrado_desenv) >= 1, ps_filtrado_desenv)

# 2.5. AGRUPAR (GLOM) OS OTUS POR GÊNERO (Rank6)
ps_genero_bruto <- tax_glom(ps_filtrado_taxa, taxrank = "Rank6", NArm = FALSE) 
cat(paste("Número de Gêneros antes do filtro de prevalência:", ntaxa(ps_genero_bruto), "\n"))

# 2.6. FILTRO DE PREVALÊNCIA (Manter Gêneros presentes em >= 3 amostras)
ps_genero_filtrado <- prune_taxa(
  taxa_sums(otu_table(ps_genero_bruto) > 0) >= 3, 
  ps_genero_bruto
)
cat(paste("Número de Gêneros após filtro de prevalência (>= 3 amostras):", ntaxa(ps_genero_filtrado), "\n"))
cat("--- (Etapa 2) Filtros e 'ps_genero_filtrado' criados ---\n")

# =================================================================
# 3. PREPARAÇÃO DOS DADOS PARA CORRELAÇÃO (com CLR)
# =================================================================
# 3.1. Transformação dos Dados (CLR - Centered Log-Ratio)
# 3.1.1. Extrair contagens brutas (Nível Gênero)
otu_bruto_genero <- as.data.frame(otu_table(ps_genero_filtrado))

# 3.1.2. Adicionar Pseudocount (para lidar com Zeros)
otu_bruto_genero_pseudo <- otu_bruto_genero + 1

# 3.1.3. Aplicar a transformação CLR
# Queremos aplicar o CLR *por amostra* (coluna). 
# Precisamos transpor (t()), aplicar o clr() [que age por linha], e transpor de volta (t()).
clr_data_matrix <- t(clr(t(otu_bruto_genero_pseudo)))

# 3.1.4. Criar o dataframe 'taxa_abund' final
# Este é o dataframe que usaremos nas correlações
taxa_abund <- as.data.frame(clr_data_matrix)
cat("\n--- (Etapa 3.1) Abundância transformada com CLR (Centered Log-Ratio) ---\n")

# 3.2. Extrair vetores de TAXAS e metadados
metadados_filtrados <- as(sample_data(ps_genero_filtrado), "data.frame")

# Calcular as TAXAS (proporções)
metadados_filtrados <- metadados_filtrados %>%
  mutate(
    Taxa_Premorula = premorula / Total_Embrioes,
    Taxa_Morula = morula / Total_Embrioes,
    Taxa_Blastocisto = blastocisto / Total_Embrioes
  )

# Vetores para a correlação
premorula_taxa <- metadados_filtrados$Taxa_Premorula
morula_taxa <- metadados_filtrados$Taxa_Morula
blastocisto_taxa <- metadados_filtrados$Taxa_Blastocisto
grupos_fator <- metadados_filtrados$Grupos 
cat("--- (Etapa 3.2) Metadados e taxas de desenvolvimento prontos ---\n")


# =================================================================
# 4. FUNÇÃO DE CORRELAÇÃO DE SPEARMAN (PARA OS 3 ESTÁGIOS)
# =================================================================
run_spearman_analysis <- function(count_vector, taxa_abund_df, ps_obj, stage_name) {
  
  cat(paste("\n--- Iniciando Correlação de Spearman para:", stage_name, "--- \n"))
  
  correlacoes <- apply(taxa_abund_df, 1, function(x) {
    if (sd(x, na.rm = TRUE) == 0) return(list(estimate = NA, p.value = NA))
    suppressWarnings(cor.test(x, count_vector, method = "spearman"))
  })
  
  cor_resultados_matrix <- t(sapply(correlacoes, function(x) {
    c(Rho = x$estimate, P_valor_Bruto = x$p.value)
  }))
  
  corr_df <- as.data.frame(cor_resultados_matrix)
  colnames(corr_df) <- c("Rho", "P_valor_Bruto")
  corr_df <- na.omit(corr_df) 
  corr_df$P_valor_BH <- p.adjust(corr_df$P_valor_Bruto, method = "BH")
  corr_df_ordenado <- corr_df %>% arrange(P_valor_Bruto)
  
  # Filtro por P-valor Bruto < 0.05 (conforme solicitado)
  corr_significativos_bruto <- corr_df_ordenado %>%
    filter(P_valor_Bruto < 0.05)
  
  if (nrow(corr_significativos_bruto) > 0) {
    cat(paste("\n--- 🎯 RESULTADOS (P-valor Bruto < 0.05) para:", stage_name, "---\n"))
    corr_significativos_bruto_taxa <- cbind(
      corr_significativos_bruto,
      as(tax_table(ps_obj)[rownames(corr_significativos_bruto), ], "matrix")
    )
    print(corr_significativos_bruto_taxa)
  } else {
    cat(paste("\n--- Nenhum Gênero significativo (P-valor Bruto < 0.05) encontrado para:", stage_name, "---\n"))
  }
  
  return(corr_df_ordenado)
}
cat("--- (Etapa 4) Função 'run_spearman_analysis' definida ---\n")

# =================================================================
# 5. EXECUÇÃO DA ANÁLISE DE CORRELAÇÃO
# =================================================================
corr_results_premorula <- run_spearman_analysis(premorula_taxa, taxa_abund, ps_genero_filtrado, "Taxa de Pré-Mórula")
corr_results_morula <- run_spearman_analysis(morula_taxa, taxa_abund, ps_genero_filtrado, "Taxa de Mórula")
corr_results_blastocisto <- run_spearman_analysis(blastocisto_taxa, taxa_abund, ps_genero_filtrado, "Taxa de Blastocisto")

cat("\n--- (Etapa 5) Análise de Correlação de Spearman concluída ---\n")

# =================================================================
# 6. FUNÇÃO DE BUSCA DE GÊNERO (Rank6) [CORRIGIDA]
# =================================================================
get_tax_name <- function(otu_id, ps_obj, rank_level = "Rank6") {
  # Converte a linha da tabela de taxonomia em uma matriz de 1 linha
  tax_info <- as(tax_table(ps_obj)[otu_id, ], "matrix")
  
  # CORREÇÃO: Acessa a [linha 1, coluna "Rank6"]
  name <- tax_info[1, rank_level] 
  
  # Verifica se o nome do Gênero (Rank6) é problemático
  if (is.na(name) || name == "" || grepl("_unclassified", name) || grepl("_group", name)) {
    
    # CORREÇÃO: Acessa a [linha 1, coluna "Rank5"]
    name_rank5 <- as.character(tax_info[1, "Rank5"]) 
    
    # Se a Família (Rank5) existir, usa ela
    if (!is.na(name_rank5) && name_rank5 != "") {
      name <- paste0(name_rank5, " (Família)")
    } else {
      # Se tudo falhar, usa o ID do OTU
      name <- otu_id
    }
  }
  return(as.character(name))
}
cat("--- (Etapa 6) Função 'get_tax_name' (CORRIGIDA) definida ---\n")

# =================================================================
# 7. GRÁFICO DE CURVA DE CORRELAÇÃO (SPEARMAN RHO) POR ESTÁGIO
# =================================================================
# Coleta os 3 Gêneros com P-valor mais baixo em cada estágio
top_premorula <- rownames(head(corr_results_premorula[corr_results_premorula$P_valor_Bruto < 0.05, ], 5))
top_morula <- rownames(head(corr_results_morula[corr_results_morula$P_valor_Bruto < 0.05, ], 5))
top_blastocisto <- rownames(head(corr_results_blastocisto[corr_results_blastocisto$P_valor_Bruto < 0.05, ], 5))

otus_interesse_curva <- unique(c(top_premorula, top_morula, top_blastocisto))

if (length(otus_interesse_curva) > 0) {
  
  df_curva_rho <- data.frame(
    OTU = character(), Estagio = character(), Rho = numeric(), stringsAsFactors = FALSE
  )
  estagios_nomes <- c("Taxa de Pré-Mórula", "Taxa de Mórula", "Taxa de Blastocisto")
  
  for (otu in otus_interesse_curva) {
    df_curva_rho <- rbind(df_curva_rho, data.frame(
      OTU = otu, Estagio = estagios_nomes[1], Rho = corr_results_premorula[otu, "Rho"]
    ))
    df_curva_rho <- rbind(df_curva_rho, data.frame(
      OTU = otu, Estagio = estagios_nomes[2], Rho = corr_results_morula[otu, "Rho"]
    ))
    df_curva_rho <- rbind(df_curva_rho, data.frame(
      OTU = otu, Estagio = estagios_nomes[3], Rho = corr_results_blastocisto[otu, "Rho"]
    ))
  }
  
  df_curva_rho$Nome_Taxonomico <- sapply(df_curva_rho$OTU, get_tax_name, ps_obj = ps_genero_filtrado)
  df_curva_rho$Estagio <- factor(df_curva_rho$Estagio, levels = estagios_nomes)
  df_curva_rho <- na.omit(df_curva_rho)
  
  plot_rho_estagio <- ggplot(df_curva_rho, aes(x = Estagio, y = Rho, group = Nome_Taxonomico, color = Nome_Taxonomico)) +
    geom_line(linewidth = 1.2) + 
    geom_point(size = 3) +  
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") + 
    labs(
      title = "Mudança na Correlação (Spearman Rho) pela Taxa de Estágio",
      subtitle = "Correlação dos Gêneros com a Taxa de Desenvolvimento",
      x = "Taxa de Estágio de Desenvolvimento",
      y = "Coeficiente de Correlação (Spearman Rho)",
      color = "Táxon (Gênero/Família)"
    ) +
    scale_color_brewer(palette = "Paired") + 
    theme_bw(base_size = 14) +
    theme(
      legend.text = element_text(face = "italic"), 
      plot.title = element_text(face = "bold")
    )
  
  cat("--- (Etapa 7) Gerando Gráfico de Curva de Correlação ---\n")
  print(plot_rho_estagio)
  
} else {
  cat("\n--- (Etapa 7) AVISO: Nenhum Gênero significativo (P<0.05) foi encontrado para plotar o gráfico de curvas.\n")
}

# =================================================================
# 10. SCATTER PLOT DE VALIDAÇÃO (Ex: Prevotellaceae)
# =================================================================

# 1. Encontrar o OTU_ID do seu táxon de interesse (ex: a Prevotellaceae com Rho -0.8)
#    (Suponha que seja o 'Otu0156' do seu diagnóstico)
print(otus_interesse_curva)
otu_id_vilao <- "Otu0032" # <-- Troque pelo ID do OTU real

nome_vilao <- get_tax_name(otu_id_vilao, ps_genero_filtrado)

# --- Gráfico 1 (Contra Pré-Mórula) ---a
df_validacao_pre <- data.frame(
  Abundancia_CLR = as.numeric(taxa_abund[otu_id_vilao, ]),
  Taxa_Estagio = metadados_filtrados[colnames(taxa_abund), "Taxa_Premorula"],
  Amostra = colnames(taxa_abund)
)
ggplot(df_validacao_pre, aes(x = Abundancia_CLR, y = Taxa_Estagio)) +
  geom_point(size = 3) + geom_smooth(method = "lm") +
  labs(title = paste(nome_vilao, "vs. Taxa de Pré-Mórula"), 
       x = paste("Abundância (CLR) de", nome_vilao), y = "Taxa de Pré-Mórula") +
  theme_bw()


# --- Gráfico 2 (Contra Mórula) ---
df_validacao_mor <- data.frame(
  Abundancia_CLR = as.numeric(taxa_abund[otu_id_vilao, ]),
  Taxa_Estagio = metadados_filtrados[colnames(taxa_abund), "Taxa_Morula"],
  Amostra = colnames(taxa_abund)
)
ggplot(df_validacao_mor, aes(x = Abundancia_CLR, y = Taxa_Estagio)) +
  geom_point(size = 3) + geom_smooth(method = "lm") +
  labs(title = paste(nome_vilao, "vs. Taxa de Mórula"), 
       x = paste("Abundância (CLR) de", nome_vilao), y = "Taxa de Mórula") +
  theme_bw()

# --- Gráfico 3 (Contra Blastocisto) ---
df_validacao_bla <- data.frame(
  Abundancia_CLR = as.numeric(taxa_abund[otu_id_vilao, ]),
  Taxa_Estagio = metadados_filtrados[colnames(taxa_abund), "Taxa_Blastocisto"],
  Amostra = colnames(taxa_abund)
)
ggplot(df_validacao_bla, aes(x = Abundancia_CLR, y = Taxa_Estagio)) +
  geom_point(size = 3) + geom_smooth(method = "lm") +
  labs(title = paste(nome_vilao, "vs. Taxa de Blastocisto"), 
       x = paste("Abundância (CLR) de", nome_vilao), y = "Taxa de Blastocisto") +
  theme_bw()

# =================================================================
# 9. VISUALIZAÇÃO COM HEATMAP (CORRIGIDO 2: O Bug do Rowname)
# =================================================================
library(pheatmap)
library(ggplot2) 

# 1. Coletar todos os resultados de Rho
otus_para_heatmap <- rownames(corr_results_premorula) 

rho_matrix <- data.frame(
  Taxa_Premorula = corr_results_premorula[otus_para_heatmap, "Rho"],
  Taxa_Morula = corr_results_morula[otus_para_heatmap, "Rho"],
  Taxa_Blastocisto = corr_results_blastocisto[otus_para_heatmap, "Rho"]
)

# --- CORREÇÃO CRÍTICA APLICADA AQUI ---
# Os vetores de Rho perdem seus 'nomes' ao entrar no data.frame.
# Precisamos reatribuir os rownames corretos (os IDs dos OTUs).
rownames(rho_matrix) <- otus_para_heatmap
# --- FIM DA CORREÇÃO ---


# 2. Encontrar os 15 gêneros com MAIOR VARIÂNCIA
# (Estes são os mais "interessantes" porque mudam muito entre as amostras)
variances <- apply(taxa_abund, 1, var)
variances_sorted <- sort(variances, decreasing = TRUE)
top_15_otus <- names(variances_sorted[1:15])

cat("\n--- Heatmap: Top 15 Gêneros Mais Variáveis ---\n")
print(top_15_otus)

# 3. Filtrar a matriz para os 15 OTUs de interesse
# Agora esta linha vai funcionar, pois os rownames da rho_matrix são os IDs dos OTUs
rho_matrix_top15 <- rho_matrix[rownames(rho_matrix) %in% top_15_otus, ]

# Limpar NAs (se houver algum)
rho_matrix_top15 <- na.omit(rho_matrix_top15)

cat(paste("\nNúmero de táxons válidos (Top 15 sem NAs) para o heatmap:", nrow(rho_matrix_top15), "\n"))


# 4. Mapear os nomes de Gênero (usando sua função da Etapa 6)
if(nrow(rho_matrix_top15) > 1) { # Verifica se sobrou algo após o na.omit
  
  nomes_taxonomicos_raw <- sapply(
    rownames(rho_matrix_top15), 
    get_tax_name, 
    ps_obj = ps_genero_filtrado
  )
  nomes_taxonomicos_vec <- as.character(nomes_taxonomicos_raw)
  rownames(rho_matrix_top15) <- make.unique(nomes_taxonomicos_vec)
  
  # 5. Gerar o Heatmap
  pheatmap(
    rho_matrix_top15,
    main = "Heatmap de Correlação (Spearman Rho)",
    subtitle = "Exibindo os Gêneros mais variáveis",
    xlab = "Taxa de Estágio de Desenvolvimento",
    ylab = "Táxon (Gênero/Família)",
    cluster_rows = TRUE,  # Agrupa gêneros com padrões similares
    cluster_cols = FALSE, # Mantém a ordem dos estágios
    display_numbers = TRUE, 
    number_format = "%.2f",
    color = colorRampPalette(c("blue", "white", "red"))(100), 
    fontsize = 10,
    border_color = "grey60"
  )
} else {
  cat("\nAVISO (Heatmap): Não há táxons suficientes (mínimo 2) para plotar após a filtragem e remoção de NAs.\n")
}

# =================================================================
# 8. ANÁLISE DE DIVERSIDADE ALFA
# =================================================================
cat("\n\n--- Iniciando Análise de Diversidade Alfa ---\n")

# 8.1. Calcular a Diversidade Alfa (Shannon, Chao1, etc.)
# IMPORTANTE: A diversidade é calculada no nível de OTU.
# Por isso, usamos 'ps_filtrado_taxa' (ANTES do tax_glom por Gênero).
div_df <- estimate_richness(
  ps_filtrado_taxa, 
  measures = c("Observed", "Chao1", "Shannon", "Simpson")
)

# 8.2. Combinar Metadados com a Diversidade
# Usamos o dataframe 'metadados_filtrados' que já criamos na Etapa 3.2,
# pois ele já contém as colunas 'Taxa_Blastocisto', 'Taxa_Morula', etc.

# Garantir que ambos dataframes estejam na mesma ordem
metadados_com_div <- metadados_filtrados[rownames(div_df), ]
# Combinar os dois
metadados_com_div <- cbind(metadados_com_div, div_df)

cat("Metadados e Índices de Diversidade combinados:\n")
print(head(metadados_com_div[, c("Taxa_Blastocisto", "Shannon", "Chao1")]))


# 8.3. Correlacionar Diversidade vs. Taxa de Blastocisto
cat("\n--- Correlação (Spearman) Diversidade vs. Taxa de Blastocisto ---\n")

# Exemplo: Shannon vs. Taxa de Blastocisto
cor_shannon <- cor.test(
  metadados_com_div$Shannon, 
  metadados_com_div$Taxa_Blastocisto, 
  method = "spearman"
)
cat("\nResultado Shannon vs. Taxa de Blastocisto:\n")
print(cor_shannon)

# Exemplo: Riqueza (Chao1) vs. Taxa de Blastocisto
cor_chao1 <- cor.test(
  metadados_com_div$Chao1, 
  metadados_com_div$Taxa_Blastocisto, 
  method = "spearman"
)
cat("\nResultado Riqueza (Chao1) vs. Taxa de Blastocisto:\n")
print(cor_chao1)


# 8.4. Visualizar: Scatterplot (Gráfico de Dispersão)
library(ggplot2)

# Gráfico para Shannon
plot_shannon_vs_blast <- ggplot(metadados_com_div, aes(x = Shannon, y = Taxa_Blastocisto)) +
  geom_point(aes(color = Grupos), size = 4) + # Pontos coloridos pelo "Grupo" original
  geom_smooth(method = "lm", color = "black") + # Linha de tendência (regressão linear)
  labs(
    title = "Diversidade (Shannon) vs. Taxa de Blastocisto",
    subtitle = paste0(
      "Correlação de Spearman: Rho = ", round(cor_shannon$estimate, 3), 
      ", P-valor = ", round(cor_shannon$p.value, 4)
    ),
    x = "Índice de Diversidade Shannon",
    y = "Taxa de Blastocisto (Proporção)"
  ) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(face = "bold"))

print(plot_shannon_vs_blast)

# Gráfico para Riqueza (Chao1)
plot_chao1_vs_blast <- ggplot(metadados_com_div, aes(x = Chao1, y = Taxa_Blastocisto)) +
  geom_point(aes(color = Grupos), size = 4) +
  geom_smooth(method = "lm", color = "black") + 
  labs(
    title = "Riqueza (Chao1) vs. Taxa de Blastocisto",
    subtitle = paste0(
      "Correlação de Spearman: Rho = ", round(cor_chao1$estimate, 3), 
      ", P-valor = ", round(cor_chao1$p.value, 4)
    ),
    x = "Índice de Riqueza Chao1",
    y = "Taxa de Blastocisto (Proporção)"
  ) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(face = "bold"))

print(plot_chao1_vs_blast)

cat("\n--- Análise de Diversidade Alfa Concluída ---\n")
# =================================================================
# 11. ANÁLISE DE DIVERSIDADE BETA (CORRIGIDA)
# =================================================================
cat("\n\n--- Iniciando Análise de Diversidade Beta (PCoA e PERMANOVA) ---")

library(vegan)

if (exists("ps_filtrado_taxa") && exists("metadados_com_div")) {
  
  # 11.1. Calcular a Matriz de Distância (Bray-Curtis)
  dist_bray <- phyloseq::distance(ps_filtrado_taxa, method = "bray")
  
  # 11.2. Teste Estatístico (PERMANOVA / adonis)
  # Usamos os metadados do dataframe da Etapa 10
  metadata_beta <- metadados_com_div
  
  cat("\n--- 📊 Resultado do Teste PERMANOVA (adonis) --- \n")
  cat("Testando a diferença entre os 'Grupos' (ATB vs NATB)...\n")
  
  adonis_result <- adonis2(dist_bray ~ Grupos, data = metadata_beta)
  print(adonis_result)
  
  # 11.3. Calcular a Ordenação (PCoA)
  ord_pcoa_bray <- ordinate(ps_filtrado_taxa, "PCoA", "bray")
  
  
  # --- CORREÇÃO APLICADA AQUI ---
  # Precisamos adicionar nossos metadados (que contêm "Grupos" e 
  # "Taxa_Blastocisto") de volta ao objeto phyloseq ANTES de plotar.
  sample_data(ps_filtrado_taxa) <- sample_data(metadata_beta)
  # --- FIM DA CORREÇÃO ---
  
  
  # 11.4. Plotar o PCoA - Colorido por GRUPO (ATB vs NATB)
  cat("\nGerando Gráfico PCoA por Grupo...\n")
  
  p_valor_adonis <- adonis_result$`Pr(>F)`[1]
  
  # Este gráfico agora vai funcionar corretamente
  plot_pcoa_grupos <- plot_ordination(
    ps_filtrado_taxa,
    ord_pcoa_bray,
    color = "Grupos" 
  ) +
    geom_point(size = 5) +
    stat_ellipse(type = "t") + 
    labs(
      title = "PCoA da Diversidade Beta (Bray-Curtis)",
      subtitle = paste0("PERMANOVA (adonis) p-valor: ", round(p_valor_adonis, 4)),
      color = "Grupos"
    ) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(face = "bold"))
  
  print(plot_pcoa_grupos)
  
  
  # 11.5. Plotar o PCoA - Colorido pela TAXA DE BLASTOCISTO
  cat("\nGerando Gráfico PCoA por Taxa de Blastocisto...\n")
  
  # Este gráfico (o que falhou) agora vai funcionar
  plot_pcoa_blastocisto <- plot_ordination(
    ps_filtrado_taxa,
    ord_pcoa_bray,
    color = "Taxa_Blastocisto" # Agora ele encontrará esta coluna!
  ) +
    geom_point(size = 5) +
    scale_color_gradient(low = "red", high = "blue") + # Vermelho (baixa) -> Azul (alta)
    labs(
      title = "PCoA da Diversidade Beta (Bray-Curtis)",
      subtitle = "Pontos coloridos pela Taxa de Blastocisto",
      color = "Taxa de Blastocisto"
    ) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(face = "bold"))
  
  print(plot_pcoa_blastocisto)
  
} else {
  cat("\nERRO (Etapa 11): Objetos 'ps_filtrado_taxa' ou 'metadados_com_div' não encontrados. 
      Execute as etapas anteriores primeiro.\n")
}

cat("\n--- Análise de Diversidade Beta Concluída ---\n")