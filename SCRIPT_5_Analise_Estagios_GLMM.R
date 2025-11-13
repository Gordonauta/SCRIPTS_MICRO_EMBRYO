# =================================================================
# SCRIPT 5: ANÁLISE DE ESTÁGIOS (GLMM)
# =================================================================

cat("--- Iniciando SCRIPT 5: Análise de Estágios (GLMM)... ---\n")

# --- 0. CARREGAR BIBLIOTECAS ---
library(phyloseq)
library(tidyverse)
library(compositions) 
library(lme4)     # Para 'glmer'
library(car)      # Para 'Anova()'
#install.packages("bayesm")

# --- 1. IMPORTAÇÃO E PREPARAÇÃO DOS DADOS ---
biom_file <- "final.opti_mcc.filter.pick.biom"
metadata_file <- "metadata.txt"
ps <- import_biom(biom_file)
sample_data_df <- read.table(metadata_file, header = FALSE, row.names = 1, sep = "", stringsAsFactors = FALSE)
colnames(sample_data_df) <- c("Grupos")
sample_data_df$Grupos <- as.factor(sample_data_df$Grupos)
sample_data(ps) <- sample_data(sample_data_df)

desenvolvimento_df <- data.frame(
  NAME = c("JC_1", "JC_10", "JC_11", "JC_12", "JC_13", "JC_14", "JC_2", "JC_3", "JC_4", "JC_5", "JC_6", "JC_7", "JC_8", "JC_9"),
  premorula = c(0, 0, 4, 0, 0, 0, 0, 0, 7, 0, 14, 0, 0, 1),
  morula = c(0, 0, 0, 15, 1, 4, 15, 11, 5, 0, 0, 0, 0, 1),
  blastocisto = c(12, 11, 5, 1, 5, 1, 0, 2, 0, 0, 0, 9, 7, 0),
  row.names = 1
)
desenvolvimento_df <- desenvolvimento_df %>%
  mutate(Total_Embrioes = premorula + morula + blastocisto)

# Combinar e salvar metadados no objeto 'ps'
novos_metadados <- sample_data(ps)
novos_metadados <- cbind(novos_metadados, desenvolvimento_df[rownames(novos_metadados), c("premorula", "morula", "blastocisto", "Total_Embrioes")])
sample_data(ps) <- sample_data(novos_metadados)

# Filtros e Glom
amostras_para_remover <- rownames(novos_metadados)[novos_metadados$Total_Embrioes == 0]
ps_filtrado_desenv <- prune_samples(!sample_names(ps) %in% amostras_para_remover, ps)
ps_filtrado_taxa <- prune_taxa(taxa_sums(ps_filtrado_desenv) >= 1, ps_filtrado_desenv)
ps_genero_bruto <- tax_glom(ps_filtrado_taxa, taxrank = "Rank6", NArm = FALSE)
ps_genero_filtrado <- prune_taxa(
  taxa_sums(otu_table(ps_genero_bruto) > 0) >= 3,
  ps_genero_bruto
)
colnames(tax_table(ps_genero_filtrado)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

# 1.4. Preparar dados CLR e Formato Longo
otu_bruto_genero_pseudo <- as.data.frame(otu_table(ps_genero_filtrado)) + 1
clr_data_matrix <- t(clr(t(otu_bruto_genero_pseudo)))
taxa_abund <- as.data.frame(clr_data_matrix)
metadados_filtrados_clr <- data.frame(sample_data(ps_genero_filtrado))

metadados_longo_base <- metadados_filtrados_clr %>%
  rownames_to_column("Amostra") %>% 
  pivot_longer(
    cols = c(premorula, morula, blastocisto), 
    names_to = "Estagio", 
    values_to = "Contagem"
  ) %>%
  mutate(
    Estagio = factor(Estagio, levels = c("premorula", "morula", "blastocisto")),
    Amostra = as.factor(Amostra)
  )
taxa_abund_t <- t(taxa_abund) %>% as.data.frame() %>% rownames_to_column("Amostra")
nomes_dos_generos <- rownames(taxa_abund)

cat("--- Dados prontos. Iniciando loop GLMM... ---\n")

# --- 2. ITERAR O MODELO GLMM EM TODOS OS GÊNEROS ---
lmm_resultados <- lapply(nomes_dos_generos, function(genero_nome) {
  abund_data_genero <- taxa_abund_t[, c("Amostra", genero_nome)]
  colnames(abund_data_genero)[2] <- "Abundancia"
  df_para_modelo <- left_join(metadados_longo_base, abund_data_genero, by = "Amostra")
  
  tryCatch({
    # Usamos glmer() com family = poisson(link="log")
    modelo_glmm <- glmer(Contagem ~ Abundancia * Estagio + Grupos + (1 | Amostra), 
                         data = df_para_modelo, 
                         family = poisson,
                         # Adicionar o controle para usar um otimizador robusto
                         control=glmerControl(optimizer="bobyqa", 
                                              optCtrl=list(maxfun=2e5)))
    # Usamos Anova() (do pacote 'car') para obter os p-valores
    anova_res <- car::Anova(modelo_glmm, type="II")
    
    # Extrair o P-valor da interação
    p_valor_interacao <- anova_res["Abundancia:Estagio", "Pr(>Chisq)"]
    return(p_valor_interacao)
    
  }, error = function(e) { return(NA) })
})
cat("--- Modelos GLMM executados ---\n")

# --- 3. COMPILAR E APRESENTAR OS RESULTADOS ---
resultados_interacao_df <- data.frame(
  OTU_ID = nomes_dos_generos, P_valor_Interacao = unlist(lmm_resultados)
)
resultados_interacao_df <- na.omit(resultados_interacao_df)
resultados_interacao_df$P_valor_BH <- p.adjust(
  resultados_interacao_df$P_valor_Interacao, method = "BH"
)
resultados_sig_interacao <- resultados_interacao_df %>%
  filter(P_valor_Interacao < 0.05) %>%
  arrange(P_valor_Interacao)

if (nrow(resultados_sig_interacao) > 0) {
  cat(paste("\n--- 🎯 GÊNEROS COM MUDANÇA DE CURVA (GLMM, P-Bruto < 0.05) ---\n"))
  resultados_sig_interacao_taxa <- cbind(
    resultados_sig_interacao,
    as(tax_table(ps_genero_filtrado)[resultados_sig_interacao$OTU_ID, ], "matrix")
  )
  print(resultados_sig_interacao_taxa)
} else {
  cat("\n--- Nenhum Gênero mostrou mudança de curva significativa (GLMM, P-Bruto < 0.05) ---\n")
}
cat("--- SCRIPT 5 Concluído ---\n")

# --- ADICIONE ESTA FUNÇÃO AO SEU SCRIPT (BASEADO EM ANALISE LMM.txt) ---
# Função para buscar nomes de Gênero
get_tax_name <- function(otu_id, ps_obj, rank_level = "Genus") {
  # CORREÇÃO: Usar Rank6 do seu SCRIPT 5 (que foi renomeado para "Genus")
  tax_info_matrix <- as(tax_table(ps_obj)[otu_id, ], "matrix")
  name <- tax_info_matrix[1, rank_level]
  if (is.na(name) || name == "" || grepl("_unclassified", name) || grepl("_group", name)) {
    # Usar a Família (Rank5)
    name_rank5 <- as.character(tax_info_matrix[1, "Family"]) 
    if (!is.na(name_rank5) && name_rank5 != "") {
      name <- paste0(name_rank5, " (Família)")
    } else {
      name <- otu_id
    }
  }
  return(as.character(name))
}

# --- SEÇÃO DE PLOTAGEM (BASEADA NA LÓGICA DO ANALISE LMM.txt) ---

# 1. VERIFICAR SE HÁ RESULTADOS SIGNIFICATIVOS PARA PLOTAR
if (nrow(resultados_sig_interacao) > 0) {
  
  cat(paste("\n--- Plotando curvas para os", nrow(resultados_sig_interacao), "gêneros significativos ---\n"))
  
  # Dataframe para armazenar os coeficientes
  df_curva_coef <- data.frame(
    OTU = character(),
    Nome_Taxonomico = character(),
    Estagio = character(),
    Coeficiente = numeric(),
    stringsAsFactors = FALSE
  )
  
  estagios_nomes <- c("premorula", "morula", "blastocisto")
  
  # Loop através dos Gêneros que foram significativos
  for (genero_nome in resultados_sig_interacao$OTU_ID) {
    
    # 2. Preparar os dados para este gênero
    abund_data_genero <- taxa_abund_t[, c("Amostra", genero_nome)]
    colnames(abund_data_genero)[2] <- "Abundancia"
    df_para_modelo <- left_join(metadados_longo_base, abund_data_genero, by = "Amostra")
    
    # 3. Rodar o GLMM novamente para obter os coeficientes
    # (Usando as mesmas configurações de controle para estabilidade)
    modelo <- tryCatch({
      glmer(Contagem ~ Abundancia * Estagio + Grupos + (1 | Amostra), 
            data = df_para_modelo, 
            family = poisson,
            control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
    }, error = function(e) { NULL })
    
    if (is.null(modelo)) next # Pula se o modelo falhar novamente
    
    # 4. Extrair os coeficientes dos Efeitos Fixos
    coefs <- summary(modelo)$coefficients
    
    # Nomes dos coeficientes de interação (podem variar dependendo do R)
    interacao_morula_nome <- "Abundancia:Estagiomorula"
    interacao_blasto_nome <- "Abundancia:Estagioblastocisto"
    
    # 5. Coeficiente Base (Impacto na Pré-Mórula)
    coef_base <- coefs["Abundancia", "Estimate"]
    
    # Coeficientes Interação (Mudança em relação à Pré-Mórula)
    coef_int_morula <- if (interacao_morula_nome %in% rownames(coefs)) coefs[interacao_morula_nome, "Estimate"] else 0
    coef_int_blasto <- if (interacao_blasto_nome %in% rownames(coefs)) coefs[interacao_blasto_nome, "Estimate"] else 0
    
    # 6. Calcular o impacto real em cada estágio
    impacto_premorula <- coef_base
    impacto_morula <- coef_base + coef_int_morula
    impacto_blastocisto <- coef_base + coef_int_blasto
    
    # 7. Obter o nome do Gênero
    nome_tax <- get_tax_name(genero_nome, ps_genero_filtrado)
    
    # 8. Adicionar ao dataframe de plotagem
    df_curva_coef <- rbind(df_curva_coef, 
                           data.frame(
                             OTU = genero_nome, Nome_Taxonomico = nome_tax, Estagio = estagios_nomes[1], Coeficiente = impacto_premorula, stringsAsFactors = FALSE
                           ),
                           data.frame(
                             OTU = genero_nome, Nome_Taxonomico = nome_tax, Estagio = estagios_nomes[2], Coeficiente = impacto_morula, stringsAsFactors = FALSE
                           ),
                           data.frame(
                             OTU = genero_nome, Nome_Taxonomico = nome_tax, Estagio = estagios_nomes[3], Coeficiente = impacto_blastocisto, stringsAsFactors = FALSE
                           )
    )
  }
  
  # 9. Reordenar os estágios
  df_curva_coef$Estagio <- factor(df_curva_coef$Estagio, levels = estagios_nomes)
  
  # 10. FILTRAR APENAS TOP 15 para um gráfico legível
  # Você tem 77, então vamos plotar os 15 mais significativos
  top_15_otus <- head(resultados_sig_interacao$OTU_ID, 15)
  df_plot_final <- df_curva_coef %>% filter(OTU %in% top_15_otus)
  
  # 11. GERAR O GRÁFICO DE CURVAS DE COEFICIENTES
  
  plot_coef_estagio_lmm <- ggplot(df_plot_final, 
                                  aes(x = Estagio, y = Coeficiente, 
                                      group = Nome_Taxonomico, 
                                      color = Nome_Taxonomico)) +
    
    geom_line(linewidth = 1.2) + 
    geom_point(size = 3) +     
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") + 
    labs(
      title = "Curva de Impacto da Microbiota no Desenv. Embrionário (GLMM)",
      subtitle = "Coeficiente (Impacto) dos Top 15 Gêneros-Chave por Estágio",
      x = "Estágio de Desenvolvimento",
      y = "Coeficiente (Impacto na Contagem de Embriões)",
      color = "Táxon (Gênero/Família)"
    ) +
    scale_color_brewer(palette = "Paired") + 
    scale_x_discrete(labels=c("Pré-Mórula", "Mórula", "Blastocisto")) +
    theme_bw(base_size = 14) +
    theme(
      legend.text = element_text(face = "italic"), 
      plot.title = element_text(face = "bold")
    )
  
  print(plot_coef_estagio_lmm)
  
} else {
  cat("\nAVISO: Nenhum Gênero com interação significativa (P<0.05) foi encontrado para plotar.\n")
}