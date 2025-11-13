# =================================================================
# SCRIPT DE GRÁFICOS PARA AS PROVAS FINAIS (LOG2FC E INTERAÇÃO)
# =================================================================
# Este script assume que os objetos 'res_final_genero', 
# 'taxa_abund_clr' e 'metadados_filtrados' já foram criados
# pelos seus scripts de análise (DESeq2 e Análise Cruzada).

cat("--- Iniciando Geração de Gráficos de Prova ---\n")

# --- 0. CARREGAR BIBLIOTECAS ---
library(phyloseq)
library(tidyverse)
library(DESeq2)
library(ggrepel)

# =================================================================
# GRÁFICO 1: BAR PLOT DE LOG2FC (PROVA DE SUPERIORIDADE DO ATB)
# =================================================================

if (exists("res_final_genero")) {
  cat("\n--- Gerando Bar Plot de Log2FC (DESeq2, padj < 0.05) ---\n")
  
  # 1. Preparar a tabela (Filtrar por padj < 0.05)
  plot_data_barplot <- res_final_genero %>%
    filter(padj < 0.05) %>%
    # Criar a coluna de cor (Log2FC > 0 = NATB/Maléfico; Log2FC < 0 = ATB/Benéfico)
    mutate(
      grupo_enriquecido = if_else(log2FoldChange > 0, 
                                  "NATB (Pior Desfecho - Maléfico)", 
                                  "ATB (Melhor Desfecho - Benéfico)")
    ) %>%
    distinct(Genus, .keep_all = TRUE)
  
  # 2. Gerar o Gráfico
  if (nrow(plot_data_barplot) > 0) {
    
    g_barplot_log2fc <- ggplot(plot_data_barplot, 
                               aes(x = log2FoldChange, 
                                   # Ordena as barras pelo Log2FC
                                   y = reorder(Genus, log2FoldChange), 
                                   fill = grupo_enriquecido)) +
      geom_col() +
      scale_fill_manual(values = c("NATB (Pior Desfecho - Maléfico)" = "red", 
                                   "ATB (Melhor Desfecho - Benéfico)" = "blue")) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
      labs(title = "Gêneros Diferencialmente Abundantes (DESeq2, padj < 0.05)",
           subtitle = "Log2FC > 0: Mais abundante no Pior Desfecho (NATB)",
           x = "Log2 Fold Change (NATB vs. ATB)",
           y = "Gênero",
           fill = "Grupo Enriquecido (Impacto)") +
      theme_minimal() +
      theme(panel.grid.major.y = element_blank(),
            plot.title = element_text(face = "bold"))
    
    print(g_barplot_log2fc)
    
  } else {
    cat("!!! NENHUM GÊNERO SIGNIFICATIVO (padj < 0.05) encontrado para o Bar Plot de Log2FC.\n")
  }
} else {
  cat("AVISO: Objeto 'res_final_genero' não encontrado. Rode SCRIPT_2_Abundancia_Diferencial.R primeiro.\n")
}

# =================================================================
# GRÁFICO 2: INTERAÇÃO SCATTER PLOT (PROVA DE EFEITO CONDICIONAL)
# =================================================================
# Visualiza o gênero mais maléfico (maior Log2FC positivo)
# para demonstrar graficamente o efeito de interação.

if (exists("taxa_abund_clr") && exists("metadados_filtrados") && exists("res_final_genero")) {
  
  cat("\n--- Gerando Scatter Plot de Interação Condicional ---\n")
  
  # 1. Encontrar o Gênero mais Maléfico (maior Log2FC positivo, padj < 0.05)
  top_malefico <- res_final_genero %>%
    filter(padj < 0.05 & log2FoldChange > 0) %>%
    arrange(desc(log2FoldChange)) %>%
    head(1)
  
  if (nrow(top_malefico) > 0) {
    otu_id_malefico <- top_malefico$feature_id
    nome_malefico <- top_malefico$Genus
    
    # Adicionar taxa de Blastocisto ao metadados (se ainda não existir)
    if (!"Taxa_Blastocisto" %in% colnames(metadados_filtrados)) {
      metadados_filtrados <- metadados_filtrados %>%
        mutate(Taxa_Blastocisto = blastocisto / Total_Embrioes)
    }
    
    # 2. Criar o Dataframe para Plotagem
    df_plot_interacao <- data.frame(
      Abundancia_CLR = as.numeric(taxa_abund_clr[otu_id_malefico, ]),
      Taxa_Blastocisto = metadados_filtrados[colnames(taxa_abund_clr), "Taxa_Blastocisto"],
      Grupos = metadados_filtrados[colnames(taxa_abund_clr), "Grupos"]
    )
    
    # 3. Gerar o Gráfico
    g_scatter_interacao <- ggplot(df_plot_interacao, 
                                  aes(x = Abundancia_CLR, y = Taxa_Blastocisto, 
                                      color = Grupos)) +
      
      geom_point(size = 4, alpha = 0.8) +
      # A CHAVE: Geom_smooth separado por cor, provando a interação
      geom_smooth(method = "lm", se = TRUE, fullrange = TRUE) +
      
      labs(
        title = paste("Interação Condicional: Impacto de", nome_malefico, "no Desfecho"),
        subtitle = "Relação entre Abundância (CLR) e Taxa de Blastocisto, colorida por grupo",
        x = paste("Abundância do Gênero", nome_malefico, "(CLR)"),
        y = "Taxa de Blastocisto (Proporção)",
        color = "Grupo"
      ) +
      scale_color_manual(values = c("NATB" = "red", "ATB" = "blue")) +
      theme_bw(base_size = 14) +
      theme(plot.title = element_text(face = "bold"),
            legend.position = "bottom")
    
    print(g_scatter_interacao)
    
  } else {
    cat("AVISO: Nenhum gênero maléfico significativo (Log2FC > 0, padj < 0.05) encontrado para o Scatter Plot de Interação.\n")
  }
} else {
  cat("AVISO: Objetos necessários (taxa_abund_clr, metadados_filtrados ou res_final_genero) não encontrados.\n")
}

cat("\n--- Geração de Gráficos de Prova Concluída ---\n")