# =================================================================
# SCRIPT DE SALVAMENTO MESTRE (ATUALIZADO COM "MOVER")
# =================================================================
# (Rode este script POR ÚLTIMO, após executar todo o fluxo de análise)

# --- 0. CARREGAR BIBLIOTECAS ---
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(fs) # ✅ Carregar a biblioteca 'fs' para mover arquivos

cat("--- Iniciando Salvamento Mestre de Todos os Resultados ---\n")

# --- 1. CRIAR PASTAS DE SAÍDA ORGANIZADAS ---
dir_out <- "TODOS_RESULTADOS_FINAIS"
dir_base_tables <- file.path(dir_out, "1_Tabelas_Base_e_Cinetica")
dir_network_tables <- file.path(dir_out, "2_Tabelas_Rede_Complexa")
dir_main_plots <- file.path(dir_out, "3_Graficos_Principais") # 📂 Pasta de destino

dir.create(dir_out, showWarnings = FALSE)
dir.create(dir_base_tables, showWarnings = FALSE)
dir.create(dir_network_tables, showWarnings = FALSE)
dir.create(dir_main_plots, showWarnings = FALSE)

cat("--- Pastas de saída criadas em 'TODOS_RESULTADOS_FINAIS/' ---\n")

# =================================================================
# 2. SALVAR FASE 1: ANÁLISE BASE (Correlação e DESeq2)
# =================================================================
cat("Salvando resultados da Fase 1 (Base, DESeq2)...\n")

# --- Tabelas (DESeq2 - SCRIPT_2_Abundancia_Diferencial.R) ---
if (exists("res_final_genero")) {
  write.csv(res_final_genero, 
            file = file.path(dir_base_tables, "deseq2_abundancia_diferencial_completa.csv"), 
            row.names = FALSE)
}

# --- Tabelas (Spearman - analise-taxa-xlr.R) ---
if (exists("corr_results_blastocisto")) {
  write.csv(corr_results_blastocisto, 
            file = file.path(dir_base_tables, "spearman_corr_geral_blastocisto.csv"))
}

# --- Gráficos Principais (GRAFICOS_PROVAS_FINAIS.R) ---
if (exists("g_barplot_log2fc")) {
  ggsave(file.path(dir_main_plots, "PROVA_1_barplot_log2fc_atb_vs_natb.png"), 
         plot = g_barplot_log2fc, width = 10, height = 8, dpi = 300)
}

# =================================================================
# 3. SALVAR FASE 2: CINÉTICA E DISCREPÂNCIA
# =================================================================
cat("Salvando resultados da Fase 2 (Cinética GLMM e Discrepância)...\n")

# --- Tabelas (GLMM - SCRIPT_5_Analise_Estagios_GLMM.R) ---
if (exists("resultados_sig_interacao")) {
  write.csv(resultados_sig_interacao, 
            file = file.path(dir_base_tables, "cinetica_glmm_generos_interacao_significativa.csv"), 
            row.names = FALSE)
}

# --- Gráficos Principais (Gráficos de Discrepância) ---
if (exists("plot_discrepancia")) {
  ggsave(file.path(dir_main_plots, "PROVA_2_discrepancia_promotores_atb.png"), 
         plot = plot_discrepancia, width = 12, height = 10, dpi = 300)
}
if (exists("g_cinetica_rho_comparativa")) {
  ggsave(file.path(dir_main_plots, "PROVA_2B_discrepancia_cinetica_rho_comparativa.png"), 
         plot = g_cinetica_rho_comparativa, width = 12, height = 10, dpi = 300)
}

# =================================================================
# 4. SALVAR FASE 3: MECANISMOS DE REDE (Tabelas)
# =================================================================
cat("Salvando resultados da Fase 3 (Tabelas de Rede)...\n")

# --- Tabela (ANÁLISE DE REDE DE CORRELAÇÃO.R) ---
if (exists("rede_comparativa")) {
  write.csv(rede_comparativa, 
            file = file.path(dir_network_tables, "rede_geral_comparativa_atb_vs_natb.csv"), 
            row.names = FALSE)
}

# --- Tabela (SCRIPT_ANALISE_ANTAGONISTAS_COMPLETA.R) ---
if (exists("tabela_antagonistas")) {
  write.csv(tabela_antagonistas, 
            file = file.path(dir_network_tables, "rede_tabela_antagonistas_natb.csv"), 
            row.names = FALSE)
}

# --- Tabela (SCRIPT_ANALISE_AUXILIARES_COMPLETA.R) ---
if (exists("tabela_auxiliares")) {
  write.csv(tabela_auxiliares, 
            file = file.path(dir_network_tables, "rede_tabela_auxiliares_atb.csv"), 
            row.names = FALSE)
}

# --- Tabelas (SCRIPT_REDE_COMPLEXA_INDIVIDUAL.R) ---
if (exists("mediadores_encontrados")) {
  write.csv(mediadores_encontrados, 
            file = file.path(dir_network_tables, "rede_tabela_mediadores_natb.csv"), 
            row.names = FALSE)
}
if (exists("guardas_encontrados")) {
  write.csv(guardas_encontrados, 
            file = file.path(dir_network_tables, "rede_tabela_guardacostas_atb.csv"), 
            row.names = FALSE)
}

# =================================================================
# 5. (NOVO) MOVER GRÁFICOS SOLTOS (.PNG)
# =================================================================
cat("\n--- 5. Movendo gráficos .png soltos para a pasta '3_Graficos_Principais'... ---\n")

tryCatch({
  # Identificar todos os arquivos .png no diretório de trabalho principal
  arquivos_png_soltos <- dir_ls(path = ".", # "." significa o diretório atual
                                glob = "*.png") # glob = "*.png" seleciona apenas PNGs
  
  if (length(arquivos_png_soltos) > 0) {
    # Mover os arquivos para a pasta de destino
    file_move(arquivos_png_soltos, dir_main_plots)
    
    cat(paste("--- Sucesso:", length(arquivos_png_soltos), "gráficos .png movidos para '3_Graficos_Principais' ---\n"))
  } else {
    cat("--- Nenhum gráfico .png solto encontrado no diretório principal para mover. ---\n")
  }
}, error = function(e) {
  cat("--- ERRO AO MOVER GRÁFICOS: Verifique se a biblioteca 'fs' está instalada (install.packages('fs')) ---\n")
  print(e)
})


cat("\n--- SALVAMENTO MESTRE CONCLUÍDO --- \n")
cat("Verifique a pasta 'TODOS_RESULTADOS_FINAIS' no seu diretório de trabalho.\n")