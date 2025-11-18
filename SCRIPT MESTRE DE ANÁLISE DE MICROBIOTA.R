# =================================================================
# SCRIPT MESTRE DE ANÁLISE DE MICROBIOTA (VERSÃO .R)
# =================================================================
# OBJETIVO:
# 1. Instalar e carregar todas as dependências.
# 2. Executar todos os 16 scripts de análise essenciais na ordem correta,
#    assumindo que TODOS terminam com a extensão .R.
# 3. Registrar todo o output (incluindo erros) em um arquivo de log.
# 4. Salvar todos os resultados finais (tabelas e gráficos).
# =================================================================
# --- 0. PREPARAÇÃO DO AMBIENTE ---

# Inicia o registro de log. Todo o output do console será salvo aqui.
sink("mestre_script.log", split = TRUE)

cat("--- INICIANDO SCRIPT MESTRE ---\n")
cat(paste("Data e Hora:", Sys.time(), "\n\n"))

# Lista de pacotes necessários para TODOS os scripts
pacotes_necessarios <- c(
  "phyloseq", "tidyverse", "dplyr", "vegan", "compositions", 
  "pheatmap", "ggplot2", "DESeq2", "ggrepel", "lme4", "car", 
  "Hmisc", "pls", "caret", "pROC", "e1071", "fs"
)

cat("--- 0.1 Instalando/Carregando Pacotes ---\n")
for (pacote in pacotes_necessarios) {
  if (!require(pacote, character.only = TRUE)) {
    cat(paste("Instalando pacote:", pacote, "\n"))
    install.packages(pacote)
  }
  library(pacote, character.only = TRUE)
}
cat("--- Todos os pacotes foram carregados com sucesso ---\n\n")


# =================================================================
# ATUALIZAÇÃO: Todos os arquivos .txt foram alterados para .R
# =================================================================
scripts_para_rodar <- c(
  
  # --- FASE 1: BASE DE DADOS (Setup, CLR, Diversidade) ---
  "analise-taxa-xlr.R",
  
  # --- FASE 2: DESCOBERTAS CENTRAIS (Abundância e Cinética) ---
  "SCRIPT_2_Abundancia_Diferencial.R",
  "SCRIPT_5_Analise_Estagios_GLMM.R",
  
  # --- FASE 3: VALIDAÇÃO (Machine Learning) ---
  "machine_learning.R",
  
  # --- FASE 4: MECANISMOS (Redes de Correlação) ---
  "ANÁLISE DE REDE DE CORRELAÇÃO.R",
  "SCRIPT_ANALISE_ANTAGONISTAS_COMPLETA.R",
  "SCRIPT_ANALISE_AUXILIARES_COMPLETA.R",
  "SCRIPT_REDE_COMPLEXA_INDIVIDUAL.R",
  
  # --- FASE 5: SÍNTESE (Cruzamento de Dados) ---
  "A HIPÓTESE CENTRAL (INTERAÇÃO + REDE + ABUNDÂNCIA).R",
  "ANÁLISE CRUZADA (INTERAÇÃO vs. ABUNDÂNCIA).R",
  "SCRIPT_CURVA_CINÉTICA_RHO_IMPACTO_GERAL.R", # Dependência para o próximo
  "CRUZAMENTO ML (VIP) vs. IMPACTO.R",
  
  # --- FASE 6: GRÁFICOS FINAIS (Visualização) ---
  "GRAFICOS_PROVAS_FINAIS.R",
  "CURVA_DE_CINÉTICA_COMPARATIVA_(ATB VS. NATB).R",
  "GRAFICO_PROMOTORES_ATB_DISCREPANCIA.R",
  
  # --- FASE 7: SALVAMENTO FINAL (Deve ser o último) ---
  "salvar resultados.R"
)

cat("--- 1. INICIANDO EXECUÇÃO SEQUENCIAL DOS SCRIPTS (.R) ---\n\n")

# --- 2. LOOP DE EXECUÇÃO ---

for (script in scripts_para_rodar) {
  cat("=================================================================\n")
  cat(paste("--- Executando:", script, "---\n"))
  cat("=================================================================\n\n")
  
  tryCatch({
    
    # Verifica se o arquivo existe
    if (!file.exists(script)) {
      stop(paste("Arquivo não encontrado:", script))
    }
    
    # Executa o script
    source(script, echo = TRUE, max.deparse.length = 1000)
    
    cat(paste("\n--- SUCESSO:", script, "concluído ---\n\n"))
    
  }, error = function(e) {
    
    # Captura e registra o erro
    cat("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
    cat(paste("--- ERRO FATAL AO EXECUTAR:", script, "---\n"))
    cat("MENSAGEM DE ERRO:\n")
    print(e)
    cat("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n")
    
  }) # Fim do tryCatch
  
} # Fim do loop

cat("=================================================================\n")
cat("--- EXECUÇÃO MESTRE CONCLUÍDA --- \n")
cat(paste("Data e Hora:", Sys.time(), "\n"))
cat("Verifique o arquivo 'mestre_script.log' para detalhes e erros.\n")
cat("Verifique a pasta 'TODOS_RESULTADOS_FINAIS' para os outputs.\n")
cat("=================================================================\n")

# Para o registro de log
sink()