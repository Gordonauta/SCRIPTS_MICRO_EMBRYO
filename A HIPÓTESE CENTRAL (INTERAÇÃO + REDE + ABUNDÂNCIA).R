# =================================================================
# SCRIPT 3: A HIPÓTESE CENTRAL (INTERAÇÃO + REDE + ABUNDÂNCIA)
# =================================================================
# OBJETIVO: Encontrar gêneros (Parceiros) que:
# 1. Se correlacionam com um Gênero-Chave (GC) no grupo NATB (Pior Desfecho)
# 2. E que são, eles mesmos, diferencialmente abundantes entre ATB e NATB
# =================================================================

cat("--- Iniciando Script da Hipótese Central ---\n")

# --- 0. CARREGAR BIBLIOTECAS ---
library(phyloseq)
library(tidyverse)
library(compositions) # Para CLR
library(DESeq2)       # Para Abundância Diferencial
library(Hmisc)        # Para rcorr (Rede de Correlação)

# --- 1. IMPORTAÇÃO E PREPARAÇÃO TOTAL DOS DADOS ---

# 1.1. Importação Básica
biom_file <- "final.opti_mcc.filter.pick.biom"
metadata_file <- "metadata.txt"
ps <- import_biom(biom_file)
sample_data_df <- read.table(metadata_file, header = FALSE, row.names = 1, sep = "", stringsAsFactors = FALSE)
colnames(sample_data_df) <- c("Grupos")
sample_data_df$Grupos <- as.factor(sample_data_df$Grupos)
sample_data(ps) <- sample_data(sample_data_df)

# 1.2. Integrar Dados de Embriões
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

# 1.3. Filtros e Agrupamento (GLOM)
amostras_para_remover <- rownames(sample_data(ps))[sample_data(ps)$Total_Embrioes == 0]
ps_filtrado_desenv <- prune_samples(!sample_names(ps) %in% amostras_para_remover, ps)
ps_filtrado_taxa <- prune_taxa(taxa_sums(ps_filtrado_desenv) >= 1, ps_filtrado_desenv)
ps_genero_bruto <- tax_glom(ps_filtrado_taxa, taxrank = "Rank6", NArm = FALSE)
ps_genero_filtrado <- prune_taxa(
  taxa_sums(otu_table(ps_genero_bruto) > 0) >= 3,
  ps_genero_bruto
)

cat(paste("--- Dados base prontos com", ntaxa(ps_genero_filtrado), "gêneros ---\n"))

# 1.4. [cite_start]Função para buscar nomes de Gênero [cite: 55-57]
get_tax_name <- function(otu_id, ps_obj, rank_level = "Rank6") {
  tax_info <- as(tax_table(ps_obj)[otu_id, ], "matrix")
  name <- tax_info[1, rank_level]
  if (is.na(name) || name == ""  || grepl("_group", name)) {
    name_rank5 <- as.character(tax_info[1, "Rank5"])
    if (!is.na(name_rank5) && name_rank5 != "") {
      name <- paste0(name_rank5, " (Família)")
    } else {
      name <- otu_id
    }
  }
  return(as.character(name))
}


# --- ETAPA 1: ENCONTRAR OS "GÊNEROS-CHAVE" (GCs) ---
# (Lógica do Teste de Interação LM) [cite_start][cite: 75-81]

cat("--- Etapa 1: Encontrando Gêneros-Chave (GCs)... ---\n")
otu_bruto_genero_pseudo <- as.data.frame(otu_table(ps_genero_filtrado)) + 1
clr_data_matrix <- t(clr(t(otu_bruto_genero_pseudo)))
taxa_abund_clr <- as.data.frame(clr_data_matrix)
metadados_filtrados <- as(sample_data(ps_genero_filtrado), "data.frame") %>%
  mutate(Taxa_Blastocisto = blastocisto / Total_Embrioes)

generos_para_testar <- rownames(taxa_abund_clr)
resultados_lm_interacao <- data.frame(OTU_ID = character(), P_valor_Interacao = numeric())

for (genero_id in generos_para_testar) {
  df_modelo <- data.frame(
    Taxa_Blastocisto = metadados_filtrados$Taxa_Blastocisto,
    Grupos = metadados_filtrados$Grupos,
    Abundancia_CLR = as.numeric(taxa_abund_clr[genero_id, ])
  )
  modelo_lm <- try(lm(Taxa_Blastocisto ~ Abundancia_CLR * Grupos, data = df_modelo), silent = TRUE)
  if (!inherits(modelo_lm, "try-error")) {
    sumario_modelo <- summary(modelo_lm)
    if ("Abundancia_CLR:GruposNATB" %in% rownames(sumario_modelo$coefficients)) {
      p_valor_interacao <- sumario_modelo$coefficients["Abundancia_CLR:GruposNATB", "Pr(>|t|)"]
      resultados_lm_interacao <- rbind(resultados_lm_interacao, data.frame(
        OTU_ID = genero_id, P_valor_Interacao = p_valor_interacao
      ))
    }
  }
}

df_interacao_sig <- resultados_lm_interacao %>%
  na.omit() %>%
  filter(P_valor_Interacao < 0.05)
lista_generos_chave_ids <- df_interacao_sig$OTU_ID

cat(paste("--- Etapa 1 Concluída:", length(lista_generos_chave_ids), "Gêneros-Chave (GCs) encontrados ---\n"))
print(df_interacao_sig)


# --- ETAPA 2: ENCONTRAR GÊNEROS DIFERENCIALMENTE ABUNDANTES (GDAs) ---
# (Lógica do Teste DESeq2) [cite_start][cite: 140-141]

cat("\n--- Etapa 2: Encontrando Gêneros Diferencialmente Abundantes (GDAs)... ---\n")
ps_genero_filtrado_plus1 <- ps_genero_filtrado
otu_table(ps_genero_filtrado_plus1) <- otu_table(ps_genero_filtrado_plus1) + 1
dds_obj <- phyloseq_to_deseq2(ps_genero_filtrado_plus1, design = ~ Grupos)
dds_result <- DESeq(dds_obj)
res_deseq <- results(dds_result, contrast = c("Grupos", "NATB", "ATB"))

# Usamos P-Bruto (pvalue) para uma rede mais ampla
df_gda <- as.data.frame(res_deseq) %>%
  rownames_to_column(var = "OTU_ID") %>%
  filter(pvalue < 0.05) %>%
  select(OTU_ID, log2FoldChange, pvalue, padj) %>%
  arrange(pvalue)

cat(paste("--- Etapa 2 Concluída:", nrow(df_gda), "GDAs encontrados (P-Bruto < 0.05) ---\n"))
print(df_gda)


# --- ETAPA 3: ENCONTRAR "PARCEIROS" DOS GCs (na rede NATB) ---
# (Lógica do Script 2 de Rede)

cat("\n--- Etapa 3: Encontrando 'Parceiros' dos GCs no grupo NATB (Pior Desfecho)... ---\n")

# 3.1. Preparar a função de rede (do Script 2)
run_and_flatten_network <- function(taxa_abund_df, ps_obj) {
  taxa_abund_t <- t(taxa_abund_df)
  cor_matrix <- Hmisc::rcorr(as.matrix(taxa_abund_t), type = "spearman")
  df_long <- data.frame()
  generos <- rownames(cor_matrix$r)
  for (i in 1:(length(generos) - 1)) {
    for (j in (i + 1):length(generos)) {
      otu_1 <- generos[i]
      otu_2 <- generos[j]
      df_temp <- data.frame(
        OTU_1 = otu_1, OTU_2 = otu_2,
        Rho = cor_matrix$r[otu_1, otu_2], P_valor_Bruto = cor_matrix$P[otu_1, otu_2]
      )
      df_long <- rbind(df_long, df_temp)
    }
  }
  df_long <- df_long %>% na.omit() %>%
    mutate(P_valor_BH = p.adjust(P_valor_Bruto, method = "BH"))
  return(df_long)
}

# 3.2. Rodar a rede APENAS no NATB (Pior Desfecho)
amostras_natb <- rownames(metadados_filtrados[metadados_filtrados$Grupos == "NATB", ])
taxa_abund_natb <- taxa_abund_clr[, amostras_natb]
rede_natb_df <- run_and_flatten_network(taxa_abund_natb, ps_genero_filtrado)

# 3.3. Filtrar apenas correlações fortes e significativas
rede_natb_sig <- rede_natb_df %>% 
  filter(P_valor_Bruto < 0.05)

# 3.4. Encontrar os Parceiros (o passo chave)
# Pares onde o GC é o Gênero 1
parceiros_1 <- rede_natb_sig %>%
  filter(OTU_1 %in% lista_generos_chave_ids) %>%
  select(Genero_Chave_ID = OTU_1, Parceiro_ID = OTU_2, Rho_NATB = Rho)

# Pares onde o GC é o Gênero 2
parceiros_2 <- rede_natb_sig %>%
  filter(OTU_2 %in% lista_generos_chave_ids) %>%
  select(Genero_Chave_ID = OTU_2, Parceiro_ID = OTU_1, Rho_NATB = Rho)

# Combinar e obter lista única de Parceiros
lista_parceiros_final <- rbind(parceiros_1, parceiros_2) %>%
  distinct(Parceiro_ID, .keep_all = TRUE) # Focamos nos Parceiros

cat(paste("--- Etapa 3 Concluída:", nrow(lista_parceiros_final), "'Parceiros' únicos de GCs encontrados na rede NATB ---\n"))
print(lista_parceiros_final)


# --- ETAPA 4: A DESCOBERTA - CRUZAR PARCEIROS E GDAs ---

cat("\n--- Etapa 4: Cruzando 'Parceiros' (da Rede NATB) com 'GDAs' (do DESeq2) ---\n")

# Hipótese: Os "Parceiros" (Etapa 3) são os "GDAs" (Etapa 2)?
# Vamos cruzar as duas listas.
candidatos_finais <- inner_join(
  lista_parceiros_final, # Lista de Parceiros (Etapa 3)
  df_gda,                  # Lista de GDAs (Etapa 2)
  by = c("Parceiro_ID" = "OTU_ID") # O "Join"
)

# --- 5. RESULTADO FINAL ---

if (nrow(candidatos_finais) > 0) {
  
  # Adicionar taxonomia para leitura
  candidatos_finais$Genero_Chave_Nome <- sapply(candidatos_finais$Genero_Chave_ID, get_tax_name, ps_obj = ps_genero_filtrado)
  candidatos_finais$Parceiro_Nome <- sapply(candidatos_finais$Parceiro_ID, get_tax_name, ps_obj = ps_genero_filtrado)
  
  cat("\n\n--- 🎯 SUCESSO! HIPÓTESE CONFIRMADA: ---")
  cat("\nOs seguintes gêneros são os 'elos perdidos':\n")
  cat("1. Eles são 'Parceiros' de um Gênero-Chave (no grupo NATB 'Doente').\n")
  cat("2. Eles são 'GDAs' (sua abundância mudou significativamente com o tratamento ATB).\n\n")
  
  print(candidatos_finais %>% 
          select(Parceiro_ID, Parceiro_Nome, 
                 Genero_Chave_Nome, Rho_NATB, 
                 log2FoldChange, pvalue, padj))
  
} else {
  cat("\n\n--- HIPÓTESE NÃO CONFIRMADA ---")
  cat("\nNenhum 'Parceiro' (da rede NATB) também foi um 'Gênero Diferencialmente Abundante' (do DESeq2).\n")
}

cat("\n--- Script da Hipótese Central Concluído ---\n")


# =================================================================
# SCRIPT DE VISUALIZAÇÃO: BOXPLOT DOS PARCEIROS (VST - CORRIGIDO)
# =================================================================
# OBJETIVO: Focar em 1 Gênero-Chave e plotar os boxplots de
# abundância VST (a "real", normalizada pelo DESeq2)
# dos seus "Parceiros".
# =================================================================

library(phyloseq)
library(tidyverse)
library(ggplot2)
library(DESeq2)

cat("--- Iniciando Script de Boxplot VST dos Parceiros (Corrigido) ---\n")

# --- 1. IMPORTAÇÃO E PREPARAÇÃO DOS DADOS ---
# (Nenhuma mudança aqui)

biom_file <- "final.opti_mcc.filter.pick.biom"
metadata_file <- "metadata.txt"
ps <- import_biom(biom_file)
sample_data_df <- read.table(metadata_file, header = FALSE, row.names = 1, sep = "", stringsAsFactors = FALSE)
colnames(sample_data_df) <- c("Grupos")
sample_data_df$Grupos <- as.factor(sample_data_df$Grupos)
sample_data(ps) <- sample_data(sample_data_df)

ps_filtrado_taxa <- prune_taxa(taxa_sums(ps) > 0, ps)
ps_genero_bruto <- tax_glom(ps_filtrado_taxa, taxrank = "Rank6", NArm = FALSE)
ps_genero_filtrado <- prune_taxa(
  taxa_sums(otu_table(ps_genero_bruto) > 0) >= 3,
  ps_genero_bruto
)

cat("--- Dados base prontos ---\n")

# --- 2. DEFINIR O ALVO (Baseado no seu output do Script 3) ---
# (Nenhuma mudança aqui)

genero_chave_alvo_nome <- "Streptococcus"
parceiros_alvo_nomes <- c(
  "UCG-010_genus", 
  "Phascolarctobacterium", 
  "Prevotellaceae_UCG-004"
)
cat(paste("--- Focando nos parceiros de:", genero_chave_alvo_nome, "---\n"))

# --- 3. RODAR DESEQ2 E OBTER DADOS VST ---

cat("--- Rodando DESeq2 e VST para obter abundância normalizada... ---\n")
ps_genero_filtrado_plus1 <- ps_genero_filtrado
otu_table(ps_genero_filtrado_plus1) <- otu_table(ps_genero_filtrado_plus1) + 1
dds_obj <- phyloseq_to_deseq2(ps_genero_filtrado_plus1, design = ~ Grupos)

# ***** MUDANÇA 1 AQUI: Corrigindo o 'note' do DESeq2 *****
# Adicionamos fitType='local' como o aviso sugeriu
dds_result <- DESeq(dds_obj, fitType='local')

vst_obj <- varianceStabilizingTransformation(dds_result, blind = TRUE)
vst_counts_matrix <- assay(vst_obj)

cat("--- Abundância VST (Normalizada) calculada ---\n")

# --- 4. PREPARAR DADOS PARA O GGPLOT (VST) ---

# 4.1. Função get_tax_name (Nenhuma mudança)
get_tax_name <- function(otu_id, ps_obj, rank_level = "Rank6") {
  tax_info <- as(tax_table(ps_obj)[otu_id, ], "matrix")
  name <- tax_info[1, rank_level]
  if (is.na(name) || name == "" || grepl("_unclassified", name) || grepl("_group", name)) {
    name_rank5 <- as.character(tax_info[1, "Rank5"])
    if (!is.na(name_rank5) && name_rank5 != "") {
      name <- paste0(name_rank5, " (Família)")
    } else {
      name <- otu_id
    }
  }
  return(as.character(name))
}

# 4.2. Tabela de Consulta (Nenhuma mudança)
taxa_nomes_df <- data.frame(
  OTU_ID = taxa_names(ps_genero_filtrado),
  Nome_Taxonomico = sapply(taxa_names(ps_genero_filtrado), get_tax_name, ps_obj = ps_genero_filtrado)
)

# 4.3. Encontrar IDs (Nenhuma mudança)
parceiros_alvo_ids <- taxa_nomes_df %>%
  filter(Nome_Taxonomico %in% parceiros_alvo_nomes) %>%
  pull(OTU_ID)

# 4.4. "Derreter" (Melt) a matriz VST

# ***** MUDANÇA 2 AQUI: Corrigindo o erro do 'left_join' *****
# Criamos o 'metadados_para_join' de forma mais robusta primeiro

metadados_para_join <- as(sample_data(ps_genero_filtrado), "data.frame") %>%
  rownames_to_column("Sample") %>%
  select(Sample, Grupos)

plot_data_full_vst <- vst_counts_matrix %>%
  t() %>% # Transpor (amostras = linhas)
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  # Adicionar os Grupos (ATB/NATB)
  left_join(metadados_para_join, by = "Sample") %>% # <-- Join corrigido
  
  # Converter para formato longo
  pivot_longer(
    cols = starts_with("Otu"), # Seleciona todas as colunas de OTU
    names_to = "OTU_ID",
    values_to = "Abundancia_VST"
  )

# 4.5. Filtrar (Nenhuma mudança)
plot_data_filtrado_vst <- plot_data_full_vst %>%
  filter(OTU_ID %in% parceiros_alvo_ids) %>%
  left_join(taxa_nomes_df, by = "OTU_ID")

cat("--- Dados VST prontos para plotagem ---\n")

# --- 5. GERAR OS BOX PLOTS (VST) ---
# (Nenhuma mudança aqui)

if (nrow(plot_data_filtrado_vst) > 0) {
  
  g_boxplot_parceiros_vst <- ggplot(plot_data_filtrado_vst, 
                                    aes(x = Grupos, y = Abundancia_VST, fill = Grupos)) +
    
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 2) +
    
    facet_wrap(~ Nome_Taxonomico, scales = "free_y") + 
    
    labs(
      title = paste("Abundância Normalizada (VST) dos 'Parceiros' de", genero_chave_alvo_nome),
      subtitle = "Esta é a abundância 'real' (corrigida) que o DESeq2 analisou",
      x = "Grupo de Tratamento (e Desfecho)",
      y = "Abundância Normalizada (VST)"
    ) +
    
    theme_bw() +
    theme(legend.position = "none", 
          plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "italic")) 
  
  print(g_boxplot_parceiros_vst)
  
} else {
  cat("!!! ERRO: Não foi possível encontrar os gêneros 'Parceiros' nos dados.\n")
  cat("Verifique os nomes no Passo 2.\n")
}

cat("\n--- Script de Boxplot VST dos Parceiros Concluído ---\n")