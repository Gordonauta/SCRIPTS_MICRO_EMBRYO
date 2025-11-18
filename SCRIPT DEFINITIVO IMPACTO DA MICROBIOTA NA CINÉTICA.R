# =================================================================
# SCRIPT DEFINITIVO: IMPACTO DA MICROBIOTA NA CINÉTICA (DIA 5)
# Foco: Identificar bactérias que causam "Bloqueio" (Pré-Mórula)
# =================================================================

# --- 0. CARREGAR BIBLIOTECAS ---
library(phyloseq)
library(tidyverse)
library(vegan)
library(compositions)
library(lme4)      # Para modelagem cinética robusta (GLMM)
library(car)       # Para Anova
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

# Dados de Embriões (D5)
desenvolvimento_df <- data.frame(
  NAME = c("JC_1", "JC_10", "JC_11", "JC_12", "JC_13", "JC_14", "JC_2", "JC_3", "JC_4", "JC_5", "JC_6", "JC_7", "JC_8", "JC_9"),
  premorula = c(0, 0, 4, 0, 0, 0, 0, 0, 7, 0, 14, 0, 0, 1),
  morula = c(12, 0, 0, 15, 1, 4, 15, 11, 5, 0, 0, 0, 0, 1),
  blastocisto = c(0, 11, 5, 1, 5, 1, 0, 2, 0, 0, 0, 9, 7, 0),
  row.names = 1
)

# Calcular Variáveis de Cinética D5
desenvolvimento_df <- desenvolvimento_df %>%
  mutate(
    Total_Embrioes = premorula + morula + blastocisto,
    # Variável Crítica: Taxa de Atraso (Bloqueio)
    Taxa_Atraso = premorula / Total_Embrioes, 
    # Variável de Sucesso (Transição Completa)
    Taxa_Sucesso = (morula + blastocisto) / Total_Embrioes,
    # Contagens para o modelo GLMM
    Count_Pre = premorula,
    Count_Sucesso = morula + blastocisto
  )

# Integrar e Filtrar
novos_metadados <- sample_data(ps)
novos_metadados <- cbind(novos_metadados, desenvolvimento_df[rownames(novos_metadados), ])
sample_data(ps) <- sample_data(novos_metadados)
ps <- prune_samples(sample_data(ps)$Total_Embrioes > 0, ps)

# Filtros e Transformação CLR
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
# PERGUNTA 1: O PERFIL DA MICROBIOTA AFETA A CINÉTICA?
# Teste: PERMANOVA usando a 'Taxa de Atraso' como explicadora
# =================================================================
cat("\n--- RESPONDENDO PERGUNTA 1: O perfil geral afeta o atraso? ---\n")

dist_bray <- phyloseq::distance(ps_filtrado_taxa, method = "bray")

# Testamos: A matriz de distâncias é explicada pela Taxa de Atraso?
permanova_atraso <- adonis2(dist_bray ~ Taxa_Atraso, data = meta_final)
print(permanova_atraso)

cat("\nINTERPRETAÇÃO:\n")
cat("Se Pr(>F) < 0.05, SIM, a composição da microbiota dita quem atrasa e quem avança.\n")

# Visualização PCoA (Colorido por Atraso)
ord_pcoa <- ordinate(ps_filtrado_taxa, "PCoA", "bray")
sample_data(ps_filtrado_taxa)$Taxa_Atraso <- meta_final$Taxa_Atraso

p_pcoa_atraso <- plot_ordination(ps_filtrado_taxa, ord_pcoa, color = "Taxa_Atraso") +
  geom_point(size = 6) +
  scale_color_viridis_c(option = "magma", direction = -1, name = "Taxa Atraso") + 
  labs(title = "PCoA: O Perfil da Microbiota e o Bloqueio Embrionário",
       subtitle = "Pontos Claros/Amarelos = Alto índice de Pré-Mórulas (Bloqueio)") +
  theme_bw(base_size = 14) +
  theme(panel.background = element_rect(fill = "gray90"))
print(p_pcoa_atraso)


# =================================================================
# PERGUNTA 2: QUAIS GÊNEROS TÊM PAPEL NO BLOQUEIO?
# Abordagem Combinada: Correlação (Triagem) + Regressão (Impacto)
# =================================================================
cat("\n--- RESPONDENDO PERGUNTA 2: Quem são os 'Bloqueadores'? ---\n")

# 1. TRIAGEM: Quem se correlaciona POSITIVAMENTE com Pré-Mórula?
# (Rho positivo aqui é ruim, pois significa associação com atraso)
cor_atraso <- apply(taxa_abund_clr, 1, function(x) cor(x, meta_final$Taxa_Atraso, method="spearman"))

# Selecionar Top 10 Candidatos a "Bloqueadores" (Rho > 0.4)
bloqueadores_candidatos <- names(sort(cor_atraso, decreasing = TRUE)[1:10])
# Selecionar Top 5 Candidatos a "Promotores" (Rho < -0.4)
promotores_candidatos <- names(sort(cor_atraso, decreasing = FALSE)[1:5])

# 2. CONFIRMAÇÃO CINÉTICA (GLMM Simplificado)
# Vamos ver se esses gêneros realmente alteram a contagem entre Estágios
cat("\n--- Validando Cinética com GLMM (Modelagem de Contagem) ---\n")

# Preparar dados longos para o modelo
meta_longo <- meta_final %>%
  rownames_to_column("Amostra") %>%
  pivot_longer(cols = c(premorula, morula, blastocisto), names_to = "Estagio", values_to = "Contagem") %>%
  mutate(Estagio = factor(Estagio, levels = c("premorula", "morula", "blastocisto")))

validacao_glmm <- data.frame()

for(otu in c(bloqueadores_candidatos, promotores_candidatos)) {
  # Dados de abundância do gênero
  abund <- data.frame(Amostra = colnames(taxa_abund_clr), Abundancia = as.numeric(taxa_abund_clr[otu, ]))
  df_mod <- left_join(meta_longo, abund, by="Amostra")
  
  # Modelo: A abundância interage com o estágio? 
  # (Ou seja, a abundância alta faz o embrião 'ficar' na pré-mórula?)
  try({
    mod <- glmer(Contagem ~ Abundancia * Estagio + (1|Amostra), data=df_mod, family=poisson, 
                 control=glmerControl(optimizer="bobyqa"))
    
    # Coeficiente de interação com Mórula (Transição chave)
    # Se for negativo, significa que ter a bactéria REDUZ a chance de virar mórula (Bloqueio)
    coefs <- summary(mod)$coefficients
    if("Abundancia:Estagiomorula" %in% rownames(coefs)) {
      est <- coefs["Abundancia:Estagiomorula", "Estimate"]
      pval <- coefs["Abundancia:Estagiomorula", "Pr(>|z|)"]
      
      validacao_glmm <- rbind(validacao_glmm, data.frame(
        OTU = otu, 
        Nome = get_tax_name(otu, ps_genero_filtrado),
        Papel = ifelse(otu %in% bloqueadores_candidatos, "Possível Bloqueador", "Possível Promotor"),
        Impacto_Transicao = est,
        P_Valor = pval
      ))
    }
  }, silent=TRUE)
}

# 3. RESULTADO FINAL
cat("\n--- 🎯 GÊNEROS COM IMPACTO NA CINÉTICA (D5) ---\n")
# Ordenar: Bloqueadores mais fortes terão Impacto de Transição mais NEGATIVO (impedem a mórula)
tabela_final <- validacao_glmm %>% 
  filter(P_Valor < 0.1) %>% # Filtro suave para visualizar tendências com N pequeno
  arrange(Impacto_Transicao)

print(tabela_final)

cat("\nNOTA DE INTERPRETAÇÃO:\n")
cat("- 'Impacto_Transicao' NEGATIVO: A bactéria dificulta a evolução para Mórula (Bloqueador).\n")
cat("- 'Impacto_Transicao' POSITIVO: A bactéria favorece a evolução para Mórula (Promotor).\n")


# =================================================================
# VISUALIZAÇÃO FINAL: CURVA CINÉTICA DOS BLOQUEADORES (CORRIGIDO)
# =================================================================

if(exists("tabela_final") && nrow(tabela_final) > 0) {
  
  cat("\n--- Gerando Gráfico de Curva Cinética (Correção de Dados) ---\n")
  
  # 1. Garantir que temos os metadados corretos com as taxas
  # Recalcula as taxas para garantir que não estão NULL
  meta_plot <- as(sample_data(ps_genero_filtrado), "data.frame") %>%
    mutate(
      Taxa_Premorula = premorula / Total_Embrioes,
      Taxa_Morula = morula / Total_Embrioes,
      Taxa_Blastocisto = blastocisto / Total_Embrioes
    )
  
  # 2. Pegar os Top 3 Bloqueadores (Mais Negativos) e Top 3 Promotores (Mais Positivos)
  top_bloqueadores <- head(tabela_final$OTU, 5) # Os 3 primeiros da lista ordenada
  top_promotores <- tail(tabela_final$OTU, 5)   # Os 3 últimos da lista ordenada
  otus_plot <- c(top_bloqueadores, top_promotores)
  
  # 3. Preparar dados para plotar curvas
  df_plot_curva <- data.frame()
  estagios_nomes <- c("1. Pré-Mórula", "2. Mórula", "3. Blastocisto")
  
  for(otu in otus_plot) {
    # Nome seguro
    nm <- get_tax_name(otu, ps_genero_filtrado)
    
    # Verificar se OTU existe na tabela de abundância
    if(otu %in% rownames(taxa_abund_clr)) {
      abund_vec <- as.numeric(taxa_abund_clr[otu, rownames(meta_plot)])
      
      # Spearman com cada estágio (usando use="complete.obs" para evitar erros de NA)
      rho1 <- cor(abund_vec, meta_plot$Taxa_Premorula, method="spearman", use="complete.obs")
      rho2 <- cor(abund_vec, meta_plot$Taxa_Morula, method="spearman", use="complete.obs")
      rho3 <- cor(abund_vec, meta_plot$Taxa_Blastocisto, method="spearman", use="complete.obs")
      
      # Define o tipo para colorir diferente no gráfico
      tipo <- ifelse(otu %in% top_bloqueadores, "Bloqueador (Vilão)", "Promotor (Mocinho)")
      
      df_plot_curva <- rbind(df_plot_curva, 
                             data.frame(Nome=nm, Tipo=tipo, Estagio=estagios_nomes[1], Rho=rho1),
                             data.frame(Nome=nm, Tipo=tipo, Estagio=estagios_nomes[2], Rho=rho2),
                             data.frame(Nome=nm, Tipo=tipo, Estagio=estagios_nomes[3], Rho=rho3))
    }
  }
  
  # 4. Plotar
  p_curva_final <- ggplot(df_plot_curva, aes(x=Estagio, y=Rho, group=Nome, color=Tipo)) +
    geom_line(linewidth=1.5, alpha=0.8) + 
    geom_point(size=4) +
    geom_hline(yintercept=0, linetype="dashed", color="gray50") +
    labs(title="Cinética: Bloqueadores vs. Promotores (D5)",
         subtitle="Bloqueadores (Vermelho): Alta correlação com atraso, cai no sucesso.\nPromotores (Azul): Baixa no atraso, sobe no sucesso.",
         y="Correlação (Spearman Rho)", x="Estágio de Desenvolvimento") +
    scale_color_manual(values=c("Bloqueador (Vilão)"="red", "Promotor (Mocinho)"="blue")) +
    theme_bw(base_size=14) +
    theme(legend.position="bottom")
  
  print(p_curva_final)
  
  # Salvar
  ggsave("GRAFICO_FINAL_CINETICA_BLOQUEADORES.png", p_curva_final, width=10, height=7, dpi=300)
  cat("\n--- Gráfico Final Salvo! ---\n")
  
} else {
  cat("\nAVISO: Tabela final vazia ou não encontrada. Rode a análise GLMM primeiro.\n")
}