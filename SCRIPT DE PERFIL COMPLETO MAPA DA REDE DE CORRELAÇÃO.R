# =================================================================
# SCRIPT DE PERFIL COMPLETO (D5): MAPA DA REDE DE CORRELAÇÃO
# =================================================================
# OBJETIVO:
# 1. Identificar os Top Gêneros associados ao "Atraso" e ao "Sucesso".
# 2. Calcular a rede de correlação (Hmisc) entre eles.
# 3. Desenhar um mapa visual (ggraph) dos perfis.
# =================================================================

# --- 0. CARREGAR BIBLIOTECAS ---
# install.packages("igraph")
# install.packages("ggraph")
library(phyloseq)
library(tidyverse)
library(compositions)
library(Hmisc)     # Para rcorr
library(ggplot2)
library(igraph)    # Para criar o objeto de rede
library(ggraph)    # Para plotar a rede com ggplot

cat("--- Bibliotecas Carregadas ---\an")

# --- 1. SETUP E DADOS (COM LÓGICA D5) ---
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


# Calcular Taxas D5
desenvolvimento_df <- desenvolvimento_df %>%
  mutate(
    Total_Embrioes = premorula + morula + blastocisto,
    Taxa_Atraso = premorula / Total_Embrioes, 
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


cat("--- Dados D5 Prontos ---\n")


# =================================================================
# PASSO 2: IDENTIFICAR OS GÊNEROS "CHAVE" (VILÕES E MOCINHOS)
# =================================================================

cat("\n--- Passo 2: Identificando Perfis de 'Vilões' e 'Mocinhos' ---\n")

# Calcular correlações de Spearman
cor_atraso <- apply(taxa_abund_clr, 1, function(x) cor(x, meta_final$Taxa_Atraso, method="spearman"))
cor_sucesso <- apply(taxa_abund_clr, 1, function(x) cor(x, meta_final$Taxa_Sucesso, method="spearman"))

# Top 10 Vilões (Maior Rho com Atraso)
viloes_ids <- names(sort(cor_atraso, decreasing = TRUE)[1:10])

# Top 10 Mocinhos (Maior Rho com Sucesso)
mocinhos_ids <- names(sort(cor_sucesso, decreasing = TRUE)[1:10])

# Lista de todos os gêneros-chave para o mapa
generos_chave_ids <- unique(c(viloes_ids, mocinhos_ids))

# Criar a Tabela de "Nós" (Nodes) para o gráfico
nodes_df <- data.frame(
  OTU_ID = generos_chave_ids,
  Nome = sapply(generos_chave_ids, get_tax_name, ps_obj = ps_genero_filtrado),
  Perfil = ifelse(generos_chave_ids %in% viloes_ids, "Vilão (Atraso)", "Mocinho (Sucesso)")
)

cat("--- Gêneros-Chave Identificados ---\n")


# =================================================================
# PASSO 3: CONSTRUIR A REDE (ARESTAS / EDGES)
# =================================================================

cat("\n--- Passo 3: Construindo a Rede de Correlação (Edges) ---\n")

# 1. Calcular a matriz de correlação COMPLETA (todos vs todos)
cor_matrix_geral <- Hmisc::rcorr(as.matrix(t(taxa_abund_clr)), type = "spearman")

# 2. "Aplanar" a matriz para uma lista de pares
flatten_cor_matrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    de = rownames(cormat)[row(cormat)[ut]],
    para = rownames(cormat)[col(cormat)[ut]],
    Rho = cormat[ut],
    P_valor = pmat[ut]
  )
}

rede_completa <- flatten_cor_matrix(cor_matrix_geral$r, cor_matrix_geral$P)

# 3. Filtrar a rede:

edges_df <- rede_completa %>%
  filter(de %in% generos_chave_ids & para %in% generos_chave_ids) %>%
  filter(P_valor < 0.05 & abs(Rho) > 0.5) %>%
  mutate(Tipo = ifelse(Rho > 0, "Positiva (Amizade)", "Negativa (Inimizade)"))

cat("--- Rede Filtrada Pronta (Edges) ---\n")


# =================================================================
# PASSO 4: DESENHAR O MAPA DE REDE (ggraph)
# =================================================================

cat("\n--- Passo 4: Gerando o Mapa de Rede Completo ---\n")

# 1. Criar o objeto de rede (igraph)
# (Usamos 'nodes_df' para os nós e 'edges_df' para as conexões)
rede_igraph <- graph_from_data_frame(d = edges_df, vertices = nodes_df, directed = FALSE)

# 2. Plotar com ggraph
p_rede_completa <- ggraph(rede_igraph, layout = 'kk') + # 'kk' é um layout "force-directed"
  
  # Desenhar as LINHAS (Edges)
  geom_edge_link(aes(color = Tipo, width = abs(Rho)), alpha = 0.8) +
  scale_edge_color_manual(values = c("Positiva (Amizade)" = "blue", "Negativa (Inimizade)" = "red")) +
  scale_edge_width(range = c(0.5, 3), name = "Força (Abs Rho)") +
  
  # Desenhar os PONTOS (Nodes)
  geom_node_point(aes(fill = Perfil), size = 12, shape = 21, color = "black") +
  scale_fill_manual(values = c("Vilão (Atraso)" = "#E41A1C", "Mocinho (Sucesso)" = "#377EB8")) +
  
  # Adicionar os NOMES
  geom_node_text(aes(label = Nome), repel = TRUE, size = 3, fontface = "bold") +
  
  theme_void() +
  labs(title = "Mapa do Perfil Completo: O Ecossistema do Atraso vs. Sucesso (D5)",
       subtitle = "Nós Vermelhos = Associados ao Atraso | Nós Azuis = Associados ao Sucesso\nLinhas = Conexões (P<0.05, |Rho|>0.5)")

print(p_rede_completa)

cat("\n--- Análise de Rede Completa Concluída ---\n")