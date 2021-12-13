
###########    Distância Molecular codominante  ################

## -- BIBLIOTECAS;

# Genind: Conversão de vcfR e análises;
library(adegenet)

# Utilidades para arquivos FSTAT e estatística genética;
library(hierfstat)  

# Para usar %>% e map
library(purrr)   

# Ler arquivos .vcf
library(vcfR)

# Para usar a função 'melt'
library(reshape2)

# Plots;
library(ggplot2)

# Criação de árvores filogenéticas;
library(ape)

# Comparação de clusters;
library(dendextend)

rm(list = ls())

# Salvar ou exibir os plots?
SAVE_FIGURES <- TRUE

# Análise baseada no seguinte tutorial:
# https://tomjenkins.netlify.app/2020/09/21/r-popgen-getting-started/#2

# Definir diretório base como o diretório onde situa-se este script;
rstudioapi::getActiveDocumentContext()$path %>% dirname %>% setwd

# Variantes localizadas dentro do gene NOD2;
vcf_file <- "1000G_chr22_19926475-19960318.vcf"

# Mapa de indivíduos e respectivas populações;
ped_file <- "20130606_g1k.ped"

# PARTE 1 - CARREGAR DADOS; ############################################

# Carregar o VCF obtido do 1000 genomes;
vcf <- read.vcfR(vcf_file, verbose = FALSE )

# Converter o VCF para genind;
GENIND <- vcfR2genind(vcf)

# Carregar o mapa indivíduo -> população;
# Arquivos .ped são basicamente .csv separados por tabs;
# Podemos carregar como data.frame sem depender de bibliotecas especializadas.
ped <- ped_file %>% read.csv(sep="\t")

# MAPEAR A POPULAÇÃO DE CADA INDIVÍDUO;
# Arquivos .vcf não possuem essa informação;
# Mas nossas análises exigem um objeto GENIND que contenha.
get_ind_pop <- function(ind) {
  ped[ped$Individual.ID == ind,]$Population[1]
}

# Vetor de população de cada indivíduo.
# (Adicionado no slot correspondente do objeto GENIND)
GENIND$pop <- GENIND$tab %>% rownames %>% map(~ get_ind_pop(.)) %>% unlist %>% as.factor

# Remover loci que não são polimórficos;
# (Existem polimorfismos que não são observados 
# na amostragem populacional do 1KG)

GENIND$loc.fac %>% length

# Extrai lista de loci que são polimórficos (isPoly verifica isso);
poly_loci = names(which(isPoly(GENIND) == TRUE))

# Filtra a lista atual
GENIND = GENIND[loc = poly_loci]

GENIND$loc.fac %>% length
# -- O arquivo genind está preparado.


# Ver os nomes das populações;
GENIND %>% popNames
# Ver os nomes das variantes;
GENIND %>% locNames

## -- PARTE 2: Análise de Heterozigosidade #########################################

# Calculate heterozygosity per SNP;
basic <- basic.stats(GENIND, diploid = TRUE)

#s, heterozigosidade esperada
#(He) e observada (Ho)
HO = apply(basic$Ho, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 3)

HE = apply(basic$Hs, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 3)


Het_df = data.frame(Site=names(HO), Ho = HO, He = HE) %>% melt(id.vars = "Site")

Het_df$value %>% sd
    
# -- Barplot de heterozigosidade;
custom_theme = theme(
  axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, face = "bold"),
  axis.text.y = element_text(size = 10),
  axis.title.y = element_text(size = 12),
  axis.title.x = element_blank(),
  axis.line.y = element_line(size = 0.5),
  legend.title = element_blank(),
  legend.text = element_text(size = 12),
  panel.grid = element_blank(),
  panel.background = element_blank(),
  plot.title = element_text(hjust = 0.5, size = 15, face="bold")
)

hetlab.o = expression(italic("H")[o])
hetlab.e = expression(italic("H")[e])

if (SAVE_FIGURES) {
  phylo_size <- 700
  png("Figures/het.png", width = phylo_size, height = phylo_size)
}  


ggplot(data = Het_df, aes(x = Site, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.6), colour = "black")+
  scale_y_continuous(expand = c(0,0), limits = c(0, 0.3)) +
  scale_fill_manual(values = c("royalblue", "#bdbdbd"), labels = c(hetlab.o, hetlab.e))+
  ylab("Heterozygosity") +
  ggtitle("NOD2 polymorphisms") +
  custom_theme

if (SAVE_FIGURES) {dev.off()}

## -- Correlação entre variação He vs Ho 
## -- contra o número de indivíduos em cada população.;
deltaH <- (HO - HE) %>% abs
n_ind <- GENIND %>% pop %>% table
ggplot(NULL, aes(x = deltaH, y=n_ind)) + geom_point() +
  ggtitle(paste("Populações: Número de Individuos vs deltaH. Corr:", cor(deltaH, n_ind) %>% round(2)))


## PARTE 3: Análise PCA ########################################################3

# Executar a PCA
# (Discriminant Analysis of Principal Components)
pca <- dapc.genind(GENIND, n.pca=20, n.da=5)


  if (SAVE_FIGURES) {
    phylo_size <- 700
    png("Figures/eig.png", width = phylo_size, height = phylo_size)
  }  

# Plotar contribuição relativa de cada componente;
percent = pca$eig / sum(pca$eig) * 100
barplot(percent, 
        ylab = "Genetic variance explained by eigenvectors (%)", 
        ylim = c(0, percent %>% max),
        names.arg = round(percent, 1),
        xlim = NULL
)

if (SAVE_FIGURES) dev.off()

if (SAVE_FIGURES) {
  phylo_size <- 700
  png("Figures/pca.png", width = phylo_size, height = phylo_size)
}  

# Plotar o scatter bidimensional da PCA;
scatter(pca, scree.da=F, bg="white", pch=20, cell=1, cstar=0, solid=1,
       cex=1.5, clab=0.7, legend = F) 

if (SAVE_FIGURES) dev.off()


# PARTE 4: DISTÂNCIAS ENTRE POPULAÇÕES, ÁRVORE E CLUSTERIZAÇÃO;

# -- Calcular distâncias entre as populações (método Cavalli-Sforza).
distances <- genet.dist(GENIND, diploid=TRUE, method="Dch")

## -- CLUSTERIZAÇÃO;
cluster1 <- hclust(distances, method = "average") %>% as.dendrogram
cluster2 <- hclust(distances, method = "complete") %>% as.dendrogram

if (SAVE_FIGURES) {
  phylo_size <- 700
  png("Figures/clusters.png", width = phylo_size, height = phylo_size)
}

dend_list <- dendlist(cluster1, cluster2)
dend_list %>%
  # Encontra o melhor layout.
  untangle(method = "step1side") %>%
  # Plot o comparativo de clusters;
  tanglegram()

if (SAVE_FIGURES) dev.off()


# -- ANÁLISES FUTURAS;
# 1. Baixar diversos .vcf: Regiões de genes diversos;
# 2. Executar toda essa análise para cada região, armazenando os objetos "cluster";
# 3. Comparar os clusters: 
#    - Para analisar marcadores de dois genes diferentes, 
#    - ou dois métodos de computar distância genômica, etc...