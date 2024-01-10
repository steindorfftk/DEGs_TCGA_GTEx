library(UCSCXenaTools)
library(data.table)
library(R.utils)
library(dplyr)

#01 - Gerar objeto para rastrear datasets no Xena
data(XenaData)
write.csv(XenaData,"00_tblXenaHubInfo.csv")

#02.1 - Baixar expected counts
GeneExpectedCnt_toil = XenaGenerate(subset = XenaHostNames == "toilHub") %>%
  XenaFilter(filterCohorts = "TCGA TARGET GTEx") %>%
  XenaFilter(filterDatasets = "TcgaTargetGtex_gene_expected_count");
XenaQuery(GeneExpectedCnt_toil) %>%
  XenaDownload (destdir = "/home/local.hcpa.ufrgs.br/tkruger/Desktop/projetos/Glioblastoma/DEGglobal_TCGA_GTEx")

#02.2 - Baixar dados clínicos
paraCohort = 'TCGA Colon Cancer'; #PRECISA MUDAR POR CÂNCER
paraDatasets = 'TCGA.COAD.samplemap/COAD_clinicalMatrix'; #PRECISA MUDAR POR CÂNCER

Clin_TCGA = XenaGenerate(subset = XenaHostNames == 'tcgaHub') %>%
  XenaFilter(filterCohorts = paraCohort) %>%
  XenaFilter(filterDatasets = paraDatasets);
XenaQuery(Clin_TCGA) %>%
  XenaDownload(destdir = "/home/local.hcpa.ufrgs.br/tkruger/Desktop/projetos/Glioblastoma/DEGglobal_TCGA_GTEx")

#02.3 - Baixar dados de sobrevida
Surv_TCGA = XenaGenerate(subset = XenaHostNames == 'toilHub') %>%
  XenaFilter(filterCohorts = 'TCGA TARGET GTEx') %>%
  XenaFilter(filterDatasets = 'TCGA_survival_data');
XenaQuery(Surv_TCGA) %>%
  XenaDownload(destdir = "/home/local.hcpa.ufrgs.br/tkruger/Desktop/projetos/Glioblastoma/DEGglobal_TCGA_GTEx")

#02.4 - Baixar dados de fenótipo
Pheno_GTEx = XenaGenerate(subset = XenaHostNames == 'toilHub') %>%
  XenaFilter(filterCohorts = 'TCGA TARGET GTEx') %>%
  XenaFilter(filterDatasets = 'TcgaTargetGTEX_phenotype');
XenaQuery(Pheno_GTEx) %>%
  XenaDownload(destdir = "/home/local.hcpa.ufrgs.br/tkruger/Desktop/projetos/Glioblastoma/DEGglobal_TCGA_GTEx")

#03 - Subamostrar dados de expressão para incluir apenas amostras desejadas
#03.1 - Recuperar IDs para amostras normais do GTEx do tecido desejado
filterGTEx01 = fread('TcgaTargetGTEX_phenotype.txt.gz');
names(filterGTEx01) = gsub('\\_','',names(filterGTEx01));

paraStudy = 'GTEX';
paraPrimarySiteGTEx = 'Colon'; #PRECISA MUDAR POR CÂNCER
paraPrimaryTissueGTEx = '^Colon'; #PRECISA MUDAR POR CÂNCER

filterGTEx02 = subset(filterGTEx01,
                      study==paraStudy &
                        primarysite==paraPrimarySiteGTEx &
                        grepl(paraPrimaryTissueGTEx, filterGTEx01$'primary disease or tissue'))

#03.2 - Recuperar IDs para amostras de tumor primário do TCGA do tecido desejado
filterTCGA01 = fread('COAD_clinicalMatrix');
names(filterTCGA01) = gsub('\\_','',names(filterTCGA01));

paraSampleType = 'Primary Tumor';
paraPrimarySiteTCGA = 'Colon'; #PRECISA MUDAR POR CÂNCER
paraHistologicalType = 'Colon Adenocarcinoma'; #PRECISA MUDAR POR CÂNCER

filterTCGA02 = subset(filterTCGA01,
                      sampletype == paraSampleType &
                        primarysite == paraPrimarySiteTCGA &
                        grepl(paraHistologicalType, filterTCGA01$histologicaltype))


#03.3 - Juntar listas do TCGA e do GTEx e coletar os perfis de expressão do arquivo .gz por ID
filterExpr = c(filterGTEx02$sample, filterTCGA02$sampleID, 'sample');
ExprSubsetBySamp = fread('TcgaTargetGtex_gene_expected_count_complete.gz',
                         select = filterExpr)

#04 - Subamostrar dados de expressão para incluir apenas genes codificantes de proteínas
probemap = fread('probeMap_gencode.v23.annotation.gene.probemap',select=c(1,2));
exprALL = merge(probemap,ExprSubsetBySamp,by.x = 'id',by.y='sample');
genesPC = fread('zz_gene.protein.coding.csv');
exprPC = subset(exprALL, gene %in% genesPC$Gene_Symbol);

#04.1 - Remover gene symbols duplicados
exprFinalA = exprPC[!duplicated(exprPC$gene) |
                     duplicated(exprPC$gene, fromLast = TRUE), ]

exprFinal = exprFinalA[!duplicated(exprFinalA$gene) |
                     duplicated(exprFinalA$gene, fromLast = TRUE), ]

#05 - Salvar dados do perfil de expressão para análises posteriores
write.csv(exprFinal, '00_ExpectedCnt.csv')

#06 - Preparar dados clínicos para incluir apenas variáveis desejadas
names(filterTCGA02);

#06.1 - Manter variável 'Lymphatic Invasion'
varClinKeep = c('sampleID','lymphaticinvasion');
clinDF01 = as.data.frame(do.call(cbind,filterTCGA02));
clinFinal = clinDF01[varClinKeep];

colSums(clinFinal=='');
colSums(is.na(clinFinal));
NA -> clinFinal[clinFinal == ''];
colSums(is.na(clinFinal));

table(clinFinal$lymphaticinvasion);

clinFinal$lymphaticinvasion=if_else(clinFinal$lymphaticinvasion=="YES",
                                    1, 0, missing = NULL);

table(clinFinal$lymphaticinvasion)

#07 - Salvar dados clínicos para análises posteriores
write.csv(clinFinal,'00_ClinTraits.csv')

#08 - Carregar pacotes para análise de expressão diferencial
library(dplyr)
library(limma)
library(edgeR)

#09 - Reverter transformação de counts esperados em log
exprFinal = read.csv('00_ExpectedCnt.csv');
exprBT = exprFinal[,-c(1:3)];
rownames(exprBT) = exprFinal$gene;
exprBT = round(((2^exprBT)-1),0);
write.csv(exprBT,'01_ExpectedCntBT.csv')

#10 - Converter dados de count em objeto DGEList
expLIMMA = exprBT[,-which(names(exprBT) %in%
                            c('GTEX.WFON.1426.SM.5CHT1',
                              'GTEX.SUCS.1026.SM.5CHTC',
                              'GTEX.WFG8.1726.SM.5CHTE',
                              'GTEX.WFG8.1526.SM.5CHSI'))] #REMOVENDO 4 AMOSTRAS DETECTADAS COMO ANOMALIAS PELO LIMMA. TALVEZ PRECISE FAZER PARA OUTROS CÂNCERES

x = DGEList(expLIMMA)

#11 - Agrupar amostras por condição
snames = colnames(x);
group = substr(snames, 1, 4);
x$samples$group = group

#12 - Computar counts por milhão
CPM = cpm(x);
lcpm = cpm(x,log = TRUE);
L = mean(x$samples$lib.size) * 1e-6;
L;

M = median(x$samples$lib.size) * 1e-6;
M;

table(x$samples$group)

#13 - Remover genes com expressão baixa
keep.exprs = filterByExpr(x, group = group);
x = x[keep.exprs, , keep.lib.sizes = FALSE];
dim(x)

#13.1 - Gerar plot de log-CPM para QC
lcpm.cutoff = log2(10/M + 2/L);
nsamples = ncol(x);
par(mfrow=c(1,2));
plot(density(lcpm[,1]),lwd=2,ylim=c(0,0.26),las=2,main='',xlab='');
title(main="A. Expected Count",xlab = 'Log-CPM');
abline(v=lcpm.cutoff,lty=3);
for (i in 2:nsamples){
  den = density(lcpm[,i])
  lines(den$x, den$y, lwd=2)
}
lcpm = cpm(x,log=TRUE);
plot(density(lcpm[,1]),lwd=2,ylim=c(0,0.26),las=2,main='',xlab='');
title(main='B. Filtered Expected Count',xlab='Log-CPM');
abline(v=lcpm.cutoff,lty=3);
for (i in 2:nsamples){
  den = density(lcpm[,i])
  lines(den$x, den$y, lwd=2)
}

#14 - Computar fatores de escala para converter bibliotecas observadas em bibliotecas efetivas
x = calcNormFactors(x,method='upperquartile')
head(x$samples$norm.factors)

#15 - Gerar matriz de design
design = model.matrix(~0 + group)

#16 - Configurar contraste para comparação
colnames(design) = gsub('group','',colnames(design));
contr.matrix = makeContrasts(TCGAvsGTEX = TCGA - GTEX,
                             levels = colnames(design))

#17 - Transformar dados RNA-seq para modelo linear
v = voom(x,design,plot=TRUE)

#18 - Ajustar um modelo linear usando mínimos quadrados ponderados para cada gene
vfit = lmFit(v,design)
vfit = contrasts.fit(vfit,contrasts=contr.matrix)

#19 - Realizar suavização Bayes empírica de erros padrão
efit = eBayes(vfit);
plotSA(efit,main='Final model: Mean-variance trend')

#20 - Examinar o número de DEGs
summary(decideTests(efit))

#20.1 - Como o número de DEGs é muito grande, treat pode ser usado para computar o empirical Bayes moderated-t p-values relativo a um log-FC de 0.58
tfit = treat(vfit,lfc=0.58);
summary(decideTests(tfit))

#21 - Salvar os DEGs
DEGsTreat = topTreat(efit, n = Inf);
write.csv(DEGsTreat, 'COAD_DEGs.csv') #MUDAR POR CÂNCER




























