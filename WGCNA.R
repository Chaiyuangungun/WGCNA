library(WGCNA)
options(stringsAsFactors = FALSE)#不将字符串转换为因子
allowWGCNAThreads()#多线程
dataExpr=read.csv("catechin.FPKMs", sep='\t', header=T,check.names=F,row.names=1)
dim(dataExpr)
head(dataExpr)[,1:8]##查看数据
####数据筛选
m.mad <- apply(dataExpr,1,mad)#mad绝对中位差
dataExprVar <- dataExpr[which(m.mad >#quantile计算百分位数
                 max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]#判断绝对中位差是否大于0.25位或者大于0.01，判断绝对中位差0.25是否大于0.01
test= (t(dataExprVar))#转置行列
dataExpr = as.data.frame(test)#类型转换
gsg = goodSamplesGenes(dataExpr, verbose = 3)#检测缺失值
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:",
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:",
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)#查看列数与行数，及基因数和样本数
head(dataExpr)[,1:8]##查看数据
#####查看是否有离群样品
sampleTree = hclust(dist(dataExpr));#计算距离，dist()函数；hclust(系谱聚类)(参数：距离矩阵；聚类算法) average(类平均法)method："euclidean"表示欧氏距离, "maximum"表示最大距离, "manhattan"表示绝对值距离, "canberra"表示兰氏距离, "binary"或 "minkowski"表示闵可夫斯基距离，默认值为"euclidean"。
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()
#####导入代谢组数据
traitData = read.csv("catechin_content", sep = "\t", row.names=1,header=T,comment.char = "",check.names=F)
dim(traitData)
names(traitData)#查看代谢物名称
allTraits = traitData
fpkmSamples = rownames(dataExpr)
traitSamples =rownames(allTraits)
traitRows = match(fpkmSamples, traitSamples)
datTraits = allTraits[traitRows,]
rownames(datTraits) 
collectGarbage()
sampleTree2 = hclust(dist(dataExpr))
traitColors = numbers2colors(datTraits, signed = FALSE)
pdf(file="2.Sample_dendrogram_and_trait_heatmap.pdf",width=12,height=12)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()###样本聚类加对应热图
save(dataExpr, file = "fpkm_forAnalysis.RData")
save(datTraits, file="trait_forAnalysis.RData")
####软阈值筛选
powers = c(1:30)#幂值，在计算相关性是如果仅选择一个数值的话，不如0.8以上相关那么0.79就会出现问题。因此会对相关性进行幂次计算，这样整体相关性会下降，但相关性低的会更低
sft = pickSoftThreshold(dataExpr, powerVector = powers, verbose = 5)#powerVector
pdf(file="3_Scale_independence.pdf",width=14,height=9)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.9,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
#####基于 TOM 的差异的基因聚类
softPower =sft$powerEstimate
adjacency = adjacency(dataExpr, power = softPower)#获得临近矩阵:基因与基因相关性系数矩阵
TOM = TOMsimilarity(adjacency);#将邻近矩阵转为TOM矩阵：利用临近矩阵计算的新的临近矩阵，WGCNA任务简单的相关性不足以计算基因的共表达
# 计算基因之间的相异度
dissTOM = 1-TOM
#使用相异度来聚类为gene tree(聚类树)
geneTree = hclust(as.dist(dissTOM));
pdf(file="4.Gene_clustering_on_TOM-based_dissimilarity.pdf",width=12,height=9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()
#####识别基因集：基于基因加权相关性得到的TOM矩阵，进行层级聚类分析，并根据设定标准切分聚类结果，获得不同的基因模块
# 使用动态剪切树挖掘模块：
minModuleSize = 30#设置每个module最少的基因
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 4, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)#给每个模块赋予颜色
table(dynamicColors)
pdf(file="5.Dynamic_Tree_Cut.pdf",width=8,height=6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()
#####绘制模块之间相关性图
MEList = moduleEigengenes(dataExpr, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss))
pdf(file="6.Clustering_of_module_eigengenes.pdf",width=30,height=20)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.2######剪切高度根据图片做出修改
abline(h=MEDissThres, col = "red")
dev.off()
###聚类树状图
merge = mergeCloseModules(dataExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
pdf(file="7.merged_dynamic.pdf", width = 9, height = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
######模块特征关系
moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
save(MEs, TOM, dissTOM,  moduleLabels, moduleColors, geneTree, sft, file = "networkConstruction-stepByStep.RData")
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
moduleTraitCor = cor(MEs, datTraits, use = "p", method = 'pearson')
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
pdf(file="8.Module_trait_relationships.pdf",width=15,height=10)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(8, 8, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50), #colors = greenWhiteRed(50), 
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()
######
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(dataExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
traitNames=names(datTraits)
geneTraitSignificance = as.data.frame(cor(dataExpr, datTraits, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep="")

for (trait in traitNames){
  traitColumn=match(trait,traitNames)
  for (module in modNames){
    column = match(module, modNames)
    moduleGenes = moduleColors==module
    if (nrow(geneModuleMembership[moduleGenes,]) > 1){
      ####进行这部分计算必须每个模块内基因数量大于2，由于前面设置了最小数量是30，这里可以不做这个判断，但是grey有可能会出现1个gene,它会导致代码运行的时候中断，故设置
      #sizeGrWindow(7, 7)
      pdf(file=paste("9_", trait, "_", module,"_Module membership vs gene significance.pdf",sep=""),width=7,height=7)
      par(mfrow = c(1,1))
      verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                         abs(geneTraitSignificance[moduleGenes, traitColumn]),
                         xlab = paste("Module Membership in", module, "module"),
                         ylab = paste("Gene significance for ",trait),
                         main = paste("Module membership vs. gene significance\n"),
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
      dev.off()
    }
  }
}
########
names(dataExpr)
probes = names(dataExpr)
geneInfo0 = data.frame(probes= probes,
                       moduleColor = moduleColors)

for (Tra in 1:ncol(geneTraitSignificance))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneTraitSignificance[,Tra],
                         GSPvalue[, Tra])
  names(geneInfo0) = c(oldNames,names(geneTraitSignificance)[Tra],
                       names(GSPvalue)[Tra])
}

for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[,mod],
                         MMPvalue[, mod])
  names(geneInfo0) = c(oldNames,names(geneModuleMembership)[mod],
                       names(MMPvalue)[mod])
}
geneOrder =order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]
write.table(geneInfo, file = "10.GS_and_MM.xls",sep="\t",row.names=F)

nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
plotTOM = dissTOM^7
diag(plotTOM) = NA
pdf(file="11.Network_heatmap_plot_all_gene.pdf",width=9, height=9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
dev.off()

nSelect = 400
set.seed(10)
select = sample(nGenes, size = nSelect)
selectTOM = dissTOM[select, select]
selectTree = hclust(as.dist(selectTOM))
selectColors = moduleColors[select]
plotDiss = selectTOM^7
diag(plotDiss) = NA
pdf(file="12.Network_heatmap_plot_selected_genes.pdf",width=9, height=9)
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
dev.off()

pdf(file="13.Eigengene_dendrogram_and_Eigengene_adjacency_heatmap.pdf", width=5, height=7.5)
par(cex = 0.9)
plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle= 90)
dev.off()

for (mod in 1:nrow(table(moduleColors)))
{ 
  modules = names(table(moduleColors))[mod]
  probes = names(dataExpr)
  inModule = (moduleColors == modules)
  modProbes = probes[inModule]
  modGenes = modProbes
  modTOM = TOM[inModule, inModule]
  dimnames(modTOM) = list(modProbes, modProbes)
 cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("CytoscapeInput-edges-", modules , ".txt", sep=""),
                                 nodeFile = paste("CytoscapeInput-nodes-", modules, ".txt", sep=""),
                                 weighted = TRUE,
                                 threshold = 0.2,
                                 nodeNames = modProbes,
                                 altNodeNames = modGenes,
                                 nodeAttr = moduleColors[inModule])
}
