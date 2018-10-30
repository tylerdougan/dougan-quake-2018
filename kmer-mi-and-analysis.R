library("seqinr")
allviruses = read.fasta("allviruses.fsa")
uniqueAccessions <- readLines("accession-list-2017.txt")
uniqueAccessions[3266] = "NC_018570.2"
uniqueAccessions[4279] = "NC_025349.3"
uniqueAccessions[4695] = "NC_027348.2"
orderedviruses = list()
for(i in 1:length(uniqueAccessions)) {
	orderedviruses[[i]] = allviruses[[uniqueAccessions[i]]]
}
rm(allviruses)

fracCounts = list()
totalCount = count(orderedviruses[[1]], 4) - count(orderedviruses[[1]], 4)
for(i in 1:length(orderedviruses)) {
	fracCounts[[i]] = count(orderedviruses[[i]], 4)
	totalCount = totalCount + fracCounts[[i]]
	fracCounts[[i]] = fracCounts[[i]]/sum(fracCounts[[i]])
}
totalCount = totalCount/sum(totalCount)
vnum = length(orderedviruses)
mimat = matrix(data = rep(0, length.out = vnum^2), nrow = vnum, ncol = vnum)
distmat = mimat
dotprodmat = mimat
for(va in 1:vnum) {
	for(vb in 1:va) {
		counta = fracCounts[[va]]
		countb = fracCounts[[vb]]
		dotprodmat[va, vb] = counta %*% countb / norm(counta, "2") / norm(countb, "2")
		dotprodmat[vb, va] = dotprodmat[va, vb]
		for(overa in c(TRUE, FALSE)) {
			for(overb in c(TRUE, FALSE)) {
				pab = sum(totalCount[((counta > totalCount) == overa) & ((countb > totalCount) == overb)])
				if(is.finite(pab) && pab > 0) {
					pa = sum(totalCount[((counta > totalCount) == overa)])
					pb = sum(totalCount[((countb > totalCount) == overb)])
					mimat[va, vb] = mimat[va, vb] + pab * log(pab/(pa * pb), base = 2)
				}
			}
		}
		mimat[vb, va] = mimat[va, vb]
		distmat[va, vb] = mimat[va, va] + mimat[vb, vb] - 2 * mimat[va, vb]
		distmat[vb, va] = distmat[va, vb]
	}
}
plot(dotprodmat, distmat, pch = 20, cex = 0.3)
rm(orderedviruses)
library(smacof)
library(R.matlab)
inverseMatrix = readMat("blast-distance-2017-3.mat")$inverseMatrix
fixedistmat = distmat
fixedistmat[fixedistmat <= 0] = 0.01
inverseMatrix = readMat("blast-distance-2017-3.mat")$inverseMatrix
fixedistmat[fixedistmat < 0.1] = 0.1
inverseMatrix = inverseMatrix + min(inverseMatrix[inverseMatrix > 0])/fixedistmat
range(inverseMatrix)
distanceMatrix = 1/inverseMatrix
weightMatrix = inverseMatrix

rm(fracCounts)
diag(inverseMatrix) = 0
diag(distanceMatrix) = 0
pc50 = mds(distanceMatrix, ndim = 50, weightmat = weightMatrix)
library("Rtsne")
allTSNE = list()
for(i in c(1:10*5)) {
	allTSNE[[i]] = Rtsne(pc50$conf, is_distance = FALSE, verbose = TRUE, theta = 0.5, perplexity = i)
	errors = ((allTSNE[[i]])$itercosts)[length((allTSNE[[i]])$itercosts)]
}
plot(genesTSNE$Y[,1][abs(grepvec) < 0.1], genesTSNE$Y[,2][abs(grepvec) < 0.1], pch = 20, col = "gray50", cex = 0.4, main = "t-SNE plot of viruses, colored by host kingdom", xlab = "", ylab = "", xlim=range(genesTSNE$Y[,1]), ylim=range(genesTSNE$Y[,2]))
for(i in 1:length(hostKingdoms)) {
	points(genesTSNE$Y[,1][grepvec == 10^(i - 1)], genesTSNE$Y[,2][grepvec == 10^(i - 1)], pch = 20, col = colorsOfKingdoms[i], cex = 0.4)
}
setwd("/Users/tylerdougan/Box Sync/Research")
baltimore <- as.list(c(1:5))
# dsDNA Viruses
baltimore[[1]] <- readLines("1-dsdna.txt")
# ssDNA Viruses
baltimore[[2]] <- readLines("2-ssdna.txt")
# dsRNA Viruses
baltimore[[3]] <- readLines("3-dsrna.txt")
# ssRNA Viruses
baltimore[[4]] <- readLines("4-ssrna+.txt")
baltimore[[5]] <- readLines("5-ssrna-.txt")

greplist <- list(rep(0,length.out=length(uniqueAccessions)),rep(0,length.out=length(uniqueAccessions)),rep(0,length.out=length(uniqueAccessions)),rep(0,length.out=length(uniqueAccessions)),rep(0,length.out=length(uniqueAccessions)))
for(i in 1:5){
	for(j in 1:length(baltimore[[i]])){
		greplist[[i]] <- greplist[[i]] + grepl((baltimore[[i]])[j],uniqueAccessions)
	}
}
grepvec <- greplist[[1]] + 10*greplist[[2]] + 100*greplist[[3]] + 1000*greplist[[4]] + 10000*greplist[[5]]
for(i in 1:length(grepvec)) {
	if((grepvec[i] != 1) && (grepvec[i] != 10) && (grepvec[i] != 100) && (grepvec[i] != 1000) && (grepvec[i] != 10000)) {
		grepvec[i] = 0
	}
}
plot(genesTSNE$Y[,1][abs(grepvec) < 0.1], genesTSNE$Y[,2][abs(grepvec) < 0.1], pch = 20, col = "gray50", cex = 0.5, main = "t-SNE Plot of Viruses by Genes in Common", xlab = "", ylab = "", xlim=range(genesTSNE$Y[,1]), ylim=range(genesTSNE$Y[,2]))
points(genesTSNE$Y[,1][abs(grepvec-1) < 0.1], genesTSNE$Y[,2][abs(grepvec-1) < 0.1], pch = 20, col = "darkcyan", cex = 0.5)
points(genesTSNE$Y[,1][abs(grepvec-1000) < 0.1], genesTSNE$Y[,2][abs(grepvec-1000) < 0.1], pch = 20, col = "darkslateblue", cex = 0.5)
points(genesTSNE$Y[,1][abs(grepvec-10) < 0.1], genesTSNE$Y[,2][abs(grepvec-10) < 0.1], pch = 20, col = "darkorange1", cex = 0.5)
points(genesTSNE$Y[,1][abs(grepvec-100) < 0.1], genesTSNE$Y[,2][abs(grepvec-100) < 0.1], pch = 20, col = "deeppink", cex = 0.5)
points(genesTSNE$Y[,1][abs(grepvec-10000) < 0.1], genesTSNE$Y[,2][abs(grepvec-10000) < 0.1], pch = 20, col = "red", cex = 0.5)

plot(genesTSNE$Y[,1][abs(grepvec) < 0.1], genesTSNE$Y[,2][abs(grepvec) < 0.1], pch = 20, col = "gray50", cex = 0.5, main = "t-SNE Plot of Viruses, colored by Baltimore classification", xlab = "", ylab = "", xlim=range(genesTSNE$Y[,1]), ylim=range(genesTSNE$Y[,2]))
points(genesTSNE$Y[,1][abs(grepvec-1) < 0.1], genesTSNE$Y[,2][abs(grepvec-1) < 0.1], pch = 20, col = "darkcyan", cex = 0.5)
points(genesTSNE$Y[,1][abs(grepvec-1000) < 0.1], genesTSNE$Y[,2][abs(grepvec-1000) < 0.1], pch = 20, col = "darkslateblue", cex = 0.5)
points(genesTSNE$Y[,1][abs(grepvec-10) < 0.1], genesTSNE$Y[,2][abs(grepvec-10) < 0.1], pch = 20, col = "darkorange1", cex = 0.5)
points(genesTSNE$Y[,1][abs(grepvec-100) < 0.1], genesTSNE$Y[,2][abs(grepvec-100) < 0.1], pch = 20, col = "deeppink", cex = 0.5)
points(genesTSNE$Y[,1][abs(grepvec-10000) < 0.1], genesTSNE$Y[,2][abs(grepvec-10000) < 0.1], pch = 20, col = "red", cex = 0.5)
pckmer = cmdscale(distmat, k = 50)
pckmer = pckmer + rnorm(length(pckmer), sd = 1e-10)
kmreTSNE = Rtsne(pckmer, is_distance = FALSE, verbose = TRUE, theta = 0.5, perplexity = 40)
plot(kmreTSNE$Y[,1][abs(grepvec) < 0.1], kmreTSNE$Y[,2][abs(grepvec) < 0.1], pch = 20, col = "gray50", cex = 0.5, main = "t-SNE Plot of Viruses by K-Mers, colored by Baltimore classification", xlab = "", ylab = "", xlim=range(kmreTSNE$Y[,1]), ylim=range(kmreTSNE$Y[,2]))
points(kmreTSNE$Y[,1][abs(grepvec-1) < 0.1], kmreTSNE$Y[,2][abs(grepvec-1) < 0.1], pch = 20, col = "darkcyan", cex = 0.5)
points(kmreTSNE$Y[,1][abs(grepvec-1000) < 0.1], kmreTSNE$Y[,2][abs(grepvec-1000) < 0.1], pch = 20, col = "darkslateblue", cex = 0.5)
points(kmreTSNE$Y[,1][abs(grepvec-10) < 0.1], kmreTSNE$Y[,2][abs(grepvec-10) < 0.1], pch = 20, col = "darkorange1", cex = 0.5)
points(kmreTSNE$Y[,1][abs(grepvec-100) < 0.1], kmreTSNE$Y[,2][abs(grepvec-100) < 0.1], pch = 20, col = "deeppink", cex = 0.5)
points(kmreTSNE$Y[,1][abs(grepvec-10000) < 0.1], kmreTSNE$Y[,2][abs(grepvec-10000) < 0.1], pch = 20, col = "red", cex = 0.5)
