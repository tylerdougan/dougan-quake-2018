# O(1)
start_time = Sys.time()
library("seqinr")
library(MASS)
library("Rtsne")

setwd("~/Documents/Quake/New 2017/Write-Up")
load ("the-end-workspace")
setwd("~/Documents/Quake/New 2017/")

# O(N) - about 4 minutes
allviruses = read.fasta("allviruses.fsa")
uniqueAccessions <- readLines("accession-list-2017.txt")
uniqueAccessions[3266] = "NC_018570.2"
uniqueAccessions[4279] = "NC_025349.3"
uniqueAccessions[4695] = "NC_027348.2"
oldviruses = list()
for(i in 1:length(uniqueAccessions)) {
	oldviruses[[i]] = allviruses[[uniqueAccessions[i]]]
}
rm(allviruses)
onum = length(oldviruses)
landmarkidx = 1:onum
t1closed = Sys.time()

# O(n)
newviruses = read.fasta("five-thousand.fsa", strip.desc = TRUE)
newNames = c()
newAccessions = c()
nnum = length(newviruses)
for(i in 1:nnum) {
	newNames[i] = substring(attr(newviruses[[i]], "Annot"), first = 10)
	newAccessions[i] = attr(newviruses[[i]], "name")
}
orderedviruses = append(oldviruses, newviruses)
rm(oldviruses)
rm(newviruses)
orderedNames = c(fullnames, newNames)
orderedAccessions = c(uniqueAccessions, newAccessions)

# O(N + n) - about 8 minutes
t2open = Sys.time()
fracCounts = list()
vnum = length(orderedviruses)
for(i in 1:vnum) {
	fracCounts[[i]] = count(orderedviruses[[i]], 4)
	fracCounts[[i]] = fracCounts[[i]]/sum(fracCounts[[i]])
}
isover = list()
ppositive = c()
for(i in 1:vnum) {
	isover[[i]] = fracCounts[[i]] > totalCount
	ppositive[i] = sum(isover[[i]] * totalCount)
}
t2closed = Sys.time()

# O((N+n)^2) or, at best, O(N*n + n^2) - 1.146519 hours
# should work for ~14k new points if optimized or ~11k if not
mimat2 = matrix(data = rep(0, length.out = vnum^2), nrow = vnum, ncol = vnum)
distmat2 = matrix(data = rep(0, length.out = vnum^2), nrow = vnum, ncol = vnum)
for(va in 1:vnum) {
	for(vb in 1:va) {
		abool = isover[[va]]
		bbool = isover[[vb]]
		for(overa in c(TRUE, FALSE)) {
			for(overb in c(TRUE, FALSE)) {
				pab = sum(totalCount[(abool == overa) & (bbool == overb)])
				if(is.finite(pab) && pab > 0) {
					mimat2[va, vb] = mimat2[va, vb] + pab * log(
						pab/(sum(totalCount[isover[[va]] == overa]) * 
							sum(totalCount[isover[[vb]] == overb])), base = 2)
				}
			}
		}
		mimat2[vb, va] = mimat2[va, vb]
	}
}

for(va in 1:vnum) {
	for(vb in 1:va) {
		distmat2[va, vb] = mimat2[va, va] + mimat2[vb, vb] - 2 * mimat2[va, vb]
		distmat2[vb, va] = distmat2[va, vb]
	}
}
rm(orderedviruses)

newdist = distmat2
diag(newdist) = 0

perplex = 30
contenders2 = c()
smallcluster2 = c()
agreement2 = c()
winnerCount = c()
for(i in 1:nnum) {
	knn = order(distmat2[i+onum,1:onum])[1:perplex]
	weights = 1/sqrt(clusterSizes[groups[knn]])
	options = unique(groups[knn])
	contenders2[i] = length(options)
	countingWeights = c()
	for(j in 1:contenders2[i]) {
		countingWeights[j] = sum(weights[groups[knn] == options[j]])
	}
	agreement2[i] = max(countingWeights)/sum(countingWeights)
	smallcluster2[i] = options[which.max(countingWeights)]
	winnerCount[i] = length(which(groups[knn] == smallcluster2[i]))
}