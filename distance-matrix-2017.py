import numpy as np
import numpy.matlib
import scipy.io as sio
import os
from numpy import inf

# File dependencies: 'accession-list-2017.txt', 'genes-fasta-2017.fsa', 'blast-results-2017-1.txt',

genesToIndicies = dict()
listOfGenes = list() # list of all the genes in order
listOfViruses = list()
virusesToIndicies = dict() # dictionary from virus names to their corresponding index 0...5816
entropiesOfGenes = dict()
with open('accession-list-2017.txt') as uniqueAccessions:
	for j, line in enumerate(uniqueAccessions):
		virusesToIndicies[line.strip()] = j # give each virus its index 0...5816
		listOfViruses.append(line.strip())
with open('genes-fasta-2017.fsa') as geneNames:
	for j, line in enumerate(geneNames):
		if line[0] == '>':
			currentGeneName = line.split(',')[0].split(' ')[0][5:]
			listOfGenes.append(currentGeneName) # add each gene to the list of genes
			if currentGeneName[1] == 'C':
				genesToIndicies[currentGeneName] = virusesToIndicies[currentGeneName[0:11]]
			elif currentGeneName[1] == 'Z':
				genesToIndicies[currentGeneName] = virusesToIndicies[currentGeneName[0:13]]
			else:
				print('error reading file, line ' + str(j) + 'currentGeneName[1] is ' + currentGeneName[1])

if os.path.isfile('entropies-2017.txt'):
	with open('entropies-2017.txt') as entropiesOutput:
		for j, line in enumerate(entropiesOutput):
			separatedLine = line.split(',')
			entropiesOfGenes[separatedLine[0]] = float(separatedLine[1])
else:
	# read BLAST file to calculate entropies
	# Read the file once and make a map of the bitscores of genes against themselves
	with open('blast-results-2017-1.txt') as blastOutput:
		for j, line in enumerate(blastOutput):
			separatedLine = line.split(',')
			if separatedLine[0][4:] == separatedLine[1]:
				entropiesOfGenes[separatedLine[1]] = float(separatedLine[11])
	print('read file, found self-matches for ' + str(len(entropiesOfGenes)) + ' genes')
	# find entropies for genes without self-matches
	for i in range(0, len(listOfGenes)):
		if listOfGenes[i] not in entropiesOfGenes:
			highestBitScore = 20.01
			# print('getting entropy for ' + listOfGenes[i])
			with open('blast-results-2017-1.txt') as blastOutput:
				for j, line in enumerate(blastOutput):
					separatedLine = line.split(',')
					if (separatedLine[0][4:] == listOfGenes[i] or separatedLine[1] == listOfGenes[i]) and float(separatedLine[11]) > highestBitScore:
						highestBitScore = float(separatedLine[11])
			entropiesOfGenes[listOfGenes[i]] = highestBitScore + 100
	print('found all entropies.')
	with open('entropies-2017.txt', 'w') as entropiesOutput:
		for gene, entropy in entropiesOfGenes.iteritems():
			entropiesOutput.write(gene + ',' + str(entropy) + ',\n')

genePairScores = dict()
inverseDistanceMatrix = np.zeros((5817, 5817), dtype = float)
with open('blast-results-2017-1.txt') as blastOutput:
	for j, line in enumerate(blastOutput):
		separatedLine = line.split(',')
		# if the two species are different
		if separatedLine[0][4:] != separatedLine[1]:
			# We're looking at how X matches Y. Do we already have a score for how Y matches X?
			if (separatedLine[1] + ',' + separatedLine[0][4:]) in genePairScores:
				# Does Bxx + Byy - Bxy - Byx > 0? If it's less than zero, make it equal to zero.
				if entropiesOfGenes[separatedLine[0][4:]] + entropiesOfGenes[separatedLine[1]] \
				- float(separatedLine[11]) - genePairScores[(separatedLine[1] + ',' + separatedLine[0][4:])] < 0:
					# print('how did you get a number less than zero?')
					# print(separatedLine[0][4:] + ', ' + separatedLine[1])
					separatedLine[11] = str(entropiesOfGenes[separatedLine[0][4:]] + entropiesOfGenes[separatedLine[1]]
						- genePairScores[(separatedLine[1] + ',' + separatedLine[0][4:])])
				inverseDistanceMatrix[genesToIndicies[separatedLine[0][4:]], genesToIndicies[separatedLine[1]]] += \
				1. / (entropiesOfGenes[separatedLine[0][4:]] + entropiesOfGenes[separatedLine[1]]
					- float(separatedLine[11]) - genePairScores[(separatedLine[1] + ',' + separatedLine[0][4:])] + 1)
				inverseDistanceMatrix[genesToIndicies[separatedLine[1]], genesToIndicies[separatedLine[0][4:]]] += \
				1. / (entropiesOfGenes[separatedLine[0][4:]] + entropiesOfGenes[separatedLine[1]]
					- float(separatedLine[11]) - genePairScores[(separatedLine[1] + ',' + separatedLine[0][4:])] + 1)
				del genePairScores[(separatedLine[1] + ',' + separatedLine[0][4:])]
			# If we don't already have a score, then put this one into our running list
			else:
				if (separatedLine[0][4:] + ',' + separatedLine[1]) in genePairScores:
					print('you done goofed')
				genePairScores[(separatedLine[0][4:] + ',' + separatedLine[1])] = float(separatedLine[11])
for genePair, pairScore in genePairScores.iteritems():
	separatedPair = genePair.split(',')
	inverseDistanceMatrix[genesToIndicies[separatedPair[0]], genesToIndicies[separatedPair[1]]] += \
	1. / (entropiesOfGenes[separatedPair[0]] + entropiesOfGenes[separatedPair[1]] - genePairScores[genePair] + 1)
	inverseDistanceMatrix[genesToIndicies[separatedPair[1]], genesToIndicies[separatedPair[0]]] += \
	1. / (entropiesOfGenes[separatedPair[0]] + entropiesOfGenes[separatedPair[1]] - genePairScores[genePair] + 1)
del genePairScores

# Lowest inverse distance, which corresponds to greatest actual distance
minInverseDistance = np.amin(inverseDistanceMatrix[np.nonzero(inverseDistanceMatrix)])
print('Lowest inverse distance is ' + str(minInverseDistance) + ', which corresponds to a greatest actual distance of ' + str(1/minInverseDistance))
# Greatest inverse distance, which corresponds to lowest actual distance
maxInverseDistance = np.amax(inverseDistanceMatrix)
print('Greatest inverse distance is ' + str(maxInverseDistance) + ', which corresponds to a lowest actual distance of ' + str(1/maxInverseDistance))
distanceMatrix = 1. / (inverseDistanceMatrix + minInverseDistance)

#triangleMatrixOne = distanceMatrix
#triangleMatrixTwo = distanceMatrix
#print('beginning triangle inequality reducing...')
#inverseIntermediateTriangleMatrix = np.zeros((5817, 5817), dtype = float)
#for i in range(1,5817):
#	for j in range(0, i):
#		if inverseDistanceMatrix[i,j] == 0:
#			for k in range(0, 5817):
#				# if d(IK) and d(JK) are defined, and d(IK) + d(JK) is the smallest value we can find for d(IJ)
#				if inverseDistanceMatrix[i,k] > 0 and inverseDistanceMatrix[j,k] > 0 and 1. / (1. / inverseDistanceMatrix[i,k] + \
#				1. / inverseDistanceMatrix[j,k]) > inverseIntermediateTriangleMatrix[i,j]:
#					inverseIntermediateTriangleMatrix[i,j] = 1. / (1. / inverseDistanceMatrix[i,k] + 1. / inverseDistanceMatrix[j,k])
#					inverseIntermediateTriangleMatrix[j,i] = 1. / (1. / inverseDistanceMatrix[i,k] + 1. / inverseDistanceMatrix[j,k])
# we want to add the maximum value of the previous stage to all of the distances in the new stage
#intermediateTriangleMatrix = np.zeros((5817, 5817), dtype = float)
#intermediateTriangleMatrix += 1./(inverseIntermediateTriangleMatrix + minInverseDistance)
#intermediateTriangleMatrix += np.amax(triangleMatrixOne)
#intermediateTriangleMatrix[intermediateTriangleMatrix == inf] = 0
#intermediateTriangleMatrix[intermediateTriangleMatrix == -inf] = 0
# triangleMatrixOne adds the maximum previous distance to all of the new distances
#triangleMatrixOne = triangleMatrixOne + intermediateTriangleMatrix
#intermediateTriangleMatrix += np.amax(triangleMatrixTwo)
#intermediateTriangleMatrix[intermediateTriangleMatrix == np.amax(triangleMatrixTwo)] = 0
# triangleMatrixTwo adds twice the maximum previous distance to all of the new distances
#triangleMatrixTwo = triangleMatrixTwo + intermediateTriangleMatrix

#sio.savemat('blast-distance-2017-1.mat', {'distanceMatrix':distanceMatrix, 'triangleMatrixOne':triangleMatrixOne, 'triangleMatrixTwo':triangleMatrixTwo})

sio.savemat('blast-distance-2017-2.mat', {'distanceMatrix':distanceMatrix})
