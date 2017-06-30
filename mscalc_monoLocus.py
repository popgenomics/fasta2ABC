#!/usr/bin/pypy
import os
import sys
import sys

datasetName = sys.argv[1]

def cr_sqrt(x):
	# returns the square root of a variable x
	if x == 0.0:
		return 0.0
	else:
		M = 1.0
		xN = x 
		while xN >= 2.0:
			xN = 0.25*xN
			M = 2.0*M
		while xN<0.5:
			xN = 4.0*xN
			M = 0.5*M
		A = xN
		B = 1.0-xN
		while 1==1:
			A = A*(1.0+0.5*B)
			B = 0.25*(3.0+B)*B*B
			if B<1.0E-15:
				return A*M

def cr_mean(x):
	# returns the mean of a list
	nElement = len(x)
	if nElement == 0:
		return(0)
	else:
		return(sum(x)/(1.0 * nElement))

def cr_std(x, exp_X):
	# returns the standard variation of a list
	nElement = len(x)
	if nElement == 0:
		return(0)
	else:
		A = sum( [ i**2 for i in x ] )
		A = A/(1.0 * nElement)
		return(cr_sqrt(A-exp_X**2))


def cr_pearsonR(x, y):
	# computes correlation between arrays x and y
	sumXi = 0.0
	sumYi = 0.0
	sumSquareX = 0.0
	sumSquareY = 0.0
	sumXiYi = 0.0
	
	nX = len(x)
	
	for i in range(nX):
		sumXi += x[i];
		sumYi += y[i];
		sumSquareX += x[i]*x[i];
		sumSquareY += y[i]*y[i];
		sumXiYi += x[i]*y[i];
	
	numerator = nX*sumXiYi - sumXi * sumYi
	denom1 = cr_sqrt(nX*sumSquareX - sumXi*sumXi)
	denom2 = cr_sqrt(nX*sumSquareY - sumYi*sumYi)
	return(numerator/(denom1*denom2))


def compFreq(sequences, segsites):
	# returns derived allele frequency for a number of 'segsites' positions
	nDerAll = []
	nInd = len(sequences)
	nPair = nInd*(nInd-1)/2.0
	for i in range(segsites):
		nDerAll.append(0)
		for j in sequences:
			if j[i] == "1":
				nDerAll[i] += 1.0
	pi = [ i * (nInd-i) / nPair for i in nDerAll ]
	freq = [ i/(1.0 * nInd) for i in nDerAll ]
	res = {}
	res['nDer'] = nDerAll
	res['pi_SNPs'] = pi
	res['pi'] = sum(pi)
	res['freq'] = freq
#			if i == 0:		# uncomment to print sequence j #1
#				print(j)	# uncomment to print sequence j #2
	return(res)

def piTot(nDerA, nDerB, nSamA, nSamB, segsites):
	# returns pi tot for the pooled populations from two vectors of allele count per locus
	piT = []
	nTot = nSamA + nSamB
	for i in range(segsites):
		tmp = nDerA[i] + nDerB[i]
		tmp = (nTot - tmp) * tmp / (nTot*(nTot-1)/2.0) # "n ancesral" x "n derived" / C(2,k)
		piT.append(tmp)
	return(piT)


def Fst(piA, piB, piT, segsites):
	# returns Fst
	if piT==0:
		res = 0.0
	else:
		res = 1.0 - (piA+piB)/(2.0*piT)
	return(res)


def sites(freqA, freqB, segsites):
	sxA, sxB, ss, sf = 0, 0, 0, 0
	successive_ss = [0]
	previous_site = ""
	for i in range(segsites):
		if freqA[i] == 0:
			if freqB[i] == 1:
				sf += 1
				if previous_site == "ss" or previous_site == "":
					successive_ss.append(0)
				previous_site = "sf"
			if freqB[i] < 1:
				sxB += 1
				if previous_site == "ss" or previous_site == "":
					successive_ss.append(0)
				previous_site = "sxB"
			continue
		if freqA[i] == 1:
			if freqB[i] == 0:
				sf += 1
				if previous_site == "ss" or previous_site == "":
					successive_ss.append(0)
				previous_site = "sf"
			if freqB[i] > 0:
				sxB += 1
				if previous_site == "ss" or previous_site == "":
					successive_ss.append(0)
				previous_site = "sxB"
			continue
		else:
			if freqB[i] == 0 or freqB[i] == 1:
				sxA += 1
				if previous_site == "ss" or previous_site == "":
					successive_ss.append(0)
				previous_site = "sxA"
			else:
				ss += 1
				if previous_site == "ss" or previous_site == "":
					successive_ss[len(successive_ss)-1] += 1
				previous_site = "ss"
			continue
	res = {'sxA':sxA, 'sxB':sxB, 'sf':sf, 'ss':ss, 'successive_ss':max(successive_ss)}
	return(res)


def tajD(pi, theta, n, S, a1, a2):
	# returns Tajima D
	# pi = pi for a locus
	# theta = theta for a locus
	# n = sample size
	# S = number of polymorphic positions
	# a1 = sum(1/i)
	# a2 = sum(1/i**2)
	b1 = (n + 1.0) / (3*(n-1))
	b2 = 2.0 * (n**2 + n + 3) / (9*n*(n-1.0))
	c1 = b1 - 1.0/a1
	c2 = b2 - (n + 2.0)/(a1*n) + a2/(a1**2)
	e1 = c1/a1
	e2 = c2/(a1**2 + a2)
	denom = cr_sqrt(e1*S + e2*S*(S-1))
	if denom != 0:
		D = (pi - theta) / denom
	else:
		D = 0.0
	return(D)

def compDiv(spA, spB, segsites):
	div = [] # vector of divergence between A and B
	nPair = 0	
	for i in spA:
		for j in spB:
			nPair += 1
			div.append(0)
			for k in range(segsites):
				if i[k]!=j[k]:
					div[nPair - 1] += 1
	res = {}
	res['divAB'] = cr_mean(div)
	res['minDivAB'] = min(div)
	res['maxDivAB'] = max(div)
	if res['divAB'] > 0:
		res['Gmin'] = res['minDivAB']/res['divAB']
		res['Gmax'] = res['maxDivAB']/res['divAB']
	else:
		res['Gmin'] = 0.0
		res['Gmax'] = 0.0
	return(res)	


def nSs_nSf(ss, sf):
	# returns the number of loci with ss and sf, with ss but no sf, etc ...
	nLoci = len(ss)
	ss_sf, noSs_sf, ss_noSf, noSs_noSf = 0, 0, 0, 0
	for i in range(nLoci):
		if ss[i] == 0:
			if sf[i] == 0:
				noSs_noSf += 1
			if sf[i] > 0:
				noSs_sf += 1
		else:
			if sf[i] == 0:
				ss_noSf += 1
			if sf[i] > 0:
				ss_sf += 1
	res = {'ss_sf': ss_sf, 'ss_noSf': ss_noSf, 'noSs_sf': noSs_sf, 'noSs_noSf': noSs_noSf}
	return(res)


# spinput.txt
if os.path.isfile("spinput.txt") == False:
	sys.exit("\n\tERROR: spinput.txt was not found")

infile = open("spinput.txt", "r")

tmp = infile.readline() # first empty line
nLoci = int(infile.readline())

nSamA, nSamB, L = [], [], []
for i in range(nLoci):
	tmp = infile.readline().strip()
	nSamA.append(int(tmp))
	tmp = infile.readline().strip()
	nSamB.append(int(tmp))
	tmp = infile.readline().strip()
	L.append(float(tmp))
tmp = infile.readline().strip()
nSim = int(tmp)
msfile = infile.readline().strip()
infile.close()

a1_spA, a1_spB, a2_spA, a2_spB= [], [], [], []
for nsam in nSamA:
	a1_spA.append(sum([ 1.0/i for i in range(1, nsam) ]))
	a2_spA.append(sum([ 1.0/(i**2) for i in range(1, nsam) ]))
for nsam in nSamB:
	a1_spB.append(sum([ 1.0/i for i in range(1, nsam) ]))
	a2_spB.append(sum([ 1.0/(i**2) for i in range(1, nsam) ]))

# ms' output file
res = "dataset\tbialsites_avg\t"
res += "sf_avg\t"
res += "sxA_avg\t"
res += "sxB_avg\t"
res += "ss_avg\t"
res += "successive_ss_avg\t"
res += "piA_avg\t"
res += "piB_avg\t"
#res += "pearson_r_pi\t"
res += "thetaA_avg\t"
res += "thetaB_avg\t"
#res += "pearson_r_theta\t"
res += "DtajA_avg\t"
res += "DtajB_avg\t"
res += "divAB_avg\t"
res += "netdivAB_avg\t"
res += "minDivAB_avg\t"
res += "maxDivAB_avg\t"
res += "Gmin_avg\t"
res += "Gmax_avg\t"
res += "FST_avg\t"
#res += "pearson_r_divAB_netDivAB\t"
#res += "pearson_r_divAB_FST\t"
#res += "pearson_r_netDivAB_FST\t"
res += "ss_sf\t" # number of loci with both ss and sf
res += "ss_noSf\t" # number of loci with ss but no sf
res += "noSs_sf\t" # number of loci without ss but with sf
res += "noSs_noSf\n" # number of loci with no Ss or Sf

infile = open(msfile, "r")

test = 0 
nSim_cnt = 0 # count the number of treated multilocus simulations
nLoci_cnt = 0 # count the number of treated loci within a simulation
for line in infile:
	line = line.strip()
	if "segsites" in line:
		if nLoci_cnt == 0:
			ss_sf, noSs_sf, ss_noSf, noSs_noSf = 0, 0, 0, 0
			bialsites = []
			sf = []
			sxA = []
			sxB = []
			ss = []
			successive_ss = []
			piA = []
			piB = []
			thetaA = []
			thetaB = []
			DtajA = []
			DtajB = []
			divAB = []
			netdivAB = []
			minDivAB = []
			maxDivAB = []
			Gmin = []
			Gmax = []
			FST = []
		nLoci_cnt += 1
		nSam_cnt = 0 # count the number of treated individuals within a locus
		test = 1
		segsites = int(line.split(":")[1])
		bialsites.append(segsites)
		spA, spB = [], []
		continue
	if test == 1:
		if segsites == 0:
			test = 0
			bialsites.append(0)
			sf.append(0)
			sxA.append(0)
			sxB.append(0)
			ss.append(0)
			successive_ss.append(0)
			piA.append(0)
			piB.append(0)
			thetaA.append(0)
			thetaB.append(0)
			DtajA.append(0)
			DtajB.append(0)
			divAB.append(0)
			netdivAB.append(0)
			minDivAB.append(0)
			maxDivAB.append(0)
			Gmin.append(0)
			Gmax.append(0)
			FST.append(0)
		if segsites != 0:
			if "positions" not in line and line!="\n":
				nSam_cnt += 1
				if nSam_cnt <= nSamA[nLoci_cnt - 1]:
					spA.append(line.strip())
				if nSam_cnt > nSamA[nLoci_cnt - 1] and nSam_cnt <= (nSamA[nLoci_cnt - 1] + nSamB[nLoci_cnt - 1]):
					spB.append(line.strip())
				if nSam_cnt == (nSamA[nLoci_cnt - 1] + nSamB[nLoci_cnt - 1]):
					tmpA = compFreq(spA, segsites)
					freqA = tmpA['freq']
					piA.append(tmpA['pi']/L[nLoci_cnt - 1])
					
					tmpB = compFreq(spB, segsites)
					freqB = tmpB['freq']
					piB.append(tmpB['pi']/L[nLoci_cnt - 1])
					
					tmp = sites(freqA, freqB, segsites) # determines the # of sxA, sxB, ss, sf
					sf.append(tmp['sf'])
					sxA.append(tmp['sxA'])
					sxB.append(tmp['sxB'])
					ss.append(tmp['ss'])
					successive_ss.append(tmp['successive_ss'])
					# test if loci has both ss, sf, or only one of two, etc ...
					if tmp['sf'] > 0:
						if tmp['ss'] > 0:
							ss_sf += 1
						if tmp['ss'] == 0:
							noSs_sf += 1
					if tmp['sf'] == 0:
						if tmp['ss'] > 0:
							ss_noSf += 1
						if tmp['ss'] == 0:
							noSs_noSf += 1
					
					# theta = sx + ss
					thetaA_locus = (tmp['sxA']+tmp['ss'])/a1_spA[nLoci_cnt - 1]
					thetaB_locus = (tmp['sxB']+tmp['ss'])/a1_spB[nLoci_cnt - 1]
					thetaA.append(thetaA_locus/L[nLoci_cnt - 1 ])
					thetaB.append(thetaB_locus/L[nLoci_cnt - 1 ])
					
					# Taj
					DtajA.append(tajD(tmpA['pi'], thetaA_locus, nSamA[nLoci_cnt - 1], tmp['sxA']+tmp['ss'], a1_spA[nLoci_cnt-1], a2_spA[nLoci_cnt-1]))
					DtajB.append(tajD(tmpB['pi'], thetaB_locus, nSamB[nLoci_cnt - 1], tmp['sxB']+tmp['ss'], a1_spB[nLoci_cnt-1], a2_spB[nLoci_cnt-1]))
					
					# divAB
					div = compDiv(spA, spB, segsites)
					divAB.append(div['divAB']/L[nLoci_cnt - 1 ])
					netdivAB.append(div['divAB']/L[nLoci_cnt - 1 ] - (tmpA['pi']/L[nLoci_cnt - 1] + tmpB['pi']/L[nLoci_cnt - 1]) / 2.0)
					minDivAB.append(div['minDivAB']/L[nLoci_cnt - 1 ])
					maxDivAB.append(div['maxDivAB']/L[nLoci_cnt - 1 ])
					Gmin.append(div['Gmin'])
					Gmax.append(div['Gmax'])
					
					# Fst
					# vector of piT over segsites
					piT = piTot(tmpA['nDer'], tmpB['nDer'], nSamA[nLoci_cnt - 1], nSamB[nLoci_cnt - 1], segsites)
					# mean Fst
					FST.append(Fst(sum(tmpA['pi_SNPs']), sum(tmpB['pi_SNPs']), sum(piT), segsites))
					
	# compute average and std over of statistics over loci
	if nLoci_cnt != 0 and len(ss) == nLoci:
		test = 0
		nSim_cnt += 1
		nLoci_cnt = 0
		
		# statistics
		bialsites_avg = cr_mean(bialsites)
#		bialsites_std = cr_std(bialsites, bialsites_avg)
		sf_avg = cr_mean(sf)
#		sf_std = cr_std(sf, sf_avg)
		sxA_avg = cr_mean(sxA)
#		sxA_std = cr_std(sxA, sxA_avg)
		sxB_avg = cr_mean(sxB)
#		sxB_std = cr_std(sxB, sxB_avg)
		ss_avg = cr_mean(ss)
#		ss_std = cr_std(ss, ss_avg)
		successive_ss_avg = cr_mean(successive_ss)
#		successive_ss_std = cr_std(successive_ss, successive_ss_avg)
		piA_avg = cr_mean(piA)
#		piA_std = cr_std(piA, piA_avg)
		piB_avg = cr_mean(piB)
#		piB_std = cr_std(piB, piB_avg)
#		pearson_r_pi = cr_pearsonR(piA, piB)
		thetaA_avg = cr_mean(thetaA)
#		thetaA_std = cr_std(thetaA, thetaA_avg)
		thetaB_avg = cr_mean(thetaB)
#		thetaB_std = cr_std(thetaB, thetaB_avg)
#		pearson_r_theta = cr_pearsonR(thetaA, thetaB)
		DtajA_avg = cr_mean(DtajA)
#		DtajA_std = cr_std(DtajA, DtajA_avg)
		DtajB_avg = cr_mean(DtajB)
#		DtajB_std = cr_std(DtajB, DtajB_avg)
		divAB_avg = cr_mean(divAB)
#		divAB_std = cr_std(divAB, divAB_avg)
		netdivAB_avg = cr_mean(netdivAB)
#		netdivAB_std = cr_std(netdivAB, netdivAB_avg)
		minDivAB_avg = cr_mean(minDivAB)
#		minDivAB_std = cr_std(minDivAB, minDivAB_avg)
		maxDivAB_avg = cr_mean(maxDivAB)
#		maxDivAB_std = cr_std(maxDivAB, maxDivAB_avg)
		Gmin_avg = cr_mean(Gmin)
#		Gmin_std = cr_std(Gmin, Gmin_avg)
		Gmax_avg = cr_mean(Gmax)
#		Gmax_std = cr_std(Gmax, Gmax_avg)
		FST_avg = cr_mean(FST)
#		FST_std = cr_std(FST, FST_avg)
#		pearson_r_div_netDiv = cr_pearsonR(divAB, netdivAB)
#		pearson_r_div_FST = cr_pearsonR(divAB, FST)
#		pearson_r_netDiv_FST = cr_pearsonR(netdivAB, FST)
		
		#print("dataset {0}: {1} loci".format(nSim_cnt-1, len(ss)))
		res += "{0}\t{1:.5f}\t".format(datasetName, bialsites_avg)
		res += "{0:.5f}\t".format(sf_avg)
		res += "{0:.5f}\t".format(sxA_avg)
		res += "{0:.5f}\t".format(sxB_avg)
		res += "{0:.5f}\t".format(ss_avg)
		res += "{0:.5f}\t".format(successive_ss_avg)
		res += "{0:.5f}\t".format(piA_avg)
		res += "{0:.5f}\t".format(piB_avg)
#		res += "{0:.5f}\t".format(pearson_r_pi)
		res += "{0:.5f}\t".format(thetaA_avg)
		res += "{0:.5f}\t".format(thetaB_avg)
#		res += "{0:.5f}\t".format(pearson_r_theta)
		res += "{0:.5f}\t".format(DtajA_avg)
		res += "{0:.5f}\t".format(DtajB_avg)
		res += "{0:.5f}\t".format(divAB_avg)
		res += "{0:.5f}\t".format(netdivAB_avg)
		res += "{0:.5f}\t".format(minDivAB_avg)
		res += "{0:.5f}\t".format(maxDivAB_avg)
		res += "{0:.5f}\t".format(Gmin_avg)
		res += "{0:.5f}\t".format(Gmax_avg)
		res += "{0:.5f}\t".format(FST_avg)
#		res += "{0:.5f}\t".format(pearson_r_div_netDiv)
#		res += "{0:.5f}\t".format(pearson_r_div_FST)
#		res += "{0:.5f}\t".format(pearson_r_netDiv_FST)
		res += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}".format(ss_sf, ss_noSf, noSs_sf, noSs_noSf)
		res += "\n"
		if nSim_cnt == nSim:
			outfile = open("ABCstat_{0}_py.txt".format(datasetName), "w")
			outfile.write(res)
			outfile.close()
infile.close()
