'''
Generates all 2n-node k-color valise adinkras and computes all their gadgets.
Adam Artymowicz
'''
import math
import numpy as np
import time
import matplotlib.pyplot as plt
start_time = time.time()
data_filename = "adinkraData.npz"

def perm_parity(lst):
    '''\
    Given a permutation of the digits 0..N in order as a list, 
    returns its parity (or sign): +1 for even parity; -1 for odd.
    '''
    parity = 1
    for i in range(0,len(lst)-1):
        if lst[i] != i:
            parity *= -1
            mn = min(range(i,len(lst)), key=lst.__getitem__)
            lst[i],lst[mn] = lst[mn],lst[i]
    return parity    
    
#Saves the adinkras, vTildes, and the gadget matrix to a numpy file
def save(fname, data, printing = False):
    if printing:
        print('Writing to ' + fname)
    np.savez_compressed(fname, **data)

#Loads adinkras, vTildes, and gadget matrix from numpy file
def load(fname):
    print('Loading from ' + fname)
    return dict(np.load(fname))

#generates all permutations of length a
def permutations(a):
	if type(a) == int:
		assert a >= 0
		return permutations(range(a))

	if a == []:
		return [[]]

	out = []
	for i in range(len(a)):
		out = out + [[a[i]] + x for x in permutations(a[:i] + a[i+1:])]
	return out

def evenPermutations(a):
        def perm_parity2(a):
                a = list(a)
                b = sorted(a)
                inversions = 0
                while a:
                    first = a.pop(0)
                    inversions += b.index(first)
                    b.remove(first)
                return -1 if inversions % 2 else 1
    
        def permutations(a):
	       if type(a) == int:
		      assert a >= 0
		      return permutations(range(a))

	       if a == []:
		      return [[]]

	       out = []
	       for i in range(len(a)):
		      out = out + [[a[i]] + x for x in permutations(a[:i] + a[i+1:])]
	       return out
	
        z=[]
        for p in permutations(a):
            if perm_parity2(p)==1:
                z.append(p)
        return z

def oddPermutations(a):
         def perm_parity2(a):
                a = list(a)
                b = sorted(a)
                inversions = 0
                while a:
                    first = a.pop(0)
                    inversions += b.index(first)
                    b.remove(first)
                return -1 if inversions % 2 else 1
    
         def permutations(a):
	       if type(a) == int:
		      assert a >= 0
		      return permutations(range(a))

	       if a == []:
		      return [[]]

	       out = []
	       for i in range(len(a)):
		      out = out + [[a[i]] + x for x in permutations(a[:i] + a[i+1:])]
	       return out
	
         z=[]
         for p in permutations(a):
             if perm_parity2(p)==-1:
                 z.append(p)
         return z

#turns a permutation (as a valid list of integers) into a permutation matrix
def permutationMatrix(a):
	n = len(a)
	out = np.zeros((n,n))
	for i in range(n):
		out[i,a[i]] = 1
	return np.matrix(np.transpose(out))

#generates all matrices corresponding to color flips
def colorFlips(n):
	out = []
	for i in range(2**n):
		x = np.matrix(np.identity(n))
		for j in range(n):
			if i&2**j:
				x[j,j] = -1
		out.append(x)
	return out

#generates propagator for permuting fermions for 2n-vertex adinkra with k colors
def generateLeftPermutationPropagator(n, k):
	def permutationPropagator(seed):
		return [[permutationMatrix(p)*seed[i] for i in range(k)] 
				for p in permutations(n)]
	return permutationPropagator

def generateLeftEvenPermutationPropagator(n, k):
	def permutationPropagator(seed):
		return [[permutationMatrix(p)*seed[i] for i in range(k)] 
				for p in evenPermutations(n)]
	return permutationPropagator

def generateLeftOddPermutationPropagator(n, k):
	def permutationPropagator(seed):
		return [[permutationMatrix(p)*seed[i] for i in range(k)] 
				for p in oddPermutations(n)]
	return permutationPropagator
	
#generates propagator for permuting fermions for 2n-vertex adinkra with k colors
def generateRightPermutationPropagator(n, k):
	def permutationPropagator(seed):
		return [[seed[i]*permutationMatrix(p) for i in range(k)] 
				for p in permutations(n)]
	return permutationPropagator

def generateRightEvenPermutationPropagator(n, k):
	def permutationPropagator(seed):
		return [[seed[i]*permutationMatrix(p) for i in range(k)] 
				for p in evenPermutations(n)]
	return permutationPropagator

def generateRightOddPermutationPropagator(n, k):
	def permutationPropagator(seed):
		return [[seed[i]*permutationMatrix(p) for i in range(k)] 
				for p in oddPermutations(n)]
	return permutationPropagator
	
#generates propagator for flipping fermions for 2n-vertex adinkra with k colors
def generateLeftFlipPropagator(n,k):
	def flipPropagator(seed):
		return [[flip*seed[i] for i in range(k)] for flip in colorFlips(n)]
	return flipPropagator

#generates propagator for flipping fermions for 2n-vertex adinkra with k colors
def generateRightFlipPropagator(n,k):
	def flipPropagator(seed):
		return [[seed[i]*flip for i in range(k)] for flip in colorFlips(n)]
	return flipPropagator

#initializes seed L-matrices for 2n-node, k-color adinkras
def initializeSeedMatrices(n,k):
	if (n,k) == (4,4):
		seed_1 = [np.matrix(np.identity(n), dtype = 'int64')]
		seed_1.append(np.matrix([[0,0,1,0],[0,0,0,1],[-1,0,0,0],[0,-1,0,0]]))
		seed_1.append(np.matrix([[0,1,0,0],[-1,0,0,0],[0,0,0,-1],[0,0,1,0]]))
		seed_1.append(np.matrix([[0,0,0,1],[0,0,-1,0],[0,1,0,0],[-1,0,0,0]]))

		seed_2 = list(seed_1)
		seed_2[3] = -1*seed_2[3]
		return [seed_1,seed_2]

	if (n,k) == (4,3):
		seed = [np.matrix(np.identity(n), dtype = 'int64')]
		seed.append(np.matrix([[0,0,1,0],[0,0,0,1],[-1,0,0,0],[0,-1,0,0]]))
		seed.append(np.matrix([[0,1,0,0],[-1,0,0,0],[0,0,0,-1],[0,0,1,0]]))
		return [seed]
	
	if (n,k) == (2,2):
	       seed = [np.matrix([[0, -1], [1, 0]], dtype = 'int64')]
	       seed.append(np.matrix([[-1, 0], [0, -1]]))
	       return [seed]
        if (n,k) == (8,4):
               seed = [np.matrix(np.identity(n), dtype = 'int64')]
               seed.append(np.matrix([[0,1,0,0,0,0,0,0], [-1,0,0,0,0,0,0,0], [0,0,0,-1,0,0,0,0], [0,0,1,0,0,0,0,0], [0,0,0,0,0,-1,0,0], [0,0,0,0,1,0,0,0], [0,0,0,0,0,0,0,1], [0,0,0,0,0,0,-1,0]]))
               seed.append(np.matrix([[0,0,1,0,0,0,0,0], [0,0,0,1,0,0,0,0], [-1,0,0,0,0,0,0,0], [0,-1,0,0,0,0,0,0], [0,0,0,0,0,0,-1,0], [0,0,0,0,0,0,0,-1], [0,0,0,0,1,0,0,0], [0,0,0,0,0,1,0,0]]))
               seed.append(np.matrix([[0,0,0,0,1,0,0,0], [0,0,0,0,0,1,0,0], [0,0,0,0,0,0,1,0], [0,0,0,0,0,0,0,1], [-1,0,0,0,0,0,0,0], [0,-1,0,0,0,0,0,0], [0,0,-1,0,0,0,0,0], [0,0,0,-1,0,0,0,0]]))
               return [seed]
#generates the orbit of seeds under propagators
def generateAdinkras(seeds, propagators, n,k, removeDuplicates = False):
	#print "number of propagators: " + str(len(propagators))
	if propagators == []:
		return seeds

	if not removeDuplicates:
		out = []
		for seed in seeds:
			out = out + propagators[0](seed)
	else:
		out = set()
		for i in range(len(seeds)):
			if i%10 == 0:
				print("propagating %i / %i" % (i,len(seeds)))
			toAdd = propagators[0](seeds[i])
			for adinkra in toAdd:
				out.add(npToTuple(adinkra))


		out = list(out)
		out = [tupleToNp(adinkra,n,k) for adinkra in out]

	return generateAdinkras(out, propagators[1:], n,k, removeDuplicates)

#computes an ordered set of vTilde matrices given an ordered set of L-matrices
def vTilde(adinkra):
	k = len(adinkra)
	n = len(adinkra[0])
	W = [[None]*k for i in range(k)]
	for i in range(k):
		for j in range(i+1,k):
			x = 0.5*(np.transpose(adinkra[i])*adinkra[j] - 
				np.transpose(adinkra[j])*adinkra[i])
			W[i][j] = x
			W[j][i] = -x
		W[i][i] = np.zeros((n,n))
	return W

#computes the gadget given two sets of V tilde matrices
def gadget(W_1, W_2):
	n = len(W_1)
	s = 0

	for i in range(n):
		for j in range(i+1,n):
			s += np.trace(W_1[i][j]*W_2[i][j])
	return 2*s

#given a matrix M, returns all the values that matrix takes on
def valuesInMatrix(M):
	s = set()
	for row in M:
		for element in row:
			if element not in s:
				s.add(element)
	return list(s)

#Checks if a list of adinkras contains duplicates

def areDistinct(adinkras, breakIfNondistinct = True, alertIfDistinct = True, printing = True):
	n = len(adinkras)
	#print 'total number of adinkras ' + str(n)
	m = len(adinkras[0])
	for i in range(n):
		for j in range(i+1, n):
			if j == n-1:
				print 'checking' + str((i,j))
			adinkra_comparison = [(np.matrix(adinkras[i][k]) == 
				np.matrix(adinkras[j][k])).all() for k in range(m)]
			if all(adinkra_comparison):
				if not alertIfDistinct:
					print('adinkras %i and %i are the same!' % (i,j))
				if breakIfNondistinct:
					return False
			else:
				if alertIfDistinct:
					print('adinkras %i and %i are not the same!' % (i,j))


#checks if two adinkras are the same
def sameAdinkra(adinkra_1, adinkra_2, k):
	adinkra_comparison = [(np.matrix(adinkra_1[i]) == 
							np.matrix(adinkra_2[i])).all() for i in range(k)]
	return all(adinkra_comparison)

#converting between adinkra representations
def npToTuple(adinkra):
	n,m = adinkra[1].shape
	out = []
	for i in range(len(adinkra)):
		out = out + adinkra[i].reshape((n*n,)).tolist()[0]
	return tuple(out)

#converting between adinkra representations
def tupleToNp(adinkra_tuple, n, k):
	assert len(adinkra_tuple) == n*n*k
	out = []
	for i in range(0,n*n*k,n*n):
		L = list(adinkra_tuple[i:i+n*n])
		out.append(np.matrix(np.reshape(L, (n,n))))
	return out

#gives a histogram of gadget values given a gadget matrix
def gadgetSummary(gadget_matrix):
	gadget_values = {}
	for entry in np.nditer(gadget_matrix):
		try:
			gadget_values[float(entry)] += 1
		except KeyError:
			gadget_values[float(entry)] = 1
	return gadget_values

#verifies that all color permutations of the seeds lie in the orbit
def verifyRainbowInvariance(adinkras, seeds, n,k):
	print 'verifying \'rainbow-invariance\''
	permuted_adinkras = []
	for seed in seeds:
		for p in permutations(k):
			permuted_adinkra = [seed[p[i]] for i in range(k)]
			permuted_adinkras.append(permuted_adinkra)
	
	for p_adinkra in permuted_adinkras:
		distinct = True
		for adinkra in adinkras:
			if sameAdinkra(adinkra, p_adinkra, k):
				distinct = False
				break
		if distinct:
			print 'This color-permuted adinkra does not lie in the orbit:'
			print p_adinkra
			return False

	print "Adinkra set is rainbow-invariant"
	return True

def miniMain(n,k, seeds = None):
    seed = initializeSeedMatrices(n,k)[0]
    m=math.factorial(n)*(2**n)
    gadgets = np.zeros((m))
    j=0
    print len(permutations(n))*len(colorFlips(n))
    for p in permutations(n):
        for flip in colorFlips(n):
           if j%50 == 0:
               print j
           a=[seed[i]*permutationMatrix(p)*flip for i in range(k)]
           gadgets[j] = gadget(vTilde(a), vTilde(seed))
           j=j+1
    return gadgetSummary(gadgets)

def miniMain2(n,k, seeds = None):
    seed = initializeSeedMatrices(n,k)[0]
    m=math.factorial(n)*(2**n)
    gadgets = np.zeros((m))
    j=0
    print len(permutations(n))*len(colorFlips(n))
    for p in permutations(n):
        for flip in colorFlips(n):
           if j%50 == 0:
               print j
           a=[seed[i]*permutationMatrix(p)*flip for i in range(k)]
           gadgets[j] = gadget(vTilde(a), vTilde(seed))
           if gadgets[j] == -24 or gadgets[j] == 24:
               print j
               print a
               print gadgets[j]
               return gadgetSummary(gadgets)   
           j=j+1
    
#computes gadgets between all valise adinkras with 2n vertices and k colors
def main(n,k, propagators, mode = "parallel", seeds = None):
	d = {'leftPermutations':generateLeftPermutationPropagator(n,k), 
	        'leftEvenPermutations':generateLeftEvenPermutationPropagator(n,k),
	        'leftOddPermutations':generateLeftOddPermutationPropagator(n,k),
		'rightPermutations':generateRightPermutationPropagator(n,k),
		'rightEvenPermutations':generateRightEvenPermutationPropagator(n,k),
		'rightOddPermutations':generateRightOddPermutationPropagator(n,k),
		'leftFlips':generateLeftFlipPropagator(n,k),
		'rightFlips':generateRightFlipPropagator(n,k)}
        
	#if no seed matrices were given, use the standard ones for n,k
	if seeds == None and mode != "serial":
		seeds = initializeSeedMatrices(n,k)

	#generate adinkras
	if mode == "parallel":
		propagator_list = [d[p] for p in propagators]
 		adinkras = generateAdinkras(seeds, propagator_list,n,k, 
 									removeDuplicates = True)

	#so far this mode supports only 2 propagators
	elif mode == "serial":
		assert len(propagators) == 2
		gadget_matrices = []
		assert seeds == None
		propagator = d[propagators[0]]
		seeds = propagator(initializeSeedMatrices(n,k)[0])
		for seed in seeds:
			gadget_matrices.append(main(n,k, propagators[1:], 
						mode = "parallel", seeds = [seed]))

		gadget_comparison = [(gadget_matrices[i] == gadget_matrices[i+1]).all() 
			for i in range(len(gadget_matrices) - 1)]
		if all(gadget_comparison):
			print('All of the %i x %i gadget matrices are the same' 
				% gadget_matrices[0].shape)
		return

	vTildes = [vTilde(a) for a in adinkras]
	m = len(vTildes)
	gadgets = np.zeros((m))

	print

	print('calculating gadget values')
	for i in range(m):
		if i%10 == 0:
			print('computing gadgets: %i / %i' % (i, len(vTildes)))
		gadgets[i] = gadget(vTildes[i], vTildes[1])
	print

	#verifyRainbowInvariance(adinkras, seeds, n,k)
	print

	print("--- %s seconds ---" % (time.time() - start_time))
	print

	data = {'adinkras': np.asarray(adinkras), 
			'vTildes': np.asarray(vTildes), 
			'gadgets': gadgets}
	save(data_filename, data, True)
	print
	#areDistinct(adinkras, breakIfNondistinct = True, alertIfDistinct = False, printing = True)
	return gadgetSummary(gadgets)
	           
#print "even" 
#print len(evenPermutations(4))
#print "odd" 
#print len(oddPermutations(4))
#print main(4,4, propagators = ["rightPermutations", "rightFlips"], mode = "parallel")
#print main(4,3, propagators = ["rightPermutations","rightFlips"], mode = "parallel")
print miniMain2(8, 4)
#print initializeSeedMatrices(8,4)[0]
#print main(2, 2, propagators = ["rightPermutations", "rightFlips"], mode = "parallel")
#seeds = initializeSeedMatrices(4,4)
#print sameAdinkra(seeds[0], seeds[1], 4)		