import sys
import numpy as np

def globalAlignment(s, t, mismatch, indel):
    #matrices
    m = len(s) #length of s string
    n = len(t) #length of t string
    V = np.zeros((m+1, n+1), dtype = int)
    V[0, 0] = 0
    #
    for j in range(1, n+1):
        V[0, j] = j * indel
    for i in range(1, m+1):
        V[i, 0] = i * indel
    S = np.empty((m+1, n+1), dtype = str)
    for i in range(1, m+1):
        for j in range(1, n+1):

            if i == 0 or j == 0:
                #V[i,j] = 0
                S[i,j] = '#'
            elif s[i-1] is t[j-1]: #if the characters match
                S[i,j] = 'F' #S[i-1,j-1] + s[i-1]
                V[i,j] = V[i-1,j-1] + 1

            elif (V[i-1,j-1]+mismatch) > (max(V[i,j-1], V[i-1,j]) + indel): #if mismatch
                V[i,j] = V[i-1,j-1]+mismatch
                S[i,j] = 'M' #S[i-1,j-1]

            elif V[i-1,j] > V[i,j-1]: #if i shift
                V[i,j] = V[i-1,j] + indel
                S[i,j] = 'R' #S[i-1,j] #+ s[i-1] # same here with s[i-1]
            else: #else j shift
                V[i,j] = V[i,j-1] + indel
                S[i,j] += 'L' #S[i,j-1] #+ t[j-1] #check if t[j-1] should be here

    c = V[i,j]
    a,b = getGlobalAlignments(S, s, t, m, n)
    return c,a,b
        


def getGlobalAlignments(S, s, t, i, j):
    #if i>0 and j>0:
    if S[i,j] == '#':
        return '',''
    elif 'F'== S[i,j]: #if match
        a,b = getGlobalAlignments(S, s, t, i-1, j-1)
        return a + s[i-1], b + t[j-1]
    
    elif 'M' == S[i,j]:
        a,b = getGlobalAlignments(S, s, t, i-1, j-1)
        return a + s[i-1], b + t[j-1]
    
    elif i>0 and S[i,j] == 'R':#if i shift
        a,b = getGlobalAlignments(S, s, t, i-1, j)
        return a + s[i-1] , b + '-'
    
    elif j>0 and S[i,j] == 'L':#if j shift
        a,b = getGlobalAlignments(S, s, t, i, j-1)
        return a + '-', b + t[j-1]
    
    else:
        return '',''

    


def localAlignment(s, t, mismatch, indel):
    #matrices
    m = len(s)
    n = len(t)
    V = np.zeros((m+1, n+1), dtype = int)
    S = np.empty((m+1, n+1), dtype = str)
    for i in range(1,m+1):
        for j in range(1,n+1):

            if i == 0 or j == 0: #zero case
                V[i,j] = 0
                S[i,j] = '#'

            elif (0 >= max(V[i-1,j-1] + 1, V[i,j-1] + indel, V[i-1,j] + indel, V[i-1,j-1] + mismatch)): #backtrdacking end marker
                V[i,j] = 0
                S[i,j] = '#'

            elif (V[i-1,j-1]+mismatch) > (max(V[i,j-1] + indel, V[i-1,j] + indel, V[i-1,j-1]+1)): #if mismatch
                V[i,j] = V[i-1,j-1] + mismatch
                S[i,j] = 'M' 
        
            elif s[i-1] is t[j-1]: #match
                S[i,j] = 'F'
                V[i,j] = V[i-1,j-1] + 1

            elif V[i,j-1] + indel > V[i-1,j] + indel: #j shift
                V[i,j] = V[i,j-1] + indel
                S[i,j] = 'L'
            else: #i shift
                V[i,j] = V[i-1,j] + indel
                S[i,j] = 'R'

            #refactor this algo, same as global but 0>max at the end in its own if statemnet 
    
    c = 0
    i_max = 0
    j_max = 0
    for x in range(1,m+1):
        for y in range(1,n+1):
            if V[x,y] > c:
                c = V[x,y]
                i_max = x
                j_max = y


    a,b = getLocalAlignment(S, s, t, i_max, j_max)
    return c,a,b




def getLocalAlignment(S, s, t, i, j):
    #if i>0 and j>0:
    if S[i,j] == '#':
        return '',''
    elif 'F'== S[i,j]: #if match
        a,b = getLocalAlignment(S, s, t, i-1, j-1)
        return a + s[i-1], b + t[j-1]
    
    elif 'M' == S[i,j]:
        a,b = getLocalAlignment(S, s, t, i-1, j-1)
        return a + s[i-1], b + t[j-1]
    
    elif i>0 and S[i,j] == 'R':#if i shift
        a,b = getLocalAlignment(S, s, t, i-1, j)
        return a + s[i-1] , b + '-'
    
    elif j>0 and S[i,j] == 'L':#if j shift
        a,b = getLocalAlignment(S, s, t, i, j-1)
        return a + '-', b + t[j-1]
    
    else:
        return '',''


def readFASTA(filename):
    """Read a FASTA file containing a single sequence and return
       the sequence as a string.
    Parameter:
        filename: a string representing the name of a FASTA file
    Return value: a string containing the sequence in the FASTA file
    """
    inputFile = open(filename, 'r', encoding = 'utf-8')
    header = inputFile.readline()
    dna = ''
    for line in inputFile:
        dna = dna + line[:-1]
    inputFile.close()
    return dna


def main():
    sys.setrecursionlimit(1500)
    print("Reading sequence from FASTA file...")
    sequence1 = readFASTA('eyeless_genomic_sequence.fasta')
    sequence2 = readFASTA('PAX6_genomic_sequence.fasta')


    s = sequence1[:1000]
    t = sequence2[:1000]


    mismatch = -2
    indel = -3


    print("Performing global alignment...")
    global_score, global_alignment1, global_alignment2 = globalAlignment(s, t, mismatch, indel)
    print("Global Alignment Score:", global_score)
    print("Global Alignment S:", global_alignment1)
    print("Global Alignment T:", global_alignment2)

    print("Performing local alignment...")
    local_score, local_alignment1, local_alignment2 = localAlignment(s, t, mismatch, indel)
    print("Local Alignment Score:", local_score)
    print("Local Alignment S:", local_alignment1)
    print("Local Alignment T:", local_alignment2)



main()