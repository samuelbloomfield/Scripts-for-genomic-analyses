###
#	Read_searcher_python.py
#	Script to quantitate the number of SNP variants in processed read files
#	Author: Samuel Bloomfield
###



import re
import fileinput

isolate_count = 0

#Command for obtaining reverse complement of a sequence

def reverse_complement(seq):    
    for k,v in alt_map.iteritems():
        seq = seq.replace(k,v)
    bases = list(seq) 
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    for k,v in alt_map.iteritems():
        bases = bases.replace(v,k)
    return bases

for line in fileinput.input():
    print line
    o = re.findall('\S+\.\w+', line)
    p = o[0]
    q = o[1]
    r = o[2]
    s = o[3]

    count = 0

    A_list = []
    C_list = []
    G_list = []
    T_list = []
    SNP_list = []

    alt_map = {'ins':'0'}
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 

    with open (p) as e:
        for line in e:

            #Forms sequence variants for each kmer

            m = re.match('\w+\.\w+', line)
            
            SNP = m.group(0)
            
            SNP_list.append(SNP)

            R_SNP = reverse_complement(SNP)
            
            number_A = 0
            number_C = 0
            number_G = 0
            number_T = 0
            
            SNP_A =SNP.replace(".", "A")
            SNP_C =SNP.replace(".", "C")
            SNP_G =SNP.replace(".", "G")
            SNP_T =SNP.replace(".", "T")
            
            R_SNP_A = R_SNP.replace(".", "T")
            R_SNP_C = R_SNP.replace(".", "G")
            R_SNP_G = R_SNP.replace(".", "C")
            R_SNP_T = R_SNP.replace(".", "A")

            #Counts the number of each SNP variant in each processed read file

            with open (q) as f:
                for line in f:
                    if SNP_A in line: number_A +=1
                    if SNP_C in line: number_C +=1
                    if SNP_G in line: number_G +=1
                    if SNP_T in line: number_T +=1
                    
                    if R_SNP_A in line: number_A +=1
                    if R_SNP_C in line: number_C +=1
                    if R_SNP_G in line: number_G +=1
                    if R_SNP_T in line: number_T +=1

            with open (r) as g:
                for line in g:
                    if SNP_A in line: number_A +=1
                    if SNP_C in line: number_C +=1
                    if SNP_G in line: number_G +=1
                    if SNP_T in line: number_T +=1
                    
                    if R_SNP_A in line: number_A +=1
                    if R_SNP_C in line: number_C +=1
                    if R_SNP_G in line: number_G +=1
                    if R_SNP_T in line: number_T +=1
                                
            A_list.append(number_A)
            C_list.append(number_C)
            G_list.append(number_G)
            T_list.append(number_T)
            count +=1
            print count
            
    f = open(s, 'a')

    print>>f, ('\t'.join([str(x) for x in SNP_list]))
    print>>f, 'A_variant'
    print>>f, ('\t'.join([str(x) for x in A_list]))
    print>>f, 'C_variant'
    print>>f, ('\t'.join([str(x) for x in C_list]))
    print>>f, 'G_variant'
    print>>f, ('\t'.join([str(x) for x in G_list]))
    print>>f, 'T_variant'
    print>>f, ('\t'.join([str(x) for x in T_list]))
