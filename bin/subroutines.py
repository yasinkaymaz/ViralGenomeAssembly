import numpy
from Bio.Align import AlignInfo
from scipy import stats
import pandas as pd
import sys
import collections
from Bio import AlignIO
import Bio.Align



class GeneticDist:
    def __init__(self, seq1, seq2):
        self.seq1 = seq1
        self.seq2 = seq2

    def estimate_nucleotide_frequencies(self):
        seq = self.seq1.replace('-','').upper()
        A = seq.count('A')
        C = seq.count('C')
        G = seq.count('G')
        T = seq.count('T')
        length = float(len(seq))
        return [ x/length for x in [A,C,G,T] ]

    def pdistance(self):
        p = 0
        pairs = []
        for x in zip(self.seq1,self.seq2):
            if '-' not in x: pairs.append(x)
        #for (x,y) in zip(seq1,seq2):
        for (x,y) in pairs:
            if x != y:
                p += 1
        #length = (len(seq1) + len(seq2)) / 2
        length = len(pairs)
        return float(p) / length

    def JCdistance(self):
        """
        distance = -b log(1 - p / b)
        where:
        b = 3/4
        and p = p-distance, i.e. uncorrected distance between seq1 and seq2
        """
        from math import log
        b = 0.75
        p = pdistance(self.seq1,self.seq2)
        try: d = -b * log( 1-p/b )
        except ValueError:
            print "Tried to take log of a negative number"
            return None
        return d

    def TNdistance(self):
        """
        Tajima-Nei distance = -b log(1 - p / b)
        where:
        b = 0.5 * [ 1 - Sum i from A to T(Gi^2+p^2/h) ]
        h = Sum i from A to G( Sum j from C to T (Xij^2/2*Gi*Gj))
        p = p-distance, i.e. uncorrected distance between seq1 and seq2
        Xij = frequency of pair (i,j) in seq1 and seq2, with gaps removed
        Gi = frequency of base i over seq1 and seq2 """
        from math import log

        ns = ['A','C','G','T']
        G = estimate_nucleotide_frequencies(self.seq1 + self.seq2)
        p = pdistance(self.seq1,self.seq2)
        pairs = []
        h = 0

        #collect ungapped pairs
        for x in zip(self.seq1,self.seq2):
            if '-' not in x: pairs.append(x)

        #pair frequencies are calculated for AC, AG, AT, CG, CT, GT (and reverse order)
        for i in range(len(ns)-1):
            for j in range(i+1,len(ns)):
                if i != j: paircount = pairs.count( (ns[i], ns[j]) ) + pairs.count( (ns[j], ns[i]) )
                Xij_sq = (float(paircount)/len(pairs))**2
                GiGj = G[i]*G[j]
                h += 0.5*Xij_sq/GiGj  #h value used to calculate b

        b = 0.5*(1-sum([x**2 for x in G])+p**2/h)
        try: d = -b * log(1 - p/b)
        except ValueError:
            print "Tried to take log of a negative number"
            return None
        return d

    def K2Pdistance(self):
        """
        Kimura 2-Parameter distance = -0.5 log( (1 - 2p -q) * sqrt( 1 - 2q ) )
        where:
        p = transition frequency
        q = transversion frequency
        """
        from math import log, sqrt
        pairs = []

        #collect ungapped pairs
        for x in zip(self.seq1,self.seq2):
            if '-' not in x: pairs.append(x)

        ts_count=0
        tv_count=0
        length = len(pairs)

        transitions = [ "AG", "GA", "CT", "TC"]
        transversions = [ "AC", "CA", "AT", "TA",
                          "GC", "CG", "GT", "TG" ]

        for (x,y) in pairs:
            if x+y in transitions: ts_count += 1
            elif x+y in transversions: tv_count += 1

        p = float(ts_count) / length
        q = float(tv_count) / length
        try: d = -0.5 * log( (1 - 2*p - q) * sqrt( 1 - 2*q ) )
        except ValueError:
            print "Tried to take log of a negative number"
            return None
        return d

    def Tamuradistance(self):
        """
        Tamura distance = -C log( 1 - P/C - Q ) - 0.5( 1 - C )log( 1 - 2Q )
        where:
        P = transition frequency
        Q = transversion frequency
        C = GC1 + GC2 - 2 * GC1 * GC2
        GC1 = GC-content of sequence 1
        GC2 = GC-coontent of sequence 2
        """
        from math import log
        pairs = []

        #collect ungapped pairs
        for x in zip(self.seq1,self.seq2):
            if '-' not in x: pairs.append(x)

        ts_count=0
        tv_count=0
        length = len(pairs)

        transitions = [ "AG", "GA", "CT", "TC"]
        transversions = [ "AC", "CA", "AT", "TA",
                          "GC", "CG", "GT", "TG" ]

        for (x,y) in pairs:
            if x+y in transitions: ts_count += 1
            elif x+y in transversions: tv_count += 1

        p = float(ts_count) / length
        q = float(tv_count) / length
        gc1 = sum(estimate_nucleotide_frequencies(self.seq1)[1:3])
        gc2 = sum(estimate_nucleotide_frequencies(self.seq2)[1:3])
        c = gc1 + gc2 - 2 * gc1 * gc2

        try: d = -c * log( 1 - p/c - q) - 0.5 * ( 1 - c ) * log ( 1 - 2*q )
        except ValueError:
            print "Tried to take log of a negative number"
            return None
        return d



# FUNCTION DEFS:
def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    """ Check if two values are almost equal """
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def geneticCode(name):

    """ Dictionary that maps codons to amino acids """
    gc = {  'AAA':'K', 'AAC':'N', 'AAG':'K', 'AAT':'N', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 'AGA':'R', 'AGC':'S', 'AGG':'R', \
            'AGT':'S','ATA':'I','ATC':'I','ATG':'M','ATT':'I','CAA':'Q','CAC':'H','CAG':'Q','CAT':'H','CCA':'P','CCC':'P','CCG':'P', \
            'CCT':'P','CGA':'R','CGC':'R','CGG':'R','CGT':'R','CTA':'L','CTC':'L','CTG':'L','CTT':'L','GAA':'E','GAC':'D','GAG':'E', \
            'GAT':'D','GCA':'A','GCC':'A','GCG':'A','GCT':'A','GGA':'G','GGC':'G','GGG':'G','GGT':'G','GTA':'V','GTC':'V','GTG':'V', \
            'GTT':'V','TAA':'*','TAC':'Y','TAG':'*','TAT':'Y','TCA':'S','TCC':'S','TCG':'S','TCT':'S','TGA':'*','TGC':'C','TGG':'W', \
            'TGT':'C','TTA':'L','TTC':'F','TTG':'L','TTT':'F'  }
    if name == 'standard':
        table = gc

    elif name == 'withmissing':
        table = gc
        for codon in table.keys():
            table['N'+codon[1:3]] = "*"
            table['NN'+codon[2:3]] = "*"
            table[codon[0]+'N'+codon[2]] = "*"
            table['N'+codon[1]+'N'] = "*"
            table[codon[0]+'NN'] = "*"
            table[codon[0:2]+'N'] = "*"
            table['NNN'] = "*"
    return table


def potential_changes_dict(nt_to_aa):

    """ Generate a dictionary, with S and N pre-calculated for all
    possible pairs of codons (key: pair of codons, value: (S,N).
    ARGS:
        nt_to_aa, a dict mapping codons (keys), e.g. 'TTA', to
            amino-acid letter (values), e.g. 'L'

            e.g. geneticCode("standard")
    Notes:
        Sources of formulae:
        http://www.megasoftware.net/mega4/WebHelp/part_iv___evolutionary_analysis/computing_evolutionary_distances/distance_models/synonymouse_and_nonsynonymous_substitution_models/hc_nei_gojobori_method.htm
    """
    # IMPORTS:
    import numpy as np
    import pickle
    from  itertools import permutations
    import copy # @todo: remove the deepcopying?
    import pdb

    potential_changes = {   'S': {  'AAA':0.0,'AAC':0.0,'AAG':0.0,'AAT':0.0, 'ACA':0.0, 'ACC':0.0, 'ACG':0.0, 'ACT':0.0, 'AGA':0.0, 'AGC':0.0, \
                                    'AGG':0.0,'AGT':0.0, 'ATA':0.0, 'ATC':0.0, 'ATG':0.0, 'ATT':0.0, 'CAA':0.0, 'CAC':0.0, 'CAG':0.0, 'CAT':0.0, \
                                    'CCA':0.0,'CCC':0.0,'CCG':0.0,'CCT':0.0,'CGA':0.0,'CGC':0.0,'CGG':0.0,'CGT':0.0,'CTA':0.0,'CTC':0.0,'CTG':0.0, \
                                    'CTT':0.0,'GAA':0.0,'GAC':0.0,'GAG':0.0,'GAT':0.0,'GCA':0.0,'GCC':0.0,'GCG':0.0,'GCT':0.0,'GGA':0.0,'GGC':0.0, \
                                    'GGG':0.0,'GGT':0.0,'GTA':0.0,'GTC':0.0,'GTG':0.0,'GTT':0.0,'TAA':0.0,'TAC':0.0,'TAG':0.0,'TAT':0.0,'TCA':0.0, \
                                    'TCC':0.0,'TCG':0.0,'TCT':0.0,'TGA':0.0,'TGC':0.0,'TGG':0.0,'TGT':0.0,'TTA':0.0,'TTC':0.0,'TTG':0.0,'TTT':0.0},

                            'N': {  'AAA':0.0,'AAC':0.0,'AAG':0.0,'AAT':0.0,'ACA':0.0,'ACC':0.0,'ACG':0.0,'ACT':0.0,'AGA':0.0,'AGC':0.0,'AGG':0.0, \
                                    'AGT':0.0,'ATA':0.0,'ATC':0.0,'ATG':0.0,'ATT':0.0,'CAA':0.0,'CAC':0.0,'CAG':0.0,'CAT':0.0,'CCA':0.0,'CCC':0.0,'CCG':0.0, \
                                    'CCT':0.0,'CGA':0.0,'CGC':0.0,'CGG':0.0,'CGT':0.0,'CTA':0.0,'CTC':0.0,'CTG':0.0,'CTT':0.0,'GAA':0.0,'GAC':0.0,'GAG':0.0, \
                                    'GAT':0.0,'GCA':0.0,'GCC':0.0,'GCG':0.0,'GCT':0.0,'GGA':0.0,'GGC':0.0,'GGG':0.0,'GGT':0.0,'GTA':0.0,'GTC':0.0,'GTG':0.0, \
                                    'GTT':0.0,'TAA':0.0,'TAC':0.0,'TAG':0.0,'TAT':0.0,'TCA':0.0,'TCC':0.0,'TCG':0.0,'TCT':0.0,'TGA':0.0,'TGC':0.0,'TGG':0.0, \
                                    'TGT':0.0,'TTA':0.0,'TTC':0.0,'TTG':0.0,'TTT':0.0}}


    # Mutate (substitutions) all possible codons in the given genetic code, and count proportions of mutations that are synonymous and non-synonmyous
    for codon in nt_to_aa.keys():

        # assert (codon not in codon_record)  @DONE: no duplicate entries
        # codon_record.append(codon)

        # Calculate S and N (i.e. potential synonymous and potential
        # non-synonymous sites) ()

        # i.e. Of possible mutations that can occur on the codon,
        # what proportion are synonymous (no aa change)?

        # for each bp position in the codon...
        for codon_p in range(0,2+1):

            nts = ['A','G','T','C']  # @DONE: refactor away, A: we can't, since the next line

            nts.remove(codon[codon_p]) # we do not consider self substitutions, e.g. A->A

            # ...and for each nucleotide that the bp can change
            # into...
            for nt in nts:

                codon_mutated = list(copy.deepcopy(codon))
                #codon_mutated = codon
                codon_mutated[codon_p] = nt  # mutate the basepair
                codon_mutated = ''.join(codon_mutated)

                # ...count how many of them are synonymous.
                if nt_to_aa[codon]==nt_to_aa[codon_mutated]:
                    potential_changes['S'][codon]+=1.0/3.0 #@DONE: Q: but why 1/3? to me it should be 1/(3*4), A: since it represents 1/3 of a "site"
                else:
                    potential_changes['N'][codon]+=1.0/3.0 #@DONE: Q: but why 1/3? to me it should be 1/(3*4), A: since it represents 1/3 of a "site"

    codons      = nt_to_aa.keys()
    codonPairs  = list(permutations(codons,2))
    selfies     = [(i,i) for i in codons]
    codonPairs  = codonPairs + selfies

    codonPair_to_potential = {}

    for pair in codonPairs:

        codon1 = pair[0]
        codon2 = pair[1]
        pn1 = potential_changes['N'][codon1]
        pn2 = potential_changes['N'][codon2]
        ps1 = potential_changes['S'][codon1]
        ps2 = potential_changes['S'][codon2]

        codonPair_to_potential[pair] = {'N':(pn1+pn2)/2.,'S':(ps1+ps2)/2.}

    return codonPair_to_potential

def observed_changes_dict(nt_to_aa):

    """ Generate a dictionary, with Sd and Nd pre-calculated for all
    possible pairs of codons (key: pair of codons, value: (Sd,Nd).
    ARGS:
        nt_to_aa, a dict mapping codons (keys), e.g. 'TTA', to
            amino-acid letter (values), e.g. 'L'

            e.g. geneticCode("standard")
    Notes:
        Sources of formulae:
        http://www.megasoftware.net/mega4/WebHelp/part_iv___evolutionary_analysis/computing_evolutionary_distances/distance_models/synonymouse_and_nonsynonymous_substitution_models/hc_nei_gojobori_method.html
    """
    # IMPORTS:
    import numpy as np
    import pickle
    from  itertools import permutations
    import copy # @todo: remove the deepcopying?
    import pdb

    codons      = nt_to_aa.keys()
    codonPairs  = list(permutations(codons,2))
    selfies     = [(i,i) for i in codons]
    codonPairs  = codonPairs + selfies

    codonPair_to_observed = {}

    for pair in codonPairs:
        codon1 = pair[0]
        codon2 = pair[1]
        indices_to_permute = []

        # Collect the position of the letters (1, 2, 3) where the two sequences differ...
        for letter_i in range(0,3):
            if not codon1[letter_i] == codon2[letter_i]:
                indices_to_permute.append(letter_i)

        # We now have all the possible mutational pathways, represented as indices
        permuted_indices = list(permutations(indices_to_permute))
        syn = []
        non = []

        for i,path in enumerate(permuted_indices):
            syn.append(int())
            non.append(int())

            codon1_path1 = list(codon1) # copies of seqs for 'mutating'
            codon2_path1 = list(codon2)

            for site in path:

                codon1_past         = ''.join(codon1_path1)
                codon1_path1[site]  = codon2_path1[site]        # s1 = 'TTT' , s2 = 'ATA'  ==> 'TTT' --> 'ATT'
                codon1_path1        = ''.join(codon1_path1)
                #@comparison-step mutants successively
                if nt_to_aa[codon1_path1] == nt_to_aa[codon1_past]:  # 'TTT --> 'ATT'
                    syn[i] = syn[i] + 1
                    non[i] = non[i] + 0
                else:
                    syn[i] = syn[i] + 0
                    non[i] = non[i] + 1
                codon1_path1 = list(codon1_path1)

        try:
            assert isclose(np.mean(syn)+np.mean(non),float(len(path)))
        except AssertionError:
            # pdb.set_trace()
            raise ValueError("Calculations are incorrect, mutation pathways calculation failed...")


        codonPair_to_observed[pair] = {'S':np.mean(syn),'N':np.mean(non)}

    return codonPair_to_observed

def align_query_vs_reference( qry_seq, ref_seq, aln_gap_open = -10, aln_gap_extend = -0.5 ):

   """ Globally align two input sequences using NW-algorithm, return highest-scoring alignment.
   ARGS:

      qry_seq, ref_seq,    two different orthologous CDS (exon-only) DNA
                           sequences to be aligned (see: "** Input sequences").
         e.g. qry_seq = s1_AGAP010815_RA = "ATGTGGCAGTTCATAAGGTCACGAATATTAACGGTGATAATCTTCATAGGTGCTGCTCATGGGCTACTGGTTGTGGGTCCGAAATTTATACGGGCCAACCAGGAATACACTCTGGTGATCAGCAACTTTAACTCACAGCTAAGCAAAGTGGACCTGCTGTTAAAACTGGAAGGCGAAACTGATAATGGTTTAAGCGTTCTGAACGTTACCAAGATGGTTGACGTGCGACGTAATATGAACCGAATGATCAACTTCAATATGCCTGAGGATCTGACGGCTGGAAACTACAAAATAACTATCGATGGACAGCGTGGCTTCAGCTTTCACAAGGAGGCAGAGCTGGTGTATCTCAGCAAATCGATATCGGGGCTAATACAGGTCGATAAGCCCGTATTTAAACCTGGGGATACGGTGAACTTCCGTGTGATCGTGCTGGACACGGAGCTGAAACCGCCGGCGAGGGTCAAGTCGGTTTATGTAACTATACGAGATCCTCAGCGCAATGTGATTCGCAAATGGTCCACGGCAAAACTGTATGCCGGTGTGTTCGAGAGCGATCTACAGATAGCGCCTACTCCAATGCTCGGGGTCTGGAATATCTCGGTGGAGGTGGAAGGAGAAGAGCTTGTGTCAAAGACGTTTGAGGTGAAGGAGTACGTGTTGTCAACGTTCGACGTGCAGGTCATGCCATCGGTGATTCCACTGGAAGAGCATCAAGCTGTGAATCTTACAATCGAAGCGAACTATCACTTTGGTAAGCCAGTGCAAGGAGTGGCCAAGGTGGAGCTGTACCTAGACGACGATAAGCTAAAACTGAAAAAAGAGCTGACTGTGTACGGAAAGGGCCAGGTAGAGTTGCGCTTTGACAATTTTGCAATGGATGCGGATCAGCAGGATGTACCAGTGAAGGTGTCGTTCGTCGAGCAGTACACAAATCGTACGGTGGTCAAACAGTCACAAATCACGGTATATAGGTATGCGTACCGAGTAGAGTTGATAAAAGAGAGTCCACAGTTTCGTCCGGGACTCCCGTTCAAATGTGCGCTTCAGTTTACACACCATGATGGAACACCGGCTAAAGGCATTAGCGGTAAGGTAGAGGTATCCGATGTACGATTCGAAACGACAACAACGAGTGATAACGATGGATTGATTAAGCTCGAGCTGCAACCAAGTGAGGGTACTGAACAACTCAGTATTCACTTCAATGCTGTTGATGGATTCTTTTTTTATGAAGATGTGAATAAGGTAGAAACGGTTACGGATGCGTATATTAAACTGGAGCTGAAATCACCGCATCAAACGGAACAAATTGATGCGTTTCATGGTGACGTGCACGGAGCGCATGACATTCTTCGTGTACTATGTCATGTCAAAGGGCAATATCATCGATGCAGGATTCATGCGACCCAACAAGCAACCGAAGTACCTGTTGCAGCTGAACGCAACAGAAAAGATGATTCCGAGGGCGAAAATTCTCATCGCTACCGTAGCGGGCCGCACGGTGGTGTACGACTTCGCAGACCTCGATTTCCAAGAGCTTCGCAATAATTTTGATTTAAGCATTGACGAGCAAGAGATCAAGCCGGGACGACAAATCGAGCTGAGCATGTCTGGACGCCCAGGAGCGTACGTTGGGCTGGCCGCGTATGACAAAGCCTTGCTGCTTTTCAACAAGAACCACGACCTGTTCTGGGAGGACATTGGGCAGGTGTTTGATGGGTTCCATGCAATCAATGAGAACGAGTTTGACATATTCCACAGCTTGGGTCTGTTCGCCAGGACATTGGACGATATCTTGTTCGACAGTGCAAATGAAAAGACGGGGCGTAATGCACTGCAGTCAGGCAAGCCGATCGGCAAGCTGGTGTCGTATCGGACGAACTTCCAGGAATCGTGGTTGTGGAAAAATGTTTCCATCGGACGATCGGGAAGTCGCAAGTTGATCGAGGTAGTACCGGACACGACCACCTCCTGGTATCTGACGGGCTTCTCGATCGATCCCGTGTACGGGTTGGGTATCATCAAGAAGCCAATCCAGTTCACAACAGTCCAGCCGTTCTACATCGTAGAGAACTTACCATATTCAATCAAACGAGGCGAAGCGGTTGTGTTGCAGTTTACGCTGTTCAACAACCTTGGAGCGGAGTATATAGCCGATGTGACGCTGTACAATGTGGCCAACCAGACCGAGTTCGTCGGACGTCCAAATACGGATCTCAGCTACACCAAATCCGTGAGCGTTCCTCCAAAAGTTGGTGTGCCAATCTCGTTCCTCATCAAGGCCCGCAAGCTCGGCGAGATGGCGGTTCGTGTAAAGGCTTCGATAATGCTGGGACACGAAACGGACGCCCTGGAAAAGGTAATACGGGTGATGCCTGAAAGTTTGGTGCAGCCGAGAATGGATACACGCTTTTTCTGCTTCGACGATCACAAAAATCAAACGTTTCCGATCAACTTGGACATCAACAAGAAGGCCGACAGTGGATCGACAAAGATTGAGTTTCGACTAAATCCCAATTTGTTGACCACGGTCATCAAGAACCTGGACCATCTTCTCGGCGTTCCGACGGGATGTGGTGAGCAGAATATGGTCAAATTTGTTCCCAACATTTTGGTACTGGATTATTTGCATGCCATCGGGTCGAAAGAACAGCATCTAATCGACAAAGCTACGAATTTGTTGCGTCAAGGATATCAAAACCAGATGCGCTACCGTCAGACGGATGGTTCATTTGGTTTGTGGGAGACTACTAATGGTAGCGTGTTTCTCACCGCGTTCGTTGGCACATCGATGCAAACTGCAGTAAAATACATAAGCGATATTGATGCAGCAATGGTGGAGAAGGCATTGGATTGGTTAGCCTCGAAGCAGCATTTCTCGGGACGGTTTGACAAGGCCGGTGCAGAGTATCACAAAGAAATGCAAGGAGGGTTGCGCAATGGTGTGGCCCTCACATCATATGTGTTGATGGCATTGCTGGAGAATGACATTGCCAAAGCAAAGCACGCAGAGGTGATTCAAAAAGGAATGACCTATCTGAGCAATCAGTTTGGATCCATCAACAATGCATACGACCTATCGATAGCAACCTACGCGATGATGTTGAACGGACACACCATGAAGGAGGAGGCACTCAATAAGCTGATTGATATGTCTTTCATTGATGCTGATAAAAACGAACGGTTCTGGAACACAACGAATCCAATAGAAACCACCGCATATGCTCTGCTGTCGTTTGTGATGGCCGAGAAGTACACAGACGGTATACCGGTCATGAATTGGTTGGTGAATCAACGTTACGTTACCGGTAGCTTTCCGAGCACGCAAGACACGTTTGTGGGGCTGAAAGCGCTGACCAAAATGGCGGAAAAGATATCTCCGTCCCGAAACGACTACACCGTTCAACTGAAGTACAAGAAGAGTGCAAAATACTTCAAAATAAACTCGGAGCAAATTGATGTGGAAAACTTCGTGGATATACCGGAGGACACAAAAAAGCTCGAGATCAATGTGGGGGGCATTGGATTTGGGTTGTTAGAGGTGGTTTATCAATTTAATTTGAATCTCGTCAACTTTGAGAATAGATTCCAACTAGACCTGGAGAAACAGAACACAGGCTCTGACTACGAGCTGAGGCTGAAGGTCTGTGCCAGCTACATACCCCAGCTGACCGACAGACGATCGAACATGGCACTGATTGAGGTAACCTTACCGAGCGGTTACGTGGTTGATCGCAATCCGATCAGCGAGCAGACGAAGGTGAATCCGATTCAGAAAACTGAAATCCGTTACGGTGGCACTTCAGTCGTTTTATACTACGACAATATGGGCAGCGAGCGTAACTGTTTCACCCTGACCGCGTACAGACGCTTTAAGGTCGCATTGAAGCGTCCAGCGTATGTGGTTGTGTATGATTATTATAATACAAATCTGAACGCCATCAAAGTGTACGAAGTGGACAAGCAGAATTTGTGCGAAATCTGTGACGAAGAAGACTGTCCTGCAGAGTGCAAAAAATAG"
         e.g. ref_seq = s2_AAEL001802_RA = "ATGTCGGTATTCATACAAACGGACAAACCGGTGTATACCCCGGGAGATCTGATACGTTTTCGGGTAATCGTGGTGGATGCTGACACTAGACCTGTGACTAGTATTAAAACGGTAAATATAGCGATCGACGATTCTGCAAAAAATTCCATTCGAAAGTGGCCTTATGCCAAGTTGTTAAACGGCATCTTTGAGTCACAAGTGCAATTAGCTTCTTCGCCTGTTCTTGGCACCTGGATTATCAACGTAACAGCTTCCGACGACATCATTGTCACCAAACAGATAGAAGTTAAGGAATATGTGTTGCCAAAATTTTTCGTGAAAGTTTACCCTTCGGAGGTTCTATTGGGGAAAAATAAGAAGGTTTCTCTTACCTTAGATGCCTATTACACGTTCAAAGAACCCGTCGACGGCAATTACAAAGTTGAGTTATTTTTGGACCATACCAAGAGAAAGCCTGACTTCATAAAAAGTGATCGAATCACCGGTAAAACATCACTTGAGTTTCAATTGAAAAATGAAGTAGACATTGATGGCGACGAGCAGTACACTGATGTCACGGTTGAAGTTGAAGTTGTCGAGGCATTTTCTAATCGCACAGTTAGTATAACTGAGAATATTCCGATTTATCGTCAGCCTTATACCGTGACCCTTCTTCCATCTGCACCATCATTTCGACCAGGAGTTCCATTCAATGTACAAATAGTTGTGAAAGATCAGCTTGGACACCCTCCTGCCGAAGAAAAAGCGGCATCAATTGACCTTACTGTAGAGTTCCATTTGCCCATTGACAGTGACACCAAATCTATCACTGTAGATCTGGACGAGAAAGGAACAGGTCAGCTCACATTAGAGCCCCGCCCAGACGCCCAAGAACTGAAAGTGAACGCTACATATGACTCTCAACAATACGATGTAATTCACGATCCGATACATGGTTTCAGTTCGCAAAGTAAGCAGTACATCACAGTAACTCTGAATCCAAAATACTATAACAACATTAAAGTCGATAAGGACATCGTACTGGACATCTCCTGCACTGAAACAATGACGCACTTCTCGTACATCGTTGTCACCAGAGGAAACATAGTGGAAGCATCGAACGTTCCTGTCAGGATAAAAAAGAAACATTCTCTGAGATTGAAAATGACTTCAAAAATGTCTCCGGAGTCGAGGCTTCTAGTGTACTATACAAACAGGGAGTATCTCATCTTTGATGATATTGAGCTGAAGTTCGATTCGTTCAACAACGACTTCAAATTCGATTTGAACGATGATGAGTATTTTCCAGGGCAATCAGTTTATATCGATGTATACGCTTCAAAGGATTCATACGTTGCGTTCAGTGGAATCGATGAAAGTGTACTCCTGGTAGGCAAAGAGCGCCATGACTTCAACAAAGGAGATGTGCTCAAGGAACTCGCTCTTTACGGAGCAACAAATGATGCCGAGTTTGACTTGTTCCACGTAAGTTTCATGTCAAATGGTTTGATTATTCCAGTTAATGTATCTGTAACTCGCTCACAGAATGCACGATTTGGTACTCTACTAGGAAGGACTAGGCAGCAAGCGATTGAAATTCGAACTCAATTCCTAGAATCCTGGTTATGGAAATCCTTTTCCATGGATGGTCGAAACAACTTCAAAGCAATAGAAGACTCGGTTCCGGATACTATTACAACGTATCACGTGTCAGGATTTGCTTTAAGTCCAACACTAGGTCTTGGAGTAATCCAACAACCAGTGAGTTTCACCGTTCGTAAAAAATTCTACTTGGTTGCAAATTTGCCTTACTCGATCAAACGGGGTGAAGTGGCGTTGATTCAGGTTACCGTCTTCAACTTCCTAGGAAGCAGCATAACAACCGATGTGACGCTGTTCAATAAACGCGATGAAATTGAGTTTGTCGAGAATGCATCCACTAATAATACACATCGAACAAAGGCGGTAATTGTCCCGAATAACAATGGAAAATCTGTATCATTTATGGTGAAAGCAAAGAAATTAGGACAGATTGCGATCAAATTCCAGGCGGTAAACCTGCTGGAAACGGATGCATTGGAGCACATGTTACGAGTAACCCCAGAGAGCCATCGCTATGAGAAAAATGTAGCTCGATTCGTTGAGCTACCAAAGTTTGAGACGCAAACTTTCGATGTGAAGCTGGACATTCCCAAAAATATCGACGAGGGTTCTGCTCAAATCAAATTCACGTTAGACCCGGACATTTTGGGAACAGCCATCAGCAACCTAGACGGGTTGATCCGGAAACCCTTTGGATGTGGCGAACAAAATATGCTCCATTTTGTGCCAAATATAGTCGTTTTGGATTATCTTAACGAAACCAACACAGCGGCAGAAGATGTGAGGACCAAAGCGATAAATTTTCTTAGCAGCGGATATCAAAACCAGCTACGCTACAAACGTTCGGATGGGGCCTTCAGTGTCTGGGGACAATCGTATGCTGGCAGTACATTTTTGACGGCCTTTGTGGCGAAATCATTCAAAATAGCAGCCAAATACATTCAGGTGGATAAGTCTATAGTAGACGCGGCATTCGACTGGTTAGTGAAACAACAACAATCAGATGGGCGGTTCCCAGAAGTGGGGCAAGTATTCCAAGCAGATATGCAGGGTGGGCTTCGTAATAACGGTTTTGCGCTTACCGCGTATGTTCTGATCGCTTTTGCTGAAAATAAGGAAGTATACAGAAAATACCAATCACAACTGAACAAAACTACTAACTTCATAGCAGATAGACTTGCTAATATGGAGAATCCATACGACCTCTCGCTGTCCACTTATGCGTTGATGCTAACAAATCATGGCAAGCGCACCGAGTTTCTTCACAAATTAGTCGAAAAGTCGATATTTGACCGCAATCAAACTGAGAGATATTGGGACAGCAAACCAGTTGATATTGAAGTTGCTGGATATGCTCTATTGTCATACGTAGCTGCCGGTAAATTATTGGATGCAACGCCTATCATGCGGTGGCTCAACAAGCAGCGTTATGGTCTCGGAGGCTATCCTGGAACTCAGGAAACATTCGTTGGATTGAAAGCATTGGCAACGTTCGCTGCAAATGTAACTAGTAGGAGAAACGAATATACTGTAAGGATATTCTACGAACCAAATGGTCGACGAACATTCGACGTACACATGCACAATTCGTTTAATATTCAAGAGCTTGACATTCCTAATAACATCAGAAAAATGAAGGTGGAAGTTGAAGGCATCGGCAGAGGCTTCTTCCAAGTGGCATATCAGTACTATCAAAATATGCAGGTGGCTAAGCCCAGTTTCAGCATTACAATTAATCAGCTTAACACCACGACGGAACACATGCAGCAATTGGACGTGTGTGTGAAATACATACCAAAAGAGGCTTATCAAAAATCGAATATGGCTTTGGTGGAAATATTCTTGCCTAGTGGGCTTGTAGCAGACTCAGATGCCATTACGGACAAGACTGGAGGAATTCGAAGAATTGAAAGACGTTTTTCGGACACCTCAGTAGTTATATATTATGATAATTTGGACCCCGAAGACAAGTGCTTCCGAGTGACTGCTTATCGTCGGTATAAAATTGCATTGCATTTGCCATCATATATTATAGTTTATGATTATTATAATTTTGAGCGCTTTGCCATTCAAAAGTACGAAGGAAAGGTGCTGCAGCTCTGCGATATTTGTGAAGACGAGGACTGCGAAACTTTATCATGTCAAAATAGCTCGAAATTGGCAATAATGTAA"
      aln_gap_open, aln_gap_extend -   float/integer specifying how to score alignments, i.e. fewest gaps and smallest gaps get higher score

         e.g. aln_gap_open = -10
         e.g. aln_gap_extend = -0.5

   RETURNS:
      @todo
   """
   import numpy as np
   import pdb

   from Bio import pairwise2

   alns = pairwise2.align.globalxs(qry_seq, ref_seq, aln_gap_open, aln_gap_extend)  # @todo: make sure open and extend are in right order

   top_aln = alns[0] # sorte list of alternative alignments, highest-scoring alignent is alns[0]

   qry_seq_aligned = top_aln[0]
   ref_seq_aligned = top_aln[1]

   #print ref_seq_aligned
   #print qry_seq_aligned

   return qry_seq_aligned, ref_seq_aligned

def trim_gaps_from_aligned_seqs( qry_seq_aln, ref_seq_aln ):

   """ Trim the gaps away from aligned sequences generated by align_query_vs_reference.
   Aligned then gap-trimmed sequenes are thus in a format ready for input to dnds.py.
      e.g.  Aligned qry and ref seqs:
               ATGTGC----TAA  qry_seq
               ATG--A--TTTGA  ref_seq
                  ^^ ^^^^     gap positions****
            After trimming:
               ATGCTAA        qry_seq (after trimming)
               ATGATGA        ref_seq (after trimming)
   ARGS:
      qry_seq_aln, ref_seq_aln -    @todo
         e.g. qry_seq_aln = 'ATGTGC----TAA'
         e.g. ref_seq_aln = 'ATG--A--TTTGA'
   RETURNS:
      qry_seq_trimmed, ref_seq_trimmed -    @todo
         e.g. qry_seq_trimmed = 'ATGCTAA'
         e.g. ref_seq_trimmed = 'ATGATGA'
   """
   import numpy as np
   import pdb

   # Triming gaps and adjusting the qry indices
   qry_R=np.zeros(len(qry_seq_aln))

   j=0
   for i in range(len(qry_seq_aln)):
      if("-"!=qry_seq_aln[i]):
         qry_R[i]=j
         j=j+1

   qry_and_ref_arr = [(j,ref_seq_aln[i],qry_R[i]) for i,j in enumerate(qry_seq_aln) if (("-" != ref_seq_aln[i]) and ("-" != j))]

   qry_and_ref_arr = np.array(qry_and_ref_arr)

   qry_trimmed= "".join(list(qry_and_ref_arr[:,0]))
   ref_trimmed= "".join(list(qry_and_ref_arr[:,1]))
   indeces_pro= list(qry_and_ref_arr[:,2])
   indeces_pro1= indeces_pro[0::3]
   qry_indices= ",".join(indeces_pro1)

   #return qry_trimmed, ref_trimmed, qry_indices
   return qry_trimmed, ref_trimmed


def dnds( seq1, seq2, changes_potential, changes_observed, msCorrect='approximate', sliding=False, windowLength=3, stepLength=1):
    """ Perform dN/dS analysis, using the 'NG' algoritm, includes both whole sequence or sliding window, and either an approximate or exact multiple-substiution correction method. (@todo: make sure it actually is exact... it could be
             something else)

    ARGS:
        seq1,  a DNA sequence as string of letters, AGTC. Seq1 must be equal in length
            to, and aligned with, seq2, with gaps trimmed away. @todo: how on earth can
            we reliably make this work on the web service?

            e.g. seq1 = 'ATGCGCAAATACTCCCCCTTCCGAAATGGATACATGGAACCCACCCTTGGGCAGCACCTCCCAACCCTGTCTTTTCCAGACCCCGGACTCCGGCCCCAAAACCTGTACACCCTCTGGGGAGGCTCCGTTGTCTGCATGTACCTCTACCAGCTTTCCCCCCCCATCACCTGGCCCCTCCTGCCCCATGTGATTTTTTGCCACCCCGGCCAGCTCGGGGCCTTCCTCACCAATGTTCCCTACAAACGAATAGAAAAACTCCTCTATAAAATTTCCCTTACCACAGGGGCCCTAATAATTCTACCCGAGGACTGTTTGCCCACCACCCTTTTCCAGCCTGCTAGGGCACCCGTCACGCTGACAGCCTGGCAAAACGGCCTCCTTCCGTTCCACTCAACCCTCACCACTCCAGGCCTTATTTGGACATTTACCGATGGCACGCCTATGATTTCCGGGCCCTGCCCTAAAGATGGCCAGCCATCTTTAGTACTACAGTCCTCCTCCTTTATATTTCACAAATTTCAAACCAAGGCCTACCACCCCTCATTTCTACTCTCACACGGCCTCATACAGTACTCTTCCTTTCATAATTTGCATCTCCTATTTGAAGAATACACCAACATCCCCATTTCTCTACTTTTTAACGAAAAAGAGGCAGATGACAATGACCATGAGCCCCAAATATCCCCCGGGGGCTTAGAGCCTCTCAGTGAAAAACATTTCCGTGAAACAGAAGTC'

        seq2,  a DNA sequence similar to seq1 but with differences (substitutions),
            representing a CDS orthologue of seq1 from a different species. Read
            description of seq1 for other required similarities to avoid errors.

            e.g. seq2 = 'ATGCGCAAGTACTCCCCCTTCCGAAACGGATACATGGAACCCACCCTTGGGCAACACCTCCCAACCCTGTCTTTTCCAGACCCCGGCCTCCGGCCCCAAAACCTGTACACCCTCTGGGGAGACTCTGTTGTCTGCCTGTACCTCTACCAGCTCTCCCCCCCCATCACCTGGCCCCTCCCGCCCCATGTGATTTTTTGCCACCCCGGCCAGCTCGGGGCCTTCCTCACCAATGTTCCCTACAAGCGTATGGAAGAACTCCTCTATAAAATTTCCCTTACCACAGGGGCCCTAATAATTCTACCCGAGGACTGTTTACCAACCACCCTTTTCCAGCCTGCTAGGGCCCCCGTCACGTTGACCGCCTGGCAGAACGGCCTCCTTCCGTTCCACTCAACCCTCACCACTCCAGGCCTTATTTGGACATTTACCGATGGCACGCCTATGGTTTCCGGACCCTGCCCCAAAGATGGCCAGCCATCTTTAGTACTACAGTCCTCCTCATTTATATTTCACAAATTTCAAACCAAGGCCTACCACCCTTCATTTCTACTCTCACACGGCCTCATACAGTACTCCTCCTTTCACAATTTACATCTCCTTTTTGAAGAATACACCAACATCCCCGTTTCTCTACTTTTTAACGAAAAAGAGGCAAATGACACTGACCATGAGCCCCAAATATCCCCCGGGGGCTTAGAGCCTCCCGCTGAAAAACATTTCCGCGAAACAGAAGTC'
        changes_potential, a dict, with key=pair of codons tuple, e.g. ('ATG','ATG'), and value=['S':<S>,'N':<N>]. Where <S> is the number of potential synonmyous sites for each codon (averaged between the two codons), and <N> is the same but for non-synonymous sites.
            e.g. changes.potential_changes_dict(...)  (see: ./changes.py)
        changes_observed, @todo
        msCorrect, a string to toggle between multiple-substitution correction methods:
            "approximate", "exact" (@todo: make sure it actually is exact... it could be
             something else)

            e.g. msCorrect = 'approximate'
        sliding, a boolean to toggle between sliding window analysis (vector of dN/dS values at successive chunks of sequence) or whole sequence analysis (a single
            dN/dS value for the given pair of input sequences), either: True, False

            e.g. sliding = False
        windowLength, an integer specifying the width of the sliding window, measured in no. of codons in the window to measure dN/dS over, from 1-to-length(seq1)

            e.g. windowLength = 50
        stepLength, an integer specifying no. of codons to shift the sliding window with each iteration. If stepLength < windowLength then windows will overlap, overlapping is dealt with prior to plotting (acts to smooth values, averages along the overlaps are taken as a dN/dS value for any codon).
            e.g. stepLength = 1
    NOTES:
        Sources of formulae:
            http://www.megasoftware.net/mega4/WebHelp/part_iv___evolutionary_analysis/computing_evolutionary_distances/distance_models/synonymouse_and_nonsynonymous_substitution_models/hc_nei_gojobori_method.html
    """
    import sys
    import pickle
    import math
    import warnings
    import numpy as np
    import pdb
    import time

    def chunks(l, n):
        """ Yield successive n-sized chunks from l. """
        for i in xrange(0, len(l), n):
            yield l[i:i+n]

    # todo: stop codons to deal with, reject
    # todo: ambiguous bases to deal with:
        # gaps,  @done
        # Ns,    @todo
        # Xs     @todo

    warning_count = 0

    # STATS per CODON-PAIR:
    codons_seq1   = [codon for codon in chunks(str(seq1),3)]  #splits
    codons_seq2   = [codon for codon in chunks(str(seq2),3)]
    #print codons_seq1
    #print codons_seq2
    codons_paired = [pair for pair in zip(codons_seq1,codons_seq2) if (len(pair[0])+len(pair[1]))== 6 and 'N' not in pair[0]+pair[1] and '-' not in pair[0]+pair[1]] # aligned codons are paired into tuples, excess codons are truncated, @todo: in main example, we lose 5 bps of data
    # @done: the next for loop is extremely innefficient, I should set the structure of the changes_potential and changes_observed dicts to how I want it to look, a priori, i.e. when it's instantiated in changes.py.
    # OR just remove this chunk and access the data as they come in the args
    #print codons_paired
    changes_all = {'observed':{'S':[],'N':[]},'potential':{'S':[],'N':[]}}

    for pair in codons_paired:
        changes_all['potential']['S'].append(changes_potential[pair]['S'])
        changes_all['potential']['N'].append(changes_potential[pair]['N'])
        changes_all['observed']['S'].append(changes_observed[pair]['S'])
        changes_all['observed']['N'].append(changes_observed[pair]['N'])

    list_S  = changes_all['potential']['S']
    list_Sd = changes_all['observed']['S']

    list_N  = changes_all['potential']['N']
    list_Nd = changes_all['observed']['N']

    if sliding:
        # STATS for each WINDOW seq
        intervals    = range(0,len(codons_paired)-windowLength+1,stepLength)
        windows      = zip(intervals,[i + windowLength - 1 for i in intervals])

        window_stats = {}

        #window_stats_list = []

        # @done: test against matlab's sliding window, also @todo: find out what stepLength does, @todo: try to plot the sliding window version

        for window_i,window in enumerate(windows):

            start = window[0]
            end   = window[1]+1

            window_stats[window] = {    'S':sum(list_S[start:end]),
                                        'Sd':sum(list_Sd[start:end]),
                                        'N': sum(list_N[start:end]),
                                        'Nd':sum(list_Nd[start:end])    }


            pS = window_stats[window]['Sd']/window_stats[window]['S']
            pN = window_stats[window]['Nd']/window_stats[window]['N']

            try:
                if msCorrect=='approximate':
                    dN = -(3./4.)*math.log(1.-(4./3.)*pN)
                    dS = -(3./4.)*math.log(1.-(4./3.)*pS)

                # @todo: what is this commented code? I don't remember...
                # elif msCorrect=='exact':
                #     d=ln(1-p*4/3)/ln(1-3/(4*N))

                else: # msCorrect=='????'  # @todo: is this the exact one? Or something else?
                    dN = pN
                    dS = pS
                window_stats[window]['dNdS'] = dN/dS
            # @todo: I'm not sure I'm treating the following exceptions in the right way...
            # technically it woud be best to exclude these from downstream analyses?
            # e.g. missing value/datapoint on a plot of dN/dS (y-axis) vs. window interval (x-axis)
            except ZeroDivisionError:
                warning_count += 1
                #warn_msg = "Query and Reference sequences are too divergent. Approximate multiple-substitutions correction cannot be achieved, for window: "+str(window_i)+", dS is zero, leading to a division error when trying dN/dS... try alternative value for argument: msCorrect (e.g. 'exact') OR alternative value for argument: windowLength (e.g. "+str(windowLength+20)+") ...\n"   # # @TODO: uncomment for verbose warning message prints // @ANDY-2017-01-30
                #warnings.warn(warn_msg)  # @TODO: uncomment for verbose warning message prints // @ANDY-2017-01-30
                window_stats[window]['dNdS'] = float('Inf')
            except ValueError:
                warning_count += 1
                #warn_msg="Query and Reference sequences are too divergent. Approximate multiple-substitutions correction cannot be achieved, for window: "+str(window_i)+",  SYNONYMOUS changes per synonymous site, pS>=3/4, log() operation will yeild return undefined... try alternative value for argument: msCorrect (e.g. 'exact') OR alternative value for argument: windowLength (e.g. "+str(windowLength+20)+") ...\n"  # # @TODO: uncomment for verbose warning message prints // @ANDY-2017-01-30
                #warnings.warn(warn_msg)  # @TODO: uncomment for verbose warning message prints // @ANDY-2017-01-30
                window_stats[window]['dNdS'] = float('nan')

        return window_stats,warning_count  # list of dnds per window interval // dict of dnds, key=(<from #base pair>,<to #base pair>), value=<dN/dS of the window specified in the key>
    else:
        # STATS for WHOLE SEQ
        S   = sum(list_S)
    #    print "S:", S
        Sd  = sum(list_Sd)
    #    print "Sd:", Sd
        pS  = Sd/S
        N   = sum(list_N)
    #    print "N:", N
        Nd  = sum(list_Nd)
    #    print "Nd:", Nd
        pN  = Nd/N

        try:
            if msCorrect=='approximate':

                if (pS>=3./4.):
                    raise ValueError("Query and reference sequences are too divergent. Approximate multiple-substitutions correction cannot be achieved, SYNONYMOUS changes per synonymous site, pS>=3/4, log() operation will yeild return undefined. Try alternative value for argument: msCorrect (e.g. 'exact')...")

                if (pN>=3./4.):
                    raise ValueError("Query and reference sequences are too divergent. Approximate multiple-substitutions correction cannot be achieved, NON-SYNONYMOUS changes per synonymous site, pN>=3/4, log() operation will yeild return undefined. Try alternative value for argument: msCorrect (e.g. 'exact')...")

                dS  = -(3./4.)*math.log(1.-((4./3.)*pS))
                dN  = -(3./4.)*math.log(1.-((4./3.)*pN))
                dN_dS = dN/dS

            else: # @todo: is this the exact one? Or something else?

                # @DONE: one day the following three lines of code will error, giving a ZeroDivisionError, this needs to be handled with try
                dS = pS  # i.e. dS = Sd/S
                dN = pN
                dN_dS = dN/dS
        except ValueError:
            warning_count += 1
            warnings.warn("Query and reference sequencea are too divergent. ValueError: Approximate multiple-substitutions correction cannot be achieved: UNKNOWN reason, probably due to illegal numbers in a log() function...\n")
            dN_dS = float("nan")
        except ZeroDivisionError:
            warning_count += 1
            warnings.warn("Query and reference sequences are too divergent. ZeroDiviSionError: Approximate multiple-substitutions correction cannot be achieved: UNKNOWN reason, probably due to illegal numbers in a log() function...\n")
            dN_dS = float('Inf')

        return dN, dS, dN_dS, warning_count  # i.e. omega = dN/dS = (Nd/N)/(Sd/S)


def plot_dnds_sliding(dnds_slide_dict):

    """ Plots sliding dN/dS values (y-axis) along the input sequence's aligned codons (x-axis). If the sliding windows overlap (see below), then plot_dnds_sliding will also average the overlapping dN/dS values for each codon.
        ----     window 1
         ----    window 2
          ----   ...
        ^^^^^^^  take average along each column (codon)
    ARGS:
        dnds_slide_dict,    the output of dnds() if the 'sliding' optional argument is set to True
            e.g. dnds_slide_dict = dnds( s1, s2, potential_changes, observed_changes, msCorrect='approximate', sliding=True, windowLength=50, stepLength=1 )
    RETURNS:
        None,   a plot is generated and written to: py/data/dnds_sliding_test.png
    """
    import sys
    import pickle
    import math
    import warnings
    import numpy as np
    import changes as codon_pair_data
    import align as align_then_trim
    import pdb
    import time

    window_intervals = dnds_slide_dict.keys() # @todo: I dont think sorting the windows will make a difference to final result, but it will make it slower, @todo: test this just in case
    max_window_position = np.amax(window_intervals) # e.g. 243 @done: a whole order of magnitude faster than: max_window_interval = max(window_intervals, key=lambda x: x[1])
    overlap_matrix      = np.empty((len(window_intervals),max_window_position+1))  # @todo: are you sure it's +1? initialize empty matrix, note: entries are not actually NaN, just near-zero
    overlap_matrix[:]   = np.NAN # initiate empty np array with NaN, so later we can mask

    for window_i,window in enumerate(window_intervals):

        start = window[0] # e.g. 0
        end   = window[1] # e.g. 49

        # in the i-th row, fill all elements from the window[0]-th to window[1]-th with the dN/dS value for this window
        overlap_matrix[window_i,start:end+1] = dnds_slide_dict[window]['dNdS'] # @todo: are you sure it's +1? test, keep in mind for these indices it does -1 for the "to" part

    nan_masker              = ~np.isfinite(overlap_matrix) # boolean matrix, True if element is finite, False if element is Inf or NaN
    overlap_matrix_masked   = np.ma.masked_array(overlap_matrix,mask=nan_masker)
    overlap_matrix_avg      = overlap_matrix_masked.mean(axis=0)

    return list(overlap_matrix_avg), overlap_matrix_avg.mean()



class DNAmsaDist:

    def __init__(self, file, genename):

        self.file = file
        self.genename = genename


    def MeanKimuraDist2Cons(self):
        alignment = AlignIO.read(open(self.file), "fasta")
        summary_align = AlignInfo.SummaryInfo(alignment)
        cons_seq = summary_align.dumb_consensus()
        K2Plist=[]
        for record in alignment:
            pct_missing = 100*(record.seq.count('-')+record.seq.count('N'))/len(record.seq)
            if record.id == "NC_007605" or pct_missing > 20:
                pass
            else:
                dist = GeneticDist(cons_seq, record.seq)
                K2Plist.append(dist.K2Pdistance())

        #gene name, mean distance, standard error of mean
        return self.genename, numpy.mean(K2Plist), stats.sem(K2Plist)

    def MeanPairwiseKimuraDist(self):
        alignment = AlignIO.read(open(self.file), "fasta")
        K2Plist=[]
        for i in range(0,len(alignment)):
            for j in range(i+1,len(alignment)):
                p1_pct_missing = 100*(alignment[i].seq.count('-')+alignment[i].seq.count('N'))/len(alignment[i].seq)
                p2_pct_missing = 100*(alignment[j].seq.count('-')+alignment[j].seq.count('N'))/len(alignment[j].seq)

                if alignment[i].id == "NC_007605" or alignment[j].id == "NC_007605" or p1_pct_missing > 20 or p2_pct_missing > 20:
                    pass
                else:
                    dist = GeneticDist(alignment[i].seq, alignment[j].seq)
                    K2Plist.append(dist.K2Pdistance())

        output=str(self.genename)+"\t"+str(numpy.mean(K2Plist))+"\t"+str(stats.sem(K2Plist))+"\t"+str(len(K2Plist))
        return output

    def MeandNdS(self):
        """
        functions are from https://github.com/a1ultima/hpcleap_dnds/
        """

        alignment = AlignIO.read(open(self.file), "fasta")
        dNlist=[]
        dSlist=[]
        dNdSlist=[]
        nt_to_aa_dict = geneticCode("standard")

        observed_changes  = observed_changes_dict(nt_to_aa_dict)
        potential_changes = potential_changes_dict(nt_to_aa_dict)


        for i in range(0,len(alignment)):
            for j in range(i+1,len(alignment)):
                p1_pct_missing = 100*(alignment[i].seq.count('-')+alignment[i].seq.count('N'))/len(alignment[i].seq)
                p2_pct_missing = 100*(alignment[j].seq.count('-')+alignment[j].seq.count('N'))/len(alignment[j].seq)

                if alignment[i].id == "NC_007605" or alignment[j].id == "NC_007605" or p1_pct_missing > 20 or p2_pct_missing > 20:
                    pass
                else:
                    dN, dS, dNdS, warning_count = dnds( alignment[i].seq, alignment[j].seq, potential_changes, observed_changes, msCorrect='approximate', sliding=False, windowLength=50, stepLength=1 )
                    dNlist.append(dN)
                    dSlist.append(dS)
                    dNdSlist.append(dNdS)

        dNlist = numpy.array(dNlist)
        dNlist = dNlist[dNlist < 1E308]
        dSlist = numpy.array(dSlist)
        dSlist = dSlist[dSlist < 1E308]
        dNdSlist = numpy.array(dNdSlist)
        dNdSlist = dNdSlist[dNdSlist < 1E308]

        output=str(self.genename)+"\t"+str(numpy.mean(dNlist))+"\t"+ str(stats.sem(dNlist))+"\t"+str(numpy.mean(dSlist))+"\t"+str(stats.sem(dSlist))+"\t"+str(numpy.mean(dNdSlist))+"\t"+str(stats.sem(dNdSlist))+"\t"+str(len(dNdS))

        return output
