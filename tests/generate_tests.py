"""Generate test cases for merge-mates
"""

from random import randint

# excerpt HIV HXB2 reference sequence
seq = ''.join("""
GACTTACTCTTGATTGTAACGAGGATTGTGGAACTTCTGGGACGCAGGGGGTGGGAAGCCCTCAAATATT
GGTGGAATCTCCTACAGTATTGGAGTCAGGAACTAAAGAATAGTGCTGTTAGCTTGCTCAATGCCACAGC
CATAGCAGTAGCTGAGGGGACAGATAGGGTTATAGAAGTAGTACAAGGAGCTTGTAGAGCTATTCGCCAC
ATACCTAGAAGAATAAGACAGGGCTTGGAAAGGATTTTGCTATAAGATGGGTGGCAAGTGGTCAAAAAGT
AGTGTGATTGGATGGCCTACTGTAAGGGAAAGAATGAGACGAGCTGAGCCAGCAGCAGATAGGGTGGGAG
CAGCATCTCGAGACCTGGAAAAACATGGAGCAATCACAAGTAGCAATACAGCAGCTACCAATGCTGCTTG
TGCCTGGCTAGAAGCACAAGAGGAGGAGGAGGTGGGTTTTCCAGTCACACCTCAGGTACCTTTAAGACCA
ATGACTTACAAGGCAGCTGTAGATCTTAGCCACTTTTTAAAAGAAAAGGGGGGACTGGAAGGGCTAATTC
ACTCCCAAAGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAAGGCTACTTCCCTGATTAGC
""".split())
read_length = 251

def fill(l):
    """generate a random DNA sequence of length l"""
    s = ''
    for _ in range(0, l):
        s += "ACGT"[randint(0,3)]
    return s

def rc(seq):
    """reverse complement sequence"""
    comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join(map(lambda x: comp[x], reversed(seq)))

def read(seq, serial, r=1, data=0):
    """print fastq record"""
    seqid = "@M00000-{} {}:0:{}".format(serial, r, data)
    qual = 'N' * len(seq)
    return "{}\n{}\n+\n{}".format(seqid, seq, qual)

def write_reads(frag):
    pass

if __name__ == "__main__":
    forward = open('r1.fq', 'w')
    reverse = open('r2.fq', 'w')

    for i in range(0, read_length * 2):
        fragment = seq[100:100 + i]

        r1 = None
        r2_rc = None

        if i < read_length:
            r1 = fragment + fill(read_length - i)
            r2_rc = fragment + fill(read_length - i)
        else:
            r1 = fragment[0:read_length] 
            r2_rc = fragment[-read_length:]
 
        print(read(r1, i, r=1, data=i), file=forward)
        print(read(rc(r2_rc), i, r=2, data=i), file=reverse)
