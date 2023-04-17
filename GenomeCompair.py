from Bio import SeqIO
from Bio import pairwise2
from Bio.Seq import Seq

# read in the two genome sequences from file
seq1 = SeqIO.read("genome1.fasta", "fasta")
seq2 = SeqIO.read("genome2.fasta", "fasta")

# calculate sequence identity
identity = seq1.seq == seq2.seq
percent_identity = identity.count(True) / len(identity) * 100

print(f"The two genomes have {percent_identity:.2f}% sequence identity using a simple comparison.")

# perform pairwise sequence alignment
alignments = pairwise2.align.globalxx(seq1.seq, seq2.seq)
best_alignment = alignments[0]
aligned_seq1 = best_alignment.seqA
aligned_seq2 = best_alignment.seqB
alignment_score = best_alignment.score
percent_alignment_identity = identity.count(True) / len(identity) * 100

print(f"The two genomes have {percent_alignment_identity:.2f}% identity using pairwise sequence alignment.")

# translate the sequences and compare the resulting proteins
proteins1 = seq1.seq.translate()
proteins2 = seq2.seq.translate()

identity = proteins1 == proteins2
percent_identity = identity.count(True) / len(identity) * 100

print(f"The two genomes have {percent_identity:.2f}% protein sequence identity.")

# reverse complement one of the sequences and compare again
rc_seq2 = seq2.seq.reverse_complement()

identity = seq1.seq == rc_seq2
percent_identity = identity.count(True) / len(identity) * 100

print(f"The two genomes have {percent_identity:.2f}% sequence identity after reverse complementing one sequence.")
