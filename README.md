[README.txt](https://github.com/user-attachments/files/23760205/README.txt)
k-mer primers
C++
Curtis Wood
11.01.25

Upon compiling and running k-mer primers, the user is prompted for a filename containing a genome sequence.  Any FASTA sequence is acceptable,
so long as it is a text file without header info.  Only sequences of standard nucleotides (G,A,T,C) are supported, no combo nucleotides (W,S,etc).
Program has been validated with genome sequences up to 6 Mbp.

User provides an integer K of base pairs defining the k-mer length that start and end in a number P of nucleotides that will be used as primers
with melting temperature (Tm) between TL, the lower bound, and TU, the upper bound.  Temperature bounds are to be provided in celsius.  With the parameters
given, the program identifies unique k-mers of length K in the provided sequence that also have forward/reverse primer sequences of length P with
a calculated Tm falling between the desired temperature bounds.

For example, for K = 100, P = 20, TL = 58, TU = 60,
the provided genome will be searched for 100 bp sequences that occur only once in the entire genome that also return a Tm between TL and TU
calculated for nucleotides 1-20 and nucleotides 81-100 of the 100 bp sequence.  The standard Tm calculation is used, being:

Tm = 64.9 + ((41.0 * (g + c -16.4)) / (g + c +a +t))
g = # of guanine in primer
c = # of cytosine in primer
a = ## of adenine in primer
t = # of thymine in primer

Identified k-mers matching the criteria are output to the console and written to a file named "supplied genome filename".kmer with the following 
format for each k-mer found:

Tm_f : Start position - f_primer sequence ... r_primer sequence - Final position : Tm_r

Example (not calculated, just to illustrate formatting of output):
58.9 : 16789 - AATGCGTACCCGTAA ... AATGCGTACCCGTAA - 16989 : 60.1

The positions identified are simply counted from the first bp in the sequence to be searched.  The actual k-mer is not output as that would become
unreadable for long sequences.  Using the positions, one can easily find the identified k-mers by opening the sequence in any standard genome editor.
The program's utility lies primarily in quickly identifying valid targets in bacterial genomes that can act as single copy amplicons for qPCR.
