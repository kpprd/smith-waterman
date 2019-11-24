# smith-waterman
Python version: Python 3.7
Testrun: Macbook Pro OS 10.15.1. 


This program uses the  Smith-Waterman algorithm to compute the optimal local alignment(s) of two sequences. It prompts the
user to input file name/paths for sequence files and scoring matrix, and, optionally, gap penalties and whether the user wants
to find all optimal local alignments or just one. Then it reads the sequence files and matrix file and stores the sequence as
instances of the class Sequence, and the scoring matrix as a nested scoring dictionary. The scoring dictionary is formatted so
that scoring_dict[“A”][“B”] gives the score for a (mis)match between A and B. The function compute_alignments then does two
things: firstly it computes two matrices, one containing alignments scores and one containing information indicating insertion
(and how large an insertion gap), deletion (and how large a deletion gap) or match-mismatch. Seondly, it computes the alignment(s) based on the information stored in the two matrices. All alignments are
stored as instances of the class Alignment while they are being computed. This class stores the sequence strings so far
computed (including gaps) and the current matrix indices and names of the sequences etc. The program loops over all the stored
alignments (if there are more than one and the program is set to find all) and updates the sequence strings and indices for
each round (or, more precisely, the program initializes new Alignment instances that “inherit” the information from previous
alignments.) When the program encounters an ‘s’ in the direction matrix, indicating that the beginning of the alignment is
found, the alignment is finalized in a presentable format and stored in the found_alignments list, that is ultimately
returned. The alignment(s) is/are then printed, and optionally saved to a file.
