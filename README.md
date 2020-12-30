# Band-Align

Global alignment of two sequences obtained through a variation of Smith-Waterman's algorithm: this algorithm uses the concept of band and BLOSUM62 matrix in order to calculate the score of the alignment.  

## Input (from standard input):

- **l1**: Length of the first string to align
- **s1**: First sequence to align
- **l2**: Length of the second string to align
- **s2**: Second sequence to align

## Output (to standard output):

- **score**: Score obtained by the global alignment of the two sequences
- **s1_align**: Sequence obtained from *s1* adding indels in order to get the maximum score for the alignment 
- **s2_align**: Sequence obtained from *s2* adding indels in order to get the maximum score for the alignment
- **res_align**: for each position of the alignment, it denotes whether there is an indel or whether the two characters in the sequences are equal or not
  - **I** = Indel
  - **E** = Equal
  - **U** = Unequal
  
## Sample:

In 'samples/' there are some examples of possible inputs and outputs: the 'input[i].txt' contains the i-th input sample and 'output[i].txt' contains the i-th result of the execution using the corresponding i-th input.
