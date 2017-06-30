# fasta2ABC  Produces input files for ABC 
 
## Compilation:  
g++ fasta2ABC.cpp -std=c++17 -O3 -o fasta2ABC
 
## Exemple: 
./fasta2ABC locA_ama.fas locA_txn.fas 0.25 10 locusA  
**arg1** name of the fasta file for species A (string).  
**arg2** name of the fasta file for species B (string).  
**arg3** maximum proportion of N allowed to keep an individual (float).  
**arg4** minimum number of individuals to keep a species (integer).  
**arg5** gene's name used for output files names (string).  
  
