# GapFiller

Since the genome length of an organism is very long, it is not possible for any current sequencer machine to sequence the whole genome at once. 
When the sequencer stitches the small reads from genome sequences together, there exist some gaps where the nucleotide sequence is unknown. 
To solve this problem, we have designed a probabilistic method for filling the gaps, which unlike other methods takes into account the gap length. 
Our experiment on simulated data withour sequencing error shows that this novel approach can fill up gaps of length with `< 400` with quick convergence alongside significant accuracy. 
We are working on real data sets to handle sequencing errors