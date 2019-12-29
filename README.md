# GapProber

GapProber is a software to fill up gaps in scaffolds. 

## Background
The task of determining the entire genome sequence has diverse applications ranging from performing experiments to know about an organism to understanding how individuals are related to each other. However because the genome length is very long, it is not possible for any sequencer to sequence the whole genome at once using current sequencing machines. When the sequencer stitches the small reads from genome sequences together, there exist some gaps where the nucleotide sequence is unknown.

## Results

To solve this problem, we have present a probabilistic method for filling the gaps, which unlike other methods takes into account the gap length and formulates the problem according to an iterative algorithm similar to expectation-maximization (EM) algorithm. Our experiment on simulated reads without errors shows that this novel approach can fill up gaps of length up-to 500 perfectly with quick convergence. We are working on filling gaps using reads with errors from [GAGE]({http://gage.cbcb.umd.edu/data/index.html}) dataset and eventually compare our result with state of the art gap filling tools `Gap2Seq`, `GapFiller`, `GapCloser`,  and `MindTheGap`.
