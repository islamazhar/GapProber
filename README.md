# GapProber

GapProber is a software to fill up gaps in genome assemblies. 

## Background
The task of determining the entire genome sequence has diverse applications ranging from performing experiments to know about an organism to understanding how individuals are related to each other. However because the genome length is very long, it is not possible for any sequencer to sequence the whole genome at once using current sequencing machines. When the sequencer stitches the small reads from genome sequences together, there exist some gaps where the nucleotide sequence is unknown.

## Results

To solve this problem, we have present a probabilistic method for filling the gaps, which unlike other methods takes into account the gap length and formulates the problem according to an iterative algorithm similar to expectation-maximization (EM) algorithm. Our experiment on simulated reads from [wgsim](https://github.com/lh3/wgsim) without errors shows that this novel approach can fill up gaps of length up-to `500` perfectly with quick convergence. Our experiments on [GAGE](http://gage.cbcb.umd.edu/data/index.html) dataset reveals the GapProber extensively outperfroms state of the art gap filling tools `Gap2Seq`, `GapFiller`, `GapCloser` on  `Rhodobacter sphaeroides` genome assemblies while keeping competitive performance on genome assemblies of other datasets including `Human Chromosome 14:` and `Staphylococcus aureus`.   
<!--
We are working on filling gaps using reads with errors from [GAGE]({http://gage.cbcb.umd.edu/data/index.html}) dataset and eventually compare our result with state of the art gap fi. lling tools `Gap2Seq`, `GapFiller`, `GapCloser`,  and `MindTheGap`.
-->
## Running GapProber

GapProber depends on 

- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.html).
- python 3.7
- C++ 

GapProber needs to map the read paired to the scaffolds. Essentially any mapper which reports multiple alignment reads and allows insertions/deletions can be used. We used  `Bowtie2` which is a popular mapping tool in research community to do this.

`Bowtie2` is available on deafult package manager for most linux distributions. For example, if your Ubuntu linux machine do not have bowtie2, you can install it simply by running for  

`sudo apt-get install bowtie2`

You can  get `Bowtie2`  via `Bioconda` as well (https://anaconda.org/bioconda/bowtie2). But for that you will need to install `Bioconda`. 
 
You also need to have `python3.7` and `C++`. To install them run

`sudo apt-get install build-essential` 

Once you have `Bowtie2`, `Python 3.7`, and `C++` install,  to execute GapProber cd to the cloned repository and run the following command. 

 `Python GapProber.py --output <path to the output filled contig> --scaffolds <path to the scaffold file> --reads <fragment paried reads seperated by commas>`
 
Running `Python GapProber.py --help` will show you a list of available options.

 ## Notice
This is an ongoing work. We are always looking ways to improve our code.

If you want to use any part of the code / have any suggestion / report bugs please shoot an email at [mazhar.buet11@gmail.com](mailto:mazharbuet11@gmail.com). 

   

