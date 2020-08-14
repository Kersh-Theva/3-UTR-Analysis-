# 3' UTR Analysis using CATS-seq data
Developed a program to analyze 3'UTRs using raw CATS-seq data and Integrated Genome Viewer (IGV).   

Main results: <br>
<img align ="center" src="https://github.com/Kersh-Theva/3-UTR-Analysis-/blob/master/TermLengthHistogram.png" width="600" height="400"> <br>
**Figure 1. 3'UTR terminator read count vs. length for RNA from *Methanococcus maripaludis***. 3'UTRs were obtained from analysis of processed CATS-seq data. Terminators were defined as contiguous 3'UTR reads. Primary, secondary and tertiary terminators were defined by the position of the 3' boundary of the terminator with respect to the end of the gene, where primary terminators ended most closely to the end of the gene and tertiary terminators were the furthest. 

<img src="https://github.com/Kersh-Theva/3-UTR-Analysis-/blob/master/IntraGene.png" width="950" height="400"> <br>
**Figure 2. Count of intragenic and intergenic terminators for primary, secondary and tertiary terminators from *Methanococcus maripaludis*.** Count of intragenic and intergenic terminators suggests that 30% or less of terminators were found within the downstream gene for primary, secondary and tertiary terminators (Term. Number = 1.0, 2.0  or 3.0 respectively). This resembles the observed overlap of 3' terminators and downstream genes in other methanogenic species. 

## Code and Resources Used
Python Version: 3.8.3 <br>
Packages: pandas, numpy, seaborn, matplotlib <br>

Bowtie Read Aligner Version: 1.2.3 (For more information about Bowtie, check out this [link](http://bowtie-bio.sourceforge.net/index.shtml) <br>
Samtools Version: 1.7-1 (For more information about Samtools, check out this [link](http://samtools.sourceforge.net/)<br>
Integrated Genome Viewer Version: 2.8.0 (For more information about IGV, check out this [link](http://software.broadinstitute.org/software/igv/)) <br>

Relevant publications:
[Brenneis *et al.*.(2007). PLoS Genetics **3**: 229](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2151090/) <br>
[Dar *et al.*. (2016). Nature Microbiology **1**:1-9](https://www.nature.com/articles/nmicrobiol2016143#Sec15). <br>
[Turchinovich *et al.*. (2014). RNA Biology **11**:817-828](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4179956/) <br>
Thank you to UC Davis Genome Center for preparing CATS-seq libraries and performing NGS using Illumina HiSeq. <br>

## What is CATS-Seq? 
Capture and amplification by tailing and switching (CATS)-seq is a RNA library preparation method that involves adding polyadenylated tails to the 3' end of RNA in order to then reverse transcribe the RNA with  primers complementary to the poly-A tail. Using this approach, one can obtain sequence data for the 3' untranslated regions (3'-UTR) for mRNA molecules from bacterial or archaeal species that do not typically have polyadenylated messenger RNA. 

![CATS-Seq](https://github.com/Kersh-Theva/3-UTR-Analysis-/blob/master/CATS-Seq%20illustration.png) <br>
**Figure 3. CATS-Seq Workflow** Process of RNA library preparation for the CATS-seq protocol. After polynucleotide tailing of the RNA with either poly-A polymerase or terminal deoxytransferase, followed by the use of a template-switching Moloney murine leukemia virus reverse transcriptase to extend the 3' end of a complementary poly-T containing primer and generate the first cDNA strand for downstream NGS library preparation.  (credit to Diagenode for image) 

## Data Processing Pipeline

**Convert Raw FASTA -> SAM file**: Use Bowtie read mapper and the command line tool to convert raw FASTA/FASTQ CATS-Seq files to reads aligned to your genome of interest. Remember to build your own index if required. Sample code is below: 
~~~
bowtie-build reference.fna output_file #build indexed reference genome
bowtie -S your_genome_of_interest inputreads.fastq  output.sam #Map reads and store in SAM file 
~~~

**Convert SAM file to BAM file for viewing with IGV**: In general, use Samtools and the command line tool to convert Bowtie SAM files into BAM files for use with IGV. Sample code is below: 
~~~
samtools view -bS -o output.bam input.sam
samtools sort input.bam -o input.sorted.bam
~~~
For this project, strand specific data was required to get 3'UTR data for genes on the (+) and (-) strand of the genome. Therefore, the following code was used to perform SAM -> BAM conversion: 
~~~
samtools view alignments.bam | gawk '(and(16, $2))' > forwardStrandReads.sam #forward strand read data
samtools view alignments.bam | gawk '(! and(16, $2))' > reverseStrandReads.sam #reverse strand read data
~~~

Before sorting the BAM file, I also found the need to add a Samheader manually to the SAM files.  Example: <br>
@HD	VN:1.0	SO:unsorted <br>
@SQ	SN:NC_005791.1	LN:1661137<br>
@PG	ID:Bowtie	VN:1.2.3	CL:"/opt/anaconda2/bin/bowtie-align-s --wrapper basic-0 -S MM_S2_genome <br> /Users/kershtheva/Desktop/Thesis/Chapter_4/forwardStrand_3_UTR_1.sam" <br>

Afterwards, I was able to make a sorted bam file using the following command line code: 
~~~
samtools view -bS -o forwardStrand_3_UTR_1.bam forwardStrand_3_UTR_1.sam
samtools sort forwardStrand_3_UTR_1.bam -o forwardStrand_3_UTR_1.sorted.bam
~~~
**Use IGV to convert BAM file into csv file for analysis**: Using IGVtools, I made an index file for the BAM file and then converted it into a .wig file, which could then be read by Microsoft Excel and saved as a csv file.

## Analysis of forward and reverse strand CATS-seq data
Processed CATS-seq data was integrated with gene tables for the organism to answer the two following  questions: <br>
1. What are the length distributions for primary, secondary and tertiary terminators? 
2. What % of terminators are found intragenically in this species? 

Answers to these two questions were plotted and are displayed at the top of this readme. To group terminators with genes, CATS seq reads were first averaged across triplicates and an algorithm was written to iterate through the genome and identify the 5' and 3' boundaries of CATS seq reads associated with each annotated gene. These boundaries were defined at edges where reads counts went from 0 to a positive number for 5' boundaries or vice versa for 3' boundaries. 
