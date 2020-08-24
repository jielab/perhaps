# PERHAPS (Paired-End short Reads-based HAPlotyping from next-generation Sequencing data), is a new and simple approach to directly call haplotypes from short-read, paired-end Whole Exome Sequencing (WES) data. 
Author: Jie Huang, MD, PhD, Department of Global Health, Peking University School of Public Health


The technical bottleneck in direct haplotype calling from short-read sequencing lies in the length of sequenced DNA fragments, often too short to include two or multiple variable nucleotide positions that define the haplotype of interest. Indeed, while sequencing reads length in UKBB WES data is 76bp, the two APOE SNPs (rs7412 and rs429358), defining the common APOE polymorphism, are located 138 bp apart. We pieced short reads by utilizing their labels to generate a composite haplotype longer than 138bp. For illustraton purpose, we downloaded the WES data of two samples: HG01173 NA20525.


# Steps:

# #1. Download 1000 genomes sequencing data
start from 1000 genomes project main page https://www.internationalgenome.org. 
Then Click the "EBI FTP site" link under the "Alignments" section, and click "1000_genomes_project" link on the next page.
Users will directed to http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/. 
The "1000genomes.exome.GRCh38DH.alignment.index" file listed the FTP URL for 2,692 samples. 
The New York Genome Center (NYGC) released high-coverage (30x) data for a total of 3,202 samples. 
Users could download the aligned sequencing data for any set of samples from this link https://www.internationalgenome.org/data-portal/data-collection/30x-grch38, aligned to the GRCh38 reference genome. Once the CRAM file is downloaded, users could use samtools to extract certain regions of the genome to created a much smaller dataset, by using scripts such as below:

```
# create a 2genes.bed file with the following two rows, tab separated.
chr1    159204012       159206500       ACKR1
chr19   44905781        44909393        APOE

# run SAMTOOLS to extract the two genomic regions and create a new 2gene.bam file
samtools view -L 2genes.bed -O BAM -o 2genes.bam [raw-cram-file]
```

# #2. run the following code to piece together haplotypes from WES

```
## only the first 3 lines need to be changed
IID=HG01173 ## sample ID
rawfile=../BAM/$IID.bam ## the location of the BAM or CRAM file, indexed
SNPs=1:159205564-159205704-159205737 ## the chr and positions of SNPs for directy haplotype detection.

chr=${SNPs/:*/} # extract "chr" from the "SNPs" defined above 
pos=${SNPs/*:/} # extract "positions" of the "SNPs" defined above 
samtools view -O SAM -o $IID.sam $rawfile # convert the BAM/CRAM file to SAM format (txt format)   
readlen=`awk 'NR==1 {printf length($10)}' $IID.sam`  # find the read length of the sequencing data

# remove reads with soft sequencing, extract first 10 fields
cut -f 1-10 $IID.sam | awk '$6 !~/S/ {if ($1 in reads) print reads[$1]" "$0; reads[$1]=$0}' > $IID.sam.paired 

# sanity check of haplotype size range formed by paired reads
awk '{if ($9<0) print -$9; else print $9}' $IID.sam.paired | uniq | sort -n | uniq  > hap.len 

# this is the core script for PERHAPS
awk -v readlen=$readlen -v c=$chr -v pos=$pos '$3=="chr"c {
	printf NR" "$1 # print row number and the read name in the first column
	split(pos,pa,"-"); # split the input positions into an array
	for (i in pa) { # for each input SNP position
		pos1=pa[i]-$4+1; pos2=pa[i]-$14+1; # the starting position of two read pairs are located in the 4th and 14th column, respectively
		if (pos1>=1 && pos1<=readlen) { split($10,seq1,""); printf " SNP"i"-left|"seq1[pos1]}; # call the genotype intercepted by the read upstream (on the left) 
		if (pos2>=1 && pos2<=readlen) { split($20,seq2,""); printf " SNP"i"-right|"seq2[pos2]}; # call the genotype intercepted by the read downstream (on the right)
	};
	print "" # write a newline
}' $IID.sam.paired  > $IID.hap

# report the number of haplotypes that meet certain critia, for example, those read pairs that overlap with SNP1 + SNP2, or more
awk '($0 ~/SNP1/ && $0~/SNP2/) {$1=$2=""; print $0}' $IID.hap | sort | uniq -c 

```

# #3. Visualize and validate the directly called haplotypes.

Researchers could then open IGV (http://www.igv.org/) to visualize the genomic region in study and also visualize the directly called haplotype
 
![Figure 1](./Pictures/Figure1.PNG)



# #4. extract phased haplotypes from PGEN file (for UKB dataset)
```
gendir=XXX # the directory for the UKB haplotypes file
snps="rs429358 rs7412"
chr=19; begin=44908684; end=44908822 # GRCh38 positions for two SNPs that define the APOE haplotype

###. extract haplotype from phased data ###
echo $snps | tr ' ' '\n' > snps.txt

plink2 --pfile $gendir/hap/chr$chr --extract snps.txt --export vcf id-paste=iid bgz --out hap; tabix hap.vcf.gz
zcat hap.vcf.gz | awk '$1 !~/##/' | datamash -W transpose > hap.tmp
sed '1,9d; s/|/ /g' hap.tmp | awk '{print $1, $2$4, $3$5}' > hap.txt
sed 's/ 00/ e1/g; s/ 10/ e2/g; s/ 11/ e3/g; s/ 01/ e4/g; s/ //2' hap.txt > apoe.hap.txt

```

# #5. run Perhaps.R and more analyses to explore the haplotypes

![Figure 2](./Pictures/Figure2.jpg)

