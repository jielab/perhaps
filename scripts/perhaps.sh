#!/bin/bash


# APOE: C-T;  e2: T-T;  e3: T-C;  e4: C-C
gendir=/mnt/d/projects/001UKB # the master directory that holds UKB genotic data

snps="rs429358 rs7412"
chr=19; begin=44908684; end=44908822 # GRCh38 positions 


### extract haplotype from phased data ###

echo $snps | tr ' ' '\n' > snps.txt

plink2 --pfile $gendir/hap/chr$chr --extract snps.txt --export vcf id-paste=iid bgz --out hap; tabix hap.vcf.gz

zcat hap.vcf.gz | awk '$1 !~/##/' | datamash -W transpose > hap.tmp

sed '1,9d; s/|/ /g' hap.tmp | awk '{print $1, $2$4, $3$5}' > hap.txt

sed 's/ 00/ e1/g; s/ 10/ e2/g; s/ 11/ e3/g; s/ 01/ e4/g; s/ //2' hap.txt > apoe.hap.txt


### piece haplotype based on paired end reads ###

echo -e "$chr\t$begin\t$end" > loc.bed

for dat in 2244305; do # 1466576
  for chr in 20; do # {1..22}; do

    samtools view ${dat}_23183_0_0.cram chr$chr > $dat.chr$chr.sam # sometimes without "chr"
	  awk '{print $9}' $dat.chr$chr.sam | sort -n > mate.len 
	  awk '$6=="151M" && NF==18 {$11="QQ"; if ($1 in reads) print reads[$1]" "$0; reads[$1]=$0}' $dat.chr$chr.sam > $dat.chr$chr.tmp.sam
	  awk 'NF !=36 || $4 > $22 {print NR, $4, $22}' $dat.chr$chr.tmp.sam | head # sanity check 
	  awk -v d=$dat '{print d, $1, $3, $4, $22, $22-$4}' $dat.chr$chr.tmp.sam | sort -k 6n > $dat.chr$chr.pairs.len
	  awk -v d=$dat -v b=$begin -v e=$end '{pos1=b-$4+1; pos2=e-$22+1; if (pos1>=1 && pos1<=76 && pos2>=1 && pos2<=76) { split($10,seq1,""); split($28,seq2,""); print d, seq1[pos1] "-" seq2[pos2]}}' $dat.tmp.sam > $dat.hap

  done
done




id=NA19146 ## sample ID
bamfile=$id.bam ## the location of the BAM or CRAM file, indexed
SNPs=19:44908684-44908822 ## the chr and positions of SNPs for phasing.

chr=${SNPs/:*/} # extract "chr" from the "SNPs" defined above 
pos=${SNPs/*:/} # extract "positions" of the "SNPs" defined above 
samtools view -O SAM -o $id.sam $bamfile chr$chr # extract the specified CHR and convert to SAM format   
readlen=`awk 'NR==1 {printf length($10)}' $id.sam`  # find the read length of the sequencing data

# remove reads with soft sequencing, extract first 10 fields
cut -f 1-10 $id.sam | awk '$6 !~/S/ {if ($1 in reads) print reads[$1]" "$0; reads[$1]=$0}' > $id.sam.paired 

# sanity check of haplotype size range formed by paired reads
awk '{if ($9<0) print -$9; else print $9}' $id.sam.paired | uniq | sort -n | uniq  > hap.len 

# this is the core script for PERHAPS
awk -v readlen=$readlen -v pos=$pos '{
	printf NR" "$1" ";
	split(pos,pa,"-");
	cnt=0; hap="";
	for (i in pa) {
		pos1=pa[i]-$4+1; 
		pos2=pa[i]-$14+1;
		if (pos1>=1 && pos1<=readlen) { split($10,seq1,""); printf "-SNP"i"-left("seq1[pos1]")" };
		if (pos2>=1 && pos2<=readlen) { split($20,seq2,""); printf "-SNP"i"-right("seq2[pos2]")" };
		if ((pos1>=1 && pos1<=readlen) || (pos2>=1 && pos2<=readlen)) cnt++; else printf "<>NA"
	};
	print " "cnt
}' $id.sam.paired | sed -e 's/-left//g' -e 's/-right//g' | sort -k 4,4nr -k 1,1n > $id.hap


# report the haplotypes that include all input SNPs
num=`echo $SNPs |  awk '{print gsub("-","") +1}'`
awk -v num=$num '$NF==num {$1=$2=""; print $0}' $id.hap | sort | uniq -c

# create a subset SAM file that only includes the reads that form the haplotype mentioned above, for IGV visualization
awk -v num=$num '$NF==num {print $2}' $id.hap > $id.subset.reads
fgrep -wf $id.subset.reads $id.sam > $id.subset.sam