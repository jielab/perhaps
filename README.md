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


# #2. LINUX version of PERHAPS

```

## only the first 3 lines need to be changed
IID=NA20525 ## sample ID
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
	printf NR" "$1" "
	split(pos,pa,"-");
	cnt=0; hap="";
	for (i in pa) {
		pos1=pa[i]-$4+1; pos2=pa[i]-$14+1;
		if (pos1>=1 && pos1<=readlen) { split($10,seq1,""); printf "-SNP"i"-left("seq1[pos1]")" };
		if (pos2>=1 && pos2<=readlen) { split($20,seq2,""); printf "-SNP"i"-right("seq2[pos2]")" };
		if ((pos1>=1 && pos1<=readlen) || (pos2>=1 && pos2<=readlen)) cnt++; else printf "<>NA"
	};
	print " "cnt
}' $IID.sam.paired | sed -e 's/-left//g' -e 's/-right//g' | sort -k 4,4nr -k 1,1n > $IID.hap


# report the haplotypes that include all input SNPs
num=`echo $SNPs |  awk '{print gsub("-","") +1}'`
awk -v num=$num '$NF==num {$1=$2=""; print $0}' $IID.hap | sort | uniq -c

# create a subset SAM file that only includes the reads that form the haplotype mentioned above, for IGV visualization
awk -v num=$num '$NF==num {print $2}' $IID.hap > $IID.subset.reads
fgrep -wf $IID.subset.reads $IID.sam > $IID.subset.sam

```


# #3. Windows Version of PERHAPS

We also developed a script that could be run on Windows OS. The tools including awk.exe, cut.exe, samtools.exe, sort.exe, uniq.exe are needed, which are put under the windows-tools folder.
To run a test, you can input `python perhaps.py -i NA20525 -d .\test-data -s 1:159205564-159205704-159205737` on Windows cmd.

```
import os

def cmd_res(cmd: str):
    with os.popen(cmd) as p:
        return p.read()

def find_file(IID: str, dir: str):
    file = '\\'.join([dir, IID])
    if os.path.isfile(file + '.bam'):
        os.system(rf".\tools\samtools view {file + '.bam'} > {IID}.sam")
    elif os.path.isfile(file + '.sam'):
        os.system(rf"copy {file + '.sam'} {IID}.sam")
    else:
        raise Exception('File not find!')

def process(IID: str, rawfile: str, SNPs: str):
    [Chr, pos] = SNPs.split(':')
    os.system(rf'.\tools\samtools view {rawfile} > {IID}.sam')
    readlen = cmd_res(r'.\tools\awk "NR==1 {printf length($10)}" %s.sam' % IID).strip()

    os.system(r'.\tools\cut -f 1-10 %s.sam | '
              r'.\tools\awk "$6 !~/S/ {if ($1 in reads) print reads[$1]\" \"$0; reads[$1]=$0}" '
              '> %s.sam.paired' % (IID, IID))

    os.system(r'.\tools\awk "{if ($9<0) print -$9; else print $9}" %s.sam.paired'
              r' | .\tools\uniq | .\tools\sort -n | .\tools\uniq  > hap.len' % IID)

    os.system(r'.\tools\awk -v readlen=%s -v c=%s -v pos=%s -f main_process.awk'
              r' %s.sam.paired  > %s.hap' % (readlen, Chr, pos, IID, IID))

    os.system(r'.\tools\awk "($0 ~/SNP1/ && $0~/SNP2/) {$1=$2=\"\"; print $0}" %s.hap '
              r'| .\tools\sort | .\tools\uniq -c' % IID)

if __name__ == '__main__':
    IID = 'NA20525'  ## sample ID
    rawfile = rf'.\test-data\{IID}.bam'  ## the location of the BAM or CRAM file, indexed
    SNPs = '1:159205564-159205704-159205737'  ## the chr and positions of SNPs for directy haplotype detection.
    process(IID, rawfile, SNPs)

```


# #4. GUI version of PERHAPS.

After putting windows_tools directory in the same directory, then click perhaps.gui.exe to run the GUI version.
The default value is pre-filled, and users only need to click the "submit" button to get the same results as above.
Below are the screenshots of the GUI version.
 
![Figure 1](./Pictures/windows_1.png)

![Figure 2](./Pictures/windows_1.png)



# #5. Visualize and validate the directly called haplotypes.

Researchers could then open IGV (http://www.igv.org/) to visualize the genomic region in study and also visualize the directly called haplotype
 
![Figure 3](./Pictures/Figure1S.JPG)



# #6. extract phased haplotypes from PGEN file (for UKB dataset)

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


# #7. run Perhaps.R and more analyses to explore the haplotypes

![Figure 4](./Pictures/figure2.png)



# #8. download and extract ~50,000 WES data (for UKB dataset)
the UKB server allows no more than 10 jobs to download the WES data simultaneously for each approved project. 
Therefore, to download ~50,000 WES samples, we designed a strategy to put create ~500 list files, each containing links for 100 WES files.
Then we use LSF to run 9 ukbgene jobs in a batch, and put the other jobs in batches of 9 jobs, and waiting on the queue.
To save disk space, we delete the downloaded raw genotypic data after extracting the target regions.

```

dir=XXXX # the root directory for data
pwd=`pwd` # find the current dictory in a variable
awk '{print $1,"23164_0_0"}' $dir/data/ukb/wes/fe.fam | split -d -a 3 -l 100 - list  # extract the sample list from the WES plink file.


cnt=0 # this counting variable will be used in the following loop

for dat in `ls list*`; do  # to loop through all the list* files generated by the "split" command above
	
	let "cnt=$cnt+1"  # increase the count by 1
	let "gp=($cnt-1)/9"  # group every 9 jobs into a batch
	let "gp2=$gp +1"
	if [[ $gp2 == 1 ]]; then
		qsub_str="-N g$gp2.$cnt" 
	else
		qsub_str="-N g$gp2.$cnt -hold_jid g$gp.*"  # wait till the previous batch of jobs finished
	fi

	raw=${dat/list/raw}  
	outdir=$dir/data/ukb/wes/BAM/$raw; mkdir -p $outdir  # create a new folder for each group of downloading files
	mv $pwd/$dat $outdir
	
	echo "#!/bin/bash -l
	module load samtools
	sed 's/23164/23163/' $dat > $dat.2  # use the new data-field 23163 to download the associated cram.crai file.
	ukbfetch -b$dat  # !!! this line to actually download the .cram file
	ukbfetch -b$dat.2 # !!! this line to actually downoad the .cram.crai file
	for d in \`awk '{printf \" \"\$1}' $dat\`; do
		mv \${d}_23163_0_0.cram \$d.cram  # rename the file to a simpler name
		mv \${d}_23164_0_0.cram.crai \$d.cram.crai 
		samtools view -L $dir/files/apoe.b38.bed -O BAM -o \$d.bam \$d.cram # !!! this line extract the target region from the downloaded raw file 
		samtools index \$d.bam; rm \$d.cram*
		rm \$d.bam \$d.bam.bai # remove the downloaded raw file to save disk space.
	done
	" > $outdir/$dat.cmd
	
	cd $outdir
	qsub $qsub_str -o $dat.LOG -e $dat.ERR < $dat.cmd

done

```

