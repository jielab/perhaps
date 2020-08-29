$3=="chr"c {
	printf NR" "$1 # print row number and the read name in the first column
	split(pos,pa,"-"); # split the input positions into an array
	for (i in pa) { # for each input SNP position
		pos1=pa[i]-$4+1; pos2=pa[i]-$14+1; # the starting position of two read pairs are located in the 4th and 14th column, respectively
		if (pos1>=1 && pos1<=readlen) { split($10,seq1,""); printf " SNP"i"-left|"seq1[pos1]}; # call the genotype intercepted by the read upstream (on the left)
		if (pos2>=1 && pos2<=readlen) { split($20,seq2,""); printf " SNP"i"-right|"seq2[pos2]}; # call the genotype intercepted by the read downstream (on the right)
	};
	print "" # write a newline
}