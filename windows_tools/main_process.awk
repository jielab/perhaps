$3=="chr"c {
	printf NR" "$1" "
	split(pos,pa,"-");
	cnt=0;
	for (i in pa) {
		pos1=pa[i]-$4+1; pos2=pa[i]-$14+1;
		if (pos1>=1 && pos1<=readlen) { split($10,seq1,""); if (cnt!=0) printf "--"; printf "SNP"i"("seq1[pos1]")"};
		if (pos2>=1 && pos2<=readlen) { split($20,seq2,""); if (cnt!=0) printf "--"; printf "SNP"i"("seq2[pos2]")" };
		if ((pos1>=1 && pos1<=readlen) || (pos2>=1 && pos2<=readlen)) cnt++
	};
	print " "cnt
}