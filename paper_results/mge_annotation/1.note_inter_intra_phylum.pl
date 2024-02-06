open LS,"representive_genome.metadata.tsv";
<LS>;
while(<LS>){
	chomp;
	@l = split /\t/;
	$tax = $l[18];
	$phy = (split /;/,$tax)[1];
	$phy{$l[0]} = $phy;
}

open IN,"identified_event.merge_bp.csv";
<IN>;
while(<IN>){
	chomp;
	@l = split /,/;
	
	($g1,$g2) = ($l[2],$l[4]);
	$g1 =~ s/_\d+$//;
	$g2 =~ s/_\d+$//;
	if($phy{$g1} eq $phy{$g2}){$type = 'intra'}else{$type = 'inter'}
	print "$_,$phy{$g1},$phy{$g2},$type\n";
}

