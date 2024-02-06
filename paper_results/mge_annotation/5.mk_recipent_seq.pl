open LS,"scaffold_len";
while(<LS>){
	chomp;
	($g,$len) = split /\t/;
	$len{$g} = $len;
	
}

$type = $ARGV[0]; # inter or intra


open LS,"cat identified_event.merge_bp.anno_type |grep $type|";
$len = 5000;
while(<LS>){
	chomp;
	@l = split /,/;
	$id = "event$l[0]_recip";
	($g1,$p0,$g2,$p1,$p2) = ($l[2],$l[3],$l[4],$l[5],$l[6]);

	if($p0-$len-1 <=0 ){$pp1 = 1}else{$pp1 = $p0-$len-1}
	if($p0+$len > $len{$g1}){$pp2 = $len{$g1}}else{$pp2 = $p0+$len}
	
	$bed = join("-",($g1,$pp1,$pp2,"False"));

	$out = ">$id $_\n";

	foreach $bed($bed){
		($tg,$tp1,$tp2,$rev) = split /-/,$bed;
		if($rev eq 'False'){
			$seq = `samtools faidx scaffold/$tg.fa $tg:$tp1-$tp2|sed 1d`;
		}else{
			$seq = `samtools faidx scaffold/$tg.fa -i $tg:$tp1-$tp2|sed 1d`;
		}
		$seq =~ s/\n//g;
		$out .= "$seq";
	}
	print "$out\n";
}
