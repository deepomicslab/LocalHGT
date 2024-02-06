open LS,"identified_event.csv";
<LS>;
while(<LS>){
	chomp;
	@l = split /,/;
	$id1 = "$l[2]-$l[3]";
	$id2 = "$l[4]-$l[5]";
	$id3 = "$l[4]-$l[6]";
	foreach $id($id1,$id2,$id3){
		($g,$p) = split /-/,$id;
		$pos{$g}{$p} = 1;
	}
}

foreach $g(keys %pos){
	%tmp = %{$pos{$g}};
	@pos = sort {$a<=>$b} keys %tmp;
	$p0 = $pos[0];
	@tmp = ($p0);
	for($i=1;$i<@pos;$i++){
		$p = $pos[$i];
		if($p - $p0 < 5000){ # origin 20
			push (@tmp,$p);
			$p0 = &aver(@tmp);
		}else{
			foreach $in_p(@tmp){
				$pos{$g}{$in_p} = $p0;
			}
			$p0 = $p;
			@tmp = ($p0);
		}
	}
	foreach $in_p(@tmp){
		$pos{$g}{$in_p} = $p0;
	}
}

open LS,"identified_event.csv";
$h = <LS>;
print $h;
while(<LS>){
	chomp;
	@l = split /,/;
	$l[3] = $pos{$l[2]}{$l[3]};
	$l[5] = $pos{$l[4]}{$l[5]};
	$l[6] = $pos{$l[4]}{$l[6]};
	$new = join(",",@l);
	print "$new\n";

}

sub aver {
	$n = scalar(@_);
	$sum = 0;
	foreach $item (@_) {
		$sum += $item;
	}
	$average = int($sum / $n) + 1;
	return $average;
}
