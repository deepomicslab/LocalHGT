




($input,$dbtype) = @ARGV;
open LS,"cat ../ref/$dbtype.seq_len|";
while(<LS>){
	chomp;
	@l = split /\t/;
	$len{$l[0]} = $l[1];
}


open IN,"$input.fa";
while($id = <IN>){
	$id = (split / /,$id)[0];
	$id =~ s/>//;
	$seq = <IN>;
	chomp $seq;
	$len{$id} = length($seq);
}

open IN,"$input.$dbtype.m6";
while(<IN>){
	chomp;
	@l = split /\t/;
	($eid,$dbid,$idt,$len,$evalue) = ($l[0],$l[1],$l[2],$l[3],$l[10]);
	if($evalue >= 0.00001){next}
	$full_len = $len{$dbid};
	if($len/$full_len > 0.9 ){
		$type = "A";
	}else{
		if($len > 100){
			$type = 'B';
		}else{
			$type = 'C';
		}
	}
	print "$type\t$len{$eid}\t$len{$dbid}\t$_\n";
}
