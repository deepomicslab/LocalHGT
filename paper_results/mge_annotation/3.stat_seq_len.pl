@file = glob("scaffold/*.fa");
foreach $file(@file){
	($genome) = $file =~ /\/(.*)\.fa/;
	$seq = `cat $file|sed 1d`;
	$seq =~ s/\s//g;
	$len = length($seq);
	print "$genome\t$len\n";
}
