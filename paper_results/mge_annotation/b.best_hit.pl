$input = "intra_recip_seq.IS.m6.anno";
$input = $ARGV[0];
open IN,$input;
while(<IN>){
	chomp;
	@l = split /\t/;
	$eid = $l[3];
	$all{$eid} = 1;;
}


foreach $eid(sort keys %all){
	$sort = `cat $input|grep $eid|sort -k 1,1 -k 15,15gr `;
	$select = `cat $input|grep $eid|sort -k 1,1 -k 15,15gr |head -1`;
	
	print $select;
}
