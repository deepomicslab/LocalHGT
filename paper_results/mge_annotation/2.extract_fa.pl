#$id = 'GUT_GENOME000055_5';
$id = $ARGV[0];

($g) = $id =~ /(GUT_.*)_/;

#genome folder contains the gff files download from UHGG database
open IN,"gunzip -dc genome/$g.gff.gz|";

$/ = ">";
<IN>;
while(<IN>){
	s/>$//;
	@line = split /\n/;
	if($line[0] eq $id){
		open OU,">scaffold/$id.fa";
		print OU ">$_";
		close OU;
	}
}

