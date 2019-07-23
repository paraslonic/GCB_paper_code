#!/srv/common/bin/perl

use Bio::SeqIO;
use File::Basename;

my $fin = shift or die "[genbank file]";

$name = fileparse($fin, qr/\.[^.]*/);

$in = Bio::SeqIO->new(-file => $fin, -format => 'Genbank');
$id = 1;
while (my $seq = $in->next_seq()) {
	$contig = $seq->display_id();
	my @cds = grep { $_->primary_tag eq 'CDS' } $seq->get_SeqFeatures();
	foreach my $feature (@cds) {
		if($feature->has_tag('gene') ){
			my @gene = $feature->get_tag_values('gene');
			my @prod = $feature->get_tag_values('product');
			my @locus = $feature->location;
			my $prod = $prod[0];
			$prod =~ s/\s+/_/g; 
			$start  = $feature->start - 1;  # bed base is 0, while genbank is 1
			$end    = $feature->end;
			#if(@gene > 0) { $gene = $gene[1] } else { $gene = "" }
			#($name) = $feature->get_tag_values('locus_tag');
			$strand = $feature->strand;
			$strand = $strand < 0 ? '-' : '+';
			$gene = $gene[0];
			print join("\t", $contig, $start, $end, "${name}|$id|$prod","${gene}", 0, $strand), "\n";
			$id++;
		}
	}
}
