#!/usr/bin/perl
use Data::Dumper;
use List::Util qw(sum);
use FindBin;
use lib "$FindBin::Bin/../lib";
use List::MoreUtils 'pairwise';


my @filetypes=("matrix.txt","filter.nompval.0.05.matrix.txt","filter.fdr.0.25.matrix.txt");
#my @filetypes=("matrix.txt");
my @gene_set_types=("c1.all","c2.cgp","c2.cp.biocarta","c2.cp.kegg","c2.cp.reactome","c3.mir","c3.tft","c4.cgn","c4.cm","c5.bp","c5.cc","c5.mf","c6.all","c7.all");
my @sexes=("na_pos");
my %matrix;

foreach my $gene_set_type (@gene_set_types) {
	foreach my $sex (@sexes) {
		foreach my $filetype (@filetypes) {
#c7.all.top20_gene_sets.males.filter.nompval.0.05.matrix.txt
			my $file=join('.',$gene_set_type,'top20_gene_sets',$sex,$filetype);
			open(my $fh,"<",$file);
			$header=<$fh>;
			while(<$fh>) {
				@values=split(" ",$_);
				my $gene_set=shift(@values);
				#unless(exists($matrix{$sex}{$gene_set_type}{$gene_set})) { $matrix{$sex}{$gene_set_type}{$gene_set}=[qw(0 0 0 0)]};
				unless(exists($matrix{$sex}{$gene_set_type}{$gene_set})) { $matrix{$sex}{$gene_set_type}{$gene_set}=[qw(0 0 0 0 0 0 0 0)]};
				if($filetype eq 'filter.fdr.0.25.matrix.txt') {@values = map { $_ * 2 } @values;}
				$matrix{$sex}{$gene_set_type}{$gene_set} = [pairwise { $a + $b } @{$matrix{$sex}{$gene_set_type}{$gene_set}}, @values];
			}
			close($fh);
		}
	
		my $fileout=join('.',$gene_set_type,'top20_gene_sets',$sex,'total.matrix.txt');
		open ($out,'>',$fileout);
		print $out "Gene_set immvar.CD14 immvar.CD4 Fairfax GenCord MESA.M MESA.T meta.CD14 meta.CD4\n";
		for my $key ( sort { sum(@{$matrix{$sex}{$gene_set_type}{$b}}) <=> sum(@{$matrix{$sex}{$gene_set_type}{$a}})  } keys %{$matrix{$sex}{$gene_set_type}} )
        	{
			my $string = join(" ", $key,@{$matrix{$sex}{$gene_set_type}{$key}});
                	print $out $string."\n"; 
        	}
		close($out);
	}
}
#print Dumper(%matrix);

