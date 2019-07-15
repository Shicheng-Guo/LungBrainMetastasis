use strict;
use List::Uniq ':all';
my @file=glob("*hg19_multianno.txt");

my %data;
my %gene;
my %sam;

foreach my $file(@file){
open F,$file;
my($sam)=split/.hg19_multianno.txt/,$file;
while(<F>){
my @line=split/\s+/;
$gene{$line[6]}=$line[6];
$sam{$sam}=$sam;
$data{$line[6]}{$sam}++;
}
}

my $sam=join("\t",sort keys %sam);
print "\t$sam\n";
foreach my $gene(sort keys %gene){
print "$gene";
foreach my $sam(sort keys %sam){
	if(defined $data{$gene}{$sam}){
		print "\t$data{$gene}{$sam}";
	}else{
		print "\t0";
	}
}
print "\n";
}
