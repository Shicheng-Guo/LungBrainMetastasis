#!/usr/bin/perl
use strict;
chdir "/home/guosa/hpc/project/LungBrainMetastasis/vaf";
my @file=glob("*vaf");

my %data;
my %gene;
my %sam;
my %pos;

foreach my $file(@file){
open F,$file;
my($sam)=split/.vafs/,$file;
while(<F>){
my @line=split/\s+/;
my $pos="$line[0]:$line[1]";
$pos{$pos}=$pos;
$sam{$sam}=$sam;
$data{$sam}{$pos}=$line[2];
}
}

my $sam=join("\t",sort keys %sam);
print "\t$sam\n";
foreach my $pos(sort keys %pos){
print "$pos";
foreach my $sam(sort keys %sam){
if(defined $data{$sam}{$pos}){
print "\t$data{$sam}{$pos}";
}else{
print "\t0";
}
}
print "\n";
}
