#!/usr/bin/perl

# DESCRIPTION:	obtendo o gene_name depositado no Uniprot das proteínas vindas do stringDB.

use strict;
use LWP::Simple; 

my $usage="\nCommand Line:\n$0 [Arquivo de PPI do stringDB] \n\n";

my $file=$ARGV[0] || die "$usage";

#system("sed -i 's/ /\t/g' $file");
open (FILE, "$file");
my $line; 
my @lstproteins; 
my $lstprotein; 
while ($line=<FILE>){
	chomp ($line);
	push (@lstproteins,$line)
}
close FILE;

foreach $lstprotein (@lstproteins){
    my @ids_name= split('\t', $lstprotein);
    my $id=@ids_name[0];
    my $specie=@ids_name[1];
    open (SAIDA, ">>$specie.fasta");
	my $url1 = "https://www.uniprot.org/uniprot/?query=proteome:$id&format=fasta";
	my $output1 = get($url1);
	#$output1=~s/Entry.*\n//g;
	#$output1=~s/\n//g;	
	print SAIDA ("$output1");
    close SAIDA;
}
