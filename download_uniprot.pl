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
open (SAIDA, ">>$file\.GO.tab");
foreach $lstprotein (@lstproteins){
	my $url1 = "https://www.uniprot.org/uniprot/?query=proteome:$lstprotein&format=fasta"
	#my $url1 = "http://www.uniprot.org/uniprot/?query=$lstprotein&columns=id,entry,genes,go-id,go(biological%20process),go(molecular%20function),go(cellular%20component)&format=tab";
	my $output1 = get($url1);
	#$output1=~s/Entry.*\n//g;
	#$output1=~s/\n//g;	
	print SAIDA ("$output1");
}
close SAIDA;