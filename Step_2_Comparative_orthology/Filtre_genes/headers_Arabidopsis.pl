#!/usr/bin/perl -w
use strict;
 
my $File1 = "Lista_Arabidopsisv2.txt";
open(FILE1, "$File1");
my @archivo1 = <FILE1>;
chomp @archivo1;
my $c1 = scalar @archivo1;

my $File2 = "Arabidopsis_thaliana.TAIR10.pep.all.fa";
open(FILE2, "$File2");
my @archivo2 = <FILE2>;
chomp @archivo2;
my $c2 = scalar @archivo2;

my $Count2 = "headers_Arabidopsis_andy_OK.txt";
open(FILE3, ">$Count2");
#print FILE3 "$h\n";

for (my $i=0; $i < $c1; $i++){
    my $renglon1 = $archivo1[$i];
    chomp $renglon1;
    
        for (my $ii=0; $ii < $c2; $ii++){
        my $renglon2 = $archivo2[$ii];
        chomp $renglon2;
        
        if(($renglon2 =~ $renglon1) and ($renglon2 =~ ">") ){
        #if($ensemble2 =~ $renglon2   ){
            #print "$ensemble2\n";
            print "$renglon2\n";
            print FILE3 "$renglon2\n";
            last;
        }
        
    }
    
}