#!/usr/bin/perl -w
use strict;

my $reporte = "Lista_RE.txt";
open (FILE1, "$reporte");
my @Count1 = <FILE1>;
chomp @Count1;
my $header = shift @Count1;
my $c1 = scalar @Count1;
my @segmentos4 = split ("\t", $header);
chomp @segmentos4;
my $header1corregido = join(",", @segmentos4);

my $reporte2 = "TablaMaestraVanilla_Limpio.csv";
open (FILE2, "$reporte2");
my @Count2 = <FILE2>;
chomp @Count2;
my $header2 = shift @Count2;

my $c2 = scalar @Count2;

my $datafinal = "Anotada." . $reporte . ".csv";
open (OUTPUT, ">$datafinal");
my $headerF = "ID" . "," . $header1corregido . "," . $header2;
print OUTPUT "$headerF\n";

for (my $i = 0; $i < $c1; $i++) {
    my $renglon = $Count1[$i];
    chomp $renglon;
    my @segmentos = split ("\t", $renglon);
    chomp @segmentos;
    my $mew = join(",", @segmentos);
    my $gene = $segmentos[0];
    $gene =~ s/"//g;
    $gene =~ s/\s//g;

    my $match_found = 0;

    for (my $ii = 0; $ii < $c2; $ii++) {
        my $renglon2 = $Count2[$ii];
        chomp $renglon2;
        my @segmentos1 = split (",", $renglon2);
        chomp @segmentos1;
        my $gene2 = $segmentos1[0];
        $gene2 =~ s/\s//g;

        if ($gene eq $gene2) {
            my $linea = $mew . "," . $renglon2;
            print OUTPUT "$linea\n";
            $match_found = 1;
        }
    }

    # Si no encontró coincidencias, de todas formas imprime el ID con campos vacíos (opcional)
    if (!$match_found) {
        my $linea = $mew . ",";
        $linea .= join(",", ("") x (scalar(split(",", $header2))));
        print OUTPUT "$linea\n";
    }
}

close OUTPUT;
