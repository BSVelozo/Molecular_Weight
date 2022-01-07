#!/usr/bin/perl

## Laboratorio de Bioinformatica da UFRJ
## Aluno: Bernardo Velozo e Yuri
## Uso: perl .\Molecular_Weigth.pl .\File_with_ptns.fasta


use strict;
use warnings;
use Encode;

for my $file (@ARGV) {
    open my $fh, '<:encoding(UTF-8)', $file;
    my $input = join q{}, <$fh>; 
    close $fh;

    my $filename = 'proteins.fasta';

    open(FC, ">", $filename);  
    
    while ( $input =~ /^(>.*?)$([^>]*)/smxg ) {
        my $name = $1;
        my $seq = $2;
        my $outfile;
        my @list;
        my $n=0;
        $seq =~ s/\n//smxg;

        my $mass = calc_mass($seq);

        ## Só alterar aqui para fazer a comparação de massa
        if ($mass >= 35 and $mass <= 45) {
          my $selectSeq = $name . $mass . "Kda" . "\n" . $seq . "\n" ;
          print FC $selectSeq;
	    }
    }

    close(FC);
}

sub calc_mass {
  my $a = shift;
  my @a = ();
  my $x = length $a;
  @a = split q{}, $a;
  my $b = 0;

  my %data = (
      A=>0.089,  R=>0.174,  D=>0.133,  N=>0.132,
      C=>0.121,  E=>0.147,  Q=>0.146,  G=>0.075,
      H=>0.155,  I=>0.131,  L=>0.131,  K=>0.146,
      M=>0.149,  F=>0.165,  P=>0.115,  S=>0.105,
      T=>0.119,  W=>0.204,  Y=>0.181,  V=>0.117
  );

  for my $i( @a ) {
      $b += $data{$i};
  }
 $b = $b - (0.018 * ($x - 1));
 $b = int($b + 0.5);
 return $b;
}