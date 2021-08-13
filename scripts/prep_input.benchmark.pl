#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;

open (INDEL, "-|", "tail -n +2 $ARGV[0]") || die "Cannot open file: $!";
my ($pname, $path, $suffix) = fileparse($ARGV[0], qr/\.[^.]*/);

foreach (<INDEL>) {

	chomp $_;

	my @data = split("\t", $_);

   	my $chr = $data[0];
	my $pos = $data[1];
	my $mapq = $data[18];
	my $ddg2p_annote = $data[24];
	my $perct_disc = $data[20];
	my $blast_hit = $data[32];
	my $sr_total = $data[5];
	my $cov = $data[2] + 1;
	my $exonic = $data[27];
	my $mum_sr = $data[38];
	my $dad_sr = $data[39];
	my $prob_y = $data[22];
	my $maf = $data[29];
	my $size = $data[33];
	my $transcript = $data[28];
	
	my $mum_tot;
	my $dad_tot;
	if ($mum_sr eq "NA" && $dad_sr eq "NA") {
		$mum_tot = 0;
		$dad_tot = 0;
	} elsif ($mum_sr eq "NA" & $dad_sr ne "NA") {
		$mum_tot = 0;
		$dad_tot = $dad_sr;
	}elsif ($mum_sr ne "NA" & $dad_sr eq "NA") {
		$mum_tot = $mum_sr;
		$dad_tot = 0;
	} else {
		$mum_tot = $mum_sr;
		$dad_tot = $dad_sr;
	}

   	my $depth = $sr_total / $cov;

	if ($mapq >= 20) {
		if ($perct_disc < 0.5 && (($perct_disc > 0.1 && $blast_hit ne "UNK" && $blast_hit ne "TRANSSEGDUP") || ($perct_disc <= 0.1))) {
			if ($depth >= 0.1) {
				if ($mum_tot < 2 && $dad_tot < 2) {
					print $_ . "\n";
				}
			}
		}
	}
}
