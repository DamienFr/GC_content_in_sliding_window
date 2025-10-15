#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;

=head1 Description

    This script will first calculate, for each nucleotide position, the %GC of the surrounding nucleotid sequence according to the window size.
    For the first nucleotid, half of the window size will be taken from the end of the sequence, because this script has been developped for cicular genoms.
    Then, i will print a value of % each "step" nucleotides
	
=head1 Usage

    perl GC_content.pl -fasta genome.fasta [-window 1000] [-step 100] [-log] [-help] 

    Arguments details :
    -fasta	Fasta file, (REQUIRED) Can contain multiple sequences: script will produce as many output files as input sequences. Each sequence can be multiline.
    -window 	Window size, (optional, default 1000) Even number ONLY. Sets the number of nucleotides used to calculate the %GC value of each position.
    -step	Step size, (optional, default 100) The output will contain a sliding GC% value every "step" nucleotides. The numbers of values you get is therefore (length genome)/(step).
    -log	(optional) For debugging purposes only
    -help	(optional) Shows this help
   
=head1 Output files

    [input_name].GC_content
    
    column 1: "chr1" for convenience when the program is used to create an input for Circos (www.circos.ca/) 
    column 2: start position of interval
    column 3: end position of interval 
    column 4: GC% of the interval 
    
    example:
    chr1 1 100 0.499
    chr1 101 200 0.5
    chr1 201 300 0.5
    chr1 301 400 0.5
    chr1 401 500 0.5
    chr1 501 600 0.5
    chr1 601 700 0.499
    chr1 701 800 0.5

    [input_name].GC_deviation

    column 1: "chr1" for convenience when the program is used to create an input for Circos (www.circos.ca/) 
    column 2: start position of interval
    column 3: end position of interval 
    column 4: GC deviation of the interval 
    
    example:
    chr1 1 100 0.162
    chr1 101 200 0.16
    chr1 201 300 0.16
    chr1 301 400 0.16
    chr1 401 500 0.16
    chr1 501 600 0.16
    chr1 601 700 0.162
    
=head1 AUTHOR

	Damien Richard, 2019

=cut

my ($fasta,$help,$log,%h_pcent, $size_fasta, $sequence, $not_sequence,%h_gcdev,%hhh,$title_provisoire,$last_end,$last_end2);

GetOptions(
	"fasta=s" => \$fasta,
	"window=s" => \(my $window = 1000),
	"step=s" => \(my $step = 100),
	"log" => \$log,
	"help" => \$help
  );

die `pod2text $0` unless $fasta;
die `pod2text $0` if $help;

$fasta =~ s/\r?\n//g;
$window =~ s/\r?\n//g;
$step =~ s/\r?\n//g;

if ($window % 2 != 0){print "\n!!! Warning : window must be an even number... !!!\n\n"; die `pod2text $0`}

if($log){ open (LOG, ">", $fasta . ".log") or die "can't open $!";}

# reading input fasta file WITHOUT using perl module Bio::SeqIO;
#open (IN, "<", $fasta) or die "can't open  $!";
#$not_sequence = 0;
#while (my $line = <IN>) {
#	$line =~ s/\r?\n//g;
#	$line =~ s/\s//g;
#	if($line =~ /^>/){ $line =~ s/^>//; if(exists($hhh{$line})){ print "\nTwo sequence titles have the exact same name : $line, i'll die now ...\n"; die } ;
#	$title_provisoire = $line}else{ $hhh{$title_provisoire} .= $line }
#}
#close IN;

# reading input fasta file WITH Bio::SeqIO;
my $fasta_in  = Bio::SeqIO->new( -file => $fasta, -format => 'Fasta' );

while ( my $seq = $fasta_in->next_seq() ) {
$title_provisoire = $seq->id;
    if ( exists( $hhh{$seq->id} ) ) { print "\nTwo sequence titles have the exact same name : $title_provisoire, i'll die now ...\n"; die ;}
$hhh{$seq->id} = $seq->seq();
}

#### loop on all sequences of the input fasta file ####
foreach my $kkk ( keys %hhh ){
my $kkk2 = $kkk;

print "\n\tWorking on sequence $kkk\n\n";

open (GCPC, ">", $fasta . "_" . $kkk2 . ".GC_content" ) or die "can't open $!";
open (GCDEV, ">", $fasta . "_" . $kkk2 . ".GC_deviation") or die "can't open $!";

my $sequence = $hhh{$kkk};
$size_fasta = length($sequence);
if ($sequence =~ m/[^ATGC]/){print "Your fasta sequence contains other caracters than [ATGC], i'll keep working normally, it's just so you know ...\n ";}

#### Calculating mean GC content ####
my $count1 = 0;
while ($sequence =~ /[GC]/g) { $count1++ }
my $pcent =  $count1 / $size_fasta * 100 ;
print "\tSequence\tLength\tMean GC content\n\t$kkk\t$size_fasta\t$pcent % \n";
my $line_mod = (substr $sequence, -$window/2 ) . $sequence . (substr $sequence, 0, ($window / 2));
print LOG $line_mod . "\n" if $log;

my $start = 0;

until ($start > (length($sequence)))
{
	my $count = 0;
	my $countG = 0; my $countC = 0;  #GC_deviation
	my $tmp = substr $line_mod, $start, $window;
	my $position = $start +1 ;
	while ($tmp =~ /[GC]/g) {$count++}
	while ($tmp =~ /[C]/g) {$countC++} #GC_deviation
	while ($tmp =~ /[G]/g) {$countG++} #GC_deviation
	my $gcdev;
	if(($countG+$countC) > 0){$gcdev = ($countG-$countC)/($countG+$countC)}else{ $gcdev = 0 } ; #GC_deviation
	#print LOG $count . " " . $tmp . "\n";
	my $pcent =  $count / $window ;
	$h_pcent{$position} = $pcent;
	$h_gcdev{$position} = $gcdev; #GC_deviation
	print LOG $position . " " . $tmp . " $pcent\n" if $log;
	$start ++;  # $start += $window; if i choose to calculate the mean on the windw instead of making it for each nucleotid
}

my $k = 0;
my $acc = 0;

foreach my $key2 ( sort {$a<=>$b} keys %h_pcent)
{
	if ($key2 % $step == 0 || $key2  == ($size_fasta + 1)) {
		my $start_window2 = $key2 - ($step - 1);
		my $end_window2 = $key2; if($end_window2 >= $size_fasta){$start_window2 = $last_end + 1; $end_window2 = $size_fasta }
		print GCPC "chr1 " . $start_window2 . " " . $end_window2 . " " . $h_pcent{$key2} . "\n";
		$last_end = $end_window2;
	}
}

foreach my $key3 ( sort {$a<=>$b} keys %h_gcdev)
{
	if ($key3 % $step == 0 || $key3  == ($size_fasta + 1)) {
		my $start_window2 = $key3 - ($step - 1);
		my $end_window2 = $key3; if($end_window2 >= $size_fasta){$start_window2 = $last_end2 + 1; $end_window2 = $size_fasta }
		print GCDEV "chr1 " . $start_window2 . " " . $end_window2 . " " . $h_gcdev{$key3} . "\n";
		$last_end2 = $end_window2;
	}
}

close GCDEV;
close GCPC;
print "\n\tOutput files are :\n\t" . $fasta . "_" . $kkk2 . ".GC_content\n\t" . $fasta . "_" . $kkk2 . ".GC_deviation\n";
}

if($log){close LOG;}

unlink $fasta . ".log" if !$log;

