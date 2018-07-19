#!/usr/bin/perl

# Script that retrieves the reads for a particular read group from a BAM file, and
# does a calculation of the average insert size

my $infile;
my $rg;
# my $tempfile;
# my $outfile;

if (scalar(@ARGV) != 2) {
	print "Usage: perl bam-rg-insert-size-calc.pl [infile] [read group]\n";
	exit(1);
}

$infile = shift(@ARGV);
$rg = shift(@ARGV);
# $tempfile = shift(@ARGV);
# $outfile = shift(@ARGV);

# my @insert_sizes = ();
my $n = 0;
my $sum = 0;
my @insert_sizes = ();

# my $samtools_cmd = "samtools view -r $rg -f 0x2 $infile | head -n 1000000 > $tempfile";
# system($samtools_cmd);

my @lines = `samtools view -r $rg -f 0x2 $infile | head -n 1000000`;

# open INFILE, "<$tempfile" or die "Can't open $tempfile: $!\n";
foreach my $line (@lines) {
	chomp($line);
	
	my @pieces = split(/\s+/, $line);
	
	# DEBUG
# 	print $pieces[2]."\t".$pieces[6]."\n";
# 	print $pieces[3]."\n";
	
	if ($pieces[6] eq "=") {
		my $insert_size = abs($pieces[7] - $pieces[3]);
		my $length = length($pieces[9]);
		# my $length = abs($pieces[8]);
		
		# DEBUG
		# print "Length: ".$length."\n";
		
		if ($insert_size > 0 && $insert_size <= 5000) {
			# push(@insert_sizes, $insert_size);
			
			# DEBUG
# 			if ($insert_size > 150000000) {
# 				next;
# 				# print $line."\n";
# 			}
			# print $insert_size."\n";
			
			push(@insert_sizes, ($insert_size+$length));
			$sum += ($insert_size + $length);
			$n++;
			if ($n >= 1000000) {
				last;
			}
		}
	}
}

# close(INFILE);
# open OUTFILE, ">$outfile" or die "Can't open $outfile: $!\n";
my $avg;
# print $n."\n";
if ($n == 0) {
	print "No insert sizes above zero\n";
} else {
	$avg = $sum/$n;
# 	open OUTFILE, ">$outfile" or die "Can't open $outfile: $!\n";
	print "Average: ".$avg."\n";
# 	close(OUTFILE);
}

my $median;
if ($n == 0) {
	# print "No insert sizes above zero\n";
} else {
	@insert_sizes = sort {$a <=> $b} @insert_sizes;
	if ($n % 2) { # Odd case
		$median = $insert_sizes[int($n/2)];
	} else { # Even case
		$median = ($insert_sizes[($n/2)-1]+$insert_sizes[($n/2)])/2;
	}
	# open OUTFILE, ">$outfile" or die "Can't open $outfile: $!\n";
	print "Median: ".$median."\n";
}

# close(OUTFILE);
exit();
