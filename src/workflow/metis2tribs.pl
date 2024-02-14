#!/usr/bin/perl

# This program converts a metis partitioning to a tRIBS partition file.
#
# Usage: metis2tribs.pl <input file> <conn_file> <output file>
#

use strict;

sub Usage
{
   print("Usage: $0 <input file> <conn_file> <output file>\n");
}

my($inputFile, $connFile, $outputFile, $open, $date, $line, @fields, $line, $cnt);
my($i, $j, $rid, $pid, $maxpart, $p, $pcount, $totalWeight, @actualP, %actualP, @part, @weight, $array);

# If not 3 args, print the usage and exit.
if ($#ARGV != 2) {
   Usage();
   exit 1;
}

$inputFile = $ARGV[0];
$connFile = $ARGV[1];
$outputFile = $ARGV[2];

# Open an input filehandle. 
$open = $inputFile;
open(INPUT, $open) || die("Unable to open input file $inputFile");

# Open a connection filehandle.
$open = $connFile;
open(CONN, $open) || die("Unable to open input file $connFile");

# Open an output filehandle.
$open = "> $outputFile";
open(OUTPUT, $open) || die("Unable to open output file $outputFile");

# Collect all reaches per partition
$maxpart = 0;
$rid = 0;
while (<INPUT>)  {
    chomp;
    @fields = split(" ", $_);
    $pid = $fields[0];
    # Check fpr max partition
    if ($pid > $maxpart) { $maxpart = $pid; }
    # Check for existing partitions
    if (not exists($actualP{$pid})) {
        $actualP{$pid} = $pid;
    }
    push (@{$part[$pid]}, $rid);
    $rid++;
}
close(INPUT);

# Read in weights per node
while (<CONN>) {
    chomp;
    @fields = split(" ", $_);
    $pid = $fields[0];
    $weight[$pid] = $fields[1];
}
close(CONN);

$pcount = 0;
foreach $p (keys %actualP) { $pcount++; }
print "Actual number of partitions = $pcount\n";
print "Maximum partition = $maxpart\n";
print "Number of reaches = $rid\n";

# Write the tRIBS partition file.
# One partition per line.
# Write partition number followed by reach ids in
# that partition.
for  ($i = 0; $i <= $maxpart; $i++) {
   $array = $part[$i];
   $totalWeight = 0;
   for ($j = 0; $j <= $#{$array}; $j++) {
      # Print the index number instead of label.
      print(OUTPUT "$i $part[$i][$j]\n");
      $totalWeight += $weight[$part[$i][$j]];
   }
   print "Partition $i weight = $totalWeight\n"
}
close(OUTPUT);
