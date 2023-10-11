#!/usr/bin/perl

# This program converts tRIBS connectivity table format to that read in by METIS.
#
# Usage: connectivity2metis.pl <input file> <output file> <node cnt> <extra weight flag> <nflux weight> <big_small_flag> <edge type>
#
# The input and output files can end in .gz, in which case they'll be
# piped through gzip or gunzip appropriately.
#

use strict;

sub Usage
{
   print("Usage: $0 <input file> <output file> <node cnt> <extra weight flag> <flux weight flag> <big_small_flag> <edge type>\n");
   print("       The files can be .gz files, which are handled transparently.\n");
}

my($inputFile, $outputFile, $open, $date, $line, @fields, $line, $cnt);
my($esize, $ine, $k, $jstart, $item, $nflux, $unode, %seen, @dup);
my($lcnt, $i, $j, $id, $part, $ndown, $jend, $index, $jindex, @headid, @outid, @edges, @vweight, @bweight, $array);
my($exWeightFlag, $bigWeightFlag, $nodeWeightFlag, $fluxWeightFlag, $edgeType, $nweights, @iweight, @fweight);

# If not three args, print the usage and exit.
if ($#ARGV != 6) {
   Usage();
   exit 1;
}

$inputFile = $ARGV[0];
$outputFile = $ARGV[1];
$nodeWeightFlag = $ARGV[2]; # Use reach number of nodes as weight (1) or not (0)
$exWeightFlag = $ARGV[3];   # Use reach inside/outside (1) or not (0)
$fluxWeightFlag = $ARGV[4]; # Use reach's number of flux reaches as a weight (1) or not
$bigWeightFlag = $ARGV[5];  # Use whether reach has > 2000 nodes
$edgeType = $ARGV[6];       # 0 = flow only, 1 = flux only, 2 = both

# Open an input filehandle.  If a compressed file is specified, pipe it
# through gunzip.
if ($inputFile =~ /.gz$/) {
   $open = "gunzip -c $inputFile |";
} else {
   $open = $inputFile;
}
open(INPUT, $open) || die("Unable to open input file $inputFile");

# Open an output filehandle.  If a compressed file is specified, pipe
# it through gzip.
if ($outputFile =~ /.gz$/) {
   $open = "| gzip -c > $outputFile";
} else {
   $open = "> $outputFile";
}
open(OUTPUT, $open) || die("Unable to open output file $outputFile");

# Write comment at the top.
chomp($date = `date`);
print(OUTPUT "% METIS input file generated from $inputFile on $date.\n");
print(OUTPUT "%\n");

# DEBUG: Simply copy the input to the output.
# while (<INPUT>) {
#    print OUTPUT;
# }

$line = <INPUT>; # skip header
chomp( $line = <INPUT>);
@fields = split(" ", $line);
$cnt = $fields[0]; # Number of reaches
#print "$cnt reaches\n";
# Init array for outside weights
for ($i = 0; $i < $cnt; $i++) {
    $iweight[$i] = 0;
    $fweight[$i] = 0;
    $bweight[$i] = 0;
}
$lcnt = 0;
for ($i = 0; $i < $cnt; $i++) {
   chomp($line = <INPUT>);
   @fields = split(" ", $line);
   $id = $fields[0];
   $vweight[$id] = $fields[1];
   if ($vweight[$id] > 2000) { $bweight[$id] = 1; }
   $headid[$id] = $fields[2];
   $outid[$id] = $fields[3];
   $ndown = $fields[4];
   # Read in downstream vertices
   $jstart = 5;
   $jend = $jstart;
   if ($ndown > 0) {
      $jend = $ndown + $jstart;
      for ($j = $jstart; $j < $jend; $j++) {
          # Save edges if edge type not "flux only"
          if ($edgeType != 1) {
	      push(@{$edges[$id]}, $fields[$j]);
	      push(@{$edges[$fields[$j]]}, $id);
 	      $lcnt++;
          }
              # Set as inside
              $iweight[$fields[$j]] = 1;
      }
   }
   # Read in flux vertices
   $nflux = $fields[$jend];
   # Set as flux weight
   $fweight[$i] = $nflux; 
   if ($nflux > 0) {
      $jstart = $jend + 1;
      $jend = $jstart + $nflux;
      for ($j = $jstart; $j < $jend; $j++) {
         # Check for duplicates before adding
         # Save edge if edge type is not "flow only"
         if ($edgeType > 0) {
             $ine = 0;
             if ($edgeType == 2) {
             $esize = scalar @{$edges[$id]};
             for ($k = 0; $k < $esize; $k++) {
                 if ($edges[$id][$k] == $fields[$j]) {
                     $ine = 1;
                 }
             }
             }
             if ($ine != 1) {
                 push(@{$edges[$id]}, $fields[$j]);
                 push(@{$edges[$fields[$j]]}, $id);
                 $lcnt++;
             }
         }
      }
   }
}
close(INPUT);

# Write the number of vertices and edges.  Currently, only vertices are
# weighted.
if ($edgeType == 1) { $lcnt /= 2; }
print(OUTPUT "% Next line lists # vertices and # edges.\n");
print(OUTPUT "$cnt $lcnt");
# Add number of weights
$nweights = 0;
if ($nodeWeightFlag == 1) { $nweights++; }
if ($exWeightFlag == 1) { $nweights++; }
if ($fluxWeightFlag == 1)  {$nweights++; }
if ($bigWeightFlag == 1)  {$nweights++; }
if ($nweights > 0) { print(OUTPUT " 10 $nweights\n"); }
else { print(OUTPUT "\n"); }

print(OUTPUT "%\n");

print(OUTPUT "% Format of vertex/edge entries:\n");
print(OUTPUT "%----------------------------------------------------------------------\n");
print(OUTPUT "% % <index of current vertex> <reach id> <weight> <head node id> <outlet node id>\n");
print(OUTPUT "% <vertex weight> <indices of neighbors of current vertex, separated by spaces>\n");
print(OUTPUT "%----------------------------------------------------------------------\n");
print(OUTPUT "%\n");

# For each vertex, write out the adjacent vertices.  The i'th line
# contains all neighbors of vertex i.
$index = 1;
for  ($i = 0; $i < $cnt; $i++) {
   print(OUTPUT "% $index $i $vweight[$i] $headid[$i] $outid[$i]\n");
   if ($nodeWeightFlag == 1) {
       print(OUTPUT "$vweight[$i] ");
   }

   # Add extra inside/outside weight if requested
   if ($exWeightFlag == 1) {
      print(OUTPUT "$iweight[$i] ");
   }

   # Add extra flux weight if requested
   if ($fluxWeightFlag == 1) {
      print(OUTPUT "$fweight[$i] ");
   }

   # Add big weight flag if requested   
   if ($bigWeightFlag == 1) {      
       print(OUTPUT "$bweight[$i] ");   
   }

   # Eliminate duplicates
   %seen = ();
   foreach $item (@{$edges[$i]}) {
      $seen{$item}++;
   }
   foreach $unode (keys %seen) {
      # Print the index number instead of label. 
      $jindex = $unode + 1;
      print(OUTPUT "$jindex ");
   }

   $index++;
   print(OUTPUT "\n");
}
