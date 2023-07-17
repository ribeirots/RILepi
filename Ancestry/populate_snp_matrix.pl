use strict ;
use warnings ;

## minimum distance between alleles shoudl be set to the read length
my $min_dist = 150 ;

### load parental SNPs
my %alleles ;
foreach (`cat $ARGV[0]`) {
chomp ;
my @split = split ( /\t/, $_ ) ;
$alleles{$split[1]} = $_ ;
}

my $last = 0 ;

### read from mpileup file, do not give samtools a reference when you make this
while (<STDIN>) {
chomp $_ ;
my @split = split ( /\t/, $_ ) ;
if ( !exists( $alleles{$split[1]} ) || $split[1] - $last < $min_dist ) {
next ;
}

my @a = split ( /\t/, $alleles{$split[1]} ) ;

### this uses a uniform recombination rate, 1e-8 per site, or 1cm/mb
print $a[0], "\t", $a[1], "\t", $a[4], "\t", $a[5], "\t", $a[6], "\t",
$a[7], "\t", ( $a[1] - $last ) * 1e-8 ;
$last = $a[1] ;

for ( my $i = 4 ; $i < $#split + 1 ; $i += 3 ) {

my $c1 = 0 ; my $c2 = 0 ;

for ( my $p = 0 ; $p < length( $split[$i] ) ; $p ++ ) {
if ( substr( $split[$i], $p, 1 ) =~ m/$a[2]/i ) {
$c1 ++ ;
}
elsif ( substr( $split[$i], $p, 1 ) =~ m/$a[3]/i ) {
                       $c2 ++ ;
}
}

print "\t", $c1, "\t", $c2 ;

}

print "\n" ;

}