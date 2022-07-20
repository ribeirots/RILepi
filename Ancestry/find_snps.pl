use strict ; 
use warnings ; 

my %a ;
@{$a{"Y"}} = ("C","T") ;
@{$a{"R"}} = ("A","G") ;
@{$a{"W"}} = ("A","T") ;
@{$a{"S"}} = ("C","G") ;
@{$a{"K"}} = ("T","G") ;
@{$a{"M"}} = ("A","C") ;

my $ef_seq = `cat $ARGV[0]` ; 
my $zi_seq = `cat $ARGV[1]` ; 
my $chrom = $ARGV[2] ;

foreach ( my $i = 0 ; $i < length($ef_seq) ; $i ++ ) { 
	if ( substr($ef_seq,$i,1) eq "N" || substr($zi_seq,$i,1) eq "N" ) { 
		next ;
	}
	elsif ( substr($ef_seq,$i,1) ne substr($zi_seq,$i,1) ) { 

		my %al ; 		

		if ( exists( $a{substr($ef_seq,$i,1)} ) ) { 
			$al{${ $a{substr($ef_seq,$i,1)} }[0]} ++ ;
			$al{${ $a{substr($ef_seq,$i,1)} }[1]} ++ ;
		}
		else {
			$al{substr($ef_seq,$i,1)} ++ ;
		}

		if ( exists( $a{substr($zi_seq,$i,1)} )	) { 
                        $al{${ $a{substr($zi_seq,$i,1)} }[0]} ++ ;
                        $al{${ $a{substr($zi_seq,$i,1)} }[1]} ++ ;
                } 
                else {
                        $al{substr($zi_seq,$i,1)} ++ ;
                }

		if ( scalar keys %al != 2 ) {
			next ;
		}

		my @alleles = sort keys %al ; 

		print $chrom, "\t", $i + 1, "\t", $alleles[0], "\t", $alleles[1] ;
		
		if ( substr($ef_seq,$i,1) eq $alleles[0] ) { 
			print "\t", 1, "\t", 0 ; 
		}
		elsif ( substr($ef_seq,$i,1) eq $alleles[1] ) { 
			print "\t", 0, "\t", 1 ; 
		}
		else {
			print "\t", 1, "\t", 1 ; 
		}

		if ( substr($zi_seq,$i,1) eq $alleles[0] ) { 
                        print "\t", 1, "\t", 0 ; 
		}
                elsif (	substr($zi_seq,$i,1) eq $alleles[1] ) {	
                        print "\t", 0, "\t", 1 ; 
                }
                else {
                        print "\t", 1, "\t", 1 ; 
                }
		print "\n" ;
	}
}
