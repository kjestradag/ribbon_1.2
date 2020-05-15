#!/usr/bin/env perl
use List::Util qw( min max );
use strict;

my $mast= shift; # genome_windows.mast (mast result)
my $org= shift; # $dir
my $genome= shift; #$genome
my ($seq, $names, %data, %fildata);

ext_data_from_mast();
fil_parts();
extract_16S();
print_structure();

sub ext_data_from_mast{
    open RESUME, ">${org}_windows.mast.resume";
    print RESUME "$org\n";
    ($seq, $names)= read_fasta("$genome", '^>(\S+)');

    open IN, $mast or die "Cant read $mast\n";
    while(<IN>){
        chomp;
        if( /^(\S+)\.(\d+)\s+([\+\-])(\d+)\s+\d+\s+\w+\s+(\d+)\s+(\d+)/ ){ # mast -v == 5.0.1
            my($contig, $part, $strand, $motif, $start, $end)= ($1, $2, $3, $4, $5, $6);
            push @{$data{$contig}{$part}{motif}}, $motif;
            $data{$contig}{$part}{strand}= $strand;

            my ($sm4, $sm3, $em2, $em3);
            $start-4 >= 0 ? ($sm4= $start-4) : ($sm4= 0);
            $start-3 >=0 ? ($sm3= $start-3) : ($sm3= 0);
            $end+2 <= length($seq->{$contig}) ? ($em2= $end+2) : ($em2= length($seq->{$contig}));
            $end+3 <= length($seq->{$contig}) ? ($em3= $end+3) : ($em3= length($seq->{$contig}));

            $strand eq '+' ? (push @{$data{$contig}{$part}{start}}, $sm4) : (push @{$data{$contig}{$part}{start}}, $sm3);
            $strand eq '+' ? (push @{$data{$contig}{$part}{end}}, $em2) : (push @{$data{$contig}{$part}{end}}, $em3);
        }
    }
}

sub fil_parts{
    foreach my $contig ( sort keys %data ){
        my $i= 0;
        my %ya;
        foreach my $part ( sort keys %{$data{$contig}} ){
            my $first= $part;
            my $from= ($part-1)*2000+ min(@{$data{$contig}{$part}{start}});
            my $to= ($part-1)*2000+ max(@{$data{$contig}{$part}{end}});
            $i++;
            if( $data{$contig}{$part-1} || $data{$contig}{$part+1} ){
                foreach my $part2 ( sort keys %{$data{$contig}} ){
                    if ($part2 == $first && keys(%{$data{$contig}})> 1){next;}
                    my $from2= ($part2-1)*2000+ min(@{$data{$contig}{$part2}{start}});
                    my $to2= ($part2-1)*2000+ max(@{$data{$contig}{$part2}{end}});
                    next if $ya{$from}{$to} || $ya{$from}{$to2};
                    if( $from2 >= $from && $from2 <= $to ){
                        $fildata{$contig}{$i}{start}= $from;
                        $fildata{$contig}{$i}{strand}= $data{$contig}{$part}{strand};
                        $fildata{$contig}{$i}{oristart}= $part;
                        if( $to2 >= $from && $to2 <= $to ){
                            $fildata{$contig}{$i}{stop}= $to;
                            $fildata{$contig}{$i}{orito}= $part;
                            foreach my $partmotif ( @{$data{$contig}{$part}{motif}} ){
                                push @{$fildata{$contig}{$i}{motif}}, $partmotif;
                            }
                        }elsif( $to2 >= $from && $to2 >= $to ){
                            $fildata{$contig}{$i}{stop}= $to2;
                            $fildata{$contig}{$i}{orito}= $part2;
                            foreach my $partmotif ( @{$data{$contig}{$part}{motif}} ){
                                push @{$fildata{$contig}{$i}{motif}}, $partmotif;
                            }
                            foreach my $partmotif ( @{$data{$contig}{$part2}{motif}} ){
                                push @{$fildata{$contig}{$i}{motif}}, $partmotif;
                            }

                        }
                    $ya{$fildata{$contig}{$i}{start}}{$fildata{$contig}{$i}{stop}}++;
                    }
                }
            }else{
                $fildata{$contig}{$i}{start}= $from;
                $fildata{$contig}{$i}{strand}= $data{$contig}{$part}{strand};
                $fildata{$contig}{$i}{oristart}= $part;
                $fildata{$contig}{$i}{stop}= $to;
                $fildata{$contig}{$i}{orito}= $part;
                @{$fildata{$contig}{$i}{motif}}= @{$data{$contig}{$part}{motif}};
            }
        }
    }

    foreach my $filcont ( sort keys %fildata ){
        print RESUME "$filcont\n";
        foreach my $filpart ( sort keys %{$fildata{$filcont}} ){
            my $doubleend= 0;
            my $pos= 0;
            my $bad= 1;
            if( ($fildata{$filcont}{$filpart}{stop} - $fildata{$filcont}{$filpart}{start})>= 1000 && ($fildata{$filcont}{$filpart}{stop} - $fildata{$filcont}{$filpart}{start})<= 3000 && @{$fildata{$filcont}{$filpart}{motif}}> 3 ){
                foreach my $filpart_motif ( @{$fildata{$filcont}{$filpart}{motif}} ){
                    if( $fildata{$filcont}{$filpart}{strand} eq '+' ){
                        $doubleend++ if ( $filpart_motif ==  '3' || $filpart_motif ==  '4' );

                    }else{
                        $doubleend++ if ( $filpart_motif ==  '6' || $filpart_motif ==  '12' );
                    }
                }
                foreach my $filpart_motif ( @{$fildata{$filcont}{$filpart}{motif}} ){
                    $pos++;
                    if( $fildata{$filcont}{$filpart}{strand} eq '+' ){
                        $bad= 0 if ( ($filpart_motif ==  '3' || $filpart_motif ==  '4') && (@{$fildata{$filcont}{$filpart}{motif}} - $pos <= 4) );

                    }else{
                        $bad= 0 if ( ($filpart_motif ==  '6' || $filpart_motif ==  '12') && (@{$fildata{$filcont}{$filpart}{motif}} - $pos <= 4) );
                    }
                }
                print RESUME "$fildata{$filcont}{$filpart}{oristart}\t$fildata{$filcont}{$filpart}{orito}\t$fildata{$filcont}{$filpart}{start}\t$fildata{$filcont}{$filpart}{stop}\t$fildata{$filcont}{$filpart}{strand}\n" if ($doubleend <= 4 && $bad);
            }
        }
    }
}


sub extract_16S{
    open RIB, ">$org.rib_complete_tail.fna"; # aqui imprime los 16S que cumplen: >= 100b, <=2000b, al menos 3 motivos de los 21 que tiene que tener el 16S, tiene el motivo 3 en el que esta contenido el aSD
    open RIB2, ">$org.rib_nocomplete_tail.fna"; # aqui imprime los 16S que cumplen: >= 200b, <=2000b, al menos 3 motivos de los 21 que tiene que tener el 16S, no tiene el motivo 3 en el que esta contenido el aSD
    foreach my $filcont ( sort keys %fildata ){
        my $i= 0;
        foreach my $filpart ( sort keys %{$fildata{$filcont}} ){
            $i++;
            my $doubleend= 0;
            my $pos= 0;
            my $bad= 1;
            my $from= $fildata{$filcont}{$filpart}{start};
            my $to= $fildata{$filcont}{$filpart}{stop};
            my $strand= $fildata{$filcont}{$filpart}{strand};
            my $length16s= $to - $from;
            foreach my $filpart_motif ( @{$fildata{$filcont}{$filpart}{motif}} ){
                if( $strand eq '+' ){
                    $doubleend++ if ( $filpart_motif ==  '3' || $filpart_motif ==  '4' );

                }else{
                    $doubleend++ if ( $filpart_motif ==  '6' || $filpart_motif ==  '12' );
                }

            }
           foreach my $filpart_motif ( @{$fildata{$filcont}{$filpart}{motif}} ){
                $pos++;
                $bad= 0 if ($filpart_motif ==  '3');
            }

            $fildata{$filcont}{$filpart}{strand} eq "+" ? (print RIB ">${filcont}_$i\t$from..$to\t$strand\n", substr($seq->{$filcont}, $from, $length16s),"\n") : (print RIB ">${filcont}_$i\t$from..$to\t$strand\n", reverse_seq(substr($seq->{$filcont}, $from, $length16s)),"\n") if ($length16s >= 100 && $length16s <= 2000 && @{$fildata{$filcont}{$filpart}{motif}}> 10 && !$bad);

            $fildata{$filcont}{$filpart}{strand} eq "+" ? (print RIB2 ">${filcont}_$i\t$from..$to\t$strand\n", substr($seq->{$filcont}, $from, $length16s),"\n") : (print RIB2 ">${filcont}_$i\t$from..$to\t$strand\n", reverse_seq(substr($seq->{$filcont}, $from, $length16s)),"\n") if ($length16s >= 200 && $length16s <= 2000 && @{$fildata{$filcont}{$filpart}{motif}}> 10 && $bad);
        }
    }
}

sub print_structure{
    open STR, ">$org.str.txt";
    foreach my $contig ( sort keys %data ){
        print STR "$contig\n";
        foreach my $part ( sort keys %{$data{$contig}} ){
            print STR "\t$part\t$data{$contig}{$part}{strand}\n";
            foreach my $motif ( @{$data{$contig}{$part}{motif}} ){
                print STR "\t\t$motif";
            }
            print STR "\n";
            print STR "\t\t",($part-1)*2000+min(@{$data{$contig}{$part}{start}}), "\t", ($part-1)*2000+max(@{$data{$contig}{$part}{end}}), "\n";
            print STR "\n";
        }
    }
}

sub reverse_seq{
# reverse_seq($seq|\$seq) returns the complementary sequence of the DNA $seq
    my $nucs= ref($_[0]) ? scalar(reverse(${$_[0]} ) ) : scalar(reverse($_[0]));
    $nucs=~ tr/ACGTacgt/TGCAtgca/;
    return $nucs;
}

sub read_fasta{
	my ($name, %seq);
    open IN, $genome or die "Cant read $genome\n";
    while(<IN>){
    	chomp;
        if( /^>(\S+)/ ){
            $name= $1; next
        }
        $seq{$name}.= $_;
    }
    return \%seq;
}

__END__
Contig.#pos# == (#pos# - 1)*2000 # caso de ventana igual 4000
ex: NZ_LT601384.1.10 == 18000


Bacteria motifs order

+ 6* 16 20 12* 5 13 21 7 14 11 8 15 2 18 1 17 9 10 4* 19 3*

- 3* 19 4* 10 9 17 1 18 2 15 8 11 14 7 21 13 5 12* 20 16 6*

Archaeas motifs order

+ 9* 10* 20 13 6 18 8 16 11 5 14 1 15 4 17 7 12 3* 19 2*
