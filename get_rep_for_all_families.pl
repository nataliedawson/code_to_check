#!/usr/bin/env perl

# Script to choose a representative domain sequence for all functional families
# Families with structural domains - find domain with highest cumulative SSAP score
# Families with no structural dom  - find domain with highest no. of deepest GO (F) terms - TO DO

# note that 'family' always means 'functional family'

# PA: I couldn't run the program successfully, due to a problem with the SAP score module.
# The error I got was:
# "Command 'SSAP' not found in /opt/local/perls/build-trunk/bin, /usr/lib64/qt-3.3/bin, 
# /usr/lib64/mpich/bin, /usr/lib64/ccache, /usr/local/bin, /usr/bin, /usr/local/sbin, 
# /usr/sbin at /usr/local/svn/trunk/lib/perl5/Cath/Role/Cmd.pm line 37."

use strict;
use warnings;

use Path::Class;

use File::Basename;

use Cath::SsapList::Matrix;
use Cath::SsapList::File;
#use Cath::Ssap;
use Getopt::Std;
use File::stat;

#OPTIONS######################
my %arg;
getopts( 'hi:o:d:v:g:m:', \%arg );
help() if ( defined $arg{'h'} );
die "Pass directory path to functional family FASTA file with -i\n" unless ( defined $arg{'i'} );
die "Pass directory path to store family representative with -o\n" unless ( defined $arg{'o'} );
die "Pass directory path to SSAP data cache with -d\n" unless ( defined $arg{'d'} );
die "Pass CATH version with -v\n" unless ( defined $arg{'v'} );
#die "Pass parent directory path to family GO terms found with -g\n" unless ( defined $arg{'g'} );
#die "Pass parent directory path to superfamily GO to MD5 mappings with -m\n" unless ( defined $arg{'m'} );

#Open the family FASTA file and extract all (if any) structural domains ids#####
my $family_fasta = file($arg{'i'})->absolute;
my $ssap_directory_path = dir( $arg{'d'} )->absolute;

my $superfamily = basename( $family_fasta->dir );
my $family_id = basename($family_fasta, ".faa"); #extract family id

# check whether rep file already exists
my $rep_file = file("$ssap_directory_path/$family_id.rep"); # define rep file path
if ( -e $rep_file ) {
    my $filestat = stat($rep_file);
    my $filesize = $filestat->size;
    print "Filesize = $filesize\n";
    if ( $filesize > 0 ) {
        print "Found non-empty rep file. Skipping family...\n";
        exit 0; # exit safely
    }
}

info("Processing functional family ($family_id) belonging to superfamily: $superfamily\n");
my $family_fh = $family_fasta->openr or die "Cannot open $family_fasta for reading: $!\n";
info(" Opened FASTA file: $family_fasta\n");


# get sequence header format
my $cath_version = $arg{'v'};

# get parent dir of all family GO term files
#my $go_term_dir = dir( $arg{'g'} );

# get parent dir of all superfamily GO to MD5 mappings
#my $go_md5_dir = dir( $arg{'m'} );

# define array to hold non-redundant list of all domains in family 
my @domains;
my @md5s;

while ( <$family_fh> ) {
 
    chomp $_;
    
    # check cath version and extract CATH structural domains accordingly
    if ( $cath_version eq '4.0' ) {
        if ($_ =~ /^>cath/) {
            my $domain = substr $_, 12, 7; # extract domain from header
            push (@domains, $domain);      # push onto array
        }
        elsif ( $_ =~ /^>\w{32}/ ) {
            my $md5 = substr $_, 1, 32;    # extract MD5 from header
            push( @md5s, $md5 );           # push onto array
        }
    }
    else {
        if ($_ =~ /^>domain/) {
            my $domain = substr $_, 8, 7; # extract domain from header
            push (@domains, $domain);     # push onto array
        }
        elsif ( $_ =~ /^>\w{32}/ ) {
            my $md5 = substr $_, 1, 32;    # extract MD5 from header
            push( @md5s, $md5 );           # push onto array
        }
    }

}

#print "\@domain ids: @domains\n";

# create output directory path to hold reps
my ( $out_path ) = dir( $arg{'o'} )->absolute;
if ( !-e $out_path ) {
	info( "output path '$out_path' does not exist, trying to create..." );
	$out_path->mkpath;
}

# count the number of domains/md5s found in the functional family
my $ssap_count_all = scalar @domains;
my $md5_count_all = scalar @md5s;

info(" Found $ssap_count_all CATH structural domains\n");
info(" Found $md5_count_all MD5 sequences\n");



# choose which method of assignment to use based upon the number of CATH domains
if ( $ssap_count_all == 1) {
    info(" Found one CATH structural domain. Defining @domains as the representative.\n");
    # assign rep
    my $rep_file = file( "$out_path/${family_id}.rep" );
    my $rep_fh = $rep_file->openw or die "Cannot open $rep_file for writing; $!";
    print $rep_fh @domains; # print the only domain present out as rep
} elsif ( $ssap_count_all == 0 ) {
    info(" Found zero structural domain sequences in this family. Using GO method...\n");
    info("TEMP: SKIPPING GO METHOD THIS TIME.\n");
    # map all GO terms in superfamily to their MD5 sequence IDs
    #my $go_md5_map = map_go_to_md5( $go_md5_dir, $superfamily, \@md5s );
    
    # GO method
    #run_GO_method( $out_path, $family_id, $go_term_dir, $go_md5_dir, $superfamily, $go_md5_map, \@md5s ); # to finish
} elsif ( $ssap_count_all > 1 ) {
    info(" At least two CATH structural domains found. Using SSAP method...\n");
    # SSAP method
    run_ssap_method( $ssap_directory_path, $out_path, $family_id, $cath_version, \@domains );
}




####################SUBROUTINES#################################################


sub help {
	print <<"END_MESSAGE";

  Usage: $0 [Options]

  Options:

    -h                   Print this message.

    -i string            Input path to family FASTA file.
                         Required.

    -o /path/to/dir      Path to output directory to hold the 
                         family representative.
                         Required.
                         
    -d /path/to/dir      Path to output directory to store the SSAP output.
                         Required.
                         
    -v string            CATH version. E.g. 3.5.
                         Required.
                         
    -g /path/to/dir      Path to parent directory of family GO terms identified.
                         Not in use.
                         
    -m /path/to/dir      Path to parent directory of GO to MD5 mappings.
                         Not in use.
                         
  
  Chooses a representative domain sequence for each functional family supplied.
  
  Uses structural information where possible to calculate the most structurally 
  similar domain in a family. Otherwise finds the domain with the highest 
  number of "deepest GO terms".
  
  Outputs the family representative, the SSAP files, and the list of family 
  structural IDs with the family representative at the top of the list.

END_MESSAGE

	exit;
}

################################################################################

sub create_output_directory_path {
    
    my ( $out_path ) = dir( shift ); # will this syntax work?!
     
    # Create the output directory path if it does not already exist
    my $out_dir_parent_path = $out_path->parent;
    
    if ( ! -e $out_dir_parent_path ) {
        info(" $out_dir_parent_path does not exist, creating folder...\n");
        mkdir( $out_dir_parent_path ) or die "Cannot create $out_dir_parent_path";
    }
    else {
        info(" Cannot create folder $out_dir_parent_path, it already exists...\n");
    }
    
    if ( ! -e $out_path )
    {
        # change to get ../ folder so can mkdir on that first
        info(" $out_path does not exist, creating folder...\n");
        mkdir( $out_path ) or die "Cannot create $out_path";
    }
    else {
        info(" Cannot create folder for $out_path, it already exists...\n");
    }
    return( $out_path );
} 

################################################################################

sub run_ssap_method {
 
    my ( $ssap_directory_path, $out_path, $family_id, $cath_version, $domains ) = @_;
#    print "$ssap_directory_path, $out_path, $family_id, $cath_version\n";

    info( " Changing directory to $ssap_directory_path" ); # change to ssap dir
    chdir "$ssap_directory_path"
    	or die "Error: failed to chdir to '$ssap_directory_path': $!";
    	
    my $rep_file = file("$out_path/$family_id.rep"); # define rep file path
    
    # create file to store SSAP output
    my $ssap_out_file = file("$out_path/$family_id.ssap");
#    my $ssap_out_fh = $ssap_out_file->openr or die "Cannot open $ssap_out_file for reading: $!";

    my $ssaplist = ensure_complete_ssap_file( "$ssap_out_file", $domains, $cath_version );
    
    #get all ssap entries in ssaplist
    my @ssaps = $ssaplist->all_ssaps;

    info( "After write file..." );
    
    #write the file out
#    $ssaplist->to_file( $ssap_out_file );
    
    my $max_ssap_score = 0;
    my $max_ssap_query;
    my %results; # define hash to total up SSAP scores
    my $r_results; # reference to results hash
    
    my %accumulated_ssap_score_by_domain_id;
    
    foreach my $ssap (@ssaps) {
        $accumulated_ssap_score_by_domain_id{ $ssap->query_id } += $ssap->ssap_score;
        $accumulated_ssap_score_by_domain_id{ $ssap->match_id } += $ssap->ssap_score;        
    }
    
    
    
    # calculate the highest cumulative SSAP score and the corresponding domain
    my ($max_score, $best_query) = &get_best_ssap_total(\%accumulated_ssap_score_by_domain_id);
    
    info( " max SSAP score: $max_score\tbest query: $best_query\n");

    my $rep_file_fh = $rep_file->openw or die "Cannot open $rep_file for writing: $!";
    
    print $rep_file_fh "$best_query\n"; #print out the determined rep to file
    
    # print out list of structural domain ids with the family representative id at the top
    # this is a file needed later on for the mapping of CSA residues
    my $domain_ids_file = file("$out_path/$family_id.ids");
    
    my $domain_ids_fh = $domain_ids_file->openw or die "Cannot open file for writing: $!";
    
    print $domain_ids_fh "$best_query\n";
    
    foreach my $domain (@domains)
    {
    	next if ($domain eq $best_query);
    	print $domain_ids_fh "$domain\n";
    }
}

################################################################################

sub ensure_complete_ssap_file {
 my ($ssap_file, $domains, $cath_version) = @_;

 my $ssaplist;
 if ( -e $ssap_file ) { # check the SSAP matrix in a file
  $ssaplist = Cath::SsapList::File->new( version => $cath_version, file => $ssap_file ); 
 }
 else { # check SSAP matrix of all family domains
  $ssaplist = Cath::SsapList::Matrix->new( domain_ids => $domains, version => $cath_version );
 }

 info( " Checking SSAP matrix for completeness of file $ssap_file...\n");

 my $count_ssaps_before = $ssaplist->count_ssaps;
 $ssaplist->check_matrix( domains => $domains, on_error => 'RUN_SSAPS' );
# $ssaplist->check_matrix( domains => $domains );
 my $count_ssaps_after = $ssaplist->count_ssaps;
 
 info( " $count_ssaps_before SSAPs in original. Should be $count_ssaps_after.\n");
 
 if ( ! -e $ssap_file || $count_ssaps_before != $count_ssaps_after ){
  info( "Writing SSAP file: $ssap_file" );
#  my $ssap_fh = $ssap_file->openw or die "Cannot open $ssap_file for writing: $!";
#  $ssaplist->run_ssaps;
  $ssaplist->to_file( $ssap_file );
 }

 return $ssaplist;
}


################################################################################

sub info {
	my $msg = "@_";
	chomp($msg);
	print "[$$] $msg\n";
}

################################################################################

#sub accumulate_ssap_score {
#
#    # Cath::Ssap::Entry
#	my ( $ssap, $results ) = @_;
#
#	#$results{$query}++;
#
#	$results->{$ssap->query_id}{'total_ssap_score'} += $ssap->ssap_score; #accumulates the ssap score for each domain query
#	$results->{$ssap->match_id}{'total_ssap_score'} += $ssap->ssap_score; #accumulates the ssap score for each domain *match*
#
##    print "query: ", $ssap->query_id, " ", $results->{$ssap->query_id}{'total_ssap_score'},  " += " , $ssap->ssap_score, "\n";
##    print "match: ", $ssap->match_id, " ", $results->{$ssap->match_id}{'total_ssap_score'},  " += " , $ssap->ssap_score, "\n";
#
#	return ($results); 
#}

################################################################################
sub get_best_ssap_total {

	my $r_results = shift;

	my $max_score = 0;
	my $best_score = 0;
	my $best_query;

	for my $domain_id (sort keys %$r_results)
	{
	    my $score = $r_results->{$domain_id};
	    if ( $score > $best_score ) {
	         $best_score = $score;
	         $best_query = $domain_id;
	    }
	}
	
#	my @sorted_domain_ids = sort { $results->{$a} <=> $results->{$b} } keys %$results;
	
#	my $top_domain_id = $sorted_domain_ids[0];

	
	return ($best_score, $best_query);
}

################################################################################

sub map_go_to_md5 { # not currently used
 
 # mapping family MD5 sequences to GO terms. 
 # N.B. Mapping not present for all sequences. 
    
    my ( $go_md5_dir, $superfamily, $md5s ) = @_;
    my %go_md5; # hash to store mappings for the superfamily
    
    my $superfamily_go_md5_file = file( "$go_md5_dir/$superfamily.anno" );
    my $superfamily_go_md5_fh = $superfamily_go_md5_file->openr 
      or die "Cannot open $superfamily_go_md5_file for reading: $!";
    
    print " Opened $superfamily_go_md5_file\n";
    
    my @lines = <$superfamily_go_md5_fh>;
    
    # store mapping of MD5 ID to its GO term as long the MD5 is present in family
    foreach my $line ( @lines ) {
        chomp $line;
        
        my ( $md5, $go_term_list, $ec ) = split( /\t/, $line );
        
#        print " Searching for $md5 in @$md5s\n";
        if ( grep( /$md5/, @$md5s) ) {
#            print " Found $md5 in family MD5 list\n";
            
            # split go terms variable into each go term
            my @go_terms = split( /;/, $go_term_list );
            
            # push on GO term list for each family sequence
            push( @{ $go_md5{$md5} }, @go_terms );
        }

    }
   return( \%go_md5 ); # return mappings of family MD5s to GO terms
}


################################################################################

sub run_GO_method { # not currently used
    # get rep for functional families with no CATH structural domains
    
    my ( $out_path, $family_id, $go_term_dir, $go_md5_dir, $superfamily, $go_md5_map, $md5s ) = @_; # output dir and family id
    
    # array to store GO terms
    my @go_terms;
    
    # construct filename for file of predicted family GO terms (and open file)
    my $family_go_file = file( "$go_term_dir/$superfamily/${family_id}.anno" );
    my $family_go_fh = $family_go_file->openr or die "Cannot open $family_go_file for reading: $!";
    my @family_go_lines = <$family_go_fh>;
    
    # get list of most specific GO terms in that family
    my ( $go_terms ) = get_most_specific_GO_terms( \@family_go_lines );
    
    print " List of most-specific GO terms for $family_id: @$go_terms\n";
    
    # count the number of most-specific GO terms mapped to each MD5
    my %md5_go_count;
    
    # count the number of MD5s with the most specific GO terms
    foreach my $go_term ( @$go_terms ) { # each specific GO term
        foreach my $md5 ( @$md5s ) {       # each family MD5
            if ( grep( /$go_term/, @{ $go_md5_map->{$md5}} ) ) { # find if GO term mapped to MD5
                print " Most specific GO term $go_term maps to $md5\n";
                $md5_go_count{$md5}{'COUNT'}++;
            }
        }
    }
    
    my $max_score = 0;
    my $max_md5;
    
    # determine the MD5 with the highest number of most-specific GO terms
    foreach my $md5 ( sort keys %md5_go_count ) {
        if ( $md5_go_count{$md5}{'COUNT'} > $max_score ) {
            print " Specific GO count for $md5 is higher than max score $max_score\n",
                  " Updating max score to $md5_go_count{$md5}{'COUNT'}\n";
            $max_score = $md5_go_count{$md5}{'COUNT'}; # updating max count
            $max_md5 = $md5;                           # updating MD5 ID
        }
    } 
    
    my $rep = $max_md5;
    if ( ! defined $rep ) {
        print "DEBUG: Representative (MD5) not defined for $superfamily $family_id.\n";
    }
    print " Choosing $rep as the representative\n";

    my $rep_file = file( "$out_path/${family_id}.rep");
    my $rep_fh = $rep_file->openw or die "Cannot open $rep_file for writing: $!";
    
    print $rep_fh "$rep"; # print rep to file
}

################################################################################

sub get_most_specific_GO_terms { # not currently used
    
    my ( $family_go_lines ) = @_;
    my @go_terms; # to hold list of most-specific GO terms
    
    # go through each line and find the most specific GO terms (col 7 == yes)
    foreach my $go_line ( @$family_go_lines ) {
        chomp $go_line;
        
        my @parts = split( /\t/, $go_line );
        
        my $go_term = $parts[0];             # GO term
        my $ontology_category = $parts[1];   # GO category (F, P or C)
        my $is_most_specific = $parts[6];    # yes / no
        
        # find lines with most specific GO term from the F category
        if ( ( $is_most_specific eq "yes" ) && ( $ontology_category eq "F" ) ) { 
            # extract GO term
            push( @go_terms, $go_term );  
        }
    }
    # return list of most-specific GO terms for the functional family
    return( \@go_terms );
}
