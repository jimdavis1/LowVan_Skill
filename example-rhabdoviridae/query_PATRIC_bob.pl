#! /usr/bin/env perl
use strict;
use gjoseqlib; 
use P3DataAPI;
use Data::Dumper;
use Getopt::Long;
use URI::Escape;
$| = 1;

my $usage = 'query_PATRIC.pl -c solr core -s sources -r return fields <list of valid identifiers
		
		Given a declared solr core, input type and a list of identifiers, this will return 
		The appropriately matching data from solr.
		
		So for instance:
		
		query_PATRIC.pl -c genome -i genome_id -s sources -r "product organism" <list of genome ids
		
		
		
        
        Options
		-c = solr core 
				e.g., 
				      sp_gene
				      genome 
				      etc.
				      
		-i input type:
		         default  = genome_id,
		         other valid types:
		         plfam_id
		         pgfam_id
		         sequence_id
		         etc.
		     
		
		-r = return fields
			 space-delimited text wrapped in double quotes
			 valid return fields include but are not limited to:
			 genome_id
			 patric_id
			 product
			 source
			 organism
			 query_coverage
			 subject_coverage
			 identity
			 e_value

	
		-s = sources, optional, for processing things like speciality proteins
		  (returns only enumerated sources)
			 space-delimited text wrapped in double quotes
			 valid sources include:
			     DrugBank
			     Victors
			     PATRIC_VF
			     VFDB
			     CARD
			     Human
			     TTD
			     ARDB
		
		-p = print header    			     	
		-h = help				
';

my $help;
my $input_type = "genome_id";
my ($sources, $return, $print_header, $core,); 
my $opts = GetOptions( 'h'   => \$help,
                       'c=s' => \$core,
                       'i=s' => \$input_type,
                       's=s' => \$sources,
                       'r=s' => \$return,
                       'p'   => \$print_header);
if ($help){die "$usage\n";}
my @sources;
my %sources;

if ($sources)
{
	@sources = split (" ", $sources);
	print Dumper %sources; 
	%sources = map{$_, 0}(@sources);
}

#unless ($return){die "must declare -r return fields\n";}
$return =~ s/ /,/g;
my @return = split (",", $return);

if ($print_header)
{
	print join ("\t", @return), "\n";
}


while (<>)
{	
	chomp;
	print STDERR "$_\n";
    my $q = $_;
	#my $api = P3DataAPI->new("http://cherry:3001");
    my $api = P3DataAPI->new();

	#$api->query_cb($core, \&handle_data, ["eq", $input_type, $q], ["select", $return]);
	$api->query_cb($core, \&handle_data, ["eq", $input_type, uri_escape($q)], ["select", $return]);
}


sub handle_data
{
    my($data) = @_;
    foreach my $ent (@$data)
    {
		for my $i (0..$#return)
		{		
			print "$ent->{$return[$i]}\t";
		}
	print "\n";
    }
	return 1;
}


















