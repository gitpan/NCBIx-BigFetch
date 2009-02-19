package NCBIx::BigFetch;
use warnings;
use strict;
use Class::Std;
use Class::Std::Utils;
use Carp;
use LWP::Simple;
use YAML qw(DumpFile LoadFile);
use Time::HiRes qw(usleep);
use Bio::SeqIO;

use version; our $VERSION = qv('0.5.5');

our $config_file   = 'efetch_N.yml';
our $esearch_file  = 'esearch_N.txt';
our $data_file     = 'sequences_N_M.txt';
our $clean_file    = 'sequences_N.txt';
our $sleep_policy  = 2_750_000;

# Silence Perl::Critic from complaining about $config hashref
## no critic (ProhibitAccessOfPrivateData)

{
	my %project_id_of         :ATTR( :get<asdf> :set<asdf> );
	my %base_url_of           :ATTR( :get<asdf> :set<asdf> );
	my %base_dir_of           :ATTR( :get<asdf> :set<asdf> );
	my %db_of                 :ATTR( :get<asdf> :set<asdf> );
	my %query_of              :ATTR( :get<asdf> :set<asdf> );
	my %querykey_of           :ATTR( :get<asdf> :set<asdf> );
	my %webenv_of             :ATTR( :get<asdf> :set<asdf> );
	my %count_of              :ATTR( :get<asdf> :set<asdf> );
	my %index_of              :ATTR( :get<asdf> :set<asdf> );
	my %start_date_of         :ATTR( :get<asdf> :set<asdf> );
	my %start_time_of         :ATTR( :get<asdf> :set<asdf> );
	my %return_max_of         :ATTR( :get<asdf> :set<asdf> );
	my %return_type_of        :ATTR( :get<asdf> :set<asdf> );
	my %missing_of            :ATTR( :get<asdf> :set<asdf> );

	sub next_index           { my ($self) = @_; my $ident = ident $self; $index_of{$ident} += $return_max_of{$ident}; $self->_save(); }

	sub get_config_filename  { my ($self)  = @_; my $project_id = $self->get_project_id(); $config_file =~ s/N/$project_id/; return $self->get_base_dir() . '/' . $config_file; } 
	sub get_esearch_filename { my ($self)  = @_; my $project_id = $self->get_project_id(); $esearch_file =~ s/N/$project_id/; return $self->get_base_dir() . '/' . $esearch_file; } 
	sub get_clean_filename   { my ($self)  = @_; my $project_id = $self->get_project_id(); $clean_file =~ s/N/$project_id/; return $self->get_base_dir() . '/' . $clean_file; } 
	sub get_data_filename    { my ($self, $index)  = @_; my $project_id = $self->get_project_id(); $index = defined($index) ? $index : $self->get_index(); my $filename = $data_file; $filename =~ s/N/$project_id/g; $filename =~ s/M/$index/g; return $self->get_base_dir() . '/' . $filename; } 

	sub BUILD {      
		my ($self, $ident, $arg_ref) = @_;

		# Set environment
		$config_of{$ident} = $self->_init( $arg_ref );

		# Check for existing project
		if (-e $self->get_config_filename()) { 
			$self->_status("Loading existing project");

			# Get existing config
			$config_of{$ident} = $self->_load();
		} else {
			$self->_status("Starting new project");

			# Set start date and time
			$config_of{$ident} = $self->_set_date( $config_of{$ident} );
	
			# Submit search and save config
			$config_of{$ident} = $self->_search( $config_of{$ident} );
			$self->_save();
		}

		return;
	}

	sub results_waiting {
		my ( $self )   = @_;
		if ( $self->get_index() < $self->get_count() ) { 
			return 1; 
		} else { 
			$self->_status("Found " . $self->_commify( scalar(@{ $self->get_missing() }) ) . " missing batches." );
			return 0; 
		}
	}

	sub missing_batches {
		my ( $self )   = @_;
		if ( @{ $self->get_missing() } ) { return 1; } else { return 0; }
	}

	sub get_next_batch {
		my ( $self )   = @_;
		my $index      = $self->get_index();

		# Get the batch
		$self->get_batch( $index );

		# Update the index
		$self->next_index();

		return;
	}

	sub get_batch {
		my ( $self, $index )   = @_;
		my $return_max = $self->get_return_max();
		my $return_type = $self->get_return_type();

		$self->_status("Starting with index " . $self->_commify( $index ) );
		
		# Ethics requires we wait sleep_policy microseconds before retrieving
		$self->_sleep();
		
		# Define a batch through URL
		my $efetch_url  = $self->get_base_url() . 'efetch.fcgi?db=' . $self->get_db();
		   $efetch_url .= '&WebEnv=' . $self->get_webenv() . '&query_key=' . $self->get_querykey() . "&rettype=$return_type";
		   $efetch_url .= "&retstart=$index&retmax=$return_max";
		   $efetch_url .= '&tool=ncbix_bigfetch&email=roger@iosea.com';

		# Get the batch using LWP::Simple (get)
		my $results  = get($efetch_url);
		
		# Check results # TODO: capture expired WebEnv and restart query
		if ( $results =~ m/resource is temporarily unavailable/i ) { $self->note_missing_batch( $index ); }
		if ( $results =~ m/NCBI C\+\+ Exception/i )                { $self->note_missing_batch( $index ); }
		if ( $results eq '' )                                      { $self->note_missing_batch( $index ); }

		# Save the sequences
		$self->_set_file_text( $self->get_data_filename( $index ), $results );

		return;
	}

	sub get_missing_batch {
		my ( $self )   = @_;

		# Get the next missing batch index
		my @missing    = @{ $self->get_missing() };
		my $index      = shift @missing;

		# Update the missing batch list
		$self->set_missing( \@missing );

		# Get the batch
		$self->get_batch( $index );

		return;
	}

	sub note_missing_batch { 
		my ( $self, $index )    = @_; 
		my @missing;
		my $missing             = $self->get_missing();
		if ( defined $missing ) { 
			@missing = @{ $missing }; 
		} else { 
			@missing = (); 
		}
		push @missing, $index;
		$self->set_missing( \@missing );
		$self->_save();
	}

	sub get_sequence {
		my ( $self, $id )   = @_;
		my $return_type = $self->get_return_type();

		$self->_status("Fetching sequence $id");
		
		# Ethics requires we wait sleep_policy microseconds before retrieving
		$self->_sleep();
		
		# Define a batch through URL
		my $efetch_url  = $self->get_base_url() . 'efetch.fcgi?db=' . $self->get_db();
		   $efetch_url .= '&id=' . $id . "&rettype=$return_type";
		   $efetch_url .= '&tool=ncbix_bigfetch&email=roger@iosea.com';

		# Get the sequence
		my $results  = get($efetch_url);
		
		# Save the sequences in missing file
		$self->_add_file_text( $self->get_data_filename( 0 ), $results );

		return;
	}

	sub unavailable_ids {
		my ( $self )     = @_;
		my $count        = $self->get_index();
		my $return_max   = $self->get_return_max();
		my $index        = 1;
		my @unavailables = ();

		while ( $index < $count ) {
			$self->_status("Checking " . $self->_commify( $index ) . " through " . $self->_commify( $index + $return_max - 1 ) );
		
			# Get the sequences
			my $text = $self->_get_file_text( $self->get_data_filename( $index ) );

			while ( $text =~ m/Error:\s(\d+)\sis\snot\savailable\sat\sthis\stime/g ) { push @unavailables, $1; }

			# Update the index
			$index += $return_max;
		}

		$self->_status("Found " . $self->_commify( scalar(@unavailables) ) . " unavailable ids." );

		return \@unavailables;
	}

	sub clean_sequences {
		my ( $self )     = @_;
		my $count        = $self->get_index();
		my $return_max   = $self->get_return_max();
		my $index        = 1;
		my $clean_file   = $self->get_clean_filename();
		my $out          = Bio::SeqIO->new(-file => ">$clean_file", -format => 'Fasta');

		while ( $index < $count ) {
			my @good_sequences = $self->_check_ambiguous( $index );
			foreach my $seq (@good_sequences) { $out->write_seq($seq); }

			# Update the index
			$index += $return_max;
		}

		my @good_sequences = $self->_check_ambiguous( 0 );
		foreach my $seq (@good_sequences) { $out->write_seq($seq); }

		$self->_status("Clean complete.");

		return;
	}

	sub _check_ambiguous {
		my ( $self, $index ) = @_;
		my @good_sequences   = ();

		my $in_file = $self->get_data_filename( $index );

		if (-e $in_file) {
			$self->_status("Cleaning $in_file");
			
			# Get the sequences
			my $in  = Bio::SeqIO->new(-file => "<$in_file",  -format => 'Fasta');
	
			# Process each sequence from the input file (create $seq object)
			while ( my $seq = $in->next_seq() ) {
			        # Get the sequence
				my $sequence = $seq->seq();
				
				# Filter for ambiguous characters
				my $good_sequence = 1;
				if ($sequence =~ m/J/)     { $good_sequence = 0; } 
				if ($sequence =~ m/Error/) { $good_sequence = 0; } 
					
				if ($good_sequence) { push @good_sequences, $seq; }
			}
		}

		return @good_sequences;
	}

	sub _init {
		my ( $self, $arg_ref ) = @_;

		$arg_ref->{base_url}     = $arg_ref->{base_url}   ? $arg_ref->{base_url}   : "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/";    
		$arg_ref->{base_dir}     = $arg_ref->{base_dir}   ? $arg_ref->{base_dir}   : $self->_get_base_dir();    
		$arg_ref->{query}        = $arg_ref->{query}      ? $arg_ref->{query}      : "";    
		$arg_ref->{db}           = $arg_ref->{db}         ? $arg_ref->{db}         : "protein";    
		$arg_ref->{return_max}   = $arg_ref->{return_max} ? $arg_ref->{return_max} : "500";    
		$arg_ref->{index}        = $arg_ref->{index}      ? $arg_ref->{index}      : "1";    
		$arg_ref->{project_id}   = $arg_ref->{project_id} ? $arg_ref->{project_id} : "1";    
		$arg_ref->{missing}      = $arg_ref->{missing}    ? $arg_ref->{missing}    : [];
		$arg_ref->{return_type}  = $arg_ref->{return_type} ? $arg_ref->{return_type} : "fasta";	
		return $arg_ref;
	}

	sub _set_date {
		my ( $self, $arg_ref ) = @_;

		my @time  = localtime;
		my $year  = 1900 + $time[5];
		my $month = $time[4] + 1; $month =~ s/^(\d)$/0$1/;
		my $day   = $time[3];     $day   =~ s/^(\d)$/0$1/;
		my $hour  = $time[2];     $hour  =~ s/^(\d)$/0$1/;
		my $min   = $time[1];     $min   =~ s/^(\d)$/0$1/;
		my $sec   = $time[0];     $sec   =~ s/^(\d)$/0$1/;
		
		$arg_ref->{start_date} = "$year-$month-$day";
		$arg_ref->{start_time} = "$hour:$min:$sec";

		return $arg_ref;
	}

	sub _search {
		my ( $self, $arg_ref ) = @_;

		# Get search result ticket
		my $esearch_url          = $self->get_base_url() . 'esearch.fcgi?db=' . $self->get_db();
		   $esearch_url         .= '&term=' . $self->get_query() . '&usehistory=y';
		   $esearch_url         .= '&tool=ncbix_bigfetch&email=roger@iosea.com';
		my $esearch_result       = get($esearch_url);

		# Save search result
		$self->_set_file_text( $self->get_esearch_filename(), $esearch_result );
		
		# Parse the relevant keys
		$esearch_result =~ m/<Count>([0-9]*)<\/Count>/g;               $arg_ref->{count}    = $1;
		$esearch_result =~ m/<QueryKey>([0-9]*)<\/QueryKey>/g;         $arg_ref->{querykey} = $1;
		$esearch_result =~ m/<WebEnv>([\.a-zA-Z0-9_@\-]*)<\/WebEnv>/g; $arg_ref->{webenv}   = $1;

		return $arg_ref;
	}

	sub _load {
		my ( $self ) = @_;
		my %config = %{ LoadFile( $self->get_config_filename() ) };
		$self->set_project_id( $config{project_id} );
		$self->set_base_url( $config{base_url} );
		$self->set_base_dir( $config{base_dir} );
		$self->set_db( $config{db} );
		$self->set_query( $config{query} );
		$self->set_querykey( $config{querykey} );
		$self->set_webenv( $config{webenv} );
		$self->set_count( $config{count} );
		$self->set_index( $config{index} );
		$self->set_start_date( $config{start_date} );
		$self->set_start_time( $config{start_time} );
		$self->set_return_max( $config{return_max} );
		$self->set_return_type( $config{return_type} );
		$self->set_missing( $config{missing} );
	}

	sub _save {
		my ( $self ) = @_;
		my $config   = {  project_id   => $project_id_of,
				  base_url     => $base_url_of,
				  base_dir     => $base_dir_of,
				  db           => $db_of,
				  query        => $query_of,
				  querykey     => $querykey_of,
				  webenv       => $webenv_of,
				  count        => $count_of,
				  index        => $index_of,
				  start_date   => $start_date_of,
				  start_time   => $start_time_of,
				  return_max   => $return_max_of,
				  return_type  => $return_type_of,
				  missing      => $missing_of };
		DumpFile( $self->get_config_filename(), $config );
		return;
	}

	sub _get_base_dir {
		my ( $self, $base_dir ) = @_;
		chomp( my $id = `id -nu`);
		if ($id eq 'root') { $base_dir = '/root'; } else { $base_dir = '/home/' . $id; }
		return $base_dir;
	}

	sub _status {
		my ( $self, $msg ) = @_;
		print STDOUT "  STATUS: $msg \n";
		return;
	}

	sub _sleep {
		my ( $self ) = @_;
		usleep($sleep_policy);
		return;
	}

	sub _get_file_text {
		my ( $self, $path_file_name ) = @_;
		my ($text, $line);
		if (-e $path_file_name) {
			open  (my $IN, '<', $path_file_name) || croak( "Cannot open $path_file_name: $!" );
			while ($line = <$IN>) { $text .= $line; }
			close ($IN)                          || croak( "Cannot close $path_file_name: $!" );
		}
		return $text;
	}
	
	sub _set_file_text {
		my ( $self, $path_file_name, $text ) = @_;
		open  (my $OUT, '>', $path_file_name)        || croak( "Cannot open $path_file_name: $!" );
		print $OUT $text                             || croak( "Cannot write $path_file_name: $!" );
		close ($OUT)                                 || croak( "Cannot close $path_file_name: $!" );
	}
	
	sub _add_file_text {
		my ( $self, $path_file_name, $text ) = @_;
		open  (my $OUT, '>>', $path_file_name)       || croak( "Cannot open $path_file_name: $!" );
		print $OUT $text                             || croak( "Cannot write $path_file_name: $!" );
		close ($OUT)                                 || croak( "Cannot close $path_file_name: $!" );
	}

	sub _commify { # Perl Cookbook 2.17
		my ( $self, $string ) = @_;
		my $text = reverse $string;
		$text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
		return scalar reverse $text;
	}
	
	sub authors { return 'Roger Hall <roger@iosea.com>, Michael Bauer <mbkodos@gmail.org>, Kamakshi Duvvuru <kduvvuru@gmail.com>'; }
}

1; # Magic true value required at end of module
__END__

=head1 NAME

NCBIx::BigFetch - Retrieve very large NCBI sequence result sets based on keyword search

=head1 SYNOPSIS

  use NCBIx::BigFetch;
  
  # Parameters
  my $params = { project_id => "1", 
                 base_dir   => "/home/user/data", 
  	         db         => "protein",
  	         query      => "apoptosis",
                 return_max => "500" };
  
  # Start project
  my $project = NCBIx::BigFetch->new( $params );
  
  # Love the one you're with
  print " AUTHORS: " . $project->authors() . "\n";
  
  # Attempt all batches of sequences
  while ( $project->results_waiting() ) { $project->get_next_batch(); }
  
  # Get missing batches 
  while ( $project->missing_batches() ) { $project->get_missing_batch(); }
  
  # Find unavailable ids
  my $ids = $project->unavailable_ids();
  
  # Retrieve unavailable ids
  foreach my $id ( @$ids ) { $project->get_sequence( $id ); }

=head1 DESCRIPTION

NCBIx::BigFetch uses the esearch and efetch services of NCBI 
to retrieve sequences by keyword. It was designed for very large 
result sets; the first project it was used on had over 11,000,000 
sequences. 

Downloaded data is organized by "project id" and "base directory" 
and saved in text files. Each file includes the project id in 
its name. Besides the data files, two other files are saved: 
1) the initial search result, which includes the WebEnv key, and 
2) a configuration file, which saves the parsed data and is used 
to pick-up the download and recover missing batches or sequences. 

Results are retrived in batches depending on the "retmax" size. 

=head2 METHODS

=over 4

=item * new()

  my $project = NCBIx::BigFetch->new( $params );

The parameters hash reference should include the following minimum 
keys: project_id, base_dir, db, and query. 

=item * results_waiting()

  while ( $project->results_waiting() ) { ... }

This method is used to determine if all of the batches have been 
attempted. It compares the current index to the total count, and 
is TRUE if the index is less than the count.

=item * get_next_batch()

  $project->get_next_batch();

Attempts to retrieve the next batch of "retmax" sequences, starting 
with the current index, which is updated every time a batch is 
downloaded. When used as in the Synopsis above, the index is both 
kept in memory and updated in the configuration file. If the 
download is interrupted and restarted, the correct index will be 
used and no data will be lost.

=item * note_missing_batch()

  $project->note_missing_batch( $index );

Adds the batch index to the list of missing batches.

=item * missing_batches()

  while ( $project->missing_batches() ) { ... }

This method is used to determine if any batches have been noted 
as "missing". It measures the "missing" list (which is stored 
in the configuration file) and returns TRUE when at leat one batch 
is listed. The batches are listed by starting index, which 
together with the return_max setting is used to describe a batch.

=item * get_missing_batch()

  $project->get_missing_batch();

Warning: do not kill the script during this phase. 

Gets a single batch, using the first index on the "missing" list. 
The index is shifted off the list and then attempted, so if you 
break during this phase you may actually lose track of the batch.

Recovery: edit the configuration file and add the index back to the 
missing list. The index will be reported to STDOUT in the status 
message.

=item * get_batch()

  $project->get_batch( $index );

Gets a single batch using the index parameter. This routine may be 
called on its own, but it is intended to only be used by get_next_batch() 
and get_missing_batch().

=item * unavailable_ids()

  my $ids = $project->unavailable_ids();

Notice that this method depends on a loaded (or started) project. It 
reads through all data files and creates a list of individual 
sequences that were unavailable when a batch was reported. The list 
is returned as a perl list reference.

=item * get_sequence()

  $project->get_sequence( $id );

Notice that this method depends on a loaded (or started) project. It 
retrieves the sequence by id and saves it to a special data file 
which uses "0" as an index. All unavailable sequences retrieved 
this way are saved to this file, so it could potentially be larger 
than the rest.

=item * clean_sequences()

  $project->clean_sequences();

Removes non-sequence text from sequence files and optionally removes 
sequences with ambiguous characters.

=item * authors()

  $project->authors();

Surely you can stand a few bytes of vanity for the price of free software!

=item * BUILD()

  my $project = NCBIx::BigFetch->new( $params );

This method is *not* called directly, but rather included in the new() 
method thanks to Class::Std.

=back

=head2 PROPERTIES

These get/set functions manage the modules properties.

=over 4

=item * get_base_dir()

  $project->get_base_dir();

Gets the base directory for project data.

=item * get_base_url()

  $project->get_base_url();

Gets the base URL for NCBI eUtils.

=item * get_clean_filename()

  $project->get_clean_filename();

Creates a filename for a project's sequences to store "cleaned" sequences.

=item * get_config_filename()

  $project->get_config_filename();

Creates a filename for the configuration file based on the project_id.

=item * get_count()

  $project->get_count();

Returns the count of results for the query.

=item * get_data_filename()

  $project->get_data_filename();

Creates a filename for a given batch based on project_id and result index.

=item * get_db()

  $project->get_db();

Gets the eSearch database setting.

=item * get_esearch_filename()

  $project->get_esearch_filename();

Creates a filename for saving the intial search request.

=item * get_index()

  $project->get_index();

Gets the current result index. The index is reset after every attempted batch by 
retmax amount.

=item * get_missing()

  $project->get_missing();

Gets the list of missing batch indices.

=item * get_project_id()

  $project->get_project_id();

Gets the project_id for the loaded project.

=item * get_query()

  $project->get_query();

Gets the query string used for eSearch.

=item * get_querykey()

  $project->get_querykey();

Gets the querykey setting from the eSearch results.

=item * get_return_max()

  $project->get_return_max();

Gets the retmax setting used to limit the batch size.

=item * get_return_type()

  $project->get_return_type();

Gets the rettype setting used to determine the format 
of fetched sequences.

=item * get_start_date()

  $project->get_start_date();

Calculates the start date for the project.

=item * get_start_time()

  $project->get_start_time();

Calculates the start time for the project.

=item * get_webenv()

  $project->get_webenv();

Gets the WebEnv key returned from eSearch. It is used 
to build the eFetch URL for retrieving batches of 
sequences.

=item * next_index()

  $project->next_index();

Gets the next result index, which defines the batch id.

=item * set_index()

  $project->set_index();

Sets the result index.

=item * set_missing()

  $project->set_missing();

Sets the list of missing batches.

=back

=head2 EXPORT

None

=head1 SEE ALSO

http://bioinformatics.ualr.edu/

http://www.ncbi.nlm.nih.gov/entrez/query/static/efetch_help.html

http://eutils.ncbi.nlm.nih.gov/entrez/query/static/efetchseq_help.html

http://www.ncbi.nlm.nih.gov/entrez/query/static/eutils_example.pl

=head1 AUTHORS

Roger Hall (roger@iosea.com), (rahall2@ualr.edu)

Michael Bauer (mbkodos@gmail.com), (mabauer@ualr.edu) 

Kamakshi Duvvuru (kduvvuru@gmail.com) 

=head1 COPYRIGHT AND LICENSE

Copyleft (C) 2009 by the Authors

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.5 or,
at your option, any later version of Perl 5 you may have available.

=cut
