package AminoAcidMatrix;

#===============================================================================
#  Routines to efficiently provide position-by-position scoring of alignments.
#  BLOSUM and PAM matrices are built-in, and can be accessed by name.  It is
#  also possible to initialize with a user-supplied matrix in the PAM matrix
#  style, or as a hash of values, or an array of values.
#
#  Amino acids J and U are supported.
#  Symbols X, * and - are supported.
#  All other characters are converted to X.
#  Afine gap penalties are not supported, they do not localize to specific sites.
#-------------------------------------------------------------------------------
#  Creation:
#
#    $MatrixObj = new(                  @Options )   #  Default matrix is BLOSUM62
#    $MatrixObj = new(           $name, @Options )   #  Look up existing matrix by name
#    $MatrixObj = new( NAME  =>  $name, @Options )   #  Look up existing matrix by name
#    $MatrixObj = new( FILE  =>  $file, @Options )   #  Read from file
#    $MatrixObj = new( FILE  => \*FH,   @Options )   #  Read from open file handle 
#    $MatrixObj = new( FILE  =>  '',    @Options )   #  Read from STDIN 
#    $MatrixObj = new( FILE  => \$text, @Options )   #  Read from reference to string
#    $MatrixObj = new( TEXT  =>  $text, @Options )   #  Read from text string
#    $MatrixObj = new( TEXT  => \$text, @Options )   #  Read from reference to string
#    $MatrixObj = new( ARRAY => [ \@residues, \@scr_rows ], @options )
#                                              #  List of residue names and array of scores
#    $MatrixObj = new( HASH  => { symb_pair => $scr, ... },          @options )
#    $MatrixObj = new( HASH  => { symb1 => { symb2 => $scr, ... } }, @options )
#
#  Creation options:
#
#    COMMENTS => \@comments    #  Add comment lines to the matrix
#    GAP      =>  $int         #  Raw score for a gap; overwrites existing values
#    NAME     =>  $name        #  Name the supplied matrix, or lookup matrix by name
#    SCALE    =>  $float       #  Bits per unit score (usually 1/2, 1/3, 1/4 or 1/5);
#                              #      overwrites any value divined from other data
#
#  Creation notes:
#
#      By far, the most common use is expected to be requesting a matrix by
#      name (e.g. BLOSUM45 or PAM120).  The numerical part of the name will
#      be matched as closely as possible from the available versions.  You can
#      see what you got with the name method.  Including a request for a given
#      scale factor (e.g., BLOSUM45_5, for 1/5 bit scoring units) is rarely
#      useful, but attempts will be made to match the request.  A matrix is
#      never rescaled to the requested value.  Get a list of available matrices
#      with the function AminoAcidMatrix::matrix_list().
#
#      When a matrix object is created, the values might be rescaled to either
#      compress the range of scores to -127 <= score <= 127, or to expand the
#      range of scores to better represent fractional values that might be
#      supplied.
#
#      When supplying raw matrix scores, including the SCALE option is highly
#      encouraged.  Providing it later might not do what you intend.
#
#      In the case of accessing an existing matrix by name, the scale is
#      is already associated with the matrix; supplying a conflicting value
#      will override this, which you probably do not wish to do.
#
#-------------------------------------------------------------------------------
#  Methods:
#
#  Set or retrieve the current gap score.  The value should be < 0.
#
#      $gap_score = $MatrixObj->gap_score( $new_score )
#      $gap_score = $MatrixObj->gap_score( )
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  Set or retrieve the name of the matrix.  This does not change the data.
#
#      $name = $MatrixObj->name( $name )
#      $name = $MatrixObj->name( )
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  Set or retrieve the scale used to interpret the matrix.  This does not
#  change the data matrix, but the evaluation functions need it.  Setting
#  the value of scale when the object is created is much safer.
#
#      $scale                  = $MatrixObj->scale( $scale )
#    ( $scale, $score_offset ) = $MatrixObj->scale( $scale )
#      $scale                  = $MatrixObj->scale( )
#    ( $scale, $score_offset ) = $MatrixObj->scale( )
#
#  Rescale all scores by a factor.  This changes the score matrix values, and
#  scale value, so reported bit scores should remain unchanged.
#
#      $scale                  = $MatrixObj->rescale( $factor )
#    ( $scale, $score_offset ) = $MatrixObj->rescale( $factor )
#
#  The score offset it the value added to the scores when encoding them as a
#  character string.
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  Add to or retrieve the comments associated with the matrix.  This does not
#  change the data matrix.
#
#    @comments = comments(  @comments )
#   \@comments = comments(  @comments )
#    @comments = comments( \@comments )
#   \@comments = comments( \@comments )
#    @comments = comments( )
#   \@comments = comments( )
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  These are the cool functions
#
#  There are 3 alignment evaluation methods, for each of 3 input data types.
#
#  The evaluation methods are:
#
#      1. Return alignment position scores, as
#         1a. (in scalar context) a string of unnormalized score values, or
#         1b. (in list context) a list of numerical bit scores;
#      2. Return the normalized bit score; and
#      3. Return the fraction of positions with a positive score.
#
#  The input data types are:
#
#      1. A pair of sequences,
#      2. A pair of sequence entries (i.e., [ $id, $def, $seq ] triples), and
#      3. A BlastInterface HSP.
#
#  The method invocations are:
#
#     $raw_score_string = $MatrixObj->score_seqs(  $seq1, $seq2 )
#     @bit_score_list   = $MatrixObj->score_seqs(  $seq1, $seq2 )
#     $norm_bit_score   = $MatrixObj->nbs_of_seqs( $seq1, $seq2 )
#     $fract_positives  = $MatrixObj->pos_of_seqs( $seq1, $seq2 )
#
#     $raw_score_string = $MatrixObj->score_entries(  $id_def_seq1, $id_def_seq2 )
#     @bit_score_list   = $MatrixObj->score_entries(  $id_def_seq1, $id_def_seq2 )
#     $norm_bit_score   = $MatrixObj->nbs_of_entries( $id_def_seq1, $id_def_seq2 )
#     $fract_positives  = $MatrixObj->pos_of_entries( $id_def_seq1, $id_def_seq2 )
#
#     $raw_score_string = $MatrixObj->score_hsp(  $hsp )
#     @bit_score_list   = $MatrixObj->score_hsp(  $hsp )
#     $norm_bit_score   = $MatrixObj->nbs_of_hsp( $hsp )
#     $fract_positives  = $MatrixObj->pos_of_hsp( $hsp )
#
#  Esch of these three scoring methods can be output as a code reference to a
#  subroutine, which can then be used without the object.  The subroutines
#  behave identically to the corresponding object methods described above.
#
#    \&score_func                         = $MatrixObj->seq_scr_func();
#  ( \&func, $mat_name, $scale, $offset ) = $MatrixObj->seq_scr_func();
#    \&nbs_func                           = $MatrixObj->seq_nbs_func();
#    \&pos_func                           = $MatrixObj->seq_pos_func();
#
#    \&score_func                         = $MatrixObj->entry_scr_func();
#  ( \&func, $mat_name, $scale, $offset ) = $MatrixObj->entry_scr_func();
#    \&nbs_func                           = $MatrixObj->entry_nbs_func();
#    \&pos_func                           = $MatrixObj->entry_pos_func();
#
#    \&score_func                         = $MatrixObj->hsp_scr_func();
#  ( \&func, $mat_name, $scale, $offset ) = $MatrixObj->hsp_scr_func();
#    \&nbs_func                           = $MatrixObj->hsp_nbs_func();
#    \&pos_func                           = $MatrixObj->hsp_pos_func();
#
#
#  The score offset it the value added to the scores when encoding them as a
#  character string.  The scale value is then used to convert to bit score.
#  That is:
#
#     @bit_score_list = map { $scale * ( ord($_) - $offset ) } split //, $raw_score_string;
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  Print the matrix in PAM format:
#
#        print_matrix( )            #  Print to STDOUT
#        print_matrix( \*FH )       #  Print to open file handle
#        print_matrix(  $file )     #  Print to named file
#        print_matrix( \$string )   #  Print to text string
#
#-------------------------------------------------------------------------------
#  Get a list of the available matrices.
#
#     @matrix_data = AminoAcidMatrix::matrix_list()
#    \@matrix_data = AminoAcidMatrix::matrix_list()
#
#  For each matrix one gets
#
#        [ name,
#          matrix_family,  # BLOSUM, PAM, etc.
#          matrix_number,  # meaning varies by matrix family
#          matrix_scale,   # bits per score unit
#          inv_scale       # usually units per bit
#        ]
#
#===============================================================================

use strict;
use ScoringMatrices qw( @aa_as_text );
use Data::Dumper;

#  offset in list   0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
#  index for score  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
my @aa_order  = qw( R  K  Q  E  N  D  H  G  S  T  A  U  C  V  I  J  L  M  F  Y  W  P  X  *  - );

#  The residue symbols that are supported
my %our_symb = map { $_ => 1 } @aa_order;

my @canonical_aa = qw( A C D E F G H I K L M N P Q R S T V W Y );
my %canonical_aa = map { $_ => 1, lc $_ => 1 } @canonical_aa;


sub new
{
    my $opts = ref( $_[-1] ) eq 'HASH' ? pop : {};

    my %opt_key = map { $_ => 1 } qw( ARRAY COMMENTS FILE GAP HASH NAME SCALE TEXT );

    my @source;

    while ( @_ > 1 && $opt_key{ $_[-2] } )
    {
        my ( $key, $val ) = splice @_, @_-2, 2;
        if    ( $key eq 'ARRAY' ) { push @source, [ $key, $val ] }
        elsif ( $key eq 'FILE' )  { push @source, [ $key, $val ] }
        elsif ( $key eq 'HASH' )  { push @source, [ $key, $val ] }
        elsif ( $key eq 'TEXT' )  { push @source, [ $key, $val ] }

        elsif ( $key eq 'COMMENTS' )
        {
            $opts->{ comments } = $val if $val & ref( $val ) eq 'ARRAY';
        }
        elsif ( $key eq 'GAP' )   { $opts->{ gap_score } = $val }
        elsif ( $key eq 'NAME' )  { $opts->{ name }      = $val }
        elsif ( $key eq 'SCALE' ) { $opts->{ scale }     = $val }
    }

    if ( @_ > 1 )
    {
        print STDERR "AminoAcidMatrix::new(): Bad option '$_[-2]'.";
        return undef;
    }

    push @source, [ 'MAT_NAME', shift ]  if @_;

    if ( ! @source )
    {
        if ( $opts->{ name } )
        {
            push @source, [ 'MAT_NAME', $opts->{ name } ];
            delete $opts->{ name };
        }
        else
        {
            push @source, [ 'MAT_NAME', 'BLOSUM62' ];
        }
    }

    if ( @source > 1 )
    {
        print STDERR "AminoAcidMatrix::new(): Multiple data sources supplied:\n";
        foreach ( map { $_->[0] } @source ) { print "   $_\n" }
        return undef;
    }

    my ( $type, $val ) = @{ $source[0] };
    my $self;

    if ( $type eq 'FILE' )
    {
        $self = _ingest_file( $val, $opts )
            or return undef;
        $self = _fill_missing( $self, $opts );
    }

    elsif ( $type eq 'TEXT' )
    {
        $val && ( ! ref( $val ) || ref( $val ) eq 'SCALAR' )
            or print STDERR "AminoAcidMatrix::new(): TEXT option supplied with invalid data.\n"
                and return undef;

        $self = _ingest_file( ref( $val ) ? $val : \$val, $opts )
            or return undef;
        $self = _fill_missing( $self, $opts );
    }

    elsif ( $type eq 'ARRAY' )
    {
        $self = _ingest_array( $val, $opts )
            or return undef;
        $self = _fill_missing( $self, $opts );
    }

    elsif ( $type eq 'HASH' )
    {
        $self = _ingest_hash( $val, $opts )
            or return undef;
        $self = _fill_missing( $self, $opts );
    }

    elsif ( $type eq 'MAT_NAME' )
    {
        $val && ! ref( $val )
            or print STDERR "AminoAcidMatrix::new(): MAT_NAME value supplied with invalid data.\n"
                and return undef;

        my ( $id, $scale, $text ) = matrix_by_name( $val );
        $id && $text or return undef;

        $self = _ingest_file( \$text, $opts )
            or return undef;

        $self->{ _name }  = $id;
        $self->{ _scale } = $scale;

        $self = _fill_missing( $self, $opts );
    }

    return undef if ! $self;

    $self->{ _name }     = $opts->{ name }      if $opts->{ name };
    $self->{ _scale }    = $opts->{ scale }     if $opts->{ scale };
    $self->{ _mat_name } = $opts->{ mat_name }  if $opts->{ mat_name };

    push @{ $self->{ _comments } }, @{ $opts->{comments} } if $opts->{comments};

    bless( $self, 'AminoAcidMatrix' );
}


#-------------------------------------------------------------------------------
#  Find a scoring matrix by its name.  Imperfect matching is allowed for the
#  number, and the scale factor.
#
#    ( $id, $scale, $matrix_text ) = matrix_by_name( $name, \%opts )
#
#  The name is expected to have the format:  FAMILYnumber_inverseScale 
#  For example, BLOSUM35_5 is the BLOSUM matrix at 35% clustering, with
#  a scale of 1/5 bits per score unit.  PAM250_3 is the PAM matrix at 250
#  replacements per 100 sequence positions, with 1/3 bit per score unit.
#  For most matrices, only one scale factor is available.
#
#  Available matrices can be obtained with the function:
#
#      AminoAcidMatrix::matrix_list
#
#-------------------------------------------------------------------------------
sub matrix_by_name
{
    my ( $name, $opts ) = @_;
    $opts ||= {};

    my ( $family, $number, $inv_scale ) = _parse_mat_name( $name );
    $family
        or print STDERR "AminoAcidMatrix::new(): Invalid matrix name '$name'.\n"
            and return ();

    my %default = ( BLOSUM =>  62,
                    GONNET =>   0,
                    PAM    => 120,
                    VTML   => 200
                  );
    $number ||= $default{ $family } || 0;
    my $scale = $inv_scale ? 1/$inv_scale : 1/10;  #  Go for hi res if not given

    my @mats;
    foreach ( @ScoringMatrices::aa_as_text )
    {
        my ( $id, $s, $text ) = @$_;
        my ( $f, $n, $i ) = _parse_mat_name( $id );
        next unless $f && $f eq $family;
        $n ||= 0;
        push @mats, [ abs( $number - $n ), abs( $scale - $s ), $id, $s, $text ];
    }

    @mats
        or print STDERR "AminoAcidMatrix::new(): Failed to find matrix '$name'.\n"
            and return ();

    my ( $mat ) = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @mats;

    #  ( $id, $scale, $text )
    @$mat[2,3,4];
}


sub _parse_mat_name
{
    my ( $name, $opts ) = @_;
    my ( $family, $number, $inv_scale ) = $name =~ /^([^\d]+)(\d*)(?:_(\d+))?$/;
    $family
        or print STDERR "AminoAcidMatrix::new(): Invalid matrix name '$name'.\n"
            and return ();

    ( uc $family, $number, $inv_scale );
}


#-------------------------------------------------------------------------------
#  Get a list of the available matrices.
#
#     @matrix_data = matrix_list()
#    \@matrix_data = matrix_list()
#
#  For each matrix one gets
#
#    [ name,
#      matrix_family,
#      matrix_number,
#      matrix_scale,   # bits per score unit
#      inv_scale       # usually units per bit
#    ]
#
#-------------------------------------------------------------------------------
sub matrix_list
{
    my @mats;
    foreach ( @ScoringMatrices::aa_as_text )
    {
        my ( $id, $s ) = @$_;
        my ( $f, $n, $i ) = _parse_mat_name( $id );
        push @mats, [ $id,
                      $f || '',
                      $n || '',
                      sprintf( "%0.3f", $s || 0 ),
                      $i || ''
                    ];
    }

    @mats = sort { lc $a->[1] cmp lc $b->[1]
                ||    $a->[2] <=>    $b->[2]
                ||    $b->[3] <=>    $a->[3]
                ||    $a->[4] cmp    $b->[4]
                ||    $a->[0] cmp    $b->[0]
                 }
            @mats;

    wantarray ? @mats : \@mats;
}


#-------------------------------------------------------------------------------
#  Read matrix from a text string.
#
#      \%matrix = _ingest_text(  $text, \%opts )
#
#  This just calls _ingest_file with a reference to the text string, so it
#  equivalent to:
#
#      \%matrix = _ingest_file( \$text, \%opts )
#
#-------------------------------------------------------------------------------
sub _ingest_text { $_[0] ? _ingest_file( \$_[0], $_[1] ) : undef }


#-------------------------------------------------------------------------------
#  Read matrix from a file, a file handle, a reference to a text string, or
#  STDIN (the default with an empty string or undef).
#
#      \%matrix = _ingest_file(  $file, \%opts )
#      \%matrix = _ingest_file( \*FH,   \%opts )
#      \%matrix = _ingest_file(  '',    \%opts )
#      \%matrix = _ingest_file(  undef, \%opts )
#      \%matrix = _ingest_file( \$text, \%opts )
#
#  An attempt is made to read the scale (bits per unit score) from comment
#  lines.
#-------------------------------------------------------------------------------
sub _ingest_file
{
    my $opts = ref( $_[-1] ) eq 'HASH' ? pop : {};
    my ( $fh, $close ) = input_file_handle( @_ );

    my @comments;
    my @res_list;
    my %pair_scr;
    my $scale = 0;     #  Bits per unit score
    my $n_symb;        #  The number of residue symbols in the file
    my @res_list2;     #  A temporary copy of the residues

    local $_;
    while ( <$fh> )
    {
        chomp;
        next unless /\S/;

        if ( /^\s*#/ )
        {
            push @comments, $_;
            $scale = 1 / $1 if m#1/(\d+) Bit Units#i;
            $scale = 1 / $1 if m#scale = ln\(2\)/(\d+) #i;
            $scale = 1 / 2  if m#Units = Half-Bits#i;
            $scale = 1 / 3  if m#Units = Third-Bits#i;
            $scale = 1 / 4  if m#Units = Quarter-Bits#i;
            $scale = 1 / 5  if m#Units = Fifth-Bits#i;
        }
        elsif ( ! @res_list )
        {
            s/\s*#.*$//;
            @res_list = @res_list2 = split ' ', uc $_;
            $n_symb = @res_list;
        }
        else
        {
            #  This is the symbol we expect on the line
            my $res1 = shift @res_list2
                or print STDERR "Score matrix framing error: found extra data line '$_'."
                    and return undef;

            s/\s*#.*$//;
            my @val  = split ' ', uc $_;

            #  Data line starts with the expected residue (canonical format)
            if ( ( $val[0] eq $res1 ) && ( @val == $n_symb + 1 ) )
            {
                shift @val;
            }
            #  Data line does not include the residue (we can try)
            elsif ( @val == $n_symb )
            {
            }
            else
            {
                print STDERR "Score matrix framing error: Looking for data for '$res1', but found '$_'.";
                return undef;
            }

            next unless $our_symb{ $res1 };

            for ( my $j = 0; $j < $n_symb; $j++ )
            {
                my $res2 = $res_list[$j];
                next unless $our_symb{ $res2 };

                my $res1_res2 = "$res1,$res2";
                $pair_scr{ $res1_res2 } = $val[$j];
            }
        }
    }

    close $fh if $close;

    if ( @res_list2 )
    {
        my $bad = list_as_text( @res_list2 );
        print STDERR "Score matrix framing error: matrix rows not found for $bad.";
        return undef;
    }

    if ( my $missing = missing_aa( @res_list ) )
    {
        print STDERR "Score matrix error: matrix is missing data for $missing.";
        return undef;
    }

    my $self = { _comments => \@comments,
                 _res_list => \@res_list,
                 _pair_scr => \%pair_scr,
                 _aa_order => \@aa_order,
                 _join_chr => ','        # hash keys are join(',', $res1, $res2 )
               };

    $self->{ _scale } = $scale  if $scale;

    $self;
}


#-------------------------------------------------------------------------------
#  Read scores from a list of residue identities, and a matrix of scores in
#  that order.
#
#      \%matrix = _ingest_array(   \@res_list, \@rows,   \%opts )
#      \%matrix = _ingest_array( [ \@res_list, \@rows ], \%opts )
#
#-------------------------------------------------------------------------------
sub _ingest_array
{
    my $opts = ref( $_[-1] ) eq 'HASH' ? pop : {};

    my ( $res_list, $array );
    if    ( @_ == 2 && ref( $_[0] )      eq 'ARRAY'
                    && ! ref( $_[0]->[0] )
                    && ref( $_[1] )      eq 'ARRAY'
                    && ref( $_[1]->[0] ) eq 'ARRAY'
          )
    {
        ( $res_list, $array ) = @_;
    }

    elsif ( @_ == 1 && ref( $_[0] )      eq 'ARRAY'
                    && @{$_[0]} == 2
                    && ref( $_[0]->[0] ) eq 'ARRAY'
                    && ref( $_[0]->[1] ) eq 'ARRAY'
          )
    {
        ( $res_list, $array ) = @{ $_[0] };
    }

    else
    {
        print STDERR "AminoAcidMatrix::new(): Invalid parameters passed.\n";
        return undef;
    }

    @$res_list >= 20
        or print STDERR "AminoAcidMatrix::new(): Bad residue list passed.\n"
            and return undef;

    @$array == @$res_list
        or print STDERR "AminoAcidMatrix::new(): Mismatch between residue list and array size.\n"
            and return undef;

    ref( $array->[0] ) eq 'ARRAY' && @{$array->[0]} == @$res_list
        or print STDERR "AminoAcidMatrix::new(): Array of scores is not square.\n"
            and return undef;

    if ( my $missing  = missing_aa( $$res_list ) )
    {
        print STDERR "AminoAcidMatrix::new(): Residue list is missing data for $missing.";
        return undef;
    }

    my @res_list = map { uc } @$res_list;
    my %pair_scr;
    my $n_symb = @res_list;

    foreach my $res1 ( @res_list )
    {
        my @val = @{ shift @$array };

        #  Data line does not include the residue (we can try)
        if ( @val != $n_symb )
        {
            print STDERR "Score matrix framing error: Looking for data for '$res1', but found '$_'.";
            return undef;
        }

        next unless $our_symb{ $res1 };

        for ( my $j = 0; $j < $n_symb; $j++ )
        {
            my $res2 = $res_list[$j];
            next unless $our_symb{ $res2 };

            my $res1_res2 = "$res1,$res2";
            $pair_scr{ $res1_res2 } = $val[$j];
        }
    }

    my $self = { _res_list => \$res_list,
                 _pair_scr => \%pair_scr,
                 _aa_order => \@aa_order,
                 _join_chr => ','        # hash keys are join(',', $res1, $res2 )
               };

    $self;
}


#-------------------------------------------------------------------------------
#  Read scores from a hash, or nested hashes.
#
#      \%matrix = _ingest_hash( \%hash, \%opts )
#
#  Input hashes can have either one level, with the keys being pairs of
#  residues (keys are parsed for their first and last characters):
#
#      $hash = { AA => 1, AC => 0, ...,
#                CA => 0, CC => 3, ...,
#                ...
#              }
#
#  or a hash of hashes, with the keys being the residues per se.
#
#      $hash = { A => { A => 1, C => 0, ... },
#                C => { A => 0, C => 3, ... },
#                ...
#              }
#
#-------------------------------------------------------------------------------
sub _ingest_hash
{
    my ( $data, $opts ) = @_;
    $opts ||= {};

    ref( $data ) eq 'HASH' && keys %$data >= 20
        or print STDERR "AminoAcidMatrix::new(): Bad hash data passed.\n"
            and return undef;

    my @res_list;
    my %pair_scr;
    my ( $val0 ) = values %$data;

    if ( ! ref( $val0 ) )
    {
        my @data = grep { $our_symb{ $_->[0] } && $our_symb{ $_->[1] } }
                   map  { /^(.).*(.)$/ ? [ uc $1, uc $2, $data->{ $_ } ] : () }
                   keys %$data;

        my %res_keys1 = map { $_->[0] => 1 } @data;
        my %res_keys2 = map { $_->[1] => 1 } @data;
        @res_list     = grep { $res_keys2{ $_ } } keys %res_keys1;  # in both sets

        @res_list == keys %res_keys1
            && @res_list == keys %res_keys2
            && @data == @res_list ** 2
                or print STDERR "AminoAcidMatrix::new(): Problem with hash keys.\n"
                    and return undef;

        if ( my $missing = missing_aa( @res_list ) )
        {
            print STDERR "AminoAcidMatrix::new(): Hash is missing data for $missing.\n";
            return undef;
        }

        foreach ( @data )
        {
            my ( $res1, $res2, $val ) = @$_;
            my $res1_res2 = "$res1,$res2";
            $pair_scr{ $res1_res2 } = $val;
        }
    }

    elsif ( ref( $val0 ) eq 'HASH' )
    {
        @res_list = map { $our_symb{ uc $_ } ? uc $_ : () } keys %$data;
        if ( my $missing = missing_aa( @res_list ) )
        {
            print STDERR "AminoAcidMatrix::new(): Hash is missing data for $missing.\n";
            return undef;
        }

        foreach my $key1 ( keys %$data )
        {
            my $res1 = uc $key1;
            next unless $our_symb{ $res1 };

            my $data1 = $data->{ $key1 };
            if ( ref( $data1 ) ne 'HASH' )
            {
                print STDERR "AminoAcidMatrix::new(): Datum for '$key1' is not a hash.\n";
                return undef;
            }

            my @keys2 = grep { $data->{$_} } keys %$data1;
            @keys2 == keys %$data
                or print STDERR "AminoAcidMatrix::new(): Nested hash for '$key1' has a different set of keys.\n"
                    and return undef;

            foreach my $key2 ( @keys2 )
            {
                my $res2 = uc $key2;
                next unless $our_symb{ $res2 };

                my $res1_res2 = "$res1,$res2";
                $pair_scr{ $res1_res2 } = $data1->{ $key2 };
            }
        }
    }

    else
    {
        print STDERR "AminoAcidMatrix::new(): Invalid data hash format.\n";
        return undef;
    }

    my $self = { _res_list => \@res_list,
                 _pair_scr => \%pair_scr,
                 _aa_order => \@aa_order,
                 _join_chr => ','        # hash keys are join(',', $res1, $res2 )
               };

    $self;
}


#-------------------------------------------------------------------------------
#  Fill in missing values.
#
#  In the following:
#
#     J, U, X, * and - are literal symbols in the alignment.
#     x2 is the symbol in the second sequence.
#     S(x1,x2) is the score of aligning symbol x1 with symbol x2.
#     . indicates an operation over {A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y}
#     Although shown in one direction, all scores are symmetrical.
#
#    S(X,x2) = int( average( S(.,x2) ) ), except
#    S(X,X)  = 0
#    S(X,J)  = 1 bit
#    S(X,U)  = 1 bit
#
#    S(*,x2) = min S(.,x2), except
#    S(*,*)  = 1 bit
#    S(*,J)  = 1 bit
#    S(*,U)  = 1 bit
#
#    S(-,x2) = -1 - max( max S(.,x) - 2 * min S(.,x) ), for any x, except
#    S(-,-)  = 0
#
#    S(J,x2) = int( ( S(I,x2) + S(L,x2) ) / 2 ), except
#    S(J,X)  = 1 bit
#    S(J,*)  = 1 bit
#
#  Rewrite the values of U, they are horrible in some of the matrices:
#
#    S(U,x2) = int( ( 2 * S(C,x2) + S(S,x2) ) / 3 ), except
#    S(U,X)  = 1 bit
#    S(U,*)  = 1 bit
#
#  We might consider increasing the scores for U and J with X and * since they
#  are often shown as the X a sequence, and appear as * in a naive translation.
#-------------------------------------------------------------------------------
sub _fill_missing
{
    my ( $self, $opts ) = @_; 
    $opts ||= {};

    my $res_list = $self->{ _res_list };
    my $pair_scr = $self->{ _pair_scr };
    my $scale    = $self->{ _scale };
    my $one_bit  = $scale ? round( 1 / $scale ) : 1;
    my $join_chr = $self->{ _join_chr };

    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #  Some values that will be useful below
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    my $have_res = { map { $_ => 1 } @$res_list };

    my $max_gain   = 0;  #  Biggest gain converting 2 bad scores to 1 good score
    my $max_scr    = 0;
    my $min_scr    = 0;
    foreach my $res1 ( @canonical_aa )
    {
        my $row_max = 0;
        my $row_min = 0;
        foreach my $res2 ( @canonical_aa )
        {
            my $res1_res2 = "$res1,$res2";
            my $val = $pair_scr->{ $res1_res2 };
            $row_max = $val if $val > $row_max;
            $row_min = $val if $val < $row_min;
        }

        $max_scr = $row_max if $row_max > $max_scr;
        $min_scr = $row_min if $row_min < $min_scr;

        #  Adding 2 gaps can convert two bad scores to one good score, so
        my $gap_gain = $row_max - 2 * $row_min;
        $max_gain = $gap_gain if $gap_gain > $max_gain;
    }

    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    #  S(X,.)
    #
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if ( ! $have_res->{ X } )
    {
        foreach my $res1 ( @canonical_aa )
        {
            my $scr = 0;
            foreach ( @canonical_aa )
            {
                $scr += $pair_scr->{ "$res1,$_" }
            }
            $pair_scr->{ "$res1,X" } = $pair_scr->{ "X,$res1" } = int( $scr/20 );
        }
        $pair_scr->{ "X,X" } = 0;
        $pair_scr->{ "*,X" } = $pair_scr->{ "X,*" } = $min_scr;
        #  See note at top of file.
        $pair_scr->{ "U,X" } = $pair_scr->{ "X,U" } = $one_bit;
        $pair_scr->{ "J,X" } = $pair_scr->{ "X,J" } = $one_bit;

        push @$res_list, 'X';
        $have_res->{ X } = 1;
    }

    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    #  S(U,.)
    #
    #  For now, we are going to unconditionally rewrite the values for U, they
    #  are horrible in some of the matrices.
    #
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if ( 1 || ! $have_res->{ U } )
    {
        my $scr;
        foreach my $res1 ( @canonical_aa )
        {
            $scr = int( ( 2 * $pair_scr->{ "$res1,C" } + $pair_scr->{ "$res1,S" } ) / 3 );
            $pair_scr->{ "$res1,U" } = $pair_scr->{ "U,$res1" } = $scr;
        }
        $pair_scr->{ "U,U" } = $pair_scr->{ "C,C" };
        #  See note at top of file.
        $pair_scr->{ "X,U" } = $pair_scr->{ "U,X" } = $one_bit;
        $pair_scr->{ "*,U" } = $pair_scr->{ "U,*" } = $one_bit;

        push @$res_list, 'U';
        $have_res->{ U } = 1;
    }

    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    #  S(J,.)
    #
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if ( ! $have_res->{ J } )
    {
        my $scr;
        foreach my $res1 ( @canonical_aa, qw( U ) )
        {
            $scr = int( ( $pair_scr->{ "$res1,I" } + $pair_scr->{ "$res1,L" } ) / 2 );
            $pair_scr->{ "$res1,J" } = $pair_scr->{ "J,$res1" } = $scr;
        }
        $pair_scr->{ "J,J" } = int( ( $pair_scr->{ "I,I" } + $pair_scr->{ "L,L" } ) / 2 );
        #  See note at top of file.
        $pair_scr->{ "X,J" } = $pair_scr->{ "J,X" } = $one_bit;
        $pair_scr->{ "*,J" } = $pair_scr->{ "J,*" } = $one_bit;

        push @$res_list, 'J';
        $have_res->{ J } = 1;
    }
    elsif ( ! exists $pair_scr->{ 'J,U' } )
    {
        my $scr = int( ( 2 * $pair_scr->{ "J,C" } + $pair_scr->{ "J,S" } ) / 3 );
        $pair_scr->{ "J,U" } = $pair_scr->{ "U,J" } = $scr;
    }

    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    #  S(*,.)
    #
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if ( ! $have_res->{ '*' } )
    {
        foreach my $res1 ( @canonical_aa, qw( X ) )
        {
            $pair_scr->{ "$res1,*" } = $pair_scr->{ "*,$res1" } = $min_scr;
        }
        $pair_scr->{ "*,*" } = $one_bit;
        #  See note at top of file.
        $pair_scr->{ "*,U" } = $pair_scr->{ "U,*" } = $one_bit;
        $pair_scr->{ "*,J" } = $pair_scr->{ "J,*" } = $one_bit;

        push @$res_list, '*';
        $have_res->{ '*' } = 1;
    }

    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    #  S(-,.)
    #
    #  It should always be the case that eliminating 2 bed scores to get
    #  one good score should not justify 2 gaps.  For example,
    #
    #    ...CY...  should always be better than  ...-CY...
    #    ...YC...                                ...YC-...
    #
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if ( ! $have_res->{ '-' } )
    {
        my $gap_score = $self->{ _gap_scr }
                    ||= $opts->{  gap_score }
                    ||= -1 - int( $max_gain / 2 );

        gap_score( $self, $gap_score );
        push @$res_list, '-';
        $have_res->{ '-' } = 1;
    }


    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #  Do we want to rescale the scores for range, or to better fit non-integer
    #  values?  Regardless, we will round everthing to the nearest int.
    #
    #  If the range is too large for our scores, we scale down.  This range
    #  is 24 bits, and is mine, not due to perl.
    #
    #    $max_scr >  8388607
    #    $min_scr < -8388607
    #
    #  Are there non-integer scores that would round better if we expanded
    #  the scale?  We will only rescale by an integer multiple, so there
    #  is no relative shift of existing integers.  We could use 32,000 as
    #  the default max magnitude, but it is not clear that we need more
    #  resolution than a 2000 value scale.
    #
    #    $max_nonint * int( 2047 / max( $max_scr, -$min_scr ) ) > 0.5
    #
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    my $max_mag    = 0;
    my $max_nonint = 0;
    foreach ( map { abs $_ } values %$pair_scr )
    {
        $max_mag = $_ if $_ > $max_mag;
        my $nonint = nonint( $_ );
        $max_nonint = $nonint if $nonint > $max_nonint;
    }

    my $rescale = ( $max_mag > 8388607 ) ? 8388607 / $max_mag
                : ( $max_nonint * int( 2047 / $max_mag ) > 0.5 ) ? int( 2047 / $max_mag )
                : 1;

    rescale( $self, $rescale )
        or return undef;

    $self;
}


#-------------------------------------------------------------------------------
#  Rescale all scores by a factor
#
#      $scale                  = $MatrixObj->rescale( $factor )
#    ( $scale, $score_offset ) = $MatrixObj->rescale( $factor )
#
#  The score offset it the value added to the scores when encoding them as a
#  character string.
#-------------------------------------------------------------------------------
sub rescale
{
    my ( $self, $rescale ) = @_;
    $self && $rescale && $rescale > 0
        or print STDERR "AminoAcidMatrix::rescale: Called with bad parameters.\n"
            and return undef;

    my $pair_scr = $self->{ _pair_scr }
        or print STDERR "AminoAcidMatrix::rescale: Improperly initialized object.\n"
            and return undef;

    my $scr;
    my $max_mag = 0;
    foreach ( keys %$pair_scr )
    {
        $pair_scr->{ $_ } = $scr = round( $rescale * $pair_scr->{ $_ } );
        $max_mag = abs( $scr ) if abs( $scr ) > $max_mag;
    }

    $self->{ _max_mag } = $max_mag;

    #  Check that we do not have an out of bounds value (the limit is mine,
    #  not perl's):

    if ( $max_mag > 2147483647 )
    {
        print STDERR "AminoAcidMatrix::rescale: Scores rescaled out of range.\n";
        return undef;
    }

    my $gap_score = $self->{ _gap_scr };
    $self->{ _gap_scr } = round( $rescale * $gap_score )  if $gap_score;

    #  Mark derived data as invalid.

    _invalidate_scores( $self );

    my $scale = ( $self->{ _scale } || 1 ) / $rescale;

    $self->{ _scale } = $scale  if $self->{ _scale };

    wantarray ? ( $scale, scalar _score_offset( $self ) ) : $scale;
}


sub nonint { local $_ = abs( $_[0] || 0 ); abs( $_ - int( $_ + 0.5 ) ) }


sub round { $_[0] ? int( $_[0] + 0.5 * ( $_[0] / abs( $_[0] ) ) ) : 0 }


#-------------------------------------------------------------------------------
#  Set the gap score, or return current gap score.  The value should be < 0.
#
#      $gap_score = $MatrixObj->gap_score( $new_score )
#      $gap_score = $MatrixObj->gap_score( )
#
#-------------------------------------------------------------------------------
sub gap_score
{
    my ( $self, $gap_score ) = @_;
    return undef               unless $self;
    return $self->{ _gap_scr } unless defined $gap_score;

    #  Gap scores must be negative
    if ( $gap_score > 0 )
    {
        print STDERR "AminoAcidMatrix::gap_score: Postive gap score converted to negative.\n";
        $gap_score *= -1;
    }

    my $pair_scr =    $self->{ _pair_scr };
    my @aa_order = @{ $self->{ _aa_order } || [] };
    $pair_scr && @aa_order
        or print STDERR "AminoAcidMatrix::gap_score: Improperly initialized object.\n"
            and return undef;

    foreach my $res1 ( @aa_order )
    {
        next if $res1 eq '-';
        $pair_scr->{ "$res1,-" } = $pair_scr->{ "-,$res1" } = $gap_score;
    }
    $pair_scr->{ "-,-" } = 0;

    my $max_mag = $self->{ _max_mag };
    $self->{ _max_mag } = abs( $gap_score )  if $max_mag && abs( $gap_score ) > $max_mag;

    #  Mark derived data as invalid.

    _invalidate_scores( $self );

    $self->{ _gap_scr } = $gap_score;
}


#-------------------------------------------------------------------------------
#  Set or retrieve the name of the matrix.  This does not change the data.
#
#      $name = $MatrixObj->name( $name )
#      $name = $MatrixObj->name( )
#
#-------------------------------------------------------------------------------
sub name
{
    my ( $self, $name ) = @_;
    return undef unless $self;

    defined $name ? ( $self->{ _name } = $name ) : $self->{ _name };
}


#-------------------------------------------------------------------------------
#  Set or retrieve the scale used to interpret the matrix.  This does not
#  change the data matrix, but the evaluation functions include it.
#
#      $scale                  = $MatrixObj->scale( $scale )
#    ( $scale, $score_offset ) = $MatrixObj->scale( $scale )
#      $scale                  = $MatrixObj->scale( )
#    ( $scale, $score_offset ) = $MatrixObj->scale( )
#
#  The score offset it the value added to the scores when encoding them as a
#  character string.
#-------------------------------------------------------------------------------
sub scale
{
    my ( $self, $scale ) = @_;
    return wantarray ? () : undef  unless $self;

    if ( $scale )
    {
        delete $self->{ _seq_scr_func };
        delete $self->{ _entry_scr_func };
        delete $self->{ _hsp_scr_func };

        $self->{ _scale } = $scale;
    }
    else
    {
        $scale = $self->{ _scale };
    }

    wantarray ? ( $scale, _score_offset( $self ) ) : $scale;
}


#-------------------------------------------------------------------------------
#  Add to or retrieve the comments associated with the matrix.  This does not
#  change the data matrix.
#
#      @comments = $MatrixObj->comments(  @comments )
#     \@comments = $MatrixObj->comments(  @comments )
#      @comments = $MatrixObj->comments( \@comments )
#     \@comments = $MatrixObj->comments( \@comments )
#      @comments = $MatrixObj->comments( )
#     \@comments = $MatrixObj->comments( )
#
#-------------------------------------------------------------------------------
sub comments
{
    my $self = shift;
    return wantarray ? () : undef unless $self;

    my $com;
    if ( @_ )
    {
        $com = $self->{ _comments } ||= [];
        push @$com, map { s/^\s*$/#/;
                          s/^(\s*[^#])/#  $1/;
                          $_;
                        }
                    ref( $_[0] ) eq "ARRAY" ? @{ $_[0] } : @_;
    }
    else
    {
        $com = $self->{ _comments } || [];
    }        

    wantarray ? @$com : $com;
}


#-------------------------------------------------------------------------------
#  Character transliteration from amino acid sequence to offsets used in the
#  scoring matrices.
#
#  If we map a first sequence into the range 1-31 and a second sequence into
#  the range 32-992 (step 32), then the logical or of their bits will uniquely
#  identify the combination of residues, allowing a direct score lookup.
#
#      $s1_index    = seq_to_index( $seq1 );
#      $s2_index    = seq2_to_index( $seq2 );
#      $s1_s2_index = $s1_index | $s2_index;
#  or
#      $s1_s2_index = seq_to_index( $seq1 ) | seq2_to_index( $seq2 );
#
#  One liners to make the translation strings:
#
#   perl -e 'foreach ( 23, 25 ) { for ($i=1;$i<=$_;$i++) {printf "\\x%02X",      $i} } print "\n"'
#   perl -e 'foreach ( 23, 25 ) { for ($i=1;$i<=$_;$i++) {printf "\\x{%03X}", 32*$i} } print "\n"'
#
#-------------------------------------------------------------------------------
sub seq_to_index
{
    local $_ = shift;
    tr/RKQENDHGSTAUCVIJLMFYWPXrkqendhgstaucvijlmfywpx*-/X/c;
    tr/RKQENDHGSTAUCVIJLMFYWPXrkqendhgstaucvijlmfywpx*-/\x01\x02\x03\x04\x05\x06\x07\x08\x09\x0A\x0B\x0C\x0D\x0E\x0F\x10\x11\x12\x13\x14\x15\x16\x17\x01\x02\x03\x04\x05\x06\x07\x08\x09\x0A\x0B\x0C\x0D\x0E\x0F\x10\x11\x12\x13\x14\x15\x16\x17\x18\x19/;
    $_;
}


sub seq2_to_index
{
    local $_ = shift;
    tr/RKQENDHGSTAUCVIJLMFYWPXrkqendhgstaucvijlmfywpx*-/X/c;
    tr/RKQENDHGSTAUCVIJLMFYWPXrkqendhgstaucvijlmfywpx*-/\x{020}\x{040}\x{060}\x{080}\x{0A0}\x{0C0}\x{0E0}\x{100}\x{120}\x{140}\x{160}\x{180}\x{1A0}\x{1C0}\x{1E0}\x{200}\x{220}\x{240}\x{260}\x{280}\x{2A0}\x{2C0}\x{2E0}\x{020}\x{040}\x{060}\x{080}\x{0A0}\x{0C0}\x{0E0}\x{100}\x{120}\x{140}\x{160}\x{180}\x{1A0}\x{1C0}\x{1E0}\x{200}\x{220}\x{240}\x{260}\x{280}\x{2A0}\x{2C0}\x{2E0}\x{300}\x{320}/;
    $_;
}


sub seq_align_index { seq_to_index( $_[0] ) | seq2_to_index( $_[1] ) }


#-------------------------------------------------------------------------------
#  Retrieve (building if necessary) a 26 x 26 scoring matrix from residue index
#  values to corresponding scores.  The index 0 is not used.
#
#     \@matrix                = $MatrixObj->scoring_matrix( \%opts )
#   ( \@matrix, $id, $scale ) = $MatrixObj->scoring_matrix( \%opts )
#
#  @aa_order  = qw( R  K  Q  E  N  D  H  G  S  T  A  U  C  V  I  J  L  M  F  Y  W  P  X  *  - );
#  offset in list   0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
#  index for score  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
#
#-------------------------------------------------------------------------------
sub scoring_matrix
{
    my ( $self, $opts ) = @_;
    $opts ||= {};

    my $mat = $self->{ _matrix }
           || _build_scoring_matrix( $self, $opts )
        or return undef;

    my $id    = $self->{ _name } || 'Unnamed';
    my $scale = $self->{ _scale };

    wantarray ? ( $mat, $id, $scale ) : $mat;
}


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  Convert the hash of scores into a matrix.
#
#    \Smatrix_rows = _build_scoring_matrix( $self, \%opts )
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub _build_scoring_matrix
{
    my ( $self, $opts ) = @_;
    $opts ||= {};

    my $pair_scr = $self->{ _pair_scr };
    my @aa_order = @{ $self->{ _aa_order } || [] };
    my $n_symb   = @aa_order;

    $pair_scr && $n_symb
        or print STDERR "AminoAcidMatrix::_build_scoring_matrix: Improperly initialized object.\n"
            and return undef;

    my @mat = ( [ (0) x ($n_symb+1) ] );  #  First row is all zeros (it is not used)
    for ( my $i = 1; $i <= $n_symb; $i++ )
    {
        my $res1 = $aa_order[$i-1];

        my @row = ( 0 );         #  First column is all zeros (it is not used)
        for ( my $j = 1; $j <= $n_symb; $j++ )
        {
            push @row, $pair_scr->{ "$res1,$aa_order[$j-1]" } || 0;
        }
        push @mat, \@row;
    }

    # print_matrix( \@mat ); exit;

    $self->{ _matrix } = \@mat;
}



#-------------------------------------------------------------------------------
#  It is possible to change one or more scores after object creation, so
#  this is a way to delete derived data.
#
#      _invalidate_scores( $self )
#
#-------------------------------------------------------------------------------
sub _invalidate_scores
{
    my $self = shift
        or print STDERR "AminoAcidMatrix::_invalidate_scores: Improperly initialized object.\n"
            and return undef;

    delete $self->{ _matrix };
    delete $self->{ _scr_offset };
    delete $self->{ _score_str };
    delete $self->{ _seq_scr_func };
    delete $self->{ _entry_scr_func };
    delete $self->{ _hsp_scr_func };

    1;
}


#-------------------------------------------------------------------------------
#  Make a subroutine to do alignment scoring, based on the template:
#
#  The score string can incorporate scores that use characters larger than
#  one byte.  They dynamically adjust based on the _max_mag value in the object,
#  changing the score offset and the number of hex characters used to represent
#  the character values.
#-------------------------------------------------------------------------------
sub _index_and_score_strings
{
    my ( $self, $opts ) = @_;
    $opts ||= {};

    my $indices = $self->{ _index_str };
    my $scores  = $self->{ _score_str };

    if ( ! ( $indices && $scores ) )
    {
        my $mat = $self->{ _matrix }
               || _build_scoring_matrix( $self, $opts )
            or return undef;

        my $n_symb = @$mat - 1;

        #
        #  Handle larger score ranges, if we ever want to use them.
        #  $self->{ _max_mag } should be set when pair_scr hash is done.
        #
        my ( $scr_offset, $scr_hex_digits ) = _score_offset( $self );
        my $scr_format = "\\x{%0${scr_hex_digits}X}";

        my @index_pair = ();
        my @score      = ();
        for ( my $i = 1; $i <= $n_symb; $i++ )
        {
            my $mat_row = $mat->[$i];
            for ( my $j = 1; $j <= $n_symb; $j++ )
            {
                push @index_pair, 32 * $i + $j;
                push @score, $mat_row->[$j] + $scr_offset;
            }
        }

        $self->{ _index_str } = $indices = join('', map { sprintf "\\x{%03X}", $_ } @index_pair );
        $self->{ _score_str } = $scores  = join('', map { sprintf $scr_format, $_ } @score );
    }

    ( $indices, $scores );
}


#-------------------------------------------------------------------------------
#  Find the number of hex digits required to represent the scores, and the
#  score offset to make all values positive.
#
#      $score_offset                = _score_offset( $self )
#    ( $score_offest, $hex_digits ) = _score_offset( $self )
#
#-------------------------------------------------------------------------------
sub _score_offset
{
    my $self = shift
        or print STDERR "AminoAcidMatrix::_score_offset: Called without object.\n"
            and return undef;

    if ( $self->{ _scr_offset } )
    {
        return wantarray ? ( $self->{ _scr_offset }, $self->{ _scr_hex_digits } )
                         :   $self->{ _scr_offset };
    }

    my $max_mag = $self->{ _max_mag } || 127;  #  set when pair hash is done
    my $scr_hex_digits = 2;
    while ( ( 16 ** $scr_hex_digits ) / 2 <= $max_mag ) { $scr_hex_digits++ }
    my $scr_offset = ( 16 ** $scr_hex_digits ) / 2;
    $self->{ _scr_offset }     = $scr_offset;
    $self->{ _scr_hex_digits } = $scr_hex_digits;

    wantarray ? ( $scr_offset, $scr_hex_digits ) : $scr_offset;
}


#-------------------------------------------------------------------------------
#  These are the cool functions
#
#  There are 3 alignment evaluation methods, for each of 3 input data types.
#
#  The evaluation methods are:
#
#      1. Return alignment position scores, as
#         1a. (in scalar context) a string of unnormalized score values, or
#         1b. (in list context) a list of numerical bit scores;
#      2. Return the normalized bit score; and
#      3. Return the fraction of positions with a positive score.
#
#  The input data types are:
#
#      1. A pair of sequences,
#      2. A pair of sequence entries (i.e., [ $id, $def, $seq ] triples), and
#      3. A BlastInterface HSP.
#
#  The method invocations are:
#
#     $raw_score_string = $MatrixObj->score_seqs(  $seq1, $seq2 )
#     @bit_score_list   = $MatrixObj->score_seqs(  $seq1, $seq2 )
#     $norm_bit_score   = $MatrixObj->nbs_of_seqs( $seq1, $seq2 )
#     $fract_positives  = $MatrixObj->pos_of_seqs( $seq1, $seq2 )
#
#     $raw_score_string = $MatrixObj->score_entries(  $id_def_seq1, $id_def_seq2 )
#     @bit_score_list   = $MatrixObj->score_entries(  $id_def_seq1, $id_def_seq2 )
#     $norm_bit_score   = $MatrixObj->nbs_of_entries( $id_def_seq1, $id_def_seq2 )
#     $fract_positives  = $MatrixObj->pos_of_entries( $id_def_seq1, $id_def_seq2 )
#
#     $raw_score_string = $MatrixObj->score_hsp(  $hsp )
#     @bit_score_list   = $MatrixObj->score_hsp(  $hsp )
#     $norm_bit_score   = $MatrixObj->nbs_of_hsp( $hsp )
#     $fract_positives  = $MatrixObj->pos_of_hsp( $hsp )
#
#  Esch of these three scoring methods can be output as a code reference to a
#  subroutine, which can then be used without the object.  The subroutines
#  behave identically to the corresponding object methods described above.
#
#    \&score_func                         = $MatrixObj->seq_scr_func();
#  ( \&func, $mat_name, $scale, $offset ) = $MatrixObj->seq_scr_func();
#    \&nbs_func                           = $MatrixObj->seq_nbs_func();
#    \&pos_func                           = $MatrixObj->seq_pos_func();
#
#    \&score_func                         = $MatrixObj->entry_scr_func();
#  ( \&func, $mat_name, $scale, $offset ) = $MatrixObj->entry_scr_func();
#    \&nbs_func                           = $MatrixObj->entry_nbs_func();
#    \&pos_func                           = $MatrixObj->entry_pos_func();
#
#    \&score_func                         = $MatrixObj->hsp_scr_func();
#  ( \&func, $mat_name, $scale, $offset ) = $MatrixObj->hsp_scr_func();
#    \&nbs_func                           = $MatrixObj->hsp_nbs_func();
#    \&pos_func                           = $MatrixObj->hsp_pos_func();
#
#
#  The score offset it the value added to the scores when encoding them as a
#  character string.  The scale value is then used to convert to bit score.
#  That is:
#
#     @bit_score_list = map { $scale * ( ord($_) - $offset ) } split //, $raw_score_string;
#
#-------------------------------------------------------------------------------
#  Usages:
#
#      print Dumper( [ AminoAcidMatrix::new()->score_seqs("ACGTVW","CGGSLY") ] );
#      $VAR1 = [
#               '-0.25',
#               '-1.25',
#                '2.75',
#                '0.75',
#                '0.5',
#                '1'
#              ];
#
#      print Dumper( scalar AminoAcidMatrix::new()->score_seqs("ACGTVW","CGGSLY") );
#      $VAR1 = "{\x{8b}\x{83}\x{82}\x{84}";
#
#  or creating a function and using it:
#
#      my   $sub                = $AminoAcidMatrix::new->seq_scr_func();
#      my ( $sub, $id, $scale ) = $AminoAcidMatrix::new->seq_scr_func();
#
#  Then, get a list of bit scores:
#
#      my @bit_scores = &$sub( $seq1, $seq2 );
#
#  or get a string with ordinal values of characters equal to raw score + 128:
#
#      my $raw_scores_plus_128 = &$sub( $seq1, $seq2 );
#
if ( 0 )
{
my $x = <<'End_of_Code'
set seq1=MRYISTRGQAPALNFEDVLLAGLASDGGLYVPENLPRFTLEEIASWVGLPYHELAFRVMR
set seq2=MKLYNLKDHNEQVSFAQAVTQGLGKNQGLFFPHDLPEFSLTEIDEMLKLDFVTRSAKILS

#  Run function in array context
perl -e 'use strict; \
use AminoAcidMatrix; \
my ($func,$mat,$scale)=AminoAcidMatrix::new("BLOSUM62")->seq_scr_func; \
print "Matrix = $mat; \
Scale = $scale bits per unit score\n"; \
print "List of per position bit scores:\n", join(", ", &$func(@ARGV)), "\n"'  $seq1 $seq2 

#  Run function in scaler context
perl -e 'use strict; \
use AminoAcidMatrix; \
my $func=AminoAcidMatrix::new("BLOSUM62")->seq_scr_func; \
use Data::Dumper; \
print Dumper( "Character string with ordinal values that are raw score + 128", scalar &$func(@ARGV) )'  $seq1 $seq2 

#  Speed test
time perl -e 'use strict; \
use AminoAcidMatrix; \
my $func=AminoAcidMatrix::new("BLOSUM62")->seq_scr_func; \
for (my $i=0; \
$i<1000000; \
$i++) { my $a = &$func(@ARGV) }'  $seq1 $seq2 
7.667u 0.004s 0:07.67 99.8%	0+0k 0+0io 0pf+0w

perl -e 'use strict; \
use AminoAcidMatrix; \
my $mat = AminoAcidMatrix::new(); \
$mat->rescale(100); \
my $func = $mat->seq_scr_func; \
print $mat->scale(), "\n", join(", ", &{$func}(@ARGV)), "\n";'  $seq1 $seq2 

perl -e 'use strict; \
use AminoAcidMatrix; \
use Data::Dumper; \
my $mat = AminoAcidMatrix::new(); \
$mat->rescale(2); \
my ($func,$name,$scale,$offset) = $mat->seq_scr_func; \
print Dumper( $name, $scale, $offset, scalar( &{$func}(@ARGV)) );'  $seq1 $seq2 

set seq1=GKGLGTKLVRALVELLFNDPEVTKIQTDPSPSNLRAIRCYEKAGFERQGT
set seq2=NRGVASALMRTMIDMCDNWLRVERIELTVFADNAPAIAVYKKYGFEIEGT

perl -e 'use strict; \
use AminoAcidMatrix; \
my $mat = AminoAcidMatrix::new("BLOSUM62_5"); \
my $func1 = $mat->seq_nbs_func; \
my $func2 = $mat->seq_pos_func; \
printf "%.3f, %.3f\n", &$func1(@ARGV), &$func2(@ARGV)'  $seq1 $seq2

End_of_Code
}
#
#-------------------------------------------------------------------------------
#
#  Evaluation methods for a pair of sequences
#
sub score_seqs
{
    my $self = shift
        or return wantarray ? () : undef;

    my $func = $self->{ _seq_scr_func }
            || ( _build_score_funcs( $self ) && $self->{ _seq_scr_func } )
        or return wantarray ? () : undef;

    &$func( @_ );
}

sub nbs_of_seqs
{
    my $self = shift
        or return wantarray ? () : undef;

    my $func = $self->{ _seq_nbs_func }
            || ( _build_score_funcs( $self ) && $self->{ _seq_nbs_func } )
        or return wantarray ? () : undef;

    &$func( @_ );
}

sub pos_of_seqs
{
    my $self = shift
        or return wantarray ? () : undef;

    my $func = $self->{ _seq_pos_func }
            || ( _build_score_funcs( $self ) && $self->{ _seq_pos_func } )
        or return wantarray ? () : undef;

    &$func( @_ );
}

#
#  Evaluation methods for a pair of sequence entries (triples)
#
sub scr_entries
{
    my $self = shift
        or return wantarray ? () : undef;

    my $func = $self->{ _entry_scr_func }
            || ( _build_score_funcs( $self ) && $self->{ _entry_scr_func } )
        or return wantarray ? () : undef;

    &$func( @_ );
}

sub nbs_of_entries
{
    my $self = shift
        or return wantarray ? () : undef;

    my $func = $self->{ _entry_nbs_func }
            || ( _build_score_funcs( $self ) && $self->{ _entry_nbs_func } )
        or return wantarray ? () : undef;

    &$func( @_ );
}

sub pos_of_entries
{
    my $self = shift
        or return wantarray ? () : undef;

    my $func = $self->{ _entry_pos_func }
            || ( _build_score_funcs( $self ) && $self->{ _entry_pos_func } )
        or return wantarray ? () : undef;

    &$func( @_ );
}

#
#  Evaluation methods for a BlastInterface HSP
#
sub scr_hsp
{
    my $self = shift
        or return wantarray ? () : undef;

    my $func = $self->{ _hsp_scr_func }
            || ( _build_score_funcs( $self ) && $self->{ _hsp_scr_func } )
        or return wantarray ? () : undef;

    &$func( @_ );
}

sub nbs_of_hsp
{
    my $self = shift
        or return wantarray ? () : undef;

    my $func = $self->{ _nbs_of_hsp }
            || ( _build_score_funcs( $self ) && $self->{ _nbs_of_hsp } )
        or return wantarray ? () : undef;

    &$func( @_ );
}

sub pos_of_hsp
{
    my $self = shift
        or return wantarray ? () : undef;

    my $func = $self->{ _pos_of_hsp }
            || ( _build_score_funcs( $self ) && $self->{ _pos_of_hsp } )
        or return wantarray ? () : undef;

    &$func( @_ );
}

#-------------------------------------------------------------------------------
#  Methods that return an evaluation function for each of the above
#-------------------------------------------------------------------------------
#
#  Functions for a pair of sequences
#
sub seq_scr_func
{
    my ( $self, $opts ) = @_;
    $self or return wantarray ? () : undef;
    $opts ||= {};

    my $func = $self->{ _seq_scr_func }
            || ( $self->_build_score_funcs( $opts ) && $self->{ _seq_scr_func } )
        or return wantarray ? () : undef;

    wantarray ? ( $func, $self->{_name}, $self->{_scale}, $self->{_scr_offset} )
              : $func;
}

sub seq_nbs_func
{
    my ( $self, $opts ) = @_;
    $self or return undef;
    $opts ||= {};

    $self->{ _seq_nbs_func }
        || ( $self->_build_score_funcs( $opts ) && $self->{ _seq_nbs_func } );
}

sub seq_pos_func
{
    my ( $self, $opts ) = @_;
    $self or return undef;
    $opts ||= {};

    $self->{ _seq_pos_func }
        || ( $self->_build_score_funcs( $opts ) && $self->{ _seq_pos_func } );
}

#
#  Functions for a pair of sequence entries (triples)
#
sub entry_scr_func
{
    my ( $self, $opts ) = @_;
    $self or return wantarray ? () : undef;
    $opts ||= {};

    my $func = $self->{ _entry_scr_func }
            || ( $self->_build_score_funcs( $opts ) && $self->{ _entry_scr_func } )
        or return wantarray ? () : undef;

    wantarray ? ( $func, $self->{_name}, $self->{_scale}, $self->{_scr_offset} )
              : $func;
}

sub seq_nbs_func
{
    my ( $self, $opts ) = @_;
    $self or return undef;
    $opts ||= {};

    $self->{ _seq_nbs_func }
        || ( $self->_build_score_funcs( $opts ) && $self->{ _seq_nbs_func } );
}

sub seq_pos_func
{
    my ( $self, $opts ) = @_;
    $self or return undef;
    $opts ||= {};

    $self->{ _seq_pos_func }
        || ( $self->_build_score_funcs( $opts ) && $self->{ _seq_pos_func } );
}

#
#  Functions for a BlastInterface HSP
#
sub hsp_scr_func
{
    my ( $self, $opts ) = @_;
    $self or return wantarray ? () : undef;
    $opts ||= {};

    my $func = $self->{ _hsp_scr_func }
            || ( $self->_build_score_funcs( $opts ) && $self->{ _hsp_scr_func } )
        or return wantarray ? () : undef;

    wantarray ? ( $func, $self->{_name}, $self->{_scale}, $self->{_scr_offset} )
              : $func;
}

sub hsp_nbs_func
{
    my ( $self, $opts ) = @_;
    $self or return undef;
    $opts ||= {};

    $self->{ _hsp_nbs_func }
        || ( $self->_build_score_funcs( $opts ) && $self->{ _hsp_nbs_func } );
}

sub hsp_pos_func
{
    my ( $self, $opts ) = @_;
    $self or return undef;
    $opts ||= {};

    $self->{ _hsp_pos_func }
        || ( $self->_build_score_funcs( $opts ) && $self->{ _hsp_pos_func } );
}


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  The scoring functions are very similar, they just extract the sequences from
#  the parameters in different ways.  It makes the most sense to create them
#  all at once.
#
#       $okay = $MatrixObj->_build_score_funcs( )
#
#  The functions themselves are placed in $self.
#
#  The fraction positives and nbs scoring functions remove columns with shared
#  gap characters.
#
#  These functions can incorporate scores that use characters larger than
#  one byte.  Only the score offset changes.
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  Prototype for creating function references that do transliteration with
#  runtime data in them.
#
#  use strict;
#  my $a = "1234";
#  my $sub;
#  eval "\$sub = sub { local \$_ = shift; tr/abcd/$a/; wantarray ? map { ord(\$_)-128 } split // : \$_ }";
#  die $@ if $@;
#
#  See the string:
#  print scalar &$sub("abcdABCD"), "\n";
#
#  See the character value list:
#  print join( ", ", &$sub("abcdABCD")), "\n";
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub _build_score_funcs
{
    my ( $self, $opts ) = @_;
    $opts ||= {};

    my ( $indices, $scores ) = _index_and_score_strings( $self, $opts );
    $indices && $scores
        or return undef;

    my $scale      = $self->{ _scale }      ||   1;
    my $scr_offset = $self->{ _scr_offset } || 128;

    #  Scores from two input sequences:

    my $seq_scr_func;
    eval "\$seq_scr_func = sub { local \$_ = seq_align_index( \$_[0], \$_[1] );
                                 tr/$indices/$scores/;
                                 wantarray ? map { $scale * ( ord(\$_)-$scr_offset ) } split // : \$_;
                               }";
    print STDERR $@ and return undef if $@;
    $self->{ _seq_scr_func } = $seq_scr_func;

    my $seq_nbs_func;
    eval "\$seq_nbs_func = sub { local \$_ = seq_align_index( \$_[0], \$_[1] );
                                 tr/\x{339}//d;
                                 return undef if ! length( \$_ );
                                 tr/$indices/$scores/;
                                 my \$ttl = 0;
                                 foreach ( split // ) { \$ttl += ord( \$_ ) }
                                 $scale * ( \$ttl / length( \$_ ) - $scr_offset );
                               }";
    print STDERR $@ and return undef if $@;
    $self->{ _seq_nbs_func } = $seq_nbs_func;

    my $seq_pos_func;
    eval "\$seq_pos_func = sub { local \$_ = seq_align_index( \$_[0], \$_[1] );
                                 tr/\x{339}//d;
                                 return undef if ! length( \$_ );
                                 tr/$indices/$scores/;
                                 my \$ttl = 0;
                                 foreach ( split // ) { \$ttl++  if ord( \$_ ) > $scr_offset }
                                 \$ttl / length( \$_ );
                               }";
    print STDERR $@ and return undef if $@;
    $self->{ _seq_pos_func } = $seq_pos_func;

    #  Scores from two input sequence entries (triples):

    my $entry_scr_func;
    eval "\$entry_scr_func = sub { local \$_ = seq_align_index( \$_[0]->[2], \$_[1]->[2] );
                                   tr/$indices/$scores/;
                                   wantarray ? map { $scale * ( ord(\$_)-$scr_offset ) } split // : \$_;
                                 }";
    print STDERR $@ and return undef if $@;
    $self->{ _entry_scr_func } = $entry_scr_func;

    my $entry_nbs_func;
    eval "\$entry_nbs_func = sub { local \$_ = seq_align_index( \$_[0]->[2], \$_[1]->[2] );
                                   tr/\x{339}//d;
                                   return undef if ! length( \$_ );
                                   tr/$indices/$scores/;
                                   my \$ttl = 0;
                                   foreach ( split // ) { \$ttl += ord( \$_ ) }
                                   $scale * ( \$ttl / length( \$_ ) - $scr_offset );
                                 }";
    print STDERR $@ and return undef if $@;
    $self->{ _entry_nbs_func } = $entry_nbs_func;

    my $entry_pos_func;
    eval "\$entry_pos_func = sub { local \$_ = seq_align_index( \$_[0]->[2], \$_[1]->[2] );
                                   tr/\x{339}//d;
                                   return undef if ! length( \$_ );
                                   tr/$indices/$scores/;
                                   my \$ttl = 0;
                                   foreach ( split // ) { \$ttl++  if ord( \$_ ) > $scr_offset }
                                   \$ttl / length( \$_ );
                                 }";
    print STDERR $@ and return undef if $@;
    $self->{ _entry_pos_func } = $entry_pos_func;

    #  Scores from a BlastInterface HSP:

    my $hsp_scr_func;
    eval "\$hsp_scr_func = sub { local \$_ = seq_align_index( \$_[0]->[17], \$_[0]->[20] );
                                 tr/$indices/$scores/;
                                 wantarray ? map { $scale * ( ord(\$_)-$scr_offset ) } split // : \$_;
                               }";
    print STDERR $@ and return undef if $@;
    $self->{ _hsp_scr_func } = $hsp_scr_func;

    my $hsp_nbs_func;
    eval "\$hsp_nbs_func = sub { local \$_ = seq_align_index( \$_[0]->[17], \$_[0]->[20] );
                                 tr/\x{339}//d;
                                 return undef if ! length( \$_ );
                                 tr/$indices/$scores/;
                                 my \$ttl = 0;
                                 foreach ( split // ) { \$ttl += ord( \$_ ) }
                                 $scale * ( \$ttl / length( \$_ ) - $scr_offset );
                               }";
    print STDERR $@ and return undef if $@;
    $self->{ _hsp_nbs_func } = $hsp_nbs_func;

    my $hsp_pos_func;
    eval "\$hsp_pos_func = sub { local \$_ = seq_align_index( \$_[0]->[17], \$_[0]->[20] );
                                 tr/\x{339}//d;
                                 return undef if ! length( \$_ );
                                 tr/$indices/$scores/;
                                 my \$ttl = 0;
                                 foreach ( split // ) { \$ttl++  if ord( \$_ ) > $scr_offset }
                                 \$ttl / length( \$_ );
                               }";
    print STDERR $@ and return undef if $@;
    $self->{ _hsp_pos_func } = $hsp_pos_func;

    1;
}


#-------------------------------------------------------------------------------
#  Print the matrix in PAM matrix format
#
#       $self->print_matrix( )
#       $self->print_matrix( $file )
#       $self->print_matrix( \*FH )
#       $self->print_matrix( \$string )
#
#-------------------------------------------------------------------------------
sub print_matrix
{
    my $self = shift
        or print STDERR "AminoAcidMatrix::print_matrix: Improperly initialized object.\n"
            and return undef;
    my @aa_order = @{ $self->{ _aa_order } || [] }
        or print STDERR "AminoAcidMatrix::print_matrix: Improperly initialized object.\n"
            and return undef;
    my $mat = scoring_matrix( $self )
        or print STDERR "AminoAcidMatrix::print_matrix: Improperly initialized object.\n"
            and return undef;

    my ( $fh, $close ) = output_file_handle( $_[0] );
    $fh or print STDERR "AminoAcidMatrix::print_matrix: Problem creating output file handle.\n"
            and return undef;
    
    my $n_symb = @aa_order;

    foreach ( @{ $self->{ _comments } || [] } ) { print $fh "$_\n" }

    print $fh '  ', ( map { "   $_" } @aa_order ), "\n";
    for ( my $i = 1; $i <= $n_symb; $i++ )
    {
        my $res1 = $aa_order[$i-1];
        print $fh "$res1 ", ( map { sprintf ' %3d', $mat->[$i]->[$_] } (1 .. @$mat-1) ), "\n";
    }
}


#-------------------------------------------------------------------------------
#  Identify the canonical amino acids missing from a list or hash
#
#    $missing_as_text = missing_aa(  @residues )
#    $missing_as_text = missing_aa( \@residues )
#    $missing_as_text = missing_aa( \%res_keys )
#
#    @missing_as_list = missing_aa(  @residues )
#    @missing_as_list = missing_aa( \@residues )
#    @missing_as_list = missing_aa( \%res_keys )
#
#-------------------------------------------------------------------------------

sub missing_aa
{
    my %have = map { $canonical_aa{ $_ } ? ( uc $_ => 1 ) : () }
               ref( $_[0] ) eq 'ARRAY' ? @{$_[0]}      :
               ref( $_[0] ) eq 'HASH'  ? keys %{$_[0]} :
                                         @_;
    my @missing = grep { ! $have{ $_ } } @canonical_aa;
    wantarray ? @missing : list_as_text( @missing );
}


sub list_as_text
{
    local $_ = pop @_;
    ! defined( $_ ) ? '' : @_ ? join( ', ', @_ ) . " and $_" : $_;
}


#-------------------------------------------------------------------------------
#  Get an input file handle, and boolean on whether to close or not:
#
#  ( \*FH, $close ) = input_file_handle(  $filename );
#  ( \*FH, $close ) = input_file_handle( \*FH );
#  ( \*FH, $close ) = input_file_handle( );                   # D = STDIN
#  ( \*FH, $close ) = input_file_handle( \$text );
#
#-------------------------------------------------------------------------------

sub input_file_handle
{
    my ( $file ) = @_;

    my ( $fh, $close );
    if ( defined $file )
    {
        if ( ref $file eq 'GLOB' )
        {
            $fh = $file;
            $close = 0;
        }
        elsif ( ref $file eq 'SCALAR' )
        {
            open( $fh, "<", $file) || die "input_file_handle could not open scalar reference.\n";
            $close = 1;
        }
        elsif ( -f $file )
        {
            open( $fh, "<", $file) || die "input_file_handle could not open '$file'.\n";
            $close = 1;
        }
        else
        {
            die "input_file_handle could not find file '$file'.\n";
        }
    }
    else
    {
        $fh = \*STDIN;
        $close = 0;
    }

    wantarray ? ( $fh, $close ) : $fh;
}


#-------------------------------------------------------------------------------
#  Get an output file handle, and boolean on whether to close or not:
#
#  ( \*FH, $close ) = output_file_handle(  $filename );
#  ( \*FH, $close ) = output_file_handle( \*FH );
#  ( \*FH, $close ) = output_file_handle( );                   # D = STDOUT
#  ( \*FH, $close ) = output_file_handle( \$text );
#
#-------------------------------------------------------------------------------

sub output_file_handle
{
    my ( $file, $umask ) = @_;

    my ( $fh, $close );
    if ( defined $file )
    {
        if ( ref $file eq 'GLOB' )
        {
            $fh = $file;
            $close = 0;
        }
        elsif ( ref $file eq 'SCALAR' )
        {
            open( $fh, ">", $file) || die "output_file_handle could not open scalar reference.\n";
            $close = 1;
        }
        else
        {
            open( $fh, ">", $file) || die "output_file_handle could not open '$file'.\n";
            $umask ||= 0664;
            chmod $umask, $file;  #  Seems to work on open file!
            $close = 1;
        }
    }
    else
    {
        $fh = \*STDOUT;
        $close = 0;
    }

    wantarray ? ( $fh, $close ) : $fh;
}


1;
