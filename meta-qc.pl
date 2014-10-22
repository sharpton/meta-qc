#!/usr/bin/perl -w

# Inspiration: http://www.hmpdacc.org/doc/ReadProcessing_SOP.pdf
# And: http://bioinformatics.oxfordjournals.org/content/suppl/2011/12/12/btr669.DC1/SupplementaryFile1.pdf

use strict;
use File::Basename;
use File::Path qw( make_path remove_tree );
use File::Spec;
use File::Spec::Functions;
use Parallel::ForkManager;
use Carp;
use Data::Dumper;
use Getopt::Long;

my $masterdir = "/home/micro/sharptot/projects/ibdmouse/testdata";
my $logdir;
my $nprocs    = 80;
my $bmtagger  = "/home/micro/sharptot/bin/bmtagger.sh";
my $fastqc    = "/home/micro/sharptot/bin/fastqc";
my $prinseq   = "/home/micro/sharptot/bin/prinseq-lite.pl";
my $seqret    = "/local/cluster/bin/seqret";
my $run_fastqc = 0;
my $run_bmtagger = 0;
my $run_prinseq  = 0;
my $cat_reads    = 0;
my $derep        = 0;
my $check_qc     = 0;
my $make_fasta   = 0;
my $compress     = 0;
my $paired_end   = 0;

GetOptions(
    "i=s" => \$masterdir,
    "log-dir:s"   => \$logdir,
    "nprocs=i"    => \$nprocs,
    "fastqc"      => \$run_fastqc,
    "host-filter" => \$run_bmtagger,
    "qc"          => \$run_prinseq,
    "cat"         => \$cat_reads,
    "derep"       => \$derep,
    "check-qc"    => \$check_qc,
    "make-fasta"  => \$make_fasta,
    "compress"    => \$compress,
    "paired-end"  => \$paired_end,
    );

if( ! defined( $logdir ) ){
    $logdir    = $masterdir . "/logs/";
}

#initialize log dir
make_path( $logdir );

#get a list of fastq files from masterdir
my @reads = @{ _get_fastq_file_names( $masterdir, 0 ) };

#set nprocs
if( scalar( @reads ) < $nprocs ){
    $nprocs = scalar( @reads );
}

print scalar(localtime()) . "\n";

# FastQC
my $fastq_result_dir = File::Spec->catdir( $masterdir, "/fastqc_raw/" );
make_path( $fastq_result_dir );
my $fastq_log_dir    = File::Spec->catdir( $logdir, "/fastqc_raw/" );
make_path( $fastq_log_dir );
if( $run_fastqc ){
    print "FASTQC: " . scalar(localtime()) . "\n";
    _run_fastqc( {
	in_dir      => $masterdir, 
	file_names  => \@reads, 
	result_dir  => $fastq_result_dir, 
	log_dir     => $fastq_log_dir, 
	nprocs      => $nprocs,
	fastqc      => $fastqc,
	paired_end  => $paired_end,
		 });
}

# BMTagger
my $bmtagger_result_dir = File::Spec->catdir( $masterdir, "/bmtagger/" );
make_path( $bmtagger_result_dir );
my $bmtagger_log_dir    = File::Spec->catdir( $logdir, "/bmtagger/" );
make_path( $bmtagger_log_dir );
my $bitmask_dir      = "/scratch/data/databases/mouse_genome_bmtagger/"; #files should end in .bitmask                                                                                       
my $srprism_dir      = "/scratch/data/databases/mouse_genome_bmtagger/";
my $tmp_dir          = "/tmp/";
my @db_names         = ( "Mus_musculus.NCBI.GRCm38.dna" ); #appropriate extensions (i.e., .fa, .bitmask) added below
my $extract          = 1; #YOU PROBABLY WANT THIS, prints only non-host sequences in output file
my $overwrite        = 0;
if( $run_bmtagger ){
    print "BMTAGGER: " . scalar(localtime()) . "\n";
    _run_bmtagger( {
	bmtagger    => $bmtagger,
	in_dir      => $masterdir, 
	file_names  => \@reads, 
	result_dir  => $bmtagger_result_dir, 
	log_dir     => $bmtagger_log_dir, 
	nprocs      => $nprocs,
	bitmask_dir => $bitmask_dir,
	srprism_dir => $srprism_dir,
	tmp_dir     => $tmp_dir,
	db_names    => \@db_names,
	extract     => $extract,
	overwrite   => $overwrite,
	paired_end  => $paired_end,
		   });
}

# Prinseq round 1: trim and filter low quality reads
my $prinseq_trim_result_dir = File::Spec->catdir( $masterdir, "/prinseq_trim/" );
make_path( $prinseq_trim_result_dir );
my $prinseq_trim_log_dir        = File::Spec->catdir( $logdir, "/prinseq_trim/" );
make_path( $prinseq_trim_log_dir );

if( $run_prinseq ){
    print "PRINSEQ, TRIM&FILTER: " . scalar(localtime()) . "\n";
    _run_prinseq( {
	prinseq    => $prinseq,
	in_dir     => $bmtagger_result_dir,
	file_names => \@reads,
	log_dir    => $prinseq_trim_log_dir,
	result_dir => $prinseq_trim_result_dir,
	nprocs     => $nprocs,
	overwrite  => $overwrite,
	derep      => 0,
	paired_end => $paired_end,
		  });
}

# Merge results
my $cat_reads_dir = File::Spec->catdir( $masterdir, "/catted_reads/" );
make_path( $cat_reads_dir );
my $cat_reads_log = File::Spec->catdir( $logdir, "/catted_reads/" );
make_path( $cat_reads_log );

if( $cat_reads ){
    print "CATTING READS: " . scalar(localtime()) . "\n";
    _cat_reads( {
	in_dir     => $prinseq_trim_result_dir,
	file_names => \@reads,
	log_dir    => $cat_reads_log,
	result_dir => $cat_reads_dir,
	nprocs     => $nprocs,	
	paired_end => $paired_end,
		});
}

# Prinseq derep
my $prinseq_derep_result_dir = File::Spec->catdir( $masterdir, "/prinseq_derep/" );
make_path( $prinseq_derep_result_dir );
my $prinseq_derep_log_dir        = File::Spec->catdir( $logdir, "/prinseq_derep/" );
make_path( $prinseq_derep_log_dir );

if( $derep ){
    print "PRINSEQ, DEREP: " . scalar(localtime()) . "\n";
    _run_prinseq( {
	prinseq    => $prinseq,
	in_dir     => $cat_reads_dir,
	file_names => \@reads,
	log_dir    => $prinseq_derep_log_dir,
	result_dir => $prinseq_derep_result_dir,
	nprocs     => $nprocs,
	overwrite  => $overwrite,
	derep      => 1,
		  });
}

my @qc_reads = @{ _get_fastq_file_names( $prinseq_derep_result_dir, 1 ) };
if( $check_qc ){
    my $fastq_clean_result_dir = File::Spec->catdir( $masterdir, "/fastqc_clean/" );
    make_path( $fastq_clean_result_dir );
    my $fastq_clean_log_dir    = File::Spec->catdir( $logdir, "/fastqc_clean/" );
    make_path( $fastq_clean_log_dir );

    print "FASTQC, CLEAN: " . scalar(localtime()) . "\n";
    _run_fastqc( {
	in_dir      => $prinseq_derep_result_dir, 
	file_names  => \@qc_reads, 
	result_dir  => $fastq_clean_result_dir, 
	log_dir     => $fastq_clean_log_dir, 
	nprocs      => 1,
	fastqc      => $fastqc,
	is_clean_check => 1,
	paired_end  => $paired_end,
		 });
}

my $fasta_result_dir = File::Spec->catdir( $masterdir, "/fasta_clean/" );
make_path( $fasta_result_dir );
my $fasta_log_dir    = File::Spec->catdir( $logdir, "/fasta_clean/" );
make_path( $fasta_log_dir );

if( $make_fasta ){    
    print "SEQRET: " . scalar(localtime()) . "\n";
    _run_seqret( {
	in_dir      => $prinseq_derep_result_dir, 
	file_names  => \@qc_reads, 
	result_dir  => $fasta_result_dir, 
	log_dir     => $fasta_log_dir, 
	nprocs      => 1,
	seqret      => $seqret,
		 });
}

#need to build this function. Also, might just want to remove, not compress, intermediate results...
if( $compress ){
    print "COMPRESSING DATA: " . scalar(localtime()) . "\n";
    _compress_results( { in_dir => $fasta_result_dir } );
    _compress_results( { in_dir => $prinseq_derep_result_dir } );
    _compress_results( { 
	in_dir => $cat_reads_dir,
	delete => 1,
		       });
    _compress_results( { 
	in_dir => $prinseq_trim_result_dir, 
	delete => 1,
		       });
    _compress_results( { 
	in_dir => $bmtagger_result_dir,
	delete => 1,
		       });
}

print scalar(localtime()) . "\n";

sub _run_seqret{
    my $args = shift;
    
    my $in_dir = $args->{"in_dir"};
    my @reads  = @{ $args->{"file_names"} };
    my $result_dir = $args->{"result_dir"};
    my $log_dir    = $args->{"log_dir"};
    my $nprocs     = $args->{"nprocs"};
    my $seqret     = $args->{"seqret"};

    #forward reads
    my $pm = Parallel::ForkManager->new($nprocs);
    for( my $i=1; $i<=$nprocs; $i++ ){
	my $pid = $pm->start and next;
	my $cmd;
	#do some housekeeping
	my $read = $reads[$i-1];
	my $f_mate = $read;
	my $f_in  = File::Spec->catfile( $in_dir, $f_mate . ".fastq" );
	my $f_log = File::Spec->catfile( $log_dir, $f_mate . ".log" );
	my @suffixlist = ( ".fastq.fastq", ".fastq", ".fastq.gz" );
	my ($name,$path,$suffix) = fileparse($f_in ,@suffixlist);
	my $outseq = File::Spec->catfile( $result_dir, $name . ".fa" );	
	#Start threads
	$cmd = "$seqret -sequence $f_in -outseq $outseq -sformat1 fastq -osformat2 fasta &> $f_log";
	print "$cmd\n";
	system( $cmd );
	$pm->finish; # Terminates the child process
    }
    $pm->wait_all_children;    
}


sub _compress_results{
    my $args = shift;
    
    my $in_dir = $args->{"in_dir"};
    my $delete = $args->{"delete"};

    opendir( DIR, $in_dir ) || die "Can't opendir $in_dir for read: $!\n";
    my @files = readdir( DIR );
    closedir DIR;
    foreach my $file( @files ){
	next if ( $file =~ m/^\./ );
	my $path = $in_dir . "/" . $file;
	if( defined( $delete ) && $delete ){
	    print "Deleting $path\n";
	    unlink( $path )
	}
	else{
	    print "Compressing $path\n";
	    system( "gzip $path" );
	}
    }
}

sub _get_fastq_file_names{
    my $masterdir = shift;
    my $is_clean_check = shift; #are we running this on dereped data rather than raw reads?
    opendir( MASTER, $masterdir) || die "Can't open $masterdir for read: $!\n";
    my @files = readdir( MASTER );
    closedir MASTER;
    my @reads = (); 
    foreach my $read( @files ){
	#only grab the R1 reads; will add R2 later
	if( !$is_clean_check ){
	    next if( $read !~ m/\.fastq/ || $read !~ m/\_R1\_/ );
	} else {
	    next if( $read !~ m/\.fastq/ );
	} 
	my @suffixlist = ( ".fastq", ".fastq.gz" );
	my( $name, $path, $suffix ) = fileparse( $read, @suffixlist );
	push( @reads, $name );
    }
    return \@reads;
}

sub _run_fastqc{
    my( $args ) = @_;
 
    #args->{in_dir} is directory containing input files
    #args->{file_names} is arrayref of filenames in in_dir to process
    #args->{result_dir} is output directory
    #args->{log_dir} is directory containing log files for each task
    #args->{nprocs} is number of tasks to run in parallel
    
    my $in_dir     = $args->{"in_dir"};
    my @reads      = @{ $args->{"file_names"} };
    my $result_dir = $args->{"result_dir"};
    my $log_dir    = $args->{"log_dir"};
    my $nprocs     = $args->{"nprocs"};
    my $fastqc     = $args->{"fastqc"};
    my $is_clean   = $args->{"is_clean_check"};
    my $paired_end = $args->{"paired_end"};

    #forward reads
    my $pm = Parallel::ForkManager->new($nprocs);
    for( my $i=1; $i<=$nprocs; $i++ ){
	my $pid = $pm->start and next;
	my $cmd;
	#do some housekeeping
	my $read = $reads[$i-1];
	my $f_mate = $read;
	my $r_mate = $read;
	$r_mate =~ s/\_R1\_/\_R2\_/;
	my( $f_in, $r_in, $f_log, $r_log );
	if( defined( $is_clean ) && $is_clean ){
	    $f_in  = File::Spec->catfile( $in_dir, $f_mate . ".fastq" );
	    $r_in  = File::Spec->catfile( $in_dir, $r_mate . ".fastq" );	    
	} else {
	    $f_in  = File::Spec->catfile( $in_dir, $f_mate . ".fastq.gz" );
	    $r_in  = File::Spec->catfile( $in_dir, $r_mate . ".fastq.gz" );
	}
	$f_log = File::Spec->catfile( $log_dir, $f_mate . ".log" );
	$r_log = File::Spec->catfile( $log_dir, $r_mate . ".log" );
	#Start threads
	$cmd = "$fastqc -o=$result_dir $f_in > $f_log 2>&1";
	print "$cmd\n";
	system( $cmd );
	$pm->finish; # Terminates the child process
    }
    $pm->wait_all_children;    
    #reverse reads
    return if ( defined( $is_clean ) && $is_clean );
    if( $paired_end ){
	my $pm2 = Parallel::ForkManager->new($nprocs);
	for( my $i=1; $i<=$nprocs; $i++ ){
	    my $pid = $pm2->start and next;
	    my $cmd;
	    #do some housekeeping
	    my $read = $reads[$i-1];
	    my $f_mate = $read;
	    my $r_mate = $read;
	    $r_mate =~ s/\_R1\_/\_R2\_/;
	    my $f_in  = File::Spec->catfile( $in_dir, $f_mate . ".fastq.gz" );
	    my $r_in  = File::Spec->catfile( $in_dir, $r_mate . ".fastq.gz" );
	    my $f_log = File::Spec->catfile( $log_dir, $f_mate . ".log" );
	    my $r_log = File::Spec->catfile( $log_dir, $r_mate . ".log" );
	    #Start threads
	    $cmd = "$fastqc -o=$result_dir $r_in > $r_log 2>&1";
	    print "$cmd\n";
	    system( $cmd );
	    $pm2->finish; # Terminates the child process
	}
	$pm2->wait_all_children;    
    }
}

sub _run_bmtagger{
    my( $args ) = @_;

    my $in_dir      = $args->{"in_dir"};
    my @reads       = @{ $args->{"file_names"} };
    my $result_dir  = $args->{"result_dir"};
    my $log_dir     = $args->{"log_dir"};
    my $nprocs      = $args->{"nprocs"};
    my $bitmask_dir = $args->{"bitmask_dir"};
    my $sprism_dir  = $args->{"srprism_dir"};
    my @db_names    = @{ $args->{"db_names"} };
    my $extract     = $args->{"extract"};
    my $overwrite   = $args->{"overwrite"};
    my $tmp_dir     = $args->{"tmp_dir"};
    my $bmtagger    = $args->{"bmtagger"};
    my $paired_end  = $args->{"paired_end"};

    my $pm = Parallel::ForkManager->new($nprocs);
    for( my $i=1; $i<=$nprocs; $i++ ){
	my $pid = $pm->start and next;
	my $cmd;
	#do some housekeeping
	my $read = $reads[$i-1];          
	my $f_mate = $read . ".fastq.gz";
	my $r_mate = $read . ".fastq.gz";
	$r_mate =~ s/\_R1\_/\_R2\_/;
	my $f_in  = File::Spec->catfile( $in_dir, $f_mate  );
	my $r_in  = File::Spec->catfile( $in_dir, $r_mate );
	my $f_log = File::Spec->catfile( $log_dir, $f_mate . ".log" );
	my $r_log = File::Spec->catfile( $log_dir, $r_mate . ".log" );
	#prepare the output
        #might need to consider if we want a single output file or not.  
	my $out_stem = $read;
	$out_stem =~ s/\_R1\_/\_RX\_/; #bmtagger appends it's own mate pair id, we no longer need this one
	my $out_path = File::Spec->catfile( $result_dir, $out_stem . ".bmtagged" );
	#if( -e $out_path . "_1.fastq" && -e $out_path . "_2.fastq" ){
	#    $pm->finish unless( $overwrite );
	#}
	#loop over dbs and run bmtagger
	foreach my $db( @db_names ){
	    my $bitmask  = File::Spec->catfile( $bitmask_dir, $db . ".bitmask" );
	    my $srprism  = File::Spec->catfile( $sprism_dir, $db . ".srprism"  );
	    my $database = File::Spec->catfile( $sprism_dir, $db . ".fa" );
	    my $f_string = "-1 " . $f_in;
	    my $r_string = "-2 " . $r_in;
	    if( $extract ){
		#bmtagger.sh -X -b $BITMASK -x $SPRISM -T $TMP -q1 $FREAD $RREAD -o $OUTPUT  >> $LOGS/bmtagger/${JOB_ID}.all 2>&1
#		$cmd = "run_bmtagger.sh -X -b $bitmask -x $sprism -T $tmp_dir -q1 $f_string $r_string -o $out_path > $f_log 2>&1";
		$cmd =  "$bmtagger ";
		$cmd .= "-X ";
		#$cmd .= "-b $bitmask -x $srprism -T $tmp_dir -d $database ";
		$cmd .= "-b $bitmask -x $srprism -d $database ";
		if( $paired_end ){
		    $cmd .= "-q 1 $f_string $r_string -o $out_path"; 
		} else {
		    $cmd .= "-q 1 $f_string -o $out_path" ;
		}
	    }
	    else{
                #bmtagger.sh -b $BITMASK -x $SPRISM -T $TMP -q1 $FREAD $RREAD -o $OUTPUT  >> $LOGS/bmtagger/${JOB_ID}.all 2>&1
#		$cmd = "run_bmtagger.sh -b $bitmask -x $sprism -T $tmp_dir -q1 $f_string $r_string -o $out_path > $f_log 2>&1";
		$cmd =  "$bmtagger ";
		$cmd .= "-b $bitmask -x $srprism -T $tmp_dir -d $database ";
		if( $paired_end ){
		    $cmd .= "-q 1 $f_string $r_string -o $out_path"; 
		} else {
		    $cmd .= "-q 1 $f_string -o $out_path"; 
		}
	    }
	    $cmd .= " &> $f_log";
	    print "$cmd\n";
	    system( $cmd );	    
	}
	$pm->finish; # Terminates the child process
    }
    $pm->wait_all_children;    
}

sub _run_prinseq{
    my( $args ) = @_;

    my $in_dir      = $args->{"in_dir"};
    my @reads       = @{ $args->{"file_names"} };
    my $result_dir  = $args->{"result_dir"};
    my $log_dir     = $args->{"log_dir"};
    my $nprocs      = $args->{"nprocs"};
    my $overwrite   = $args->{"overwrite"};
    my $prinseq     = $args->{"prinseq"};
    my $derep       = $args->{"derep"};
    my $paired_end  = $args->{"paired_end"};

    if( ! $derep ){
	my $pm = Parallel::ForkManager->new($nprocs);
	for( my $i=1; $i<=$nprocs; $i++ ){
	    my $pid = $pm->start and next;
	    my $cmd;
	    #do some housekeeping
	    my $read = $reads[$i-1];          
	    $read =~ s/\_R1\_/\_RX\_/; #bmtagger appends it's own mate pair id, we no longer need this one
	    my $f_mate = $read . ".bmtagged_1.fastq";
	    my $r_mate = $read . ".bmtagged_2.fastq";
	    $r_mate =~ s/\_bmtagged_1\_/\_bmtagged_2\_/;
	    my $f_in  = File::Spec->catfile( $in_dir, $f_mate );
	    my $r_in  = File::Spec->catfile( $in_dir, $r_mate );
	    my $f_log = File::Spec->catfile( $log_dir, $f_mate . ".log" );
	    my $r_log = File::Spec->catfile( $log_dir, $r_mate . ".log" );
	    #prepare the output
	    #might need to consider if we want a single output file or not.  
	    my $out_path = File::Spec->catfile( $result_dir, $read . ".trim_filter" );
	    #if( -e $out_path . "_1.fastq" && -e $out_path . "_2.fastq" ){
		#next unless( $overwrite );
	    #}
	    $cmd =  "$prinseq -verbose -derep 14 -derep_min 2 -no_qual_header "; #do we want -exact_only?
	    $cmd .= "-min_len 60 -max_len 200 -min_qual_mean 25 -ns_max_n 0 ";
	    $cmd .= "-lc_method entropy -lc_threshold 60 -trim_qual_left 20 -trim_qual_right 20 ";
	    if( $paired_end ){
		$cmd .= "-out_good $out_path -fastq $f_in -fastq2 $r_in -log $f_log ";
	    } else {
		$cmd .= "-out_good $out_path -fastq $f_in -log $f_log ";
	    }
	    $cmd .= "-out_bad null ";
	    print "$cmd\n";
	    system( $cmd );
	    $pm->finish; # Terminates the child process
	}
	$pm->wait_all_children;    
    } else {
	#only working with a single file from catted_reads dir
	my $in_file   = $reads[0];
	$in_file =~ s/\_R1.*$/\.fastq/; #R1, not RX, because this is read, not catted file.
	my $in_path = File::Spec->catfile( $in_dir, $in_file  );
	my $log_path = File::Spec->catfile( $log_dir, $reads[0] . ".log" );
	my $out_file  = $reads[0]; #again, R1, not RX
	$out_file =~ s/\_R1.*$/\.cleaned\.fastq/; #bmtagger appends it's own mate pair id, we no longer need this one
	my $out_path = File::Spec->catfile( $result_dir, $out_file );
	my $cmd =  "$prinseq -verbose -derep 14 -derep_min 2 -no_qual_header "; #do we want -exact_only?
	$cmd    .= "-out_good $out_path -fastq $in_path -log $log_path ";
	$cmd    .= "-out_bad null ";
	print "$cmd\n";
	system( $cmd );	    
    }
}

#Currently ignores sequences in the prinseq singletons file! This is the conservative move, it would seem.
#Also, F and R reads are pushed into the same file.
sub _cat_reads{
    my( $args ) = shift;
    my $in_dir  = $args->{"in_dir"};
    my $paired_end = $args->{"paired_end"};
    my @reads   = @{ $args->{"file_names"} };
    @reads      = map { s/\_R1\_/\_RX\_/; $_ } @reads;
    my $log_dir = $args->{"log_dir"};
    my $result_dir = $args->{"result_dir"};
    my @sorted     = sort( @reads );    
    @sorted        = map { $in_dir . "/" . $_ } @sorted;
    my @f_sorted   = map { $_ . ".trim_filter_1.fastq" } @sorted;
    my @r_sorted   = map { $_ . ".trim_filter_2.fastq" } @sorted;
    my $out_stem   = $reads[0];
    $out_stem =~ s/\_RX.*$/\.fastq/; #bmtagger appends it's own mate pair id, we no longer need this one
    my $out_path = File::Spec->catfile( $result_dir, $out_stem );
    if( $paired_end ){
	print("cat @f_sorted @r_sorted > $out_path\n");
	system("cat @f_sorted @r_sorted > $out_path");
    } else {
	print("cat @f_sorted > $out_path\n");
	system("cat @f_sorted > $out_path");
    }
}
