#!/usr/bin/perl -w

# Inspiration: http://www.hmpdacc.org/doc/ReadProcessing_SOP.pdf
# And: http://bioinformatics.oxfordjournals.org/content/suppl/2011/12/12/btr669.DC1/SupplementaryFile1.pdf

use strict;
use File::Basename;
use File::Copy;
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
my $deconseq  = "/home/micro/sharptot/bin/deconseq.pl";
my $tmp_dir   = "/tmp/";
my $run_fastqc   = 0;
my $run_bmtagger = 0;
my $run_prinseq  = 0;
my $run_deconseq = 0;
my $cat_reads    = 0;
my $derep        = 0;
my $check_qc     = 0;
my $make_fasta   = 0;
my $compress     = 0;
my $paired_end   = 0;
my $entire_pipe  = 0;
my $overwrite    = 0; #overwrite the old results?
my $run_type;

GetOptions(
    "i=s"         => \$masterdir,
    "log-dir|l:s" => \$logdir,
    "nprocs=i"    => \$nprocs,
    "r=s"         => \$run_type,
    "compress"    => \$compress,
    #the following options are not needed and invoked internally
    #but advanced users can use them at command line to alter
    #specific behavior
    "fastqc"      => \$run_fastqc,
    "bmtagger"    => \$run_bmtagger,
    "deconseq"    => \$run_deconseq,
    "qc"          => \$run_prinseq,
    "cat"         => \$cat_reads,
    "derep"       => \$derep,
    "check-qc"    => \$check_qc,
    "make-fasta"  => \$make_fasta,
    "paired-end"  => \$paired_end,
    #this is obsolete. No longer does anything.
    "complete"    => \$entire_pipe, #note: does not invoke --compress
    );

my $settings = _set_settings( $run_type );
( $run_prinseq, $run_deconseq, $run_bmtagger,
  $cat_reads, $derep, $check_qc, $make_fasta ) = @{ $settings->{"parameters"} };
 
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

# Prinseq round 1: trim and filter low quality reads
if( $run_prinseq ){
    my $prinseq_trim_in_dir      = File::Spec->catdir( $masterdir, $settings->{"run_prinseq"}->{"input"} );
    my $prinseq_trim_results_dir = File::Spec->catdir( $masterdir, $settings->{"run_prinseq"}->{"output"} );
    make_path( $prinseq_trim_results_dir );
    my $prinseq_trim_log_dir        = File::Spec->catdir( $logdir, $settings->{"run_prinseq"}->{"log"} );
    make_path( $prinseq_trim_log_dir );
    
    print "PRINSEQ, TRIM&FILTER: " . scalar(localtime()) . "\n";
    _run_prinseq( {
	prinseq    => $prinseq,
	in_dir     => $prinseq_trim_in_dir,
	file_names => \@reads,
	log_dir    => $prinseq_trim_log_dir,
	result_dir => $prinseq_trim_results_dir,
	nprocs     => $nprocs,
	overwrite  => $overwrite,
	derep      => 0,
	paired_end => $paired_end,
		  });
}

# Deconseq
if( $run_deconseq ){
    my $deconseq_in_dir     = File::Spec->catdir( $masterdir, $settings->{"run_deconseq"}->{"input"} );
    my $deconseq_results_dir = File::Spec->catdir( $masterdir, $settings->{"run_deconseq"}->{"output"} );
    make_path( $deconseq_results_dir );
    my $deconseq_log_dir    = File::Spec->catdir( $logdir, $settings->{"run_deconseq"}->{"log"} );
    make_path( $deconseq_log_dir );
#note that decon_db_names MUST point to a key in the deconseq perl module
#where database locations are specified. A change the array below
#may require revising said perl module. You can probably find the module
#in some place like /src/deconseq-standalone-0.4.3/DeconSeqConfig.pm
    my @decon_db_names  = ( "mouse" ); 
    
    
    print "DECONSEQ: " . scalar(localtime()) . "\n";
    _run_deconseq( {
	deconseq    => $deconseq,
	in_dir      => $deconseq_in_dir, 
	file_names  => \@reads, 
	result_dir  => $deconseq_results_dir,
	log_dir     => $deconseq_log_dir, 
	nprocs      => $nprocs,
	tmp_dir     => $tmp_dir,
	db_names    => \@decon_db_names,
	overwrite   => $overwrite,
	paired_end  => $paired_end,
		   });
}


# BMTagger
if( $run_bmtagger ){
    my $bmtagger_in_dir      = File::Spec->catdir( $masterdir, $settings->{"run_bmtagger"}->{"input"});
    my $bmtagger_results_dir = File::Spec->catdir( $masterdir, $settings->{"run_bmtagger"}->{"output"});
    make_path( $bmtagger_results_dir );
    my $bmtagger_log_dir     = File::Spec->catdir( $logdir, $settings->{"run_bmtagger"}->{"log"} );
    make_path( $bmtagger_log_dir );
    my $bitmask_dir      = "/scratch/data/databases/mouse_genome_bmtagger/"; #files should end in .bitmask                                                                                       
    my $srprism_dir      = "/scratch/data/databases/mouse_genome_bmtagger/";
    my @db_names         = ( "Mus_musculus.NCBI.GRCm38.dna" ); #appropriate extensions (i.e., .fa, .bitmask) added below
    my $extract          = 1; #YOU PROBABLY WANT THIS, prints only non-host sequences in output file


    print "BMTAGGER: " . scalar(localtime()) . "\n";
    _run_bmtagger( {
	bmtagger    => $bmtagger,
	in_dir      => $bmtagger_in_dir, 
	file_names  => \@reads, 
	result_dir  => $bmtagger_results_dir, 
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

# Merge results
if( $cat_reads ){
    my $cat_reads_in_dir      = File::Spec->catdir( $masterdir, $settings->{"cat_reads"}->{"input"} );
    my $cat_reads_results_dir = File::Spec->catdir( $masterdir, $settings->{"cat_reads"}->{"output"});
    make_path( $cat_reads_results_dir );
    my $cat_reads_log = File::Spec->catdir( $logdir, $settings->{"cat_reads"}->{"log"} );
    make_path( $cat_reads_log );
        
    print "CATTING READS: " . scalar(localtime()) . "\n";
    _cat_reads( {
	in_dir     => $cat_reads_in_dir,
	file_names => \@reads,
	log_dir    => $cat_reads_log,
	result_dir => $cat_reads_results_dir,
	nprocs     => $nprocs,	
	paired_end => $paired_end,
		});
}

# Prinseq derep
if( $derep ){
    my $prinseq_derep_in_dir     = File::Spec->catdir( $masterdir, $settings->{"derep"}->{"input"} );
    my $prinseq_derep_results_dir = File::Spec->catdir( $masterdir, $settings->{"derep"}->{"output"} );
    make_path( $prinseq_derep_results_dir );
    my $prinseq_derep_log_dir        = File::Spec->catdir( $logdir, $settings->{"derep"}->{"log"});
    make_path( $prinseq_derep_log_dir );
        
    print "PRINSEQ, DEREP: " . scalar(localtime()) . "\n";
    _run_prinseq( {
	prinseq    => $prinseq,
	in_dir     => $prinseq_derep_in_dir,
	file_names => \@reads,
	log_dir    => $prinseq_derep_log_dir,
	result_dir => $prinseq_derep_results_dir,
	nprocs     => $nprocs,
	overwrite  => $overwrite,
	derep      => 1,
		  });
}


my @qc_reads = ();
if( $run_type eq "complete"){
    @qc_reads = @{ _get_fastq_file_names( File::Spec->catdir( $masterdir, $settings->{"derep"}->{"output"} ), 1  ) };
} elsif( $run_type eq "fast" ){
    @qc_reads = @{ _get_fastq_file_names( File::Spec->catdir( $masterdir, $settings->{"cat_reads"}->{"output"} ), 1 ) };
}

if( $check_qc ){
    my $fastq_clean_in_dir      = File::Spec->catdir( $masterdir, $settings->{"check_qc"}->{"input"}  );
    my $fastq_clean_results_dir = File::Spec->catdir( $masterdir, $settings->{"check_qc"}->{"output"});
    make_path( $fastq_clean_results_dir );
    my $fastq_clean_log_dir     = File::Spec->catdir( $logdir, $settings->{"check_qc"}->{"log"} );
    make_path( $fastq_clean_log_dir );

    print "FASTQC, CLEAN: " . scalar(localtime()) . "\n";
    _run_fastqc( {
	in_dir      => $fastq_clean_in_dir, 
	file_names  => \@qc_reads, 
	result_dir  => $fastq_clean_results_dir, 
	log_dir     => $fastq_clean_log_dir, 
	nprocs      => 1,
	fastqc      => $fastqc,
	is_clean_check => 1,
	paired_end  => $paired_end,
		 });
}

if( $make_fasta ){    
    my $fasta_in_dir      = File::Spec->catdir( $masterdir, $settings->{"make_fasta"}->{"input"}  );
    my $fasta_results_dir = File::Spec->catdir( $masterdir, $settings->{"make_fasta"}->{"output"} );
    make_path( $fasta_results_dir );
    my $fasta_log_dir     = File::Spec->catdir( $logdir, $settings->{"make_fasta"}->{"log"} );
    make_path( $fasta_log_dir );
    
    print "SEQRET: " . scalar(localtime()) . "\n";
    _run_seqret( {
	in_dir      => $fasta_in_dir,
	file_names  => \@qc_reads, 
	result_dir  => $fasta_results_dir, 
	log_dir     => $fasta_log_dir, 
	nprocs      => 1,
	seqret      => $seqret,
		 });
}

if( $compress ){
    print "COMPRESSING DATA: " . scalar(localtime()) . "\n";
    _compress_results( { 
	in_dir => File::Spec->catdir( $masterdir, $settings->{"make_fasta"}->{"output"}  )
		       });
    if( $run_type eq "complete" ){
	_compress_results( { 
	    in_dir => File::Spec->catdir( $masterdir, $settings->{"derep"}->{"output"}  )
			   });
    }
    _compress_results( { 
	in_dir => File::Spec->catdir( $masterdir, $settings->{"cat_reads"}->{"output"}  ),
	delete => 1,
		       });
    _compress_results( { 
	in_dir => File::Spec->catdir( $masterdir, $settings->{"run_prinseq"}->{"output"}  ),
	delete => 1,
		       });
    if( $run_type eq "complete" ){
	_compress_results( { 
	    in_dir => File::Spec->catdir( $masterdir, $settings->{"run_bmtagger"}->{"output"}  ),
	    delete => 1,
			   });
    }
    _compress_results( { 
	in_dir => File::Spec->catdir( $masterdir, $settings->{"run_deconseq"}->{"output"}  ),
	delete => 1,
		       });
}

print scalar(localtime()) . "\n";


########
#
# SUBROUTINES
#
#######

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
	#my $f_in  = File::Spec->catfile( $in_dir, $f_mate . ".fastq" );
	my $f_in  = File::Spec->catfile( $in_dir, $f_mate );
	my $f_log = File::Spec->catfile( $log_dir, $f_mate . ".log" );
	my @suffixlist = ( ".fastq.fastq", ".fastq", ".fastq.gz", ".fq", ".fq.gz" );
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
	    next if( -d $path );
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
	    next if( $read !~ m/\.fastq/ &&
		     $read !~ m/\.fq/    && 
		     $read !~ m/\_R1\_/ );
	} else {
	    next if( $read !~ m/\.fastq/ &&
		     $read !~ m/\.fq/    );
	} 
	print "Grabbing: $read\n";
	#my @suffixlist = ( ".fastq", ".fastq.gz", ".fq", ".fq.gz" );
	my @suffixlist = ( ".gz" );
	my( $name, $path, $suffix ) = fileparse( $read, @suffixlist );
	#my( $name, $path, $suffix ) = fileparse( $read ) ; 
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
	#$r_mate =~ s/\_R1\_/\_R2\_/;
	my( $f_in, $r_in, $f_log, $r_log );
	if( defined( $is_clean ) && $is_clean ){	    
	    #$f_in  = File::Spec->catfile( $in_dir, $f_mate . ".fastq" );
	    $f_in  = File::Spec->catfile( $in_dir, $f_mate );
	    #$r_in  = File::Spec->catfile( $in_dir, $r_mate . ".fastq" );	    
	    $r_in  = File::Spec->catfile( $in_dir, $r_mate);	    
	} else {
	    #$f_in  = File::Spec->catfile( $in_dir, $f_mate . ".fastq.gz" );
	    $f_in  = File::Spec->catfile( $in_dir, $f_mate );
	    #$r_in  = File::Spec->catfile( $in_dir, $r_mate . ".fastq.gz" );
	    $r_in  = File::Spec->catfile( $in_dir, $r_mate);
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
	#my $f_mate = $read . ".fastq.gz";
	my $f_mate = $read;
	#my $r_mate = $read . ".fastq.gz";
	my $r_mate = $read;
	#$r_mate =~ s/\_R1\_/\_R2\_/;
	my $f_in  = File::Spec->catfile( $in_dir, $f_mate  );
	my $r_in  = File::Spec->catfile( $in_dir, $r_mate );
	my $f_log = File::Spec->catfile( $log_dir, $f_mate . ".log" );
	my $r_log = File::Spec->catfile( $log_dir, $r_mate . ".log" );
	#prepare the output
        #might need to consider if we want a single output file or not.  
	my $out_stem = $read;
	#$out_stem =~ s/\_R1\_/\_RX\_/; #bmtagger appends it's own mate pair id, we no longer need this one
	#my $out_path = File::Spec->catfile( $result_dir, $out_stem . ".bmtagged" );
	my $out_path = File::Spec->catfile( $result_dir, $out_stem );
	#if( -e $out_path . "_1.fastq" && -e $out_path . "_2.fastq" ){
	#    $pm->finish unless( $overwrite );
	#}
	#loop over dbs and run bmtagger
	foreach my $db( @db_names ){
	    my $bitmask  = File::Spec->catfile( $bitmask_dir, $db . ".bitmask" );
	    my $srprism  = File::Spec->catfile( $sprism_dir, $db . ".srprism"  );
	    my $database = File::Spec->catfile( $sprism_dir, $db . ".fa" );
	    #my $f_string = "-1 " . $f_in;
	    my $f_string = $f_in;
	    #my $r_string = "-2 " . $r_in;
	    my $r_string = $r_in;
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

sub _run_deconseq{
    my( $args )  = @_;

    my $deconseq    = $args->{ "deconseq" };
    my $in_dir      = $args->{ "in_dir" };
    my @reads       = @{ $args->{ "file_names" } }; 
    my $result_dir  = $args->{ "result_dir" };
    my $log_dir     = $args->{"log_dir"}; 
    my $nprocs      = $args->{"nprocs"};
    my $tmp_dir     = $args->{"tmp_dir"};
    my @db_names    = @{ $args->{"db_names"} };
    my $overwrite   = $args->{"overwrite"};
    my $paired_end  = $args->{"paired_end"};

    my $db_string   = join(",", @db_names );

    #forward reads
    my $pm = Parallel::ForkManager->new($nprocs);
    for( my $i=1; $i<=$nprocs; $i++ ){
	my $pid = $pm->start and next;
	my $cmd;
	#do some housekeeping
	my $read = $reads[$i-1];
	#$read =~ s/\_R1\_/\_RX\_/; #bmtagger appends it's own mate pair id, we no longer need this one
	my ( $f_mate, $f_out_name );
	if( $paired_end ){
	    #$f_mate     = $read . ".bmtagged_1.fastq";
	    $f_mate     = $read;
	    #$f_out_name = $read . ".deconseq_1.fastq";
	    $f_out_name = $read;
	} else {
	    #$f_mate     = $read . ".bmtagged.fastq";
	    $f_mate     = $read;
	    #$f_out_name = $read . ".deconseq.fastq";
	    $f_out_name = $read;
	}
	my( $f_in, $r_in, $f_log, $r_log );
	$f_in  = File::Spec->catfile( $in_dir, $f_mate );
	$f_log = File::Spec->catfile( $log_dir, $f_mate . ".log" );
	#Start threads
	$cmd = "$deconseq -id $f_out_name -out_dir $result_dir -f $f_in -dbs $db_string";
	print "$cmd\n";
	system( $cmd );
	my $result_stem = File::Spec->catfile( $result_dir, $f_out_name );
	move( $result_stem . "_clean.fq", $result_stem ); #test.fq_clean.fq
	$pm->finish; # Terminates the child process
    }
    $pm->wait_all_children;    
    #reverse reads
    if( $paired_end ){
	my $pm2 = Parallel::ForkManager->new($nprocs);
	for( my $i=1; $i<=$nprocs; $i++ ){
	    my $pid = $pm2->start and next;
	    my $cmd;
	    #do some housekeeping
	    my $read = $reads[$i-1];
	    #$read =~ s/\_R1\_/\_RX\_/; #bmtagger appends it's own mate pair id, we no longer need this one
	    #my $r_mate = $read . ".bmtagged_2.fastq";
	    my $r_mate = $read;
	    #my $r_out_name = $read . ".deconseq_2.fastq";
	    my $r_out_name = $read;
	    #$r_mate =~ s/\_bmtagged_1\_/\_bmtagged_2\_/; #I don't think this is needed
	    my $r_in  = File::Spec->catfile( $in_dir, $r_mate );
	    my $r_log = File::Spec->catfile( $log_dir, $r_mate . ".log" );
	    #Start threads
	    $cmd = "$deconseq -id $r_out_name -out_dir $result_dir -f $r_in -dbs $db_string";
	    print "$cmd\n";
	    system( $cmd );
	    $pm2->finish; # Terminates the child process
	}
	$pm2->wait_all_children;    
    }
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
	    #do some housekeeping
	    my $read = $reads[$i-1];          
	    #$read =~ s/\_R1\_/\_RX\_/; #bmtagger appends it's own mate pair id, we no longer need this one
	    my ( $f_mate, $r_mate );
	    if( $paired_end ){
		#$f_mate = $read . ".deconseq_1.fastq_clean.fq";
		$f_mate = $read;
	        #$r_mate = $read . ".deconseq_2.fastq_clean.fq";
	        $r_mate = $read;
		#$r_mate =~ s/\_deconseq_1\_/\_deconseq_2\_/;
	    } else {
		#$f_mate = $read . ".deconseq.fastq";
	        #$r_mate = $read . ".deconseq.fastq"; #we don't actually call this below
		$f_mate = $read;
	    }
	    my $f_in  = File::Spec->catfile( $in_dir, $f_mate );
	    #my $r_in  = File::Spec->catfile( $in_dir, $r_mate );
	    my $r_in;
	    my $f_log = File::Spec->catfile( $log_dir, $f_mate . ".log" );
	    #my $r_log = File::Spec->catfile( $log_dir, $r_mate . ".log" );
	    my $r_log;
	    #prepare the output
	    #might need to consider if we want a single output file or not.  
	    #my $out_path = File::Spec->catfile( $result_dir, $read . ".trim_filter" );
	    my $out_path = File::Spec->catfile( $result_dir, $read );
	    #if( -e $out_path . "_1.fastq" && -e $out_path . "_2.fastq" ){
		#next unless( $overwrite );
	    #}
	    my $compressed = 0;	   
	    my $gz_file = File::Spec->catfile( $f_in . ".gz" );
	    if( -e $gz_file ){
		$compressed = 1;
	    }
	    my $cmd = "";
	    if( $compressed ){
		$cmd .= "zcat ${f_in}.gz | ";
	    }
	    #you might want to turn derep back on here for downstream efficiency, but we turned off for
	    #courtney's analysis
	    #$cmd =  "$prinseq -verbose -derep 14 -derep_min 2 -no_qual_header "; #do we want -exact_only?
	    $cmd .=  "$prinseq -verbose -no_qual_header "; #do we want -exact_only?
	    $cmd .= "-min_len 60 -max_len 200 -min_qual_mean 25 -ns_max_n 0 ";
	    $cmd .= "-lc_method entropy -lc_threshold 60 -trim_qual_left 20 -trim_qual_right 20 ";
	    if( $paired_end ){
		$cmd .= "-out_good $out_path -fastq $f_in -fastq2 $r_in -log $f_log ";
	    } else {
		if( $compressed ){
		    $cmd .= "-fastq stdin ";
		} else {
		    $cmd .= "-fastq $f_in ";
		}
		$cmd .= "-out_good $out_path -log $f_log ";
	    }
	    $cmd .= "-out_bad null ";
	    print "$cmd\n";
	    system( $cmd );
	    #trim the name of the output file to standard to ease next steps. Use a move
	    move( $out_path . ".fastq", $out_path );
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
    #@reads      = map { s/\_R1\_/\_RX\_/; $_ } @reads;
    my $log_dir = $args->{"log_dir"};
    my $result_dir = $args->{"result_dir"};
    my @sorted     = sort( @reads );    
    @sorted        = map { $in_dir . "/" . $_ } @sorted;
    my @f_sorted   = ();
    my @r_sorted   = ();
    if( $paired_end ) {
	#@f_sorted   = map { $_ . ".trim_filter_1.fastq" } @sorted;	
	@f_sorted    = @sorted;
	#@r_sorted   = map { $_ . ".trim_filter_2.fastq" } @sorted;
	@r_sorted   = @sorted;
    } else {
	#@f_sorted   = map { $_ . ".trim_filter.fastq" } @sorted;
	@f_sorted    = @sorted;
	#@r_sorted   = map { $_ . ".trim_filter.fastq" } @sorted; #this never actually gets used
	@r_sorted    = @sorted
    }
    my $out_stem   = $reads[0];
    #$out_stem =~ s/\_RX.*$/\.fastq/; #bmtagger appends it's own mate pair id, we no longer need this one
    my $out_path = File::Spec->catfile( $result_dir, $out_stem );
    if( $paired_end ){
	print("cat @f_sorted @r_sorted > $out_path\n");
	system("cat @f_sorted @r_sorted > $out_path");
    } else {
	print("cat @f_sorted > $out_path\n");
	system("cat @f_sorted > $out_path");
    }
}

# settings is a hashref with the following pointers:
# setting->method->data_type
# where method is a function name (e.g., bmtagger)
# data_type is input, output, log, or something custom to function
# for now, there are only two settings configured: fast and complete
sub _set_settings{
    my $setting_type = shift;   
    my ( $run_prinseq, $run_deconseq, $run_bmtagger,
	 $cat_reads, $derep, $check_qc, $make_fasta );
    if( $setting_type eq "fast" ){
	$run_prinseq  = 1;
	$run_deconseq = 1;
	$run_bmtagger = 0;
	$cat_reads    = 1;
	$derep        = 0;
	$check_qc     = 1;
	$make_fasta   = 1;
    } elsif( $setting_type eq "complete" ){
	$run_prinseq  = 1;
	$run_deconseq = 1;
	$run_bmtagger = 1;
	$cat_reads    = 1;
	$derep        = 1;
	$check_qc     = 1;
	$make_fasta   = 1;
    } else{
	die( "I don't know how to process -r $setting_type\n") ;
    }
    my $settings     = _build_settings( $setting_type );
    my @params = ( $run_prinseq, $run_deconseq, $run_bmtagger,
		   $cat_reads, $derep, $check_qc, $make_fasta );
    $settings->{"parameters"} = \@params;
    return $settings;
}

sub _build_settings{
    my $setting_type = shift;
    my $settings     = ();
    if( $setting_type eq "fast" ){
	my @keys = qw( run_prinseq run_deconseq cat_reads check_qc make_fasta );
	foreach my $key( @keys ){
	    if( $key eq "run_prinseq" ){
		$settings->{$key}->{"input"}  = "";
		$settings->{$key}->{"output"} = "prinseq_trim";
		$settings->{$key}->{"log"}    = "prinseq_log";
	    } elsif( $key eq "run_deconseq" ){
		$settings->{$key}->{"input"}  = "prinseq_trim";	       
		$settings->{$key}->{"output"} = "deconseq";
		$settings->{$key}->{"log"}    = "deconseq_log";
	    } elsif( $key eq "cat_reads" ){
		$settings->{$key}->{"input"}  = "deconseq";  
		$settings->{$key}->{"output"} = "cat_reads";
		$settings->{$key}->{"log"}    = "cat_reads_log";
	    } elsif( $key eq "check_qc" ){
		$settings->{$key}->{"input"}  = "cat_reads"; 
		$settings->{$key}->{"output"} = "fastqc_clean";     
		$settings->{$key}->{"log"}    = "fastqc_clean_log";
	    } elsif( $key eq "make_fasta" ){
		$settings->{$key}->{"input"}  = "cat_reads";
		$settings->{$key}->{"output"} = "fasta_clean";     
		$settings->{$key}->{"log"}    = "fasta_log";
	    } else {
		die( "I don't know how to process the setting key $key\n" );
	    }
	}
    } elsif( $setting_type eq "complete" ){
	my @keys = qw( run_prinseq run_deconseq run_bmtagger cat_reads derep check_qc make_fasta );
	foreach my $key( @keys ){
	    if( $key eq "run_prinseq" ){
		$settings->{$key}->{"input"}  = "";
		$settings->{$key}->{"output"} = "prinseq_trim";
		$settings->{$key}->{"log"}    = "prinseq_trim_log";
	    } elsif( $key eq "run_deconseq" ){
		$settings->{$key}->{"input"}  = "prinseq_trim";	       
		$settings->{$key}->{"output"} = "deconseq";
		$settings->{$key}->{"log"}    = "deconseq_log";
	    } elsif( $key eq "run_bmtagger" ){
		$settings->{$key}->{"input"}  = "deconseq";  
		$settings->{$key}->{"output"} = "bmtagger";
		$settings->{$key}->{"log"}    = "bmtagger_log";
	    } elsif( $key eq "cat_reads" ){
		$settings->{$key}->{"input"}  = "bmtagger";  
		$settings->{$key}->{"output"} = "cat_reads";
		$settings->{$key}->{"log"}    = "cat_reads_log";
	    } elsif( $key eq "derep" ){
		$settings->{$key}->{"input"}  = "cat_reads";  
		$settings->{$key}->{"output"} = "prinseq_derep";
		$settings->{$key}->{"log"}    = "prinseq_derep_log";
	    } elsif( $key eq "check_qc" ){
		$settings->{$key}->{"input"}  = "prinseq_derep"; 
		$settings->{$key}->{"output"} = "fastqc_clean";     
		$settings->{$key}->{"log"}    = "fastqc_clean_log";
	    } elsif( $key eq "make_fasta" ){
		$settings->{$key}->{"input"}  = "cat_reads";
		$settings->{$key}->{"output"} = "fasta_clean";     
		$settings->{$key}->{"log"}    = "fasta_log";
	    } else {
		die( "I don't know how to process the setting key $key\n" );
	    }
	}
    }
    return $settings;
}
