#!/usr/bin/perl

use strict;
my $MAXNJOBS = 12;

my $STARTJOB = 0;
my $ENDJOB = 800;

my $SCRIPTNAME = "CourtshipBowlTrack";

my $nargs = $#ARGV + 1;
if($nargs < 1){
    print "Usage: qsub_$SCRIPTNAME.pl <expdirlist.txt>\n";
    exit(1);
}

my $KEYWORD = "WING_GO";

my $ANALYSIS_PROTOCOL = "20130212_galit_rejection_wingclipped_20130207";

my $BASEDIR = "/groups/branson/bransonlab/projects/CourtshipBowls/CourtshipBowlAnalysis";

if($nargs >= 2){
    $ANALYSIS_PROTOCOL = $ARGV[1];
}

my $SCRIPT = "$BASEDIR/$SCRIPTNAME/distrib/run_$SCRIPTNAME.sh";

my $PARAMS = "analysis_protocol $ANALYSIS_PROTOCOL";

my $MCR = "/groups/branson/bransonlab/projects/olympiad/MCR/v717";

my $TMP_ROOT_DIR = "/tmp";
my $MCR_CACHE_ROOT = "$TMP_ROOT_DIR/mcr_cache_root";
`mkdir -p $MCR_CACHE_ROOT`;

# read in expdirs
my $expdirfile = $ARGV[0];
open(FILE,$expdirfile) or die("Could not open file $expdirfile for reading");

my $njobs = 0;
my $nskipped = 0;
my $njobssubmitted = 0;

# loop through each experiment directory
my @childs = ();
while(my $expdir = <FILE>){
    chomp $expdir;

    if($nskipped < $STARTJOB){
	print "Skipping job $nskipped: $expdir\n";
	$nskipped++;
	next;
    }
    if($njobssubmitted+$STARTJOB > $ENDJOB){
	last;
    }
    
    if(!$expdir){
	next;
    }
    if($expdir =~ /^\s*#/){
	next;
    }

    $njobssubmitted++;

    if(! -d $expdir){
	print "Directory $expdir does not exist\n";
	next;
    }

    #if($isfirst){
    #$isfirst = 0;
    #next;
    #}

    $expdir =~ /^(.*)\/([^\/]+)$/;
    my $rootdir = $1;
    my $basename = $2;

    print "*** $basename\n";
    
    # make a name for this job
    my $sgeid = "$KEYWORD" . "_$basename";
    $sgeid =~ s/\//_/g;
    $sgeid =~ s/\./_/g;
    $sgeid =~ s/\;/_/g;

    # names for temporary script and log file
    my $shfilename = "$expdir/$SCRIPTNAME"."_$sgeid" . ".sh";
    my $outfilename = "$expdir/$SCRIPTNAME"."_$sgeid" . ".log";

    # create temporary script to be submitted
    write_qsub_sh($shfilename,$expdir,$sgeid);

    # submit command
    my $cmd = "$shfilename 2>&1 > $outfilename";
    #my $cmd = qq~qsub -N $sgeid -j y -b y -o $outfilename -cwd -l limit50=1 $shfilename~;
    my $pid = fork();
    if($pid){
	$njobs++;
	push(@childs, $pid);
    }
    else{
	print "job $njobs: $cmd\n";
	#sleep(10);
	`$cmd`;
	$njobs--;
	exit(0);
    }

    if($njobs >= $MAXNJOBS){
	print "njobs = $njobs\n";
	my $child = wait();
	print "child $child finished\n";
	$njobs--;
    }

}

close(FILE);

foreach (@childs) {
    waitpid($_, 0);
    print "child $_ finished\n";
}

my $tmp = $njobssubmitted - $STARTJOB;
print "$tmp jobs queued.";

sub write_qsub_sh {
	my ($shfilename,$expdir,$jobid) = @_;
	
	open(SHFILE,">$shfilename") || die "Cannot write $shfilename";

	print SHFILE qq~#!/bin/bash
# $SCRIPTNAME test script.
# this script will be qsubed

export MCR_CACHE_ROOT=$MCR_CACHE_ROOT.$jobid

echo \$MCR_CACHE_ROOT
$SCRIPT $MCR $expdir $PARAMS

~;
	
#	print SHFILE qq~#delete itself
#rm -f \$0
#~;	
	
	close(SHFILE);
	
	chmod(0755, $shfilename);
}
