#!/usr/bin/env perl 
#/////////////////////////////////////////////////////////////////////////////
# Author: Matthias Gerstl 
# Email: matthias.gerstl@acib.at 
# Company: Austrian Centre of Industrial Biotechnology (ACIB) 
# Web: http://www.acib.at Copyright
# (C) 2015 Published unter GNU Public License V3
#/////////////////////////////////////////////////////////////////////////////
#Basic Permissions.
# 
# All rights granted under this License are granted for the term of copyright
# on the Program, and are irrevocable provided the stated conditions are met.
# This License explicitly affirms your unlimited permission to run the
# unmodified Program. The output from running a covered work is covered by
# this License only if the output, given its content, constitutes a covered
# work. This License acknowledges your rights of fair use or other equivalent,
# as provided by copyright law.
# 
# You may make, run and propagate covered works that you do not convey,
# without conditions so long as your license otherwise remains in force. You
# may convey covered works to others for the sole purpose of having them make
# modifications exclusively for you, or provide you with facilities for
# running those works, provided that you comply with the terms of this License
# in conveying all material for which you do not control copyright. Those thus
# making or running the covered works for you must do so exclusively on your
# behalf, under your direction and control, on terms that prohibit them from
# making any copies of your copyrighted material outside their relationship
# with you.
# 
# Disclaimer of Warranty.
# 
# THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE
# LAW. EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR
# OTHER PARTIES PROVIDE THE PROGRAM “AS IS” WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE
# ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU.
# SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY
# SERVICING, REPAIR OR CORRECTION.
# 
# Limitation of Liability.
# 
# IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING WILL
# ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MODIFIES AND/OR CONVEYS THE
# PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES, INCLUDING ANY
# GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE
# OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF DATA
# OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD
# PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS),
# EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
# SUCH DAMAGES.
#/////////////////////////////////////////////////////////////////////////////

use strict;
use warnings;
use utf8;
use XML::Simple;
use Math::CPLEX::OP;
use Math::CPLEX::Env;
use Getopt::Std;
use Data::Dumper;

# constants
##################################################
use constant OPTI_MIP_STATE   => Math::CPLEX::Base::CPXMIP_OPTIMAL();
use constant OPTTOL_MIP_STATE => Math::CPLEX::Base::CPXMIP_OPTIMAL_TOL();
use constant TIME_LIM         => Math::CPLEX::Base::CPX_STAT_ABORT_TIME_LIM();
use constant TIME_LIM_FEAS    => Math::CPLEX::Base::CPXMIP_TIME_LIM_FEAS();
use constant TIME_LIM_INFEAS  => Math::CPLEX::Base::CPXMIP_TIME_LIM_INFEAS();
use constant TIME_LIMIT_UNLIM => 31536000; # 1 year

# parameters
##################################################
my $opt_string = 'hi:s:t:l:f:n:x:d';
my %opt;
getopts( "$opt_string", \%opt ) or usage();
if ( $opt{h} or !$opt{i} ){
	usage();
}
my $input_lp       = $opt{i};
my $threads        = $opt{t} ? $opt{t} : 1;
my $lpfile         = $opt{l} ? $opt{l} : 'temp_';
my $solfile        = $opt{f} ? $opt{f} : 'sol_';
my $tightbounds    = $opt{b} ? 1 : 0;
my $debug          = $opt{d} ? 1 : 0;
my $start_val      = $opt{s} ? $opt{s} : 0;
my $time_limit     = $opt{x} ? $opt{x} : 600; # ten minutes
my $max_intervall  = $opt{n} ? $opt{n} : 0;

if (!$start_val) {
    $start_val = `grep OBJ_MAX $input_lp | awk -F '<= ' '{print \$2}'`;
    chomp($start_val);
}
my $efm_nr         = getEfmNr($input_lp);

# print log
##################################################
print "input lp file:         $input_lp\n";
print "lp outputfile:         $lpfile\n";
print "lp solution file:      $solfile\n";
print "threads:               $threads\n";
print "intervall:             clever\n";
print "max intervall:         ";
if ($max_intervall)
{
    print $max_intervall."\n";
}
else
{
    print "not defined\n";
}
print "time limit:            $time_limit\n";
print "tightened bounds:      ";
if ($tightbounds)
{
    print "yes\n";
}
else
{
    print "no\n";
}

# set cplex parameters
##################################################
my $cplex_env      = Math::CPLEX::Env->new();
$cplex_env->setintparam(&Math::CPLEX::Base::CPX_PARAM_THREADS, $threads);
$cplex_env->setdblparam(&Math::CPLEX::Base::CPX_PARAM_TILIM, $time_limit);
$cplex_env->writeparams('cplex_params.txt');

# Load LP
#########
my $lp = $cplex_env->createOP();
$lp->readcopyprob($input_lp);

my $sol_count = `grep -c SOL_ $input_lp`;

findAllSolutions($lp, $start_val, $lpfile, $solfile, $sol_count,
    $efm_nr, $max_intervall, $debug);

##################################################
# FUNCTIONS
##################################################

sub getEfmNr {
    my ($f) = @_;
    open(IN, "<$f") or die ("Could not open $f: $!\n");
    my @in = <IN>;
    close(IN);
    my $inBound = 0;
    my $efm_nr = 0;
    for (my $i = 0; $i < @in; $i++) {
        if ($inBound){
            if ($in[$i] =~ /^Binaries/){
                $i = @in;
            } elsif ($in[$i] =~ /0 <= l\d+ <= 1/) {
                my @spl = split(/l/, $in[$i]);
                my @spl2 = split(/ <=/, $spl[1]);
                my $t = $spl2[0];
                chomp($t);
                if ($t > $efm_nr){
                    $efm_nr = $t;
                }
            }
        } else {
            if ($in[$i] =~ /^Bounds/){
                $inBound = 1;
            }
        }
    }
    return $efm_nr;
}

sub findAllSolutions {
    my ($lp, $start_val, $lpfile, $solfile, $sol_count, $efm_nr,
        $max_intervall, $debug) = @_;

	my $sol_status;
    my $mip_obj_val;
	my @vals;

    my $obj_min_ix = $lp->getrowindex("OBJ_MIN");

    my @match;

    my $first = 0;
    my $last = 0;
    for (my $i=1; $i <= $efm_nr; $i++) {
        my $temp = $lp->getcolindex("l".$i);
        $match[$i-1] = $temp;
        if ($temp < $first or $first < 1){
            $first = $temp;
        }
        if ($temp > $last) {
            $last = $temp;
        }
    }

    my $max = $start_val;
    my $intervall = getIntervall($max, $max_intervall);
    my $min = changeBounds($lp, $max, $intervall);

	while (1){
        printDebug($debug, $intervall, $min, $max);

		# solve mip
		############################
		$lp->mipopt();
		my $get_stat = $lp->getstat();

        my $printProblem = 0;
        if( defined $get_stat )
        {
            if ($get_stat == TIME_LIM or $get_stat == TIME_LIM_FEAS or $get_stat == TIME_LIM_INFEAS)
            {
                if ($intervall > 0) {
                    my $temp = length($intervall - 1) - 1;
                    $intervall = 10**$temp;
                    if ($intervall == 1) {
                        $intervall = 0;
                    }
                } else {
                    die ("Could not find solution within time limit\n");
                }
            }
            elsif ($get_stat == OPTI_MIP_STATE or $get_stat == OPTTOL_MIP_STATE or $get_stat == 118)
            {
                $sol_count++;
                ($sol_status, $mip_obj_val, @vals) = $lp->solution();
                my $solution_file = $solfile.$sol_count.'.csv';
                printSolution2file($solution_file, \@vals, \@match);
                addSolution2lp($lp, $sol_count, \@vals, $first, $last);
                $mip_obj_val = int($mip_obj_val + 0.5);
                print "Cardinality: $mip_obj_val \t saved in $solution_file \n";
                $max = $mip_obj_val;
                $intervall = getIntervall($max, $max_intervall);
                $printProblem = 1;
            } 
            else 
            {
                if ($min > 1)
                {
                    $max = $min - 1;
                    $intervall = getIntervall($max, $max_intervall);
                } 
                else 
                {
                    last;
                }
            }
        } 
        else 
        {
            if ($min > 1)
            {
                $max = $min - 1;
                $intervall = getIntervall($max, $max_intervall);
            } 
            else 
            {
                last;
            }
        }
        $min = changeBounds($lp, $max, $intervall);
        setTimeLimit($intervall);
        if ($printProblem)
        {
            $lp->writeprob($lpfile.$sol_count.".lp", "LP");
        }
    }
}

sub getIntervall
{
    my ($max, $max_intervall) = @_;
    my $intervall = $max - 1;
    if ($max_intervall)
    {
        if ($max_intervall < $intervall)
        {
            $intervall = $max_intervall;
        }
    }
    return $intervall;
}

sub changeBounds {
	my ($lp, $max, $intervall) = @_;
    my $min = 1;
    if ($intervall > -1){
		my $new_lwr = $max - $intervall;
		if ($new_lwr < 1){
            $new_lwr = 1;
		}
        $min = $new_lwr;
	}
    my $obj_min_ix = $lp->getrowindex("OBJ_MIN");
    my $obj_max_ix = $lp->getrowindex("OBJ_MAX");
    my $newrhs->[$obj_min_ix] = $min;
    $newrhs->[$obj_max_ix] = $max;
    die "ERROR: change new bounds of objective failed\n" unless $lp->chgrhs($newrhs);
    return $min;
}

sub addSolution2lp {
	my ($lp, $count, $vals, $start, $end) = @_;
	my $new_row_coef;
	my $new_row_name  = ['SOL_'.$count];
	my $new_row_rhs   = [1];
	my $new_row_sense = ['G'];
	for (my $i = $start; $i < $end; $i++) {
		$new_row_coef->[0][$i] = $vals->[$i] ? 0 : 1;
	}
	my $new_row = {
		num_rows  => 1,
		rhs       => $new_row_rhs,
		sense     => $new_row_sense,
		row_names => $new_row_name,
		row_coefs => $new_row_coef};
	die "ERROR: addrows() failed\n" unless $lp->addrows($new_row);
}

sub printSolution2file {
	my ($file, $vals, $match) = @_;
	open(OUT, ">$file") or die ("Could not open $file: $!\n");
	for (my $i = 0; $i < @{$match}; $i++) {
		print OUT ',' if ($i > 0);
		print OUT $vals->[$match->[$i]];
	}
	print OUT "\n";
	close(OUT);
}

sub printDebug {
    my ($debug, $intervall, $min, $max) = @_;
    if ($debug){
        my $dtstring = localtime();
        print "$dtstring: Intervall = $intervall,\tobjective minimum = $min,\tmaximum = $max\n";
    }
}

sub setTimeLimit {
    my $intervall = shift;
    if ($intervall > 0){
        $cplex_env->setdblparam(&Math::CPLEX::Base::CPX_PARAM_TILIM, $time_limit);
    } else {
        $cplex_env->setdblparam(&Math::CPLEX::Base::CPX_PARAM_TILIM, TIME_LIMIT_UNLIM);
    }
}


sub usage {
    print <<"END_TEXT";

finds the largest feasible set of EFMs when given a LP file
produced by ltcstool.pl

Usage: $0 [OPTIONS]


   -h    show this help
   -i    input LP file, created by ltcstool.pl
   -s    start maximum value [optional]
   -l    output prefix of lp files (default: temp_)
   -f    output prefix of solution files (default: sol_)
   -t    number of threads (default: 1)
   -n    search within intervalls (number between 1 and number of EFMs) [optional]
         increases efficiency of the solver for bigger problems and therfore it
         is recommended to use an appropriate interval 
   -x    time out intervall in seconds [default: 600]
   -d    debug output

   Examples:
   $0 -i step3.lp -s 32000 -l lp_ -f sol_ -t 8 -n 1000 -d
   $0 -i step3.lp -s 32000 -l lp_ -f sol_ -t 8 -n 1000 -x 60

END_TEXT

   exit( 0 );
}

