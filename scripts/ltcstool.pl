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
use constant R                => 8.31451;
use constant EPSILON          => 0; # -1e-4;
use constant CPLEX_INF        => Math::CPLEX::Base::CPX_INFBOUND();
use constant OPTI_MIP_STATE   => Math::CPLEX::Base::CPXMIP_OPTIMAL();
use constant OPTTOL_MIP_STATE => Math::CPLEX::Base::CPXMIP_OPTIMAL_TOL();

# parameters
##################################################
my $opt_string = 'hs:r:e:c:g:k:p:i:t:l:f:n:o:y:bd';
my %opt;
getopts( "$opt_string", \%opt ) or usage();
if ( $opt{h} or !$opt{s} or !$opt{r} or !$opt{e} or !$opt{c} or !$opt{g} ){
	usage();
}
my $sbml_file      = $opt{s};
my $r_file         = $opt{r};
my $modes_file     = $opt{e};
my $conc_file      = $opt{c};
my $gibbs_file     = $opt{g};
my $yield_file     = $opt{y} ? $opt{y} : undef;
my $temperature    = $opt{k} ? $opt{k} : 310.15;
my $ionic_strength = $opt{i} ? $opt{i} : 0.15;
my $ph             = $opt{p} ? $opt{p} : 7;
my $threads        = $opt{t} ? $opt{t} : 1;
my $lpfile         = $opt{l} ? $opt{l} : 'temp_';
my $solfile        = $opt{f} ? $opt{f} : 'sol_';
my $proton         = $opt{o} ? $opt{o} : '';
my $intervall      = -1;
if ($opt{n}) 
{
    $intervall      = $opt{n} > -1 ? $opt{n} : -1;
}
my $tightbounds    = $opt{b} ? 1 : 0;
my $debug          = $opt{d} ? 1 : 0;

# print log
##################################################
print "sbml file:             $sbml_file\n";
print "r-file:                $r_file\n";
print "EFM file:              $modes_file\n";
print "concentration file:    $conc_file\n";
print "energy file:           $gibbs_file\n";
print "yields file:           ";
if (defined($yield_file))
{
    print "$yield_file\n";
}
else
{
    print "not used\n";
}
print "temperature:           $temperature\n";
print "ionic strength:        $ionic_strength\n";
print "pH-value:              $ph\n";
print "proton:                $proton\n";
print "lp outputfile:         $lpfile\n";
print "lp solution file:      $solfile\n";
print "threads:               $threads\n";
print "intervall:             ";
if ($intervall > 0)
{
    print "$intervall\n";
}
else
{
    print "not used\n";
}
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
$cplex_env->writeparams('cplex_params.txt');

# parse input files
##################################################
my $rx_order      = parseRfile($r_file);
my ($spec, $reac) = parseSbmlFile($sbml_file, $rx_order);
my $gHash         = parseGibbsFile($gibbs_file);
my $cHash         = parseConcentrationFile($conc_file);
my $dfgs          = calcDfg($temperature, $ionic_strength, $ph, $gHash);
my $useRx         = getUseableReactions($reac, $gHash, $cHash, $proton);
my $rt            = R * $temperature / 1000;
my $yields;
if (defined($yield_file))
{
    my $yf  = readFile($yield_file);
    $yields = parseYieldsFile($yf, $rx_order);
}

# CREATE LP
###########
my $lp_info = createLpProblem($reac, $spec, $cHash, $rt, $dfgs, $modes_file,
		$useRx, $tightbounds, $proton, $yields); 

findAllSolutions($lp_info, $intervall, $lpfile, $solfile, $debug);

##################################################
# FUNCTIONS
##################################################

sub findAllSolutions {
	my ($lp_info, $intervall, $lpfile, $solfile, $debug) = @_;
	my $nr_efms       = $lp_info->{'info'}{'nr_efms'};
	my $lda_col_start = $lp_info->{'info'}{'lda_col_start'};
	my $obj_ix        = $lp_info->{'info'}{'objective_index'};
	my $lp            = $lp_info->{'lp'};

	my $sol_status;
	my @vals;
	my $mip_obj_val = $nr_efms;
	my $sol_count   = 0;
	my $first       = $lda_col_start;
	my $last        = $lda_col_start + $nr_efms;

	# write lp file without solutions
	#################################
	$lp->writeprob($lpfile.$sol_count.".lp", "LP");

    my $obj_col = $lp_info->{'info'}{'objective_col'};
    my $min = 1;
    if ($intervall > -1){
        $min = $nr_efms - $intervall;
        if ($min < 1){
            $min = 1;
        }
    }
    my $newRows;
    $newRows->[0][$obj_col] = 1;
    $newRows->[1][$obj_col] = 1;
    my $temp_rows = {
        num_rows => 2,
        rhs => [$min, $nr_efms],
        sense => ['G', 'L'],
        row_names => ['OBJ_MIN', 'OBJ_MAX'],
        row_coefs => $newRows
    };
    $lp_info->{'info'}{'obj_min'} = $min;
    $lp_info->{'info'}{'obj_max'} = "max";
    die "ERROR: addrows() failed for setting bounds of objective
    intervall\n" unless $lp->addrows($temp_rows);

	while (1){
		if ($debug){
            my $dtstring = localtime();
			print "$dtstring: objective minimum = $lp_info->{'info'}{'obj_min'}, maximum = $lp_info->{'info'}{'obj_max'}\n";
		}

		# solve mip
		############################
		$lp->mipopt();
		my $get_stat = $lp->getstat();

		if( defined $get_stat ){
			if ($get_stat == OPTI_MIP_STATE or $get_stat == OPTTOL_MIP_STATE or $get_stat == 118){
				$sol_count++;
				($sol_status, $mip_obj_val, @vals) = $lp->solution();
				my $solution_file = $solfile.$sol_count.'.csv';
				printSolution2file($solution_file, \@vals, $first, $last);
				addSolution2lp($lp_info, $sol_count, \@vals, $first, $last);
				changeBounds($lp_info, $mip_obj_val, $intervall);
				print "Cardinality: $mip_obj_val \t saved in $solution_file \n";
				$lp->writeprob($lpfile.$sol_count.".lp", "LP");
			} else {
				if ($intervall > -1 and $lp_info->{'info'}{'obj_min'} > 1){
					changeBounds($lp_info, $lp_info->{'info'}{'obj_min'}-1, $intervall);
				} else {
					last;
				}
			}
		} else {
			if ($intervall > -1 and $lp_info->{'info'}{'obj_min'} > 1){
				changeBounds($lp_info, $lp_info->{'info'}{'obj_min'}-1, $intervall);
			} else {
				last;
			}
		}
	}
}

sub changeBounds {
	my ($lp_info, $max, $intervall) = @_;
	my $lp = $lp_info->{'lp'};
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
    $lp_info->{'info'}{'obj_min'} = $min;
    $lp_info->{'info'}{'obj_max'} = $max;
    return $min;
}

sub addSolution2lp {
	my ($lp_info, $count, $vals, $start, $end) = @_;
	my $lp = $lp_info->{'lp'};
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
	my ($file, $vals, $start, $end) = @_;
	open(OUT, ">$file") or die ("Could not open $file: $!\n");
	for (my $i = $start; $i < $end; $i++) {
		print OUT ',' if ($i > $start);
		print OUT $vals->[$i];
	}
	print OUT "\n";
	close(OUT);
}

sub createLpProblem {
    my ($rx, $sp, $cH, $rt, $dfgs, $efms, $useRx, $tightbounds, $proton,
        $yields) = @_;

	my $lp = $cplex_env->createOP();
	my %lp_hash;
	$lp_hash{'col_count'} = -1;
	$lp_hash{'row_count'} = -1;

	# define dfG
	##############
	my $dfg_col_start = 0;
	$lp_hash{'dfg_col_start'} = $dfg_col_start;
	defineDfG($cH, $dfgs, $tightbounds, \%lp_hash, $proton); 

	# define drG
	##############
	my $drg_col_start = $lp_hash{'col_count'} + 1;
	$lp_hash{'drg_col_start'} = $drg_col_start;
	defineDrG($useRx, $rx, $cH, $tightbounds, \%lp_hash, $proton);

	# define rows for lambda sum
	#############################
	my $sum_rx_col_start = $lp_hash{'col_count'} + 1;
	my $sum_rx_row_start = $lp_hash{'row_count'} + 1;
	$lp_hash{'sum_rx_col_start'} = $sum_rx_col_start;
	$lp_hash{'sum_rx_row_start'} = $sum_rx_row_start;
	defineRowsForLambdaSum($useRx, $rx, \%lp_hash);
	$lp_hash{'sum_rx_col_end'} = $lp_hash{'col_count'};

	# define active EFMs (lambda)
	##############################
	my $lda_col_start = $lp_hash{'col_count'} + 1;
	$lp_hash{'lda_col_start'} = $lda_col_start;
	my $nr_efms = defineEfms($efms, $useRx, $rx, $yields, \%lp_hash);
	$lp_hash{'nr_efms'} = $nr_efms;

	# set upper bound
	##################
	setUpperBoundsForLambdaSums($useRx, \%lp_hash);

	# define columns for indicators
	#################################
	my $ind_col_start = $lp_hash{'col_count'} + 1;	
	$lp_hash{'ind_col_start'} = $ind_col_start;
	defineColsForIndicators($useRx, $rx, \%lp_hash);

	# define objective 
	##################
	my $objective_col = $lp_hash{'col_count'} + 1;
	$lp_hash{'objective_col'} = $objective_col;
	$lp_hash{'obj_min'} = 1;
	defineObjectiveInfo(\%lp_hash);
	$lp->maximize();

	# add cols and rows to lp
	###########################
	addRowsAndCols2lp($lp, \%lp_hash);
    addYieldConstraints($lp, $yields);

	# add indicators to lp
	#######################
	addIndicators2lp($lp, $useRx, $rx, \%lp_hash);
	addIndicatorsForReversibleReactions($lp, $useRx, $rx, \%lp_hash);

	my $ret;
	$ret->{'lp'} = $lp;
	$ret->{'info'} = \%lp_hash;
	return $ret; 
}

sub defineDfG {
	my ($cH, $dfgs, $tightbounds, $lp_hash, $proton) = @_;
	my $col_count = $lp_hash->{'col_count'};
	my $row_count = $lp_hash->{'row_count'};
    foreach my $x (keys(%{$cH})) {
        if ($x ne $proton) {
            my $metabolite = $cH->{$x};
            my $dfgName = $metabolite->{'dfgname'};
            if (defined($dfgs->{$dfgName})){
                my ($min, $max);
                if ($tightbounds){
                    ($min, $max) = getMinMaxDfg($metabolite, $dfgs, $dfgName);
                } else {
                    $min = -1 * CPLEX_INF;
                    $max = CPLEX_INF;
                }

                # dfg variables
                ###############
                $col_count++;
                $cH->{$x}{'index'}  = $col_count;
                $cH->{$x}{'dfgmin'} = $min;
                $cH->{$x}{'dfgmax'} = $max;
                $lp_hash->{'lwr_bound'}[$col_count] = $min;
                $lp_hash->{'upr_bound'}[$col_count] = $max;
                $lp_hash->{'col_types'}[$col_count] = 'C';
                $lp_hash->{'col_names'}[$col_count] = "dfg_".$x;

                # concentration variables
                #########################
                $col_count++;
                $lp_hash->{'lwr_bound'}[$col_count] = $metabolite->{'min'};
                $lp_hash->{'upr_bound'}[$col_count] = $metabolite->{'max'};
                $lp_hash->{'col_types'}[$col_count] = 'C';
                $lp_hash->{'col_names'}[$col_count] = "c_".$x;

                # dfg equation
                ##############
                $row_count++;
                $lp_hash->{'row_rhs'}[$row_count] = $dfgs->{$dfgName}{'dfg0'};
                $lp_hash->{'row_sense'}[$row_count]               = 'E';
                $lp_hash->{'row_names'}[$row_count]               = "DfG_$x";
                $lp_hash->{'row_coefs'}[$row_count][$col_count-1] = 1;
                $lp_hash->{'row_coefs'}[$row_count][$col_count]   = -$rt;
            }
        }
    }
	$lp_hash->{'col_count'} = $col_count;
	$lp_hash->{'row_count'} = $row_count;
}

sub getMinMaxDfg {
	my ($met, $dfgs, $dfgName) = @_;
	my $rhs = $dfgs->{$dfgName}{'dfg0'};
	my $min = $rhs + $rt * $met->{'min'};
	my $max = $rhs + $rt * $met->{'max'};
	if ($min > $max){
		my $t = $min;
		$min = $max;
		$max = $t;
	}
	$dfgs->{$dfgName}{'min'} = $min;
	$dfgs->{$dfgName}{'max'} = $max;
	return ($min, $max);
}

sub defineDrG {
	my ($useRx, $rx, $cH, $tightbounds, $lp_hash) = @_;
	my $col_count = $lp_hash->{'col_count'};
	my $row_count = $lp_hash->{'row_count'};
	my ($fct, $pre);
	my ($min, $max);
	foreach my $x (@{$useRx}) {
		if ($x->{'use'}){
			my $temp_rx = $rx->[$x->{'id'}];
			my $id  = $temp_rx->{'id'};
			my $rct = $temp_rx->{'reactant'};
			my $prd = $temp_rx->{'product'};
			for (my $i = 0; $i < 2; $i++) {
				if ($i == 0) {
					$fct = 1;
					$pre = '';
				} elsif ($i > 0 and $x->{'rev'}){
					$fct = -1;
					$pre = 'r';
				} else {
					$i = 2;
				}
				if ($tightbounds){
					$min = 0;
					$max = 0;
				} else {
					$min = -1 * CPLEX_INF;
					$max = CPLEX_INF;
				}
				if ($i < 2){
					$col_count++;
					$row_count++;
					$lp_hash->{'col_names'}[$col_count] = $pre.'drg_'.$id;
					$lp_hash->{'col_types'}[$col_count] = 'C';
					$lp_hash->{'row_rhs'}[$row_count]   = 0;
					$lp_hash->{'row_sense'}[$row_count] = 'E';
					$lp_hash->{'row_names'}[$row_count] = $pre.'DrG_'.$id;
					foreach my $temp (@{$rct}) {
                        if ($temp->{'species_id'} ne $proton)
                        {
                            my $tix = $cH->{$temp->{'species_id'}}{'index'};
                            my $stoich = $temp->{'stoichiometry'};
                            $lp_hash->{'row_coefs'}[$row_count][$tix] = $fct * $stoich;
                            if ($tightbounds){
                                my $val = ($i < 1) ?
                                -$cH->{$temp->{'species_id'}}{'dfgmax'}:
                                $cH->{$temp->{'species_id'}}{'dfgmin'};
                                $min += $stoich * $val;
                                $val = ($i < 1) ?
                                -$cH->{$temp->{'species_id'}}{'dfgmin'}:
                                $cH->{$temp->{'species_id'}}{'dfgmax'};
                                $max += $stoich * $val;
                            }
                        }
                    }
                    foreach my $temp (@{$prd}) {
                        if ($temp->{'species_id'} ne $proton)
                        {
                            my $tix = $cH->{$temp->{'species_id'}}{'index'};
                            my $stoich = $temp->{'stoichiometry'};
                            $lp_hash->{'row_coefs'}[$row_count][$tix] = -$fct * $stoich;
                            if ($tightbounds){
                                my $val = ($i < 1) ?
                                $cH->{$temp->{'species_id'}}{'dfgmin'}:
                                -$cH->{$temp->{'species_id'}}{'dfgmax'};
                                $min += $stoich * $val;
                                $val = ($i < 1) ?
                                $cH->{$temp->{'species_id'}}{'dfgmax'}:
                                -$cH->{$temp->{'species_id'}}{'dfgmin'};
                                $max += $stoich * $val;
                            }
                        }
                    }
					$lp_hash->{'lwr_bound'}[$col_count] = $min;
					$lp_hash->{'upr_bound'}[$col_count] = $max;
					$lp_hash->{'row_coefs'}[$row_count][$col_count] = 1;
				}
			}
		}
	}
	$lp_hash->{'col_count'} = $col_count;
	$lp_hash->{'row_count'} = $row_count;
}

sub defineRowsForLambdaSum {
	my ($useRx, $rx, $lp_hash) = @_;
	my $col_count = $lp_hash->{'col_count'};
	my $row_count = $lp_hash->{'row_count'};
	my $pre;
	foreach my $x (@{$useRx}) {
		my $id = $rx->[$x->{'id'}]{'id'};
		for (my $i = 0; $i < 2; $i++) {
			if ($i == 0) {
				$pre = '';
			} elsif ($i > 0 and $x->{'rev'}){
				$pre = 'r_';
			} else {
				$i = 2;
			}
			if ($i < 2){
				$row_count++;
				$col_count++;
				$lp_hash->{'lwr_bound'}[$col_count] = 0;
				$lp_hash->{'upr_bound'}[$col_count] = 0;
				$lp_hash->{'col_names'}[$col_count] = 's_'.$pre.'drg_'.$id;
				$lp_hash->{'col_types'}[$col_count] = 'C';
				$lp_hash->{'row_rhs'}[$row_count]               = 0;
				$lp_hash->{'row_sense'}[$row_count]             = 'E';
				$lp_hash->{'row_names'}[$row_count]             = 'SUM_'.$pre.'DrG_'.$id;
				$lp_hash->{'row_coefs'}[$row_count][$col_count] = -1;
			}
		}
	}
	$lp_hash->{'col_count'} = $col_count;
	$lp_hash->{'row_count'} = $row_count;
}

sub defineEfms {
	my ($efms, $useRx, $rx, $yields, $lp_hash) = @_;
	my $col_count = $lp_hash->{'col_count'};
	my $row_count = $lp_hash->{'row_count'};
	my $nr_efms = 0;
	open(IN, "<$efms") or die ("Could not open $efms: $!\n");
	while (my $line = <IN>){
        chomp($line);
        $line =~ s/^\s+//;
        $line =~ s/\s+$//;
        $line =~ s/Inf/1000/g;
        $line =~ s/NaN/0/g;
        $line =~ s/NA/0/g;
		if ($line =~ /\d/){
			$nr_efms++;
			$col_count++;
			$lp_hash->{'lwr_bound'}[$col_count] = 0;
			$lp_hash->{'upr_bound'}[$col_count] = 1;
			$lp_hash->{'col_names'}[$col_count] = "l$nr_efms";
			$lp_hash->{'col_types'}[$col_count] = 'B';
			my @spl = split(/\s+/, $line);
			my $act_rx_ix = $lp_hash->{'sum_rx_row_start'} - 1;
			foreach my $x (@{$useRx}){
				$act_rx_ix++;
				if ($spl[$x->{'id'}] > 0){
					$lp_hash->{'row_coefs'}[$act_rx_ix][$col_count] = 1;
				}
				if ($x->{'rev'}){
					$act_rx_ix++;
					if ($spl[$x->{'id'}] < 0){
						$lp_hash->{'row_coefs'}[$act_rx_ix][$col_count] = 1;
					}
				}
			}
            foreach my $x (keys(%{$yields})) {
                my $ix = $yields->{$x}{ix};
                my $val = $yields->{$x}{yield};
                my $comp = $yields->{$x}{comp};
                if ($comp eq "l"){
                    if ($spl[$ix] <= $val){
                        push(@{$yields->{$x}{efms}}, $col_count);
                    }
                } elsif ($comp eq "g"){
                    if ($spl[$ix] >= $val){
                        push(@{$yields->{$x}{efms}}, $col_count);
                    }
                } elsif ($comp eq "e"){
                    if ($spl[$ix] == $val){
                        push(@{$yields->{$x}{efms}}, $col_count);
                    }
                }
            }
		}
	}
	$lp_hash->{'col_count'} = $col_count;
	$lp_hash->{'row_count'} = $row_count;
	return $nr_efms;
}

sub setUpperBoundsForLambdaSums {
	my ($useRx, $lp_hash) = @_;
	my $nr_efms = $lp_hash->{'nr_efms'};
	my $index = $lp_hash->{'sum_rx_col_start'} - 1;
	foreach my $x (@{$useRx}) {
		$index++;
		$lp_hash->{'upr_bound'}[$index] = $nr_efms;
		if ($x->{'rev'}){
			$index++;
			$lp_hash->{'upr_bound'}[$index] = $nr_efms;
		}
	}
}

sub defineColsForIndicators {
	my ($useRx, $rx, $lp_hash) = @_;
	my $col_count = $lp_hash->{'col_count'};
	my $pre;
	foreach my $x (@{$useRx}) {
		my $id = $rx->[$x->{'id'}]{'id'};	
		for (my $i = 0; $i < 2; $i++) {
			if ($i == 0) {
				$pre = 'ind_';
			} elsif ($i > 0 and $x->{'rev'}){
				$pre = 'ind_r';
			} else {
				$i = 2;
			}
			if ($i < 2){
				$col_count++;
				$lp_hash->{'lwr_bound'}[$col_count] = 0;
				$lp_hash->{'upr_bound'}[$col_count] = 1;
				$lp_hash->{'col_names'}[$col_count] = $pre."drg_".$id;
				$lp_hash->{'col_types'}[$col_count] = 'B';
			}
		}
	}
	$lp_hash->{'col_count'} = $col_count;
}

sub defineObjectiveInfo {
	my ($lp_hash) = @_;
	my $nr_efms = $lp_hash->{'nr_efms'};
	my $col_count = $lp_hash->{'col_count'};
	my $row_count = $lp_hash->{'row_count'};
	my $lda_col_start = $lp_hash->{'lda_col_start'};
	$col_count++;
	$lp_hash->{'lwr_bound'}[$col_count] = 1;
	$lp_hash->{'upr_bound'}[$col_count] = $nr_efms;
	$lp_hash->{'col_names'}[$col_count] = "obj_col";
	$lp_hash->{'col_types'}[$col_count] = 'C';
	$lp_hash->{'obj_coefs'}[$col_count] = 1;

	$row_count++;
	$lp_hash->{'row_rhs'}[$row_count]   = 0;
	$lp_hash->{'row_sense'}[$row_count] = 'E';
	$lp_hash->{'row_names'}[$row_count] = 'obj_sum';
	for (my $i = 0; $i < $nr_efms; $i++) {
		my $ix = $lda_col_start + $i;
		$lp_hash->{'row_coefs'}[$row_count][$ix] = 1;
	}
	$lp_hash->{'row_coefs'}[$row_count][$col_count] = -1;

	$lp_hash->{'col_count'} = $col_count;
	$lp_hash->{'row_count'} = $row_count;
}

sub addRowsAndCols2lp {
	my ($lp, $lp_hash) = @_;
	my $cols = { num_cols => $lp_hash->{'col_count'} + 1,
		obj_coefs => $lp_hash->{'obj_coefs'},
		lower_bnd => $lp_hash->{'lwr_bound'},
		upper_bnd => $lp_hash->{'upr_bound'},
		col_names => $lp_hash->{'col_names'},
		col_types => $lp_hash->{'col_types'}
	};
	my $rows = { num_rows => $lp_hash->{'row_count'} + 1,
		rhs       => $lp_hash->{'row_rhs'},
		sense     => $lp_hash->{'row_sense'},
		row_names => $lp_hash->{'row_names'},
		row_coefs => $lp_hash->{'row_coefs'}
	};
	$lp->newcols($cols);
	$lp->addrows($rows);
}

sub addYieldConstraints {
    my ($lp, $yields) = @_;
    foreach my $x (keys(%{$yields})) {
        my $efms = $yields->{$x}{efms};
        if (defined($efms)) {
            my $row_coef;
            foreach my $y (@{$efms}) {
                $row_coef->[0][$y] = 1;
            }
            my $rhs->[0] = 1;
            my $sense->[0] = 'G';
            my $row_name->[0] = 'Yield_'.$x;
            my $new_row = {
                num_rows  => 1,
                rhs       => $rhs,
                sense     => $sense,
                row_names => $row_name,
                row_coefs => $row_coef};
            die "ERROR: addrows() failed\n" unless $lp->addrows($new_row);
        } else {
            print "not defined for $x\n";
        }
    }
}

sub	addIndicatorsForReversibleReactions {
	my ($lp, $useRx, $rx, $lp_hash) = @_;
	my $t_ind_ix = $lp_hash->{'ind_col_start'} - 1;
	my $new_row_coef;
	my $new_row_name;
	my $new_row_rhs;
	my $new_row_sense;
	my $num_rows = 0;
	foreach my $x (@{$useRx}) {
		if ($x->{'rev'}){
			$t_ind_ix++;
			my $id = $rx->[$x->{'id'}]{'id'};
			$new_row_coef->[$num_rows][$t_ind_ix] = 1;
			$t_ind_ix++;
			$new_row_coef->[$num_rows][$t_ind_ix] = 1;
			$new_row_rhs->[$num_rows] = 1;
			$new_row_sense->[$num_rows] = 'L';
			$new_row_name->[$num_rows] = 'DIR_'.$id;
			$num_rows++;
		} else {
			$t_ind_ix++;
		}
	}
	my $new_rows = {
		num_rows  => $num_rows,
		rhs       => $new_row_rhs,
		sense     => $new_row_sense,
		row_names => $new_row_name,
		row_coefs => $new_row_coef};
	die "ERROR: addrows() failed\n" unless $lp->addrows($new_rows);
}

sub	addIndicators2lp {
	my ($lp, $useRx, $rx, $lp_hash) = @_;
	my $t_sum_ix = $lp_hash->{'sum_rx_col_start'} - 1;
	my $t_drg_ix = $lp_hash->{'drg_col_start'} - 1;
	my $t_ind_ix = $lp_hash->{'ind_col_start'} - 1;
	my $pre;
	foreach my $x (@{$useRx}) {
		my $id = $rx->[$x->{'id'}]{'id'};
		for (my $i = 0; $i < 2; $i++) {
			if ($i == 0) {
				$pre = 'IND_';
			} elsif ($i > 0 and $x->{'rev'}){
				$pre = 'IND_r';
			} else {
				$i = 2;
			}
			if ($i < 2){
				$t_ind_ix++;
				$t_sum_ix++;
				my $sum_val;
				$sum_val->[$t_sum_ix] = 1;
				my $ind1 = {
					indvar       => $t_ind_ix,
					complemented => 1,
					rhs          => 0,
					sense        => 'E',
					val          => $sum_val,
					name         => $pre.'SUM_z_'.$id
				};
				die "ERROR: addindconstr() failed\n" unless
					$lp->addindconstr($ind1);
				my $ind1a = {
					indvar       => $t_ind_ix,
					complemented => 0,
					rhs          => 1,
					sense        => 'G',
					val          => $sum_val,
					name         => $pre.'SUM_'.$id
				};
				die "ERROR: addindconstr() failed\n" unless
					$lp->addindconstr($ind1a);
				if ($x->{'use'}){
					$t_drg_ix++;
					my $drg_val;
					$drg_val->[$t_drg_ix] = 1;
					my $ind2 = {
						indvar       => $t_ind_ix,
						complemented => 0,
						rhs          => EPSILON,
						sense        => 'L',
						val          => $drg_val,
						name         => $pre.'DrG_'.$id
					};
					die "ERROR: addindconstr() failed\n" unless
						$lp->addindconstr($ind2);
				}
			}
		}
	}
}

sub getUseableReactions {
	my ($rx, $gHash, $cHash, $proton) = @_;
	my @useRxs;
	for (my $i = 0; $i < @{$rx}; $i++) {
		my $act = $rx->[$i];
		my $rev = ($act->{'reversible'} eq 'true') ? 1 : 0;
		my %temp;
		my $use = 1;
		$use = isUsableReaction($act->{'product'}, $gHash, $cHash, $proton);
		if ($use){
			$use = isUsableReaction($act->{'reactant'}, $gHash, $cHash, $proton);
		}
		$temp{'id'} = $i;
		$temp{'rev'} = $rev;
		$temp{'use'} = $use ? 1 : 0;
		if ($use or $rev){
			push(@useRxs, \%temp);
		}
	}
	return (\@useRxs);
}

sub isUsableReaction {
    my ($array, $gHash, $cHash, $proton) = @_;
    my $use = 1;
    foreach my $x (@{$array}) {
        my $id = $x->{'species_id'};
        if ($id ne $proton)
        {
            if (!$cHash->{$id}){
                $use = 0;
                last;
            } elsif (!$gHash->{$cHash->{$id}{'dfgname'}}) {
                $use = 0;
                last;
            }
        }
    }
	return $use;
}

sub parseYieldsFile {
    my ($yf, $rf) = @_;
    my $yields;
    foreach my $x (@{$yf}) {
        if ($x =~ /[<=>]/)
        {
            $x =~ s/\s+/ /g;
            chomp($x);
            my @spl = split(/ /, $x);
            my $name = $spl[0];
            my $ix = -1;
            for (my $i = 0; $i < @{$rf}; $i++) {
                if ($rf->[$i] eq $name)
                {
                    $ix = $i;
                    last;
                }
            }
            if ($ix < 0) {
                die ("Could not find $name in list of reactions. Stopped at parsing yields file!\n");
            }
            $yields->{$name}{ix} = $ix;
            if ($spl[1] eq "<") {
                $yields->{$name}{comp} = "l";
            } elsif ($spl[1] eq ">") {
                $yields->{$name}{comp} = "g";
            } elsif ($spl[1] eq "=") {
                $yields->{$name}{comp} = "e";
            } else {
                die ("Comparison value of $name is not <, = or > ( $spl[1] ). Stopped at parsing yields file!\n");
            }
            if ($spl[2] !~ /^[+-]*\d+\.*\d*$/) {
                die ("Yield value of $name is not a number ( $spl[2] ). Stopped at parsing yields file!\n");
            }
            $yields->{$name}{yield} = $spl[2];
            $yields->{$name}{efms} = ();
        }
    }
    return $yields;
}

sub parseRfile {
	my $f = shift;
	my $rf = readFile($f);
	my @rx;
	foreach my $x (@{$rf}) {
		if ($x =~ /\D/){
			chomp($x);
			$x =~ s/"//g;
			$x =~ s/^\s+//;
			$x =~ s/\s+$//;
			@rx = split(/\s+/, $x);
			return \@rx;
		}
	}
	return undef;
}

sub parseSbmlFile {
	my ($sbml_file, $rx_order) = @_;
	my $sbml_ref = XMLin($sbml_file, KeyAttr => 'id' );
	my $spec     = getSpecies($sbml_ref);
	my $reac     = getReactions($sbml_ref, $spec, $rx_order);
	return ($spec, $reac);
}

sub getSpecies {
	my ($sr) = @_;
	my $spc;
	foreach my $id (keys %{$sr->{'model'}{'listOfSpecies'}{'species'}}){
		$spc->{$id}{'id'}          = $id;      
		$spc->{$id}{'compartment'} = $sr->{'model'}{'listOfSpecies'}{'species'}{$id}{'compartment'};
	}   
	return $spc;
}

sub getReactions {
	my ($sr, $sp, $rx_order) = @_;
	my $reacs;
	my $nr  = -1;
	my $rev = 0;
	my $reac_ref = $sr->{'model'}{'listOfReactions'}{'reaction'};
	foreach my $reaction (keys %{$reac_ref}) {
		$nr = -1;
		for (my $rxi = 0; $rxi < @{$rx_order}; $rxi++) {
			if ($reaction eq $rx_order->[$rxi]){
				$nr = $rxi;
				last;
			}
		}
		if ($nr < 0){
			die ("$reaction not found in getReactions\n");
		}
		$reacs->[$nr]{'id'}         = $reaction;
		$reacs->[$nr]{'reversible'} = $reac_ref->{$reaction}{'reversible'} || 'true';

		# list of reactant species
		##########################
		my $nr_rct = 0;
		my $inpSpec_ref = $reac_ref->{$reaction}{'listOfReactants'}{'speciesReference'};
		if ($inpSpec_ref) {
			my $var_type = ref $inpSpec_ref;
			if( $var_type eq 'ARRAY' ) {
				foreach my $species (@{$inpSpec_ref} ) {
					my $specID = $species->{'species'};
					$reacs->[$nr]{'reactant'}[$nr_rct]{'species_id'}    = $specID;
					$reacs->[$nr]{'reactant'}[$nr_rct]{'stoichiometry'} = $species->{'stoichiometry'} || 1.0;
					$nr_rct++;
				}
			} elsif( $var_type eq 'HASH' ) {
				my $specID = $inpSpec_ref->{'species'};
				$reacs->[$nr]{'reactant'}[$nr_rct]{'species_id'}    = $specID;
				$reacs->[$nr]{'reactant'}[$nr_rct]{'stoichiometry'} = $inpSpec_ref->{'stoichiometry'} || 1.0;
				$nr_rct++;
			} else {
				warn "ERROR: unexpected variable type '$var_type' for listOfReactants of reaction '$reaction'\n";
			}
		}

		# list of product species
		#########################
		my $nr_prd = 0;
		my $outSpec_ref = $reac_ref->{$reaction}{'listOfProducts'}{'speciesReference'};
		if ($outSpec_ref) {
			my $var_type = ref $outSpec_ref;
			if( $var_type eq 'ARRAY' ) {
				foreach my $species (@{$outSpec_ref} ) {
					my $specID = $species->{'species'};
					$reacs->[$nr]{'product'}[$nr_prd]{'species_id'}    = $species->{'species'};
					$reacs->[$nr]{'product'}[$nr_prd]{'stoichiometry'} = $species->{'stoichiometry'} || 1.0;
					$nr_prd++;
				}
			}
			elsif( $var_type eq 'HASH' ) {
				my $specID = $outSpec_ref->{'species'};
				$reacs->[$nr]{'product'}[$nr_prd]{'species_id'}    = $outSpec_ref->{'species'};
				$reacs->[$nr]{'product'}[$nr_prd]{'stoichiometry'} = $outSpec_ref->{'stoichiometry'} || 1.0;
				$nr_prd++;
			} else {
				warn "ERROR: unexpected variable type '$var_type' for listOfReactants of reaction '$reaction'\n";
			}
		}
	}
	return ($reacs);
}

sub parseGibbsFile {
	my $f = shift;
	my $file = readFile($f);
	my %hash;
	foreach my $line (@{$file}) {
		if ($line =~ /=/){
			$line =~ s/\s+//g;
			chomp($line);
			my @spl  = split(/=/, $line);
			my $name = $spl[0];
			my $temp = $spl[1];
			my @spl2 = split(/\)/, $temp);
			my @dfg;
			my @nh;
			my @charge;
			for (my $i = 0; $i < @spl2; $i++) {
				if ($spl2[$i] =~ /\(/){	
					my @spl3 = split(/\(/, $spl2[$i]);
					my @spl4 = split(/,/, $spl3[1]);
					if (@spl4 != 3){
						die("Execution aborted:\nGibbs file is not correct in
								following line:\n$line\n");
					}
					$dfg[$i]    = $spl4[0];
					$charge[$i] = $spl4[1];
					$nh[$i]     = $spl4[2];
				} else {
					die("Execution aborted:\nGibbs file is not correct in
					following line:\n$line\n");
				}
			}
			$hash{$name}{'dfg'}    = \@dfg;
			$hash{$name}{'charge'} = \@charge;
			$hash{$name}{'nh'}     = \@nh;
			$hash{$name}{'nsp'}    = scalar(@spl2);
		}
	}
	return \%hash;
}

sub parseConcentrationFile {
	my $file = shift;
	my $f = readFile($file);
	my %hash;
	foreach my $x (@{$f}) {
		if ($x =~ /;/){
			$x =~ s/\s+//g;
			chomp($x);
			my @spl = split(/;/, $x);
			$hash{$spl[0]}{'dfgname'} = $spl[1];
			$hash{$spl[0]}{'min'}     = log($spl[2]);
			$hash{$spl[0]}{'max'}     = log($spl[3]);
		}
	}
	return \%hash;
}

sub calcDfg {
	my ($tm, $is, $ph, $gh) = @_;
	my $rt      = R * $tm / 1000;
	my $pH_part = $rt * log( 10 ** -$ph );
	my $is_pow  = $is ** 0.5;
	my $is_part = 2.91482 * $is_pow / ( 1 + ( 1.6 * $is_pow ) );
	my %hash;
	foreach my $k (keys(%{$gh})) {
		my $sp = $gh->{$k};
		my @exp;
		for (my $i = 0; $i < $sp->{'nsp'}; $i++) {
			my $ph_term = $sp->{'nh'}[$i] * $pH_part;
			my $is_term = ( $sp->{'charge'}[$i] ** 2 - $sp->{'nh'}[$i] ) * $is_part;
			my $gpfnsp  = $sp->{'dfg'}[$i] - $ph_term - $is_term;
			my $temp    = -$gpfnsp / $rt;
			push(@exp, $temp);
		}
		my $dfg   = -$rt * logSumOfExponentials(\@exp);
		$hash{$k}{'dfg0'} = $dfg;
	}
	return \%hash;
}

sub logSumOfExponentials {
	my $exp = shift;
	if (@{$exp} == 1){
		return $exp->[0];
	}
	my $max = $exp->[0];
	for (my $i = 1; $i < @{$exp}; $i++) {
		if ($exp->[$i] > $max){
			$max = $exp->[$i];
		}
	}
	my $sum = 0;
	foreach my $x (@{$exp}) {
		$sum += exp( $x - $max );
	}
	my $logSum = $max + log($sum);
	return $logSum;
}

sub readFile {
    my ($file) = @_;
    open(IN, "<$file") or die ("Could not open $file for reading: $!\n");
    my @in = <IN>;
    close(IN);
    return \@in;
}

sub usage {
    print <<"END_TEXT";

finds the largest feasible set of EFMs

Usage: $0 [OPTIONS]


   -h    show this help
   -s    sbml file of the model
   -r    rfile of the model
   -e    file with EFMs in double-text format as printed by efmtool
   -c    file with concentration of metabolites
   -g    file with gibbs energies
   -o    proton name in model
   -k    temperature in K (default: 310.15 K)
   -p    pH value (default: 7)
   -i    ionic strength (default: 0.15)
   -l    output prefix of lp files (default: temp_)
   -f    output prefix of solution files (default: sol_)
   -t    number of threads (default: 1)
   -b    bounds tightening [optional]
         this option was implemented for testing, but it slows down the solver 
         dramatically and therefore it is recommended not be used
   -n    search within intervalls (number between 1 and number of EFMs) [optional]
         increases efficiency of the solver for bigger problems and therfore it
         is recommended to use an appropriate interval 
   -y    yield file [optional]
         format:
           R_BIOMASS > 0.009528785
           R_EX_co2_e > 0.2288889
           R_EX_o2_e < -0.2301605
           R_ATPS4r > 0.1653086
   -d    debug output

   Examples:
   $0 -s e_coli.xml -r e_coli.rfile -e efm.modes -c glc.conc -g group_contribution.txt -t 8 -o M_h_c
   $0 -s e_coli.xml -r e_coli.rfile -e efm.modes -c glc.conc -g group_contribution.txt -k 273 -i 0.2 -p 7.4 -n 8 -proton M_h_c
   $0 -s e_coli.xml -r e_coli.rfile -e efm.modes -c glc.conc -g group_contribution.txt -n 1000 -b -proton M_h_c

END_TEXT

   exit( 0 );
}

