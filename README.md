# ltcsCalculator

ltcsCalculator provides tools to find largest thermodynamic consistent sets of
EFMs.

If you use this software please cite 

> Matthias P. Gerstl, Christian Jungreuthmayer, Stefan Müller, and Jürgen
> Zanghellini. (2016), Which sets of elementary flux modes form
> thermodynamically feasible flux distributions?. FEBS J, 283: 1782–1794.
> doi:10.1111/febs.13702

## Table of Contents

[Installation](#Installation)

[Examples](#Examples)

[LTCS calculation](#ltcs calculation)

* ltcstool.pl
* ltcstool_clever.pl
* ltcstool_hotstart.pl
* ltcstool_hotstart_clever.pl
* calcLtcs

## <a name="Installation"></a>Installation

1. Download ltcsCalculator
2. The generic calcLtcs is written in C and can be compiled using gcc compiler
   by performing:
   ```
   cd ltcsCalculator
   make
   ```
   calcLtcs is then located in ltcsCalculator/bin
3. Perl scripts are located in folder scripts and can be executed without
   compilation. All perl scripts can be started with -h to see the help page.

## <a name="Examples"></a>Examples

Examples for all programs and scripts are located in the examples folder. To
remove all produced files run `clean_examples.sh`.

## <a name="ltcs calculation"></a>LTCS calculation

This section describes tools that calculates LTCS. All scripts result in the
same LTCS. However, I would recomment to use the _clever scripts, as it seems
that they have a better performance.

**ltcstool.pl**

This perl script calculates LTCS by efmtool input files. To make it
easier for the solver an intervall can be given, which is used to narrow the
space of allowed objective value. Additional yield information can be given.
But be careful to use an already normalized EFM set. This script DOES NOT check
if the EFMs are normalized! For detailed information run `ltcstool.pl -h`.

**ltcstool_clever.pl** 

This perl script calculates LTCS as well as ltcstool.pl, but it reduces the
search intervall itself after timeout of cplex. An alternative timeout for
cplex can be given as argument. Additional yield information can be given. See
ltcstool.pl for or run `ltcstool_clever.pl -h` for more information.

**ltcstool_hotstart.pl**

This perl script calculates LTCS when given a LP file created by ltcstool.pl or
ltcstool_clever.pl. The calculation method is the same as for ltcstool.pl. Run
`ltcstool_hotstart.pl -h` for detailed information.

**ltcstool_hotstart_clever.pl**

This perl script further calculates LTCS when given a LP file created by
ltcstool.pl or ltcstool_clever.pl. The calculation method is the same as for
ltcstool_clever.pl. Run `ltcstool_hotstart_clever.pl -h` for detailed
information.

**calcLtcs**

This C program calculates LTCS without considering concentrations. This tool
only uses the reversibility information of the reactions in EFMs. For detailed
information run `calcLtcs` without any argument
