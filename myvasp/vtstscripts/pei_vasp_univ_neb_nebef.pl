#!/usr/bin/env perl
#;-*- Perl -*-

# Script prints the force, energy etc of OUTCAR's in immediate subdir
# of present working directory.

# Taken from vtst, only closed the zip part, which is not needed in my case.
# Please move to vtst scripts, if you want to use it. change name to pei_vasp_univ_neb_nebef.pl
# he need to stay with other original vtst scripts, so that it can be used by vtst scripts.


use Cwd;

$dir = cwd;
$outlist =`ls -1 "$dir"/[0-9][0-9]/OUTCAR`;   #specifies location of OUTCARs
@outlist = split("\n",$outlist);

#printf "%4s %16s %16s %16s\n", "ID", "Max Force (eV/Å)", "Total Energy (eV)", "Relative E (eV)";

$i = 0;
foreach $outfile (@outlist) {
    # revised by J. Pei
    $energy = `grep 'energy  without' "$outfile" | tail -1 | awk '{print \$NF}'`;
    $force  = `grep 'max at' "$outfile" | tail -1 | awk '{print \$NF}'`;
    # revision closed
    if(!$i) { $e0 = $energy; }
    $rel = $energy - $e0;
    @f4 = ($i,$force,$energy,$rel);
    printf "%4i %16.8f %16.8f %16.8f \n",@f4;
    $i++;
}

