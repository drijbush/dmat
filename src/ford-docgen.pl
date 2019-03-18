#!/usr/bin/env perl
## This script is to filter Fortran files and add Ford compliant comments at key places

## Key assumptions:
# Comments on variables are *inline* and not on the next line
# Does *not* filter any existing !! out of main body
# Things cannot be predocumented.
# We are not using Ford's metadata system
# Fixes end block statements to have a space in them!
# Does not cover headers in interface blocks

## Written by J. Wilkins 2018

use v5.10;

# For bad code -- Good docs set to 1
$uglycode_nicedocs = 0;

$blank = qr/^\s*$/;
$comment = qr/^\s*!/;
$continuation = qr/&\s*$/;
$type = qr/^\s*type\s*(?!\()/;
$block = qr/(?:subroutine|function|do|module|program|if|interface|type|class)/;

while (<>) {
    # Forces one space between end and block type. If anyone doesn't want this, that's their fault, not mine.
    # Don't affect pre-processor commands
    s/(?<!#)end\s*($block)/end $1/;
    # Skip interfaces
    if (/^\s*(?<!end)(abstract\s)?\Winterface/) {
        while (!/end interface/) {print; $_ = <>;}
        print;
    }
    elsif (/^[^!]+(?<!end)\W(function|subroutine|module)\s+[a-zA-Z_]+/){
	$blocktype = $1;
        print;
        while (/$continuation/) {$_ = <>;print;}
        $_ = <>;
        #Loop through 'use's, 'implicit's and blank lines
        while(/^\s*(use|implicit.*)?$/){
            print;
            $_ = <>;
        }
        # Catch unheadered routine
        if(!/$comment/){say '!! This $blocktype has no header !'}
        # Loop through header
        while(/$comment/){
            s/!+/!!/;
            if ($uglycode_nicedocs) {s/!\s*$/<br>\n/;}
            print;
            $_ = <>;
        }
        print;
    }
    # Catch inline comments on declarations
    elsif (/^[^!]+::/){
        # Catch continuations
        while (/$continuation/) {$_ .= <>;}
        s/!+/!!/g;
        if ($uglycode_nicedocs) {s/!\s*$/<br>\n/;}
        print;
    }
    # Handle types
    elsif (/$type/) {
        print;
        $_ = <>;
        if (!/$comment/) {say '!! This type has no header !';}
        while (/$comment/) {
            s/!+/!!/;
            print;
            $_ = <>;
        }
    }
    else {print;}
}
