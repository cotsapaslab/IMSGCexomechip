#!/usr/bin/perl

# cleanplink.pl
# (c) Stephen D. Turner 2010 http://www.stephenturner.us/
# This is free open-source software.
# See http://gettinggeneticsdone.blogspot.com/p/copyright.html

my $help = "\nUsage: $0 <input whitespace file> <tab or comma>\n\n";
die $help if @ARGV<2;

$delimiter=pop(@ARGV);
die $help unless ($delimiter=~/tab/i|$delimiter=~/comma/i);
@inputfiles=@ARGV;

if ($delimiter =~ /comma/i) {
    foreach (@inputfiles) {

        open (IN,"<$_");
        open (OUT,">$_.csv");
        while (<IN>) {
            chomp;
            $_ =~ s/^\s+//;  #Trim whitespace at beginning
            $_ =~ s/\s+$//;  #Trim whitespace at end
            $_ =~ s/\s+/,/g; #Remaining whitespace into commas
            #$_ =~ s/NA/-9/g;#If you want to recode NA as -9
            print OUT "$_\n";
        }
    }
} elsif ($delimiter =~ /tab/i) {
    foreach (@inputfiles) {
        open (IN,"<$_");
        open (OUT,">$_.tab");
        while (<IN>) {
            chomp;
            $_ =~ s/^\s+//;  #Trim whitespace at beginning
            $_ =~ s/\s+$//;  #Trim whitespace at end
            $_ =~ s/\s+/\t/g;#Remaining whitespace into commas
            #$_ =~ s/NA/-9/g;#If you want to recode NA as -9
            print OUT "$_\n";
        }
    }
} else {
    die $help;
}