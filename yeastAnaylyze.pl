#Ming Creekmore
#A program to determine the TiTv ratios between 4 different yeast alignments

#open files
use Cwd;
my $dir = getcwd;
$yeastalign = $dir."/YeastAlignments";
opendir(DIR, $yeastalign) or die "can't open directory $yeastalign:$!"; print"\n";
mkdir "YeastAlignTiTvRatios";

#4 different yeast type names to look for in files
$yeastname1 = "Scer";
$yeastname2 = "Spar";
$yeastname3 = "Smik";
$yeastname4 = "Sbay";

#read each file in directory
while($filename = readdir DIR) {
    $ORFname = substr($filename, 1, 7);
    $WatsonCrick = substr ($filename, 7, 1);

    #only open if the filename is a yeast gene alignment
    $filelocation = $yeastalign."/"."$filename";
    if (length $ORFname == 7){
        open (INFILE, $filelocation) or die "Cannot open file";
    }
    else {
        next;
    }
    open (OUTFILE, ">"."./YeastAlignTiTvRatios/"."$ORFname".".txt") || die " could not open output file\n";

    #alignment sequences corresponding to 4 yeast
    $align1 = "";
    $align2 = "";
    $align3 = "";
    $align4 = "";

    #go through each line in file and read alignment,
    #concatenating to secified align sequence
    while(<INFILE>) {
        #removes return carriages
        chomp;
        #splits string based on spaces of what was read in
        @a = split(/\s+/, $_);

        #see if heading matches yeastname, if so then add to the
        #alignment string
        if($_ =~ m/^$yeastname1/){
            $align1 = $align1.$a[1];
        }
        elsif($_ =~ m/^$yeastname2/){
            $align2 = $align2.$a[1];
        }
        elsif($_ =~ m/^$yeastname3/){
            $align3 = $align3.$a[1];
        }
        elsif($_ =~ m/^$yeastname4/){
            $align4 = $align4.$a[1];
        }
    }

    #convert alignments to uppercase
    $align1 = uc($align1);
    $align2 = uc($align2);
    $align3 = uc($align3);
    $align4 = uc($align4);

    #if C then reverse alignment and get compliment
    if ($WatsonCrick eq "C"){
        my $revcom = reverse $align1; 
        $revcom =~ tr/ACGT/TGCA/; 
        $align1 = $revcom;

        $revcom = reverse $align2;
        $revcom =~ tr/ACGT/TGCA/; 
        $align2 = $revcom;

        $revcom = reverse $align3;
        $revcom =~ tr/ACGT/TGCA/; 
        $align3 = $revcom;

        $revcom = reverse $align4;
        $revcom =~ tr/ACGT/TGCA/; 
        $align4 = $revcom;
    }

    #list to keep track of number of translations or transversions
    #1st num: Compare align1 & 2
    #2nd num: Compare align1 & 3
    #3rd num: Compare align1 & 4
    #4th num: Compare align2 & 3
    #5th num: Compare align2 & 4
    #6th num: Compare align3 & 4
    @transitions = (0, 0, 0, 0, 0, 0);
    @transversions = (0, 0, 0, 0, 0, 0);

    #Add to transition or transversion for each pair in sequence
    for(my $i=0; $i<length($align1); $i++) {

        #for 1 & 2
        $state = gettrans(substr($align1, $i, 1), substr($align2, $i, 1));
        if($state eq "S") {
            @transitions[0]++;
        }
        elsif($state eq "T") {
            @transversions[0]++;
        }

        #for 1 & 3
        $state = gettrans(substr($align1, $i, 1), substr($align3, $i, 1));
        if($state eq "S") {
            @transitions[1]++;
        }
        elsif($state eq "T") {
            @transversions[1]++;
        }

        #for 1 & 4
        $state = gettrans(substr($align1, $i, 1), substr($align4, $i, 1));
        if($state eq "S") {
            @transitions[2]++;
        }
        elsif($state eq "T") {
            @transversions[2]++;
        }


        #for 2 & 3
        $state = gettrans(substr($align2, $i, 1), substr($align3, $i, 1));
        if($state eq "S") {
            @transitions[3]++;
        }
        elsif($state eq "T") {
            @transversions[3]++;
        }

        #for 2 & 4
        $state = gettrans(substr($align2, $i, 1), substr($align4, $i, 1));
        if($state eq "S") {
            @transitions[4]++;
        }
        elsif($state eq "T") {
            @transversions[4]++;
        }

        #for 3 & 4
        $state = gettrans(substr($align3, $i, 1), substr($align4, $i, 1));
        if($state eq "S") {
            @transitions[5]++;
        }
        elsif($state eq "T") {
            @transversions[5]++;
        }
    }

    #print TiTv ratios to outfile
    if(@transversions[0]!=0) {
        $ratio = @transitions[0]/@transversions[0];
    }
    else {
        $ratio = 1;
    }
    print OUTFILE "$yeastname1"." and "."$yeastname2\t"."TiTv: "."@transitions[0]"."/"."@transversions[0]"."\t"."Calculated Ratio: "."$ratio\n";
    if(@transversions[1]!=0) {
        $ratio = @transitions[1]/@transversions[1];
    }
    else {
        $ratio = 1;
    }
    print OUTFILE "$yeastname1"." and "."$yeastname3\t"."TiTv: "."@transitions[1]"."/"."@transversions[1]"."\t"."Calculated Ratio: "."$ratio\n";
    if(@transversions[2]!=0) {
        $ratio = @transitions[2]/@transversions[2];
    }
    else {
        $ratio = 1;
    }
    print OUTFILE "$yeastname1"." and "."$yeastname2\t"."TiTv: "."@transitions[2]"."/"."@transversions[2]"."\t"."Calculated Ratio: "."$ratio\n";
    if(@transversions[3]!=0) {
        $ratio = @transitions[3]/@transversions[3];
    }
    else {
        $ratio = 1;
    }
    print OUTFILE "$yeastname1"." and "."$yeastname2\t"."TiTv: "."@transitions[3]"."/"."@transversions[3]"."\t"."Calculated Ratio: "."$ratio\n";
    if(@transversions[4]!=0) {
        $ratio = @transitions[4]/@transversions[4];
    }
    else {
        $ratio = 1;
    }
    print OUTFILE "$yeastname1"." and "."$yeastname2\t"."TiTv: "."@transitions[4]"."/"."@transversions[4]"."\t"."Calculated Ratio: "."$ratio\n";
    if(@transversions[5]!=0) {
        $ratio = @transitions[5]/@transversions[5];
    }
    else {
        $ratio = 1;
    }
    print OUTFILE "$yeastname1"." and "."$yeastname2\t"."TiTv: "."@transitions[5]"."/"."@transversions[5]"."\t"."Calculated Ratio: "."$ratio\n";
    close OUTFILE;
}

print "end program\n";
exit;


#Given two nucleotide arguments, determines if there is a translation or a transversion
#Returns a character
#S stands for translation
#V stands for transversion
#otherwise return X
sub gettrans {
    $nuc1 = @_[0];
    $nuc2 = @_[1];
    if($nuc1 eq $nuc2) {
        return "X";
    }
    elsif($nuc1 eq '-' or $nuc2 eq '-') {
        return "X"
    }
    elsif(($nuc1 eq 'A' and $nuc2 eq 'G') or ($nuc1 eq 'G' and $nuc2 eq 'A')
        or ($nuc1 eq 'C' and $nuc2 eq 'T') or ($nuc1 eq 'T' and $nuc2 eq 'C')) {
        return "S";
    }
    else {
        return "T";
    }
}