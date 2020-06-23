
#!/usr/bin/perl
use strict;
use warnings;

#the program uses need to be put in folder having files K12genesCDS.txt and K12.fasta to compute the basic numbers, as well as with a file ObservedMutations.txt (see teh file for how to fill it)


my %CodeGenetic=('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'_','TAG'=>'_','TGC'=>'C','TGT'=>'C','TGA'=>'_','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G');


my  (%AAlist)=("F"=>0,"L"=>0,"S"=>0,"Y"=>0, "STOP"=>0,"C"=>0,"W"=>0,"P"=>0,"H"=>0,"Q"=>0,"R"=>0,"I"=>0,"M"=>0,"T"=>0,"N"=>0,"K"=>0,"V"=>0,"A"=>0,"D"=>0,"E"=>0,"G"=>0);

my %MutTypeAction=("1_A"=>"C","1_T"=>"G","2_A"=>"G","2_T"=>"C","3_A"=>"T","3_T"=>"A","4_C"=>"G","4_G"=>"C","5_C"=>"T","5_G"=>"A","6_G"=>"T","6_C"=>"A");
my %Codonlist;
my %CodontoAA_mutType;

for my $key (keys %CodeGenetic)
{
    $Codonlist{$key}=0;
    
}

my @GeneName=(1..6000);

my %ContigShift=();
my $TransitShift=0;
my @b;
open (GENELIST, "<patient_17.gff");
open (OUTPUTCDS, ">OutputCDS.txt");
my $refsequence="";


my $line = <GENELIST>;
my $currentname="";
my $nblines=0;
while( $line = <GENELIST>)
{
    chomp $line;
    if ($line=~/##/)
    {
        if ($line!~/FASTA/)
        {
        @b=split(/ /,$line);
        $ContigShift{$b[1]}=$TransitShift;
        $TransitShift+=$b[3];
        }
    }
    else
    {
     if ($line=~/^>/)
     {
         $currentname=$line;
     }
        else
        {
        if ($line=~/NODE_/)
        {
         @b=split(/\t/,$line);
        if ($b[2]eq "CDS")
        {
            my $deb=$b[3]+$ContigShift{$b[0]};
            my $end=$b[4]+$ContigShift{$b[0]};
            my $strand="D";
            if ($b[6]eq"-")
            {
               $strand="C";
            }
            $b[8]=~/ID=(PROKKA_\d+);/;
            # print "$b[8]\n";
            my $name=$1;
            
            printf OUTPUTCDS "$deb\t$end\t$strand\t$name\tCDS\n";
        }
            
            
        }
            else
            {
                $refsequence.="$line";
                $nblines++;
                my $ll=length($refsequence);#-$nblines;
                print "$currentname\t$nblines\t$ll\n";
            }
            
        }
    
        
        
    }
    
}
close GENELIST;
close OUTPUTCDS;

open (OUTPUTFASTA, ">concat.fa");
printf OUTPUTFASTA ">ConcatenatedGenome\n$refsequence";
close OUTPUTCDS;

open (GENELIST, "<concat.fa");
my $refsequence2="_";
$line = <GENELIST>;
while( $line = <GENELIST>)
{
    chomp $line;
    $refsequence2.=$line;
    
}
close GENELIST;

my $ll=length($refsequence);
print "$ll\n";

open (GENELIST, "<OutputCDS.txt");

my $nb=0;

while( $line = <GENELIST>)
{
    chomp $line;	# drop \n at end of line
    #print "$line\n";
    @b=split(/\t/,$line);
    # print "\n$b[0]\t$b[1]\t$ll\n";
    for (my $j=$b[0];$j<=int($b[1]);$j+=3)
    {
       
        if ($b[2]eq"D" && $b[4] eq "CDS")
        {
            my $codon=substr($refsequence, $b[0]+int(($j-$b[0])/3)*3, 3);
            $Codonlist{$codon}++;
            
        }
        
        if ($b[2]eq"C" && $b[4] eq "CDS")
        {
            my $codon=substr($refsequence, $b[1]-int(($b[1]-$j)/3)*3-2, 3);
            
            # print("$codon\n");
            $codon=reverse($codon);
            $codon=~tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
            $Codonlist{$codon}++;
            
        }
        
        
        
    }
    # print ">",$b[3],"\n",substr($refsequence,$b[0],$b[1]-$b[0]+1),"\n";
    $nb++;
}

close GENELIST;


my $sum=0;
foreach my $i (keys %Codonlist)
{
    #print "$i\t",$Codonlist{$i},"\n";
    $sum+=$Codonlist{$i};
}
print "number of codons in CDS $sum\n";


my @MutType=(1,2,3,4,5,6);
my %Codon_Mut_S=();
my %Codon_Mut_N=();
my @codonpos=(0,1,2);
print "@MutType\n";

foreach my $i (@MutType)
{
    foreach my $codon (keys %Codonlist)
    {
        my $N=0;
        my $S=0;
        $Codon_Mut_S{"$codon"."_$i"}="_";
        $Codon_Mut_N{"$codon"."_$i"}="_";
        foreach my $j (@codonpos)
        {
            my $newcodon=$codon;
            my $aa=$CodeGenetic{$codon};
            my $newaa=$aa;
            
            my $baseref=substr($codon,$j,1);
            my $mutname="$i"."_$baseref";
            if (exists $MutTypeAction{$mutname})
            {
                substr($newcodon, $j, 1)=$MutTypeAction{$mutname};
                $newaa=$CodeGenetic{$newcodon};
                
                if ($newaa eq $aa)
                {
                    $S++;
                }else
                {
                    $N++;
                }
                
                #here compute the number of AA change
                
            
            }
        }
        $Codon_Mut_S{"$codon"."_$i"}=$S;
        $Codon_Mut_N{"$codon"."_$i"}=$N;
        
    }
    
    
}
print "Codon\tQuantity\tS1_AT_to_CG\tS2_AT_to_GC\tS3_AT_to_TA\tS4_CG_to_GC\tS5_CG_to_TA\tS6_GC_to_TA\tN1_AT_to_CG\tN2_AT_to_GC\tN3_AT_to_TA\tN4_CG_to_GC\tN5_CG_to_TA\tN6_GC_to_TA\n";

my @Synonymous=(0,0,0,0,0,0);
my @NonSynonymous=(0,0,0,0,0,0);

foreach my $i (keys %Codonlist)
{
    print "$i\t$Codonlist{$i}\t";
    foreach my $j (@MutType)
    {
    print $Codon_Mut_S{"$i"."_$j"},"\t";
    $Synonymous[$j-1]+=$Codon_Mut_S{"$i"."_$j"}*$Codonlist{$i};
    }
    foreach my $j (@MutType)
    {
    print $Codon_Mut_N{"$i"."_$j"},"\t";
    $NonSynonymous[$j-1]+=$Codon_Mut_N{"$i"."_$j"}*$Codonlist{$i};

    }
    print "\n";
    
}
print "Synonymous: @Synonymous\n";
print "Non-Synonymous: @NonSynonymous\n";

open (GENELIST, "<ObservedMuts.txt");

$line = <GENELIST>;
chomp $line;	# drop \n at end of line
@b=split(/\t/,$line);
close GENELIST;

#HERE loose programming to be more explicit:
my $expetedNS_muttype_1=($b[0]/$Synonymous[0])*$NonSynonymous[0];
my $expetedNS_muttype_2=($b[1]/$Synonymous[1])*$NonSynonymous[1];
my $expetedNS_muttype_3=($b[2]/$Synonymous[2])*$NonSynonymous[2];
my $expetedNS_muttype_4=($b[3]/$Synonymous[3])*$NonSynonymous[3];
my $expetedNS_muttype_5=($b[4]/$Synonymous[4])*$NonSynonymous[4];
my $expetedNS_muttype_6=($b[5]/$Synonymous[5])*$NonSynonymous[5];

my $expectedNS=($expetedNS_muttype_1+$expetedNS_muttype_2+$expetedNS_muttype_3+$expetedNS_muttype_4+$expetedNS_muttype_5+$expetedNS_muttype_6);
my $observedNS=($b[6]+$b[7]+$b[8]+$b[9]+$b[10]+$b[11]);




my $KaKs=$observedNS/$expectedNS;


print "KA/KS= $KaKs\n"



my $expetedS_muttype_1=($b[6]/$NonSynonymous[0])*$Synonymous[0];
my $expetedS_muttype_2=($b[7]/$NonSynonymous[1])*$Synonymous[1];
my $expetedS_muttype_3=($b[8]/$NonSynonymous[2])*$Synonymous[2];
my $expetedS_muttype_4=($b[9]/$NonSynonymous[3])*$Synonymous[3];
my $expetedS_muttype_5=($b[10]/$NonSynonymous[4])*$Synonymous[4];
my $expetedS_muttype_6=($b[11]/$NonSynonymous[5])*$Synonymous[5];

my $expectedS=($expetedS_muttype_1+$expetedS_muttype_2+$expetedS_muttype_3+$expetedS_muttype_4+$expetedS_muttype_5+$expetedS_muttype_6);
my $observedS=($b[0]+$b[1]+$b[2]+$b[3]+$b[4]+$b[5]);




$KaKs=1/($observedS/$expectedS);

