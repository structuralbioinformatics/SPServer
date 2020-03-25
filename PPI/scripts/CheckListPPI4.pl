#!/usr/bin/perl
#


 use Getopt::Long;
 use File::Basename;


 &GetOptions("f=s"   => \$file,
             "c=s"   => \$cutoff,
             "cac=s" => \$cutoff_acc,
             "pot=s" => \$pot,
             "fw=s"  => \$file_weight,
             "r=s"   => \$nran,
             "l=s"   => \$list,
             "o=s"   => \$outext,
             "sc=s"  => \$outscore,
             "xml"   => \$xml,
             "h"     => \$help       );

 unless (defined $file  && not defined $help && defined $nran && defined $pdb && $cutoff >= 2.0){
   printf "
          f    File with Potentials
          c    Cut-off for interacting pairs. It should be real, larger than 2 and smaller than 50
          cac  Cut-off for accumulated residue pairs (Potential independent of distances)
          pot  Potential can be by CA, CB or Minimum atom-atom distance: ca,cb,min.
               Default is min.
          fw   Force Field for Weights: File with precomputed environment weights
          r    Number of random decoys
          l    List of PDB complex structures
	  sc   Score of which you wish to have an output list for all PDB complexes
	       Select among combinations: 
	       [D, C, M, U, L] - [Score, Energy, ZEnergy, Dual, Eaa3Denv, Elocal, Ecmp, E3Denv, E3D]
	       For example: D-Score, C-Eaa3Denv, etc.
          xml  Flag to write output in XML format
          o    Output file. Used also as suffix extension for the rest of output files
          \n";
   exit;
  }
  if (not defined $pot){$pot="min";}
  if ($cutoff>50.0){$cutoff=50.0;}

  $outglobal="$outext".".global";
  open (OUTSCR,">$outglobal");
  if (defined $xml){
    $outxml="$outext".".xml";
    open (OUTXML,">$outxml");
    printf OUTXML "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n";
    printf OUTXML "<xml>\n";
 }

  $|=1;
  $boltz=8.31441e-3;
  $T=298.0;
  
$code[0]="ALA";
$code[1]="CYS";
$code[2]="ASP";
$code[3]="GLU";
$code[4]="PHE";
$code[5]="GLY";
$code[6]="HIS";
$code[7]="ILE";
$code[8]="LYS";
$code[9]="LEU";
$code[10]="MET";
$code[11]="ASN";
$code[12]="PRO";
$code[13]="GLN";
$code[14]="ARG";
$code[15]="SER";
$code[16]="THR";
$code[17]="VAL";
$code[18]="TRP";
$code[19]="TYR";
$code[20]="HSD";
$code[21]="HSE";

$lode="ACDEFGHIKLMNPQRSTVWYHHX";
$sode="acdefghiklmnpqrstvwyhhX";
$pode="nnppnnpnpnnpnppppnnpppX";

 $dssp="$ENV{SP4I}"."/bin/dssp";
 $interact="$ENV{SP4I}"."/bin/interact";
 $split="perl "."$ENV{SP4I}"."/scripts/PDBtoAddChain.pl";

# Read Potentials

 $skip=0;
 if (-e $file) {printf "Potential file $file found\n";}else{printf "$file not found\n";}
 open (INP,"<$file") || exit();  
   while(<INP>){
     if (/END/){$skip=0;}
     if (/^PAIR-POTENTIALS OF RESIDUES  \n/){$skip=1;next;}
     if (/^PAIR-POTENTIALS OF ENVIRONMENTS/){$skip=2;next;}
     if (/^LOCAL PAIR-POTENTIALS OF RESIDUES AND ENVIRONMENTS/){$skip=3;next;}
     if (/^PAIR-POTENTIALS OF RESIDUES AND ENVIRONMENTS  \n/){$skip=4;next;}
     if (/^STANDARD PAIR-POTENTIALS OF DISTANCES/){$skip=5;next;}
     if (/^SOLUTION FOR ENVIRONMENTAL WEIGHTS/){$skip=6;next;}
     if (/^PAIR-POTENTIALS OF ACCUMULATED PAIRS/){$skip=7;next;}
     if ($skip==1){ if (/PAIR/){chop;@data=split;$pair=$data[1]; for ($i=2;$i<=$#data;$i++){ $respmf{$pair}[$i-2]=$data[$i];}} }
     if ($skip==2){ if (/PAIR/){chop;@data=split;$pair=$data[1]; for ($i=2;$i<=$#data;$i++){ $epmf{$pair}[$i-2]=$data[$i];}} }
     if ($skip==3){ if (/AA-ENV/){chop;@data=split;$pair=$data[1]; $rpmf{$pair}=$data[2];} }
     if ($skip==4){ if (/PAIR/){chop;@data=split;$pair=$data[1]; for ($i=2;$i<=$#data;$i++){ $pmf{$pair}[$i-2]=$data[$i];}} }
     if ($skip==5){ if (/Distance/){chop;@data=split;$d=$data[1]; $dpmf[$d]=$data[2];} }
     if ($skip==6){ if (/Pair/){
                                chop;
                               ($pair,$factor)=split /=>/,substr($_,10,100);
                                $pair=~s/ //g;
                               ($env1,$env2)=split /:/,$pair;
                                $pair_reverse="$env2".":"."$env1";
                                $weight{$pair}=$factor/(-1*$boltz*$T);
                                $weight{$pair_reverse}=$factor/(-1*$boltz*$T);
                               } }
     if ($skip==7){ if (/PAIR/){
        chop;
        @data=split;
        $pair=$data[1];
        ($pair1,$pair2)=split /:/,$pair;
        $pair_reverse="$pair2".":"."$pair1";
        $i=2+$cutoff_acc;
        $pmfacc{$pair}=$data[$i];
        $pmfacc{$pair_reverse}=$data[$i];
     } }
   }
  close INP;

  if (-e $file_weight) {printf "Weights file $file_weight found\n";}else{printf "Weight File $file_weight not found\n";}

  open (INP,"<$file_weight");
   while(<INP>){
     if (/END/){$skip=0;}
     if (/^SOLUTION FOR ENVIRONMENTAL WEIGHTS/){$skip=6;next;}
     if ($skip==6){ if (/Pair/){
                                chop;
                               ($pair,$factor)=split /=>/,substr($_,10,100);
                                $pair=~s/ //g;
                               ($env1,$env2)=split /:/,$pair;
                                $pair_reverse="$env2".":"."$env1";
                                $weight{$pair}=$factor/(-1*$boltz*$T);
                                $weight{$pair_reverse}=$factor/(-1*$boltz*$T);
                               } }
   }
  close INP;




    undef %Mrespmf;
  foreach $pair (keys %respmf){
    $Mrespmf{$pair}=0.0;
    $n=0;
    for $d (0..$#dpmf){ $Mrespmf{$pair} += $respmf{$pair}[$d];$n++;}
    $Mrespmf{$pair}=$Mrespmf{$pair}/$n;
    for $d (0..$#dpmf){ $Mrespmf{$pair}[$d]=$Mrespmf{$pair};}
  }
    undef %Mepmf;
  foreach $pair (keys %epmf){
    $Mepmf{$pair}=0.0;
    $n=0;
    for $d (0..$#dpmf){ $Mepmf{$pair} += $epmf{$pair}[$d];$n++;}
    $Mrespmf{$pair}=$Mepmf{$pair}/$n;
    for $d (0..$#dpmf){ $Mepmf{$pair}[$d]=$Mepmf{$pair};}
  }
    undef %Mpmf;
  foreach $pair (keys %pmf){
    $Mpmf{$pair}=0.0;
    $n=0;
    for $d (0..$#dpmf){ $Mpmf{$pair} += $pmf{$pair}[$d];$n++;}
    $Mpmf{$pair}=$Mpmf{$pair}/$n;
    for $d (0..$#dpmf){ $Mpmf{$pair}[$d]=$Mpmf{$pair};}
  }

  undef %Lrespmf;
  foreach $pair (keys %respmf){
    $Lrespmf{$pair}=0.0;
    $n=0;
    for $d (0..$#dpmf){ if ($Lrespmf{$pair}>$respmf{$pair}[$d]){$Lrespmf{$pair}=$respmf{$pair}[$d];}}
    for $d (0..$#dpmf){ $Lrespmf{$pair}[$d]=$Lrespmf{$pair};}
  }
  undef %Lepmf;
  foreach $pair (keys %epmf){
    $Lepmf{$pair}=0.0;
    $n=0;
    for $d (0..$#dpmf){ if ($Lepmf{$pair}<$epmf{$pair}[$d]){$Lepmf{$pair}=$epmf{$pair}[$d];}}
    for $d (0..$#dpmf){ $Lepmf{$pair}[$d]=$Lepmf{$pair};}
  }
  undef %Lpmf;
  foreach $pair (keys %pmf){
    $Lpmf{$pair}=0.0;
    $n=0;
    for $d (0..$#dpmf){ if ($Lpmf{$pair}>$pmf{$pair}[$d]){$Lpmf{$pair}=$pmf{$pair}[$d];}}
    for $d (0..$#dpmf){ $Lpmf{$pair}[$d]=$Lpmf{$pair};}
  }


  undef %Urespmf;
  foreach $pair (keys %respmf){
    $Urespmf{$pair}=0.0;
    $n=0;
    for $d (0..$#dpmf){ if ($Urespmf{$pair}<$respmf{$pair}[$d]){$Urespmf{$pair}=$respmf{$pair}[$d];}}
    for $d (0..$#dpmf){ $Urespmf{$pair}[$d]=$Urespmf{$pair};}
  }
  undef %Uepmf;
  foreach $pair (keys %epmf){
    $Uepmf{$pair}=0.0;
    $n=0;
    for $d (0..$#dpmf){ if ($Uepmf{$pair}>$epmf{$pair}[$d]){$Uepmf{$pair}=$epmf{$pair}[$d];}}
    for $d (0..$#dpmf){ $Uepmf{$pair}[$d]=$Uepmf{$pair};}
  }
  undef %Upmf;
  foreach $pair (keys %pmf){
    $Upmf{$pair}=0.0;
    $n=0;
    for $d (0..$#dpmf){ if ($Upmf{$pair}<$pmf{$pair}[$d]){$Upmf{$pair}=$pmf{$pair}[$d];}}
    for $d (0..$#dpmf){ $Upmf{$pair}[$d]=$Upmf{$pair};}
  }

  undef %Crespmf;
  foreach $pair (keys %respmf){
    for $d (0..$cutoff_acc)         { $Crespmf{$pair}[$d]=$respmf{$pair}[$cutoff_acc];}
    for $d (($cutoff_acc+1)..$#dpmf){ $Crespmf{$pair}[$d]=0.0;}
  }
  undef %Cepmf;
  foreach $pair (keys %epmf){
    for $d (0..$cutoff_acc)         { $Cepmf{$pair}[$d]=$epmf{$pair}[$cutoff_acc];}
    for $d (($cutoff_acc+1)..$#dpmf){ $Cepmf{$pair}[$d]=0.0;}
  }
  undef %Cpmf;
  foreach $pair (keys %pmf){
    for $d (0..$cutoff_acc)         { $Cpmf{$pair}[$d]=$pmf{$pair}[$cutoff_acc];}
    for $d (($cutoff_acc+1)..$#dpmf){ $Cpmf{$pair}[$d]=0.0;}
  }
  undef @Cdpmf;
  for $d (0..$cutoff_acc)         { $Cdpmf[$d]=$dpmf[$cutoff_acc];}
  for $d (($cutoff_acc+1)..$#dpmf){ $Cdpmf[$d]=0.0;}


# Read PDB list

 undef @pdb;
 open (LIST, "<$list");
 while(<LIST>){
  chop;
  @line=split;
  push @pdb,$line[0];
 }
# Read PDB file (starting loop on @PDB files)
 for $pdb (@pdb){

 $output="$pdb"."$outext";
 $np=0; 
 undef @prta;
 undef @prtb;
 $pdb_chain="$pdb".".chain";
 $pdbA="$pdb_chain"."A.pdb";
 $pdbB="$pdb_chain"."B.pdb";
 $dsspf[0]="$pdbA".".dssp";
 $dsspf[1]="$pdbB".".dssp";
 system("$split -i $pdb -o $pdb_chain");
 system("$dssp $pdbA $dsspf[0] ");
 system("$dssp $pdbB $dsspf[1] ");
 printf "$dssp $pdbA $dsspf[0] \n";
 printf "$dssp $pdbB $dsspf[1] \n";
 undef %in;
 $in={};
 $in={
      file  => $pdb,
      chain => "A",
     };
 @prta=&ReadPDB(%in);
 $in={};
 $in={
      file  => $pdb,
      chain => "B",
     };
 @prtb=&ReadPDB(%in);

  undef @ss1;
  undef @ss2;
  undef @bur1;
  undef @bur2;
  undef %out;
 $in={};
 $in={
      file  => $dsspf[0],
      chain => "A",
     };
 %out=&ReadDSSP(%in);
 @ss1=@{$out->{ss}};
 @bur1=@{$out->{bur}};

 undef %out;
 $in={};
 $in={
      file  => $dsspf[1],
      chain => "B",
     };
 %out=&ReadDSSP(%in);
 @ss2=@{$out->{ss}};
 @bur2=@{$out->{bur}};

#TRANSFORM INTO SEQUENCE

  $natoma=$#prta+1;
  $natomb=$#prtb+1;
  undef @seq1;
  undef @seq2;
  undef @pol1;
  undef @pol2; 
  $nresa=0;
  for $i (0..$#prta)
   {
     $res=  $prta[$i]->{residue};
     $letter=22;
     for $j (0..$#code){if ($code[$j] eq $res){$letter=$j;last;}}
     if ($prta[$i]->{atom} eq "CA" || $prta[$i]->{atom} eq "CAA"){ 
         push @seq1,substr($lode,$letter,1);push @pol1,substr($pode,$letter,1); $nresa++; 
     }
   }

  if ($natoma <= 4*$nresa || (1+$#ss1) ne $nresa){
   printf "REJECT $pdb NAT= $natoma NRES= $nresa SS $#ss1 \nSEQUENCE\n @seq1\n @ss1\n";
   exit;
  }
  $nresb=0;
  for $i (0..$#prtb)
   {
     $res=  $prtb[$i]->{residue};
     $letter=22;
     for $j (0..$#code){if ($code[$j] eq $res){$letter=$j;last;}}
     if ($prtb[$i]->{atom} eq "CA" || $prtb[$i]->{atom} eq "CAA"){ 
         push @seq2,substr($lode,$letter,1);push @pol2,substr($pode,$letter,1); $nresb++; 
     }
   }

  if ($natomb <= 4*$nresb || (1+$#ss2) ne $nresb){
   printf "REJECT $pdb NAT= $natomb NRES= $nresb SS $#ss2 \nSEQUENCE\n @seq2\n @ss2\n";
   exit;
  }

# Read alignment 

  undef @querya;
  undef @queryb;
  undef @query1;
  undef @query2;
  undef @seqa;
  undef @seqb;
  undef @qola;
  undef @qolb;
  undef @tsa;
  undef @tsb;
  undef %qsa;
  undef %qsb;

  if (open (ALN,"<$aln1") ){
   undef @sequence1;
   $skip=0;
   $seq="";
   while (<ALN>){
    if (/^>/ && $skip ==0){$skip=1;$seq="";next;}
    if (/^>/ && $skip ==1){push @sequence1,$seq;$seq="";next;}
    if ($skip == 1){chop;$seq.="$_";}
   }
   push @sequence1,$seq;
   close ALN;
  }else{
   undef @sequence1;
   push @sequence1,join '',@seq1;
   $seq1=join '',@seq1;
   push @sequence1,join '',@seq1;
  }
  $found=0;
  for $seq (@sequence1){
    $seq=~s/\s//g;
    undef @string;
    undef @word;
    @word=split /-/,$seq; 
    for $word (@word){ push @string,$word; } 
    $test    =join '',@string;
    $template=join '',@seq1;
    if ($test eq $template && $found==0){
       $n=0;
       @word=split //,$seq;
       for $word (@word){ 
            push @seqa,$word;
            if ($word ne "-"){push @tsa,$n;}
            $n++;
        }
       $found=1;
    }else{
       $n=0;
       @word=split //,$seq;
       for $word (@word){ 
            push @querya,$word;
            if ($word ne "-"){
               for($ka=0;$ka<22;$ka++){if ($word eq substr($lode,$ka,1)){$nk=$ka;last;}}
               $pol=substr($pode,$nk,1);
               push @qola,$pol;
               $qsa{$n}=1;
               push @query1,$word;
            }else{
               push @qola,"-";
               $qsa{$n}=0;
            }
            $n++;
        }
    }
  }
  if (open (ALN,"<$aln2") ){
   undef @sequence2;
   $skip=0;
   while (<ALN>){
    if (/^>/ && $skip ==0){$skip=1;$seq="";next;}
    if (/^>/ && $skip ==1){push @sequence2,$seq;$seq="";next;}
    if ($skip == 1){chop;$seq.="$_";}
   }
   push @sequence2,$seq;
   close ALN;
  }else{
   undef @sequence2;
   push @sequence2,join '',@seq2;
   push @sequence2,join '',@seq2;
  }
  $found=0;
  for $seq (@sequence2){
    $seq=~s/\s//g;
    undef @string;
    undef @word;
    @word=split /-/,$seq; 
    for $word (@word){ push @string,$word; } 
    $test    =join '',@string;
    $template=join '',@seq2;
    if ($test eq $template && $found==0){
       $n=0;
       @word=split //,$seq;
       for $word (@word){ 
            push @seqb,$word;
            if ($word ne "-"){push @tsb,$n;}
            $n++;
        }
        $found=1;
    }else{
       $n=0;
       @word=split //,$seq;
       for $word (@word){ 
            push @queryb,$word;
            if ($word ne "-"){
               for($ka=0;$ka<22;$ka++){if ($word eq substr($lode,$ka,1)){$nk=$ka;last;}}
               $pol=substr($pode,$nk,1);
               push @qolb,$pol;
               $qsb{$n}=1;
               push @query2,$word;
            }else{
               push @qolb,"-";
               $qsb{$n}=0;
            }
            $n++;
        }
    }
  }

# CALCULATE DISTANCES
 undef @Energy;
 undef @Energy2;
 undef @Eaa3Denv;
 undef @E3D;
 undef @E3Denv;
 undef @Elocal;
 undef @Pscore;
 undef @Ecmp;

 undef @aEnergy;
 undef @aEnergy2;
 undef @aEaa3Denv;
 undef @aE3D;
 undef @aE3Denv;
 undef @aElocal;
 undef @aPscore;
 undef @aEcmp;

 undef @bEnergy;
 undef @bEnergy2;
 undef @bEaa3Denv;
 undef @bE3D;
 undef @bE3Denv;
 undef @bElocal;
 undef @bPscore;
 undef @bEcmp;

 undef @abEnergy2;
 undef @abEnergy;
 undef @abEaa3Denv;
 undef @abE3D;
 undef @abE3Denv;
 undef @abElocal;
 undef @abPscore;
 undef @abEcmp;

 undef @CEaa3Denv;
 undef @CE3D;
 undef @CE3Denv;
 undef @CElocal;
 undef @CPscore;
 undef @CEcmp;

 undef @aCEaa3Denv;
 undef @aCE3D;
 undef @aCE3Denv;
 undef @aCElocal;
 undef @aCPscore;
 undef @aCEcmp;

 undef @bCEaa3Denv;
 undef @bCE3D;
 undef @bCE3Denv;
 undef @bCElocal;
 undef @bCPscore;
 undef @bCEcmp;

 undef @aMEaa3Denv;
 undef @aME3Denv;
 undef @aMPscore;
 undef @bMEaa3Denv;
 undef @bME3Denv;
 undef @bMPscore;
 $MEaa3Denv=0;
 $ME3Denv=0;
 $MPscore=0;
 $MEnergy=0;
 $MEnergy2=0;

 undef @aLEaa3Denv;
 undef @aLE3Denv;
 undef @aLPscore;
 undef @bLEaa3Denv;
 undef @bLE3Denv;
 undef @bLPscore;
 $LEaa3Denv=0;
 $LE3Denv=0;
 $LPscore=0;
 $LEnergy=0;
 $LEnergy2=0;

 undef @aUEaa3Denv;
 undef @aUE3Denv;
 undef @aUPscore;
 undef @bUEaa3Denv;
 undef @bUE3Denv;
 undef @bUPscore;
 $UEaa3Denv=0;
 $UE3Denv=0;
 $UPscore=0;
 $UEnergy=0;
 $UEnergy2=0;


   $pdbDF="$pdb".".interact"."$pot";
   print "Running $interact with $pdbA and $pdbB \n";
   $cutoff2= 1 + $cutoff;
   if (-e $pdbDF){
        printf "FOUND  $pdbDF after previous running interact -a $pdbA -b $pdbB -p $pot -c $cutoff2\n";
    }else{
        system("$interact -a $pdbA -b $pdbB -p $pot -c $cutoff2 -o $pdbDF");
        printf "$interact -a $pdbA -b $pdbB -p $pot -c $cutoff2 -o $pdbDF\n";
    }

 open (DIST,"<$pdbDF") || next;
   $skip=0;
   while(<DIST>){
     if ($skip==1){last;}
     ($n,$m,$dr,$aaa,$aab)=split;
     $i=$n-1;
     $j=$m-1;
     $testa=$seq1[$i];
     $testb=$seq2[$j];
     $seqa=$querya[$tsa[$i]];
     $seqb=$queryb[$tsb[$j]];
     if ($seqa eq "-" || $seqb eq "-") {next;}
     $pola=$qola[$tsa[$i]];
     $polb=$qolb[$tsb[$j]];
     $ss1 =$ss1[$i];
     $ss2 =$ss2[$j];
     $bur1=$bur1[$i];
     $bur2=$bur2[$j];
     $pair1= $seqa."-".$pola."-".$ss1."-".$bur1;
     $pair2= $seqb."-".$polb."-".$ss2."-".$bur2;
     $env1="$pola"."-"."$ss1"."-"."$bur1";
     $env2="$polb"."-"."$ss2"."-"."$bur2";
     if (not defined $seqa || not defined $seqb){
       printf "Error(1) (N,M)=($n,$m) pair in $pdb\n";
       $skip=1;
       next;
     }
     if ($testa ne $aaa || $testb ne $aab){
       printf "Error(1) (N,M)=($n,$m) pair($aaa,$aab) not ($testa,$testb)  in $pdb\n";
       $skip=1;
       next;
     }
     if ( $pair1 eq "---" || $pair2 eq "---"){
       printf "Error(2) (N,M)=($n,$m) pair in $pdb\n";
       $skip=1;
       next;
     }
     if ( $env1 eq "--" || $env2 eq "--"){
       printf "Error(3) (N,M)=($n,$m) pair in $pdb\n";
       $skip=1;
       next;
     }
     if ( $seqa eq "X" || $seqb eq "X"){
       printf "Unknown residue X (N,M)=($n,$m) ($seqa,$seqb) pair in $pdb\n";
       next;
     }
     for($ka=0;$ka<22;$ka++){if ($seqa eq substr($lode,$ka,1)){$aa1=$ka;last;}}
     for($ka=0;$ka<22;$ka++){if ($seqb eq substr($lode,$ka,1)){$aa2=$ka;last;}}
     if ($aa1<=$aa2){
          $respair="$seqa".":"."$seqb";
          $pair = "$pair1".":"."$pair2";
          $envpair="$env1".":"."$env2";
          $envpair_reverse="$env2".":"."$env1";
         }else{
          $respair="$seqb".":"."$seqa";
          $pair = "$pair2".":"."$pair1";
          $envpair="$env2".":"."$env1";
          $envpair_reverse="$env1".":"."$env2";
         } 
     $d=int($dr);

     $aPscore[$np][$i]  +=$respmf{$respair}[$d];
     $aEaa3Denv[$np][$i]+=$pmf{$pair}[$d];
     $abPscore[$np][$i][$j]  +=$respmf{$respair}[$d];
     $abEaa3Denv[$np][$i][$j]+=$pmf{$pair}[$d];

     $aMPscore[$i]  +=$Mrespmf{$respair}[$d];
     $aMEaa3Denv[$i]+=$Mpmf{$pair}[$d];
     $aLPscore[$i]  +=$Lrespmf{$respair}[$d];
     $aLEaa3Denv[$i]+=$Lpmf{$pair}[$d];
     $aUPscore[$i]  +=$Urespmf{$respair}[$d];
     $aUEaa3Denv[$i]+=$Upmf{$pair}[$d];

     $aE3D[$np][$i]     +=$dpmf[$d];
     $abE3D[$np][$i][$j]     +=$dpmf[$d];

     if     (defined $epmf{$envpair}[$d])         {$aE3Denv[$np][$i]  +=$epmf{$envpair}[$d];}
     elsif  (defined $epmf{$envpair_reverse}[$d]) {$aE3Denv[$np][$i]  +=$epmf{$envpair_reverse}[$d];}
 
     if     (defined $epmf{$envpair}[$d])         {$abE3Denv[$np][$i][$j]  +=$epmf{$envpair}[$d];}
     elsif  (defined $epmf{$envpair_reverse}[$d]) {$abE3Denv[$np][$i][$j]   +=$epmf{$envpair_reverse}[$d];}

     if     (defined $Cepmf{$envpair}[$d])         {$aCE3Denv[$np][$i]  +=$Cepmf{$envpair}[$d];}
     elsif  (defined $Cepmf{$envpair_reverse}[$d]) {$aCE3Denv[$np][$i]  +=$Cepmf{$envpair_reverse}[$d];}

     if     (defined $Mepmf{$envpair}[$d])         {$aME3Denv[$i]  +=$Mepmf{$envpair}[$d];}
     elsif  (defined $Mepmf{$envpair_reverse}[$d]) {$aME3Denv[$i]  +=$Mepmf{$envpair_reverse}[$d];}

     if     (defined $Uepmf{$envpair}[$d])         {$aUE3Denv[$i]  +=$Uepmf{$envpair}[$d];}
     elsif  (defined $Uepmf{$envpair_reverse}[$d]) {$aUE3Denv[$i]  +=$Uepmf{$envpair_reverse}[$d];}

     if     (defined $Lepmf{$envpair}[$d])         {$aLE3Denv[$i]  +=$Lepmf{$envpair}[$d];}
     elsif  (defined $Lepmf{$envpair_reverse}[$d]) {$aLE3Denv[$i]  +=$Lepmf{$envpair_reverse}[$d];}



#     $aElocal[$np][$i]  +=$rpmf{$pair1}+$rpmf{$pair2};
     $aElocal[$np][$i]  +=$rpmf{$pair1};
     $aEcmp[$np][$i]    +=$boltz*$T*$weight{$envpair};
     $abElocal[$np][$i][$j]  +=$rpmf{$pair1};
     $abEcmp[$np][$i][$j]    +=$boltz*$T*$weight{$envpair};


#     $aEnergy[$np][$i]   =$aEaa3Denv[$np][$i] + $aE3D[$np][$i] - $aElocal[$np][$i] - $aE3Denv[$np][$i] - $aEcmp[$np][$i];

     if     (defined $epmf{$envpair}[$d]){
     $aEnergy[$np][$i]  +=$pmf{$pair}[$d] + $dpmf[$d] - $rpmf{$pair1} - $rpmf{$pair2} - $epmf{$envpair}[$d] - $boltz*$T*$weight{$envpair};
     $abEnergy[$np][$i][$j]  +=$pmf{$pair}[$d] + $dpmf[$d] - 0.5*$rpmf{$pair1} - 0.5*$rpmf{$pair2} - $epmf{$envpair}[$d] - $boltz*$T*$weight{$envpair};
     $aEnergy2[$np][$i]  +=$pmf{$pair}[$d]  - $epmf{$envpair}[$d] ;
     $abEnergy2[$np][$i][$j]  +=$pmf{$pair}[$d]  - $epmf{$envpair}[$d] ;
     }else{
     $aEnergy[$np][$i]  +=$pmf{$pair}[$d] + $dpmf[$d] - $rpmf{$pair1} - $rpmf{$pair2} - $epmf{$envpair_reverse}[$d] - $boltz*$T*$weight{$envpair};
     $abEnergy[$np][$i][$j]  +=$pmf{$pair}[$d] + $dpmf[$d] - 0.5*$rpmf{$pair1} - 0.5*$rpmf{$pair2} - $epmf{$envpair_reverse}[$d] - $boltz*$T*$weight{$envpair};
     $aEnergy2[$np][$i]  +=$pmf{$pair}[$d]  - $epmf{$envpair_reverse}[$d] ;
     $abEnergy2[$np][$i][$j]  +=$pmf{$pair}[$d]  - $epmf{$envpair_reverse}[$d] ;
     }

     if     (defined $Mepmf{$envpair}[$d]){
     $aMEnergy[$i]  +=$Mpmf{$pair}[$d] + $dpmf[$d] - $rpmf{$pair1} - $rpmf{$pair2} - $Mepmf{$envpair}[$d] - $boltz*$T*$weight{$envpair};
     }else{
     $aMEnergy[$i]  +=$Mpmf{$pair}[$d] + $dpmf[$d] - $rpmf{$pair1} - $rpmf{$pair2} - $Mepmf{$envpair_reverse}[$d] - $boltz*$T*$weight{$envpair};
     }
     if     (defined $Uepmf{$envpair}[$d]){
     $aUEnergy[$i]  +=$Upmf{$pair}[$d] + $dpmf[$d] - $rpmf{$pair1} - $rpmf{$pair2} - $Uepmf{$envpair}[$d] - $boltz*$T*$weight{$envpair};
     }else{
     $aUEnergy[$i]  +=$Upmf{$pair}[$d] + $dpmf[$d] - $rpmf{$pair1} - $rpmf{$pair2} - $Uepmf{$envpair_reverse}[$d] - $boltz*$T*$weight{$envpair};
     }
     if     (defined $Lepmf{$envpair}[$d]){
     $aLEnergy[$i]  +=$Lpmf{$pair}[$d] + $dpmf[$d] - $rpmf{$pair1} - $rpmf{$pair2} - $Lepmf{$envpair}[$d] - $boltz*$T*$weight{$envpair};
     }else{
     $aLEnergy[$i]  +=$Lpmf{$pair}[$d] + $dpmf[$d] - $rpmf{$pair1} - $rpmf{$pair2} - $Lepmf{$envpair_reverse}[$d] - $boltz*$T*$weight{$envpair};
     }

     if ($d<=$cutoff_acc) {
        $aCPscore[$np][$i]  +=$Crespmf{$respair}[$d];
        $aCEaa3Denv[$np][$i]+=$Cpmf{$pair}[$d];
        $aCE3D[$np][$i]     +=$Cdpmf[$d];
        $aCElocal[$np][$i]  +=$rpmf{$pair1};
        $aCEcmp[$np][$i]    +=$boltz*$T*$weight{$envpair};
        if     (defined $Cepmf{$envpair}[$d]){
         $aCEnergy[$np][$i]  +=$Cpmf{$pair}[$d] + $Cdpmf[$d] - $rpmf{$pair1} - $rpmf{$pair2} - $Cepmf{$envpair}[$d] - $boltz*$T*$weight{$envpair};
         $aCE3Denv[$np][$i]  +=$Cepmf{$envpair}[$d];
         }else{
         $aCEnergy[$np][$i]  +=$Cpmf{$pair}[$d] + $Cdpmf[$d] - $rpmf{$pair1} - $rpmf{$pair2} - $Cepmf{$envpair_reverse}[$d] - $boltz*$T*$weight{$envpair};
         $aCE3Denv[$np][$i]  +=$Cepmf{$envpair_reverse}[$d];
         }
     }


     $bPscore[$np][$j]  +=$respmf{$respair}[$d];
     $bEaa3Denv[$np][$j]+=$pmf{$pair}[$d];

     $bMPscore[$j]  +=$Mrespmf{$respair}[$d];
     $bMEaa3Denv[$j]+=$Mpmf{$pair}[$d];
     $bLPscore[$j]  +=$Lrespmf{$respair}[$d];
     $bLEaa3Denv[$j]+=$Lpmf{$pair}[$d];
     $bUPscore[$j]  +=$Urespmf{$respair}[$d];
     $bUEaa3Denv[$j]+=$Upmf{$pair}[$d];

     $bE3D[$np][$j]     +=$dpmf[$d];

     if     (defined $epmf{$envpair}[$d])         {$bE3Denv[$np][$j]  +=$epmf{$envpair}[$d];}
     elsif  (defined $epmf{$envpair_reverse}[$d]) {$bE3Denv[$np][$j]  +=$epmf{$envpair_reverse}[$d];}

     if     (defined $Mepmf{$envpair}[$d])         {$bME3Denv[$i]  +=$Mepmf{$envpair}[$d];}
     elsif  (defined $Mepmf{$envpair_reverse}[$d]) {$bME3Denv[$i]  +=$Mepmf{$envpair_reverse}[$d];}

     if     (defined $Uepmf{$envpair}[$d])         {$bUE3Denv[$i]  +=$Uepmf{$envpair}[$d];}
     elsif  (defined $Uepmf{$envpair_reverse}[$d]) {$bUE3Denv[$i]  +=$Uepmf{$envpair_reverse}[$d];}

     if     (defined $Lepmf{$envpair}[$d])         {$bLE3Denv[$i]  +=$Lepmf{$envpair}[$d];}
     elsif  (defined $Lepmf{$envpair_reverse}[$d]) {$bLE3Denv[$i]  +=$Lepmf{$envpair_reverse}[$d];}


#     $bElocal[$np][$j]  +=$rpmf{$pair1}+$rpmf{$pair2};
     $bElocal[$np][$j]  +=$rpmf{$pair2};
     $bEcmp[$np][$j]    +=$boltz*$T*$weight{$envpair};
     if     (defined $epmf{$envpair}[$d]){
     $bEnergy[$np][$j]  +=$pmf{$pair}[$d] + $dpmf[$d] - $rpmf{$pair1} - $rpmf{$pair2} - $epmf{$envpair}[$d] - $boltz*$T*$weight{$envpair};
     $bEnergy2[$np][$j]  +=$pmf{$pair}[$d]  - $epmf{$envpair}[$d] ;
     }else{
     $bEnergy[$np][$j]  +=$pmf{$pair}[$d] + $dpmf[$d] - $rpmf{$pair1} - $rpmf{$pair2} - $epmf{$envpair_reverse}[$d] - $boltz*$T*$weight{$envpair};
     $bEnergy2[$np][$j]  +=$pmf{$pair}[$d]  - $epmf{$envpair_reverse}[$d] ;
     }

     if     (defined $Mepmf{$envpair}[$d]){
     $bMEnergy[$j]  +=$Mpmf{$pair}[$d] + $dpmf[$d] - $rpmf{$pair1} - $rpmf{$pair2} - $Mepmf{$envpair}[$d] - $boltz*$T*$weight{$envpair};
     }else{
     $bMEnergy[$j]  +=$Mpmf{$pair}[$d] + $dpmf[$d] - $rpmf{$pair1} - $rpmf{$pair2} - $Mepmf{$envpair_reverse}[$d] - $boltz*$T*$weight{$envpair};
     }
     if     (defined $Uepmf{$envpair}[$d]){
     $bUEnergy[$j]  +=$Upmf{$pair}[$d] + $dpmf[$d] - $rpmf{$pair1} - $rpmf{$pair2} - $Uepmf{$envpair}[$d] - $boltz*$T*$weight{$envpair};
     }else{
     $bUEnergy[$j]  +=$Upmf{$pair}[$d] + $dpmf[$d] - $rpmf{$pair1} - $rpmf{$pair2} - $Uepmf{$envpair_reverse}[$d] - $boltz*$T*$weight{$envpair};
     }
     if     (defined $Lepmf{$envpair}[$d]){
     $bLEnergy[$j]  +=$Lpmf{$pair}[$d] + $dpmf[$d] - $rpmf{$pair1} - $rpmf{$pair2} - $Lepmf{$envpair}[$d] - $boltz*$T*$weight{$envpair};
     }else{
     $bLEnergy[$j]  +=$Lpmf{$pair}[$d] + $dpmf[$d] - $rpmf{$pair1} - $rpmf{$pair2} - $Lepmf{$envpair_reverse}[$d] - $boltz*$T*$weight{$envpair};
     }

    if ($d<=$cutoff_acc) {
        $bCPscore[$np][$j]  +=$Crespmf{$respair}[$d];
        $bCEaa3Denv[$np][$j]+=$Cpmf{$pair}[$d];
        $bCE3D[$np][$j]     +=$Cdpmf[$d];
        $bCElocal[$np][$j]  +=$rpmf{$pair2};
        $bCEcmp[$np][$j]    +=$boltz*$T*$weight{$envpair};
        if     (defined $Cepmf{$envpair}[$d]){
         $bCEnergy[$np][$j]  +=$Cpmf{$pair}[$d] + $Cdpmf[$d] - $rpmf{$pair1} - $rpmf{$pair2} - $Cepmf{$envpair}[$d] - $boltz*$T*$weight{$envpair};
         $bCE3Denv[$np][$j]  +=$Cepmf{$envpair}[$d];
         }else{
         $bCEnergy[$np][$j]  +=$Cpmf{$pair}[$d] + $Cdpmf[$d] - $rpmf{$pair1} - $rpmf{$pair2} - $Cepmf{$envpair_reverse}[$d] - $boltz*$T*$weight{$envpair};
         $bCE3Denv[$np][$j]  +=$Cepmf{$envpair_reverse}[$d];
         }
     }



     $Pscore[$np]  +=$respmf{$respair}[$d];
     $Eaa3Denv[$np]+=$pmf{$pair}[$d];
     $MPscore  +=$Mrespmf{$respair}[$d];
     $MEaa3Denv+=$Mpmf{$pair}[$d];
     $LPscore  +=$Lrespmf{$respair}[$d];
     $LEaa3Denv+=$Lpmf{$pair}[$d];
     $UPscore  +=$Urespmf{$respair}[$d];
     $UEaa3Denv+=$Upmf{$pair}[$d];


     $E3D[$np]     +=$dpmf[$d];

     if     (defined $epmf{$envpair}[$d])         {$E3Denv[$np]  +=$epmf{$envpair}[$d];}
     elsif  (defined $epmf{$envpair_reverse}[$d]) {$E3Denv[$np]  +=$epmf{$envpair_reverse}[$d];}

     if     (defined $Mepmf{$envpair}[$d])         {$ME3Denv  +=$Mepmf{$envpair}[$d];}
     elsif  (defined $Mepmf{$envpair_reverse}[$d]) {$ME3Denv  +=$Mepmf{$envpair_reverse}[$d];}

     if     (defined $Uepmf{$envpair}[$d])         {$UE3Denv  +=$Uepmf{$envpair}[$d];}
     elsif  (defined $Uepmf{$envpair_reverse}[$d]) {$UE3Denv  +=$Uepmf{$envpair_reverse}[$d];}

     if     (defined $Lepmf{$envpair}[$d])         {$LE3Denv  +=$Lepmf{$envpair}[$d];}
     elsif  (defined $Lepmf{$envpair_reverse}[$d]) {$LE3Denv  +=$Lepmf{$envpair_reverse}[$d];}


     $Elocal[$np]  +=$rpmf{$pair1}+$rpmf{$pair2};
     $Ecmp[$np]    +=$boltz*$T*$weight{$envpair};
#     $Energy[$np]   =$Eaa3Denv[$np] + $E3D[$np] - $Elocal[$np] - $E3Denv[$np] - $Ecmp[$np];

     if     (defined $epmf{$envpair}[$d]){
     $Energy[$np]  +=$pmf{$pair}[$d] + $dpmf[$d] - $rpmf{$pair1} - $rpmf{$pair2} - $epmf{$envpair}[$d] - $boltz*$T*$weight{$envpair};
     $Energy2[$np]  +=$pmf{$pair}[$d]  - $epmf{$envpair}[$d] ;
     }else{
     $Energy[$np]  +=$pmf{$pair}[$d] + $dpmf[$d] - $rpmf{$pair1} - $rpmf{$pair2} - $epmf{$envpair_reverse}[$d] - $boltz*$T*$weight{$envpair};
     $Energy2[$np]  +=$pmf{$pair}[$d]  - $epmf{$envpair_reverse}[$d] ;
     }
     if     (defined $Mepmf{$envpair}[$d]){
     $MEnergy  +=$Mpmf{$pair}[$d] + $dpmf[$d] - $rpmf{$pair1} - $rpmf{$pair2} - $Mepmf{$envpair}[$d] - $boltz*$T*$weight{$envpair};
     $MEnergy2  +=$Mpmf{$pair}[$d] - $Mepmf{$envpair}[$d] ;
     }else{
     $MEnergy  +=$Mpmf{$pair}[$d] + $dpmf[$d] - $rpmf{$pair1} - $rpmf{$pair2} - $Mepmf{$envpair_reverse}[$d] - $boltz*$T*$weight{$envpair};
     $MEnergy2  +=$Mpmf{$pair}[$d]  - $Mepmf{$envpair_reverse}[$d] ;
     }
     if     (defined $Uepmf{$envpair}[$d]){
     $UEnergy  +=$Upmf{$pair}[$d] + $dpmf[$d] - $rpmf{$pair1} - $rpmf{$pair2} - $Uepmf{$envpair}[$d] - $boltz*$T*$weight{$envpair};
     $UEnergy2  +=$Upmf{$pair}[$d] - $Uepmf{$envpair}[$d] ;
     }else{
     $UEnergy  +=$Upmf{$pair}[$d] + $dpmf[$d] - $rpmf{$pair1} - $rpmf{$pair2} - $Uepmf{$envpair_reverse}[$d] - $boltz*$T*$weight{$envpair};
     $UEnergy2  +=$Upmf{$pair}[$d] - $Uepmf{$envpair_reverse}[$d] ;
     }
     if     (defined $Lepmf{$envpair}[$d]){
     $LEnergy  +=$Lpmf{$pair}[$d] + $dpmf[$d] - $rpmf{$pair1} - $rpmf{$pair2} - $Lepmf{$envpair}[$d] - $boltz*$T*$weight{$envpair};
     $LEnergy2  +=$Lpmf{$pair}[$d] - $Lepmf{$envpair}[$d] ;
     }else{
     $LEnergy  +=$Lpmf{$pair}[$d] + $dpmf[$d] - $rpmf{$pair1} - $rpmf{$pair2} - $Lepmf{$envpair_reverse}[$d] - $boltz*$T*$weight{$envpair};
     $LEnergy2  +=$Lpmf{$pair}[$d]  - $Lepmf{$envpair_reverse}[$d] ;
     }


     if ($d<=$cutoff_acc) {
        $CPscore[$np]  +=$Crespmf{$respair}[$d];
        $CEaa3Denv[$np]+=$Cpmf{$pair}[$d];
        $CE3D[$np]     +=$Cdpmf[$d];
        $CElocal[$np]  +=$rpmf{$pair1}+$rpmf{$pair2};
        $CEcmp[$np]    +=$boltz*$T*$weight{$envpair};
        if     (defined $Cepmf{$envpair}[$d]){
         $CEnergy[$np]  +=$Cpmf{$pair}[$d] + $Cdpmf[$d] - $rpmf{$pair1} - $rpmf{$pair2} - $Cepmf{$envpair}[$d] - $boltz*$T*$weight{$envpair};
         $CEnergy2[$np] +=$Cpmf{$pair}[$d] - $Cepmf{$envpair}[$d];
         $CE3Denv[$np]  +=$Cepmf{$envpair}[$d];
         }else{
         $CEnergy[$np]  +=$Cpmf{$pair}[$d] + $Cdpmf[$d] - $rpmf{$pair1} - $rpmf{$pair2} - $Cepmf{$envpair_reverse}[$d] - $boltz*$T*$weight{$envpair};
         $CEnergy2[$np] +=$Cpmf{$pair}[$d] - $Cepmf{$envpair_reverse}[$d];
         $CE3Denv[$np]  +=$Cepmf{$envpair_reverse}[$d];
         }
      }


    }
   close DIST;
 
   $np++;

   printf "RANDOM DECOYS\n";
   $nnp=0;
   for ($iran=0;$iran<$nran;$iran++){
     if ($nnp == $nran/1000 && $nran>=1000){printf ".";}else{printf ".";}
     if ($nnp == $nran/100  && $nran>=100) {printf "$iran";}
     if ($nnp == $nran/10   && $nran>=10)  {#printf "\n";
                                            $nnp=0;}
     $nnp++;
     undef @seqaT;
     undef @seqaR;
     undef @polaR;
     @seqaT=@query1;
     for $nn (0..$#querya){
        if ($qsa{$nn} == 1){
          $nr=rand $#seqaT;
          $aa=splice(@seqaT,$nr,1);
          for($ka=0;$ka<22;$ka++){if ($aa eq substr($lode,$ka,1)){$nk=$ka;last;}}
          $pol=substr($pode,$nk,1);
        }else{
          $aa=$pol="-";
        }
        push @seqaR,$aa;
        push @polaR,$pol;
     }
     undef @seqbT;
     undef @seqbR;
     undef @polbR;
     @seqbT=@query2;
     for $nn (0..$#queryb){
        if ($qsb{$nn} == 1){
          $nr=rand $#seqbT;
          $aa=splice(@seqbT,$nr,1);
          for($ka=0;$ka<22;$ka++){if ($aa eq substr($lode,$ka,1)){$nk=$ka;last;}}
          $pol=substr($pode,$nk,1);
        }else{
          $aa=$pol="-";
        }
        push @seqbR,$aa;
        push @polbR,$pol;
     }
     $Eaa3Denv[$np]=0;
     $E3D[$np]=0;
     $E3Denv[$np]=0;
     $Elocal[$np]=0;
     $Pscore[$np]=0;
     $Ecmp[$np]=0;
     $CEaa3Denv[$np]=0;
     $CE3D[$np]=0;
     $CE3Denv[$np]=0;
     $CElocal[$np]=0;
     $CPscore[$np]=0;
     $CEcmp[$np]=0;

     $pdbDF="$pdb".".interact"."$pot";
     open (DIST,"<$pdbDF") || next;
       $skip=0;
       while(<DIST>){
         if ($skip==1){last;}
         ($n,$m,$dr,$aaa,$aab)=split;
         $i=$n-1;
         $j=$m-1;
         $testa=$seq1[$i];
         $testb=$seq2[$j];
         $seqa=$seqaR[$tsa[$i]];
         $seqb=$seqbR[$tsb[$j]];
         if ($seqa eq "-" || $seqb eq "-") {next;}
         $pola=$polaR[$tsa[$i]];
         $polb=$polbR[$tsb[$j]];
         $ss1 =$ss1[$i];
         $ss2 =$ss2[$j];
         $bur1=$bur1[$i];
         $bur2=$bur2[$j];
         $pair1= $seqa."-".$pola."-".$ss1."-".$bur1;
         $pair2= $seqb."-".$polb."-".$ss2."-".$bur2;
         $env1="$pola"."-"."$ss1"."-"."$bur1";
         $env2="$polb"."-"."$ss2"."-"."$bur2";
          if (not defined $seqa || not defined $seqb){
            printf "Error(1) (N,M)=($n,$m) pair in $pdb\n";
            $skip=1;
            next;
          }
          if ($testa ne $aaa || $testb ne $aab){
             printf "Error(1) (N,M)=($n,$m) pair($aaa,$aab) not ($testa,$testb)  in $pdb\n";
             $skip=1;
             next;
          }
          if ( $pair1 eq "---" || $pair2 eq "---"){
            printf "Error(2) (N,M)=($n,$m) pair in $pdb\n";
            $skip=1;
            next;
          }
          if ( $env1 eq "--" || $env2 eq "--"){
            printf "Error(3) (N,M)=($n,$m) pair in $pdb\n";
            $skip=1;
            next;
          }
          if ( $seqa eq "X" || $seqb eq "X"){
           printf "Unknown residue X (N,M)=($n,$m) ($seqa,$seqb) pair in $pdb\n";
           next;
          }
         for($ka=0;$ka<22;$ka++){if ($seqa eq substr($lode,$ka,1)){$aa1=$ka;last;}}
         for($ka=0;$ka<22;$ka++){if ($seqb eq substr($lode,$ka,1)){$aa2=$ka;last;}}
         if ($aa1<=$aa2){
          $pair = "$pair1".":"."$pair2";
          $respair="$seqa".":"."$seqb";
          $envpair="$env1".":"."$env2";
          $envpair_reverse="$env2".":"."$env1";
         }else{
          $pair = "$pair2".":"."$pair1";
          $respair="$seqb".":"."$seqa";
          $envpair="$env2".":"."$env1";
          $envpair_reverse="$env1".":"."$env2";
         } 
         $d=int($dr);

          $aPscore[$np][$i]  +=$respmf{$respair}[$d];
          $aEaa3Denv[$np][$i]+=$pmf{$pair}[$d];
          $aE3D[$np][$i]     +=$dpmf[$d];

          $abPscore[$np][$i][$j]  +=$respmf{$respair}[$d];
          $abEaa3Denv[$np][$i][$j]+=$pmf{$pair}[$d];
          $abE3D[$np][$i][$j]     +=$dpmf[$d];

          if     (defined $epmf{$envpair}[$d])         {$aE3Denv[$np][$i]  +=$epmf{$envpair}[$d];}
          elsif  (defined $epmf{$envpair_reverse}[$d]) {$aE3Denv[$np][$i]  +=$epmf{$envpair_reverse}[$d];}
          if     (defined $epmf{$envpair}[$d])         {$abE3Denv[$np][$i][$j]  +=$epmf{$envpair}[$d];}
          elsif  (defined $epmf{$envpair_reverse}[$d]) {$abE3Denv[$np][$i][$j]  +=$epmf{$envpair_reverse}[$d];}
          $aElocal[$np][$i]  +=$rpmf{$pair1};
          $aEcmp[$np][$i]    +=$boltz*$T*$weight{$envpair};
          $aEnergy[$np][$i]   =$aEaa3Denv[$np][$i] + $aE3D[$np][$i] - $aElocal[$np][$i] - $aE3Denv[$np][$i] - $aEcmp[$np][$i];
          $aEnergy2[$np][$i]   =$aEaa3Denv[$np][$i]  - $aE3Denv[$np][$i] ;

          $abElocal[$np][$i][$j]  +=$rpmf{$pair1};
          $abEcmp[$np][$i][$j]    +=$boltz*$T*$weight{$envpair};
          $abEnergy[$np][$i][$j]   =$abEaa3Denv[$np][$i][$j] + $abE3D[$np][$i][$j] - $abElocal[$np][$i][$j] - $abE3Denv[$np][$i][$j] - $abEcmp[$np][$i][$j];
          $abEnergy2[$np][$i][$j]   =$abEaa3Denv[$np][$i][$j]  - $abE3Denv[$np][$i][$j] ;


          $bPscore[$np][$j]  +=$respmf{$respair}[$d];
          $bEaa3Denv[$np][$j]+=$pmf{$pair}[$d];
          $bE3D[$np][$j]     +=$dpmf[$d];
          if     (defined $epmf{$envpair}[$d])         {$bE3Denv[$np][$j]  +=$epmf{$envpair}[$d];}
          elsif  (defined $epmf{$envpair_reverse}[$d]) {$bE3Denv[$np][$j]  +=$epmf{$envpair_reverse}[$d];}
          $bElocal[$np][$j]  +=$rpmf{$pair2};
          $bEcmp[$np][$j]    +=$boltz*$T*$weight{$envpair};
          $bEnergy[$np][$j]   =$bEaa3Denv[$np][$j] + $bE3D[$np][$j] - $bElocal[$np][$j] - $bE3Denv[$np][$j] - $bEcmp[$np][$j];
          $bEnergy2[$np][$j]   =$bEaa3Denv[$np][$j] - $bE3Denv[$np][$j] ;

          $Pscore[$np]  +=$respmf{$respair}[$d];
          $Eaa3Denv[$np]+=$pmf{$pair}[$d];
          $E3D[$np]     +=$dpmf[$d];
          if     (defined $epmf{$envpair}[$d])         {$E3Denv[$np]  +=$epmf{$envpair}[$d];}
          elsif  (defined $epmf{$envpair_reverse}[$d]) {$E3Denv[$np]  +=$epmf{$envpair_reverse}[$d];}
          $Elocal[$np]  +=$rpmf{$pair1}+$rpmf{$pair2};
          $Ecmp[$np]    +=$boltz*$T*$weight{$envpair};
          $Energy[$np]   =$Eaa3Denv[$np] + $E3D[$np] - $Elocal[$np] - $E3Denv[$np] - $Ecmp[$np];
          $Energy2[$np]   =$Eaa3Denv[$np] - $E3Denv[$np] ;

          if ($d<=$cutoff_acc) {
             $CPscore[$np]  +=$Crespmf{$respair}[$d];
             $CEaa3Denv[$np]+=$Cpmf{$pair}[$d];
             $CE3D[$np]     +=$Cdpmf[$d];
             $CElocal[$np]  +=$rpmf{$pair1}+$rpmf{$pair2};
             $CEcmp[$np]    +=$boltz*$T*$weight{$envpair};
             $aCPscore[$np][$i]  +=$Crespmf{$respair}[$d];
             $aCEaa3Denv[$np][$i]+=$Cpmf{$pair}[$d];
             $aCE3D[$np][$i]     +=$Cdpmf[$d];
             $aCElocal[$np][$i]  +=$rpmf{$pair1};
             $aCEcmp[$np][$i]    +=$boltz*$T*$weight{$envpair};
             $bCPscore[$np][$j]  +=$Crespmf{$respair}[$d];
             $bCEaa3Denv[$np][$j]+=$Cpmf{$pair}[$d];
             $bCE3D[$np][$j]     +=$Cdpmf[$d];
             $bCElocal[$np][$j]  +=$rpmf{$pair2};
             $bCEcmp[$np][$j]    +=$boltz*$T*$weight{$envpair};
             if     (defined $Cepmf{$envpair}[$d]){
              $CEnergy[$np]  +=$Cpmf{$pair}[$d] + $Cdpmf[$d] - $rpmf{$pair1} - $rpmf{$pair2} - $Cepmf{$envpair}[$d] - $boltz*$T*$weight{$envpair};
              $CE3Denv[$np]  +=$Cepmf{$envpair}[$d];
              $aCEnergy[$np][$i]  +=$Cpmf{$pair}[$d] + $Cdpmf[$d] - $rpmf{$pair1} - $rpmf{$pair2} - $Cepmf{$envpair}[$d] - $boltz*$T*$weight{$envpair};
              $aCE3Denv[$np][$i]  +=$Cepmf{$envpair}[$d];
              $bCEnergy[$np][$j]  +=$Cpmf{$pair}[$d] + $Cdpmf[$d] - $rpmf{$pair1} - $rpmf{$pair2} - $Cepmf{$envpair}[$d] - $boltz*$T*$weight{$envpair};
              $bCE3Denv[$np][$j]  +=$Cepmf{$envpair}[$d];
              $CEnergy2[$np] +=$Cpmf{$pair}[$d] - $Cepmf{$envpair}[$d];
              }else{
              $CEnergy[$np]  +=$Cpmf{$pair}[$d] + $Cdpmf[$d] - $rpmf{$pair1} - $rpmf{$pair2} - $Cepmf{$envpair_reverse}[$d] - $boltz*$T*$weight{$envpair};
              $CE3Denv[$np]  +=$Cepmf{$envpair_reverse}[$d];
              $aCEnergy[$np][$i]  +=$Cpmf{$pair}[$d] + $Cdpmf[$d] - $rpmf{$pair1} - $rpmf{$pair2} - $Cepmf{$envpair_reverse}[$d] - $boltz*$T*$weight{$envpair};
              $aCE3Denv[$np][$i]  +=$Cepmf{$envpair_reverse}[$d];
              $bCEnergy[$np][$j]  +=$Cpmf{$pair}[$d] + $Cdpmf[$d] - $rpmf{$pair1} - $rpmf{$pair2} - $Cepmf{$envpair_reverse}[$d] - $boltz*$T*$weight{$envpair};
              $bCE3Denv[$np][$j]  +=$Cepmf{$envpair_reverse}[$d];
              $CEnergy2[$np] +=$Cpmf{$pair}[$d] - $Cepmf{$envpair_reverse}[$d];
              }
           }          

        }
     close DIST;
     $np++;

   } #En random test

  printf "\n";
# Calculate Z-scores

  undef @aZEnergy;
  undef @aZEnergy2;
  undef @aZZEnergy;
  undef @aZEaa3Denv;
  undef @aZElocal;
  undef @aZE3Denv;
  undef @aZEcmp;
  undef @aZscore;


  undef @abZEnergy;
  undef @abZEnergy2;
  undef @abZZEnergy;
  undef @abZEaa3Denv;
  undef @abZElocal;
  undef @abZE3Denv;
  undef @abZEcmp;
  undef @abZscore;


  undef @bZEnergy;
  undef @bZZEnergy;
  undef @bZEaa3Denv;
  undef @bZElocal;
  undef @bZE3Denv;
  undef @bZEcmp;
  undef @bZscore;
  $ZEnergy=0;
  $ZZEnergy=0;
  $ZEnergy2=0;
  $Zaa3Denv=0;
  $ZElocal=0;
  $ZE3Denv=0;
  $ZEcmp=0;
  $Zscore=0;

  undef @aCZEnergy;
  undef @aCZZEnergy;
  undef @aCZEaa3Denv;
  undef @aCZElocal;
  undef @aCZE3Denv;
  undef @aCZEcmp;
  undef @aCZscore;
  undef @bCZEnergy;
  undef @bCZZEnergy;
  undef @bCZEaa3Denv;
  undef @bCZElocal;
  undef @bCZE3Denv;
  undef @bCZEcmp;
  undef @bCZscore;
  $CZEnergy=0;
  $CZEnergy2=0;
  $CZZEnergy=0;
  $CZaa3Denv=0;
  $CZElocal=0;
  $CZE3Denv=0;
  $CZEcmp=0;
  $CZscore=0;

  undef @aMZEnergy;
  undef @aMZZEnergy;
  undef @aMZEnergy2;
  undef @aMZEaa3Denv;
  undef @aMZE3Denv;
  undef @aMZscore;
  undef @bMZEnergy;
  undef @bMZZEnergy;
  undef @bMZEaa3Denv;
  undef @bMZE3Denv;
  undef @bMZscore;
  $MZEnergy=0;
  $MZZEnergy=0;
  $MZEnergy2=0;
  $MZaa3Denv=0;
  $MZE3Denv=0;
  $MZscore=0;
  undef @aLZEnergy;
  undef @aLZZEnergy;
  undef @aLZEnergy2;
  undef @aLZEaa3Denv;
  undef @aLZE3Denv;
  undef @aLZscore;
  undef @bLZZEnergy;
  undef @bLZEnergy;
  undef @bLZEaa3Denv;
  undef @bLZE3Denv;
  undef @bLZscore;
  $LZEnergy=0;
  $LZZEnergy=0;
  $LZEnergy2=0;
  $LZaa3Denv=0;
  $LZE3Denv=0;
  $LZscore=0;
  undef @aUZEnergy;
  undef @aUZZEnergy;
  undef @aUZEnergy2;
  undef @aUZEaa3Denv;
  undef @aUZE3Denv;
  undef @aUZscore;
  undef @bUZZEnergy;
  undef @bUZEnergy;
  undef @bUZEaa3Denv;
  undef @bUZE3Denv;
  undef @bUZscore;
  $UZEnergy=0;
  $UZZEnergy=0;
  $UZEnergy2=0;
  $UZaa3Denv=0;
  $UZE3Denv=0;
  $UZscore=0;

 
  for ($i=0;$i<$#seqa;$i++){
  for ($k=0;$k<$#seqb;$k++){

      $x=$abEnergy[0][$i][$k];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$aEnergy[$j+1][$i][$k]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $abZEnergy[$i][$k]=&Zscore(%in);

      $x=$abEnergy2[0][$i][$k];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$abEnergy2[$j+1][$i][$k]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $abZEnergy2[$i][$k]=&Zscore(%in);

      $x=$abEaa3Denv[0][$i][$k];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$abEaa3Denv[$j+1][$i][$k]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $abZEaa3Denv[$i][$k]=&Zscore(%in);

      $x=$abElocal[0][$i];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$abElocal[$j+1][$i][$k]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $abZElocal[$i][$k]=&Zscore(%in);

      $x=$abE3Denv[0][$i][$k];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$abE3Denv[$j+1][$i][$k]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $abZE3Denv[$i][$k]=&Zscore(%in);

      $x=$abEcmp[0][$i][$k];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$abEcmp[$j+1][$i][$k]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $abZEcmp[$i][$k]=&Zscore(%in);

      $x=$abPscore[0][$i][$k];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$abPscore[$j+1][$i][$k]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $abZscore[$i][$k]=&Zscore(%in);

      $abZZEnergy[$i][$k] =$abZEaa3Denv[$i][$k]   - $abZElocal[$i][$k]    - $abZE3Denv[$i][$k];

  }}

  for ($i=0;$i<$#seqa;$i++){
      $x=$aCEnergy[0][$i];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$aCEnergy[$j+1][$i]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $aCZEnergy[$i]=&Zscore(%in);
      $x=$aEnergy[0][$i];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$aEnergy[$j+1][$i]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $aZEnergy[$i]=&Zscore(%in);
      $x=$aMEnergy[$i];
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $aMZEnergy[$i]=&Zscore(%in);
      $x=$aUEnergy[$i];
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $aUZEnergy[$i]=&Zscore(%in);
      $x=$aLEnergy[$i];
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $aLZEnergy[$i]=&Zscore(%in);


      $x=$aEnergy2[0][$i];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$aEnergy2[$j+1][$i]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $aZEnergy2[$i]=&Zscore(%in);
      $x=$aMEnergy2[$i];
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $aMZEnergy2[$i]=&Zscore(%in);
      $x=$aUEnergy2[$i];
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $aUZEnergy2[$i]=&Zscore(%in);
      $x=$aLEnergy2[$i];
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $aLZEnergy2[$i]=&Zscore(%in);


      $x=$aEaa3Denv[0][$i];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$aEaa3Denv[$j+1][$i]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $aZEaa3Denv[$i]=&Zscore(%in);
      $x=$aCEaa3Denv[0][$i];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$aCEaa3Denv[$j+1][$i]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $aCZEaa3Denv[$i]=&Zscore(%in);
      $x=$aMEaa3Denv[$i];
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $aMZEaa3Denv[$i]=&Zscore(%in);
      $x=$aUEaa3Denv[$i];
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $aUZEaa3Denv[$i]=&Zscore(%in);
      $x=$aLEaa3Denv[$i];
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $aLZEaa3Denv[$i]=&Zscore(%in);

      $x=$aElocal[0][$i];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$aElocal[$j+1][$i]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $aZElocal[$i]=&Zscore(%in);
      $x=$aCElocal[0][$i];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$aCElocal[$j+1][$i]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $aCZElocal[$i]=&Zscore(%in);
      $x=$aE3Denv[0][$i];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$aE3Denv[$j+1][$i]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $aZE3Denv[$i]=&Zscore(%in);
      $x=$aCE3Denv[0][$i];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$aCE3Denv[$j+1][$i]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $aCZE3Denv[$i]=&Zscore(%in);
      $x=$aME3Denv[$i];
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $aMZE3Denv[$i]=&Zscore(%in);
      $x=$aUE3Denv[$i];
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $aUZE3Denv[$i]=&Zscore(%in);
      $x=$aLE3Denv[$i];
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $aLZE3Denv[$i]=&Zscore(%in);

      $x=$aEcmp[0][$i];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$aEcmp[$j+1][$i]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $aZEcmp[$i]=&Zscore(%in);

      $x=$aCEcmp[0][$i];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$aCEcmp[$j+1][$i]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $aCZEcmp[$i]=&Zscore(%in);

      $x=$aPscore[0][$i];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$aPscore[$j+1][$i]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $aZscore[$i]=&Zscore(%in);

      $x=$aCPscore[0][$i];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$aCPscore[$j+1][$i]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $aCZscore[$i]=&Zscore(%in);

      $x=$aMPscore[$i];
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $aMZscore[$i]=&Zscore(%in);
      $x=$aUPscore[$i];
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $aUZscore[$i]=&Zscore(%in);
      $x=$aLPscore[$i];
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $aLZscore[$i]=&Zscore(%in);
      $aZZEnergy[$i] =$aZEaa3Denv[$i]   - $aZElocal[$i]    - $aZE3Denv[$i];
      $aCZZEnergy[$i]=$aCZEaa3Denv[$i]  - $aCZElocal[$i]   - $aCZE3Denv[$i];
      $aUZZEnergy[$i]=$aUZEaa3Denv[$i]  - $aUZElocal[$i]   - $aUZE3Denv[$i];
      $aLZZEnergy[$i]=$aLZEaa3Denv[$i]  - $aLZElocal[$i]   - $aLZE3Denv[$i];
      $aMZZEnergy[$i]=$aMZEaa3Denv[$i]  - $aMZElocal[$i]   - $aMZE3Denv[$i];
   }
  for ($i=0;$i<$#seqb;$i++){
      $x=$bEnergy[0][$i];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$bEnergy[$j+1][$i]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $bZEnergy[$i]=&Zscore(%in);
      $x=$bCEnergy[0][$i];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$bCEnergy[$j+1][$i]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $bCZEnergy[$i]=&Zscore(%in);
      $x=$bMEnergy[$i];
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $bMZEnergy[$i]=&Zscore(%in);
      $x=$bUEnergy[$i];
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $bUZEnergy[$i]=&Zscore(%in);
      $x=$bLEnergy[$i];
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $bLZEnergy[$i]=&Zscore(%in);

      $x=$bEnergy2[0][$i];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$bEnergy2[$j+1][$i]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $bZEnergy2[$i]=&Zscore(%in);
      $x=$bMEnergy2[$i];
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $bMZEnergy2[$i]=&Zscore(%in);
      $x=$bUEnergy2[$i];
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $bUZEnergy2[$i]=&Zscore(%in);
      $x=$bLEnergy2[$i];
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $bLZEnergy2[$i]=&Zscore(%in);


      $x=$bEaa3Denv[0][$i];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$bEaa3Denv[$j+1][$i]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $bZEaa3Denv[$i]=&Zscore(%in);
      $x=$bCEaa3Denv[0][$i];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$bCEaa3Denv[$j+1][$i]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $bCZEaa3Denv[$i]=&Zscore(%in);

      $x=$bMEaa3Denv[$i];
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $bMZEaa3Denv[$i]=&Zscore(%in);
      $x=$bUEaa3Denv[$i];
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $bUZEaa3Denv[$i]=&Zscore(%in);
      $x=$bLEaa3Denv[$i];
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $bLZEaa3Denv[$i]=&Zscore(%in);

      $x=$bElocal[0][$i];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$bElocal[$j+1][$i]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $bZElocal[$i]=&Zscore(%in);
      $x=$bCElocal[0][$i];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$bCElocal[$j+1][$i]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $bCZElocal[$i]=&Zscore(%in);

      $x=$bE3Denv[0][$i];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$bE3Denv[$j+1][$i]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $bZE3Denv[$i]=&Zscore(%in);
      $x=$bCE3Denv[0][$i];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$bCE3Denv[$j+1][$i]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $bCZE3Denv[$i]=&Zscore(%in);

      $x=$bME3Denv[$i];
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $bMZE3Denv[$i]=&Zscore(%in);
      $x=$bUE3Denv[$i];
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $bUZE3Denv[$i]=&Zscore(%in);
      $x=$bLE3Denv[$i];
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $bLZE3Denv[$i]=&Zscore(%in);
      $x=$bEcmp[0][$i];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$bEcmp[$j+1][$i]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $bZEcmp[$i]=&Zscore(%in);
      $x=$bCEcmp[0][$i];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$bCEcmp[$j+1][$i]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $bCZEcmp[$i]=&Zscore(%in);

      $x=$bPscore[0][$i];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$bPscore[$j+1][$i]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $bZscore[$i]=&Zscore(%in);
      $x=$bCPscore[0][$i];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$bCPscore[$j+1][$i]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $bCZscore[$i]=&Zscore(%in);

      $x=$bMPscore[$i];
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $bMZscore[$i]=&Zscore(%in);
      $x=$bUPscore[$i];
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $bUZscore[$i]=&Zscore(%in);
      $x=$bLPscore[$i];
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $bLZscore[$i]=&Zscore(%in);
      $bZZEnergy[$i] =$bZEaa3Denv[$i]   - $bZElocal[$i]    - $bZE3Denv[$i];
      $bCZZEnergy[$i]=$bCZEaa3Denv[$i]  - $bCZElocal[$i]   - $bCZE3Denv[$i];
      $bUZZEnergy[$i]=$bUZEaa3Denv[$i]  - $bUZElocal[$i]   - $bUZE3Denv[$i];
      $bLZZEnergy[$i]=$bLZEaa3Denv[$i]  - $bLZElocal[$i]   - $bLZE3Denv[$i];
      $bMZZEnergy[$i]=$bMZEaa3Denv[$i]  - $bMZElocal[$i]   - $bMZE3Denv[$i];


   }

 $eneout="$output".".ene";

    if ( defined $output ) { open (OUTPUT,">$output"); select(OUTPUT);}

      $x=$Energy[0];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$Energy[$j+1]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $ZEnergy=&Zscore(%in);
      $x=$ZEnergy;
      undef @y;
      for ($j=0;$j<$nran;$j++){
          $xx=$Energy[$j+1];
          undef @yy;
          for ($jj=1;$jj<$nran;$jj++){ push @yy, $Energy[$jj];}
          $in={};
          $in={
            x => $xx,
            n => $nran,
            y => [@yy]
          };
          push @y,&Zscore(%in); 
      }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $ZZEnergy=&Zscore(%in);


      $x=$CEnergy[0];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$CEnergy[$j+1]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $CZEnergy=&Zscore(%in);
      $x=$CZEnergy;
      undef @y;
      for ($j=0;$j<$nran;$j++){
          $xx=$CEnergy[$j+1];
          undef @yy;
          for ($jj=1;$jj<$nran;$jj++){ push @yy, $CEnergy[$jj];}
          $in={};
          $in={
            x => $xx,
            n => $nran,
            y => [@yy]
          };
          push @y,&Zscore(%in); 
      }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $CZZEnergy=&Zscore(%in);


      $x=$MEnergy;
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $MZEnergy=&Zscore(%in);

      $x=$UEnergy;
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $UZEnergy=&Zscore(%in);

      $x=$LEnergy;
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $LZEnergy=&Zscore(%in);

                 


      $x=$Energy2[0];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$Energy2[$j+1]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $ZEnergy2=&Zscore(%in);

      $x=$CEnergy2[0];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$CEnergy2[$j+1]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $CZEnergy2=&Zscore(%in);

      $x=$MEnergy2;
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $MZEnergy2=&Zscore(%in);

      $x=$UEnergy2;
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $UZEnergy2=&Zscore(%in);

      $x=$LEnergy2;
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $LZEnergy2=&Zscore(%in);




                 

      $x=$Eaa3Denv[0];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$Eaa3Denv[$j+1]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $ZEaa3Denv=&Zscore(%in);

      $x=$CEaa3Denv[0];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$CEaa3Denv[$j+1]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $CZEaa3Denv=&Zscore(%in);

      $x=$MEaa3Denv;
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $MZEaa3Denv=&Zscore(%in);

      $x=$UEaa3Denv;
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $UZEaa3Denv=&Zscore(%in);

      $x=$LEaa3Denv;
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $LZEaa3Denv=&Zscore(%in);

                 


      $x=$Elocal[0];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$Elocal[$j+1]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $ZElocal=&Zscore(%in);

      $x=$CElocal[0];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$CElocal[$j+1]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $CZElocal=&Zscore(%in);

                 
      $x=$E3Denv[0];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$E3Denv[$j+1]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $ZE3Denv=&Zscore(%in);

      $x=$CE3Denv[0];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$CE3Denv[$j+1]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $CZE3Denv=&Zscore(%in);

      $x=$ME3Denv;
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $MZE3Denv=&Zscore(%in);

      $x=$UE3Denv;
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $UZE3Denv=&Zscore(%in);

      $x=$LE3Denv;
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $LZE3Denv=&Zscore(%in);

      $x=$Ecmp[0];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$Ecmp[$j+1]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $ZEcmp=&Zscore(%in);

      $x=$CEcmp[0];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$CEcmp[$j+1]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $CZEcmp=&Zscore(%in);

      $x=$Pscore[0];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$Pscore[$j+1]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $Zscore=&Zscore(%in);

      $x=$CPscore[0];
      undef @y;
      for ($j=0;$j<$nran;$j++){ push @y,$CPscore[$j+1]; }
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $CZscore=&Zscore(%in);


      $x=$MPscore;
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $MZscore=&Zscore(%in);

      $x=$UPscore;
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $UZscore=&Zscore(%in);


      $x=$LPscore;
      $in={};
      $in={
         x => $x,
         n => $nran,
         y => [@y]
       };
      $LZscore=&Zscore(%in);


   if (defined $xml){
     $pdbid= basename($pdb);
     $pdbid=~s/::/-/g;
     $pdbid=~s/#/_/g;
     printf  OUTXML " <protein>\n";
     printf  OUTXML "  <id>%s</id>\n",$pdbid;
     printf  OUTXML "  <length_a>%5d</length_a>\n",$#seqa+1;
     printf  OUTXML "  <length_b>%5d</length_b>\n",$#seqb+1;
     printf  OUTXML "  <global_energies>\n";
     printf  OUTXML "    <Epair>%10.3e</Epair>\n",$Pscore[0];
     printf  OUTXML "    <Ecomb>%10.3e</Ecomb>\n",$Energy[0];
     printf  OUTXML "    <Es3dc>%10.3e</Es3dc>\n",$Eaa3Denv[0];
     printf  OUTXML "    <Elocal>%10.3e</Elocal>\n",-$Elocal[0];
     printf  OUTXML "    <E3dc>%10.3e</E3dc>\n",-$E3Denv[0];
     printf  OUTXML "    <E3d>%10.3e</E3d>\n",$E3D[0];
     printf  OUTXML "    <ZEpair>%10.3e</ZEpair>\n",$Zscore;
     printf  OUTXML "    <ZEcomb>%10.3e</ZEcomb>\n",$ZEnergy;
     printf  OUTXML "    <ZEs3dc>%10.3e</ZEs3dc>\n",$ZEaa3Denv;
     printf  OUTXML "    <ZElocal>%10.3e</ZElocal>\n",-$ZElocal;
     printf  OUTXML "    <ZE3dc>%10.3e</ZE3dc>\n",-$ZE3Denv;
     printf  OUTXML "    <ZE3d>NA</ZE3d>\n";
     printf  OUTXML "    <Zene>%10.3e</Zene>\n",$ZZEnergy;
     printf  OUTXML "  </global_energies>\n";
     printf  OUTXML "  <residues>\n";
     for ($i=0;$i<=$#seqa;$i++){
     for ($k=0;$k<=$#seqb;$k++){
     if (abs($abEnergy[0][$i][$k])>1.0e-5 ||  abs($abPscore[0][$i][$k])>1.0e-5|| abs($abEaa3Denv[0][$i][$k])>1.0e-5|| abs($abE3Denv[0][$i][$k])>1.0e-5||abs($abElocal[0][$i][$k])>1.0e-5){
       printf  OUTXML "   <residue>\n";
       printf  OUTXML "     <number_x>%10d</number_x>\n",$i+1;
       printf  OUTXML "     <number_y>%10d</number_y>\n",$k+1;
       printf  OUTXML "     <type_x>%1s</type_x>\n",$seqa[$tsa[$i]];
       printf  OUTXML "     <type_y>%1s</type_y>\n",$seqb[$tsb[$k]];
       printf  OUTXML "     <Epair>%10.3e</Epair>\n",$abPscore[0][$i][$k];
       printf  OUTXML "     <Ecomb>%10.3e</Ecomb>\n",$abEnergy[0][$i][$k];
       printf  OUTXML "     <Es3dc>%10.3e</Es3dc>\n",$abEaa3Denv[0][$i][$k];
       printf  OUTXML "     <Elocal>%10.3e</Elocal>\n",-$abElocal[0][$i][$k];
       printf  OUTXML "     <E3dc>%10.3e</E3dc>\n",-$abE3Denv[0][$i][$k];
       printf  OUTXML "     <E3d>%10.3e</E3d>\n",$abE3D[0][$i][$k];
       printf  OUTXML "     <ZEpair>%10.3e</ZEpair>\n",$abZscore[$i][$k];
       printf  OUTXML "     <ZEcomb>%10.3e</ZEcomb>\n",$abZEnergy[$i][$k];
       printf  OUTXML "     <ZEs3dc>%10.3e</ZEs3dc>\n",$abZEaa3Denv[$i][$k];
       printf  OUTXML "     <ZElocal>%10.3e</ZElocal>\n",-$abZElocal[$i][$k];
       printf  OUTXML "     <ZE3dc>%10.3e</ZE3dc>\n",-$abZE3Denv[$i][$k];
       printf  OUTXML "     <ZE3d>NA</ZE3d>\n";
       printf  OUTXML "     <Zene>%10.3e</Zene>\n",$abZZEnergy[$i][$k];
       printf  OUTXML "   </residue>\n";
     }}}
     printf  OUTXML "  </residues>\n";
     printf  OUTXML " </protein>\n";



   }else{

   printf "ENERGIES PER RESIDUE\n"; 
   printf "PROTEIN A\n"; 
   printf "%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n","Residue","Energy","Score","Eaa3Denv","E3D","E3Denv","Elocal","Ecmp","Dual";
   for ($i=0;$i<=$#seqa;$i++){
    printf "%10d\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\n",$i,$aEnergy[0][$i],$aPscore[0][$i],$aEaa3Denv[0][$i],$aE3D[0][$i],-$aE3Denv[0][$i],-$aElocal[0][$i],$aEcmp[0][$i],$aEnergy2[0][$i]; 
   }
   printf "END\n";
   printf "ENERGIES PER RESIDUE\n"; 
   printf "PROTEIN B\n"; 
   printf "%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n","Residue","Energy","Score","Eaa3Denv","E3D","E3Denv","Elocal","Ecmp","Dual";
   for ($i=0;$i<=$#seqb;$i++){
    printf "%10d\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\n",$i,$bEnergy[0][$i],$bPscore[0][$i],$bEaa3Denv[0][$i],$bE3D[0][$i],-$bE3Denv[0][$i],-$bElocal[0][$i],$bEcmp[0][$i],$bEnergy2[0][$i]; 
   }
   printf "END\n";
   printf "ENERGIES PER RESIDUE PAIRS\n"; 
   printf "PROTEIN A & B\n"; 
   printf "%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n","Residue","Energy","Score","Eaa3Denv","E3D","E3Denv","Elocal","Ecmp","Dual";
   for ($i=0;$i<=$#seqb;$i++){
   for ($k=0;$k<=$#seqb;$k++){
    if (abs($abEnergy[0][$i][$k])>1.0e-5 ||  abs($abPscore[0][$i][$k])>1.0e-5|| abs($abEaa3Denv[0][$i][$k])>1.0e-5|| abs($abE3Denv[0][$i][$k])>1.0e-5||abs($abE3D[0][$i][$k])>1.0e-5||abs($abElocal[0][$i][$k])>1.0e-5||abs($abEcmp[0][$i][$k])>1.0e-5||abs($abEnergy2[0][$i][$k])>1.0e-5){
     printf "(%4d,%4d)\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\n",$i,$k,$abEnergy[0][$i][$k],$abPscore[0][$i][$k],$abEaa3Denv[0][$i][$k],$abE3D[0][$i][$k],-$abE3Denv[0][$i][$k],-$abElocal[0][$i][$k],$abEcmp[0][$i][$k],$abEnergy2[0][$i][$k]; 
   }}}
   printf "END\n";


   printf "\nZ-SCORES PER RESIDUE\n"; 
   printf "PROTEIN A\n"; 
   printf "%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n","Residue","Energy","Score","Eaa3Denv","E3D","E3Denv","Elocal","Ecmp","Dual","ZEnergy";
   for ($i=0;$i<=$#seqa;$i++){
    printf "%10d\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\n",$i,$aZEnergy[$i],$aZscore[$i],$aZEaa3Denv[$i],0.0,-$aZE3Denv[$i],-$aZElocal[$i],$aZEcmp[$i],$aZEnergy2[0][$i],$aZZEnergy[$i]; 
   }
   printf "END\n";
   printf "\nZ-SCORES PER RESIDUE\n"; 
   printf "PROTEIN B\n"; 
   printf "%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n","Residue","Energy","Score","Eaa3Denv","E3D","E3Denv","Elocal","Ecmp","Dual","ZEnergy";
   for ($i=0;$i<=$#seqb;$i++){
    printf "%10d\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\n",$i,$bZEnergy[$i],$bZscore[$i],$bZEaa3Denv[$i],0.0,-$bZE3Denv[$i],-$bZElocal[$i],$bZEcmp[$i],$bZEnergy2[0][$i],$bZZEnergy[$i]; 
   }
   printf "END\n";
   printf "\nZ-SCORES PER RESIDUE PAIR\n"; 
   printf "PROTEIN A & B\n"; 
   printf "%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n","Residue","Energy","Score","Eaa3Denv","E3D","E3Denv","Elocal","Ecmp","Dual","ZEnergy";
   for ($i=0;$i<=$#seqa;$i++){
   for ($k=0;$k<=$#seqb;$k++){
    if (abs($abZEnergy[$i][$k])>1.0e-5 ||  abs($abZscore[$i][$k])>1.0e-5|| abs($abZEaa3Denv[$i][$k])>1.0e-5|| abs($abZE3Denv[$i][$k])>1.0e-5||abs($abZEcmp[$i][$k])>1.0e-5||abs($abZEnergy2[$i][$k])>1.0e-5){
     printf "(%4d,%4d)\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\n",$i,$k,$abZEnergy[$i][$k],$abZscore[$i][$k],$abZEaa3Denv[$i][$k],0.0,-$abZE3Denv[$i][$k],-$abZElocal[$i][$k],$abZEcmp[$i][$k],$abZEnergy2[$i][$k],$abZZEnergy[$i][$k]; 
   }}}
   printf "END\n";



   printf "C-ENERGIES PER RESIDUE\n"; 
   printf "PROTEIN A\n"; 
   printf "%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n","Residue","Energy","Score","Eaa3Denv","E3D","E3Denv","Elocal","Ecmp","Dual";
   for ($i=0;$i<=$#seqa;$i++){
    printf "%10d\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\n",$i,$aCEnergy[0][$i],$aCPscore[0][$i],$aCEaa3Denv[0][$i],$aCE3D[0][$i],$aCE3Denv[0][$i],$aCElocal[0][$i],$aCEcmp[0][$i],$aCEnergy[0][$i]; 
   }
   printf "END\n";
   printf "C-ENERGIES PER RESIDUE\n"; 
   printf "PROTEIN B\n"; 
   printf "%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n","Residue","Energy","Score","Eaa3Denv","E3D","E3Denv","Elocal","Ecmp","Dual";
   for ($i=0;$i<=$#seqb;$i++){
    printf "%10d\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\n",$i,$bCEnergy[0][$i],$bCPscore[0][$i],$bCEaa3Denv[0][$i],$bCE3D[0][$i],$bCE3Denv[0][$i],$bCElocal[0][$i],$bCEcmp[0][$i],$bCEnergy[0][$i]; 
   }
   printf "END\n";


   printf "\nCZ-SCORES PER RESIDUE\n"; 
   printf "PROTEIN A\n"; 
   printf "%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n","Residue","Energy","Score","Eaa3Denv","E3D","E3Denv","Elocal","Ecmp","Dual","ZEnergy";
   for ($i=0;$i<=$#seqa;$i++){
    printf "%10d\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\n",$i,$aCZEnergy[$i],$aCZscore[$i],$aCZEaa3Denv[$i],0.0,$aCZE3Denv[$i],$aCZElocal[$i],$aCZEcmp[$i],$aCZEnergy[0][$i],$aCZZEnergy[$i]; 
   }
   printf "END\n";
   printf "\nCZ-SCORES PER RESIDUE\n"; 
   printf "PROTEIN B\n"; 
   printf "%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n","Residue","Energy","Score","Eaa3Denv","E3D","E3Denv","Elocal","Ecmp","Dual","ZEnergy";
   for ($i=0;$i<=$#seqb;$i++){
    printf "%10d\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\n",$i,$bCZEnergy[$i],$bCZscore[$i],$bCZEaa3Denv[$i],0.0,$bCZE3Denv[$i],$bCZElocal[$i],$bCZEcmp[$i],$bCZEnergy[0][$i],$bCZZEnergy[$i]; 
   }
   printf "END\n";



   printf "M-ENERGIES PER RESIDUE\n"; 
   printf "PROTEIN A\n"; 
   printf "%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n","Residue","Energy","Score","Eaa3Denv","E3D","E3Denv","Elocal","Ecmp";
   for ($i=0;$i<=$#seqa;$i++){
    printf "%10d\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\n",$i,$aMEnergy[$i],$aMPscore[$i],$aMEaa3Denv[$i],$aE3D[0][$i],$aME3Denv[$i],$aElocal[0][$i],$aEcmp[0][$i]; 
   }
   printf "END\n";
   printf "M-ENERGIES PER RESIDUE\n"; 
   printf "PROTEIN B\n"; 
   printf "%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n","Residue","Energy","Score","Eaa3Denv","E3D","E3Denv","Elocal","Ecmp";
   for ($i=0;$i<=$#seqb;$i++){
    printf "%10d\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\n",$i,$bMEnergy[$i],$bMPscore[$i],$bMEaa3Denv[$i],$bE3D[0][$i],$bME3Denv[$i],$bElocal[0][$i],$bEcmp[0][$i]; 
   }
   printf "END\n";


   printf "\nMZ-SCORES PER RESIDUE\n"; 
   printf "PROTEIN A\n"; 
   printf "%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n","Residue","Energy","Score","Eaa3Denv","E3D","E3Denv","Elocal","Ecmp","ZEnergy";
   for ($i=0;$i<=$#seqa;$i++){
    printf "%10d\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\n",$i,$aMZEnergy[$i],$aMZscore[$i],$aMZEaa3Denv[$i],0.0,$aMZE3Denv[$i],$aZElocal[$i],$aZEcmp[$i],$aMZZEnergy[$i]; 
   }
   printf "END\n";
   printf "\nMZ-SCORES PER RESIDUE\n"; 
   printf "PROTEIN B\n"; 
   printf "%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n","Residue","Energy","Score","Eaa3Denv","E3D","E3Denv","Elocal","Ecmp","ZEnergy";
   for ($i=0;$i<=$#seqb;$i++){
    printf "%10d\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\n",$i,$bMZEnergy[$i],$bMZscore[$i],$bMZEaa3Denv[$i],0.0,$bMZE3Denv[$i],$bZElocal[$i],$bZEcmp[$i],$bMZZEnergy[$i]; 
   }
   printf "END\n";


   printf "U-ENERGIES PER RESIDUE\n"; 
   printf "PROTEIN A\n"; 
   printf "%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n","Residue","Energy","Score","Eaa3Denv","E3D","E3Denv","Elocal","Ecmp";
   for ($i=0;$i<=$#seqa;$i++){
    printf "%10d\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\n",$i,$aUEnergy[$i],$aUPscore[$i],$aUEaa3Denv[$i],$aE3D[0][$i],$aUE3Denv[$i],$aElocal[0][$i],$aEcmp[0][$i]; 
   }
   printf "END\n";
   printf "U-ENERGIES PER RESIDUE\n"; 
   printf "PROTEIN B\n"; 
   printf "%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n","Residue","Energy","Score","Eaa3Denv","E3D","E3Denv","Elocal","Ecmp";
   for ($i=0;$i<=$#seqb;$i++){
    printf "%10d\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\n",$i,$bUEnergy[$i],$bUPscore[$i],$bUEaa3Denv[$i],$bE3D[0][$i],$bUE3Denv[$i],$bElocal[0][$i],$bEcmp[0][$i]; 
   }
   printf "END\n";


   printf "\nUZ-SCORES PER RESIDUE\n"; 
   printf "PROTEIN A\n"; 
   printf "%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n","Residue","Energy","Score","Eaa3Denv","E3D","E3Denv","Elocal","Ecmp","ZEnergy";
   for ($i=0;$i<=$#seqa;$i++){
    printf "%10d\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\n",$i,$aUZEnergy[$i],$aUZscore[$i],$aUZEaa3Denv[$i],0.0,$aUZE3Denv[$i],$aZElocal[$i],$aZEcmp[$i],$aUZZEnergy[$i]; 
   }
   printf "END\n";
   printf "\nUZ-SCORES PER RESIDUE\n"; 
   printf "PROTEIN B\n"; 
   printf "%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n","Residue","Energy","Score","Eaa3Denv","E3D","E3Denv","Elocal","Ecmp","ZEnergy";
   for ($i=0;$i<=$#seqb;$i++){
    printf "%10d\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\n",$i,$bUZEnergy[$i],$bUZscore[$i],$bUZEaa3Denv[$i],0.0,$bUZE3Denv[$i],$bZElocal[$i],$bZEcmp[$i],$bUZZEnergy[$i]; 
   }
   printf "END\n";


   printf "L-ENERGIES PER RESIDUE\n"; 
   printf "PROTEIN A\n"; 
   printf "%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n","Residue","Energy","Score","Eaa3Denv","E3D","E3Denv","Elocal","Ecmp";
   for ($i=0;$i<=$#seqa;$i++){
    printf "%10d\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\n",$i,$aLEnergy[$i],$aLPscore[$i],$aLEaa3Denv[$i],$aE3D[0][$i],$aLE3Denv[$i],$aElocal[0][$i],$aEcmp[0][$i]; 
   }
   printf "END\n";
   printf "L-ENERGIES PER RESIDUE\n"; 
   printf "PROTEIN B\n"; 
   printf "%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n","Residue","Energy","Score","Eaa3Denv","E3D","E3Denv","Elocal","Ecmp";
   for ($i=0;$i<=$#seqb;$i++){
    printf "%10d\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\n",$i,$bLEnergy[$i],$bLPscore[$i],$bLEaa3Denv[$i],$bE3D[0][$i],$bLE3Denv[$i],$bElocal[0][$i],$bEcmp[0][$i]; 
   }
   printf "END\n";


   printf "\nLZ-SCORES PER RESIDUE\n"; 
   printf "PROTEIN A\n"; 
   printf "%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n","Residue","Energy","Score","Eaa3Denv","E3D","E3Denv","Elocal","Ecmp","ZEnergy";
   for ($i=0;$i<=$#seqa;$i++){
    printf "%10d\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\n",$i,$aLZEnergy[$i],$aLZscore[$i],$aLZEaa3Denv[$i],0.0,$aLZE3Denv[$i],$aZElocal[$i],$aZEcmp[$i],$aLZZEnergy[$i]; 
   }
   printf "END\n";
   printf "\nLZ-SCORES PER RESIDUE\n"; 
   printf "PROTEIN B\n"; 
   printf "%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n","Residue","Energy","Score","Eaa3Denv","E3D","E3Denv","Elocal","Ecmp","ZEnergy";
   for ($i=0;$i<=$#seqb;$i++){
    printf "%10d\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\n",$i,$bLZEnergy[$i],$bLZscore[$i],$bLZEaa3Denv[$i],0.0,$bLZE3Denv[$i],$bZElocal[$i],$bZEcmp[$i],$bLZZEnergy[$i]; 
   }
   printf "END\n";

  }


  open (ENE,">$eneout");
  select(ENE);

  printf "\n";
  printf "%10s\t%10.3e\t%10.3e\n","D-Score   ",$Pscore[0],$Zscore;
  printf "%10s\t%10.3e\t%10.3e\n","D-Energy  ",$Energy[0],$ZEnergy;
  printf "%10s\t%10.3e\t%10.3e\n","D-ZEnergy ",$ZEnergy,$ZZEnergy;
  printf "%10s\t%10.3e\t%10.3e\n","D-Dual    ",$Energy2[0],$ZEnergy2;
  printf "%10s\t%10.3e\t%10.3e\n","D-Eaa3Denv",$Eaa3Denv[0],$ZEaa3Denv;
  printf "%10s\t%10.3e\t%10.3e\n","D-Elocal  ",$Elocal[0],$ZElocal;
  printf "%10s\t%10.3e\t%10.3e\n","D-Ecmp    ",$Ecmp[0],$ZEcmp;
  printf "%10s\t%10.3e\t%10.3e\n","D-E3Denv  ",$E3Denv[0],$ZE3Denv;
  printf "%10s\t%10.3e\t%10.3e\n","D-E3D     ",$E3D[0],0.0;
  $Score{Distance}=$Zscore;
  $Energy{Distance}=$ZEnergy;
  $Energy2{Distance}=$ZEnergy2;
  $ZEnergy{Distance}=$ZZEnergy;
  $Eaa3Denv{Distance}=$ZEaa3Denv;
  $Elocal{Distance}=$ZElocal;
  $Ecmp{Distance}=$ZEcmp;
  $E3Denv{Distance}=$ZE3Denv;
  if ($outscore eq "D-Score") {printf OUTSCR "%s\t%10.3e\t%10.3e\n",$eneout,$Pscore[0],$Zscore;}
  if ($outscore eq "D-Energy") {printf OUTSCR "%s\t%10.3e\t%10.3e\n",$eneout,$Energy[0],$ZEnergy;}
  if ($outscore eq "D-ZEnergy") {printf OUTSCR "%s\t%10.3e\t%10.3e\n",$eneout,$ZEnergy,$ZZEnergy;}
  if ($outscore eq "D-Dual") {printf OUTSCR "%s\t%10.3e\t%10.3e\n",$eneout,$Energy2[0],$ZEnergy2;}
  if ($outscore eq "D-Eaa3Denv") {printf OUTSCR "%s\t%10.3e\t%10.3e\n",$eneout,$Eaa3Denv[0],$ZEaa3Denv;}
  if ($outscore eq "D-Elocal") {printf OUTSCR "%s\t%10.3e\t%10.3e\n",$eneout,$Elocal[0],$ZElocal;}
  if ($outscore eq "D-Ecmp") {printf OUTSCR "%s\t%10.3e\t%10.3e\n",$eneout,$Ecmp[0],$ZEcmp;}
  if ($outscore eq "D-E3Denv") {printf OUTSCR "%s\t%10.3e\t%10.3e\n",$eneout,$E3Denv[0],$ZE3Denv;}
  if ($outscore eq "D-E3D") {printf OUTSCR "%s\t%10.3e\t%10.3e\n",$eneout,$E3D[0],0.0;}

  printf "\n";
  printf "%10s\t%10.3e\t%10.3e\n","C-Score   ",$CPscore[0],$CZscore;
  printf "%10s\t%10.3e\t%10.3e\n","C-Energy  ",$CEnergy[0],$CZEnergy;
  printf "%10s\t%10.3e\t%10.3e\n","C-ZEnergy ",$CZEnergy,$CZZEnergy;
  printf "%10s\t%10.3e\t%10.3e\n","C-Dual    ",$CEnergy2[0],$CZEnergy2;
  printf "%10s\t%10.3e\t%10.3e\n","C-Eaa3Denv",$CEaa3Denv[0],$CZEaa3Denv;
  printf "%10s\t%10.3e\t%10.3e\n","C-Elocal  ",$CElocal[0],$CZElocal;
  printf "%10s\t%10.3e\t%10.3e\n","C-Ecmp    ",$CEcmp[0],$CZEcmp;
  printf "%10s\t%10.3e\t%10.3e\n","C-E3Denv  ",$CE3Denv[0],$CZE3Denv;
  printf "%10s\t%10.3e\t%10.3e\n","C-E3D     ",$CE3D[0],0.0;
  $Score{CutOff}=$CZscore;
  $Energy{CutOff}=$CZEnergy;
  $Energy2{CutOff}=$CZEnergy2;
  $ZEnergy{CutOff}=$CZZEnergy;
  $Eaa3Denv{CutOff}=$CZEaa3Denv;
  $Elocal{CutOff}=$CZElocal;
  $Ecmp{CutOff}=$CZEcmp;
  $E3Denv{CutOff}=$CZE3Denv;
  if ($outscore eq "C-Score") {printf OUTSCR "%s\t%10.3e\t%10.3e\n",$eneout,$CPscore[0],$CZscore;}
  if ($outscore eq "C-Energy") {printf OUTSCR "%s\t%10.3e\t%10.3e\n",$eneout,$CEnergy[0],$CZEnergy;}
  if ($outscore eq "C-ZEnergy") {printf OUTSCR "%s\t%10.3e\t%10.3e\n",$eneout,$CZEnergy,$CZZEnergy;}
  if ($outscore eq "C-Dual") {printf OUTSCR "%s\t%10.3e\t%10.3e\n",$eneout,$CEnergy2[0],$CZEnergy2;}
  if ($outscore eq "C-Eaa3Denv") {printf OUTSCR "%s\t%10.3e\t%10.3e\n",$eneout,$CEaa3Denv[0],$CZEaa3Denv;}
  if ($outscore eq "C-Elocal") {printf OUTSCR "%s\t%10.3e\t%10.3e\n",$eneout,$CElocal[0],$CZElocal;}
  if ($outscore eq "C-Ecmp") {printf OUTSCR "%s\t%10.3e\t%10.3e\n",$eneout,$CEcmp[0],$CZEcmp;}
  if ($outscore eq "C-E3Denv") {printf OUTSCR "%s\t%10.3e\t%10.3e\n",$eneout,$CE3Denv[0],$CZE3Denv;}
  if ($outscore eq "C-E3D") {printf OUTSCR "%s\t%10.3e\t%10.3e\n",$eneout,$CE3D[0],0.0;}


  printf "\n";
  printf "%10s\t%10.3e\t%10.3e\n","M-Score   ",$MPscore,$MZscore;
  printf "%10s\t%10.3e\t%10.3e\n","M-Energy  ",$MEnergy,$MZEnergy;
  printf "%10s\t%10.3e\t%10.3e\n","M-Dual    ",$MEnergy2,$MZEnergy2;
  printf "%10s\t%10.3e\t%10.3e\n","M-Eaa3Denv",$MEaa3Denv,$MZEaa3Denv;
  printf "%10s\t%10.3e\t%10.3e\n","M-Elocal  ",$Elocal[0],$ZElocal;
  printf "%10s\t%10.3e\t%10.3e\n","M-Ecmp    ",$Ecmp[0],$ZEcmp;
  printf "%10s\t%10.3e\t%10.3e\n","M-E3Denv  ",$ME3Denv,$MZE3Denv;
  printf "%10s\t%10.3e\t%10.3e\n","M-E3D     ",$E3D[0],0.0;
  $Score{Media}=$MZscore;
  $Energy{Media}=$MZEnergy;
  $Energy2{Media}=$MZEnergy2;
  $Eaa3Denv{Media}=$MZEaa3Denv;
  $Elocal{Media}=$ZElocal;
  $Ecmp{Media}=$ZEcmp;
  $E3Denv{Media}=$MZE3Denv;
  if ($outscore eq "M-Score") {printf OUTSCR "%s\t%10.3e\t%10.3e\n",$eneout,$MPscore[0],$MZscore;}
  if ($outscore eq "M-Energy") {printf OUTSCR "%s\t%10.3e\t%10.3e\n",$eneout,$MEnergy[0],$MZEnergy;}
  if ($outscore eq "M-Dual") {printf OUTSCR "%s\t%10.3e\t%10.3e\n",$eneout,$MEnergy2[0],$MZEnergy2;}
  if ($outscore eq "M-Eaa3Denv") {printf OUTSCR "%s\t%10.3e\t%10.3e\n",$eneout,$MEaa3Denv[0],$MZEaa3Denv;}
  if ($outscore eq "M-Elocal") {printf OUTSCR "%s\t%10.3e\t%10.3e\n",$eneout,$MElocal[0],$MZElocal;}
  if ($outscore eq "M-Ecmp") {printf OUTSCR "%s\t%10.3e\t%10.3e\n",$eneout,$MEcmp[0],$MZEcmp;}
  if ($outscore eq "M-E3Denv") {printf OUTSCR "%s\t%10.3e\t%10.3e\n",$eneout,$ME3Denv[0],$MZE3Denv;}
  if ($outscore eq "M-E3D") {printf OUTSCR "%s\t%10.3e\t%10.3e\n",$eneout,$ME3D[0],0.0;}



  printf "\n";
  printf "%10s\t%10.3e\t%10.3e\n","U-Score   ",$UPscore,$UZscore;
  printf "%10s\t%10.3e\t%10.3e\n","U-Energy  ",$UEnergy,$UZEnergy;
  printf "%10s\t%10.3e\t%10.3e\n","U-Dual    ",$UEnergy2,$UZEnergy2;
  printf "%10s\t%10.3e\t%10.3e\n","U-Eaa3Denv",$UEaa3Denv,$UZEaa3Denv;
  printf "%10s\t%10.3e\t%10.3e\n","U-Elocal  ",$Elocal[0],$ZElocal;
  printf "%10s\t%10.3e\t%10.3e\n","U-Ecmp    ",$Ecmp[0],$ZEcmp;
  printf "%10s\t%10.3e\t%10.3e\n","U-E3Denv  ",$UE3Denv,$UZE3Denv;
  printf "%10s\t%10.3e\t%10.3e\n","U-E3D     ",$E3D[0],0.0;
  $Score{Lower}=$UZscore;
  $Energy{Lower}=$UZEnergy;
  $Energy2{Lower}=$UZEnergy2;
  $Eaa3Denv{Lower}=$UZEaa3Denv;
  $Elocal{Lower}=$ZElocal;
  $Ecmp{Lower}=$ZEcmp;
  $E3Denv{Lower}=$UZE3Denv;
  if ($outscore eq "U-Score") {printf OUTSCR "%s\t%10.3e\t%10.3e\n",$eneout,$UPscore[0],$UZscore;}
  if ($outscore eq "U-Energy") {printf OUTSCR "%s\t%10.3e\t%10.3e\n",$eneout,$UEnergy[0],$UZEnergy;}
  if ($outscore eq "U-Dual") {printf OUTSCR "%s\t%10.3e\t%10.3e\n",$eneout,$UEnergy2[0],$UZEnergy2;}
  if ($outscore eq "U-Eaa3Denv") {printf OUTSCR "%s\t%10.3e\t%10.3e\n",$eneout,$UEaa3Denv[0],$UZEaa3Denv;}
  if ($outscore eq "U-Elocal") {printf OUTSCR "%s\t%10.3e\t%10.3e\n",$eneout,$UElocal[0],$UZElocal;}
  if ($outscore eq "U-Ecmp") {printf OUTSCR "%s\t%10.3e\t%10.3e\n",$eneout,$UEcmp[0],$UZEcmp;}
  if ($outscore eq "U-E3Denv") {printf OUTSCR "%s\t%10.3e\t%10.3e\n",$eneout,$UE3Denv[0],$UZE3Denv;}
  if ($outscore eq "U-E3D") {printf OUTSCR "%s\t%10.3e\t%10.3e\n",$eneout,$UE3D[0],0.0;}


 
  printf "\n";
  printf "%10s\t%10.3e\t%10.3e\n","L-Score   ",$LPscore,$LZscore;
  printf "%10s\t%10.3e\t%10.3e\n","L-Energy  ",$LEnergy,$LZEnergy;
  printf "%10s\t%10.3e\t%10.3e\n","L-Dual    ",$LEnergy2,$LZEnergy2;
  printf "%10s\t%10.3e\t%10.3e\n","L-Eaa3Denv",$LEaa3Denv,$LZEaa3Denv;
  printf "%10s\t%10.3e\t%10.3e\n","L-Elocal  ",$Elocal[0],$ZElocal;
  printf "%10s\t%10.3e\t%10.3e\n","L-Ecmp    ",$Ecmp[0],$ZEcmp;
  printf "%10s\t%10.3e\t%10.3e\n","L-E3Denv  ",$LE3Denv,$LZE3Denv;
  printf "%10s\t%10.3e\t%10.3e\n","L-E3D     ",$E3D[0],0.0;
  $Score{Upper}=$LZscore;
  $Energy{Upper}=$LZEnergy;
  $Energy2{Upper}=$LZEnergy2;
  $Eaa3Denv{Upper}=$LZEaa3Denv;
  $Elocal{Upper}=$ZElocal;
  $Ecmp{Upper}=$ZEcmp;
  $E3Denv{Upper}=$LZE3Denv;
  if ($outscore eq "L-Score") {printf OUTSCR "%s\t%10.3e\t%10.3e\n",$eneout,$LPscore[0],$LZscore;}
  if ($outscore eq "L-Energy") {printf OUTSCR "%s\t%10.3e\t%10.3e\n",$eneout,$LEnergy[0],$LZEnergy;}
  if ($outscore eq "L-Dual") {printf OUTSCR "%s\t%10.3e\t%10.3e\n",$eneout,$LEnergy2[0],$LZEnergy2;}
  if ($outscore eq "L-Eaa3Denv") {printf OUTSCR "%s\t%10.3e\t%10.3e\n",$eneout,$LEaa3Denv[0],$LZEaa3Denv;}
  if ($outscore eq "L-Elocal") {printf OUTSCR "%s\t%10.3e\t%10.3e\n",$eneout,$LElocal[0],$LZElocal;}
  if ($outscore eq "L-Ecmp") {printf OUTSCR "%s\t%10.3e\t%10.3e\n",$eneout,$LEcmp[0],$LZEcmp;}
  if ($outscore eq "L-E3Denv") {printf OUTSCR "%s\t%10.3e\t%10.3e\n",$eneout,$LE3Denv[0],$LZE3Denv;}
  if ($outscore eq "L-E3D") {printf OUTSCR "%s\t%10.3e\t%10.3e\n",$eneout,$LE3D[0],0.0;}


  printf "\n";
  printf "\t%10s"," ";
  foreach $key (keys %Score){printf "\t%10s",$key;}
  printf "\n";
  printf "\t%10s","Score";
  foreach $key (keys %Score){printf "\t%10.3e",$Score{$key};}
  printf "\n";
  printf "\t%10s","Energy";
  foreach $key (keys %Energy){printf "\t%10.3e",$Energy{$key};}
  printf "\n";
  printf "\t%10s","Eaa3Denv";
  foreach $key (keys %Eaa3Denv){printf "\t%10.3e",$Eaa3Denv{$key};}
  printf "\n";
  printf "\t%10s","Elocal";
  foreach $key (keys %Elocal){printf "\t%10.3e",$Elocal{$key};}
  printf "\n";
  printf "\t%10s","Ecmp";
  foreach $key (keys %Ecmp){printf "\t%10.3e",$Ecmp{$key};}
  printf "\n";
  printf "\t%10s","E3Denv";
  foreach $key (keys %E3Denv){printf "\t%10.3e",$E3Denv{$key};}
  printf "\n";
  printf "\t%10s","Conformer";
  foreach $key (keys %Energy2){printf "\t%10.3e",$Energy2{$key};}
  printf "\n";
  printf "\t%10s","ZEnergy";
  foreach $key (keys %E3Denv){if (defined $ZEnergy{$key}) {printf "\t%10.3e",$ZEnergy{$key};}else{printf "\t%10s","-";}}
  printf "\n";

  close (ENE);
  select(STDOUT);

  } #Loop on @PDB files
 
  if (defined $xml) { 
     printf  OUTXML "</xml>\n";
     close(OUTXML);
  }
  printf "[SUCCESS] Execution completed!\n";

sub Zscore
{
 local (%in)=@_;
 local ($AV,$AV2,$z,$ii);

  $AV=$in->{x};
  $AV2=$AV*$AV;
  for ($ii=0;$ii<$in->{n};$ii++){
     $AV  += $in->{y}->[$ii];
     $AV2 += $in->{y}->[$ii] * $in->{y}->[$ii];
  }
  $AV =$AV/(1+$in->{n});
  $AV2=$AV2/(1+$in->{n});
  if ( ($AV2 - $AV*$AV) >1.0e-10) {$z=( ($in->{x} - $AV)/sqrt($AV2-$AV*$AV) );}else{ $z=0.0;}
#  printf "%10.3e\t%10.3e\t%10.3e\t%10.3e\t%10.3e\n",$z,$in->{x},$AV, $AV2, sqrt($AV2-$AV*$AV) ;
  return $z;
}
  

sub Graph
{
  local (%in)=@_;
  local (@g,%g,$g,$min,$max,$step,$n,$re,$ra,$ene2,$ene1,$ene,$i,$av,$av2,$sigma,$nr);

  $min=$max=0;
  $av2=0;
  $av=0;
  $nr=0;
  for ($i=0;$i<$in->{dim};$i++){
         $found=0;
         for $n (@{$in->{real}}){if ($i == $n){$found=1;last;}}
         if ($found==0){
           $av2+=$in->{ene}->[$i]*$in->{ene}->[$i]; 
           $av +=$in->{ene}->[$i];
           $nr++;
         }
      }
  if ($nr>0){
   $av2=$av2/$nr;
   $av =$av/$nr;
  }
  if (($av2-$av*$av)>0){$sigma=sqrt($av2-$av*$av);}else{$sigma=0;}
  $max=$av+10*$sigma;
  $min=$av-10*$sigma;
  $step=($max-$min)/$in->{step};
  if ($step eq 0){
     $min=-1.0;
     $max= 1.0;
     $step=($max-$min)/$in->{step};
     }
  undef @g;
  for $j (0..$in->{step}){
    $re=0;
    $ra=0;
    $ene1=$min+$j*$step;
    $ene2=$min+($j+1)*$step;
    $ene =($ene1+$ene2)/2;
    for ($i=0;$i<$in->{dim};$i++){
         if ( $in->{ene}->[$i] < $ene2 && $in->{ene}->[$i] >= $ene1){
            $found=0;
            for $n (@{$in->{real}}){if ($i == $n){$found=1;last;}}
            if ($found==1){$re++;}else{$ra++;}
         }
    }
    $g={};
    $g={
        ene    => $ene,
        real   => $re,
        random => $ra
       };
    push @g,$g;
  }
  return @g;
}


sub ReadPDB
{
 local (%in)=@_;
 local (@prt,@fields,%atom,$atom,$pdb,$check,$chain,$type,$atmnum,$atomn,$res,$resnum,$coordx,$coordy,$coordz,$occ,$bfact,$protname);

 $pdb=$in->{file};
 $check=$in->{chain};

 print "Read PDB $pdb CHAIN $check \n";
 open (PDB,"<$pdb");
 while (<PDB>)
 {  
      @fields = split;
      if ($fields[0] eq "ATOM" || $fields[0] eq "HETATM"){
        $type=$fields[0];
	$type=~s/\s//g;
 	$atmnum=substr($_,6,5);
	$atmnum=~s/\s//g;
	$atomn=substr($_,11,6);
	$atomn=~s/\s//g;
	$res=substr($_,17,3);
	$res=~s/\s//g;
	$chain=substr($_,21,2);
	$chain=~s/\s//g;
	$resnum=substr($_,23,7);
	$resnum=~s/\s//g;
	$coordx=substr($_,29,9);
	$coordx=~s/\s//g;
	$coordy=substr($_,39,7);
	$coordy=~s/\s//g;
	$coordz=substr($_,47,7);
	$coordz=~s/\s//g;
        $occ=substr($_,55,5);
        $occ=~s/\s//g;
        $bfact=substr($_,60,5);
        $bfact=~s/\s//g;
        $protname=substr($_,66,10);
        $protname=~s/\s//g;
        if ($chain eq $check){
        for $j (0..$#code){if ($code[$j] eq $res){$aa=substr($lode,$j,1);}}
        $atom = {};
        $atom = {
           resnum => $resnum,
           residue=> $res,
           aa     => $aa,
           atom   => $atomn,
           atmnum => $atmnum,
           occup  => $occ,
           bfact  => $bfact,
           type   => $protname,
           x      => $coordx,
           y      => $coordy,
           z      => $coordz,
          };
        push @prt,$atom; 
        }
      }else{next;};

 }

 close (PDB);
 
 return @prt;

}

sub ReadDSSP
{
 local (%in)=@_;
 local ($dsspf,$begin,$aa,$ss,@ss,@bur,$chain,$b,$bb,$bur,%sol,%out);

 $sol={};
 $sol={  A => 115.0, R => 263.0, N => 184.0, D => 170.0, C => 149.0, 
        Q => 208.0, E => 207.0,
        G => 86.0, H => 206.0, I => 187.0, L => 192.0, K => 222.0, 
        M => 210.0, F => 230.0,
        P => 140.0 , S => 140.0, T => 164.0, W => 269.0, Y => 257.0, V => 161.0};


 $dsspf=$in->{file};
 $check=$in->{chain};

  print "Read DSSP $dsspf CHAIN $check \n";
  undef @ss;
  undef @bur;

  open (DSSP,"<$dsspf");
  $begin=0;
  while (<DSSP>){
    if (/  #  RESIDUE AA STRUCTURE BP1 BP2  ACC   N-H-->O  O-->H-N  N-H-->O  O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA/){$begin=1; next;}
    if ($begin==1){
     $aa=substr($_,13,1);
     if ($aa ne "!"){
      $ss=substr($_,16,1);
      $b=substr($_,35,3);
      $chain=substr($_,11,1);
      if (defined $sol->{$aa}){$bb=int(10*$b/$sol->{$aa});}else{$bb=0;}
      if ($bb < 0){$bb=0;}
      if ($bb >= 10) {$bb= 9;}
      if ($bb>5){$bur="E";}else{$bur="B";}
      if ($chain eq $check){
         if ($ss ne "H" && $ss ne "E"){ push @ss,"C"; }else{ push @ss,$ss; }
         push @bur,$bur;
         }
     }
    }
  }
 close (DSSP);
 
 $out={};
 $out={
        ss  => [@ss],
        bur => [@bur],
      };
 return %out;
}

