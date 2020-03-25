#!/usr/bin/perl
#


 use Getopt::Long;



 &GetOptions("ff=s"  => \$field,
             "w=s"   => \$distance_weight,
             "sw=s"  => \$group_weight,
             "fw=s"  => \$file_weight,
             "c=s"   => \$cutoff,
             "pot=s" => \$pot,
             "r=s"   => \$nran,
             "g=s"   => \$gap,
             "p=s"   => \$pdb,
             "a=s"   => \$aln,
             "re=s"  => \$select,
             "o=s"   => \$out,
             "h"     => \$help       );

 unless (defined $field  && not defined $help && defined $nran && defined $pdb && $cutoff >= 2.0){
   printf "
          ff   Force Field File with Potentials
          c    Cut-off for interacting pairs should be real, larger than 2 and smaller than 50
          pot  Potential can be by CA, CB or Minimum atom-atom distance: ca,cb,min.
               Default is min
          w    if w=0 Fold weights do not include distance dependency, 
               if w=1 Fold weights include distance dependency, default is 0
          fw   Force Field for Weights: File with precomputed environment weights
          sw   Interval window to group weights by distance-intervals
          r    Number of random decoys
          g    Number of residues ahead which interaction is skipped, default is 3 
          p    PDB of structure
	  re   Region of residues [r1-r2,r3-r4,...] Default is all
          a    Sequence alignment
          o    Output file
          \n";
   exit;
  }

  if (not defined $gap){$gap=3;}
  if (not defined $pot){$pot="min";}
  if (not defined $distance_weight){$distance_weight=0;}
  if (not defined $group_weight) {$group_weight=1;}
  if (not defined $select){
	  printf "Select the whole protein\n";
	  $select=0;
  }else{
	  @fragments=split /,/,$select;
          for $frag (@fragments){
		  @border= split /-/,$frag;
		  for ($i=$border[0];$i<=$border[1];$i++){push @region,$i;}
	  }
          printf "Selected region @region \n";
  }

  if ($cutoff>50.0){$cutoff=50.0;}


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

$sol={};
$sol={  A => 115.0, R => 263.0, N => 184.0, D => 170.0, C => 149.0, 
        Q => 208.0, E => 207.0,
        G => 86.0, H => 206.0, I => 187.0, L => 192.0, K => 222.0, 
        M => 210.0, F => 230.0,
        P => 140.0 , S => 140.0, T => 164.0, W => 269.0, Y => 257.0, V => 161.0};


 $dssp="./bin/dssp";
 $interact="./bin/fold";
 
# Read Potentials

 $skip=0;

 if (-e $field) {printf "Potential file $file found\n";}else{printf "$field not found\n";}

 if ($distance_weight == 0){

  open (INP,"<$field") || exit();  
   while(<INP>){
     if (/END/){$skip=0;}
     if (/^PAIR-POTENTIALS OF RESIDUES  \n/){$skip=1;next;}
     if (/^PAIR-POTENTIALS OF ENVIRONMENTS/){$skip=2;next;}
     if (/^LOCAL PAIR-POTENTIALS OF RESIDUES AND ENVIRONMENTS/){$skip=3;next;}
     if (/^PAIR-POTENTIALS OF RESIDUES AND ENVIRONMENTS  \n/){$skip=4;next;}
     if (/^STANDARD PAIR-POTENTIALS OF DISTANCES/){$skip=5;next;}
     if (/^SOLUTION FOR ENVIRONMENTAL WEIGHTS/){$skip=6;next;}
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
   }
   close INP;

   if (-e $file_weight) {printf "Weights file $file_weight found\n";}else{printf "Weights file $file_weight not found\n";}

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

  }else{

  open (INP,"<$field") || exit();  
   while(<INP>){
     if (/END/){$skip=0;}
     if (/^PAIR-POTENTIALS OF RESIDUES  \n/){$skip=1;next;}
     if (/^PAIR-POTENTIALS OF ENVIRONMENTS/){$skip=2;next;}
     if (/^LOCAL PAIR-POTENTIALS OF RESIDUES AND ENVIRONMENTS/){$skip=3;next;}
     if (/^PAIR-POTENTIALS OF RESIDUES AND ENVIRONMENTS  \n/){$skip=4;next;}
     if (/^STANDARD PAIR-POTENTIALS OF DISTANCES/){$skip=5;next;}
     if (/^SOLUTION FOR ENVIRONMENTAL WEIGHTS/){$skip=6;next;}
     if ($skip==1){ if (/PAIR/){chop;@data=split;$pair=$data[1]; for ($i=2;$i<=$#data;$i++){ $respmf{$pair}[$i-2]=$data[$i];}} }
     if ($skip==2){ if (/PAIR/){chop;@data=split;$pair=$data[1]; for ($i=2;$i<=$#data;$i++){ $epmf{$pair}[$i-2]=$data[$i];}} }
     if ($skip==3){ if (/AA-ENV/){chop;@data=split;$pair=$data[1]; $rpmf{$pair}=$data[2];} }
     if ($skip==4){ if (/PAIR/){chop;@data=split;$pair=$data[1]; for ($i=2;$i<=$#data;$i++){ $pmf{$pair}[$i-2]=$data[$i];}} }
     if ($skip==5){ if (/Distance/){chop;@data=split;$d=$data[1]; $dpmf[$d]=$data[2];} }
     if ($skip==6){ if (/Pair/){
                                chop;
                               ($pair,$factor)=split /=>/,substr($_,10,100);
                                $pair=~s/ //g;
                               ($env1,$env2,$d)=split /:/,$pair;
                                $pair_reverse="$env2".":"."$env1";
                                $weight{$pair}[$d]=$factor/(-1*$boltz*$T);
                                $weight{$pair_reverse}[$d]=$factor/(-1*$boltz*$T);
                               } }
   }
   close INP;

   if (-e $file_weight) {printf "Weights file $file_weight found\n";}else{printf "Weights file $file_weight not found\n";}

   open (INP,"<$file_weight");
    while(<INP>){
     if (/END/){$skip=0;}
     if (/^SOLUTION FOR ENVIRONMENTAL WEIGHTS/){$skip=6;next;}
     if ($skip==6){ if (/Pair/){
                                chop;
                               ($pair,$factor)=split /=>/,substr($_,10,100);
                                $pair=~s/ //g;
                               ($env1,$env2,$d)=split /:/,$pair;
                                $pair_reverse="$env2".":"."$env1";
                                $weight{$pair}[$d]=$factor/(-1*$boltz*$T);
                                $weight{$pair_reverse}[$d]=$factor/(-1*$boltz*$T);
                               } }
    }
   close INP;

  }


# Read PDB 

 $np=0; 

 open (PDB,"<$pdb");
 $dsspf="$pdb".".dssp";
 system("$dssp $pdb $dsspf ");
 undef @prt;
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
      }else{next;};

 }

 close (PDB);

  undef @ss;
  undef @bur;


 open (DSSP,"<$dsspf");
 print "Reading DSSP $dsspf \n";
  $begin=0;
  while (<DSSP>){
    if (/  #  RESIDUE AA STRUCTURE BP1 BP2  ACC   N-H-->O  O-->H-N  N-H-->O  O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA/){$begin=1; next;}
    if ($begin==1){
      $aa=substr($_,13,1);
      if ($aa eq "!") {next;}
      $ss=substr($_,16,1);
      $b=substr($_,35,3);
      $chain=substr($_,11,1);
      if (defined $sol->{$aa}){$bb=int(10*$b/$sol->{$aa});}else{$bb=0;}
      if ($bb < 0){$bb=0;}
      if ($bb >= 10) {$bb= 9;}
      if ($bb>5){$bur="E";}else{$bur="B";}
      if ($ss ne "H" && $ss ne "E"){ push @ss,"C"; }else{ push @ss,$ss; }
      push @bur,$bur;
    }
  }
 close (DSSP);


 #TRANSFORM INTO SEQUENCE

  $natom=$#prt+1;
  undef @seq;
  undef @pol;
  $nres=0;
  for $i (0..$#prt)
   {
     $res=  $prt[$i]->{residue};
     $letter=22;
     for $j (0..$#code){if ($code[$j] eq $res){$letter=$j;last;}}
     if ($prt[$i]->{atom} eq "CA" ){ push @seq,substr($lode,$letter,1);push @pol,substr($pode,$letter,1); $nres++; }
   }

  if ($natom <= 4*$nres || (1+$#ss) ne $nres){
   printf "Rejected ($#ss NE $nres) $pdb\n";
   next;
  }

 printf "SEQUENCE $#prt ($#seq) =>\t @seq\n";

# Read alignment 

  undef @query;
  undef @query1;
  undef @tmp;
  undef @qol;
  undef @tsa;
  undef %qsa;

  if (defined $aln && -e $aln) {
   open (ALN,"<$aln");
   undef @sequence;
   $skip=0;
   $seq="";
   while (<ALN>){
    if (/^>/ && $skip ==0){$skip=1;$seq="";next;}
    if (/^>/ && $skip ==1){$seq=~/\s//g;push @sequence,$seq;$seq="";next;}
    if ($skip == 1){chop;$seq.="$_";}
   }
   push @sequence,$seq;
   close ALN;
  }else{
   for $i (1,2){
   $seq="";
   for $aa (@seq){$seq.="$aa";}
   push @sequence,$seq;
   }
  }
  $found_template=0;
  for $seq (@sequence){
    undef @string;
    undef @word;
    @word=split /-/,$seq; 
    for $word (@word){ push @string,$word; } 
    $test    =join @string;
    $template=join @seq;
    if ($test eq $template && $found_template==0){
       $n=0;
       @word=split //,$seq;
       for $word (@word){ 
            push @tmp,$word;
            if ($word ne "-"){push @tsa,$n;}
            $n++;
        }
        $found_template=1;
    }else{
       $n=0;
       @word=split //,$seq;
       for $word (@word){ 
            push @query,$word;
            if ($word ne "-"){
               for($ka=0;$ka<22;$ka++){if ($word eq substr($lode,$ka,1)){$nk=$ka;last;}}
               $pol=substr($pode,$nk,1);
               push @qol,$pol;
               $qsa{$n}=1;
               push @query1,$word;
            }else{
               push @qol,"-";
               $qsa{$n}=0;
            }
            $n++;
        }
    }
  }


# CALCULATE DISTANCES

 undef @Eaa3Denv;
 undef @E3D;
 undef @E3Denv;
 undef @Elocal;
 undef @Pscore;
 undef @Ecmp;

 undef @aEaa3Denv;
 undef @aE3D;
 undef @aE3Denv;
 undef @aElocal;
 undef @aPscore;
 undef @aEcmp;


   $pdbDF="$pdb".".interact"."$pot";
   print "Running $interact with $pdb \n";
   $cutoff2= 1 + $cutoff;
   system("$interact -i $pdb -p $pot -c $cutoff2 -g $gap -o $pdbDF");
   printf "$interact -i $pdb -p $pot -c $cutoff2 -g $gap -o $pdbDF\n";

 open (DIST,"<$pdbDF") || exit();
   $skip=0;
   while(<DIST>){
     if ($skip==1){last;}
     ($n,$m,$dr,$aaa,$aab)=split;
     $i=$n-1;
     $j=$m-1;
     if ($select eq 0) {$use=1;}else{$use=0;}
     for $r (@region){
      if ($n eq $r || $m eq $r){$use=1;}
      if ($use eq 1) {last;}
     }
     if ($use eq 0 && $select ne 0) {next;}
     $testa=$seq[$i];
     $testb=$seq[$j];
     $seqa=$query[$tsa[$i]];
     $seqb=$query[$tsa[$j]];
     if ($seqa eq "-" || $seqb eq "-") {next;}
     $pola=$qol[$tsa[$i]];
     $polb=$qol[$tsa[$j]];
     $ss1 =$ss[$i];
     $ss2 =$ss[$j];
     $bur1=$bur[$i];
     $bur2=$bur[$j];
     $pair1= $seqa."-".$pola."-".$ss1."-".$bur1;
     $pair2= $seqb."-".$polb."-".$ss2."-".$bur2;
     $env1="$pola"."-"."$ss1"."-"."$bur1";
     $env2="$polb"."-"."$ss2"."-"."$bur2";
     if (not defined $seqa || not defined $seqb){
       printf "Error(1) (N,M)=($n,$m) pair in $pdb\n";
       $skip=1;
       next;
     }
     if ($seqa ne $aaa || $seqb ne $aab){
       printf "Error(1) (N,M)=($n,$m) pair($aaa,$aab) not ($seqa,$seqb)  in $pdb\n";
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
     $dw=$group_weight*(1+int($d/$group_weight));
     if ($distance_weight==1){ $aEcmp_res = $boltz*$T*$weight{$envpair}[$dw];}
     else                    { $aEcmp_res = $boltz*$T*$weight{$envpair};}
 
     $aPscore[$np][$i]  +=$respmf{$respair}[$d];
     $aEaa3Denv[$np][$i]+=$pmf{$pair}[$d];
     $aE3D[$np][$i]     +=$dpmf[$d];
     if     (defined $epmf{$envpair}[$d])         {$aE3Denv[$np][$i]  -=$epmf{$envpair}[$d];}
     elsif  (defined $epmf{$envpair_reverse}[$d]) {$aE3Denv[$np][$i]  -=$epmf{$envpair_reverse}[$d];}
     $aElocal[$np][$i]  -= $rpmf{$pair1};
     $aEcmp[$np][$i]    -= $aEcmp_res;
     if     (defined $epmf{$envpair}[$d]) {
              $aEnergy[$np][$i]  +=$pmf{$pair}[$d] + $dpmf[$d] - $rpmf{$pair1}  - $aEcmp_res - $epmf{$envpair}[$d];
     }else{
              $aEnergy[$np][$i]  +=$pmf{$pair}[$d] + $dpmf[$d] - $rpmf{$pair1}  - $aEcmp_res - $epmf{$envpair_reverse}[$d];
     }

     $aPscore[$np][$j]  +=$respmf{$respair}[$d];
     $aEaa3Denv[$np][$j]+=$pmf{$pair}[$d];
     $aE3D[$np][$j]     +=$dpmf[$d];
     if     (defined $epmf{$envpair}[$d])         {$aE3Denv[$np][$j]  -=$epmf{$envpair}[$d];}
     elsif  (defined $epmf{$envpair_reverse}[$d]) {$aE3Denv[$np][$j]  -=$epmf{$envpair_reverse}[$d];}
     $aElocal[$np][$j]  -= $rpmf{$pair2};
     $aEcmp[$np][$j]    -= $aEcmp_res;
     if     (defined $epmf{$envpair}[$d]) {
              $aEnergy[$np][$j]  +=$pmf{$pair}[$d] + $dpmf[$d] - $rpmf{$pair2} - $aEcmp_res - $epmf{$envpair}[$d];
     }else{
              $aEnergy[$np][$j]  +=$pmf{$pair}[$d] + $dpmf[$d] - $rpmf{$pair2} - $aEcmp_res - $epmf{$envpair_reverse}[$d];
     }

     $Pscore[$np]  +=$respmf{$respair}[$d];
     $Eaa3Denv[$np]+=$pmf{$pair}[$d];
     $E3D[$np]     +=$dpmf[$d];
     if     (defined $epmf{$envpair}[$d])         {$E3Denv[$np]  +=$epmf{$envpair}[$d];}
     elsif  (defined $epmf{$envpair_reverse}[$d]) {$E3Denv[$np]  +=$epmf{$envpair_reverse}[$d];}
     $Elocal[$np]  +=$rpmf{$pair1}+$rpmf{$pair2};
     $Ecmp[$np]    +=$aEcmp_res;
     $Energy[$np]   =$Eaa3Denv[$np] + $E3D[$np] - $Elocal[$np] - $E3Denv[$np] - $Ecmp[$np];


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
     undef @seqT;
     undef @seqR;
     undef @polR;
     @seqT=@query1;
     for $nn (0..$#query){
       if ($qsa{$n}==1){
        $nr=rand $#seqT;
        $aa=splice(@seqT,$nr,1);
        for($ka=0;$ka<22;$ka++){if ($aa eq substr($lode,$ka,1)){$nk=$ka;last;}}
        $pol=substr($pode,$nk,1);
       }else{
        $aa=$pol="-";
       }
       push @seqR,$aa;
       push @polR,$pol;
     }
     $Eaa3Denv[$np]=0;
     $E3D[$np]=0;
     $E3Denv[$np]=0;
     $Elocal[$np]=0;
     $Pscore[$np]=0;
     $Ecmp[$np]=0;
     $pdbDF="$pdb".".interact"."$pot";
     open (DIST,"<$pdbDF") || next;
       $skip=0;
       while(<DIST>){
         if ($skip==1){last;}
         ($n,$m,$dr,$aaa,$aab)=split;
         $i=$n-1;
         $j=$m-1;
         if ($select eq 0) {$use=1;}else{$use=0;}
         for $r (@region){
          if ($n eq $r || $m eq $r){$use=1;}
          if ($use eq 1) {last;}
          }
         if ($use eq 0 && $select ne 0) {next;}	 
         $testa=$seq[$i];
         $testb=$seq[$j];
         $seqa=$seqR[$tsa[$i]];
         $seqb=$seqR[$tsa[$j]];
         if ($seqa eq "-" || $seqb eq "-") {next;}
         $pola=$polR[$tsa[$i]];
         $polb=$polR[$tsa[$j]];
         $ss1 =$ss[$i];
         $ss2 =$ss[$j];
         $bur1=$bur[$i];
         $bur2=$bur[$j];
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
         $dw=$group_weight*(1+int($d/$group_weight));
         if ($distance_weight==1){ $aEcmp_res = $boltz*$T*$weight{$envpair}[$dw];}
         else                    { $aEcmp_res = $boltz*$T*$weight{$envpair};}
 
          $aPscore[$np][$i]  +=$respmf{$respair}[$d];
          $aEaa3Denv[$np][$i]+=$pmf{$pair}[$d];
          $aE3D[$np][$i]     +=$dpmf[$d];
          if     (defined $epmf{$envpair}[$d])         {$aE3Denv[$np][$i]  -=$epmf{$envpair}[$d];}
          elsif  (defined $epmf{$envpair_reverse}[$d]) {$aE3Denv[$np][$i]  -=$epmf{$envpair_reverse}[$d];}
          $aElocal[$np][$i]  -=$rpmf{$pair1};
          $aEcmp[$np][$i]    -=$aEcmp_res;
          if     (defined $epmf{$envpair}[$d]) {
              $aEnergy[$np][$i]  +=$pmf{$pair}[$d] + $dpmf[$d] - $rpmf{$pair1}  - $aEcmp_res  - $epmf{$envpair}[$d];
          }else{
              $aEnergy[$np][$i]  +=$pmf{$pair}[$d] + $dpmf[$d] - $rpmf{$pair1}  - $aEcmp_res  - $epmf{$envpair_reverse}[$d];
          }

          $aPscore[$np][$j]  +=$respmf{$respair}[$d];
          $aEaa3Denv[$np][$j]+=$pmf{$pair}[$d];
          $aE3D[$np][$j]     +=$dpmf[$d];
          if     (defined $epmf{$envpair}[$d])         {$aE3Denv[$np][$j]  -=$epmf{$envpair}[$d];}
          elsif  (defined $epmf{$envpair_reverse}[$d]) {$aE3Denv[$np][$j]  -=$epmf{$envpair_reverse}[$d];}
          $aElocal[$np][$j]  -=$rpmf{$pair2};
          $aEcmp[$np][$j]    -=$aEcmp_res;
          if     (defined $epmf{$envpair}[$d]) {
              $aEnergy[$np][$j]  +=$pmf{$pair}[$d] + $dpmf[$d]  - $rpmf{$pair2} - $aEcmp_res - $epmf{$envpair}[$d];
          }else{
              $aEnergy[$np][$j]  +=$pmf{$pair}[$d] + $dpmf[$d]  - $rpmf{$pair2} - $aEcmp_res - $epmf{$envpair_reverse}[$d];
          }


          $Pscore[$np]  +=$respmf{$respair}[$d];
          $Eaa3Denv[$np]+=$pmf{$pair}[$d];
          $E3D[$np]     +=$dpmf[$d];
          if     (defined $epmf{$envpair}[$d])         {$E3Denv[$np]  +=$epmf{$envpair}[$d];}
          elsif  (defined $epmf{$envpair_reverse}[$d]) {$E3Denv[$np]  +=$epmf{$envpair_reverse}[$d];}
          $Elocal[$np]  +=$rpmf{$pair1}+$rpmf{$pair2};
          $Ecmp[$np]    +=$aEcmp_res;
          if     (defined $epmf{$envpair}[$d]) {
              $Energy[$np]  +=$pmf{$pair}[$d] + $dpmf[$d] - $rpmf{$pair1} - $rpmf{$pair2} - $aEcmp_res  - $epmf{$envpair}[$d];
          }else{
              $Energy[$np]  +=$pmf{$pair}[$d] + $dpmf[$d] - $rpmf{$pair1} - $rpmf{$pair2} - $aEcmp_res  - $epmf{$envpair_reverse}[$d];
          }


        }
     close DIST;
     $np++;

   } #En random test

  printf "\n";
# Calculate Z-scores

  undef @aZEnergy;
  undef @aZEaa3Denv;
  undef @aZElocal;
  undef @aZE3Denv;
  undef @aZEcmp;
  undef @aZscore;
  undef @aZZEnergy;
  $ZEnergy=0;
  $Zaa3Denv=0;
  $ZElocal=0;
  $ZE3Denv=0;
  $ZEcmp=0;
  $Zscore=0;
  for ($i=0;$i<=$#seq;$i++){
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
      $aZZEnergy[$i]=$aZEaa3Denv[$i]  + $aZElocal[$i]  + $aZE3Denv[$i];
   }


  if ( $out ) { open (OUTPUT,">$out"); select(OUTPUT)}

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

 $real[0]=0;
  undef @GZ;
  $gstep=100;
  $nlist=$nran+1;
  printf "\nGRAPH FOR ENERGY-ZSCORE\n"; 
  printf "\t%10s\t%10s\t%10s\n","ENERGY","RANDOM","REAL";
  $in={};
  $in={
       step=> $gstep,
       dim => $nlist,
       real=>[@real],
       ene =>[@y]
      };
  @GZ=&Graph(%in);
  for $i (0..$#GZ){
       printf "\t%10.3e\t%10.1e\t%10.1e\n",$GZ[$i]->{ene},$GZ[$i]->{random},$GZ[$i]->{real};
      }
  printf "END\n";


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

  $real[0]=0;
  undef @GZ;
  $gstep=100;
  $nlist=$nran+1;
  printf "\nGRAPH FOR ZEaa3DEnv\n"; 
  printf "\t%10s\t%10s\t%10s\n","ENERGY","RANDOM","REAL";
  $in={};
  $in={
       step=> $gstep,
       dim => $nlist,
       real=>[@real],
       ene =>[@y]
      };
  @GZ=&Graph(%in);
  for $i (0..$#GZ){
       printf "\t%10.3e\t%10.1e\t%10.1e\n",$GZ[$i]->{ene},$GZ[$i]->{random},$GZ[$i]->{real};
      }
  printf "END\n";
                 

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


  $real[0]=0;
  undef @GZ;
  $gstep=100;
  $nlist=$nran+1;
  printf "\nGRAPH FOR ZElocal\n"; 
  printf "\t%10s\t%10s\t%10s\n","ENERGY","RANDOM","REAL";
  $in={};
  $in={
       step=> $gstep,
       dim => $nlist,
       real=>[@real],
       ene =>[@y]
      };
  @GZ=&Graph(%in);
  for $i (0..$#GZ){
       printf "\t%10.3e\t%10.1e\t%10.1e\n",$GZ[$i]->{ene},$GZ[$i]->{random},$GZ[$i]->{real};
      }
  printf "END\n";



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

 $real[0]=0;
  undef @GZ;
  $gstep=100;
  $nlist=$nran+1;
  printf "\nGRAPH FOR ZE3DEnv\n"; 
  printf "\t%10s\t%10s\t%10s\n","ENERGY","RANDOM","REAL";
  $in={};
  $in={
       step=> $gstep,
       dim => $nlist,
       real=>[@real],
       ene =>[@y]
      };
  @GZ=&Graph(%in);
  for $i (0..$#GZ){
       printf "\t%10.3e\t%10.1e\t%10.1e\n",$GZ[$i]->{ene},$GZ[$i]->{random},$GZ[$i]->{real};
      }
  printf "END\n";
                 


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
  $real[0]=0;
  undef @GZ;
  $gstep=100;
  $nlist=$nran+1;
  printf "\nGRAPH FOR ZEcmp\n"; 
  printf "\t%10s\t%10s\t%10s\n","ENERGY","RANDOM","REAL";
  $in={};
  $in={
       step=> $gstep,
       dim => $nlist,
       real=>[@real],
       ene =>[@y]
      };
  @GZ=&Graph(%in);
  for $i (0..$#GZ){
       printf "\t%10.3e\t%10.1e\t%10.1e\n",$GZ[$i]->{ene},$GZ[$i]->{random},$GZ[$i]->{real};
      }
  printf "END\n";


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

  $real[0]=0;
  undef @GZ;
  $gstep=100;
  $nlist=$nran+1;
  printf "\nGRAPH FOR CLASSIC PAIR-POTENTIAL Z-SCORE\n"; 
  printf "\t%10s\t%10s\t%10s\n","ENERGY","RANDOM","REAL";
  $in={};
  $in={
       step=> $gstep,
       dim => $nlist,
       real=>[@real],
       ene =>[@y]
      };
  @GZ=&Graph(%in);
  for $i (0..$#GZ){
       printf "\t%10.3e\t%10.1e\t%10.1e\n",$GZ[$i]->{ene},$GZ[$i]->{random},$GZ[$i]->{real};
      }
  printf "END\n";
       

 

   printf "ENERGIES PER RESIDUE\n"; 
   printf "%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n","number","sequence","Energy","P-Score","Eaa3Denv","E3D","-E3Denv","-Elocal","Ecmp";
   for ($i=0;$i<=$#seq;$i++){
    printf "%10d\t%10s\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\n",$i,$query[$tsa[$i]],$aEnergy[0][$i],$aPscore[0][$i],$aEaa3Denv[0][$i],$aE3D[0][$i],$aE3Denv[0][$i],$aElocal[0][$i],$aEcmp[0][$i]; 
   }
   printf "END\n";

   printf "\nZ-SCORES PER RESIDUE\n"; 
   printf "%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n","number","sequence","Energy","P-Score","Eaa3Denv","E3D","-E3Denv","-Elocal","Ecmp","ZEne_$pot";
   for ($i=0;$i<=$#seq;$i++){
    printf "%10d\t%10s\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.1e\n",$i,$query[$tsa[$i]],$aZEnergy[$i],$aZscore[$i],$aZEaa3Denv[$i],0.0,$aZE3Denv[$i],$aZElocal[$i],$aZEcmp[$i],$aZZEnergy[$i]; 
   }
   printf "END\n";

  if ( $out ) { close (OUTPUT);select(STDOUT); }

  printf "%10s\t%10.3e\t%10.3e\n","Score   ",$Pscore[0],$Zscore;
  printf "%10s\t%10.3e\t%10.3e\n","Energy  ",$Energy[0],$ZEnergy;
  printf "%10s\t%10.3e\t%10.3e\n","ZEne_$pot ",$Eaa3Denv[0]-$Elocal[0]-$E3Denv[0],$ZEaa3Denv-$ZElocal-$ZE3Denv;
  printf "%10s\t%10.3e\t%10.3e\n","Eaa3Denv",$Eaa3Denv[0],$ZEaa3Denv;
  printf "%10s\t%10.3e\t%10.3e\n","Elocal  ",$Elocal[0],$ZElocal;
  printf "%10s\t%10.3e\t%10.3e\n","Ecmp    ",$Ecmp[0],$ZEcmp;
  printf "%10s\t%10.3e\t%10.3e\n","E3Denv  ",$E3Denv[0],$ZE3Denv;
  printf "%10s\t%10.3e\t%10.3e\n","E3D     ",$E3D[0],0.0;


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


