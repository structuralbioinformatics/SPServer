#!/usr/bin/perl
#

 use Getopt::Long;

 
 
 &GetOptions("i=s" => \$reference,
             "c=s" => \$addchain,
             "o=s" => \$crd       );


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

$lode="ACDEFGHIKLMNPQRSTVWY";


 
 open (REF,"<$reference")|| die ("Can't open input file:  $reference \n");

 if (length($crd) < 1) {$crd="$reference"."_split";printf "OUTPUT IN $crd \n";}

 $title=0;
 $NMR=0; 
 while (<REF>)
 {
 	@fields = split;
      if ($fields[0] eq "MODEL"){
         $chain=$fields[1];
         $NMR =1;
      }
      if ($fields[0] eq "ATOM"){
        $type=$fields[0];
	$type=~s/\s//g;
 	$atmnum=substr($_,6,5);
	$atmnum=~s/\s//g;
	$atomn=substr($_,12,4);
	$atomn=~s/\s//g;
        if ($atomn eq "OT1"){$atomn="O";}
        if ($atomn eq "OT2"){$atomn="OXT";}
	$res=substr($_,17,3);
	$res=~s/\s//g;
        if ($NMR == 0 ){
	$chain=substr($_,21,1);
	$chain=~s/\s//g;
        }
	$resnum=substr($_,22,7);
	$resnum=~s/\s//g;
	$coordx=substr($_,29,9);
	$coordx=~s/\s//g;
	$coordy=substr($_,39,7);
	$coordy=~s/\s//g;
	$coordz=substr($_,47,7);
	$coordz=~s/\s//g;
        $occ=substr($_,55,5);
        $occ=~s/\s//g;
        if ($occ eq ""){$occ=1.0;}
        $bfact=substr($_,60,5);
        $bfact=~s/\s//g;
        if ($bfact eq ""){$bfact=10.0;}
        $chain2=substr($_,72,2);
        $chain2=~s/\s//g;
        if ($chain eq "" && $chain2 ne ""){$chain=$chain2;}
        $atom = {};
        $atom = {
           resnum => $resnum,
           residue=> $res,
           chain  => $chain,
           atom   => $atomn,
           atmnum => $atmnum,
           occ    => $occ,
           bfact    => $bfact,
           x      => $coordx,
           y      => $coordy,
           z      => $coordz,
          };
        push @prot,$atom;
        if ($res eq "HOH" || $res eq "WAT" || $res eq "H2O"  ) {break;}
      }else{next;};

 }
 close (REF);

   $chains=" ";
   for $i (0..$#prot)
   { 
     if ( $prot[$i]->{chain} ne $chains) {$chains=$prot[$i]->{chain};push @sequences,$prot[$i]->{chain};}
   }
   for $i (0..$#sequences){
    push @out,"$crd"."$sequences[$i]".".fa";
    push @pdb,"$crd"."$sequences[$i]".".pdb";

   }

 for $k (0..$#out )
 {
   open (OUTPUT,">$out[$k]");;
   printf OUTPUT ">%s%s \n",$crd, $sequences[$k];
   $n=0;
   for $i (0..$#prot)
   {
     $res=  $prot[$i]->{residue};
     $letter="X";
     for $j (0..$#code){if ($code[$j] eq $res){$letter=substr($lode,$j,1);}}
     if ($prot[$i]->{atom} eq "CA" && $prot[$i]->{chain} eq $sequences[$k]){
         printf OUTPUT "%s",$letter;
         $n++;
         if ($n==80){printf OUTPUT "\n";$n=0;}
     }
   }
   printf OUTPUT "\n";
   close (OUTPUT);
   $outflush=select(PDB);
   $~ = "PDB_FORMAT";
   select($outflush);
   open (PDB,">$pdb[$k]");
   print PDB "REMARK TITLE\n";
   print PDB "REMARK TRANFORM PDB $reference TO $pdb\n";
   print PDB "REMARK TOTAL NUMBER OF ATOMS $#prot \n";
   $n=1;
   $resnumw=0;
   $resnum_0=" ";
   for $i (0..$#prot)
   {
     if ($prot[$i]->{chain} eq $sequences[$k]){
       $typew  ="ATOM";
       $chainw=$prot[$i]->{chain};
       if ($chainw eq "") {$chainw=$addchain;}
       $resnum = $prot[$i]->{resnum};
       if ($resnum ne $resnum_0){$resnumw++;$resnum_0=$resnum;} 
       $resw   = $prot[$i]->{residue};
       $atomw  = $prot[$i]->{atom};
       $atnumw = $n;
       @num=split /\./,$prot[$i]->{x};
       $x1 = $num[0];
       $x2 = substr($num[1],0,9);
       @num=split /\./,$prot[$i]->{y};
       $y1 = $num[0];
       $y2 = substr($num[1],0,9);
       @num=split /\./,$prot[$i]->{z};
       $z1 = $num[0];
       $z2 = substr($num[1],0,9);
       @num=split /\./,$prot[$i]->{occ};
       $o1 = $num[0];
       $o2 = substr($num[1],0,9);
       @num=split /\./,$prot[$i]->{bfact};
       $b1 = $num[0];
       $b2 = substr($num[1],0,9);
       write PDB;
       $n++;
     }else{$n=1;$resnum_0=" ";}
   }
   close PDB;
 }



format PDB_FORMAT=
@<<<<<<@>>>  @<<<@<< @@>>>    @>>>.@<<@>>>.@<<@>>>.@<<
$typew,$atnumw,$atomw,$resw,$chainw,$resnumw,$x1,$x2,$y1,$y2,$z1,$z2
.

  

