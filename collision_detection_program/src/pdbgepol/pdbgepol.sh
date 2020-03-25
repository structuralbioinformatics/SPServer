
PDBtoSplitChain.pl -i pdb4cpa.ent -o 4cpa
cat 4cpaA.pdb>4cpa_pci.pdb
cat 4cpaI.pdb >>4cpa_pci.pdb
echo "REMARK COMPLEX" >complex.pdb
perl /home/boliva/PROGRAMS/help/PERL/arrangeR.pl 4cpa_pci.pdb complex_0.pdb
grep ATOM complex_0.pdb|wc -l > natom
cat complex_0.pdb>>complex.pdb
\rm complex.inp
./pdbgepol <complex_pdbgepol.inp
\rm complex.esf complex.vec complex.dib  complex.out
\rm fort.7 fort.8 fort.15
ln -s complex.esf fort.7
ln -s complex.vec fort.8
ln -s complex.dib fort.15
./gepol <complex.inp>complex.out
\rm fort.7 fort.8 fort.15

echo "REMARK CPA">cpa.pdb
perl /home/boliva/PROGRAMS/help/PERL/arrangeR.pl 4cpaA.pdb cpa_0.pdb
grep ATOM cpa_0.pdb|wc -l > natom
cat cpa_0.pdb>>cpa.pdb
\rm cpa.inp
./pdbgepol <cpa_pdbgepol.inp
\rm cpa.esf cpa.dib cpa.vec cpa.out
\rm fort.7 fort.8 fort.15
ln -s cpa.esf fort.7
ln -s cpa.vec fort.8
ln -s cpa.dib fort.15
./gepol <cpa.inp>cpa.out
\rm fort.7 fort.8 fort.15

echo "REMARK PCI">pci.pdb
perl /home/boliva/PROGRAMS/help/PERL/arrangeR.pl 4cpaI.pdb pci_0.pdb
grep ATOM pci_0.pdb|wc -l > natom
cat pci_0.pdb>>pci.pdb
\rm pci.inp
./pdbgepol <pci_pdbgepol.inp
\rm pci.esf pci.dib pci.vec pci.out
\rm fort.7 fort.8 fort.15
ln -s pci.esf fort.7
ln -s pci.vec fort.8
ln -s pci.dib fort.15
./gepol <pci.inp>pci.out
\rm fort.7 fort.8 fort.15



