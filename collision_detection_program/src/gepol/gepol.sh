
rm fort.7 fort.8 fort.15

ln -s example.esf fort.7
ln -s example.vec fort.8
ln -s example.dib fort.15

./gepol < exam6.inp > example.out

rm fort.7 fort.8 fort.15
