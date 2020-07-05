./clean.sh

sed -i '1s/.*/0/' input.dat

./MolDyn_NVE.exe

sed -i '1s/.*/1/' input.dat

for i in {1..10} 
do
	./MolDyn_NVE.exe
done




