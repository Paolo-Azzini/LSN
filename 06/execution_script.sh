./clean.sh
rm graph.*

sed -i '3s/.*/50/' input.dat #Number of sites = 50
sed -i '4s/.*/1.0/' input.dat #Exchange energy J=1.0
sed -i '5s/.*/0.0/' input.dat #External field h=0.0
sed -i '7s/.*/20/' input.dat #Number of blocks = 20
sed -i '8s/.*/10000/' input.dat #Steps each bolck = 10000


for T in {0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0}
do
	sed -i '2s/.*/'$T'/' input.dat   #Setting temperature

	sed -i '1s/.*/0/' input.dat  
	./Monte_Carlo_ISING_1D.exe       #Equilibration run
	sed -i '1s/.*/1/' input.dat  

	sed -i '6s/.*/0/' input.dat      #Setting Gibbs mode
	./Monte_Carlo_ISING_1D.exe       


	sed -i '6s/.*/1/' input.dat      #Setting Metropolis mode
	./Monte_Carlo_ISING_1D.exe      
done

tail -q -n 1 output.ene.Gibbs* > graph.ene.Gibbs.dat
tail -q -n 1 output.heat.Gibbs* > graph.heat.Gibbs.dat
tail -q -n 1 output.chi.Gibbs* > graph.chi.Gibbs.dat

tail -q -n 1 output.ene.Metropolis* > graph.ene.Metropolis.dat
tail -q -n 1 output.heat.Metropolis* > graph.heat.Metropolis.dat
tail -q -n 1 output.chi.Metropolis* > graph.chi.Metropolis.dat

h=0.02
sed -i '5s/.*/'$h'/' input.dat #Setting external field h=0.02

for T in {0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0}
do
	sed -i '2s/.*/'$T'/' input.dat   #Setting temperature

	sed -i '1s/.*/0/' input.dat  
	./Monte_Carlo_ISING_1D.exe       #Equilibration run
	sed -i '1s/.*/1/' input.dat  

	sed -i '6s/.*/0/' input.dat      #Setting Gibbs mode
	./Monte_Carlo_ISING_1D.exe       


	sed -i '6s/.*/1/' input.dat      #Setting Metropolis mode
	./Monte_Carlo_ISING_1D.exe      
done

tail -q -n 1 output.mag.Gibbs* > graph.mag.Gibbs.dat

tail -q -n 1 output.mag.Metropolis* > graph.mag.Metropolis.dat

./clean.sh


