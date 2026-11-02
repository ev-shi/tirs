Commands
------------

BASELINE MODEL
# Parameter Estimation

Compile: 
  /usr/local/bin/gcc-15 -o execute/mcmc code-baseline/mainTIRS.c \
	code-baseline/functionsTIRS.c -O3 -lm -lgsl -lgslcblas \
	-I/usr/local/include -L/usr/local/lib -Wall
	
Run:
	./execute/mcmc -sm 10 -bn 0 \
	-ag input/c.global.mcmc.alpha.txt \
	-ai input/c.positives.alpha.txt \
  -ac input/c.incidence.txt
	
	Description:
		-sm: number of simulations
		-bn: burn-in simulations
		-ag: input global data
		-ai: input individual data
		-ac: daily contribution of community cases
		-ainf: input infection times
		-asc: input school characteristics
		-axns: daily index of attendance to the school (for the negative school 
		  contacts)
		-ninft: number of different infection times simulated per each positive case

# Pooling posteriors
Execute the R code in the "code-baseline" folder, for the selection of the 
posteriors parameters and to create the files needed for the reconstruction of 
the transmission chains.

# Reconstruction of transmission chains:
Compile:
  /usr/local/bin/gcc-15 -o execute/chain-s1 code-baseline/chainSchools.c \
  code-baseline/functionsSchools.c -O3 -lm -lgsl -lgslcblas \
  -I/usr/local/include -L/usr/local/lib -Wall

Run:
  ./execute/chain-s1 -sd 456 -ag input/c.global.mcmc.alpha.txt -ai \
  input/c.positives.alpha.txt -an input/c.negatives.alpha.txt -axf \
  input/index.attendance.household.alpha.txt -axs \
  input/index.attendance.school.alpha.txt -axc \
  input/index.attendance.community.alpha.txt -axcn \
  input/index.attendance.community.neg.alpha.txt -ao input/c.incidence.RE.txt \
  -ainf output/infection.times.alpha.2.txt -asc \
  input/c.schools.characteristics.alpha.txt -axns \
  input/index.attendance.school.hyponeg.alpha.txt -apar \
  output/posterior.params.alpha-s1.txt -ninft 5000 -np 5000
  
  Description:
  		-sd: seed
  		-ag: input global data
  		-ai: input individual data
  		-an: input negative household individuals' data
  		-axf: daily index of attendance to the household (for positive cases)
  		-axs: daily index of attendance to the school (for positive cases)
  		-axc: daily index of attendance to the community (for positive cases)
  		-axcn: daily index of attendance to the community (for the negative 
  		  household members)
  		-ao: daily contribution of community cases
  		-ainf: input infection times
  		-asc: input school characteristics
  		-axns: daily index of attendance to the school (for the negative school 
  		  contacts)
  		-apar: selection of posterior parameters
  		-ninft: = np
  		-np: number of posterior parameters selected # originally 50000