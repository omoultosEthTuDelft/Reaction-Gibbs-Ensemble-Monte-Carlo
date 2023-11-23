This code can be used to perform Monte Carlo (MC) simulations in the combined Reaction and Gibbs Ensemble. In these simulations, only Lennard-Jones interactions are considered. The directory "Source" contains the source code of the software written in Python. This code can be compiled using gfortran or Intel fortran compilers. The directory "Run" contains three example simulations. The subdirectories 'run_N2O4', 'run_NO2', and 'run_react' can be used to simulate phase and chemical equilibrium of pure N2O4 (at T=260K), pure NO2 (at T=240K) and the reactive mixture of NO2 and N2O4 (at T=270K). In the input files, the keywords represent:

Ncycle: total number of MC cycles
Ninit: number of initialization cycles
Linit: A logical that determines if the simulations start from a randomly generated configuration (.true.) or a configuration from a previous simulation (.false.). If started with a configuration from a previous simulation, a file called 'Coordold' needs to be present in the same directory (you can use the file 'Coordnew' generated at the end of each simulation).
Temp: Temperature
Pdisp: probability of choosing displacement trial move
Pswap: probability of choosing swap trial move
Pvol: probability of choosing volume trial move
Preact: probability of choosing reaction trial move
Deltax: Maximum displacement
Deltav: Maximum volume change
Lnpt: If the simulation is conducted in NPT ensemble (.true.) or NVT ensemble (.false.). If NPT ensemble is chosen, the software only simulates one simulation box in NPT ensemble.
Press: Pressure for NPT ensemble simulations
Box(1) and Box(2): Box sizes for the simulation boxes 1 and 2.
Ncomp: Number of components
Number of Molecs for Each Box: The lines after this line shows the number of molecules for each component in each box. If two simulation boxes are simulated, there must be two lines, each showing the number of molecules of each component.
Epsilon Sigma ln[Q] for each component: Three lines after this one show epsilon, sigma, and the natural logarithm of isolated molecule partition function for each component.
Nreact: Number of reactions. Stoichiometry of each reactions should be given after this line.
For more details, please refer to "Scaling Towards the Critical Point in the Combined Reaction/Gibbs Ensemble" from H. Mert Polat, Silvia Lasala, Frederick de Meyer, Celine Houriez, Othonas A. Moultos, and Thijs J. H. Vlugt.

