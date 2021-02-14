# Report on MPI parallelization of LJMD code


by Giulio Dondi

 
The Verlet molecular dynamics simulation was parallelised with MPI in thee different versions:

- **simple**  , which uses naive force calculation and splits the work equally among processes. Master broadcasts the atoms' positions before force calculations, calculates forces and potential energy in the naive way and later on reduces the forces vectors to the MAster processes, who is tasked with carrying out the integration steps and calculation of kinetic energy and temperature
- **full** ,  which parallelises the integration and kinetic energy calculations as well but still uses the naive force calculation. The positions are no longer merely broadcast from masterbut the collective MPI function Allgatherv is used, since at the beginning of any step each process is privy only to the updated positions of the particles that belong to him. Force reduction is skipped but, instead, the kinetic energy and temperature are reduced to master
- **simple + 3rd law** which is identical to the simple version but implements Newton's third law to halve the number of force calculations overall. To achieve proper processor balancing, dividing the array indices in P parts (P is the number of processes) is no longer adequate. A more sopisticated system represents the pairings of particles i-j as the upper triangle of a matrix NxN (N is the number of atoms) where each cell of this matrix represents an individual force calculation to perform. THis triangle is sliced in P portions of approximately equal area in order to try and minimise the difference between the largest and smallest area, since this determines the processor unbalancing.

Benchmark runs were conducted on the SISSA Ulysses cluster. One node with 2 sockets and 10 cores per socket was used. Process affinity was set as --map-by socket in order to distribute processes round-robin on both sockets and avoid a significant performance difference between P<=10 and P>10 runs .
System sizes of 108,2916 and 78732 were considered, while P varied as 1,2,4,8,12,16,20. Three runs per ceach case were performed and the results averaged.



![Total runtimes for each case](https://github.com/pogmat/ljmd/blob/master/reports/plots/MPI_runtimes.jpg?raw=true)

![Strong scalability curves, for different system sizes](https://github.com/pogmat/ljmd/blob/master/reports/plots/MPI_strong_scaling.jpg?raw=true)

![Parallel efficiency versus number of processes](https://github.com/pogmat/ljmd/blob/master/reports/plots/MPI_efficiency.jpg?raw=true)

![Detailed breakdown of the runtime for N=108](https://github.com/pogmat/ljmd/blob/master/reports/plots/MPI_108_times.jpg?raw=true)
![Detailed breakdown of the runtime for N=108](https://github.com/pogmat/ljmd/blob/master/reports/plots/MPI_2916_times.jpg?raw=true)
![Detailed breakdown of the runtime for N=108](https://github.com/pogmat/ljmd/blob/master/reports/plots/MPI_78732_times.jpg?raw=true)


Thee plots make evident several things about the problem and implementation:

-Scalability is strongly dependent with problem size, for very large problems the parallel overhead is much lessened, although still noticeable. 
- In the largest case, the efficiency bottoms off at 80% with about 12 or more processors, although runtimes still decrease slightly beyond this. Smaller problems highlight the parallel overhead and the efficiency takes a large hit.
- The **3rd law** version seems to have  marginally worse scalability, but comparing the runtimes it is still at least 2x faster than th other two as expected.
- In the breakdown plots, the blue represents time spend distributing the porisions around and calculating forces. This is by far the largest contribution to the total runtime, this explains the benefit of implementing Newton's third law.
- Although the other contributions are minor, the _ekin_ calculation time is noticeably large in the **full** implementation. This apparent delay is most likely not caused by the actual computations of the force, but rather by the MPI communication barrier right before the energy-temperature reduction step. This barrier evidences the processes imbalance during the various steps: as quicker processes need to wait for others the average time measured in this portion of the loop increases.

Based on these benchmarks, the **3rd law** version was chosen to be incorporated in the group project
