# NEEDS TO BE UPDATED!!!

# Re=20 cylinder (rotating)
Compare different magainst C++ implementation (lab and body-fixed frames).
Rotation is $\Omega = 0.1$ rad/s, which gives a theoretical prediction of $C_L = -0.157$ according to Kutta-Joukouski
There are three different approaches to treating the cylinder:
1. Keep the body fixed in the lab frame, but specify variable velocity.  The reg/interp operators don't need to be recomputed - only change is the RHS forcing term from body velocity
2. Simulate in the body-fixed frame with time-varying base flow (see Hsieh-Chen Tsai & Steve Brunton theses)
3. Full moving body in lab frame, recomputing operators and using TangentSE2 to update points.

In theory all three should give equivalent results, but this benchmark is on a pretty small domain, so there might be some differences due to boundary effects, etc.

__NOTE__: The "Runtime" benchmarks are to T=10 with dt=1e-2, not to full convergence (usually about T=100)

| Package      | Script | Runtime (secs) |  $C_L$ |
| ----------- | ----------- | ----------- | -----|
| C++ (lab frame)*  | N/A | 573   |  -0.139 |
| C++ (body-fixed)  | N/A |  78 |   |
| Julia (case 1)  | 'cyl_rot1.jl' |  57 | -0.200 |
| Julia (case 2)  |  'cyl_rot2.jl'   |  82 |  -0.172  |
| Julia (case 3)  |  'cyl_rot3.jl'   |  |  |

*Did not converge to steady state

__TO DO:__ FINISH THESE BENCHMARKS (and run in C++)

### Eldredge maneuver

Canonical pitch-hold maneuver (see folder for details)
