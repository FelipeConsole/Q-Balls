# Q-Balls

Q-balls are an example of solitons. More specifically, non-topological solitons.

In this repo I make available the code to solve the Q-balls equations described by Sidney Coleman in [here](https://inspirehep.net/literature/214529)

This repo contains:
1. The Mathematica notebooks to obtain the equations and translate them to Fortran form,
2. The fortran solver (colsys.f) and the program (bolasq.f90) that calls it to solve the equations.
3. A bash script (script.sh) to automate the increment in the parameters values, compile and run the programs.

There are some points I should clarify. For example, why use Fortran to solve the equations?
1. Well, obviously, Fortran is fast and efficient;
2. Non-linear differential equations are complicated and I was unable to use the standard techniques in Mathematica or Julia to solve them;
3. The most robust method I found is to use this solver, called [Colsys](https://dl.acm.org/doi/10.1145/355945.355951) that is still very used in the Gravitational Physics community.

Colsys uses an adaptive grid selection and a Newton method to solve the boundary value problems. Therefore, an initial guess solution is necessary (solu_init_q_ball_ipr_415_xl_1200.dat) that should be copied to solu.new before running bolasq.f90 with colsys.f.

It would be amazing if this method were implemented in Julia.
