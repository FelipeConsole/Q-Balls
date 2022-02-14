 #! /bin/csh

  # This is a simple shell script that will (hopefully) cd into a dir that
  # contains a simple fortran routine that increments an integer and passes
  # th VAR=${VAR/,/.} is is as a job to the GRID engine.

  # Be sure to cd into the right dir (not sure if needed)
  cd Desktop/Quali/Q-balls/fortran_files

  # compile the fortran code then execute it with the ./ option

for cf in $(seq 1.1 -0.002 0.00)
do
echo ${cf/,/.} > parametros.dat

cp solu_init_q_ball_ipr_415_xl_1200.dat solu.new

#gfortran bolasq.f90 colsys.f -std=legacy -w
ifort bolasq.f90 colsys.f -w
 
  ./a.out
  
  done
