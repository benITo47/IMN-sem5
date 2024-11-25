
set terminal pngcairo size 1600, 1200
set output "graphs/combined.png"

set palette defined (-10 "blue", 0 "white", 10 "red")
set size ratio -1
set xlabel "x"
set ylabel "y"
set pm3d map

set multiplot layout 2, 3 title "Poisson Equation Results"

set title "nx=ny=50"
splot [0:5][0:5] "data/n_50.dat" i 0 u 2:1:3

set title "nx=ny=100"
splot [0:10][0:10] "data/n_100.dat" i 0 u 2:1:3

set title "nx=ny=200"
splot [0:20][0:20] "data/n_200.dat" i 0 u 2:1:3

set title "nx=ny=100, e1=1, e2=1"
set cbrange [-0.8:0.8]
splot [0:10][0:10] "data/e_1_1.dat" i 0 u 2:1:3

set title "nx=ny=100, e1=1, e2=2"
set cbrange [-0.8:0.8]
splot [0:10][0:10] "data/e_1_2.dat" i 0 u 2:1:3

set title "nx=ny=100, e1=1, e2=10"
set cbrange [-0.8:0.8]
splot [0:10][0:10] "data/e_1_10.dat" i 0 u 2:1:3

unset multiplot


