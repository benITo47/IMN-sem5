

set terminal pngcairo enhanced size 800,600

system("mkdir -p graphs")

set output 'graphs/energy_combined.png'
set title "Energy vs Time"
set xlabel "Time (t)"
set ylabel "Energy (E)"
set grid
plot \
    "data/E_alpha_0.0_beta_0.0.txt" using 1:2 with lines title "Alpha=0.0, Beta=0.0", \
    "data/E_alpha_0.0_beta_0.1.txt" using 1:2 with lines title "Alpha=0.0, Beta=0.1", \
    "data/E_alpha_0.0_beta_1.0.txt" using 1:2 with lines title "Alpha=0.0, Beta=1.0"

set output 'graphs/energy_alpha_1.0_beta_1.0.png'
set title "Energy vs Time (Alpha=1.0, Beta=1.0)"
plot "data/E_alpha_1.0_beta_1.0.txt" using 1:2 with lines title "Alpha=1.0, Beta=1.0"




set terminal pngcairo enhanced size 1200,600
set xlabel "Time (t)"
set ylabel "Position (x)"
set grid
set pm3d map
set palette defined (0 "navy", 1 "blue", 2 "cyan", 3 "green", 4 "yellow", 5 "red")
unset key  

stats "data/u_alpha_0.0_beta_0.0.txt" using 3 nooutput
set cbrange [STATS_min:STATS_max]

set output 'graphs/displacement_alpha_0.0_beta_0.0.png'
set title "Displacement Map (Alpha=0.0, Beta=0.0)"
splot "data/u_alpha_0.0_beta_0.0.txt" using 1:2:3 with pm3d

stats "data/u_alpha_0.0_beta_0.1.txt" using 3 nooutput
set cbrange [STATS_min:STATS_max]

set output 'graphs/displacement_alpha_0.0_beta_0.1.png'
set title "Displacement Map (Alpha=0.0, Beta=0.1)"
splot "data/u_alpha_0.0_beta_0.1.txt" using 1:2:3 with pm3d

stats "data/u_alpha_0.0_beta_1.0.txt" using 3 nooutput
set cbrange [STATS_min:STATS_max]

set output 'graphs/displacement_alpha_0.0_beta_1.0.png'
set title "Displacement Map (Alpha=0.0, Beta=1.0)"
splot "data/u_alpha_0.0_beta_1.0.txt" using 1:2:3 with pm3d

stats "data/u_alpha_1.0_beta_1.0.txt" using 3 nooutput
set cbrange [STATS_min:STATS_max]

set output 'graphs/displacement_alpha_1.0_beta_1.0.png'
set title "Displacement Map (Alpha=1.0, Beta=1.0)"
splot "data/u_alpha_1.0_beta_1.0.txt" using 1:2:3 with pm3d

