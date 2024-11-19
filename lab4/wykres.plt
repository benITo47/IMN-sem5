
# Wykresy dla różnych wartości omega_G i omega_L

set term pngcairo size 800,600 enhanced font 'Verdana,10'
set output 'wykresy/relaksacja_globalna.png'
set xlabel "Numer iteracji"
set ylabel "Całka funkcjonalna S"
set title "Zbieżność metody relaksacji globalnej dla różnych wartości omega_G"
set grid
set logscale x
set ytics 1000

stats 'dane/globalna_s_0.6.txt' nooutput
max_iter_0_6 = STATS_records

stats 'dane/globalna_s_1.0.txt' nooutput
max_iter_1_0 = STATS_records

set label 1 sprintf("Max iteracji = %d", max_iter_0_6) at graph 0.5, 0.95
set label 2 sprintf("Max iteracji = %d", max_iter_1_0) at graph 0.5, 0.90

plot 'dane/globalna_s_0.6.txt' using 1:2 with lines title "omega_G = 0.6", \
     'dane/globalna_s_1.0.txt' using 1:2 with lines title "omega_G = 1.0"

unset output
reset


set term pngcairo size 800,600 enhanced font 'Verdana,10'
set output 'wykresy/relaksacja_lokalna.png'
set xlabel "Numer iteracji"
set ylabel "Całka funkcjonalna S"
set title "Zbieżność metody relaksacji lokalnej dla różnych wartości omega_L"
set grid
set logscale x
set ytics 1000

stats 'dane/lokalna_s_1.0.txt' nooutput
max_iter_1_0_local = STATS_records

stats 'dane/lokalna_s_1.4.txt' nooutput
max_iter_1_4 = STATS_records

stats 'dane/lokalna_s_1.8.txt' nooutput
max_iter_1_8 = STATS_records

stats 'dane/lokalna_s_1.9.txt' nooutput
max_iter_1_9 = STATS_records

set label 1 sprintf("Max iteracji = %d", max_iter_1_0_local) at graph 0.5, 0.95
set label 2 sprintf("Max iteracji = %d", max_iter_1_4) at graph 0.5, 0.90
set label 3 sprintf("Max iteracji = %d", max_iter_1_8) at graph 0.5, 0.85
set label 4 sprintf("Max iteracji = %d", max_iter_1_9) at graph 0.5, 0.80

plot 'dane/lokalna_s_1.0.txt' using 1:2 with lines title "omega_L = 1.0", \
     'dane/lokalna_s_1.4.txt' using 1:2 with lines title "omega_L = 1.4", \
     'dane/lokalna_s_1.8.txt' using 1:2 with lines title "omega_L = 1.8", \
     'dane/lokalna_s_1.9.txt' using 1:2 with lines title "omega_L = 1.9"

unset output
reset




# Wykres potencjału V(x, y) dla omega = 0.6

set term pngcairo size 800,600 enhanced font 'Verdana,10'
set output 'wykresy/potencjal_ustalony_06.png'
set xlabel "x"
set ylabel "y"
set title "Mapa potencjału V(x, y) w stanie ustalonym omega = 0.6"
set pm3d map
set palette defined ( 0.0 'blue', 0.5 'white', 1.0 'red')
set xrange [0:15]   
set yrange [0:10]   
set cbrange [0:10]  

splot 'dane/globalna_v_0.6.txt' using 1:2:3 notitle

unset output
reset


# Wykres potencjału V(x, y) dla omega = 1.0

set term pngcairo size 800,600 enhanced font 'Verdana,10'
set output 'wykresy/potencjal_ustalony_1.png'
set xlabel "x"
set ylabel "y"
set title "Mapa potencjału V(x, y) w stanie ustalonym omega = 1.0"
set pm3d map
set palette defined ( 0.0 'blue', 0.5 'white', 1.0 'red')
set xrange [0:15]
set yrange [0:10]
set cbrange [0:10]

splot 'dane/globalna_v_1.0.txt' using 1:2:3 notitle

unset output
reset

# Wykres błędu numerycznego δ(x, y) dla omega = 1.0

set term pngcairo size 800,600 enhanced font 'Verdana,10'
set output 'wykresy/blad_numeryczny_1.png'
set xlabel "x"
set ylabel "y"
set title "Mapa błędu numerycznego δ(x, y) w stanie ustalonym omega = 1.0"
set pm3d map
set palette rgbformulae 6,3,-1   # Violet to Blue palette for error
set xrange [0:15]
set yrange [0:10]
set cbrange [0:0.002]

splot 'dane/globalna_err_1.0.txt' using 1:2:3 notitle

unset output
reset

# Wykres błędu numerycznego δ(x, y) dla omega = 0.6

set term pngcairo size 800,600 enhanced font 'Verdana,10'
set output 'wykresy/blad_numeryczny_06.png'
set xlabel "x"
set ylabel "y"
set title "Mapa błędu numerycznego δ(x, y) w stanie ustalonym omega = 0.6"
set pm3d map
set palette rgbformulae 6,3,-1   # Violet to Blue palette for error
set xrange [0:15]
set yrange [0:10]
set cbrange [0:0.002]

splot 'dane/globalna_err_0.6.txt' using 1:2:3 notitle

unset output
reset

