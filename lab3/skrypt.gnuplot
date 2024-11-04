
# Skrypt Gnuplot dla generowania porównawczych wykresów z plików wynikowych

set terminal pngcairo size 800,600 enhanced font 'Arial,10'

set style line 1 lt 1 lw 2 lc rgb "blue"  # TOL = 10^-2
set style line 2 lt 1 lw 2 lc rgb "red"   # TOL = 10^-5

# Zestaw 1: Metoda RK2, wykres x(t) dla obu TOL
set output 'RK2_x_t.png'
set title "RK2: x(t) dla TOL=10^-2 i TOL=10^-5"
set xlabel "t"
set ylabel "x(t)"
plot 'RK2_tol1.txt' using 1:3 with lines linestyle 1 title "TOL=10^-2", \
     'RK2_tol2.txt' using 1:3 with lines linestyle 2 title "TOL=10^-5"

set output 'RK2_v_t.png'
set title "RK2: v(t) dla TOL=10^-2 i TOL=10^-5"
set xlabel "t"
set ylabel "v(t)"
plot 'RK2_tol1.txt' using 1:4 with lines linestyle 1 title "TOL=10^-2", \
     'RK2_tol2.txt' using 1:4 with lines linestyle 2 title "TOL=10^-5"

set output 'RK2_dt_t.png'
set title "RK2: Δt(t) dla TOL=10^-2 i TOL=10^-5"
set xlabel "t"
set ylabel "Δt(t)"
plot 'RK2_tol1.txt' using 1:2 with lines linestyle 1 title "TOL=10^-2", \
     'RK2_tol2.txt' using 1:2 with lines linestyle 2 title "TOL=10^-5"

set output 'RK2_v_x.png'
set title "RK2: v(x) dla TOL=10^-2 i TOL=10^-5"
set xlabel "x"
set ylabel "v"
plot 'RK2_tol1.txt' using 3:4 with lines linestyle 1 title "TOL=10^-2", \
     'RK2_tol2.txt' using 3:4 with lines linestyle 2 title "TOL=10^-5"

set output 'trapezy_x_t.png'
set title "Trapezy: x(t) dla TOL=10^-2 i TOL=10^-5"
set xlabel "t"
set ylabel "x(t)"
plot 'trapezy_tol1.txt' using 1:3 with lines linestyle 1 title "TOL=10^-2", \
     'trapezy_tol2.txt' using 1:3 with lines linestyle 2 title "TOL=10^-5"

set output 'trapezy_v_t.png'
set title "Trapezy: v(t) dla TOL=10^-2 i TOL=10^-5"
set xlabel "t"
set ylabel "v(t)"
plot 'trapezy_tol1.txt' using 1:4 with lines linestyle 1 title "TOL=10^-2", \
     'trapezy_tol2.txt' using 1:4 with lines linestyle 2 title "TOL=10^-5"

set output 'trapezy_dt_t.png'
set title "Trapezy: Δt(t) dla TOL=10^-2 i TOL=10^-5"
set xlabel "t"
set ylabel "Δt(t)"
plot 'trapezy_tol1.txt' using 1:2 with lines linestyle 1 title "TOL=10^-2", \
     'trapezy_tol2.txt' using 1:2 with lines linestyle 2 title "TOL=10^-5"

set output 'trapezy_v_x.png'
set title "Trapezy: v(x) dla TOL=10^-2 i TOL=10^-5"
set xlabel "x"
set ylabel "v"
plot 'trapezy_tol1.txt' using 3:4 with lines linestyle 1 title "TOL=10^-2", \
     'trapezy_tol2.txt' using 3:4 with lines linestyle 2 title "TOL=10^-5"

