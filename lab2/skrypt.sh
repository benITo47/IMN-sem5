set terminal pngcairo enhanced
set output 'picarda_plot.png'

set title "Metoda Picarda"
set xlabel "Czas (t)"
set grid
set xrange [0:100]
set yrange [0:500]

plot "wyniki_picarda.dat" using 1:2 title 'u(t)' with lines linecolor rgb "blue", \
     "wyniki_picarda.dat" using 1:3 title 'Z = N - u' with lines linestyle 2 linecolor rgb "cyan"

unset output

set output 'newton_plot.png'

set title "Metoda Newtona"
set xlabel "Czas (t)"
set grid
set xrange [0:100]
set yrange [0:500]

plot "wyniki_newtona.dat" using 1:2 title 'u(t)' with lines linecolor rgb "orange", \
     "wyniki_newtona.dat" using 1:3 title 'Z = N - u' with lines linestyle 2 linecolor rgb "green"
unset output

set output 'rk_plot.png'

set title "Niejawna metoda RK2"
set xlabel "Czas (t)"
set grid 
set xrange [0:100]
set yrange [0:500]

plot "wyniki_rk.dat" using 1:2 title 'u(t)' with lines linecolor rgb "violet", \
     "wyniki_rk.dat" using 1:3 title 'Z = N - u' with lines linestyle 2 linecolor rgb 'red'


unset output

set output 'picard_iter.png'

set title "Iteracje w danym kroku czasowym - picard"
set xlabel "Czas (t)"
set grid 
set xrange [0:100]
set yrange [0:10]

plot "picard_iter.dat" using 1:2 title 'I(t)' with lines linecolor rgb "violet", \


unset output
set output 'newton_iter.png'

set title "Iteracje w danym kroku czasowym - newton"
set xlabel "Czas (t)"
set grid 
set xrange [0:100]
set yrange [0:10]

plot "newton_iter.dat" using 1:2 title 'I(t)' with lines linecolor rgb "violet", \


unset output
set output 'rk_iter.png'

set title "Iteracje w danym kroku czasowym - rk"
set xlabel "Czas (t)"
set grid 
set xrange [0:100]
set yrange [0:20]

plot "rk_iter.dat" using 1:2 title 'I(t)' with lines linecolor rgb "violet", \


unset output
