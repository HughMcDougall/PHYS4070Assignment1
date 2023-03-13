set logscale xy
plot "potential.txt" every 2 using 1:($2*-1) with lines title "Potential"

set key top right