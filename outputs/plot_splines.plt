stats 'splines.txt' matrix nooutput
mid = floor(STATS_size_x/2)

plot "splines.txt" every 10 using 1:2 with lines title "Spline ".mid, \
	 "splines.txt" every 10 using 1:mid with lines title "Spline diff".mid,\
	 "splines.txt" every 10 using 1:STATS_size_x with lines title "Spline diff".mid
set key top right