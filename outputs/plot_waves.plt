stats 'waves.txt' matrix nooutput
mid = floor(STATS_size_x/2)

plot "waves.txt" every 10 using 1:2 with lines title "State Number 1", \
	 "waves.txt" every 10 using 1:3 with lines title "State Number 2", \
	 "waves.txt" every 10 using 1:4 with lines title "State Number 3", \
	 
set key top right