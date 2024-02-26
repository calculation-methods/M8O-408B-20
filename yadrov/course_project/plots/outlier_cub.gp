set term pdfcairo
set size ratio -1
set key off
set xrange [-0.6:6.6]
set yrange [-0.3:4.9]
set style line 1 lw 1.5 lt 4 pt 7 ps 0.5
$points << EOD
0 1
1 1
2 1
3 3.6
4 1
5 1
6 1
EOD
plot sample [0:1] - 0.3 * x**3 + 0 * x**2 + 0.3 * x**1 + 1 ls 1, \
[1:2] 1.5 * (x - 1)**3 - 0.9 * (x - 1)**2 - 0.6 * (x - 1)**1 + 1 ls 1, \
[2:3] - 3.1 * (x - 2)**3 + 3.6 * (x - 2)**2 + 2.1 * (x - 2)**1 + 1 ls 1, \
[3:4] 3.1 * (x - 3)**3 - 5.7 * (x - 3)**2 + 0 * (x - 3)**1 + 3.6 ls 1, \
[4:5] - 1.5 * (x - 4)**3 + 3.6 * (x - 4)**2 - 2.1 * (x - 4)**1 + 1 ls 1, \
[5:6] 0.3 * (x - 5)**3 - 0.9 * (x - 5)**2 + 0.6 * (x - 5)**1 + 1 ls 1, \
$points with points ls 1