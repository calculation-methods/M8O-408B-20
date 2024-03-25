set term pdfcairo
set size ratio -1
set key off
set xrange [-0.284:6.364]
set yrange [-1:2]
set style line 1 lw 1.5 lt 4 pt 7 ps 0.5
$points << EOD
0.27 0.4
1.2 0.4
2.12 0.4
2.17 0.6992
2.3 0.9466
2.5 1.1448
2.76 1.2764
3.04 1.32
3.32 1.2764
3.58 1.1448
3.78 0.9466
3.91 0.6992
3.96 0.4
4.88 0.4
5.81 0.4
EOD
plot sample [0.27:1.2] - 1.03131 * (x - 0.27)**3 + 0 * (x - 0.27)**2 + 0.891978 * (x - 0.27)**1 + 0.4 ls 1, \
[1.2:2.12] 5.23525 * (x - 1.2)**3 - 2.87735 * (x - 1.2)**2 - 1.78396 * (x - 1.2)**1 + 0.4 ls 1, \
[2.12:2.17] - 323.868 * (x - 2.12)**3 + 11.5719 * (x - 2.12)**2 + 6.21507 * (x - 2.12)**1 + 0.4 ls 1, \
[2.17:2.3] 104.786 * (x - 2.17)**3 - 37.0082 * (x - 2.17)**2 + 4.94326 * (x - 2.17)**1 + 0.6992 ls 1, \
[2.3:2.5] - 10.361 * (x - 2.3)**3 + 3.85832 * (x - 2.3)**2 + 0.633776 * (x - 2.3)**1 + 0.9466 ls 1, \
[2.5:2.76] 2.74438 * (x - 2.5)**3 - 2.35827 * (x - 2.5)**2 + 0.933785 * (x - 2.5)**1 + 1.1448 ls 1, \
[2.76:3.04] - 0.604405 * (x - 2.76)**3 - 0.217656 * (x - 2.76)**2 + 0.264043 * (x - 2.76)**1 + 1.2764 ls 1, \
[3.04:3.32] 0.604405 * (x - 3.04)**3 - 0.725356 * (x - 3.04)**2 - 7.21645e-16 * (x - 3.04)**1 + 1.32 ls 1, \
[3.32:3.58] - 2.74438 * (x - 3.32)**3 - 0.217656 * (x - 3.32)**2 - 0.264043 * (x - 3.32)**1 + 1.2764 ls 1, \
[3.58:3.78] 10.361 * (x - 3.58)**3 - 2.35827 * (x - 3.58)**2 - 0.933785 * (x - 3.58)**1 + 1.1448 ls 1, \
[3.78:3.91] - 104.786 * (x - 3.78)**3 + 3.85832 * (x - 3.78)**2 - 0.633776 * (x - 3.78)**1 + 0.9466 ls 1, \
[3.91:3.96] 323.868 * (x - 3.91)**3 - 37.0082 * (x - 3.91)**2 - 4.94326 * (x - 3.91)**1 + 0.6992 ls 1, \
[3.96:4.88] - 5.23525 * (x - 3.96)**3 + 11.5719 * (x - 3.96)**2 - 6.21507 * (x - 3.96)**1 + 0.4 ls 1, \
[4.88:5.81] 1.03131 * (x - 4.88)**3 - 2.87735 * (x - 4.88)**2 + 1.78396 * (x - 4.88)**1 + 0.4 ls 1, \
$points with points ls 1