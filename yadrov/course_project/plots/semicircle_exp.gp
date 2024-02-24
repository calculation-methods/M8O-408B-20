set term pdfcairo
set size ratio -1
set key off
set xrange [-0.284:6.364]
set yrange [-1:2]
set style line 1 lw 1.5 lt 1 pt 7 ps 0.5
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
plot sample [0.27:1.2] (0 * sinh(10 * (1.2 - x)) + -0.422448 * sinh(10 * (x - 0.27))) / 546901 + 0.4 * (1.2 - x) / 0.93 + 0.404224 * (x - 0.27) / 0.93 ls 1, \
[1.2:2.12] (-0.422448 * sinh(80 * (2.12 - x)) + 252.661 * sinh(80 * (x - 1.2))) / 2.94594e+35 + 0.400066 * (2.12 - x) / 0.92 + 0.360522 * (x - 1.2) / 0.92 ls 1, \
[2.12:2.17] (252.661 * sinh(0.1 * (2.17 - x)) + -161.077 * sinh(0.1 * (x - 2.12))) / 5.00002e-05 + -25265.7 * (2.17 - x) / 0.05 + 16108.4 * (x - 2.12) / 0.05 ls 1, \
[2.17:2.3] (-161.077 * sinh(36.2849 * (2.3 - x)) + -0.56139 * sinh(36.2849 * (x - 2.17))) / 73616.3 + 0.821543 * (2.3 - x) / 0.13 + 0.947026 * (x - 2.17) / 0.13 ls 1, \
[2.3:2.5] (-0.56139 * sinh(29.2539 * (2.5 - x)) + -3.94365 * sinh(29.2539 * (x - 2.3))) / 148693 + 0.947256 * (2.5 - x) / 0.2 + 1.14941 * (x - 2.3) / 0.2 ls 1, \
[2.5:2.76] (-3.94365 * sinh(0.05 * (2.76 - x)) + -0.649122 * sinh(0.05 * (x - 2.5))) / 3.25009e-05 + 1578.61 * (2.76 - x) / 0.26 + 260.925 * (x - 2.5) / 0.26 ls 1, \
[2.76:3.04] (-0.649122 * sinh(0.05 * (3.04 - x)) + -1.34383 * sinh(0.05 * (x - 2.76))) / 3.50011e-05 + 260.925 * (3.04 - x) / 0.28 + 538.852 * (x - 2.76) / 0.28 ls 1, \
[3.04:3.32] (-1.34383 * sinh(0.05 * (3.32 - x)) + -0.649122 * sinh(0.05 * (x - 3.04))) / 3.50011e-05 + 538.852 * (3.32 - x) / 0.28 + 260.925 * (x - 3.04) / 0.28 ls 1, \
[3.32:3.58] (-0.649122 * sinh(0.05 * (3.58 - x)) + -3.94365 * sinh(0.05 * (x - 3.32))) / 3.25009e-05 + 260.925 * (3.58 - x) / 0.26 + 1578.61 * (x - 3.32) / 0.26 ls 1, \
[3.58:3.78] (-3.94365 * sinh(29.2539 * (3.78 - x)) + -0.56139 * sinh(29.2539 * (x - 3.58))) / 148693 + 1.14941 * (3.78 - x) / 0.2 + 0.947256 * (x - 3.58) / 0.2 ls 1, \
[3.78:3.91] (-0.56139 * sinh(36.2849 * (3.91 - x)) + -161.077 * sinh(36.2849 * (x - 3.78))) / 73616.3 + 0.947026 * (3.91 - x) / 0.13 + 0.821543 * (x - 3.78) / 0.13 ls 1, \
[3.91:3.96] (-161.077 * sinh(0.1 * (3.96 - x)) + 252.661 * sinh(0.1 * (x - 3.91))) / 5.00002e-05 + 16108.4 * (3.96 - x) / 0.05 + -25265.7 * (x - 3.91) / 0.05 ls 1, \
[3.96:4.88] (252.661 * sinh(80 * (4.88 - x)) + -0.422448 * sinh(80 * (x - 3.96))) / 2.94594e+35 + 0.360522 * (4.88 - x) / 0.92 + 0.400066 * (x - 3.96) / 0.92 ls 1, \
[4.88:5.81] (-0.422448 * sinh(10 * (5.81 - x)) + 0 * sinh(10 * (x - 4.88))) / 546901 + 0.404224 * (5.81 - x) / 0.93 + 0.4 * (x - 4.88) / 0.93 ls 1, \
$points with points ls 1