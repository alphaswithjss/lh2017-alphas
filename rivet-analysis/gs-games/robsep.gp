reset
set ter pdfcairo enhanced size 11cm,9cm
set out 'robsep.pdf'
set colors classic

set xlabel '{/*1.4 {/Symbol D}_{robustness}}'
set xrange [0:0.2]

set ylabel '{/*1.4 {/Symbol D}_{sensitivity}}'
set yrange [0:0.007]

set grid

set label 1 'Z+q' at graph 0.795,0.95
set label 2 'Z+g' at graph 0.895,0.95
set key maxrows 3 width -7.1 spacing 1.1
set key at graph 1.0,0.92

set arrow 1 from 0.018,0.0062 to 0.005,0.0068 lc 7 lw 3
set label 11 '{/*1.2 Better}' at 0.013,0.0067

plot '< mergeidx.pl -f labels.res "Z\+q"' u ($3==0.5 ? $5 : 1/0):4 w p pt 4 ps 0.8 lc 3 lw 2 t '{/Symbol a}=0.5',\
     '< mergeidx.pl -f labels.res "Z\+q"' u ($3==1.0 ? $5 : 1/0):4 w p pt 6 ps 0.8 lc 3 lw 2 t '{/Symbol a}=1',\
     '< mergeidx.pl -f labels.res "Z\+q"' u ($3==2.0 ? $5 : 1/0):4 w p pt 8 ps 0.8 lc 3 lw 2 t '{/Symbol a}=2',\
     '< mergeidx.pl -f labels.res "Z\+g"' u ($3==0.5 ? $5 : 1/0):4 w p pt 5 ps 0.8 lc 1 lw 1 t ' ',\
     '< mergeidx.pl -f labels.res "Z\+g"' u ($3==1.0 ? $5 : 1/0):4 w p pt 7 ps 0.8 lc 1 lw 1 t ' ',\
     '< mergeidx.pl -f labels.res "Z\+g"' u ($3==2.0 ? $5 : 1/0):4 w p pt 9 ps 0.8 lc 1 lw 1 t ' ',\
     0.009-0.000365/(x+0.1) w l dt (8,8) lw 4 lc 7 not
     

set out
