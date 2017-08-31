# Set terminal and output
set terminal pngcairo size 1000,1000 enhanced font 'Verdana,15'
#set terminal pngcairo size 4000, 4000 enhanced color font 'Helvetica,80' linewidth 10
#set terminal postscript eps size 4,4 enhanced color font 'Helvetica,18' linewidth 3
#set terminal  postscript enhanced 18
set output 'P2.png'
# Set various features of the plot
set pm3d
unset surface  # don't need surfaces
set view map
#set contour
#set cntrparam bspline #cubicspline  # smooth out the lines
#set cntrparam levels -2    # sets the num of contour lines
set pm3d interpolate 20,20 # interpolate the color

# Set a nice color palette
set palette defined ( 0 "red", 1 "green", 2 "blue") 
#set palette defined ( 0 0 0 0, 1 1 1 1 )
#set palette functions gray, gray, gray
#set palette functions sqrt(gray), gray**3, sin(gray*2*pi)
#set title 'Free Energy Profile ({/Symbol D}F)'

# Axes
#set xlabel 'X(A)' font 'Helvetica,100'
#set ylabel 'y (A)' rotate by 90 font 'Helvetica,,100'
#set cblabel 'S' #rotate by 0 offset -9.,16.5,0 font 'Times Bold Italic,35'
#set cblabel 'S' rotate by 0 offset -7.,16.5,0 font 'Verdana,20'
set label 'x'  at  -14,-190 front font 'Helvetica,25'
set label '(A)'  at  0,-190.5 front font 'Helvetica,20'
set label 'y'  at  -193,-14 front rotate by 90 font 'Helvetica,25' 
set label '(A)'  at  -192.5,0 rotate by 90 front font 'Helvetica,20'

#unset tics
set tics nomirror
set tics out
set xtics offset 0,0.5 font 'Helvetica,18'
#set ytics #offset 0.0,0.5 rotate by 90 font 'Helvetica,80'
#set noytics #set 0.2,0 font 'Helvetica,80'
#set cbtics offset -1.,0 font 'Helvetica,80'
#set xtics 17.2 font 'Verdana,13'
#set ytics 17.2 font 'Verdana,13'
#set mxtics 2
#set mytics 2
#set style line 12 lc rgb 'blue' lt 1 lw 2
#set style line 13 lc rgb 'yellow' lt 1 lw 1
#set style line 14 lc rgb 'black' lt 1 lw 1
#set grid xtics ytics mxtics mytics #ls 12, ls 13
#set grid mxtics mytics ls 14


set xrange [-154.502:154.502]
set yrange [-154.502:154.502]
#set border 3

set cbrange [0.3:0.7]
#unset colorbox
# Now plot
#set multiplot
splot 'Output_scaler' using 1:2:3 notitle with lines lt 1, 'Output_vector' u ($1-5*$6):($2-5*$7):3:(10*$6):(10*$7):8 notitle with vector nohead lt 1 lc 0 lw 2
#splot 'Output_vector' u ($1-6*$6):($2-6*$7):3:(12*$6):(12*$7):8 notitle with vector nohead lt 1 lc 0 lw 2.5
#splot 'Output_vector' u ($1-5*$6):($2-5*$7):3:(10*$6):(10*$7):8 notitle with vector nohead lt 1 lc 0 lw 4

#r = 21
zz = 0
set parametric
set vrange [0:2*pi]
set urange [0:50]
fx(u,v) = u*cos(v)
fy(u,v) = u*sin(v)
unset pm3d
set surface
set iso 60

#splot fx(u,v),fy(u,v),zz notitle w l lt 2 lw 4 lc rgb "White"

unset multiplot
#pause -1
