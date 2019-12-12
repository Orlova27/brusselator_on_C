set terminal png 
set yrange [0.0:8.0]
do for [t=999:999999:1000] {
  outfile = sprintf('data%d.png',t)
  infile = sprintf('brus%d.dat',t)
  set output outfile
  plot infile u 2:3 w l lt rgb "#0480ad", infile u 2:4 w l lt rgb "#7c00ad"
}