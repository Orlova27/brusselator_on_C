set terminal png 
set yrange [1.3:2.4]
do for [t=0:999] {
  outfile = sprintf('data%d.png',t)
  infile = sprintf('brus%d.dat',t)
  set output outfile
  plot infile u 2:3 w lp
}
