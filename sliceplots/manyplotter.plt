do for [t=0:999] {
  outfile = sprintf('data', "%d", t, '.png')
  infile = sprintf('brus', "%d", t, '.dat')
  set output outfile
  plot infile u 1:3
}
