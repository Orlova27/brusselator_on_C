set term gif animate
set output "filename.gif"
set yrange [0:300]
do for [t=0:999000:1000] {
  outfile = sprintf('anim100a1b15data%d.png',t)
  infile = sprintf('filename%d.dat',t)
  set size 1,1; set origin 0,0
  set grid layerdefault
  set xlabel "x "
  set ylabel "y "

  set sample 11; set isosamples 11
  set pm3d map
  set palette
  set colorbox
  set lmargin 0


  set pm3d flush begin
  splot infile u 2:3:4 w pm3d
}
set output
