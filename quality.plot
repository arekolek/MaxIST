
set xlabel "Number of internal vertices"
set ylabel "Running time in ms"
set key left
set style fill transparent solid 0.25 noborder
set term png truecolor enhanced font "Verdana,15"

#internal(cons, improv) = '< sqlite3 test.db "select internal, avg(time), count(*) from quality where type=\"gnp\" and vertices=500 and parameter=0.003 and construction=\"'.cons.'\" and improvement=\"'.improv.'\" group by internal;" | tr "|" "\t"'

#internal(cons, improv) = '< sqlite3 test.db "select internal, avg(time), count(*) from quality where type=\"gnp\" and vertices=500 and parameter=0.003 and construction=\"'.cons.'\" and improvement=\"'.improv.'\" group by round(100.0*(internal-200)/300)+100*round(10.0*(time-500)/1000);" | tr "|" "\t"'

clustered(cons, improv) = '< python3 quality.py '.cons.' '.improv.' 20'

all(cons, improv) = '< sqlite3 test.db "select internal, time from quality where type=\"gnp\" and vertices=500 and parameter=0.003 and construction=\"'.cons.'\" and improvement=\"'.improv.'\";" | tr "|" "\t"'

cons = 'rdfs-sort rdfs-rand rdfs fifo dfs bfs'
imps = 'none'
do for [i in imps] {
  set title "Construction algorithms"
  set output 'quality-construction.png'
  plot for [c in cons] clustered(c, i) using 2:3:($1*0.2) title c with circles
}

cons = 'rdfs-sort rdfs-rand rdfs'
imps = 'none'
do for [i in imps] {
  set title "Construction algorithms"
  set output 'quality-construction-closeup.png'
  plot for [c in cons] clustered(c, i) using 2:3:($1*0.2) title c with circles #, \for [c in cons] all(c, i) title c with points
}

cons = 'bfs dfs fifo'
imps = 'lost-ex lost lost-light prieto none'
do for [c in cons] {
  outfile = sprintf('quality-%s.png', c)
  set title "Quality of improvement over ".c
  set output outfile
  plot for [i in imps] clustered(c, i) using 2:3:($1*0.2) title i with circles
}

cons = 'rdfs rdfs-rand rdfs-sort'
imps = 'lost-ex lost lost-light prieto none'
do for [c in cons] {
  outfile = sprintf('quality-%s.png', c)
  set title "Quality of improvement over ".c
  set output outfile
  plot for [i in imps] clustered(c, i) using 2:3:($1*0.1) title i with circles
}

cons = 'bfs'
imps = 'lost-ex lost lost-light'
do for [c in cons] {
  outfile = sprintf('quality-%s-closeup.png', c)
  set title "Quality of improvement over ".c
  set output outfile
  plot for [i in imps] clustered(c, i) using 2:3:($1*0.05) title i with circles
}

cons = 'dfs fifo'
imps = 'lost-ex lost lost-light'
do for [c in cons] {
  outfile = sprintf('quality-%s-closeup.png', c)
  set title "Quality of improvement over ".c
  set output outfile
  plot for [i in imps] clustered(c, i) using 2:3:($1*0.1) title i with circles
}
