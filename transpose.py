
def select(value, rowheaders, columnheaders, table, where = ''):
  import sqlite3
  from collections import defaultdict
  
  conn = sqlite3.connect('test.db')
  
  if where:
    where = 'where ' + where
  
  c = conn.execute(
    '''select {2}, {1}, {0} from {3} {4} group by {2}, {1}'''
      .format(value, rowheaders, columnheaders, table, where))
  
  s = set()
  d = defaultdict(dict)
  for row in c:
    s.add(row[0])
    d[row[1]][row[0]] = row[2]
  
  conn.close()

  print(end='\t')
  print(*s, sep='\t')
  
  for k in sorted(d):
    print(k, end='\t')
    for h in s:
      print(d[k][h], end='\t')
    print()
    

