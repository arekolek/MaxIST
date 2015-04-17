
def five_three(k):
  n = 6*k + 3
  m = 3 + k*8 + (k-1)*2
  
  a = n-3
  b = n-2
  c = n-1
  
  print(*(n, m))
  
  print(*(a, b))
  print(*(b, c))
  print(*(c, 4))
  
  for i in range(k):
    print(*(6*i + 0, b))

    print(*(6*i + 0, 6*i + 1))
    print(*(6*i + 1, 6*i + 2))
    print(*(6*i + 1, 6*i + 3))
    print(*(6*i + 1, 6*i + 5))
    print(*(6*i + 2, 6*i + 4))
    print(*(6*i + 3, 6*i + 4))
    print(*(6*i + 4, 6*i + 5))

    if i < k-1:
      print(*(6*i + 1, 6*(i+1) + 4))
      print(*(6*i + 5, 6*(i+1) + 5))

'''
print(10)
for k in range(2, 100, 10):
  five_three(k)
'''

'''
print(1)
five_three(4)
'''

def rank_same(*s):
  from itertools import chain
  print('{rank=same; ', end='')
  print(*chain(*s), end='')
  print('}')

def dot(k):
  n = 6*k + 3
  a = n-3
  b = n-2
  c = n-1
  rank_same([a, b])
  rank_same(*([6*i + 0, 6*i + 1] for i in range(k)))
  rank_same(*([6*i + 2, 6*i + 3] for i in range(k)))
  rank_same(6*i + 5 for i in range(k))
  rank_same(6*i + 6 for i in range(k))
  
  print('\nedge[style=invis]')
  for i in range(k-1):
    print('{} -- {};'.format(6*i + 1, 6*(i+1) + 0))
  from itertools import chain
  print( '{};'.format(' -- '.join(map(str, chain(*([6*i+2, 6*i+3] for i in range(k))) ))) )
  print('{} -- {};'.format(a, c))
  for i in range(k):
    print('{} -- {};'.format(6*i + 3, 6*i + 5))
  for i in range(k):
    print('{} -- {};'.format(6*i + 0, 6*i + 2))

dot(4)

