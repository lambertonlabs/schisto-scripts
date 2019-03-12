a = ['Keila','Lauren','Teteh','Elias','Billy','Rachel','Christina','Tristan','Suzan']

import random

random.shuffle(a)

for (i, x) in enumerate(a):

    print x, a[i-1]
