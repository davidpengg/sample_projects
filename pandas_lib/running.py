# a scrap paper for unit tests. not organized to run as its own file (except for the larger dataset test)

# trying larger? datasets
import seaborn as sns
import time
iris = sns.load_dataset('iris')
sepal_length = list(iris['sepal_length'])
sepal_width = list(iris['sepal_width'])
species = list(iris['species'])

from mypandas import *
start = time.time()
df1 = DataFrame({'sepal_length': Series(sepal_length), 'sepal_width': Series(sepal_width), 'species': Series(species)})
print("mypandas", time.time() - start)

import pandas as pd
start = time.time()
df2 = pd.DataFrame({'sepal_length': pd.Series(sepal_length), 'sepal_width': pd.Series(sepal_width), 'species': pd.Series(species)})
print("pandas", time.time() - start)
exit()

# contains unit tests for mypandas and experimental tests with python pandas

from pandas import * ; s1 = Series([True]); s2 = Series([False])
from pandas import * ; s1 = Series(['hello']); s2 = Series(['there'])
from pandas import * ; s1 = Series([1]); s2 = Series([2])
from pandas import * ; s1 = Series([1, 3]); s2 = Series([2])
from pandas import * ; s1 = Series([1.2]); s2 = Series([3.4])

from pandas import * ; s1 = Series([1]); s2 = Series(['there', 'ok'])



from mypandas import * ; s1 = Series([True]); s2 = Series([False])
from mypandas import * ; s1 = Series([True, False]); s2 = Series([False])
from mypandas import * ; s1 = Series(['hello']); s2 = Series(['there'])
from mypandas import * ; s1 = Series([1]); s2 = Series([2])
from mypandas import * ; s1 = Series([1, 3]); s2 = Series([2])
from mypandas import * ; s1 = Series([1.2]); s2 = Series([3.4])


# ==================
from pandas import * ; s1 = Series([1, 3, 1, 3]); s2 = Series([2, 3, 2, 3])
from mypandas import * ; s1 = Series([1, 3, 1, 3]); s2 = Series([2, 3, 2, 3])

from pandas import * ; s1 = Series(['a','b','c','d']) ; b = Series([True, False, True, False])
from mypandas import * ; s1 = Series(['a','b','c','d']) ; b = Series([True, False, True, False])


from mypandas import * ; s1 = Series({3:'a', 6: 'b', 9:'c', 12:'d'}) ; b = Series([True, False, True, False])
#===========
# trying stuff with None
from mypandas import * ; s1 = Series(['a','b',None,'d']) ; b = Series([True, False, True, False])

# dataframes

from pandas import * ; df = DataFrame({"SKU": Series(["X4E", "T3B", "F8D", "C7X"]), "price": Series([7.0, 3.5, 8.0, 6.0]), "sales" : Series([5, 3, 1, 10]), "taxed" : Series([False, False, True, False])})

result = df[(df["price"] + 5.0 > 10.0) & (df["sales"] > 3) & ~df["taxed"]]["SKU"]
print(result)


from pandas import * ; df = DataFrame({"SKU": Series(["X4E", "T3B", "F8D",]), "price": Series([7.0, 3.5, 8.0, 6.0]), "sales" : Series([5, 3, 1, 10]), "taxed" : Series([False, False, True, False])})

from mypandas import * ; df = DataFrame({"SKU": Series(["X4E", "T3B", "F8D", "C7X"]), "price": Series([7.0, 3.5, 8.0, 6.0]), "sales" : Series([5, 3, 1, 10]), "taxed" : Series([False, False, True, False])})

from mypandas import * ; df = DataFrame({"SKU": Series(["X4E", "T3B", "F8D", "C7X"]), "price": Series([7.0, 3.5, 8.0, 6.0]), "sales" : Series([5, 3, 1, 10]), "taxed" : Series([False, False, True, False])}); b = Series([True, False, True, False])

result = df[(df["price"] + 5.0 > 10.0) & (df["sales"] > 3) & ~df["taxed"]]["SKU"]
print(result)

from mypandas import * ; path = "experiment.csv" ; df = read_csv(path)

