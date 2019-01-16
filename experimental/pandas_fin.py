import numpy as np
import pandas as pd

# As mentioned above a DataFrame is somewhat like a spreadsheet,
# or a structure for storing the data matrix in a regression
# While a Series is one individual column of data,
# a DataFrame is all the columns


# Series
s = pd.Series(np.random.randn(4), name='daily returns')
print s
print s.describe()
s.index = ['A1', 'B1', 'C1', 'D1']
print s['B1']
s['B1'] = 0
print s
print "is C1 in s: %s" % ('C1' in s)


# Dataframe
df = pd.read_csv('test_pwt.csv')
print type(df), df

#e can select particular rows using standard Python array slicing notation

print df[2:5]
# To select columns, we can pass a list containing the names of the desired columns represented as strings
print df[['country', 'tcgdp']]
#To select a mix of both we can use the ix attribute
print df.ix[2:5, ['country', 'tcgdp']]
