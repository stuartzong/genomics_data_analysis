import plotly.plotly as py
import plotly.graph_objs as go
import pandas as pd

df = pd.read_csv('https://raw.githubusercontent.com/plotly/datasets/master/2014_orange_stock.csv')
df.head()

trace = go.Scatter( x=df['AAPL_x'], y=df['AAPL_y'] )
data = [trace]

# IPython notebook
# py.iplot(data, filename='pandas-time-series')

url = py.plot(data, filename='pandas-time-series')
