import os
import re
import datetime
import tempfile
 
import numpy as np
import pandas as pd
import scipy.optimize as spo
import matplotlib
import matplotlib.pyplot as plt
 
 
def symbol_to_path(symbol):
    return os.path.join('stocks', '{0}.csv'.format(symbol))
 
 
def get_data(symbols, dates):
    """Read stock data (adjusted close) for given symbols from CSV files."""
    df = pd.DataFrame(index=dates)
    if 'SPY' not in symbols:  # add SPY for reference, if absent
        symbols.insert(0, 'SPY')
 
    for symbol in symbols:
        df_temp = pd.read_csv(symbol_to_path(symbol), index_col='Date',
                              parse_dates=True, usecols=['Date', 'Adj Close'], na_values=['nan'])
        df_temp = df_temp.rename(columns={'Adj Close': symbol})
        df = df.join(df_temp)
        if symbol == 'SPY':  # drop dates SPY did not trade
            df = df.dropna(subset=["SPY"])
    return df
 
 
def plot_data(df, title="Stock prices"):
    """Plot stock prices with a custom title and meaningful axis labels."""
    ax = df.plot(title=title, fontsize=12)
    ax.set_xlabel("Date")
    ax.set_ylabel("Price")
    plt.show()
 
 
def get_rolling_mean(values, window):
    """Return rolling mean of given values, using specified window size."""
    return values.rolling(window=window, center=False).mean()
 
 
def get_rolling_std(values, window):
    """Return rolling standard deviation of given values, using specified window size."""
    # TODO: Compute and return rolling standard deviation
    return values.rolling(window=window, center=False).std()
 
def get_bollinger_bands(rm, rstd):
    """Return upper and lower Bollinger Bands."""
    # TODO: Compute upper_band and lower_band
    upper_band = rm + 2 * rstd
    lower_band = rm - 2 * rstd
    return upper_band, lower_band
 
def compute_daily_returns(df):
    """Compute and return the daily return values."""
    # res = df.copy()
    # res[1:] = (df[1:] / df.values[:-1]) - 1
    # alternatively:
    res = (df / df.shift(1)) - 1
    if isinstance(df, pd.DataFrame):
        res.ix[0, :] = 0
    elif isinstance(df, pd.Series):
        res[0] = 0
    else:
        raise ValueError('unknown type: {0}'.format(type(df)))
    return res
 
def compute_cumulative_return(df):
    """Compute and return the daily return values."""
    res = (df / df.values[0]) - 1
    return res
 
def fill_missing_values(df):
    df.fillna(method='ffill', inplace=True)
    df.fillna(method='bfill', inplace=True)
