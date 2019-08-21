"""Shared functions and parameters between all feature projection analyses."""

def minmax(df, feats):
    norm = df.copy()
    for col in feats:
        cmin = min(norm[col])
        cmax = max(norm[col])
        norm[col] = (norm[col] - cmin) / (cmax - cmin)
    return norm


def norm(df, featnorms):
    for feat, size in featnorms:
        df[feat] = df[feat] / size
    return df


def znorm(df):
    avg = df.mean()
    std = df.std()
    z = (df - avg) / std
    return z


n_components = 5  # number of PCA components
fsets = {'all': lambda df: df,
         '-len': lambda df: df.drop('loglen', axis=1),
         '-net_charge': lambda df: df.drop('net_charge', axis=1),
         '-len-net_charge': lambda df: df.drop(['loglen', 'net_charge'], axis=1),
         'minmax': lambda df: minmax(df, ['net_charge', 'net_charge_P', 'ED_ratio', 'RK_ratio', 'SCD', 'loglen', 'iso_point']),
         'norm': lambda df: norm(df.drop(['net_charge', 'net_charge_P', 'ED_ratio', 'RK_ratio', 'SCD', 'loglen', 'iso_point'], axis=1), [('NCPR', 2)]),
         'znorm': znorm}
labels = {(True, False): 'κ = -1, Ω ≠ -1',
          (False, True): 'κ ≠ -1, Ω = -1',
          (True, True): 'κ = -1, Ω = -1',
          (False, False): 'κ ≠ -1, Ω ≠ -1'}
