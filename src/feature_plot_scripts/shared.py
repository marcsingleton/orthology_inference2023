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


def top_vars(df, i):
    return df.var().sort_values(ascending=False).index[:i]


def znorm(df):
    avg = df.mean()
    std = df.std()
    z = (df - avg) / std
    return z


n_components = 5  # number of PCA components
fsets = {'minmax': lambda df: minmax(df, ['net_charge', 'net_charge_P', 'ED_ratio', 'RK_ratio', 'SCD', 'loglen', 'iso_point']),
         'norm': lambda df: norm(df.drop(['net_charge', 'net_charge_P', 'ED_ratio', 'RK_ratio', 'SCD', 'loglen', 'iso_point'], axis=1), [('NCPR', 2)]),
         'znorm': znorm,
         'drop_1': lambda df: df.drop(top_vars(df, 1), axis=1),
         'drop_2': lambda df: df.drop(top_vars(df, 2), axis=1),
         'drop_3': lambda df: df.drop(top_vars(df, 3), axis=1),
         'drop_4': lambda df: df.drop(top_vars(df, 4), axis=1),
         'drop_5': lambda df: df.drop(top_vars(df, 5), axis=1)}
labels = {(True, False): 'κ = -1, Ω ≠ -1',
          (False, True): 'κ ≠ -1, Ω = -1',
          (True, True): 'κ = -1, Ω = -1',
          (False, False): 'κ ≠ -1, Ω ≠ -1'}
