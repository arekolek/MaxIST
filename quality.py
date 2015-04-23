
import sqlite3
import numpy as np

conn = sqlite3.connect('test.db')
c = conn.cursor()

def get_data(construction, improvement):
    return np.array(c.execute(
        "SELECT internal, time FROM quality \
            WHERE type='gnp' \
                AND vertices=500 AND parameter=0.003 \
                AND construction=? AND improvement=?",
        (construction, improvement)).fetchall())

def scale_linear_bycolumn(rawpoints, high=100.0, low=0.0):
    mins = np.min(rawpoints, axis=0)
    maxs = np.max(rawpoints, axis=0)
    rng = maxs - mins
    return high - (((high - low) * (maxs - rawpoints)) / rng)


from sklearn.cluster import KMeans
from sys import argv

X = get_data(argv[1], argv[2])
X_scaled = scale_linear_bycolumn(X)
km = KMeans(int(argv[3]))
km.fit(X_scaled)
labels = km.labels_

for k, n in np.dstack(np.unique(labels, return_counts=True))[0]:
    print(n, *np.average(X[labels == k], axis=0), sep='\t')
