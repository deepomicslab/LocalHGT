from sknetwork.ranking import PageRank
import numpy as np

def pr(adjacency):
    pagerank = PageRank()
    scores = pagerank.fit_predict(adjacency)
    sp_idx = np.argmax(scores)
    return sp_idx, scores