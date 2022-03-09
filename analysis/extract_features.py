import numpy as np
import pickle

feature_importances = np.load('importances.npy', allow_pickle=True) 
# select_edges = np.load('features.npy', allow_pickle=True) 
with open('features.pkl', 'rb') as f:
    select_edges = pickle.load(f)

# print (feature_importances)
i = 0
feature_importances_dict = {}
for edge in select_edges.keys():
    feature_importances_dict[edge] = feature_importances[i]
    # print (i)
    i += 1
print ("sort")
sort_select_edges = sorted(feature_importances_dict.items(), key=lambda item: item[1], reverse = True)
print (sort_select_edges[:5])


select_nodes = {}
select_edges = {}
for i in range(50):

    edge = sort_select_edges[i][0]
    select_edges[edge] = i
    array = edge.split("&")
    node1 = array[0]
    node2 = array[1]
    select_nodes[node1] = 1
    select_nodes[node2] = 1
print (select_nodes, select_edges)
# return select_nodes, select_edges