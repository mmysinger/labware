"""
Copyright (C) 2015 Michael M Mysinger

Heirarchical clustering

Adapted from original code by Leo Gendelev & Garrett Gaskins

Michael Mysinger 201504 Created
"""

import numpy
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, fcluster

DEFAULT_THRESHOLD = 0.5

def hierarchical_y(xy_data, y_labels, metric="correlation",
                   threshold=DEFAULT_THRESHOLD):
    # Compute hierarchical clusters
    dist_matrix = pdist(xy_data, metric=metric)
    linkage_matrix = linkage(dist_matrix, metric=metric)
    linkage_matrix[linkage_matrix < 0] = 0
    threshold = threshold * max(dist_matrix)
    flat_clusters = fcluster(linkage_matrix, threshold, depth=4)
    # Arrange labels into the clusters
    cluster_dict = {}
    sums = xy_data.sum(axis=1)
    for index, cluster in enumerate(flat_clusters):
        score, members = cluster_dict.get(cluster, (0.0, []))
        item_score = sums[index]
        score += item_score
        members.append((item_score, y_labels[index]))
        cluster_dict[cluster] = (score, members) 
    # Reorder labels using the hierarchical clustering
    out_labels = []
    for cluster in sorted(cluster_dict, reverse=True,
                          key=lambda x: cluster_dict[x][0]):
        print cluster, cluster_dict[cluster]
        for item_score, label in sorted(cluster_dict[cluster][1], reverse=True):
            out_labels.append(label)
    return out_labels

def hierarchical_x(xy_data, x_labels, metric="correlation", 
                   threshold=DEFAULT_THRESHOLD):
    return hierarchical_y(xy_data.transpose(), x_labels, metric=metric,
                          threshold=threshold)

def hierarchical_both(xy_data, x_labels, y_labels, metric="correlation", 
                   threshold=DEFAULT_THRESHOLD):
    out_x_labels = hierarchical_x(xy_data, x_labels, metric=metric,
                                  threshold=threshold)
    out_y_labels = hierarchical_y(xy_data, y_labels, metric=metric,
                                  threshold=threshold)
    return out_x_labels, out_y_labels
