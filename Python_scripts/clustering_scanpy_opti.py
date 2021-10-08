# Imports
import sys
import os
import time
import h5py
import numpy as np
import loompy
import json
from scipy.sparse import issparse, coo_matrix, csr_matrix
from sklearn.metrics import pairwise_distances
from umap.umap_ import nearest_neighbors
from umap.umap_ import fuzzy_simplicial_set

# Arguments
loom_file = 'MISSING'
iAnnot = 'MISSING'
oAnnot = 'MISSING'
oJSON = 'MISSING'
method = 'MISSING'
n_neighbors = 'MISSING'
metric = 'euclidean'
resolution = 'MISSING'
random_state = 'MISSING'
time_idle = 0
warnings = []

# Check input args
if len(sys.argv) >= 2: loom_file = sys.argv[1]
if len(sys.argv) >= 3: iAnnot = sys.argv[2]
if len(sys.argv) >= 4: oAnnot = sys.argv[3]
if len(sys.argv) >= 5: oJSON = sys.argv[4]
if len(sys.argv) >= 6: method = sys.argv[5]
if len(sys.argv) >= 7: n_neighbors = int(sys.argv[6])
if len(sys.argv) >= 8: resolution = float(sys.argv[7])
if len(sys.argv) >= 9: random_state = int(sys.argv[8])

# Functions
# Handling errors
def error_json(message, jsonfile):
    print(message, file=sys.stderr)
    if not jsonfile.startswith('MISSING'):
        data_json = {}
        data_json['displayed_error'] = message
        with open(jsonfile, 'w') as outfile:
            json.dump(data_json, outfile)
    sys.exit()

# Open the Loom file while handling potential locking
def open_with_lock(loomfile, mode):
    global time_idle, oJSON
    while True:
        try:
            return h5py.File(loomfile, mode)
        except Exception as e:
            if not "unable to lock file" in format(e): error_json("Error opening Loom file:" + format(e), oJSON)
        print("Sleeping 1sec for file lock....")
        time_idle += 1
        time.sleep(1)
        
def _compute_connectivities_umap(knn_indices, knn_dists, n_obs, n_neighbors, set_op_mix_ratio=1.0, local_connectivity=1.0):
    X = coo_matrix(([], ([], [])), shape=(n_obs, 1))
    connectivities = fuzzy_simplicial_set(X, n_neighbors, None, None, knn_indices=knn_indices, knn_dists=knn_dists, set_op_mix_ratio=set_op_mix_ratio, local_connectivity=local_connectivity)

    if isinstance(connectivities, tuple):
        # In umap-learn 0.4, this returns (result, sigmas, rhos)
        connectivities = connectivities[0]

    distances = _get_sparse_matrix_from_indices_distances_umap(knn_indices, knn_dists, n_obs, n_neighbors)

    return distances, connectivities.tocsr()

def _get_sparse_matrix_from_indices_distances_umap(knn_indices, knn_dists, n_obs, n_neighbors):
    rows = np.zeros((n_obs * n_neighbors), dtype=np.int64)
    cols = np.zeros((n_obs * n_neighbors), dtype=np.int64)
    vals = np.zeros((n_obs * n_neighbors), dtype=np.float64)

    for i in range(knn_indices.shape[0]):
        for j in range(n_neighbors):
            if knn_indices[i, j] == -1:
                continue  # We didn't get the full knn for i
            if knn_indices[i, j] == i:
                val = 0.0
            else:
                val = knn_dists[i, j]

            rows[i * n_neighbors + j] = i
            cols[i * n_neighbors + j] = knn_indices[i, j]
            vals[i * n_neighbors + j] = val

    result = coo_matrix((vals, (rows, cols)), shape=(n_obs, n_obs))
    result.eliminate_zeros()
    return result.tocsr()

def get_igraph_from_adjacency(adjacency):
    import igraph as ig
    sources, targets = adjacency.nonzero()
    weights = adjacency[sources, targets]
    if isinstance(weights, np.matrix):
        weights = weights.A1
    g = ig.Graph(directed=True)
    g.add_vertices(adjacency.shape[0])  # this adds adjacency.shape[0] vertices
    g.add_edges(list(zip(sources, targets)))
    try:
        g.es['weight'] = weights
    except:
        pass
    if g.vcount() != adjacency.shape[0]:
        print(
            f'The constructed graph has only {g.vcount()} nodes. '
            'Your adjacency matrix contained redundant nodes.'
        )
    return g

def louvain(adjacency, resolution = 1.0, random_state = 0):
    import louvain
    print('running Louvain clustering using the "louvain" package')
    g = get_igraph_from_adjacency(adjacency)
    part = louvain.find_partition(g, louvain.RBConfigurationVertexPartition, resolution_parameter = resolution, seed = random_state)
    groups = np.array(part.membership)
    return groups + 1 # because it starts at clust0

def leiden(adjacency, resolution = 1.0, random_state = 0):
    use_weights: bool = True
    import leidenalg
    print('running Leiden clustering')
    # convert adjacency to igraph
    g = get_igraph_from_adjacency(adjacency)
    # clustering proper
    part = leidenalg.find_partition(g, partition_type = leidenalg.RBConfigurationVertexPartition, resolution_parameter = resolution, weights = np.array(g.es['weight']).astype(np.float64), n_iterations = -1, seed = random_state)
    # store output into adata.obs
    groups = np.array(part.membership)
    return groups + 1 # because it starts at clust0

# Argument list
print('Clustering optimized from scanpy - List of arguments:')
print('1. Loom file:', loom_file)
print('2. Metadata to read:', iAnnot)
print('3. Metadata to write:', oAnnot)
print('4. Output JSON file:', oJSON)
print('5. Method to use [louvain, leiden]:', method)
print('6. Number of neighbors:', n_neighbors)
print('7. Resolution:', resolution)
print('8. Random seed:', random_state)

# Check arguments
if(len(sys.argv) < 9): error_json("Some arguments are MISSING. Stopping...", oJSON)

# Open Loom file in reading mode
f = open_with_lock(loom_file, 'r')

# Open dataset
is_nn_computed = False
if not iAnnot in f: error_json("This dataset is not present in the Loom file", oJSON)
m = f[iAnnot][:,:]
if n_neighbors > m.shape[0]:  # very small datasets [0] = nb cells
    n_neighbors = 1 + int(0.5*m.shape[0])
    warnings.append(f'n_obs too small: adjusting to `n_neighbors = {n_neighbors}`')
graphAnnot = "connectivities_{}nn_graph".format(n_neighbors) # For storing the knn graph
if "/col_graphs/" + graphAnnot in f: is_nn_computed = True    
f.close() # Close the Loom file (reading mode)

if is_nn_computed: # If already computed/stored in .loom
    print('Loading ', n_neighbors,'-nearest-neighbor graph from .loom file...')
    ds = loompy.connect(loom_file)
    _connectivities = ds.col_graphs[graphAnnot].tocsr() # Need Compressed Sparse Row format, default is COO
    ds.close()
else: # If never computed/stored in .loom for this k
    # Neighbor search (approx nearest neighbors)
    print('Computing approx.', n_neighbors,'-nearest-neighbor graph...')
    if m.shape[0] < 4096: # Nb cells
        m = pairwise_distances(m, metric=metric)
        metric = 'precomputed'
    knn_indices, knn_distances, forest = nearest_neighbors(m, n_neighbors, metric, metric_kwds = {}, angular = False, random_state = random_state)
    _distances, _connectivities = _compute_connectivities_umap(knn_indices, knn_distances, m.shape[0], n_neighbors)

    # Note: knn_indices is equal to the the graph computed in R
    print('Saving graph in .loom file...')
    f = h5py.File(loom_file, 'r+') # Just for obtaining the lock? Not sure loompy handle this
    f.close()
    ds = loompy.connect(loom_file) # I use loompy for saving/loading the graph, coz it is much easier
    ds.col_graphs[graphAnnot] = _connectivities
    ds.close()

print('Now running clustering method on the computed graph...')
if method == 'louvain':
	clusters = louvain(_connectivities, resolution = resolution, random_state = random_state)
elif method == 'leiden':
	clusters = leiden(_connectivities, resolution = resolution, random_state = random_state)      
else: error_json("This method does not exist. Use one in [louvain, leiden].", oJSON)

# Writing clustering in .loom file
f = open_with_lock(loom_file, 'r+')
if oAnnot in f: del f[oAnnot] # Delete output dataset if exists
f.create_dataset(oAnnot, data=clusters, chunks=(min(clusters.shape[0], 64)), compression="gzip", compression_opts=2, dtype='i8')
f.close()

# Preparing JSON with results
clusts, counts = np.unique(clusters, return_counts=True)
data_json = {}
data_json['time_idle'] = time_idle
if len(warnings) != 0: data_json['warnings'] = warnings
data_json['nber_clusters'] = clusts.shape[0]
data_json['metadata'] = [{
    'name': oAnnot,
    'on': 'CELL',
    'type': 'DISCRETE',
    'nber_cols' : clusters.shape[0],
    'nber_rows' : 1,
    'categories' : {str(a):int(b) for a,b in zip(clusts, counts)}
}]

# Write output.json file
with open(oJSON, 'w') as outfile:
    json.dump(data_json, outfile)
