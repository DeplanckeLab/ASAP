# Imports
import sys
import h5py
import numpy
import json
from MulticoreTSNE import MulticoreTSNE as TSNE

# Arguments
loom_file = 'MISSING'
iAnnot = 'MISSING'
oAnnot = 'MISSING'
oJSON = 'MISSING'
n_comps = 'MISSING'
perplexity = 'MISSING'
early_exaggeration = 'MISSING'
learning_rate = 'MISSING'
random_state = 'MISSING'
n_threads = 'MISSING'
time_idle = 0

if len(sys.argv) >= 2: loom_file = sys.argv[1]
if len(sys.argv) >= 3: iAnnot = sys.argv[2]
if len(sys.argv) >= 4: oAnnot = sys.argv[3]
if len(sys.argv) >= 5: oJSON = sys.argv[4]
if len(sys.argv) >= 6: n_comps = int(sys.argv[5])
if len(sys.argv) >= 7: perplexity = float(sys.argv[6])
if len(sys.argv) >= 8: early_exaggeration = float(sys.argv[7])
if len(sys.argv) >= 9: learning_rate = float(sys.argv[8])
if len(sys.argv) >= 10: random_state = int(sys.argv[9])
if len(sys.argv) >= 11: n_threads = int(sys.argv[10])

print('t-SNE optimized from scanpy - List of arguments:')
print('1. Loom file:', loom_file)
print('2. Metadata to read:', iAnnot)
print('3. Metadata to write:', oAnnot)
print('4. Output JSON file:', oJSON)
print('5. Number of t-SNE dimensions to compute:', n_comps)
print('6. Perplexity:', perplexity)
print('7. Early exaggeration:', early_exaggeration)
print('8. Learning rate:', learning_rate)
print('9. Random state:', random_state)
print('10. Nb of threads to use for parallelization:', n_threads)

# Functions
# Handling errors
def error_json(message, jsonfile):
    print(message, file=sys.stderr)
    if(jsonfile != 'MISSING'):
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
            
# Check arguments
if(len(sys.argv) < 11): error_json("Some arguments are MISSING. Stopping...", oJSON)

# Open Loom file in reading mode
f = open_with_lock(loom_file, 'r')

# Open dataset
if not iAnnot in f: error_json("This dataset is not present in the Loom file", oJSON)
m = f[iAnnot][:,:] # full object is stored in RAM

# Close Loom
f.close()

# Prepare t-SNE object
tsne_ = TSNE(n_components=n_comps, n_jobs=n_threads, perplexity=perplexity, random_state=random_state, early_exaggeration=early_exaggeration, learning_rate=learning_rate)  

# Run t-SNE
X_tsne = tsne_.fit_transform(m)

# Open in writing to write output
f = open_with_lock(loom_file, 'r+')
if oAnnot in f: del f[oAnnot] # Delete output dataset if exists
f.create_dataset(oAnnot, data=X_tsne, chunks=(min(X_tsne.shape[0], 32), min(X_tsne.shape[1], 32)), compression="gzip", compression_opts=2)
f.close()

# Prepare output.json
data_json = {}
data_json['time_idle'] = time_idle
data_json['metadata'] = [{
    'name': oAnnot,
    'on': 'CELL',
    'type': 'NUMERIC',
    'nber_cols' : X_tsne.shape[0],
    'nber_rows' : X_tsne.shape[1]
}]

# Write output.json file
with open(oJSON, 'w') as outfile:
    json.dump(data_json, outfile)
