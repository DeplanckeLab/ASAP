# TODO Handle batches (available in scanpy)

# Imports
import sys
import os
import time
import h5py
import numpy
import math
import json
import pandas
import subprocess

# Arguments
loom_file = 'MISSING'
iAnnot = 'MISSING'
oDir = 'MISSING'
min_disp = 'MISSING'
min_mean = 'MISSING'
max_mean = 'MISSING'
max_disp = 'MISSING'
n_bins = 'MISSING'
n_top_genes = 'MISSING'
is_logged = 'MISSING'
time_idle = 0
warnings = []

# Check input args
if len(sys.argv) >= 2: loom_file = sys.argv[1]
if len(sys.argv) >= 3: iAnnot = sys.argv[2]
if len(sys.argv) >= 4: oDir = sys.argv[3]

# Output JSON file
oDir = oDir.replace('\\', '/')
if not oDir.endswith('/'): oDir = oDir + '/'
oJSON = oDir + 'output.json'

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
            

if len(sys.argv) >= 5: min_disp = float(sys.argv[4]) # Default = 0.5
if len(sys.argv) >= 6: min_mean = float(sys.argv[5]) # Default = 0.0125
if len(sys.argv) >= 7: max_mean = float(sys.argv[6]) # Default = 3
if len(sys.argv) >= 8: # Default = Inf / numpy.inf
   if sys.argv[7] == 'Inf': max_disp = numpy.inf
   else: max_disp = float(sys.argv[7])
if len(sys.argv) >= 9: n_bins = int(sys.argv[8]) # Default = 20
if len(sys.argv) >= 10: # Default = None
   if sys.argv[9] == 'None': n_top_genes = None
   else: n_top_genes = int(sys.argv[9])
if len(sys.argv) >= 11: # Default = true
   if sys.argv[10] == 'true': is_logged = True
   elif sys.argv[10] == 'false': is_logged = False
   else: error_json("Last argument (11th) should be [false, true], you've put " + sys.argv[10], oJSON)

# Argument list
print('Highly Variable Genes optimized from scanpy - List of arguments:')
print('1. Loom file:', loom_file)
print('2. Metadata to read:', iAnnot)
print('3. Output folder to write new Loom file:', oDir)
print('4. Minimum dispersion:', min_disp)
print('5. Minimum mean:', min_mean)
print('6. Maximum mean:', max_mean)
print('7. Maximum dispersion:', max_disp)
print('8. Number of bins:', n_bins)
print('9. Number of top HVG:', n_top_genes)
print('10. Is the data logged?', is_logged)

# Check arguments
if(len(sys.argv) < 11): error_json("Some arguments are MISSING. Stopping...", oJSON)

# Open Loom file in reading mode
f = open_with_lock(loom_file, 'r')

# Open dataset
if not iAnnot in f: error_json("This dataset is not present in the Loom file", oJSON)
m = f[iAnnot]
if hasattr(m, 'chunks'): chunksize = m.chunks[1]# math.ceil(chunksize / m.chunks[1]) * m.chunks[1]

# Check arguments
if n_top_genes is not None and not all([min_disp is None, max_disp is None, min_mean is None, max_mean is None]):
    print('If you pass `n_top_genes`, all cutoffs are ignored.')
    warnings.append('If you pass `n_top_genes`, all cutoffs are ignored.')

# Prepare objects
mean = numpy.zeros(m.shape[0], m.dtype) # Shape = (nbgenes)
var = numpy.zeros(m.shape[0], m.dtype) # Shape = (nbgenes)

# Run chunks of all genes x n cells
totalChunks = math.ceil(m.shape[0] / chunksize)

# Calculate mean and var by chunk
for x in range(1, m.shape[0], chunksize):
    currentChunk = int((x - 1) / chunksize) + 1
    start = x - 1
    end = x + chunksize - 2
    if(end > m.shape[0] - 1): end = m.shape[0] - 1
    chunk = m[start:(end+1),:] # chunksize lines (genes) and all columns (cells)
    if is_logged: 
        if numpy.max(chunk) > 30: error_json('You said your data is logged but there are values > 30, which seems unlikely. Change the log parameter and rerun the method.', oJSON)
        X = numpy.power(2, chunk) - 1 # In Seurat, the data is unlogged (log2 / base 2)
    else: X = chunk
    for row in range(0, X.shape[0]):
        mean[x + row - 1] = numpy.mean(X[row,:])
        var[x + row - 1] = numpy.var(X[row,:])

# Close the Loom file (reading mode)
f.close()
        
# now actually compute the dispersion
mean[mean == 0] = 1e-12  # set entries equal to zero to small value
dispersion = var / mean
dispersion[dispersion == 0] = numpy.nan
dispersion = numpy.log(dispersion) # logarithmized mean as in Seurat
mean = numpy.log1p(mean)

# all of the following quantities are "per-gene" here
df = pandas.DataFrame()
df['means'] = mean
df['dispersions'] = dispersion
df['mean_bin'] = pandas.cut(df['means'], bins=n_bins)
disp_grouped = df.groupby('mean_bin')['dispersions']
disp_mean_bin = disp_grouped.mean()
disp_std_bin = disp_grouped.std(ddof=1)

# retrieve those genes that have nan std, these are the ones where
# only a single gene fell in the bin and implicitly set them to have
# a normalized disperion of 1
one_gene_per_bin = disp_std_bin.isnull()
gen_indices = numpy.where(one_gene_per_bin[df['mean_bin'].values])[0].tolist()
if len(gen_indices) > 0:
    print('Some gene indices fell into a single bin: their normalized dispersion was set to 1.\nDecreasing `n_bins` will likely avoid this effect.')
    warnings.append('Some gene indices fell into a single bin: their normalized dispersion was set to 1.\nDecreasing `n_bins` will likely avoid this effect.')

# Circumvent pandas 0.23 bug. Both sides of the assignment have dtype==float32,
# but there’s still a dtype error without “.value”.
disp_std_bin[one_gene_per_bin.values] = disp_mean_bin[one_gene_per_bin.values].values
disp_mean_bin[one_gene_per_bin.values] = 0

# actually do the normalization
df['dispersions_norm'] = (
    (
        df['dispersions'].values  # use values here as index differs
        - disp_mean_bin[df['mean_bin'].values].values
    ) / disp_std_bin[df['mean_bin'].values].values
)
dispersion_norm = df['dispersions_norm'].values.astype('float32')

# Find the HVGs
if n_top_genes is not None:
    dispersion_norm = dispersion_norm[~numpy.isnan(dispersion_norm)]
    dispersion_norm[::-1].sort()  # interestingly, np.argpartition is slightly slower
    disp_cut_off = dispersion_norm[n_top_genes-1]
    gene_subset = numpy.nan_to_num(df['dispersions_norm'].values) >= disp_cut_off
    print('The ' + str(n_top_genes) + ' top genes correspond to a normalized dispersion cutoff of ' + str(disp_cut_off))
    warnings.append('The ' + str(n_top_genes) + ' top genes correspond to a normalized dispersion cutoff of ' + str(disp_cut_off))
else:
    dispersion_norm[numpy.isnan(dispersion_norm)] = 0  # similar to Seurat
    gene_subset = numpy.logical_and.reduce((mean > min_mean, mean < max_mean, dispersion_norm > min_disp, dispersion_norm < max_disp))
df['highly_variable'] = gene_subset

# If no more genes, output error JSON
if sum(df['highly_variable']) == 0: error_json("The output matrix has no more genes. You should put less stringent thresholds.", oJSON)
print(sum(df['highly_variable']), " genes are found to be highly variable")
    
# Preparing JSON for filtering using the Java software
data_json = {}
data_json['kept_genes'] = numpy.where(df['highly_variable'])[0].tolist() # These are the indexes to filter

# Write to_keep.json file
with open(oDir + 'to_keep.json', 'w') as outfile:
    json.dump(data_json, outfile)

# Run Java for filtering the rows in the main matrix and all metadata
cmd = 'java -jar ASAP.jar -T FilterRows -loom ' + loom_file + ' -o ' + oDir + ' -m keep -row_names_file ' + oDir + 'to_keep.json'
print('RUN ' + cmd)
run_output = subprocess.run(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
print("OUT", run_output.stdout.decode('utf-8'))
print("ERR", run_output.stderr.decode('utf-8'))

# Preparing JSON for plotting results
data_json = {}
for d in ['dispersions_norm', 'dispersions']:
    # Subsetting
    disp_ = df[d]
	# Create JSON object 'plot'
    data_json['data'] = [{
        'marker':  {
            'color' : 'black',
            'opacity' : 1,
            'size' : 1,
            'name' : 'highly variable genes'
        },
        'mode': 'markers',
        'type': 'scattergl',
        'x' : df['means'][df['highly_variable'] & ~disp_.isnull()].values.tolist(), # Remove NaN
        'y' : disp_[df['highly_variable'] & ~disp_.isnull()].values.tolist() # Remove NaN
    }, {
        'marker':  {
            'color' : 'grey',
            'opacity' : 1,
            'size' : 1,
            'name' : 'other genes'
        },
        'mode': 'markers',
        'type': 'scattergl',
        'x' : df['means'][~df['highly_variable'] & ~disp_.isnull()].values.tolist(), # Remove NaN
        'y' : disp_[~df['highly_variable'] & ~disp_.isnull()].values.tolist() # Remove NaN
    }]
    data_json['layout'] = {'title': 'Highly Variable Genes', 'showlegend' : 'true', 'xaxis' : {'title' : 'mean expressions of genes'}, 'yaxis' : {'title' : 'dispersions of genes' + (' (normalized)' if d == 'dispersions_norm' else ' (not normalized)')}}
    data_json['frames'] = []
    
    # Write updated output.json file
    with open(oDir + 'plot_hvg_' + d + '.json', 'w') as outfile:
        json.dump(data_json, outfile)

# Open JSON file generated by the Java program
with open(oDir + 'output.json') as json_file:
    data_json = json.load(json_file)
    
# Add stuff to it
data_json['time_idle'] = time_idle
if len(warnings) != 0: data_json['warnings'] = warnings
data_json['plots'] = ['plot_hvg_dispersions.json', 'plot_hvg_dispersions_norm.json']

# Write updated output.json file
with open(oDir + 'output.json', 'w') as outfile:
    json.dump(data_json, outfile)