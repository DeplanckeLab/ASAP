from ZIFA import ZIFA,block_ZIFA
#from pylab import *
import numpy as np
import random
import json
import sys

def runZIFA():
	random.seed(42)
	np.random.seed(42)
	
	print 'Number of arguments:', len(sys.argv), 'arguments.'
	print 'Argument List:', str(sys.argv)

	inputfilename = sys.argv[1]
	outputfolder = sys.argv[2]

	input_alldimensions = []
	with open(inputfilename, 'r') as infile:
		input = infile.readlines()
		cell_names = input.pop(0).rstrip('\r\n').split('\t')
		cell_names.pop(0)
		for line in input:
			line = line.rstrip('\r\n')
			linearray = []
			l = line.split('\t')
			l.pop(0)
			for it in l:
				number = float(it)
				if number < 0.0000001: 
					number = float(0)
				linearray.append(number)
			input_alldimensions.append(linearray)
				
	alldim = np.asarray(input_alldimensions)
	alldim = alldim.transpose()
	
	try:
		with open(outputfolder + "/log.zifa.txt", 'w') as f:
			sys.stdout = f
			Zhat, params = block_ZIFA.fitModel(alldim, min(5,len(cell_names)))
		sys.stdout = sys.__stdout__
	except Exception as err:
		f = open(outputfolder + "/log.zifa.txt", 'r')
		output_json = {}
		errorMsg = str(err[0])
		if errorMsg.startswith("Your input matrix contains no zeros"):
			output_json['displayed_error'] = "Zifa is not converging. This can be due to an input matrix which contains no zeros. Zifa input should be log read counts. You can try another filtering/normalization but this may not solve the issue for this dataset."
		else:
			output_json['displayed_error'] = errorMsg
		output_json['original_error'] = f.read()
		with open(outputfolder + "/output.json", 'w') as outfile:
			json.dump(output_json, outfile)
		raise

	output_json = {}
	output_json['PC1'] = []
	output_json['PC2'] = []
	output_json['PC3'] = []
	output_json['PC4'] = []
	output_json['PC5'] = []
	output_json['text'] = []
	
	i = 0
	for it in Zhat:
		if len(it)>=1: output_json['PC1'].append(it[0])
		if len(it)>=2: output_json['PC2'].append(it[1])
		if len(it)>=3: output_json['PC3'].append(it[2])
		if len(it)>=4: output_json['PC4'].append(it[3])
		if len(it)>=5: output_json['PC5'].append(it[4])
		output_json['text'].append(cell_names[i])
		i += 1
	
	if len(output_json['PC1']) == 0: del(output_json['PC1'])
	if len(output_json['PC2']) == 0: del(output_json['PC2'])
	if len(output_json['PC3']) == 0: del(output_json['PC3'])
	if len(output_json['PC4']) == 0: del(output_json['PC4'])
	if len(output_json['PC5']) == 0: del(output_json['PC5'])

	with open(outputfolder + "/output.json", 'w') as outfile:
		json.dump(output_json, outfile)

if __name__ == '__main__':
	runZIFA()
