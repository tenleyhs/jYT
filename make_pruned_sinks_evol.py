'''
cleans sinks_evol.dat of restarts and saves it as pruned_sinks_evol.dat

run this script inside each run directory
'''

import pandas as pd

df = pd.read_table('sinks_evol.dat', delim_whitespace=True)

for i, t in enumerate(df['[00]part_tag']):

	# the column headers are printed again in a restart
	if str(t) == '[00]part_tag':

		# start time (first time) after restart
		restart_time = str(df['[01]time'][i+1])

		# count how many rows back we have to delete
		# might have to increase this from 10000
		j = 0
		for j in range(100000):
			if str(df['[01]time'][i-j]) == restart_time:
				break
			j += 1
			if j==100000:
				print('Error: found restart but no matching time')

		# j+2 because we want to delete back to row j and j-1,
		# because there are two particle tags.
		for k in range(j+2):
			df = df.drop(i-k)

# this doesn't retain the nice columns with whitespace, but should still work
df.to_csv('pruned_sinks_evol.dat', sep=' ', index=False, header=False)
