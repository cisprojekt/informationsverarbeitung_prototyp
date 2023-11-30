call
./dataprocessing.py or dataprocessing_small.py
	it will generate a distmat.csv and a plot of original data (orig_data.pdf), also stores the coordinates of the points intest_file_d.csv
then call
./glimmer.x outputtest.csv
	it will generate new coordinates in outputtest.csv
then call
./plot.py outputtest.csv
	it will generate a plot from the new generated coordinates stored in outputtest.csv
	
compare orig_data.pdf with plot.pdf
