import pandas as pd
import urllib.request
import os
import argparse
import re
# Parsing arguments
parser = argparse.ArgumentParser(description='Downloading metagenomes by term.')
parser.add_argument('--term', '-t', required=True, help='Term, string value')
parser.add_argument('--input', '-i', required=True, help='Input .csv file')
parser.add_argument('--output_dir', '-d', required=True, help='The path to the directory where you want to save the files, without "/"')
d = vars(parser.parse_args())
term, input_file, path_to_dir = d['term'], d['input'], d['output_dir']
term = term.replace("_", " ")
# Open csv file like pandas df
df = pd.read_csv(input_file)
df = df[df["Name"]== term]
#Download
for i in df["URL_to_file"]:
     link = i
     file = os.path.basename(link) + '_genomic.fna.gz'
     urllib.request.urlretrieve(i, path_to_dir+"/"+file)


