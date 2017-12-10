#!/usr/bin/env python
#'''
#Nicholas Hathaway
#'''
import argparse, urllib, os

def parse_downloadFiles_args():
    parser = argparse.ArgumentParser(description="Take in a file where the first column holds the url of a file to be downloaded, will overwrite current files if they exist")
    parser.add_argument('-f', '--file', type=str, required = True)
    return parser.parse_args()

def download_files(urlsFile):
    with open(urlsFile, "r") as f:
        for line in f:
            lineSplit = line.split()
            print ("Downloading {url} to {file}".format(url = lineSplit[0], file = os.path.basename(lineSplit[0])))
            urllib.urlretrieve(lineSplit[0], os.path.basename(lineSplit[0]))

if __name__ == "__main__":
	args = parse_downloadFiles_args()
	download_files(args.file)
