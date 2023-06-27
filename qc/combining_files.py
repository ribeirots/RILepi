# This script merges multiple inputs that have a filename pattern in common, zips the output and remove the original files

import os, glob, zipfile, fnmatch, argparse, re

parser = argparse.ArgumentParser()
parser.add_argument('-p','--pattern', type=str, help='Filename pattern of files to be merged.', required=True)
parser.add_argument('-r','--range', type=str, default="0", help='If files and in an increasing number right before the extension, and you would like the merge to follow that order, supply first and last number of the sequence separated by comma, e.g. "1,10". Still include * in the pattern right before the extension and the * will be replaced with the number. The script will still zip all files matching the pattern regardless of numbering.', required=False)
parser.add_argument('-o','--outfile', type=str, help='Merged file output name.', required=True)
parser.add_argument('-header','--header', type=str, default="y", help='(y/n) - type "n" if the files to be merged do NOT have a header. If no information is given, it assumes that there is a header.', required=False)
parser.add_argument('-nc','--nocheck', type=int, help='If any integer greater than 0 is supplied, it will not check if the output files already exist and will overwrite any file with that name that is present. If nothing is supplied, performs the check.', required=False)
parser.add_argument('-nv','--noverb', type=int, help='If any integer greater than 0 is supplied, it will not output text to the screen. If nothing is supplied, outputs text.', required=False)
args = parser.parse_args()

# Counts how many files match the pattern - used to check whether is there any file to be merged and to check the zip files
def count_files(file_pattern):
    file_count = len(glob.glob(file_pattern))
    return file_count

# the merge file function
def combine_files(file_pattern, output_file, seq, header):
    if header == "y" or header == "Y":
        first_file = 1
    elif header == "n" or header == "N":
        first_file = 0
    else:
        print("-header argument used incorrect. Set it to 'y' if there is a header or 'n' is there is no header in the files to be merged.\nTerminating script.")
        quit()

    with open(output_file, "w") as outfile:
        if len(seq) == 2:
            for i in range(seq[0],seq[1]):
                with open(pattern[:-5]+str(i)+pattern[-4:]) as infile:
                    if first_file:
                        first_file = 0
                    else:
                        if header == "y" or header == "Y":
                            next(infile)
                    outfile.write(infile.read())
        else:
            for file_name in glob.glob(file_pattern):
                with open(file_name, "r") as infile:
                    if first_file:
                        first_file = 0
                    else:
                        if header == "y" or header == "Y":
                            next(infile)
                    outfile.write(infile.read())

# create zip file and remove original files
def compress_directory(directory_path, zip_path, file_pattern):
    with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for root, _, files in os.walk(directory_path):
            for file in files:
                file_path = os.path.join(root, file)
                if fnmatch.fnmatch(file, file_pattern):
                    zipf.write(file_path, os.path.relpath(file_path, directory_path))
                    os.remove(file_path)

# check if the zip file already exists, it will not overwrite it and create a new one with different name
# this ensures that the original files won't be lost in case they were deleted and only existed in the zip form
# output zip file starts with 'zip' instead of just mimicking the merged output filename to create a different pattern in case merged has to be deleted by the user it minimizes the chance of deleting the zip file by mistake
if os.path.exists("zip_"+args.outfile[:-4]+".zip"):
    zip_count = count_files("zip_"+args.outfile[:-4]+".zi*")
    zip_output = "zip_"+args.outfile[:-4]+".zip"+str(zip_count)
    if not args.noverb:
        print("zip_"+args.outfile[:-4]+".zip already exists!")
        print("Saving file as *.zip"+str(zip_count))
else:
    zip_output = "zip_"+args.outfile[:-4]+".zip"


pattern = args.pattern
directory_to_compress = "./"
nfiles = count_files(args.pattern)

# if no file match the pattern, the script doesn't run.
# this is a redundant safety measure to renaming the zip files to ensure original files aren't lost
if nfiles == 0:
    print("No files that match this pattern were found.\nTerminating script.")
    quit()

if not args.nocheck:
    if os.path.exists(args.outfile):
        print(args.outfile+' already exists.')
        resp = input("Would you like to replace it? (y/n): ")
        if resp.lower() != "y" and resp.lower() != "yes":
            print("Script aborted.\nQuitting.")
            quit()

n_order = args.range
if n_order != "0":
    n_order = re.split(",",n_order)
    if len(n_order) != 2:
        print("Sequence of numbers supplied incorrectly. Type it separated by 1 comma.\nTerminating script.")
        quit()
    else:
        n_order = list(map(int,n_order))
        n_order[1] += 1


# combine files and then compress then
combine_files(pattern, args.outfile, seq=n_order, header=args.header)
compress_directory(directory_to_compress, zip_output, pattern)

