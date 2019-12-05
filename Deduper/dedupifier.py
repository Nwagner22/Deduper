#!/usr/bin/env python



### DEFINE THE PROBLEM:
# PCR duplicates are a very common problem and an every day annoyance in sequencing.
# A PCR duplicate is when you have more than one identical copy of a read in your final 
# output. PCR duplicates occur during the step of library prep where we are PCR amplifying
# the fragments that have adapters ligated to them. We amplify because need to make sure we 
# have enough of the DNA molecules to fill the lane, so it is inevitable that there is going 
# to be at least some duplicate cells. Trying to reduce duplicate cells here will greatly reduce
# the amount of duplicate reads that come after the sequencing run. So no matter what, deduplication
# is highly reccomended for all PCR based sequencing output. 



import argparse
import os    # operating system
import re    # regular expressions
import deduper_functions as helpers   # import the functions I created to help with deduplication


def get_args():
    parser = argparse.ArgumentParser(description="A program to deduplicate PCR duplicates")
    parser.add_argument("-s", "--sam", help="use to specify the sam file to deduplicate", required=False, type = str)
    parser.add_argument("-o", "--output", help="Use to specify the output file name/path", default='_deduped.sam', type=str)
    parser.add_argument("-u", "--umi", help="Use to specify the umi file name/path", default='UMI96.txt', type=str)
    return parser.parse_args()

args = get_args()               # calls get_args method from above assigns the arguments to args
INPUT_SAM_FILE = args.sam          # assigning sam file path as string to global variable
UMI_FILE_PATH = args.umi
FILE_NAME = INPUT_SAM_FILE.split('/')[-1]
OUTPUT_PATH = FILE_NAME[:-4] + args.output       # building and assigning output file path as string to global variable     example: sorted_deduped.sam


SEEN_FORWARD_DICT = {}  # Dictionary to keep track of all of the UMIs seen at a specific 5' position and chromosome for forward reads  key: tuple of chromosome and 5' position    value: list of UMIs that have already been seen
SEEN_REVERSE_DICT = {}  # Dictionary to keep track of all of the UMIs seen at a specific 5' position and chromosome for reverse reads  key: tuple of chromosome and 5' position    value: list of UMIs that have already been seen
CHROM_NUMBER = 1
UMIS = []

OUT_FILE_FP = open(OUTPUT_PATH, 'w')    # filepointer associated with my output file


# Read in all of the 96 given UMIs into list data structure
with open(UMI_FILE_PATH, 'r') as umiFile:
    for line in umiFile:
        line = line.strip().split()
        UMIS.append(line[0])

### MY ALGORITHM:

# The first step in my algorithm is to make sure that the SAM file we are working 
# with is sorted. I wrote up a bash script that uses SAMtools and takes the 
# original SAM file, converts it to BAM, sorts it, and converts it back into SAM
# format. This pipeline outputs a SAM file that has been sorted by leftmost position
# which allows us to break the big problem down into smaller problems by analyzing 
# one chromosome at at time. This will drastically reduce the amount of RAM needed
# when working with large SAM files.


# The next step in the algorithm is to go through the entire sorted sam file and to check the
# following conditions:
    # 1. parse the bitwise flag to make sure the read was successfully mapped and to see if 
    #    the read was in the read was in the forward or reverse orientation
    # 2. Based on the information that is returned from the parse_bitwise function, pass off 
    #    the record to the respective function that handles that orientation
    # 3. Once the three pieces of information are returned from the respective orientation
    #    handler, the next step is to build a key out of the combination of the 5' position
    #    and the chromosome, and then append the UMI associated with the record as the value
    #    in the dictionary
    # 4A. IF the UMI had already been seen at that 5' position and it was on the same chromosome,
    #    then it is a PCR duplicate and we do not write it out the output file
    # 4B. IF the UMI had not already been seen at that 5' position and chromosome then append
    #    the UMI to the list of values and write the full record to the output file
    #




with open(INPUT_SAM_FILE, 'r') as inFile:
    duplication_counter = 0
    total_counter = 0

    for line in inFile:
        total_counter += 1
        if( line[0] == '@'):
            # write all initial headers to the output file here
            OUT_FILE_FP.write(line)
        else:
            line_list = line.strip().split()     # line[0]: header with UMIs   line[1]: bitwise flag     line[2]: chromosome    line[3]: leftmost position    line[4]: mapping quality    line[5]: CIGAR string
            temp = re.search(r'.*:([A-Z]+)',line_list[0])  # extract UMI        temp.group(1): first UMI
            umi = temp.group(1)

            # Check to see the UMIs are one of the 96 given
            if umi not in UMIS:
                # if they are not then break out of loop and skip this record
                continue

            # Each time the chromosome number changes I clear the memory from my two SEEN_DICTs
            # to prevent over-use of RAM
            if(line_list[2] != CHROM_NUMBER):
                # clear the contents of SEEN_DICT
                SEEN_FORWARD_DICT = {}
                SEEN_REVERSE_DICT = {}
                CHROM_NUMBER = line_list[2]

            
            # this is where I will do my checks and the bulk of the algorithm
            bool_list = helpers.parse_bitwise(line_list[1])    # pass bitwise flag from current record to parse_bitwise() function


            # I will now check to make sure it was mapped and also check its orientation
            if(bool_list[0] == True):  # read was successfully mapped
                if(bool_list[1] == True):  # reverse orientation
                    # Here is where I pass on the full record to the function that handles the reverse orientation records
                    output = helpers.handle_reverse(line_list)
                    tup_key = (output[1], umi)
                    if(tup_key in SEEN_REVERSE_DICT):
                        # the key (chrom_number, 5' pos) was already in the dictionary
                        # now ckeck and see if the current umi is in the list of values
                        duplication_counter += 1
                        pass
                    else:
                        # The key was not in the dictionary at all, write to file and add tup_key to dictionary and initialize value as a list with umi inside
                        OUT_FILE_FP.write(line)
                        SEEN_REVERSE_DICT[tup_key] = [1]
                        

                else:                       # forward orientation
                    # Here is where I pass on the full record to the function that handles the forward orientation records
                    output = helpers.handle_forward(line_list)
                    tup_key = (output[1], umi)
                    if(tup_key in SEEN_FORWARD_DICT):
                        # the key (chrom_number, 5' pos) was already in the dictionary
                        # now ckeck and see if the current umi is in the list of values
                        duplication_counter += 1
                        pass
                    else:
                        # The key was not in the dictionary at all, write to file and add tup_key to dictionary and initialize value as a list with umi inside
                        OUT_FILE_FP.write(line)
                        SEEN_FORWARD_DICT[tup_key] = [umi]
                        
            else:   # The read did not map, DO NOTHING
                pass


            
            #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            #
            # DEPENDS ON FORWARD OR REVERSE ORIENTATION BUT THE PROCESS BELOW WILL BE THE SAME:
            #
            # now check if (5' position, chromosome) combined as key are in dictionary, and
            # if already in dictionary check to see if the current UMI is in the list of UMIs
            # that have already been seen, 
            #
            # IF the UMI has already been added to the dictionary then it is a duplicate
            #
            # IF the UMI was not already one of the values, add it to the list of UMIs and write
            # the whole line to the file  OUT_FILE_FP.write(line)

print('Removed ' + str(duplication_counter) + ' duplicates' )

OUT_FILE_FP.close()




                        # if(umi in SEEN_FORWARD_DICT[tup_key]):
                        #     # if the umi is already in the list of values for this key then the read is a duplicate and nothing happens with this record
                        #     pass
                        # else:
                        #     # add the umi to list of values and write the record out to the file
                        #     SEEN_FORWARD_DICT[tup_key].append(umi)
                        #     OUT_FILE_FP.write(line)