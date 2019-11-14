

# FUNCTION to parse bitwise flag
def parse_bitwise(bitwise_flag):                            
    """
        The purpose of this function is to determine two things:
        1. if the read was properly mapped or not
        2. if the read was in the forward or reverse orientation

        bitwise_flag: integer flag to do bitwise comparisons to

        This function returns a list object of two booleans. The first being True if the read was mapped
        successfully and False otherwise. The second item in the list being True if the read was 
        in the reverse orientation, False if in the forward orientation.
    """ 
    mapped = True
    reverse_strand = False

    if ((bitwise_flag & 4) == 4):
        mapped = False
    
    if ((bitwise_flag & 16) == 16):
        reverse_strand = True
        
    return([mapped,reverse_strand])                           
# unit test:
assert parse_bitwise(20) == [False, True]
assert parse_bitwise(16) == [True, True]
assert parse_bitwise(4) == [False, False]


# FUNCTION that returns the 5' position and read length
def get_five_prime_pos_and_read_length_stub(CIGAR, orientation):                            
    """
        The purpose of this function is to parse the given CIGAR string based on the given
        orientation. The reason orientation is important here is because the SAM file gives us
        the leftmost position for both the forward and reverse reads, but what we really care
        about for determining PCR duplicate is the 5' position. This function also takes into
        account the inserts and delets from the reference and will return the length of the read.
        
        I will be utilizing regular expressions to pull the desired information from the CIGAR
        string.

        CIGAR: The CIGAR string is series of letters and numbers that encode information about
               the given read. I will be parsing this string and looking for soft clipping, 
               inserts and deletions.
        orientation: The orientation prameter is a string object denoting forward or reverse
                    'FW': Forward orientation
                    'RV': Reverse orientation

        This function returns a list object containing the 5' position as well as the read length:
        [five_prime_pos, read_length]
    """ 

    five_prime_pos = -1
    read_length = -1



    return([five_prime_pos, read_length])                           
# unit test:
assert len(get_five_prime_pos_and_read_length_stub('71M', 'FW')) == 2
assert get_five_prime_pos_and_read_length_stub('71M', 'FW') == [-1,-1]


# FUNCTION 
def handle_forward_stub(forward_record):                          
    """
        The purpose of this function is to make all of the function calls and analysis needed
        to return three things:
        1. The chromosome number
        2. The 5' position
        3. The length of the read

        forward_record: A list object containing a full record from a SAM file
                # forward_record[0]: header with UMIs   forward_record[1]: bitwise flag     forward_record[2]: chromosome    
                # forward_record[3]: leftmost position    forward_record[4]: mapping quality    forward_record[5]: CIGAR string

        This function returns a list object containing three things: 
        [chrom_number, five_prime_pos, read_length]
    """ 

    output = get_five_prime_pos_and_read_length_stub(forward_record[5], 'FW')

    return([forward_record[2], output[0], output[1]])                           
# unit test:
assert len(handle_forward_stub([0,0,11,0,0,0])) == 3
assert handle_forward_stub([0,0,11,0,0,0]) == [11,-1,-1]


# FUNCTION 
def handle_reverse_stub(reverse_record):                          
    """
        The purpose of this function is to make all of the function calls and analysis needed
        to return three things:
        1. The chromosome number
        2. The 5' position
        3. The length of the read

        reverse_record: A list object containing a full record from a SAM file
                # reverse_record[0]: header with UMIs   reverse_record[1]: bitwise flag     reverse_record[2]: chromosome    
                # reverse_record[3]: leftmost position    reverse_record[4]: mapping quality    reverse_record[5]: CIGAR string

        This function returns a list object containing three things: 
        [chrom_number, 5_prime_pos, read_length]
    """ 

    output = get_five_prime_pos_and_read_length_stub(reverse_record[5], 'RV')


    return([reverse_record[2], output[0], output[1]])                             
# unit test:
assert len(handle_reverse_stub([0,0,11,0,0,11])) == 3
assert handle_reverse_stub([0,0,11,0,0,0]) == [11,-1,-1]