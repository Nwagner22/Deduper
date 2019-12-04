
import string

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
    bitwise_flag = int(bitwise_flag)

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
def get_five_prime_pos_and_read_length(CIGAR, leftmost, orientation):                            
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
        leftmost: Integer value representing leftmost position for the read
        orientation: The orientation prameter is a string object denoting forward or reverse
                    'FW': Forward orientation
                    'RV': Reverse orientation

        This function returns a list object containing the 5' position as well as the read length:
        [five_prime_pos, read_length]
    """ 
    
    leftmost = int(leftmost)
    five_prime_pos = 0
    read_length = 0
    current_num_string = ''

    if( orientation == 'FW'): # Forward Orientation
        leading_softclipping = False

        for i in range(len(CIGAR)):
            if(CIGAR[i] not in string.ascii_uppercase):
                # Not an uppercase letter so add number to current_num_string
                current_num_string = current_num_string + CIGAR[i]
            else:
                # Is one of the valid uppercase letters so turn current_num_string into an int
                # evaluate it according to the current letter of the CIGAR string
                current_int = int(current_num_string)
                current_letter = CIGAR[i]

                if(current_letter == 'S' and (i == 1 or i == 2)):
                    # the current letter is leading soft-clipping    ADD to read length
                    # we only care about leading soft clipping on forward reads
                    read_length += current_int
                    five_prime_pos = leftmost - current_int
                    leading_softclipping = True
                elif(current_letter == 'S' and (i != 1 and i != 2)):
                    # taling soft clipping, add to read length
                    read_length += current_int

                if(current_letter == 'M'):
                    # Add matched to read length
                    read_length += current_int

                if(current_letter == 'I'):
                    # Insertion    ADD to read_length
                    read_length += current_int

                if(leading_softclipping == False):
                    five_prime_pos = leftmost

                current_num_string = ''
    else: # Reverse Orientation
        reference_length = 0

        for i in range(len(CIGAR)):
            if(CIGAR[i] not in string.ascii_uppercase):
                # Not an uppercase letter so add number to current_num_string
                current_num_string = current_num_string + CIGAR[i]
            else:
                # Is one of the valid uppercase letters so turn current_num_string into an int
                # evaluate it according to the current letter of the CIGAR string
                current_int = int(current_num_string)
                current_letter = CIGAR[i]

                if(current_letter == 'S' and (i == 1 or i == 2)):
                    # the current letter is leading soft-clipping
                    # we ignore leading soft-clipping for reverse
                    pass
                elif(current_letter == 'S' and (i != 1 and i != 2)):
                    # Add tailing soft clipping to reference length and read length
                    reference_length += current_int
                    read_length += current_int

                if(current_letter == 'M'):
                    # Add matched to reference length and read length
                    reference_length += current_int
                    read_length += current_int

                if(current_letter == 'N'):
                    # Add intron positions to reference length only
                    reference_length += current_int

                if(current_letter == 'I'):
                    # Insertion    ADD to read_length and reference_length
                    reference_length += current_int
                    read_length += current_int
                
                if(current_letter == 'D'):
                    # Deletion from reference    ADD to reference_length
                    reference_length += current_int

                current_num_string = ''

        # evaluate five_prime_pos for reverse here
        five_prime_pos = leftmost + reference_length

    return([five_prime_pos, read_length])                           
# unit test:
assert len(get_five_prime_pos_and_read_length('71M',200, 'FW')) == 2
assert get_five_prime_pos_and_read_length('23M500N24M1D11M',100, 'RV') == [659,58]
assert get_five_prime_pos_and_read_length('2S98M',102, 'FW') == [100,100]
assert get_five_prime_pos_and_read_length('100M',100, 'FW') == [100,100]



# FUNCTION 
def handle_forward(forward_record):                          
    """
        The purpose of this function is to make all of the function calls and analysis needed
        to return three things:
        1. The chromosome number
        2. The 5' position
        3. The length of the read

        forward_record: A list object containing a full record from a SAM file
                # forward_record[0]: header with UMIs   forward_record[1]: bitwise flag     forward_record[2]: chromosome    
                # forward_record[3]: leftmost position    forward_record[4]: mapping quality    forward_record[5]: CIGAR string

        This function returns a list object containing two things: 
        [(chrom_number, five_prime_pos), read_length]
    """ 

    output = get_five_prime_pos_and_read_length(forward_record[5], forward_record[3], 'FW')

    return([(forward_record[2], output[0]), output[1]])                           
# unit test:
assert len(handle_forward(['Header',0,11,100,0,'100M'])) == 2
assert handle_forward(['Header',0,11,100,0,'100M']) == [(11,100),100]


# FUNCTION 
def handle_reverse(reverse_record):                          
    """
        The purpose of this function is to make all of the function calls and analysis needed
        to return three things:
        1. The chromosome number
        2. The 5' position
        3. The length of the read

        reverse_record: A list object containing a full record from a SAM file
                # reverse_record[0]: header with UMIs   reverse_record[1]: bitwise flag     reverse_record[2]: chromosome    
                # reverse_record[3]: leftmost position    reverse_record[4]: mapping quality    reverse_record[5]: CIGAR string

        This function returns a list object containing two things: 
        [(chrom_number, 5_prime_pos), read_length]
    """ 

    output = get_five_prime_pos_and_read_length(reverse_record[5], reverse_record[3], 'RV')


    return([(reverse_record[2], output[0]), output[1]])                             
# unit test:
assert len(handle_reverse(['Header',0,11,100,0,'100M'])) == 2
assert handle_reverse(['Header',0,11,100,0,'23M500N24M1D11M']) == [(11,659),58]