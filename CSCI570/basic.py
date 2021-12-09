# Python Modules
import sys
import tracemalloc
from time import process_time

# Custom Modules
from extractInput import generateStrings

# Global Variables
# Mismatch Penalty matrix created using mismatch penalty values given
mismatch_penalty = {'AA': 0, 'AC':110, 'CA': 110, 
                    'AG': 48, 'GA': 48, 'CC':0, 
                    'AT': 94, 'TA': 94, 'CG':118, 
                    'GC': 118, 'GG': 0, 'TG':110, 
                    'GT':110, 'TT':0, 'CT': 48, 
                    'TC':48}

# Gap_Penalty
# default value = 30
gap_penalty = 30


def build_aligned_sequences(str1, str2, dp_table):
    """
    Generates the string by backtracking the dp_table using the gap_penalties assigned to each character in the sequence
    
    Args:
        str1 (str): String sequence - 1
        str2 (str): String sequence - 2
        dp_table (list): 2D Array containing all the similarity cost of sequences

    Returns:
        list, list: Returns the aligned sequences of both the strings
    """

    # Retrieving length of the strings 
    str1_length = len(str1)
    str2_length = len(str2)

    idx1 = str1_length
    idx2 = str2_length

    # Intializing empty strings to store aligned strings.
    aligned_str1, aligned_str2 = '', ''

    # We start traversing the dp_table from bottom-right corner
    # Thus, creating the aligned sequence of strings using BackTracking
    while idx1 and idx2:
        # If at dp_table[idx1][idx2] we had incur a gap_penalty from [idx1 -1][idx2] this means
        # we have to create a gap '_' in the aligned_str2
        # Thus, we append the str1[idx1 - 1] at the start of aligned_str1 
        # To make sure we avoid getting reversed string
        if dp_table[idx1][idx2] == gap_penalty + dp_table[idx1 - 1][idx2]:
            aligned_str1 = str1[idx1 - 1] + aligned_str1
            aligned_str2 = '_' + aligned_str2
            idx1 -= 1
        
        # If we incur the mismatch penalty at (idx1, idx2) we add both the 
        # characters into the start of the aligned_str{1,2}
        # Following the decrement of idx{1,2}
        elif dp_table[idx1][idx2] == mismatch_penalty[str1[idx1 - 1] + str2[idx2-1]] + dp_table[idx1 - 1][idx2 - 1]:
            aligned_str1 = str1[idx1 - 1] + aligned_str1
            aligned_str2 = str2[idx2 - 1] + aligned_str2
            idx1 -= 1
            idx2 -= 1
        
        # If at dp_table[idx1][idx2] we had incur a gap_penalty from [idx1][idx2 - 1] this means
        # we have to create a gap '_' in the aligned_str1
        # Thus, we append the str2[idx2 - 1] at the start of aligned_str2
        # To make sure we avoid getting reversed string
        elif dp_table[idx1][idx2] == gap_penalty + dp_table[idx1][idx2 - 1]:
            aligned_str1 = '_' + aligned_str1
            aligned_str2 = str2[idx2 - 1] + aligned_str2
            idx2 -= 1
    
    # We add the remaining characters of the sequence 2 if str2_length > str1_length
    while idx2:
        aligned_str1 = '_' + aligned_str1
        aligned_str2 = str2[idx2 - 1] + aligned_str2
        idx2 -= 1
    
    # We add the remaining characters of the sequence 1 if str1_length > str2_length
    while idx1:
        aligned_str1 = str1[idx1 - 1] + aligned_str1
        aligned_str2 = '_' + aligned_str2
        idx1 -= 1

    
    return aligned_str1, aligned_str2
    

def sequence_alignment(str1, str2):
    """
    Performs the Basic Sequence Alignment based on the Needleman-Wunsch Algorithm
    which uses Dynamic Programming approach

    Args: 
        str1 (str): String sequence - 1 
        str2 (str): String sequence - 2

    Returns:
        list: Returns the list of three objects -> 1. Alignment of string 1 
                                                    2. Alignment of string 2
                                                    3. Similarity Cost
    """
    
    # Retrieving length of the strings 
    str1_length = len(str1)
    str2_length = len(str2)

    # Creating a 2D (dp_table) DP table of size (str2_length + 1) * (str1_length + 1)
    dp_table = [[0 for i in range(str2_length + 1)] for j in range(str1_length + 1)]

     #Intialization of Base Case
    # Intializing First Column i.e. dp_table[0...(str1_length + 1)][0] = idx * gap_penalty
    for idx1 in range(str1_length + 1):
        dp_table[idx1][0] = idx1 * gap_penalty
    
     #Intialization of Base Case
    # Intializing First Row i.e. dp_table[0][0...(str2_length + 1)] = idx * gap_penalty
    for idx2 in range(str2_length + 1):
        dp_table[0][idx2] = idx2 * gap_penalty
    
    # Filling the DP Table
    for idx2 in range(1, str2_length + 1):
        for idx1 in range(1, str1_length + 1):
            # Taking the minimum value of the three
            # Using Needleman-Wunsch Algorithm
            dp_table[idx1][idx2] = min(mismatch_penalty[str1[idx1 - 1] + str2[idx2 - 1]] + dp_table[idx1 - 1][idx2 - 1], gap_penalty + 
                                    dp_table[idx1 - 1][idx2], 
                                    gap_penalty + dp_table[idx1][idx2 - 1])

    # Building the aligned sequence of str1 and str2 from similarity values stored 
    # in the dp_table (Backtracking the dp_table to retrieve the aligned sequence)
    aligned_str1, aligned_str2 = build_aligned_sequences(str1, str2, dp_table)
    
    # Returning three values
    # 1. Aligned Sequence of str1
    # 2. Aligned Sequence of str2
    # 3. Similarity cost stored at the bottom-right cell of the dp_table.
    return [aligned_str1, aligned_str2, dp_table[str1_length][str2_length]]

if __name__ == "__main__":
   
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else: 
        filename = 'input1.txt'

    # cpu_time_basic = []
    # memory_usage_basic = []

    for i in range(1,8):
        filename = 'input'+str(i)+'.txt'

        tracemalloc.start()
        str1, str2 = generateStrings(filename)
        t_start=process_time()
        result =  sequence_alignment(str1, str2)
        # print(result)

        t_end=process_time()
        current, peak = tracemalloc.get_traced_memory()
        # memory_usage_basic.append(peak)
        print(f"Current memory usage is {current / 10**6}MB; Peak was {peak / 10**6}MB")
        tracemalloc.stop()

        print("Alignment of A: ", result[0])
        print("Alignment of B: ", result[1])
        print("Similarity score: ", result[2], '\n')
        # cpu_time_basic.append(t_end-t_start)
        print(t_end-t_start)
    # print(time.process_time())
    # print(cpu_time_basic)
    # print(memory_usage_basic)
    # print("Alignment of A: ", result[0])
    # print("Alignment of B: ", result[1])
    # print("Similarity score: ", result[2], '\n')
    # write_file(result[0], result[1], result[2], runtime, memory_used)
    # print(timeit.timeit('SequenceAlignment(s1, s2)', setup=""))
    # runtime()
