# Python Modules

from time import process_time
import timeit
import sys
import tracemalloc

# Custom Modules
from basic import sequence_alignment
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

def get_prefix(str1, str2):

    # Retrieving length of the strings 
    str1_length, str2_length = len(str1), len(str2)
    
    # Creating a 2D prefix DP table of size (str2_length + 1) * 2
    dp_prefix = [[0 for idx1 in range(str2_length + 1)] for idx2 in range(2)]
    
    #Intialization of Base Case
    # Intializing First Row i.e. dp_prefix[0][0...(str2_length + 1)] = idx * gap_penalty
    for idx1 in range(str2_length + 1):
        dp_prefix[0][idx1] = idx1 * gap_penalty
    
    # Filling the DP table dp_prefix
    # Row-wise 
    for idx1 in range(1, str1_length + 1):
        dp_prefix[1][0] = dp_prefix[0][0] + gap_penalty
        for idx2 in range(1, str2_length + 1):
            # Taking the minimum value of the three
            # Using Needleman-Wunsch Algorithm
            dp_prefix[1][idx2] = min(dp_prefix[0][idx2 - 1] + mismatch_penalty[str1[idx1 - 1] + str2[idx2 - 1]],
                            dp_prefix[0][idx2] + gap_penalty,
                            dp_prefix[1][idx2 - 1] + gap_penalty)
        
        # Copying row 1 to row 0
        for idx in range(0, str2_length + 1):
            dp_prefix[0][idx] = dp_prefix[1][idx]
    

    return dp_prefix[1]

def get_suffix(str1, str2):

    # Retrieving length of the strings 
    str1_length, str2_length = len(str1), len(str2)
    
    # Creating a 2D suffix DP table of size (str2_length + 1) * 2
    dp_suffix = [[0 for idx1 in range(str2_length + 1)] for idx2 in range(2)]
    
    #Intialization of Base Case
    # Intializing First Row i.e. dp_suffix[0][0...(str2_length + 1)] = idx * gap_penalty
    for idx1 in range(str2_length + 1):
        dp_suffix[0][idx1] = idx1 * gap_penalty
    
    # Filling the DP table dp_suffix
    # Row-wise 
    for idx1 in range(1, str1_length + 1):
        dp_suffix[1][0] = dp_suffix[0][0] + gap_penalty
        for idx2 in range(1, str2_length+1):
            # Taking the minimum value of the three
            # Using Needleman-Wunsch Algorithm
            dp_suffix[1][idx2] = min(dp_suffix[0][idx2 - 1] + mismatch_penalty[str1[str1_length - idx1] + str2[str2_length - idx2]],
                            dp_suffix[0][idx2] + gap_penalty,
                            dp_suffix[1][idx2 - 1] + gap_penalty)
        
        # Copying row 1 to row 0
        for idx in range(0, str2_length + 1):
            dp_suffix[0][idx] = dp_suffix[1][idx]
    
    return dp_suffix[1]



def memory_efficient_sequence_alignment(str1, str2):
    
    """
    Performs the Memory Efficient Divide & Conquer based Sequence Alignment
    alogrithm based on the Needleman-Wunsch Algorithm

    Args: 
        str1 (str): String sequence - 1 
        str2 (str): String sequence - 2

    Returns:
        list: Returns the list of three objects -> 1. Alignment of string 1 
                                                    2. Alignment of string 2
                                                    3. Similarity Cost
    """

    # Retrieving length of the strings 
    str1_length, str2_length = len(str1), len(str2)

    # Here as we use the divide & conquer method to make memory efficient algorithm
    # length of both of the strings serves as a base case for our recrusion
    # If we reach to a certain minimum length (here being 2) we call the normal sequence_alignment algorithm
    if str1_length < 2 or str2_length<2:
        # Returns the output of the sequence_alignment function
        return sequence_alignment(str1, str2)
    
    # If the string length exceeds our base case value we perform Divide operation 
    # to divide the strings into 2 halves
    else:
        # Calling get_prefix function on two halves of the big sequence Str1
        # Prefix string for (0....str1_length // 2)
        prefix = get_prefix(str1[:str1_length // 2], str2)
        # suffix string for (str1_length // 2.....) 
        suffix = get_suffix(str1[str1_length // 2:], str2)
        
        # prefix and suffix have the forward and backward sequence strings as explained in lecture
        partition_matrix = [prefix[idx] + suffix[str2_length - idx] for idx in range(str2_length + 1)]
        partition_min = min(partition_matrix)
        partition_idx = partition_matrix.index(partition_min)

        # Dumping all the temporary variables 
        # Helps in reducing Memory Usage
        prefix = []
        suffix = []
        partition = []

        # Recursive exceution of the above steps completing the Conquering Step
        
        # Recrusively perform for the First Half
        firstHalf = memory_efficient_sequence_alignment(str1[:str1_length // 2], str2[:partition_idx])
        
        # Recrusively perform for the Second Half
        secondHalf = memory_efficient_sequence_alignment(str1[str1_length // 2:], str2[partition_idx:])
        
        # Performing the simple Combine Step (i.e. Concatenating the result of the Two halfs)
        # 1. sequence 1 of firstHalf + sequence 1 of secondHalf
        # 2. sequence 2 of firstHalf + sequence 2 of secondHalf
        # 3. similarity cost of firstHalf + similarity cost of secondHalf
        return [firstHalf[r] + secondHalf[r] for r in range(3)]



if __name__ == "__main__":
   
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else: 
        filename = 'input1.txt'
    
    # problem_size = []
    # cpu_time_efficient = []
    # memory_usage_efficient = []

    for i in range(1,8):
        filename = 'input'+str(i)+'.txt'

        tracemalloc.start()
        str1, str2 = generateStrings(filename)
        # problem_size.append(len(str1)+len(str2))

        t_start=process_time()
        result =  memory_efficient_sequence_alignment(str1, str2)
        # print(result)
        t_end=process_time()
        current, peak = tracemalloc.get_traced_memory()
        # memory_usage_efficient.append(peak)
        print(f"Current memory usage is {current / 10**6}MB; Peak was {peak / 10**6}MB")
        tracemalloc.stop()

        print("Alignment of A: ", result[0])
        print("Alignment of B: ", result[1])
        print("Similarity score: ", result[2], '\n')
        print(t_end-t_start)
        # cpu_time_efficient.append(t_end-t_start)

    # print(problem_size)
    # print(cpu_time_efficient)
    # print(memory_usage_efficient)
    # print()
    # write_file(result[0], result[1], result[2], runtime, memory_used)
    # print(timeit.timeit('sequence_alignment(s1, s2)', setup=""))
    # runtime()
