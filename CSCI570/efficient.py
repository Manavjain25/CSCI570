from time import process_time
import timeit
import sys
from basic import sequence_alignment
from extractInput import generateStrings
import tracemalloc

# Global Variables
mismatch_penalty = {'AA': 0, 'AC':110, 'CA': 110, 
                    'AG': 48, 'GA': 48, 'CC':0, 
                    'AT': 94, 'TA': 94, 'CG':118, 
                    'GC': 118, 'GG': 0, 'TG':110, 
                    'GT':110, 'TT':0, 'CT': 48, 
                    'TC':48}

gap_penalty = 30

def prefix(str1, str2):
    n, m = len(str1), len(str2)
    # mat = [[0 for i in range(2)] for j in range(n + 1)]
    dp_prefix = [[0 for idx1 in range(m + 1)] for idx2 in range(2)]
    #initialize the base cases 
    for idx1 in range(m + 1):
        dp_prefix[0][idx1] = idx1 * gap_penalty
    
    # update the matrix in row order 
    for idx1 in range(1, n + 1):
        dp_prefix[1][0] = dp_prefix[0][0] + gap_penalty
        for idx2 in range(1, m + 1):
            dp_prefix[1][idx2] = min(dp_prefix[0][idx2 - 1] + mismatch_penalty[str1[idx1 - 1] + str2[idx2 - 1]],
                            dp_prefix[0][idx2] + gap_penalty,
                            dp_prefix[1][idx2 - 1] + gap_penalty)
        
        for idx in range(0, m + 1):
            dp_prefix[0][idx] = dp_prefix[1][idx]
    
    return dp_prefix[1]

def suffix(str1, str2):
    n, m = len(str1), len(str2)
    # mat = [[0 for i in range(2)] for j in range(n + 1)]
    dp_suffix = [[0 for idx1 in range(m + 1)] for idx2 in range(2)]
    
    #initialize the base cases 
    for idx1 in range(m + 1):
        dp_suffix[0][idx1] = idx1 * gap_penalty
    
    #update the matrix in row order 
    for idx1 in range(1, n + 1):
        dp_suffix[1][0] = dp_suffix[0][0] + gap_penalty
        for idx2 in range(1, m+1):
            dp_suffix[1][idx2] = min(dp_suffix[0][idx2 - 1] + mismatch_penalty[str1[n - idx1] + str2[m - idx2]],
                            dp_suffix[0][idx2] + gap_penalty,
                            dp_suffix[1][idx2 - 1] + gap_penalty)
        
        for idx in range(0, m + 1):
            dp_suffix[0][idx] = dp_suffix[1][idx]
    
    return dp_suffix[1]



def space_efficient_alignment(str1, str2):
    # This is the main space_efficient_alignment routine.
    n, m = len(str1), len(str2)
    if n<2 or m<2:
        # In this case we just use the N-W algorithm.
        # return nw(str1, str2, alphEnum)
        return sequence_alignment(str1, str2)

    else:
        # Make partitions, call subroutines.
        # F, B = forwards(str1[:n//2], str2, alphEnum), backwards(str1[n//2:], str2, alphEnum)
        # F, B = forwards(str1[:n//2], str2, alphEnum), backwards(str1[n//2:], str2, alphEnum)
        F, B = prefix(str1[:n//2], str2), suffix(str1[n//2:], str2)
        # assert(E == F)
        # assert(A == B)
        partition = [F[j] + B[m-j] for j in range(m+1)]
        cut = partition.index(min(partition))
        # Clear all memory now, so that we don't store data during recursive calls.
        F, B, partition = [], [], []
        # Now make recursive calls.
        callLeft = space_efficient_alignment(str1[:n//2], str2[:cut])
        callRight = space_efficient_alignment(str1[n//2:], str2[cut:])
        # Now return result in format: [1st alignment, 2nd alignment, similarity]
        return [callLeft[r] + callRight[r] for r in range(3)]



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
        result =  space_efficient_alignment(str1, str2)
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
