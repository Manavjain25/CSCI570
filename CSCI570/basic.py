import sys
import tracemalloc
from extractInput import generateStrings
from time import process_time
# Global Variables
mismatch_penalty = {'AA': 0, 'AC':110, 'CA': 110, 
                    'AG': 48, 'GA': 48, 'CC':0, 
                    'AT': 94, 'TA': 94, 'CG':118, 
                    'GC': 118, 'GG': 0, 'TG':110, 
                    'GT':110, 'TT':0, 'CT': 48, 
                    'TC':48}

gap_penalty = 30

def build_aligned_sequences(m, n, str1, str2, dp_table):
    idx1 = m
    idx2 = n
    aligned_str1, aligned_str2 = '', ''
    while idx1 and idx2:
    
        if dp_table[idx1][idx2] == gap_penalty + dp_table[idx1 - 1][idx2]:
            aligned_str1 = str1[idx1 - 1] + aligned_str1
            aligned_str2 = '_' + aligned_str2
            idx1 -= 1
        elif dp_table[idx1][idx2] == mismatch_penalty[str1[idx1 - 1] + str2[idx2-1]] + dp_table[idx1 - 1][idx2 - 1]:
            aligned_str1 = str1[idx1 - 1] + aligned_str1
            aligned_str2 = str2[idx2 - 1] + aligned_str2
            idx1 -= 1
            idx2 -= 1
        elif dp_table[idx1][idx2] == gap_penalty + dp_table[idx1][idx2 - 1]:
            aligned_str1 = '_' + aligned_str1
            aligned_str2 = str2[idx2 - 1] + aligned_str2
            idx2 -= 1
    
    while idx2:
        aligned_str1 = '_' + aligned_str1
        aligned_str2 = str2[idx2 - 1] + aligned_str2
        idx2 -= 1
    while idx1:
        aligned_str1 = str1[idx1 - 1] + aligned_str1
        aligned_str2 = '_' + aligned_str2
        idx1 -= 1

    
    return aligned_str1, aligned_str2
    

def sequence_alignment(str1, str2):
    m = len(str1)
    n = len(str2)

    dp_table = [[0 for i in range(n + 1)] for j in range(m + 1)]

    #initialization Str[idx1, 0]= idx1 * gap_penalty, 
    for idx1 in range(m + 1):
        dp_table[idx1][0] = idx1 * gap_penalty
    
    #initialization Str[0, idx2]= idx2 * gap_penalty, 
    for idx2 in range(n + 1):
        dp_table[0][idx2] = idx2 * gap_penalty
    
    for idx2 in range(1, n + 1):
        for idx1 in range(1, m + 1):
            dp_table[idx1][idx2] = min(mismatch_penalty[str1[idx1 - 1] + str2[idx2 - 1]] + dp_table[idx1 - 1][idx2 - 1], gap_penalty + 
                                    dp_table[idx1 - 1][idx2], gap_penalty + dp_table[idx1][idx2 - 1])

    aligned_str1, aligned_str2 = build_aligned_sequences(m, n, str1, str2, dp_table)
    
    # alignment,_ = get_aligned(str1, str2, dp_table, gap_penalty, mismatch_penalty)
    # return dp_table[m][n], alignment
    return [aligned_str1, aligned_str2, dp_table[m][n]]

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
