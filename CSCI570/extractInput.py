import sys

def generateStrings(filename):
    with open(filename, 'r') as f:
        base_str1 = f.readline().strip()
        base_str2 = ''
        lines = f.readlines()
        
        base_str1_generation_indices = []
        lineIdx = 0
        for line in lines:
            if line.strip().isdigit():
                base_str1_generation_indices.append(int(line))
                lineIdx += 1
            else: 
                base_str2 = line.strip()
                lineIdx += 1
                break
        
        base_str2_generation_indices = []
        for i in range(lineIdx, len(lines)):
            base_str2_generation_indices.append(int(lines[i]))


    str1 = generator(base_str1, base_str1_generation_indices)  
    str2 = generator(base_str2, base_str2_generation_indices) 
    return str1, str2


def generator(base_string, base_string_generation_indices):
    for idx in base_string_generation_indices:
        base_string = base_string[:idx+1] + base_string + base_string[idx+1:]
    return base_string