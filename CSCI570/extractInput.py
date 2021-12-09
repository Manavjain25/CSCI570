import sys

def generateStrings(filename):
    """
    Performs cleaning of input and retrieving the base strings from the input file
    
    Args: 
        filename (str): Input filename
    
    Returns:
        list, list: Returns the generated input string from the two base strings and the given 
                    input parameters.
    """

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
    """
    Generates Sequence from the input base string and base string generations indices

    Args: 
        base_string (str): Base String sequence given in the input
        base_string_generation_indices (list): List of indices mentioned in the input file 
                                                following the base string sequence.
    
    Returns:
        str: Returns the generated string sequences.
    """
    for idx in base_string_generation_indices:
        base_string = base_string[:idx+1] + base_string + base_string[idx+1:]
    return base_string