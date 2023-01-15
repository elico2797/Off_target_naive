import time


def find_best_alignment(DNA1, DNA2, max_mismatches, max_bulges):
    # Create a matrix to store the alignments
    matrix = [[0 for x in range(len(DNA2) + 1)] for y in range(len(DNA1) + 1)]
    # Initialize the first row and column with gap penalties
    for i in range(1, len(DNA1) + 1):
        matrix[i][0] = matrix[i - 1][0] - 1
    for j in range(1, len(DNA2) + 1):
        matrix[0][j] = matrix[0][j - 1] - 1
    # Fill in the matrix by comparing the DNA strings
    for i in range(1, len(DNA1) + 1):
        for j in range(1, len(DNA2) + 1):
            # Determine the cost of a match or mismatch
            if DNA1[i - 1] == DNA2[j - 1]:
                cost = 0
            else:
                cost = -1
            # Find the minimum cost of an alignment by considering insertions, deletions, and matches/mismatches
            matrix[i][j] = max(matrix[i - 1][j - 1] + cost, matrix[i - 1][j] - 1, matrix[i][j - 1] - 1)
    # Check if the alignment is under the given constraints
    mismatches = 0
    bulges = 0
    i = len(DNA1)
    j = len(DNA2)
    while i > 0 and j > 0:
        # Check for a match
        if DNA1[i - 1] == DNA2[j - 1]:
            i -= 1
            j -= 1
        # Check for a deletion (bulge in DNA2)
        elif matrix[i][j] == matrix[i - 1][j] - 1 and bulges < max_bulges:
            bulges += 1
            i -= 1
        # Check for an insertion (bulge in DNA1)
        elif matrix[i][j] == matrix[i][j - 1] - 1 and bulges < max_bulges:
            bulges += 1
            j -= 1
        # Check for a mismatch
        else:
            mismatches += 1
            i -= 1
            j -= 1
    # Return True if the alignment is under the given constraints, False otherwise
    if mismatches <= max_mismatches and bulges <= max_bulges:
        return True, mismatches, bulges
    else:
        return False


# start counting time
start_time = time.time()

with open('C:/Users/eliav/CLionProjects/off_target_bitap/cmake-build-debug/chromosome.1.txt', 'r') as file:
    #contents = file.read(), ignore new lines
    contents = file.read().replace('\n', '')
    # print(contents)

string1 = "TCCCAAGGAACCAGAGTGGC"
max_mismatch = 5
max_bulge = 1



with open('C:/Users/eliav/CLionProjects/off_target_bitap/cmake-build-debug/Naive_output_NW.txt', 'a') as file:
    print("Naive output NW v1.0\n")
    file.write("Naive output NW v1.0\n")
    for i in range(0, len(contents)):
        for j in range(0, max_bulge + 1):
            string2 = contents[i:i + 20 + j]
            func = find_best_alignment(string1, string2, max_mismatch, max_bulge)
            if func:
                print("Inx: ", i, "Target:  ", string2, "distance: ", func[1]+func[2], "mismatches: ", func[1], "bulges: ", func[2])
                file.write("Inx: " + str(i) + " Target:  " + string2 + " distance: " + str(func[1]+func[2]) + " mismatches: " + str(func[1]) + " bulges: " + str(func[2]) + "\n")
                break
end_time = time.time()
print("time take:", end_time - start_time)
