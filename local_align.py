"""
Local Sequence Alignment with Affine Gap Penalty
Author: Aleah Hassabo

Implements local alignment for pairwise protein sequences using the
BLOSUM62 scoring matrix with affine gap penalty (opening = 11, extension = 1).
"""

import sys

def createMatrix(s1, s2):
	upper = [([0] * (len(s2)+1)) for _ in range(len(s1)+1)]
	middle = [([0] * (len(s2)+1)) for _ in range(len(s1)+1)]
	lower = [([0] * (len(s2)+1)) for _ in range(len(s1)+1)]

	blosum = [ 	
		[4, -2,  0, -2, -1, -2,  0, -2, -1, -1, -1, -1, -2, -1, -1, -1,  1,  0,  0, -3, -1, -2, -1, -8],
		[-2,  6, -3,  6,  2, -3, -1, -1, -3, -1, -4, -3,  1, -1,  0, -2,  0, -1, -3, -4, -1, -3,  2, -8],
		[0, -3,  9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -1, -2, -4, -8],
		[-2,  6, -3,  6,  2, -3, -1, -1, -3, -1, -4, -3,  1, -1,  0, -2,  0, -1, -3, -4, -1, -3,  2, -8],
		[-1,  2, -4,  2,  5, -3, -2,  0, -3,  1, -3, -2,  0, -1,  2,  0,  0, -1, -2, -3, -1, -2,  5, -8],
		[-2, -3, -2, -3, -3,  6, -3, -1,  0, -3,  0,  0, -3, -4, -3, -3, -2, -2, -1,  1, -1,  3, -3, -8],
		[0, -1, -3, -1, -2, -3,  6, -2, -4, -2, -4, -3,  0, -2, -2, -2,  0, -2, -3, -2, -1, -3, -2, -8],
		[-2, -1, -3, -1,  0, -1, -2,  8, -3, -1, -3, -2,  1, -2,  0,  0, -1, -2, -3, -2, -1,  2,  0, -8],
		[-1, -3, -1, -3, -3,  0, -4, -3,  4, -3,  2,  1, -3, -3, -3, -3, -2, -1,  3, -3, -1, -1, -3, -8],
		[-1, -1, -3, -1,  1, -3, -2, -1, -3,  5, -2, -1,  0, -1,  1,  2,  0, -1, -2, -3, -1, -2,  1, -8],
		[-1, -4, -1, -4, -3,  0, -4, -3,  2, -2,  4,  2, -3, -3, -2, -2, -2, -1,  1, -2, -1, -1, -3, -8],
		[-1, -3, -1, -3, -2,  0, -3, -2,  1, -1,  2,  5, -2, -2,  0, -1, -1, -1,  1, -1, -1, -1, -2, -8],
		[-2,  1, -3,  1,  0, -3,  0,  1, -3,  0, -3, -2,  6, -2,  0,  0,  1,  0, -3, -4, -1, -2,  0, -8],
		[-1, -1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2,  7, -1, -2, -1, -1, -2, -4, -1, -3, -1, -8],
		[-1,  0, -3,  0,  2, -3, -2,  0, -3,  1, -2,  0,  0, -1,  5,  1,  0, -1, -2, -2, -1, -1,  2, -8],
		[-1, -2, -3, -2,  0, -3, -2,  0, -3,  2, -2, -1,  0, -2,  1,  5, -1, -1, -3, -3, -1, -2,  0, -8],
		[1,  0, -1,  0,  0, -2,  0, -1, -2,  0, -2, -1,  1, -1,  0, -1,  4,  1, -2, -3, -1, -2,  0, -8],
		[0, -1, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1,  0, -1, -1, -1,  1,  5,  0, -2, -1, -2, -1, -8],
		[0, -3, -1, -3, -2, -1, -3, -3,  3, -2,  1,  1, -3, -2, -2, -3, -2,  0,  4, -3, -1, -1, -2, -8],
		[-3, -4, -2, -4, -3,  1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11, -1,  2, -3, -8],
		[-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -8],
		[-2, -3, -2, -3, -2,  3, -3,  2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1,  2, -1,  7, -2, -8],
		[-1,  2, -4,  2,  5, -3, -2,  0, -3,  1, -3, -2,  0, -1,  2,  0,  0, -1, -2, -3, -1, -2,  5, -8],
		[-8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8,  1]
	]
	return (upper, middle, lower, blosum)

def fill(upper, middle, lower, blosum):
    for i in range(len(middle)):
        middle[i][0] = 0
        upper[i][0] = 0
        lower[i][0] = 0
    for i in range(len(middle[0])):
        middle[0][i] = 0
        upper[0][i] = 0
        lower[0][i] = 0

    blosumChars = "ABCDEFGHIKLMNPQRSTVWXYZ*" 
    chartoInt = {char: i for i, char in enumerate(blosumChars)}

    for i in range(1, len(middle)):
        for j in range(1, len(middle[0])):
            s1char = sequence1[i-1]
            s2char = sequence2[j-1]
            s1int = chartoInt[s1char]
            s2int = chartoInt[s2char]
            diagonal = middle[i-1][j-1] + blosum[s1int][s2int]
            upper[i][j] = max(0, max(upper[i][j-1] - 1, middle[i][j-1] - 11))
            lower[i][j] = max(0, max(lower[i-1][j] - 1, middle[i-1][j] - 11))
            middle[i][j] = max(0, max(lower[i][j], diagonal, upper[i][j]))
			
    return middle


def alignment(s1, s2, middle, lower, upper, maxPos):
    alignment1 = ""
    alignment2 = ""
    i = maxPos[0]
    j = maxPos[1]

    while middle[i][j] != 0:
        if i > 0 and middle[i][j] == lower[i][j]:
            alignment1 = s1[i-1] + alignment1
            i -= 1
        elif j > 0 and middle[i][j] == upper[i][j]:
            alignment2 = s2[j-1] + alignment2
            j -= 1
        else:
            alignment1 = s1[i-1] + alignment1
            alignment2 = s2[j-1] + alignment2
            i -= 1
            j -= 1

    return(alignment1, alignment2)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        filename = sys.argv[1]
        outputName = "sol_q2_t" + filename[6] + ".txt"
        try:
            with open(filename, 'r') as file:
                firstline = file.readline().strip()
                line = file.readline().strip()
                
                sequence1 = ""
                sequence2 = ""

                while line[0] != ">":
                    sequence1 += line
                    line = file.readline().strip()
                
                line = file.readline().strip()

                while line:
                    sequence2 += line
                    line = file.readline().strip()

            upper, middle, lower, blosum = createMatrix(sequence1, sequence2)
            middle = fill(upper, middle, lower, blosum)
            maxScore = middle[0][0]
            maxPos = (0, 0)

            for i in range(1, len(middle)):
                for j in range(1, len(middle[0])):
                    if middle[i][j] > maxScore:
                        maxScore = middle[i][j]
                        maxPos = (i, j)

            align = alignment(sequence1, sequence2, middle, lower, upper, maxPos)

            answer = str(maxScore) + "\n" + str(align[0]) + "\n" + str(align[1])
            with open(outputName, 'w') as ofile:
                ofile.write(answer)

        except FileNotFoundError:
            print("Error: file not found")
        except IOError:
            print("Error: cannot read file")
