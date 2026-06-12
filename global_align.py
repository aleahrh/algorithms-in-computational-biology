"""
Global Sequence Alignment with Fixed Indel Penalty
Author: Aleah Hassabo

Implements global alignment for pairwise protein sequences using the 
BLOSUM62 scoring matrix with a fixed indel penalty of 5.
"""

import sys

def createMatrix(s1, s2):
	matrix = [([0] * (len(s2)+1)) for _ in range(len(s1)+1)]

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

	return (matrix, blosum)

def fill(matrix, blosum):
	for i in range(len(matrix)):
		matrix[i][0] = i * -5
	for i in range(len(matrix[0])):
		matrix[0][i] = i * -5

	blosumChars = "ABCDEFGHIKLMNPQRSTVWXYZ*"
	chartoInt = {char: i for i, char in enumerate(blosumChars)}

	for i in range(1, len(matrix)):
		for j in range(1, len(matrix[0])):
			left = matrix[i][j-1] - 5
			above = matrix[i-1][j] - 5
			s1char = sequence1[i-1]
			s2char = sequence2[j-1]
			s1int = chartoInt[s1char]
			s2int = chartoInt[s2char]
			diagonal = matrix[i-1][j-1] + blosum[s1int][s2int]
			matrix[i][j] = max(left, above, diagonal)

	return matrix


def alignment(s1, s2, matrix):
	alignment1 = ""
	alignment2 = ""
	i = len(s1)
	j = len(s2)

	while i > 0 or j > 0:
		if i > 0 and matrix[i][j] == matrix[i-1][j] - 5:
			alignment1 = s1[i-1] + alignment1
			alignment2 = "-" + alignment2
			i -= 1
		elif j > 0 and matrix[i][j] == matrix[i][j-1] - 5:
			alignment1 = "-" + alignment1
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
		outputName = "sol_q1_t" + filename[5] + ".txt"
		try:
			with open(filename, 'r') as file:
				firstline = file.readline().strip()
				sequence1 = file.readline().strip()
				thirdline = file.readline().strip()
				sequence2 = file.readline().strip()

			matrix, blosum = createMatrix(sequence1, sequence2)
			matrix = fill(matrix, blosum)
			score = matrix[len(sequence1)][len(sequence2)]
			align = alignment(sequence1, sequence2, matrix)

			answer = str(score) + "\n" + str(align[0]) + "\n" + str(align[1]) + "\n"
			print(answer)
			with open(outputName, 'w') as ofile:
				ofile.write(answer)

		except FileNotFoundError:
			print("Error: file not found")
		except IOError:
			print("Error: cannot read file")
