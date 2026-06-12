"""
Z-Algorithm for Exact Pattern Matching
Author: Aleah Hassabo

Implements the Z-algorithm for exact pattern matching on DNA sequences.
Given a text t and pattern p, finds all occurrences of p in t and reports
the number of character comparisons, matches, and mismatches.
"""

import sys


def compute_z(input):
    """
    Computes the Z-array for a given string.
    Z[i] = length of the longest substring starting at position i
    that matches a prefix of the string.
    """
    num = len(input)
    z = [0] * num

    compare = 0
    match = 0
    mismatch = 0

    l = 0
    r = 0

    for k in range(1, num):

        # Case 1: k is outside the current Z-box
        if k > r:
            i = 0
            while (k + i < num) and (input[i] == input[k + i]):
                compare += 1
                match += 1
                i += 1

            if (k + i) < num:
                compare += 1
                mismatch += 1

            z[k] = i

            if i > 0:
                l = k
                r = k + i - 1

        # Case 2: k is inside the current Z-box
        else:
            alpha = k - l
            beta = r - k + 1

            # Case 2a: Z[alpha] < beta, no extension needed
            if z[alpha] < beta:
                z[k] = z[alpha]

            # Case 2b: Z[alpha] >= beta, extend from r+1
            elif z[alpha] >= beta:
                i = beta
                while (k + i < num) and (input[i] == input[k + i]):
                    compare += 1
                    match += 1
                    i += 1

                if (k + i) < num:
                    compare += 1
                    mismatch += 1

                z[k] = i

                if i > 0:
                    l = k
                    r = k + i - 1

    return z, compare, match, mismatch


if __name__ == "__main__":
    if len(sys.argv) > 1:
        filename = sys.argv[1]
        try:
            with open(filename, 'r') as file:
                t = file.readline().strip()
                p = file.readline().strip()

                # Concatenate pattern, delimiter, and text
                concatenated = p + '$' + t

                z, compare, match, mismatch = compute_z(concatenated)

                # Find all pattern occurrences (1-based indexing)
                pattern_len = len(p)
                for i in range(len(z)):
                    if z[i] == pattern_len:
                        print(i - pattern_len)

                print(f"Number of comparisons:  {compare}")
                print(f"Number of matches:      {match}")
                print(f"Number of mismatches:   {mismatch}")

        except FileNotFoundError:
            print("Error: file not found")
        except IOError:
            print("Error: cannot read file")
