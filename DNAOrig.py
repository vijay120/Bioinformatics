
def findLTPatterns(dnaFile, patternLength, freq, windowLength):
    resultSet = set()
    for line in dnaFile:
        dnaText = line
        mappings = comprehensiveFreq(dnaText[0:windowLength], patternLength)
        for clump in mappings.keys():
            sampleFreq = mappings[clump]
            if (sampleFreq >= freq):
                resultSet.add(clump)

        for i in range(1, len(dnaText) - windowLength + 1):
            word = dnaText[i:i+windowLength]
            firstWord = word[0:patternLength]
            lastWord = word[windowLength-patternLength: windowLength]
            if (firstWord in mappings):
                mappings[firstWord] = mappings[firstWord] - 1

            if (lastWord in mappings):
                mappings[lastWord] = mappings[lastWord] + 1
            else:
                mappings[lastWord] = 1

            for clump in mappings.keys():
                sampleFreq = mappings[clump]
                if (sampleFreq >= freq):
                    resultSet.add(clump)

    return list(resultSet)

def patternCount(dnaText, pattern):
    count = 0
    for i in range(0, len(dnaText) - len(pattern) + 1):
        word = dnaText[i:i+len(pattern)]
        if (word == pattern):
            count = count + 1
    return count

def freqWordProblem(text, k):
    countWords = []
    for i in range(0, len(text) - k + 1):
        pattern = text[i:i+k]
        countWords.append(patternCount(text, pattern))

    maxCount = 0
    indexes = []
    for j in range(0, len(countWords)):
        count = countWords[j]
        if (count == maxCount):
            indexes.append(j)
        elif (count > maxCount):
            indexes = [j]
            maxCount = count

    result = set()
    for index in indexes:
        result.add(text[index:index+k])

    return list(result)


def comprehensiveFreq(text, k):
    mapping = {}
    for i in range(0, len(text) - k + 1):
        pattern = text[i:i+k]
        mapping[pattern] = patternCount(text, pattern)
    return mapping

mapDNA = {
    "A": "T",
    "G": "C",
    "T": "A",
    "C": "G"
}
def complimentDNA(text):
    result = ""
    for letter in text:
        result = result + mapDNA[letter]
    return result[::-1]

def patternFind(text, pattern):
    index = []
    for i in range(0,len(text)-len(pattern) + 1):
        word = text[i:i+len(pattern)]
        if word == pattern:
            index.append(i)
    return index

import math
import numpy as np
from numpy.linalg import inv

def PatternToNumber(dnaWord):
    mapping = ['A', 'C', 'G', 'T']
    result = 0
    for i in range(0, len(dnaWord)):
        exponent = len(dnaWord) - i - 1
        result = result + (math.pow(4, exponent) * mapping.index(dnaWord[i]))
    return int(result)


def NumberToPattern(value, length):
    operator = np.ones((1, length))
    operand = np.zeros((1, 1))
    operand[0] = value
    for i in range(0, length):
        operator[0][i] = math.pow(4, length-1-i)
    print(np.dot(np.transpose(operator), operator))
    return np.dot(np.dot(inv(np.dot(np.transpose(operator), operator)), np.transpose(operator)), operand)

print(NumberToPattern(5437, 8))

# print(PatternToNumber("ATGCAA"))

dnaText = ""
windowLength = 0
patternLength = 0
freq = 0
counter = 0

# import sys
#
# lines = sys.stdin.read().splitlines() # read in the input from STDIN
# dnaWord = ""
# freq = 0
# for i in range(len(lines)):
#     if i == 0:
#         dnaWord = lines[i]
#     else:
#         freq = int(lines[i])
#
# result = freqWordProblem(dnaWord, freq)
# print(" ".join(result))
#
# # for line in inputFile:
# #     if (counter == 0):
# #         dnaText = line
# #     else:
#         query = map(lambda x : int(x), line.split())
#         patternLength = query[0]
#         windowLength = query[1]
#         freq = query[2]
# #     counter = counter + 1
#
# print(findLTPatterns(open('/Users/vijaytramakrishnan/Documents/E-coli.txt', 'r'), 9, 3, 500))
