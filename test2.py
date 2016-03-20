import sys
import math

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

def comprehensiveFreq(text, k):
    mapping = {}
    for i in range(0, len(text) - k + 1):
        pattern = text[i:i+k]
        mapping[pattern] = patternCount(text, pattern)
    return mapping

def ComputingFrequencies(Text, k):
    frequencyArray = []
    for i in range(0, int(math.pow(4, k))):
        frequencyArray.append(0)

    for i in range(0, len(Text) - k + 1):
        textVal = Text[i:i+k]
        frequencyArray[PatternToNumber(textVal)] = frequencyArray[PatternToNumber(textVal)] + 1

    return frequencyArray
def patternCount(dnaText, pattern):
    count = 0
    for i in range(0, len(dnaText) - len(pattern) + 1):
        word = dnaText[i:i+len(pattern)]
        if (word == pattern):
            count = count + 1
    return count

def comprehensiveFreq(text, k):
    mapping = {}
    for i in range(0, len(text) - k + 1):
        pattern = text[i:i+k]
        mapping[pattern] = patternCount(text, pattern)
    return mapping

def findLTPatterns(dnaText, patternLength, freq, windowLength):
    resultSet = set()
    mappings = comprehensiveFreq(dnaText[0:windowLength], patternLength)
    for clump in mappings.keys():
        sampleFreq = mappings[clump]
        if (sampleFreq >= freq):
            resultSet.add(clump)

    for i in range(1, len(dnaText) - windowLength + 1):
        word = dnaText[i-1:i-1+windowLength]
        wordPrime = dnaText[i:i+windowLength]
        firstWord = word[0:patternLength]
        lastWord = wordPrime[windowLength-patternLength: windowLength]
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

lines = sys.stdin.read().splitlines() # read in the input from STDIN
dnaWord = ""
kmer = 0
for i in range(len(lines)):
    if i == 0:
        dnaWord = lines[i]

result = PatternToNumber(dnaWord)
print(result)
