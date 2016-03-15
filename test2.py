import sys

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
freq = 0
for i in range(len(lines)):
    if i == 0:
        dnaWord = lines[i]
    else:
        query = list(map(int, lines[i].split()))
        patternLength = query[0]
        windowLength = query[1]
        freq = query[2]

result = findLTPatterns(dnaWord, patternLength, freq, windowLength)
print(result)
print(" ".join(result))