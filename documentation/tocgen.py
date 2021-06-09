def processFile(inFile, outFile):
    mdFile = open(inFile, 'r')
    toc = []
    levels = [0,0,0,0,0 ]
    newFile = open(outFile, 'w')
    tempFile = []
    tocLoc = 0
    partOfToc = False
    
    for line in mdFile:
        if partOfToc and line != '\n':
            continue
        else:
            partOfToc = False
        if 'Table of Contents' in line:
            tocLoc = len(tempFile) + 1
            partOfToc = True
            line += "\n"
        elif line[0] == '#':
            secId = buildToc(line, toc, levels)
            line = addSectionTag(cleanLine(line), secId) + '\n'
        tempFile.append(line)


    for line in toc:
        tempFile.insert(tocLoc, line)
        tocLoc += 1

    newFile.write("\n")
    for line in tempFile:
        newFile.write(line)

    mdFile.close()
    newFile.close()

def addSectionTag(line, secId):
    startIndex = line.find(' ')
    line = line[:startIndex + 1] + '<a id=\'' + secId + '\' />' + line[startIndex + 1:]
    return line

def buildToc(line, toc, levels):
    line = cleanLine(line)
    secId = 's'
    if line[:4] == '####':
        
        #raise UserWarning('Header levels greater than 3 not supported')
        levels[4] += 1
        secId += str(levels[1]) + '-' + str(levels[2]) + '-' + str(levels[3]) + '-' + str(levels[4])
        toc.append('            - [' + line[5:] + '](#' + secId + ')\n')
    elif line[:3] == '###':
        levels[3] += 1
        secId += str(levels[1]) + '-' + str(levels[2]) + '-' + str(levels[3])
        toc.append('        - [' + line[4:] + '](#' + secId + ')\n')
    elif line[:2] == '##':
        levels[2] += 1
        levels[3] = 0
        secId += str(levels[1]) + '-' + str(levels[2])
        toc.append('    - [' + line[3:] + '](#' + secId + ')\n')
    elif line[:1] == '#':
        levels[1] += 1
        levels[3] = levels[2] = 0
        secId += str(levels[1])
        toc.append('- [' + line[2:] + '](#' + secId + ')\n')
    return secId

def cleanLine(text):
    text = stripNewline(text)
    text = removeAnchors(text)
    return text

def stripNewline(text):
    return text.replace('\n', '')

def removeAnchors(text):
    while ('<' in text and '>' in text):
        leftTag = text.index('<')
        rightTag = text.index('>')
        text = text[0:leftTag] + text[rightTag + 1:]
    return text

def clean_4_toc(inFile,outFile):
    num_str = [str(i) for i in range(1,11)]
    lines = open(inFile,'r').readlines()
    for i,line in enumerate(lines):
        if "version" in line.lower():
            lines[i] = lines[i].replace("#","")

        if line.strip().startswith("====="):
            lines[i-1] = "# " + lines[i-1]
            lines[i] = "\n"
        elif line.strip().startswith("----"):
            lines[i-1] = "## " + lines[i-1]
            lines[i] = "\n"
        elif "<img" in line:
            lines[i] = line.replace("#","")
        elif line.strip().startswith("** ") and i < 200: #trying to catch the stuff at the top
            lines[i] = "# " + line.replace("*","")
        elif "**" in line:
            lines[i] = line.replace("**","")
        elif line.strip().startswith("#") and len(line.split()) > 1 and line.lower().split()[1][0] not in num_str:
            lines[i] = "**"+line.split()[1] + "**"
    with open(outFile,'w') as f:

        for line in lines:
            if "1. introduction" in line.lower():
                f.write("\n**Table of Contents**\n\n")
            f.write(line)





if __name__ == "__main__":
    import sys
    clean_4_toc("file.md","temp.md")
    processFile("temp.md","filetoc.md")





