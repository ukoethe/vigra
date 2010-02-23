#!/usr/bin/env python

import re
import glob
import sys

if len(sys.argv) != 2:
    print 'usage: python makeFunctionIndex.py directory'
    sys.exit(1)

path = str(sys.argv[1])

def getNamespaceList():
    text = open(path + "/namespaces.html").read()
    return re.findall(r'<tr><td class="indexkey"><a class="el" href="([^"]+)">([^<]+)</a>', text)

def getFunctionList(namespaceList):
    functionList = []
    for namespace in namespaceList:
        text = open(path + '/' + namespace[0]).read()
        # start of function section in the namespace file
        start = re.search(r'<tr><td colspan="2">(?:<br>)?<h2>Functions</h2></td></tr>', text)
        if not start:
            continue # no functions in this namespace
        # end of function section in the namespace file
        end = re.search(r'<tr><td colspan="2">(?:<br>)?<h2>Variables</h2></td></tr>', text)
        if not end:
            end = re.search(r'<hr/?><a name="_details"></a><h2>Detailed Description</h2>', text)
        # extract the function section from the namespace file
        text = text[start.regs[0][0]:end.regs[0][0]]
        
        # split at the function signatures to get sections for each individual function
        functionPieces = re.split(r'</a> \([^)]*\)</td></tr>', text)
        for f in functionPieces:
            # the rightmost hyperlink contains the function name and link address
            f = f[f.rfind('<a class="el" href='):]
            functionList += re.findall(r'<a class="el" href="([^"]+)">([^<]+)$', f)
    return functionList

def addHeading(index, initial):    
    index = index + '<p><a name="index_' + initial + \
    '"><table class="function_index"><tr><th> ' + initial.upper() + \
    ' </th><td align="right" width="100%">VIGRA_NAVIGATOR_PLACEHOLDER</td></tr></table><p>\n'
    return index
    
def generateFunctionIndex(functionList):
    for k in ['srcImageRange', 'srcImage', 'destImageRange', 'destImage', 'maskImage']:
        functionList.append(('group__ImageIterators.html#ImageBasedArgumentObjectFactories', k))
    for k in ['srcMultiArrayRange', 'srcMultiArray', 'destMultiArrayRange', 'destMultiArray']:
        functionList.append(('group__ImageIterators.html#MultiArrayBasedArgumentObjectFactories', k))
    for k in ['srcIterRange', 'srcIter', 'destIterRange', 'destIter', 'maskIter']:
        functionList.append(('group__ImageIterators.html#IteratorBasedArgumentObjectFactories', k))

    functionList.sort(lambda a,b: cmp(a[1], b[1]))

    index = ""
    initials = []
    for i in range(len(functionList)):
        overloaded = None
        functionName = functionList[i][1]
        link = functionList[i][0]
        initial = functionName[0]
        if i > 0:
            lastInitial = functionList[i-1][1][0]
            if initial != lastInitial:
                initials.append(initial)
                index = addHeading(index, initial)

            lastFunctionName = functionList[i-1][1]
            if functionName == lastFunctionName:
                overloaded = 1
        else:
            initials.append(initial)
            index = addHeading(index, initial)
            
        if i < len(functionList) - 1:
            nextFunctionName = functionList[i+1][1]
            if functionName == nextFunctionName:
                overloaded = 1

        index = index + '<a href="'+ link + '">' + functionName + '</a>()'
        if overloaded:
            # disambiguate overloaded functions by their group or namespace
            group = re.sub(r'(group__|namespacevigra_1_1)([^\.]+)\.html.*', r'\2', link)
            index = index + ' [' + group + ']'
        index = index + '<br>\n'

    navigator = '['
    for i in range(len(initials)):
        initial = initials[i]
        if i != 0:
            navigator = navigator + '|'
        navigator = navigator + ' <a href="#index_' + initial + '">' + initial.upper() + '</a> '
    navigator = navigator + ']'
    index = re.sub('VIGRA_NAVIGATOR_PLACEHOLDER', navigator, index)

    # use file "/namespaces.html" as boiler plate for "/functionindex.html"
    text = open(path + "/namespaces.html").read()
    header = text[:text.find('</h1>')+5]
    footer = re.search(r'(?s)(<!-- footer.html -->.*)', text).group(1)

    text = re.sub(r'Namespace List', r'Function Index', header)
    text = text + '\n<p><hr>\n'
    text = text + index
    text = text + footer

    open(path + "/functionindex.html", 'w+').write(text)

namespaceList = getNamespaceList()
functionList = getFunctionList(namespaceList)
generateFunctionIndex(functionList)
#export functions list to c_api_replaces
replaces=open("../vigranumpy/docsrc/c_api_replaces","w")
for i in range(len(functionList)):
    functionName = functionList[i][1]
    link = functionList[i][0]
    replaces.write(functionName+":"+link+"\n")

