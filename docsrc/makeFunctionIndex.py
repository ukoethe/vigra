#!/usr/bin/env python

import re
import glob
import sys
import xml.etree.ElementTree as ET
from itertools import chain

if len(sys.argv) != 2:
    print 'usage: python makeFunctionIndex.py directory'
    sys.exit(1)

path = str(sys.argv[1])


def getClassList():
    tree = ET.parse("tagfile.tag")
    root = tree.getroot()

    class_list = []
    for classes in root.findall("*[@kind='class']"):
        long_name = classes.find('name').text
        # only classes from vigra namespace
        if long_name[0:5] == 'vigra':
            html = classes.find('filename').text
            short_name = long_name[long_name.rfind('::')+2:]
            namespace = long_name[:long_name.rfind('::')]
            class_list.append([html, short_name, namespace])

    class_list.sort(lambda a,b: cmp(a[1], b[1]))
    return class_list



def getFunctionList():
    tree = ET.parse("tagfile.tag")
    root = tree.getroot()

    function_list = []
    for function in chain(root.findall("*[@kind='group']/*[@kind='function']"),root.findall("*[@kind='namespace']/*[@kind='function']")):
        name = function.find('name').text
        anchorfile = function.find('anchorfile').text
        anchor = function.find('anchor').text
        html = anchorfile + '#' + anchor
        function_list.append((html, name))

    # add special documentation for argument object factories
    for k in ['srcImageRange', 'srcImage', 'destImageRange', 'destImage', 'maskImage']:
        function_list.append(('group__ImageIterators.html#ImageBasedArgumentObjectFactories', k))
    for k in ['srcMultiArrayRange', 'srcMultiArray', 'destMultiArrayRange', 'destMultiArray']:
        function_list.append(('group__ImageIterators.html#MultiArrayBasedArgumentObjectFactories', k))
    for k in ['srcIterRange', 'srcIter', 'destIterRange', 'destIter', 'maskIter']:
        function_list.append(('group__ImageIterators.html#IteratorBasedArgumentObjectFactories', k))

    unique = set()
    set_add = unique.add
    function_list = [ x for x in function_list if x not in unique and not set_add(x)]    
    function_list.sort(lambda a,b: cmp(a[1], b[1]))
    function_list = disambiguateOverloadedFunctions(function_list)
    return function_list

def addHeading(index, initial):    
    index = index + '<p><a name="index_' + initial + \
    '"><table class="function_index"><tr><th> ' + initial.upper() + \
    ' </th><td align="right" width="100%">VIGRA_NAVIGATOR_PLACEHOLDER</td></tr></table><p>\n'
    return index

def disambiguateOverloadedFunctions(functionList):
    for i in xrange(len(functionList)):
        overloaded = False
        functionName = functionList[i][1]
        if i > 0:
            lastFunctionName = functionList[i-1][1]
            if functionName == lastFunctionName:
                overloaded = True
        if i < len(functionList) - 1:
            nextFunctionName = functionList[i+1][1]
            if functionName == nextFunctionName:
                overloaded = True
        if overloaded:
            # disambiguate overloaded functions by their group or namespace
            link = functionList[i][0]
            group = re.sub(r'(group__|namespacevigra_1_1)([^\.]+)\.html.*', r'\2', link)
        else:
            group = ""
        functionList[i] = functionList[i] + (group,)
    
    return functionList


def generateFunctionIndex(functionList):
    index = ""
    initials = []
    for i in range(len(functionList)):
        functionName = functionList[i][1]
        link = functionList[i][0]
        initial = functionName[0]
        if i > 0:
            lastInitial = functionList[i-1][1][0]
            if initial != lastInitial:
                initials.append(initial)
                index = addHeading(index, initial)
        else:
            initials.append(initial)
            index = addHeading(index, initial)
            
        index = index + '<a href="'+ link + '">' + functionName + '</a>()'
        overloadDisambiguation = functionList[i][2]
        if overloadDisambiguation != "":
            index = index + ' [' + overloadDisambiguation + ']'
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
    if text.find('</h1>') > -1: # up to doxygen 1.7.1
        header = text[:text.find('</h1>')+5]
    else: # for doxygen 1.7.4 to 1.7.6.1
        header = text[:re.search(r'<div class="title">[^<]*</div>\s*</div>\s*</div>(?:<!--header-->)?\n<div class="contents">',text).end()]
    footer = re.search(r'(?s)(<!-- footer.html -->.*)', text).group(1)

    text = re.sub(r'Namespace List', r'Function Index', header)
    text = text + '\n<p><hr>\n'
    text = text + index
    text = text + footer

    open(path + "/functionindex.html", 'w+').write(text)


classList = getClassList()
functionList = getFunctionList()
generateFunctionIndex(functionList)

# Export class and function list to c_api_replaces.txt for 
# crosslinking of vigranumpy documentation.
# Note that '::' are not allowed in reStructuedText link names, 
# so we have to use '.' instead.
replaces=open("../vigranumpy/docsrc/c_api_replaces.txt","w")
for i in range(len(functionList)):
    functionName = functionList[i][1]
    overloadDisambiguation = functionList[i][2]
    if i > 0 and functionName == functionList[i-1][1] and \
                   overloadDisambiguation == functionList[i-1][2]:
        continue
    if overloadDisambiguation != "":
        functionName = overloadDisambiguation +'.' + functionName
    link = functionList[i][0]
    replaces.write(functionName+":"+link+"\n")
for i in range(len(classList)):
    className = classList[i][1]
    namespace = classList[i][2]
    if (i > 0 and className == classList[i-1][1]) or \
       (i < len(classList)-1 and className == classList[i+1][1]):
        namespace = namespace.replace('::', '.')
        className = namespace +'.' + className
    link = classList[i][0]
    replaces.write(className+":"+link+"\n")
replaces.close()
