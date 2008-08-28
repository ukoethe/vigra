def hook(ui, repo, **args):
    import os
    repoPath = repo.url()[5:] + '/' # cut-off 'file:' prefix
    testSuccessFile = repoPath + 'test/testSuccess'
    if not os.path.exists(testSuccessFile):
        print "File 'test/testSuccess' is missing. Run the test suite before committing."
        return True
    testTime = os.path.getmtime(testSuccessFile)
    stat = repo.status()
    files = stat[0]+stat[1]
    modified = []
    for file in files:
        if file[-4:] not in ['.cxx', '.hxx']:
            continue
        fileTime =  os.path.getmtime(repoPath+file)
        if fileTime > testTime:
            modified.append(file)
    if len(modified) > 0:
        print "Run the test suite before committing. The following files are untested:" 
        for file in modified:
            print '   ',file
        return True
    return False
    
