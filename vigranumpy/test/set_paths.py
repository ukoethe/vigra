import sys, os

if len(sys.argv) == 3:
    outdir = os.sep + sys.argv[-1]
else:
    outdir = ''

os.environ['PATH'] = os.pathsep.join([r'Z:\Ilastik\vigra\src\impex%s' % outdir, os.environ['PATH']])
sys.path.insert(0, r'Z:\Ilastik\vigra\vigranumpy\src\core%s' % outdir)
sys.path.insert(0, r'Z:\Ilastik\vigra\vigranumpy\test%s' % outdir)
sys.path.insert(0, r'Z:\Ilastik\vigra\vigranumpy\test\..\lib')
