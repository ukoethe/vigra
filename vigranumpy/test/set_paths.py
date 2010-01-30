import sys, os

if len(sys.argv) == 3:
    outdir = os.sep + sys.argv[-1]
else:
    outdir = ''

os.environ['PATH'] = os.pathsep.join([r'@VIGRAIMPEX_PATH@%s' % outdir, os.environ['PATH']])
sys.path.insert(0,r'../src')
#sys.path.insert(0, r'@VIGRANUMPYCMODULE_PATH@%s' % outdir)
sys.path.insert(0, r'@VIGRANUMPYTEST_PATH@%s' % outdir)
sys.path.insert(0, r'@VIGRANUMPYSCRIPTS_PATH@')
