import sys, os

os.environ['PATH'] = os.pathsep.join([r'@VIGRAIMPEX_PATH@', os.environ['PATH']])
sys.path.insert(0, r'@VIGRANUMPYCMODULE_PATH@')
sys.path.insert(0, r'@VIGRANUMPYTEST_PATH@')
sys.path.insert(0, r'@VIGRANUMPYSCRIPTS_PATH@')
