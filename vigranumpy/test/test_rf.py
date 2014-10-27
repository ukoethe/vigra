import vigra
import numpy as np


gaus1=10*np.random.randn(1000).reshape((500,2))+np.array([20,20])
gaus1=gaus1.astype(np.float32)
gaus2=10*np.random.randn(1000).reshape((500,2))+np.array([20,20])
gaus2=gaus2.astype(np.float32)
label_gaus1=np.ones(500).astype(np.uint32)
label_gaus2=np.zeros(500).astype(np.uint32)

def test_rf_learn():
    RF=vigra.learning.RandomForest(10)
    print "HHHHHHHHHH"
    fmat=np.vstack([gaus1,gaus2])
    lmat=np.vstack([label_gaus1,label_gaus2]).reshape(-1,1)
    RF.learnRF(fmat,lmat,1,100)