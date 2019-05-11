import re
import numpy as np
from WannierKit import *
import matplotlib.pyplot as plt

#wkit=WannierKit()
#wkit.readData()

#wcc=wkit.calcWannierCenter(iband_range=[0,7],startkpt=0,endkpt=7)
#print wcc

tb=TBmodel()
tb.plotWCC(algorithm='WilsonLoop')
tb.plotbands()
