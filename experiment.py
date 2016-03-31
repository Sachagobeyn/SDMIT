# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 15:55:02 2016
Description:  
@author: sacha gobeyn (sacha.gobeyn@ugent.be or sachagobeyn@gmail.com)
"""

import pandas as pd
import numpy as np
import os
import sys

if __name__ =="__main__":

    taxon = "Baetidae"
    run = "0"
    os.system("python script.py -t "+taxon+" -res "+run)
