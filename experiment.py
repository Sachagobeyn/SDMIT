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

def create_dir(res,L):
    
    for i in range(len(L)):
        if not os.path.exists(os.path.join(res,L[i])):
            os.makedirs(os.path.join(res,L[i]))     


if __name__ =="__main__":
    
    res="full_run"
    taxon = ["Baetidae","Caenidae","Chironomidae","Coenagrionidae","Corydalidae","Dytiscidae","Hydrophilidae","Hydroptilidae","Leptohyphidae","Libellulidae"]
    
    for i in taxon:
        
        taxon = "Chironomidae"
        run = os.path.join(res,i)
        create_dir("",[run])
        os.system("python SCRIPTS/script.py -t "+taxon+" -res "+run)


