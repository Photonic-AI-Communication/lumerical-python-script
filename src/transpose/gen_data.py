from initbase import init_base,change_pattern,fom
import lumapi
import numpy as np
from collections import OrderedDict
import os
from bidict import bidict

row_num = 5
column_num = 5
mode = lumapi.MODE()
init_base(mode, size=(row_num, column_num))
mode.save("fdtd_files/base"+str(row_num)+"_"+str(column_num)+".lms")

pattern_matrix = np.random.randint(0, 2, (row_num, column_num))
change_pattern(mode,pattern_matrix)
mode.run()
best_result=fom(mode)

t1 = mode.getresult("M1","expansion for input")["T_forward"].squeeze()
print(t1)
a=input()