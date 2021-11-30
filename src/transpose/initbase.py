import lumapi
import numpy as np
from collections import OrderedDict
import os
from bidict import bidict


def init_base(mode,size=(2,2)):
    mode.switchtolayout()
    # 基底生成
    mode.importmaterialdb("fdtd_files/material.mdf")
    mode.addrect(x=0, y=0.0, z=-1.11e-6,x_span=2e-5,y_span=1e-5,z_span=2e-6,name="box",material="SiO2 (Glass) - Palik")
    mode.addrect(x=0, y=0.0, z=0,x_span=2.6e-6,y_span=2.6e-6,z_span=0.22e-6,name="si",material="si")
    mode.addrect(x=-5.15e-6, y=0.0, z=0, x_span=7.7e-6, y_span=0.5e-6, z_span=0.22e-6, name="in",
                 material="si")
    mode.addrect(x=5.15e-6, y=0.75e-6, z=0, x_span=7.7e-6, y_span=0.5e-6, z_span=0.22e-6, name="out1",
                 material="si")
    mode.addrect(x=5.15e-6, y=-0.75e-6, z=0, x_span=7.7e-6, y_span=0.5e-6, z_span=0.22e-6, name="out2",
                 material="si")
    mode.addmesh(x=0, y=0, z=0, x_span=2.6e-6, y_span=2.6e-6, z_span=0.22e-6, name="mesh",
                 dx=0.02e-6,dy=0.02e-6,dz=0.02e-6,)
    mode.addvarfdtd(x=0, y=0, z=0,x_span=4e-6, y_span=4e-6, z_span=0.8e-6,simulation_time=1000e-15,x0=-1.48565e-6)

    source_props = OrderedDict([("name", "source"), ("injection axis", "x"),
                         ("x", -1.7e-6), ("y", 0),("y span", 2e-6),
                         ("selected mode number",2),
                         ("center wavelength",1.5e-6),("wavelength span",0.01e-6),
                         ("optimize for short pulse", 1)
                         ])
    mode.addmodesource(properties=source_props)

    mode.addpower(properties=OrderedDict([("name", "T1"), ("monitor type", "2D X-normal"),
                                        ("x", 1.45e-6), ("y", 0.75e-6),("z",0),
                                        ("y span", 0.75e-6),("z span",0.4e-6),
                                        ("override global monitor settings", True),("frequency points",1)]))
    mode.addmodeexpansion(properties=OrderedDict([
        ("name", "M1"),
        ("x", 1.65e-6),
        ("y", 0.75e-6),
        ("y span", 0.75e-6),
    ]))
    mode.setexpansion("input", "T1")

    mode.addpower(properties=OrderedDict([("name", "T2"), ("monitor type", "2D X-normal"),
                                        ("x", 1.45e-6), ("y", -0.75e-6),("z",0),
                                        ("y span", 0.75e-6),("z span",0.4e-6),
                                        ("override global monitor settings", True),("frequency points",1)]))

    mode.addmodeexpansion(properties=OrderedDict([
        ("name", "M2"),
        ("x", 1.65e-6),
        ("y", 0.75e-6),
        ("y span", 0.75e-6),
    ]))
    mode.setexpansion("input", "T2")

    mode.addpower(properties=OrderedDict([("name", "all"), ("monitor type", "2D Z-normal"),
                                          ("x", 0), ("y", 0), ("z", 0),
                                          ("x span", 2.6e-6), ("y span", 2.6e-6),
                                          ("override global monitor settings", True), ("frequency points", 1)]))
    # pattern生成
    rows, columns = size[0], size[1]
    mode.addstructuregroup(name="666")
    x_span = 2.6e-6 / columns
    y_span = 2.6e-6 / rows
    for row in range(rows):
        for column in range(columns):
            mode.addrect(x=-1.3e-6 + x_span / 2 + column * x_span, x_span=x_span,
                         y=1.3e-6 - y_span / 2 - row * y_span, y_span=y_span,
                         z=0, z_span=0.22e-6, name=str(row) + "r" + str(column) + "c",
                         material="etch")
            mode.addtogroup("666")

def change_pattern(mode,pattern_matrix):
    mode.switchtolayout()
    pattern_dic=bidict({
        "si":1,
        "etch":0
    })
    for row in range(pattern_matrix.shape[0]):
        for column in range(pattern_matrix.shape[1]):
            mode.select("666::"+str(row)+"r"+str(column)+"c")
            pix = mode.getObjectBySelection()
            material_name=pix.material
            if pattern_dic[material_name]!=pattern_matrix[row][column]:
                pix.material=pattern_dic.inverse[pattern_matrix[row][column]]

def fom(mode):
    t1 = mode.getresult("T1", "E")["E"].squeeze()
    t2 = mode.getresult("T2", "E")["E"].squeeze()
    a = np.sqrt(np.sum(np.power(np.abs(t1), 2), axis=1))
    b=np.sqrt(np.sum(np.power(np.abs(t2), 2), axis=1))
    return (a-b).max()

