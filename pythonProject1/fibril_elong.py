import Bio
import pandas
import numpy as np
import math
from biopandas.pdb import pandas_pdb
import read_pdb

#Find coordinate and radius of the rotate circle.
def get_cir(cor_list):
    x1,y1,x2,y2,x3,y3=cor_list[0],cor_list[1],cor_list[2],cor_list[3],cor_list[4],cor_list[5]
    A=x1*(y2-y3)-y1*(x2-x3)+x2*y3-x3*y2
    B=(x1*x1+y1*y1)*(y3-y2)+(x2*x2+y2*y2)*(y1-y3)+(x3*x3+y3*y3)*(y2-y1)
    C=(x1*x1+y1*y1)*(x2-x3)+(x2*x2+y2*y2)*(x3-x1)+(x3*x3+y3*y3)*(x1-x2)
    D=(x1*x1+y1*y1)*(x3*y2-x2*y3)+(x2*x2+y2*y2)*(x1*y3-x3*y1)+(x3*x3+y3*y3)*(x2*y1-x1*y2)
    x_cir=-B/(2*A)
    y_cir=-C/(2*A)

    r_cir=math.sqrt((B**2+C**2-4*A*D)/(4*A**2))
    return [x_cir,y_cir,r_cir]

#alist=[57.085,67.13,56.763,67.204,56.475,67.252]
#big_list=[[57.085,67.130,56.763,67.204,56.475,67.252],[56.763,67.204,56.475,67.252,56.125,67.294],[56.475,67.252,56.125,67.294,55.815,67.348]]

#Find degree by given 2D points (point b)
def degree(a,b,c):
    ba=a-b
    bc=c-b
    cosine_angle=np.dot(ba,bc)/(np.linalg.norm(ba)*np.linalg.norm(bc))
    angle=np.arccos(cosine_angle)
    return np.degrees(angle)

'''
Read base file
Find center axis
Calculate rotate angel
Extend fibril
output to new pdb
'''


