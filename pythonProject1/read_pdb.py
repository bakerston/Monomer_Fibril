import pandas as pd
from pandas import DataFrame, Series
from biopandas.pdb import PandasPdb

import numpy as np

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

def degree(a,b,c):
    ba=a-b
    bc=c-b
    cosine_angle=np.dot(ba,bc)/(np.linalg.norm(ba)*np.linalg.norm(bc))
    angle=np.arccos(cosine_angle)
    return np.degrees(angle)



base=PandasPdb()
base.read_pdb('base.pdb')

starter=pd.DataFrame(base.df['ATOM'])
coord=DataFrame(starter,columns=['atom_name','x_coord','y_coord','z_coord'])
#print(coord.head())
#print(coord.iloc[0:2])
#print(coord.iloc[1290:1292])
xy_mean=np.array(coord.mean(axis=0))
#x_mean=54.200072=54.200
#y_mean=54.131048=54.131
#z_increase=4.697

x_mean=float(54.200)
y_mean=float(54.131)
b=np.array([x_mean,y_mean])


#z_increase=(coord.iloc[1290:]['z_coord'].sum(axis=0)-coord[:1290]['z_coord'].sum(axis=0))/1290

z_increase=4.702
print("z_increase is", z_increase)


layer1=coord.iloc[:1290]
layer2=coord.iloc[1290:]


arr_1=DataFrame(layer1,columns=['x_coord','y_coord']).to_numpy()
arr_2=DataFrame(layer2,columns=['x_coord','y_coord']).to_numpy()
result=[degree(x,b,y) for x,y in zip(arr_1,arr_2)]
angel=(sum(result)/len(result))
#angel=1.41

'''
add new layer
'''
layer3=layer2.copy()
layer3['z_coord']=layer3['z_coord']+z_increase
print(layer1)
a=layer1.iloc[[0],[1,2]])





#result=[degree(np.array(x),b,np.array(y)) for x,y in zip(layer1[['x_coord','y_coord']],layer2[['x_coord','y_coord']])]


#for i in range(3):
#    a=np.array([layer1[i]['x_coord'],layer1[i]['y_coord']])
#    b=np.array([54.2,54.131])
#    c=np.array([layer2[i]['x_coord'],layer2[i]['y_coord']])
#    sum+=fibril_elong.degree(a,b,c)






