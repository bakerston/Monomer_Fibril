# coding=utf-8
import math
import numpy as np
from scipy.spatial.distance import cdist
"""
读取所有CA坐标
81 chains * 37 residues * 100 snaps
CA.1.1 CA.2.1 为fibril方向矢量

CA.81.* 算平均坐标为monomer位置
"""

#dist = numpy.linalg.norm(a-b)
#convert input string to list.
def str2list(input_str):
   return input_str.split()

def t(p, q, r):
    x = p-q
    return np.dot(r-q, x)/np.dot(x, x)

#get distance from point r to line section given by p, q
def dis(p, q, r):
    return np.linalg.norm(t(p, q, r)*(p-q)+q-r)

def get_com(data_in, data_len):
    ans=[]

    x_list=[x[0] for x in data_in]
    y_list=[x[1] for x in data_in]
    z_list=[x[2] for x in data_in]

    ans.append(sum(x_list)/data_len)
    ans.append(sum(y_list)/data_len)
    ans.append(sum(z_list)/data_len)

    return ans

def mid_p(p_1,p_2):
    return list(map(lambda x,y:(x+y)/2,p_1,p_2))

#setup
f=open("long_ca.dat")
line=f.readline()
tmp=[]
res_num,ori_mark=37,1

#read fibril vector for 1st time
while ori_mark<=2970:
    cur=str2list(line)
    tmp.append(list(map(float,cur)))

    line=f.readline()
    ori_mark+=1
f.close()
print(tmp[629:666])
m_18=get_com(tmp[629:666],res_num)
m_19=get_com(tmp[666:703],res_num)
m_21=get_com(tmp[740:777],res_num)
m_22=get_com(tmp[777:814],res_num)

end_a=mid_p(m_18,m_19)
end_b=mid_p(m_21,m_22)

center=mid_p(end_a,end_b)



"""
f=open("qv_ca_snap.dat","w+")
for i in range(len(q_list)):
    f.write(str(q_list[i]))
    f.write("\n")
f.write("max:")
f.write("\t")
f.write(str(max(q_list)))

f.close()"""