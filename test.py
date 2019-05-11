import numpy as np

x=np.array(range(1,4)).reshape((3,1))
y=np.array(range(2,5))
#print x
#print y
x1=np.array(range(10,13)).reshape((3,1))
y1=np.array(range(6,9))

M=x*y
M1=x1*y1


#eig, _ =np.linalg.eig()
V,E,W=np.linalg.svd(M)
V1,E1,W1=np.linalg.svd(M1)
#T=np.dot(V,W).transpose()
VW=np.dot(V,W)
VW1=np.dot(V1,W1)

H=np.dot(M,VW.transpose())
H1=np.dot(M1,VW1.transpose())
print H
print np.dot(np.dot(V,E*np.eye(3)),V.T)
#print H1
#print np.dot(H,H1)

eig, _=np.linalg.eig(H)
print eig
print E
#eig1,_=np.linalg.eig(H1)
#eig2,_=np.linalg.eig(np.dot(H,H1))

#print eig
#print eig1
#print eig2
#print eig*np.eye(3)
#print np.eye(3)*eig

