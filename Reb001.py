import numpy as np
from math import sqrt,pi,cos,sin
import datetime as dat
import os
import subprocess
import sys
import itertools as iter
from PIL import Image, ImageDraw
import networkx as ntx
from copy import deepcopy
import random as rnd
from pylab import loadtxt,show,plot
import math

# ---------------------
#  var auxiliares

tau=0.0
tt=np.linspace(2*np.pi,0.,31)
tt=list(tt)
tt.pop()
tt=np.array(tt)

#
# ---------------------

class DX:
	def __init__(self, image, data, curve1,curve2):
		self.image = image # image name
		self.data = data # name of 1st approximation curve
		self.curve1 =curve1 # Set of 1st approximation points
		self.curve2=curve2 # complete final set of coordinates
tlim='54999'
Ldxs=[]

imgs_id=os.listdir("imgs_ideal")
data={}
for carpeta in os.listdir("data_ini"):
	L=os.listdir("data_ini/"+carpeta)
	L.sort()
	if (tlim in L[-1]) == False:
		#data[carpeta[4:]]=L[-1]
		Ldxs.append(DX("imgs_ideal/"+carpeta[:],"data_ini/"+carpeta+"/"+L[-1],loadtxt("data_ini/"+carpeta+"/"+L[-1]), []))
#data=[os.listdir("data_ini/"+carpeta) for carpeta in os.listdir("data_ini")]
#Lerror=[label for label in data.keys() if (label in imgs) ==0]
#if len(Lerror)>0:
#	print Lerror

for dxs in Ldxs[:]:
	# ---------------------------
	# read image
	# 
	#
	im1=Image.open(dxs.image)
	im2=np.full((1800,1800),255.0)
	for i in range(600):
		for jj in range(450):
			im2[jj+675][i+600] = im1.getpixel((i,jj))  # im2==img  !!!!
	im = Image.fromarray(im2)
	wi,he=im.size   
	img=np.zeros((wi,he))
	fx=lambda x: (x-float(wi)*0.5)/100.
	fy=lambda y: (-y+float(he)*0.5)/100.
	f_x=lambda x: x*100.+float(wi)*0.5
	f_y=lambda y: -1.*(y*100.-float(he)*0.5)
	#for i in range(wi):
	#	for j in range(he):
	#		pixv=im.getpixel((i,j))# (0,0,0) ,(255,255,255)
	#		img[i][j]=pixv
	# -----------------------------
	N=len(dxs.curve1)
	M=dxs.curve1
	Ldata=["data/"+elem for elem in os.listdir("data/")]
	coords=[]
	coords2=[]
	for j in range(N):
		c= [(M[j][0]+M[(j+1)%N][0])/2.0,(M[j][1]+M[(j+1)%N][1])/2.0]
		RR=np.linalg.norm(M[(j+1)%N]-c)
		tau = math.acos((M[(j+1)%N][0]-c[0]) *(1/RR))
		#print(tau)
		M1=np.empty((N,2))
		for i in range(N):
			condicion=tt[0]+tau
			#print(condicion)
			if tau<np.pi and M[(j+1)%N][1]>M[j][1]:
				M1[i][0]=(M[j][0]+M[(j+1)%N][0])/2.0 + (np.linalg.norm(M[j]-M[(j+1)%N])/2.0)*cos(tt[i]+tau+np.pi)
				M1[i][1]=(M[j][1]+M[(j+1)%N][1])/2.0 + (np.linalg.norm(M[j]-M[(j+1)%N])/2.0)*sin(tt[i]+tau+np.pi)
			else:
				M1[i][0]=(M[j][0]+M[(j+1)%N][0])/2.0 + (np.linalg.norm(M[j]-M[(j+1)%N])/2.0)*cos(tt[i]-tau+np.pi)
				M1[i][1]=(M[j][1]+M[(j+1)%N][1])/2.0 + (np.linalg.norm(M[j]-M[(j+1)%N])/2.0)*sin(tt[i]-tau+np.pi)
		arcN=j
		if j == arcN:
			#print ((M[j][0]+M[(j+1)%N][0])/2.0, np.linalg.norm(M[j]-M[(j+1)%N])/2.0, tau+np.pi,(M[j][1]+M[(j+1)%N][1])/2.0)
	

			X1=np.array([M1[i][0] for i in range(N)])
			Y1=np.array([M1[i][1] for i in range(N)])
			maX,miX,maY,miY=max(X1),min(X1),max(Y1),min(Y1)
		
			Axy,Bxy,Cxy,Dxy=[miX,maY],[maX,maY],[miX,miY],[maX,miY]
			Aij,Bij,Cij,Dij=np.array([f_x(miX),f_y(maY)]).astype("int32"),np.array([f_x(maX),f_y(maY)]).astype("int32"),np.array([f_x(miX),f_y(miY)]).astype("int32"),np.array([f_x(maX),f_y(miY)]).astype("int32")
		
			# bounding box
			wibb=hebb=max((Bij-Aij)[0],(Cij-Aij)[1])
			#hebb=(Cij-Aij)[1]
			subI=np.full((wibb,hebb),255).astype("uint8")
			for r in range(wibb):
				for s in range(hebb):
					subI[s][r]=im.getpixel(( Aij[0]+r ,Aij[1]+s  ))
				
			img_new=Image.fromarray(subI) 
			img_new=img_new.resize((900,900), resample=Image.BICUBIC)
			#fm=fg[:-4]+"__"+lastCurve[:-4]+'__'+str(arcN)+'__.bmp'
		
			#Ly=[sum([img_new.getpixel((r,s))  for r in range(900) ]) for s in range(900)]
			#Lx=[sum([img_new.getpixel((r,s))  for s in range(900) ]) for r in range(900)]
			#xim=f_x(Lx[Lx.index(max(Lx))])
			#yim=f_y(Ly[Ly.index(max(Ly))])
		
			#fm=fg[:-4]+"__"+lastCurve[:-4]+'__'+str(arcN)+','+str(xim)+','+str(yim)+',_.bmp'
		
		
			#img_new.save("cutter/"+fm) ## Una vez cortadas ya no hace falta
			#subimgs=["imgs/"+elem for elem in os.listdir("imgs") if (dxs.image.split("/")[-1][:-4] in elem) == True]
			#subimgs.sort()
			#scaleF=0.25*(np.linalg.norm(Bxy-Axy)+np.linalg.norm(Dxy-Cxy)+np.linalg.norm(Axy-Cxy)+np.linalg.norm(Dxy-Bxy))
			#C=0.5*np.array([fx(850)+fx(450),fy(850)+fy(450)])
			axy,bxy,cxy,dxy=np.array([fx(450),fy(450)]),np.array([fx(1350),fy(450)]),np.array([fx(450),fy(1350)]),np.array([fx(1350),fy(1350)])
			#for fig in subimgs:
			#	Lfig=["data/"+fig.split("/")[1]+"/"+f for f in os.listdir("data/"+fig.split("/")[1])]
			#	Lfig.sort()
			#f="data/+"+dxs.data.split("/")[2][:-4]+"__"+dxs.data.split("/")[1][4:-4]+"__"+str(j).zfill(3)+"_.bmp'"#Lfig[-1]
			ap=dxs.data.split("/")[1]
			f=[elem for elem in Ldata if "__"+str(j)+"," in elem]#[elem for elem in Ldata if (ap[4:-4]+"__"+str(j).zfill(3) in elem)]
			if len(f)==0:
				pass
			else:
				f=f[0]
				print f
				ff=os.listdir(f)
				ff.sort()
				if "54999" not in ff[-1]:
					Q=loadtxt(f+"/"+ff[-1])
					QQ=np.zeros(np.shape(Q))
					for i in range(len(Q)):
						QQ[i][0]=Axy[0]+(Q[i][0]-axy[0])*(Bxy[0]-Axy[0])/(bxy[0]-axy[0])
						QQ[i][1]=Cxy[1]+(Q[i][1]-cxy[1])*(Axy[1]-Cxy[1])/(axy[1]-cxy[1])
						coords2.append(QQ[i])
					NN=len(Q)
					DList={}
					for k in range(NN):
						#dvec=QQ[(k+1)%N]-QQ[k]
						#dnor=np.array([-dvec[1],dvec[0]])
						#dnor=dnor/np.linalg.norm(dnor)
						t=np.dot(QQ[k]-M[j],QQ[k]-M[j])#np.dot(dnor,QQ[k]-M[j])
						DList[k]=abs(t)
					kcandidate=rnd.choice([k for k in range(NN) if DList[k]==min([DList[kk] for kk in range(NN)])])#[k for elem in range(NN) if DList[k]==min([DList[elem] for elem in DList.keys()])])
					print kcandidate
					DList2={}
					for k in range(NN):
						#dvec=QQ[(k+1)%N]-QQ[k]
						#dnor=np.array([-dvec[1],dvec[0]])
						#dnor=dnor/np.linalg.norm(dnor)
						t=np.dot(QQ[k]-M[(j+1)%N],QQ[k]-M[(j+1)%N])#np.dot(dnor,QQ[k]-M[(j+1)%N])
						DList2[k]=abs(t)
					k1candidate=rnd.choice([k for k in range(NN) if DList2[k]==min([DList2[kk] for kk in range(NN)])])#[k for elem in range(NN) if DList[k]==min([DList[elem] for elem in DList.keys()])])#DList2.index(min(DList2))
					print k1candidate
				m1=kcandidate
				m2=k1candidate
				
				if kcandidate>k1candidate:
					m2+=N
				print m1,m2,N
				#m1,m2=min([kcandidate,k1candidate]),max([kcandidate,k1candidate])
				#k1cadidate=max(aux)
				#kcandidate=min(aux)
				
				#subL=range(kcandidate,k1candidate+1)
				subL=range(m1,m2+1)
				for elem in subL:
						coords.append(QQ[elem%N])
	
				X=[elem[0] for elem in coords]
				Y=[elem[1] for elem in coords]
				X2=[elem[0] for elem in coords2]
				Y2=[elem[1] for elem in coords2]
				#plot(X2,Y2,'-ro')
				plot(X,Y,'-bo',ms=2)
plot([M[j][0] for j in range(N)], [M[j][1] for j in range(N)],'-yo')

import matplotlib.pyplot as plt
import pylab as pyl
#%matplotlib inline
im1 = Image.open("imgs_ideal/ISIC_0024370_1.tif")#read_"+fg[:-4]+".bmp")
#im1=Image.open(IMR)
im2=np.full((1800,1800),255.0)
for i in range(600):
	for j in range(450):
		im2[j+675][i+600] = im1.getpixel((i,j))
plt.figure(1)

implot = plt.imshow(im2,extent=[-9,9,-9,9] ,cmap='gray')
plt.grid(True)

show()
"""
imgs.sort()
for im in imgs[:]:
	arg1=im#gs[0]
	IMR="./imgs/"+arg1
	#+time_000000000041067__ISIC_0024316_1__000_.bmp
	#+time_000000000041067__ISIC_0024316_1__000_.bmp
	
	curvefiles=os.listdir("./data/"+im)
	curvefiles.sort()
	coordsf=curvefiles[:-1]#"./data/time_000000000041067.txt"
	# ------------------------------
	im1=Image.open(IMR)
	im2=np.full((1800,1800),255.0)
	for i in range(900):
		for j in range(900):
			im2[j+450][i+450] = im1.getpixel((i,j))
	#im=Image.open(IMR)
	im = Image.fromarray(im2)

	wi,he=im.size   
	img=np.zeros((wi,he))

	fx=lambda x: (x-float(wi)*0.5)/100.
	fy=lambda y: (-y+float(he)*0.5)/100.

	f_x=lambda x: x*100.+float(wi)*0.5
	f_y=lambda y: -1.*(y*100.-float(he)*0.5)
	for i in range(wi):
		for j in range(he):
			pixv=im.getpixel((i,j))# (0,0,0) ,(255,255,255)
			img[i][j]=pixv
	# -----------------------------
	#L=[sum([1-img[i][j] for j in range(wi)]) for i in range(he) ]
	#LX=[L.index(elem) for elem in L if elem==max(L)]
	#qx=sum(LX)/(1.0*len(LX))

	#L=[sum([1-img[i][j] for i in range(he)]) for j in range(wi) ]
	#LY=[L.index(elem) for elem in L if elem==max(L)]
	#qy=sum(LY)/(1.0*len(LY))
	
	#Q=iter.product(range(1800),range(1800))
	##LL=[(i,j) for i,j in Q if img[i][j]==0]
	#print len(LL)
	#LL2=deepcopy(LL)
	#G=ntx.Graph()
	#while len(LL2)!=0:
	#	elem=rnd.choice(LL2)
	#	for tgt in LL:
	#		if elem!= tgt:
	#			if (elem[0] in list(tgt)) or (elem[1] in list(tgt)):
	#				G.add_edge(elem, tgt)
	#	LL2.remove(elem)		
			
	# ------------------------------

	#arg2= fx(qx)   #carga-pos
	#arg3=fy(qy)
	#print qx," , ", qy
	subprocess.call("python Reb001.py "+coordsf, shell=True)
"""