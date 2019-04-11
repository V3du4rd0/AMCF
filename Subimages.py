import PIL
#Esta funcion es para determinar ancho y alto de la imagen que se va a extraer
#Los parametros estan dados en pixeles
def get_wi_he(phi1_x,phi1_y,phi2_x,phi2_y):
	wi = math.ceil(abs(phi2_x-phi1_x))
	he = math.ceil(abs(phi2_y-phi1_y))
	#Hay que tener cuidado para el caso que los puntos de la curva son paralelas al eje x o y
	#Que hacer en esos casos.. saltarse a la proxima secuencia?
	if wi==0:
		wi=he
	elif he==0:
		he=wi
	return wi,he


#A function that gets the image itself and the file with pixel coordinates of curve 
def get_subimages(im,pix_file_name,dir="data_pix"):
	#Here we create the list.
	x,y=pyl.loadtxt(os.path.join(dir,pix_file_name),unpack=True,delimiter = ",")

	X=list(x)
	X.append(X[0])
	Y=list(y)
	Y.append(Y[0])
	
	size = [get_wi_he(X[i],Y[i],X[i+1],Y[i+1]) for i in range(len(pix_file_name))]
	if os.path.isdir('subimages')==0:
		os.mkdir('subimages')
	for k in range(len(pix_file_name)):
		global wi,he
		wi,he=get_wi_he(X[k],Y[k],X[k+1],Y[k+1])
		img = np.zeros((wi,he))
		for i in range(wi):
			for j in range(he):
				pixv=im.getpixel((j+X[k],i+Y[k]))
				img[i][j]=pixv
		del(wi,he)
		img = img.astype(np.uint8)
		# POR QUE SALE LA IMAGEN AL REVES?
		img_new=Image.fromarray(img) 
		img_new=img_new.resize((900,900), resample=PIL.Image.BICUBIC)
		route = "subimages"
		title = "subimage"+str(k)+".bmp"
		img_new.save(os.path.join(route,title))
