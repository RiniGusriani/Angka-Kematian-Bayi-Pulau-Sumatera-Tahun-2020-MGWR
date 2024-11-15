#menjalankan library
library(readxl)
library(car)
library(lmtest)
library(spgwr)
library(fBasics)
library(AICcmodavg)
library(zoo)
library(spdep)
library(DescTools)
library(sf)
library(spatialreg)
library(raster)
library(sp) 
library(GWmodel)
library(foreign) 
library(lattice) 
library(Matrix) 
library(mvtnorm) 
library(carData)
library(MLmetrics)
library(psych)
#______________________________________________________________________________________________________________________________

#Input data
options(max.print=999999)
dataTA = read.csv("D:/Documents/TA RINI GUSRIANI/SIDANG/FIX_RINI GUSRIANI_20080108010012_DATA TA ASLI_X1-X8.csv",
                  sep = ",")
dataTA
View(dataTA)
attach(dataTA)
summary(dataTA)

#Analisis Regresi Linear Berganda
MReg<-lm(formula = Y~X1+X2+X3+X4+X5+X6+X7+X8, data = dataTA)
summary.lm(MReg)

#multikolinearitas
library(DescTools)
VIF(MReg)

#Uji Aspek Spasial
#uji heteroskedastisitas
bptest(MReg)

# Uji dependensi spasial (Moran's I)
sumatra = shapefile("D:/Documents/TA RINI GUSRIANI/SIDANG/DATA SHP/DATA TA_GQIS.shp")
X = cbind(sumatra$BBLR,sumatra$ASI,sumatra$IMD, sumatra$FASYANKES,
          sumatra$SANITASI,sumatra$`AIR MINUM`,sumatra$P.MISKIN,sumatra$RLS)
Y = sumatra$AKB
coords = coordinates(sumatra)
IDs = row.names(as(sumatra, "data.frame"))
IDs

sumatra_kn = knn2nb(knearneigh(coords, k=4), row.names = IDs)
sumatra_kn

# Row standardized weights matrix
sumatra_w = nb2listw(sumatra_kn)
listw = sumatra_w
listw

# Indeks moran
moran.test(sumatra$AKB, listw, alternative = "two.sided")

########Analisis GWR########## 
#Mencari bandwith 
#Fixed Gaussian
fixedgauss<-gwr.sel(Y~X1+X2+X3+X4+X5+X6+X7+X8,
                    coords=cbind(dataTA$u,dataTA$v),data=dataTA, adapt=FALSE,gweight=gwr.Gauss)
fixedgauss
#Estimasi Fixed Gaussian
gwr.fixedgauss<-gwr(Y~X1+X2+X3+X4+X5+X6+X7+X8,
                    coords=cbind(dataTA$u,dataTA$v),bandwidth=fixedgauss,
                    data=dataTA,hatmatrix=TRUE,gweight=gwr.Gauss)
gwr.fixedgauss

#Adaptive Gaussian
adaptivegauss=gwr.sel(Y~X1+X2+X3+X4+X5+X6+X7+X8,
                      coords=cbind(dataTA$u,dataTA$v),data=dataTA, adapt=TRUE,gweight=gwr.Gauss)
adaptivegauss
#Estimasi Adaptive Gaussian
gwr.adaptivegauss=gwr(Y~X1+X2+X3+X4+X5+X6+X7+X8,
                      coords=cbind(dataTA$u,dataTA$v),adapt=adaptivegauss,
                      data=dataTA,hatmatrix=TRUE,gweight=gwr.Gauss)
gwr.adaptivegauss
#_________________________________________________________________________________________________________________________
#Fixed Tricube
fixedtricube<-gwr.sel(Y~X1+X2+X3+X4+X5+X6+X7+X8,
                      coords=cbind(dataTA$u,dataTA$v),data=dataTA, adapt=FALSE,gweight=gwr.tricube)
fixedtricube
#Estimasi Fixed Tricube
gwr.fixedtricube<-gwr(Y~X1+X2+X3+X4+X5+X6+X7+X8,
                      coords=cbind(dataTA$u,dataTA$v),bandwidth=fixedtricube,
                      data=dataTA,hatmatrix=TRUE,gweight=gwr.tricube)
gwr.fixedtricube

#Adaptive Tricube
adapttricube<-gwr.sel(Y~X1+X2+X3+X4+X5+X6+X7+X8,
                      coords=cbind(dataTA$u,dataTA$v),data=dataTA, adapt=TRUE,gweight=gwr.tricube)
adapttricube

#Adaptive Tricube
gwr.adapttricube<-gwr(Y~X1+X2+X3+X4+X5+X6+X7+X8,
                      coords=cbind(dataTA$u,dataTA$v),adapt=adapttricube,
                      data=dataTA,hatmatrix=TRUE,gweight=gwr.tricube)
gwr.adapttricube

#____________________________________________________________________________________________________________________________
#Fixed Bisquare
fixedbisquare<-gwr.sel(Y~X1+X2+X3+X4+X5+X6+X7+X8,
                       coords=cbind(dataTA$u,dataTA$v),data=dataTA, adapt=FALSE,gweight=gwr.bisquare)
fixedbisquare

#Estimasi Fixed Bisquare
gwr.fixedbisquare<-gwr(Y~X1+X2+X3+X4+X5+X6+X7+X8,
                       coords=cbind(dataTA$u,dataTA$v),bandwidth=fixedbisquare,
                       data=dataTA,hatmatrix=TRUE,gweight=gwr.bisquare)
gwr.fixedbisquare

#Adaptive Bisquare
adaptbisquare<-gwr.sel(Y~X1+X2+X3+X4+X5+X6+X7+X8,
                       coords=cbind(dataTA$u,dataTA$v),data=dataTA, adapt=TRUE,gweight=gwr.bisquare)
adaptbisquare

#Estimasi Fixed Bisquare
gwr.adaptbisquare<-gwr(Y~X1+X2+X3+X4+X5+X6+X7+X8,
                       coords=cbind(dataTA$u,dataTA$v),adapt=adaptbisquare,
                       data=dataTA,hatmatrix=TRUE,gweight=gwr.bisquare)
gwr.adaptbisquare
#_____________________________________________________________________________________________________________________________

#Analisis GWR
#Adaptive Gaussian
adaptivegauss=gwr.sel(Y~X1+X2+X3+X4+X5+X6+X7+X8,
                      coords=cbind(dataTA$u,dataTA$v),data=dataTA, adapt=TRUE,gweight=gwr.Gauss)

#Estimasi Parameter Adaptive Gaussian
gwr1=gwr(Y~X1+X2+X3+X4+X5+X6+X7+X8,
         coords=cbind(dataTA$u,dataTA$v),adapt=adaptivegauss,
         data=dataTA,hatmatrix=TRUE,gweight=gwr.Gauss)
gwr1

#Menampilkan Nilai Bandwidth Adaptive Gaussian
gwr1$bandwidth

# Membentuk matriks euclidean
u = as.matrix(dataTA$u)
i = nrow(u)
v = as.matrix(dataTA$v)
j = nrow(v)
jarak<-matrix(0,154,154) 
for (i in 1:154){    
  for (j in 1:154) {      
    jarak[i,j]<-sqrt((u[i,]-u[j,])^2+(v[i,]-v[j,])^2)    
  }  
}  
jarak 

#Mencari nilai pembobot GWR
h<-as.matrix(gwr1$bandwidth)  
i<-nrow(h)  
W<-matrix(0,154,154)  
for (i in 1:154) {    
  for (j in 1:154) {      
    W[i,j]<-exp(-(1/2)*(jarak[i,j]/h[i,])^2)      
  }  
}  
W  

#membaca output GWR
names(gwr1)
names(gwr1$SDF)
#estimasi tiap parameter gwr masing-masing kabupaten/kota
#menampilkan nilai koefisien nilai beta
B0 = gwr1$SDF$`(Intercept)`#Beta 0
BX1 = gwr1$SDF$X1 #Beta 1
BX2 = gwr1$SDF$X2 #Beta 2
BX3 = gwr1$SDF$X3 #Beta 3
BX4 = gwr1$SDF$X4 #Beta 4
BX5 = gwr1$SDF$X5 #Beta 5
BX6 = gwr1$SDF$X6 #Beta 6
BX7 = gwr1$SDF$X7 #Beta 7
BX8 = gwr1$SDF$X8 #Beta 8
Kabupaten.Kota=dataTA$`Kabupaten/Kota` #Kabupaten/Kota
EstimasiParamerGWR = cbind(B0,BX1,BX2,BX3,BX4,BX5,BX6,BX7,BX8)
EstimasiParamerGWR

#nilai r2 lokal
R2=gwr1$SDF$localR2 
R2

#uji kesesuaian model GWR 
BFC02.gwr.test(gwr1)
ftabelgwr=qf(.95, df1=145.00, df2=73.162) 
ftabelgwr

#uji pengaruh pengaruh geografis setiap prediktor (UJI VARIABILITAS SPASIAL)
LMZ.F3GWR.test(gwr1)
#F-tabel for X1-X8
ftabelgwr.X1=qf(.95, df1=35.97623, df2=99.668) 
ftabelgwr.X2=qf(.95, df1=33.27874, df2=99.668) 
ftabelgwr.X3=qf(.95, df1=44.00736, df2=99.668) 
ftabelgwr.X4=qf(.95, df1=46.83307, df2=99.668) 
ftabelgwr.X5=qf(.95, df1=46.08666, df2=99.668) 
ftabelgwr.X6=qf(.95, df1=24.07542, df2=99.668) 
ftabelgwr.X7=qf(.95, df1=10.63367, df2=99.668) 
ftabelgwr.X8=qf(.95, df1=37.90651, df2=99.668) 
F.hitung_U.Variabilitas=cbind(ftabelgwr.X1,ftabelgwr.X2,ftabelgwr.X3,ftabelgwr.X4,
                              ftabelgwr.X5,ftabelgwr.X6,ftabelgwr.X7,ftabelgwr.X8)
F.hitung_U.Variabilitas

# Uji signifikansi GWR 
# Menampilkan t-hitung
t_X1=gwr1$SDF$X1/gwr1$SDF$X1_se
t_X2=gwr1$SDF$X2/gwr1$SDF$X2_se
t_X3=gwr1$SDF$X3/gwr1$SDF$X3_se
t_X4=gwr1$SDF$X4/gwr1$SDF$X4_se
t_X5=gwr1$SDF$X5/gwr1$SDF$X5_se
t_X6=gwr1$SDF$X6/gwr1$SDF$X6_se
t_X7=gwr1$SDF$X7/gwr1$SDF$X7_se
t_X8=gwr1$SDF$X8/gwr1$SDF$X8_se
Nilai.thitung.GWR = cbind(t_X1,t_X2,t_X3,t_X4,t_X5,t_X6,t_X7,t_X8)
Nilai.thitung.GWR

########Analisis MGWR########## 
#Mencari bandwidth
gwr.adaptivegauss=gwr(Y~X1+X2+X3+X4+X5+X6+X7+X8,
                coords=cbind(dataTA$u,dataTA$v),adapt=adaptivegauss,
                data=dataTA,hatmatrix=TRUE,gweight=gwr.Gauss)
bw=gwr.adaptivegauss$bandwidth

#Variabel pada data 
y=as.matrix(dataTA$Y)  
lat=as.matrix(dataTA$u) 
lon=as.matrix(dataTA$v) 
xg=as.matrix(cbind(dataTA$X1,dataTA$X2,dataTA$X7,dataTA$X8)) 
xl=as.matrix(cbind(dataTA$X3,dataTA$X4,dataTA$X5,dataTA$X6)) 
x=as.matrix(cbind(xl,xg))  
ng=ncol(xg) 
nl=ncol(xl) 
n=length(y) 
I=diag(1,n,n) 
w=matrix(0,n,n) 
d=matrix(0,n,n) 

#Mencari jarak euclidian dan pembobot MGWR 
for (i in 1:n){    
  for (j in 1:n)    
  {      
    d[i,j]=sqrt((lat[i,1]-lat[j,1])^2+(lon[i,1]-lon[j,1])^2)      
    if (d[i,j]>bw[i]){        
      w[i,j]=0      
    }else        
      W[i,j]=1-((d[i,j]/bw[i])^2)^2    
  }  
}
W

#Mencari estimasi parameter global MGWR
beta.l=matrix(0,nl,n)  
Sl=matrix(0,n,n)
for (i in 1:n){     
  Sl[i,]=((xl[i,]%*%(solve((t(xl)%*%diag(W[,i]))%*%xl)))%*%t(xl))%*%diag(W[,i])  
} 
beta.g=((((solve(((t(xg)%*%t(I-Sl))%*%(I-Sl))%*%xg))%*%t(xg))%*%t(I-Sl))%*%(I-Sl))%*%y 
beta.g

#Mencari estimasi parameter lokal MGWR  
for (i in 1:n){      
  beta.l[,i]=((solve((t(xl)%*%diag(W[,i]))%*%xl)%*%t(xl)%*%diag(W[,i]))%*%(y-(xg%*%beta.g)))  
} 

Sg=(xg%*%solve(t(xg)%*%xg))%*%t(xg)  
S=Sl+((((((I-Sl)%*%xg)%*%solve(((t(xg)%*%t(I-Sl))%*%(I-Sl))%*%xg))%*%t(xg))%*%t(I-Sl))%*%(I-Sl))  
y.hat=(Sl%*%y)+(((I-Sl)%*%xg)%*%beta.g)
residual=(I-S)%*%y  
H=(x%*%solve(t(x)%*%x))%*%t(x) 
beta.l
y.hat

#Pengujian model MGWR  
RSS.H0.F1=as.vector(((t(y)%*%(I-H))%*%y))  
RSS.H0.F2=as.vector((((t(y)%*%(t(I-Sl)))%*%(I-Sl))%*%y))  
RSS.H0.F3=as.vector((((t(y)%*%(t(I-Sg)))%*%(I-Sg))%*%y))  
RSS.H0=cbind(RSS.H0.F1,RSS.H0.F2,RSS.H0.F3)  
RSS.H1=as.vector((((t(y)%*%(t(I-S)))%*%(I-S))%*%y))  
DSS1=RSS.H0.F1-RSS.H1  
DSS2=RSS.H0.F2-RSS.H1  
DSS3=RSS.H0.F3-RSS.H1  
DSS=cbind(DSS1,DSS2,DSS3)  
v=c(0,0)  
u=c(0,0)  
r=c(0,0)  
t=c(0,0)  
for (i in 1:2)  
{    
  v[i]=tr(((I-H)-(t(I-S)%*%(I-S)))^i)    
  u[i]=tr((t(I-S)%*%(I-S))^i)    
  r[i]=tr((t(I-Sl)%*%(I-Sl)-t(I-S)%*%(I-S))^i)    
  t[i]=tr((t(I-Sg)%*%(I-Sg)-t(I-S)%*%(I-S))^i)  
}  

#Uji kesesuaian model MGWR 
F1=as.vector((((t(y)%*%((I-H)-(t(I-S)%*%(I-S))))%*%y)/v[1])/((((t(y)%*%t(I-S))%*%(I-S))%*%y)/u[1]))  
df1.1=(v[1]^2/v[2])  
df2=(u[1]^2)/u[2] 

#Uji serentak parameter global MGWR
F2=as.vector((((t(y)%*%(((t(I-Sl)%*%(I-Sl))-(t(I-S)%*%(I-S)))))%*%y)/r[1])/((((t(y)%*%t(I-S))%*%(I-S))%*%y)/u[1]))  
df1.2=(r[1]^2/r[2]) 

#Uji serentak parameter lokal MGWR 
F3=as.vector((((t(y)%*%(((t(I-Sg)%*%(I-Sg))-(t(I-S)%*%(I-S)))))%*%y)/t[1])/((((t(y)%*%t(I-S))%*%(I-S))%*%y)/u[1]))  
df1.3=(t[1]^2/t[2])

F=as.vector(rbind(F1,F2,F3))  
df1=c(df1.1,df1.2,df1.3)  
p.value=as.vector(matrix(0,3,1))  
for (i in 1:3)  
{    
  p.value[i]=1-(pf(F[i], df1=df1[i], df2=df2)) 
}  
Uji.Serentak=cbind(F,df1,df2,p.value)
Uji.Serentak

#F-tabel
Uji.Serentak = cbind(F,df1, df2, p.value)
Ftabel.F1 = df(0.95, df1=121.94264, df2=147.976)
Ftabel.F2 = df(0.95, df1=64.56458, df2=147.976)
Ftabel.F3 = df(0.95, df1=122.05946, df2=147.976)
F.tabel = cbind(Ftabel.F1,Ftabel.F2,Ftabel.F3)
F.tabel

#Pengujian parameter secara parsial 
#Pengujian parameter global secara parsial 
G=((solve(((t(xg)%*%t(I-Sl))%*%(I-Sl))%*%xg)%*%t(xg))%*%t(I-Sl))%*%(I-Sl)  
gkk=diag(G%*%t(G))  
t.g=as.vector(matrix(0,ng,1))  
p.val=as.vector(matrix(0,ng,1))  
sigma=as.vector(sqrt(((((t(y)%*%t(I-S))%*%(I-S))%*%y)/n)))  
for (i in 1:ng) 
{    
  t.g[i]=beta.g[i]/(sigma*sqrt(gkk[i]))    
  df=df2    
  p.val[i]=1-pt(abs(t.g[i]), df=df)  
}  
Uji.Parsial.Global=cbind(t.g,df,p.val)  
sigma=as.vector(sqrt(((((t(y)%*%t(I-S))%*%(I-S))%*%y)/n)))  
t.hit.l=matrix(0,nl,n)   
pvalue=matrix(0,nl,n)   
Uji.Parsial.Global
Uji.Parsial.Global=cbind(t.g,df,p.val)
ttabel=qt(.025, 147.976)
ttabel

#Pengujian parameter lokal secara parsial  
ringkasan=matrix(0,n,2*nl)  
for (i in 1:n) 
  
{        
  
  M=((((solve(((t(xl)%*%diag(W[,i]))%*%xl)))%*%t(xl))%*%diag(W[,i]))%*%(I-(xg%*%G)))    
  m=diag(M%*%t(M))    
  m=as.matrix(m)    
  for (j in 1:nl)    
  {      
    t.hit.l[j,i]=beta.l[j,i]/(sigma*(sqrt(m[j,])))      
    pvalue[j,i]=1-pt(t.hit.l[j,i],df=df2,lower.tail=TRUE)    
  }    
  ringkasan[i,]=t(cbind(t.hit.l[,i],pvalue[,i]))  
}  
ringkasan #nilai t, nilai pval, dst#  

print(beta.g)  
print(beta.l)  
Uji.Serentak  
RSS.H0  
RSS.H1  
DSS  
df1.1  
df2  
df1.2  
df1.3
print(Uji.Parsial.Global) 

