import numpy as np
from numpy import cos, sin, pi, flip
import matplotlib.pyplot as plt
import math
from math import sqrt
from scipy.interpolate import interp1d
from scipy.integrate import quad
from scipy.misc import derivative
from scipy.interpolate import CubicSpline
from scipy.interpolate import UnivariateSpline, InterpolatedUnivariateSpline
from math import hypot
from test import dbloss


def arclength(f, a, b, tol=1e-6):
    """Compute the arc length of function f(x) for a <= x <= b. Stop
    when two consecutive approximations are closer than the value of
    tol.
    """
    nsteps = 1  # number of steps to compute
    oldlength = 1.0e20
    length = 1.0e10
    while abs(oldlength - length) >= tol:
        nsteps *= 2
        fx1 = f(a)
        xdel = (b - a) / nsteps  # space between x-values
        oldlength = length
        length = 0
        for i in range(1, nsteps + 1):
            fx0 = fx1  # previous function value
            fx1 = f(a + i * (b - a) / nsteps)  # new function value
            length += hypot(xdel, fx1 - fx0)  # length of small line segment
    return length



stypex=[116,217,292,358,398,408,418,435,447,459,467,470,477]
stypey=[107,132,155,180,203,224,237,248,295,387,527,759,992]
nstypex=[61,107,254,423,453,466,475,478,478]
nstypey=[203,225,238,249,295,386,526,760,992]

stx=[]
sty=[]
nstx=[]
nsty=[]
for i in stypex:
    value=100-((i-18)*100)/460
    #print(value)
    stx.append(value)
for i in nstypex:
    value=100-((i-18)*100)/460
    #print(value)
    nstx.append(value)

for i in stypey:
    value=((i-62)*200)/932
    #print('sty',value)
    sty.append(value)
for i in nstypey:
    value=((i-62)*200)/932
    #print('stx',value)
    nsty.append(value)



fbendloss = interp1d(sty, stx)
xnewbendloss = np.linspace(min(sty), max(sty), num=141, endpoint=True)


a=71.61
b=55.34
theta=np.arange(0, pi/2, 0.1)
theta=np.flip(theta)
def deflection(a,b,theta):
    x=a*cos(theta)
    y=b*sin(theta)
    deflection=2*(sqrt(x**2+y**2))
    return deflection, x , y


B1=[]
xf1=[]
yf1=[]
defLoss1=[]
extralen1=0
for i in theta:
    d,xd,yd=deflection(a,b,i)
    B1.append(d)
    xf1.append(xd)
    yf1.append(yd)
    defLoss1.append(fbendloss(d))


#
ff1=interp1d(theta, defLoss1,kind='cubic')
xnewf1=np.linspace(min(theta), max(theta), num=141, endpoint=True)
fiber1=quad(ff1, min(theta), max(theta))

fiber1Int='Curve loss (%)= '+str('%.2f' % fiber1[0])

theta_dist=np.arange(0, pi/2, 0.25)
theta_dist=np.flip(theta_dist)
xf1_dist=[]
yf1_dist=[]
for i in theta_dist:
    d,xd,yd=deflection(a,b,i)
    #B1.append(d)
    xf1_dist.append(xd)
    yf1_dist.append(yd)
    #defLoss1.append(fbendloss(d))
#print(xf1_dist)
ff1_dist=UnivariateSpline(xf1_dist, yf1_dist)
xnewf1_dist=np.linspace(min(xf1_dist), max(xf1_dist), num=100, endpoint=True)


lenfiber1=arclength(ff1_dist, min(xf1_dist), max(xf1_dist), 1e-3)+extralen1
result1=[0,0]
totalloss1=(fiber1[0]+(result1[0])*100)
print(totalloss1,result1)
txtcurvfiber1='length Curve 1= '+str('%.2f' % (lenfiber1-extralen1))
txttotfiber1= 'extra length 1= '+str('%.2f' % extralen1)
fiber1Int='Total loss (%) = '+str('%.2f' % totalloss1)
a=102.86
b=76.82

extralen2=195.31
B2=[]
xf2=[]
yf2=[]
defLoss2=[]
for i in theta:
    d,xd,yd=deflection(a,b,i)
    B2.append(d)
    xf2.append(xd)
    yf2.append(yd)
    #print(d)
    if d>200:
        var=0
    else:
        var=fbendloss(d)
    defLoss2.append(var)
#print(xf2)
ff2=interp1d(theta, defLoss2, kind='cubic')
xnewf2=np.linspace(min(theta), max(theta), num=141, endpoint=True)
fiber2=quad(ff2, min(theta), max(theta))


xf2_dist=[]
yf2_dist=[]
for i in theta_dist:
    d,xd,yd=deflection(a,b,i)
    #B1.append(d)
    xf2_dist.append(xd)
    yf2_dist.append(yd)
    #defLoss1.append(fbendloss(d))
#print(xf1_dist)
ff2_dist=CubicSpline(xf2_dist, yf2_dist)
xnewf2_dist=np.linspace(min(xf2_dist), max(xf2_dist), num=100, endpoint=True)
#fiber2_dist=quad(ff2_dist, min(xf2_dist), max(xf2_dist))
#fiber2Int_dist='Int = '+str(fiber2_dist[0])



lenfiber2=arclength(ff2_dist, min(xf2_dist), max(xf2_dist), 1e-3)+extralen2
result2=dbloss(lenfiber2)
totalloss2=(fiber2[0]+(result2[0])*100)
print(totalloss2,result2)
txtcurvfiber2='length Curve 2= '+str('%.2f' % (lenfiber2-extralen2))
txttotfiber2= 'extra length 2= '+str('%.2f' % extralen2)
fiber2Int='Total loss (%) = '+str('%.2f' % totalloss2)
a=97.01
b=65.76


extralen3=17.57
B3=[]
xf3=[]
yf3=[]
defLoss3=[]
for i in theta:
    d,xd,yd=deflection(a,b,i)
    B3.append(d)
    xf3.append(xd)
    yf3.append(yd)
    #print(d)
    if d>200:
        var=0
    else:
        var=fbendloss(d)
    defLoss3.append(var)


ff3=interp1d(theta, defLoss3, kind='cubic')
xnewf3=np.linspace(min(theta), max(theta), num=141, endpoint=True)
fiber3=quad(ff3, min(theta), max(theta))


xf3_dist=[]
yf3_dist=[]
for i in theta_dist:
    d,xd,yd=deflection(a,b,i)
    #B1.append(d)
    xf3_dist.append(xd)
    yf3_dist.append(yd)
    #defLoss1.append(fbendloss(d))
#print(xf1_dist)
ff3_dist=CubicSpline(xf3_dist, yf3_dist)
xnewf3_dist=np.linspace(min(xf3_dist), max(xf3_dist), num=100, endpoint=True)
#fiber2_dist=quad(ff2_dist, min(xf2_dist), max(xf2_dist))
#fiber2Int_dist='Int = '+str(fiber2_dist[0])




lenfiber3=arclength(ff3_dist, min(xf3_dist), max(xf3_dist), 1e-3)+extralen3
result3=dbloss(lenfiber3)
totalloss3=(fiber3[0]+(result3[0])*100)
#print((result3[0])*100)
txtcurvfiber3='length Curve 3= '+str('%.2f' % (lenfiber3-extralen3))
txttotfiber3= 'Extra length 3= '+str('%.2f' % extralen3)
fiber3Int='Total loss (%) = '+str('%.2f' % totalloss3)


Lossx=[1,7,14]
lossy=[totalloss1,totalloss3,totalloss2]
model = np.polyfit(Lossx, lossy, 2)
predict = np.poly1d(model)
lossx_interp=[1,2,3,4,5,6,7,8,9,10,11,12,13,14]
IntTOT=quad(predict, min(lossx_interp), max(lossx_interp))

PreCalc=round(((1/(max(lossx_interp)-min(lossx_interp)))*IntTOT[0]),2)
INTTOTSTR='Average loss =' +str(PreCalc)


#print(B)
fig, axs = plt.subplots(nrows=4, ncols=3, sharex=False, sharey=False, squeeze=True, figsize=(6, 6))
axs[0,0].plot(xf1,yf1)
axs[0,0].plot(xnewf1_dist,ff1_dist(xnewf1_dist),',')
axs[0,0].text(10, 10, txtcurvfiber1, style='italic',
        bbox={'facecolor': 'red', 'alpha': 0.1, 'pad': 9})
        
axs[0,0].set_xlabel('Length (mm)')
axs[0,0].set_ylabel('Length (mm)')


axs[0,1].plot(theta,B1)
axs[0,1].set_ylabel('B. D. (mm)')
axs[0,1].set_xlabel('Angle (rad)')

axs[0,2].plot(xnewf1,ff1(xnewf1),',')
axs[0,2].plot(theta,defLoss1)
axs[0,2].text(0,2.1, fiber1Int, style='italic',
        bbox={'facecolor': 'red', 'alpha': 0.1, 'pad': 10})
#axs[0,2].text(0.1, 2.2, txttotfiber1, style='italic',
#            bbox={'facecolor': 'red', 'alpha': 0.1, 'pad': 10})
axs[0,2].set_xlabel('Loss (%)')
axs[0,2].set_ylabel('Angle (rad)')

axs[1,0].plot(xf2,yf2)
axs[1,0].plot(xnewf2_dist,ff2_dist(xnewf2_dist),',')
axs[1,0].text(10,25, txtcurvfiber2, style='italic',
        bbox={'facecolor': 'red', 'alpha': 0.1, 'pad': 10})
axs[1,0].set_xlabel('Length (mm)')
axs[1,0].set_ylabel('Length (mm)')

axs[1,1].plot(theta,B2)
axs[1,1].set_ylabel('B. D. (mm)')
axs[1,1].set_xlabel('Angle (rad)')

axs[1,2].plot(xnewf2,ff2(xnewf2),',')
axs[1,2].plot(theta,defLoss2)
axs[1,2].text(0,1.0, fiber2Int, style='italic',
        bbox={'facecolor': 'red', 'alpha': 0.1, 'pad': 10})
#axs[1,2].text(0, max(defLoss2)-0.6, txttotfiber2, style='italic',
#            bbox={'facecolor': 'red', 'alpha': 0.5, 'pad': 10})
axs[1,2].set_xlabel('Angle (rad)')
axs[1,2].set_ylabel('Loss (%)')

axs[2,0].plot(xf3,yf3)
axs[2,0].plot(xnewf3_dist,ff3_dist(xnewf3_dist),',')
axs[2,0].text(10, 30, txtcurvfiber3, style='italic',
        bbox={'facecolor': 'red', 'alpha': 0.1, 'pad': 10})
axs[2,0].set_xlabel('Length (mm)')
axs[2,0].set_ylabel('Length (mm)')

axs[2,1].plot(theta,B3)
axs[2,1].set_ylabel('B. D. (mm)')
axs[2,1].set_xlabel('Angle (rad)')

axs[2,2].plot(xnewf3,ff3(xnewf3),',')
axs[2,2].plot(theta,defLoss3)
axs[2,2].text(0,1.5, fiber3Int, style='italic',
        bbox={'facecolor': 'red', 'alpha': 0.1, 'pad': 10})
#axs[2,2].text(0, max(defLoss3)-0.6, txttotfiber3, style='italic',
#            bbox={'facecolor': 'red', 'alpha': 0.5, 'pad': 10})
axs[2,2].set_xlabel('Angle (rad)')
axs[2,2].set_ylabel('Loss (%)')

axs[3,0].plot(sty,stx,'.')
axs[3,0].plot(xnewbendloss,fbendloss(xnewbendloss),',')
axs[3,0].set_xlabel('Bending Diameter [B.D.] (mm)')
axs[3,0].set_ylabel('Loss (%)')

axs[3,1].plot(Lossx,lossy,'.')
axs[3,1].plot(Lossx,lossy,',')
axs[3,1].plot(lossx_interp,predict(lossx_interp))
axs[3,1].set_xlabel('Number of fiber')
axs[3,1].set_ylabel('Loss (%)')

axs[3,2].text(0.2, 0.5 , INTTOTSTR, style='italic',
            bbox={'facecolor': 'red', 'alpha': 0.1, 'pad': 10})
axs[3,2].set_xlabel('Not Applicable')
axs[3,2].set_ylabel('Not Applicable')
#axs[1,1].scatter(theta,B)
plt.subplot_tool()
plt.show()
