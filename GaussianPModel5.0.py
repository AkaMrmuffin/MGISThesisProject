import numpy as np
import pandas as pd
import scipy.stats as stats
import random
from random import randint

## function
# Leaks Initiation Function
def Init_leak(s):
    # Default Average leaks per facilities number -> 6
    # Maximum
    # s is sample size of normal distribution
    # xdim and ydim should equal to the shape of H
    # Formulate the nomral distribution for each cell
    mu=6
    maxl=12
    minl = 0
    l =range(minl,maxl+1,1)
    sigma = np.std(l)
    Pde = stats.truncnorm((minl-mu)/sigma,(maxl-mu)/sigma,loc = mu)
    dist= Pde.rvs(s)
    leak = random.choice(dist)
    leak = int(leak)

    return leak

#  Emission Size Generation Funciton
def findsize(leakn):
    # lam -> mean value of the EmissionSize distribution
    # Leakn -> Initial leaks number 2D array
    # sized -> FAQWS table
    size = np.random.choice(sized,leakn)
    Emsize = np.sum(size)

    return Emsize


def calsigma(ws,x):
    x = np.float32(x)
    if ws<2:
        sigy = 0.22*x*(1+0.0001*x)**(-1/2)
        sigz = 0.2*x
    elif(2<= ws < 5):
        sigy = 0.16*x*(1+0.0001*x)**(-1/2)
        sigz = 0.12*x
    elif(5<= ws < 6):
        sigy = 0.11*x*(1+0.0001*x)**(-1/2)
        sigz = 0.08*x*(1+0.0002*x)**(-1/2)
    else:
        sigy = 0.08*x*(1+0.0001*x)**(-1/2)
        sigz = 0.06*x*(1+0.0015*x)**(-1/2)

    return sigy,sigz


## Draw a value based on the mean and standard deviation
def drawvalue(mu,sigma):
    sob = stats.norm(loc = mu, scale= sigma)
    sdist = sob.rvs(1000)
    G = random.choice(sdist)
    return G

# Manually Create Effective Stack Height Distribution
aa = np.zeros(10000)
a = np.ones(5000)
aa[0:5000] = a
aa[5000:7000] = np.ones(2000)*2
aa[7000:8000] = np.ones(1000)*3
aa[8000:8500] = np.ones(500)*4
aa[8500:8750] = np.ones(250)*5
aa[8750:9000] = np.ones(250)*6
aa[9000:9250] = np.ones(250)*7
aa[9250:9500] = np.ones(250)*8
aa[9500:9750] = np.ones(250)*9
aa[9750:9900] = np.ones(150)*10
aa[9900:10000] = np.ones(100)*5
ESH = list(aa)

## Import Brooks Facility data
df = pd.read_csv('FacTable.csv',sep=',')
# # Mean Wind Speed
mws = df.iloc[:,2]
dd = df.iloc[:,3]

# Load the FWAQS tabel
Qdf = pd.read_excel(r'C:\Users\mozhou\Desktop\GaussianPlumeModel\FWAQS_distribution.xls')
sized = Qdf.iloc[:,3]

#### Main Model ####
rep = 1000
r = 0
MCC = []
while r<rep:
    i = 0
    # Concentration
    C = []
    while i<n:
        # Wind Speed & Downwind distance
        x = dd[i]
        ws = mws[i]

        if x < 100:

            x = 100
            sigy,sigz = calsigma(ws,x)
            ln = Init_leak(1000)
            if ln>0:
                Q = findsize(ln)
                c = (Q*np.exp(-1))/(np.pi*ws*sigy*sigz)
                C.append(c)
            else:
                C.append(0)

        else:
            # Calculate the lateral and vertical dispersion
            sigy,sigz = calsigma(ws,x)

            # Initialize the leaks
            ln = Init_leak(1000)

            if ln>0:
                # Random Select leaks size
                Q = findsize(ln)

                # Gasussian Plume Model
                H = random.sample(ESH, 1)[0]
                Hsq = H**2
                brak = np.exp(-Hsq/(2*sigz**2))
                c = (brak*Q)/(np.pi*ws*sigy*sigz)
                C.append(c)

            else:
                C.append(0)


        i += 1

    MCC.append(C)
    r+=1
    print (r)


GPM = pd.DataFrame(MCC)

GPM.to_csv('GPM5.csv',sep=',')