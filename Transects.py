import numpy as np
from collections import namedtuple
import matplotlib.pyplot as plt
import bisect

Transect = list()

def create_T(Lar_hi, Alt_hi_L, Alt_hi_R, Lar_low, Alt_low, Alt_bed, Abs_bed, Dist): 
    T = namedtuple('Transect', ['Lar_hi', 'Alt_hi_L', 'Alt_hi_R', 'Lar_low', 'Alt_low', 'Alt_bed', 'Abs_bed', 'Dist'])
    T.Lar_hi = Lar_hi
    T.Alt_hi_L = Alt_hi_L
    T.Alt_hi_R = Alt_hi_R
    T.Lar_low = Lar_low
    T.Alt_low = Alt_low ## Contains the first altitude point of reference 
    T.Alt_bed = Alt_bed ## Should be a list of size N 

    ratio = 5/2.3 ##meters in real life vs cm in sheet
    if Abs_bed == None :
        T.Abs_bed = Abs_bed
    else : 
        T.Abs_bed = [0]
        for x in Abs_bed : T.Abs_bed.append(x)
        T.Abs_bed.append(T.Lar_low)
        T.Abs_bed = ratio * np.diff(Abs_bed) ## Cumulative differences as to get distance between each altitude measuremed points ##N+1
    T.Dist = Dist
    return T

Transect.append(create_T(12.90, 33.72, 33.92, 10.50, 32.52, [31.80, 32.05, 31.70, 32.06, 31.96, 31.84, 31.93, 31.93, 31.48], 
                                                            [ 0.39,  0.40,  0.59,  1.19,   2.9,     3,  3.17,  3.61,  4.03],
                                                            81))
Transect.append(create_T(16.80, 33.62, 33.92, 13.40, 32.52, [31.89, 31.66, 32.00, 31.77, 31.36, 31.80, 31.40], 
                                                            [ 0.49,  0.99,  1.71,  2.81,  3.98,  4.08, 5.45],
                                                            206  ))
Transect.append(create_T(14.40, 34.81, 33.81, 11.60, 32.51, [31.94, 32.02, 31.82, 31.69, 31.81, 31.84, 31.79],
                                                            [  0.6,  0.77,   1.6,   3.2,   3.4,  4.63, 5],
                                                            290 ))
Transect.append(create_T(13.80, 33.71, 33.81, 11.70, 32.51, [31.73, 31.56, 31.70, 31.91, 31.96], 
                                                            [ 0.15,  1.63,  2.45,   4.4,  5],
                                                            483 ))
Transect.append(create_T(14.90, 34.81, 33.71, 11.90, 32.51, [31.73, 31.56, 31.70, 31.91, 31.96], 
                                                            [  0.3,   0.7,   1.1,   4.5, 5.3],
                                                            584 ))
Transect.append(create_T(18.30, 33.71, 33.71, 15.30, 32.51, [31.88, 32.05, 31.67, 31.56, 31.98, 31.79], 
                                                            [  0.9,   1.7,  2.55,  3.5 ,  6.  ,  6.1 ],
                                                            757 ))
Transect.append(create_T(16.80, 33.91, 33.71, 14.60, 32.51, [31.87, 31.81, 31.58, 31.59, 31.66], 
                                                            [ 0.3,   1.2,   1.85,  3.95,  5.9],
                                                            857 ))
Transect.append(create_T(16.70, 33.60, 38.50, 11.70, 32.50, [31.32, 31.31, 31.05, 31.26, 32.02, 31.87], 
                                                            [ 0.31,   1.7,   2.9,  3.35,   4.8, 5.1],
                                                            1010 ))
Transect.append(create_T(19.60, 33.50, 33.60, 15.10, 32.50, [31.25, 31.71, 31.83, 31.22, 31.52, 31.22, 31.65, 31.95], 
                                                            [ 0.1 ,  1.25,  2.15,  2.65,  4.05,   4.3,   4.6, 6.35],
                                                            1109 ))
Transect.append(create_T(13.60, 33.60, 33.50, 11.50, 32.50, [31.93, 31.85, 31.14, 31.54, 31.10, 31.04, 31.05, 31.09], 
                                                            [ 0.55,  1.25,   1.4,   1.65,  2.3,  3.25,   4.25,  4.9],
                                                            1196 ))
Transect.append(create_T(14.90, 32.46, 33.06, 9.10, 31.46,  [31.31, 31.21, 31.26, 31.12, 31.34], 
                                                            [ 0.08,  0.55,     1,   3.3,  4.15],
                                                            1306))


def rectangle_surface(x, y, true_y):
    '''
    x is the width of the rectangle
    y the bed water aka max(abs_bed[i], abs_bed[i+1])
    true_y is the altitude of the wated hence Alt_water
    ''' 
    S = (true_y - y) * x
    return max(S, 0)
def triangle_surface(x, y1, y2, true_y):
    new_x = x * (true_y-min(y1, y2))/max(0.01,np.linalg.norm(y2-y1))
    return 0.5 * new_x * (true_y-min(y1, y2)) 
def hypotenus_length(x, y1, y2, true_y):
    new_x = x * (true_y-min(y1, y2))/max(0.01,np.linalg.norm(y2-y1))
    #print(np.sqrt(new_x**2 + (true_y-min(y1, y2))**2))
    return np.sqrt(new_x**2 + (true_y-min(y1, y2))**2)

def volume_T(T, Alt_Water = None):
    if Alt_Water == None : Alt_Water = T.Alt_low ## Simple case

    Surface = 0 ## we iterate over each (rectangle+triangle) that constitues our definition of the lake. 
    for i in range(len(T.Abs_bed)):
        Surface += rectangle_surface(T.Abs_bed[i], max(T.Alt_bed[i], T.Alt_bed[i+1]), Alt_Water) ## Here we add the rectangle, heitgh implied by the max()

        if Alt_Water <= T.Alt_bed[i] and Alt_Water <= T.Alt_bed[i+1]: Surface+=0 ##bah y'a pas d'eau la 
        
        elif Alt_Water <= T.Alt_bed[i] or  Alt_Water <= T.Alt_bed[i+1]:
            Surface+= 0.5 * triangle_surface(T.Abs_bed[i], T.Alt_bed[i], T.Alt_bed[i+1], Alt_Water)

    if Alt_Water > T.Alt_low :
        #print("TODO : Here the altitude is over the reference altitude hence the surface calcululs is inexact, if this ever happend call Maxime or recode it cuz i'm busy")
        if Alt_Water > T.Alt_hi_L or Alt_Water > T.Alt_hi_R : pass#print("Consider floodind happening, and calculus even more ambiguous")
    return Surface
def surface_T(T, Alt_Water = None):
    if Alt_Water == None : Alt_Water = T.Alt_low
    Exchange_Surface = 0
    for i in range(len(T.Abs_bed)):
        #Exchange_Surface_Max = np.sqrt(T.Abs_bed[i]**2 + (T.Alt_bed[i+1]-T.Alt_bed[i+1])**2)
        if Alt_Water <= T.Alt_bed[i] and Alt_Water <= T.Alt_bed[i+1]: Exchange_Surface+=0 ##bah y'a pas d'eau la 
        elif Alt_Water <= T.Alt_bed[i] or  Alt_Water <= T.Alt_bed[i+1]:
            Exchange_Surface += hypotenus_length(T.Abs_bed[i], T.Alt_bed[i], T.Alt_bed[i+1], Alt_Water)
        else : 
            Exchange_Surface += hypotenus_length(T.Abs_bed[i], T.Alt_bed[i], T.Alt_bed[i+1], max(T.Alt_bed[i], T.Alt_bed[i+1]))
    return Exchange_Surface



def Vol_Sur(Alt_Water):
    X = np.linspace(0, 1306, 1307)

    Surface = [0]*len(X)
    Volume = [0]*len(X)
    Transect_Dist = []
    Alt_Water = Alt_Water
    for T in Transect:
        Transect_Dist.append(T.Dist)
        Surface[T.Dist] = surface_T(T, Alt_Water)
        Volume[T.Dist] = volume_T(T, Alt_Water) 

    for i,x in enumerate(X) :
        idx_T = bisect.bisect(Transect_Dist, x)
        T1, T2 = [Transect[idx_T-1], Transect[min(idx_T, 10)]]
        Surface[i] = (Surface[T2.Dist] - Surface[T1.Dist])*(x-T1.Dist)/(T2.Dist - T1.Dist) + Surface[T1.Dist]
        Volume[i]  =  (Volume[T2.Dist] - Volume[T1.Dist])*(x-T1.Dist)/(T2.Dist - T1.Dist) + Volume[T1.Dist]
    return Volume, Surface

Volumes  = []
Surfaces = []

for Alt_Water in [32.94, 33.13, 32.51, 31.2, 31.53, 32.38]:
    Volume, Surface = Vol_Sur(Alt_Water)
    Volumes.append(Volume)
    Surfaces.append(Surface)
    print(np.nansum(Surface))

plt.figure(1)
t_res_Qmna5 = [x/0.017 for x in Volumes[0][:81]]
for x in Volumes[0][:81] : pass
plt.plot(t_res_Qmna5, label = f'Temps de résidence pour Qmna5 ')
plt.plot([x for x in Volumes[0][81:]],  label = f'Volume pour Qmna5 ')
#print(len(Volumes[0][81:]))
#plt.plot([x/1.08 for x in Volumes[1][81:]],  label = f'Temps de résidence pour Module ')
#plt.plot([x/11.996 for x in Volumes[2][81:]],  label = f'Temps de résidence pour Q2 ')

#Volume_Qmna5=Volumes[0][81:]
#print(np.nansum(Volume_Qmna5)/(0.017*24*60*60))
#print(np.nansum(t_res_Qmna5)/(24*60*60))

plt.legend()
plt.savefig(f'Temps de résidence en amont avec clapet')



#plt.figure(2)
#for Sur in Surfaces :   
#    plt.plot(Sur[81:])
#plt.label(f"{Alt_Water}")
#plt.savefig(f'Surface over distance')


