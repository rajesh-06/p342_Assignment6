#assignmennt 6: Numerical Integration
import math as m
import random as r
import mymodule as mm

def q2(x):#Q2
        f = x/(1+x)
        return f
def q3(x):#Q3
        f = m.e**(-x**2)
        return f
def q4(x):#q4
        f = 4/(1+x**2)
        return f
    
#|f"(x)|_max=2 in 0 to 1 at 0
#|f""(x)|_max=12 in 0 to 1 at 0

#for Q2 Tabular format
print('Q2: Integration of x/(1+x)dx in [1,3]')
print('N','            ','Midpoint','                 ','Trapezoidal','                  ','Simpson')
print('6  ','   ',mm.midpoint(q2,1,3,6),'   ',mm.trapezoidal(q2,1,3,6),'   ',mm.simpson(q2,1,3,6))
print('12','   ',mm.midpoint(q2,1,3,12),'   ',mm.trapezoidal(q2,1,3,12),'   ',mm.simpson(q2,1,3,12))
print('24','   ',mm.midpoint(q2,1,3,24),'   ',mm.trapezoidal(q2,1,3,24),'   ',mm.simpson(q2,1,3,24))

#for Q3
print('\nQ3: Integration of e^(-x^2)dx in [0,1] with maximum error of 0.001')
print('With Midpoint Method: ',mm.midpoint(q3,0,1,mm.n_midpoint(2,0,1,0.001)))
print('With Trapezoidal Method: ',mm.trapezoidal(q3,0,1,mm.n_trapezoidal(2,0,1,0.001)))
print('With Simpson Method: ',mm.simpson(q3,0,1,mm.n_simpson(2,0,1,0.001)))

#for Q4 Monte Carlo method
print('\nQ4: Integration of 4/(1+x^2)dx in [0,1] with Monte Carlo method')
print('N','          ',m.pi,'                \sigma_f')
with open('Q4.txt', 'r+') as f:
        for i in range(1,10000):
                a,b=mm.monte_carlo(q4,0,1,10*i)
                print(10*i,'       ',a,'        ',b)  
                f.write(str(10*i))
                f.write('       ')
                f.write(str(a))
                f.write('\n')

'''Appending output
Q2: Integration of x/(1+x)dx in [1,3]
N              Midpoint                   Trapezoidal                    Simpson
6       1.3077156791250208     1.3051226551226551     1.306830206830207
12     1.3070695049216046     1.3064191671238379     1.306851337790899
24     1.3069070523279667     1.306744336022721     1.306852725655682

Q3: Integration of e^(-x^2)dx in [0,1] with maximum error of 0.001
With Midpoint Method:  0.7471308777479975
With Trapezoidal Method:  0.7464612610366896
With Simpson Method:  0.7471804289095103'''
