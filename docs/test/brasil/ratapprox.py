import sys
import numpy as np
import baryrat

from baryrat import brasil



alpha = sys.argv[1]

alpha = float(alpha)

print("alpha = ", alpha)

degree = sys.argv[2]

degree = int(degree)

print("degree = ", degree)

lower_limit = 10**(-(5+degree)/2)

r = brasil(lambda x: x**alpha, [lower_limit,1], degree)

poles = r.poles()
zeros = r.zeros()

np.savetxt('poles' + str(degree)+'.csv',poles,delimiter=',')
np.savetxt('zeros'+str(degree)+'.csv',zeros,delimiter=',')