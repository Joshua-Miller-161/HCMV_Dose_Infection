import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm, skewnorm

x = np.linspace(100, 1000, 10000)
y = skewnorm.pdf(x, 3.8, 173, 200)

plt.plot(x, y)
# plt.scatter(150, 0.001, marker='x', color='black')
# plt.scatter(300, 0.001, marker='x', color='black')
#plt.plot(x, np.ones_like(x)*0.001)
plt.plot(np.ones(50)*230, np.linspace(0, max(y), 50))
plt.plot(np.ones(50)*(1.5*230), np.linspace(0, max(y), 50))
plt.show()