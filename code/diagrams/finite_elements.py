import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    x0 = 0
    x1 = 1
    x2 = 2

    xL = np.linspace(x0, x1, 100)
    xR = np.linspace(x1, x2, 100)

    h_cubic_1 = 3*((xL - x0)/(x1 - x0))**2 - 2*((xL - x0)/(x1 - x0))**3
    h_cubic_2 = 3*((x2 - xR)/(x2 - x1))**2 - 2*((x2 - xR)/(x2 - x1))**3
    h_cubic_3 = (xL - x1) * ((xL - x0)/(x1 - x0))**2
    h_cubic_4 = (xR - x1) * ((xR - x2)/(x2 - x1))**2

    fig, ax = plt.subplots(1)
    ax.plot(xL, h_cubic_1, color="blue", linestyle="solid")
    ax.plot(xR, h_cubic_2, color="blue", linestyle="dashed")
    ax.plot(xL, h_cubic_3, color="green", linestyle="solid")
    ax.plot(xR, h_cubic_4, color="green", linestyle="dashed")

    plt.show()
