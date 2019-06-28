import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    x0 = -1
    x1 = 0
    x2 = 1

    xL = np.linspace(x0, x1, 100)
    xR = np.linspace(x1, x2, 100)

    h_quadratic_1 = 4*(xL - x0)*(x1 - xL) / (x1 - x0)**2
    h_quadratic_2 = np.zeros_like(h_quadratic_1)
    h_quadratic_3 = (2*xL - x1 - x0)*(xL - x0) / (x1 - x0)**2
    h_quadratic_4 = (2*xR - x2 - x1)*(xR - x2) / (x2 - x1)**2
	
    fig, ax = plt.subplots(2, 2, figsize=(12, 8))
    ax[0][0].plot(xL, h_quadratic_1, color="blue", linestyle="solid")
    ax[0][0].plot(xR, h_quadratic_2, color="blue", linestyle="dashed")
    ax[0][0].plot(xL, h_quadratic_3, color="green", linestyle="solid")
    ax[0][0].plot(xR, h_quadratic_4, color="green", linestyle="dashed")
    ax[0][0].set_title("Quadratic elements")

    h_cubic_1 = 3*((xL - x0)/(x1 - x0))**2 - 2*((xL - x0)/(x1 - x0))**3
    h_cubic_2 = 3*((x2 - xR)/(x2 - x1))**2 - 2*((x2 - xR)/(x2 - x1))**3
    h_cubic_3 = (xL - x1) * ((xL - x0)/(x1 - x0))**2
    h_cubic_4 = (xR - x1) * ((xR - x2)/(x2 - x1))**2

    ax[0][1].plot(xL, h_cubic_1, color="blue", linestyle="solid")
    ax[0][1].plot(xR, h_cubic_2, color="blue", linestyle="dashed")
    ax[0][1].plot(xL, h_cubic_3, color="green", linestyle="solid")
    ax[0][1].plot(xR, h_cubic_4, color="green", linestyle="dashed")
    ax[0][1].set_title("Cubic elements")

    dh_quadratic_1 = 4*(-2*xL + x1 + x0) / (x1 - x0)**2
    dh_quadratic_2 = np.zeros_like(dh_quadratic_1)
    dh_quadratic_3 = (4*xL - x1 - 3*x0) / (x1 - x0)**2
    dh_quadratic_4 = (4*xR - x1 - 3*x2) / (x2 - x1)**2

    ax[1][0].plot(xL, dh_quadratic_1, color="blue", linestyle="solid")
    ax[1][0].plot(xR, dh_quadratic_2, color="blue", linestyle="dashed")
    ax[1][0].plot(xL, dh_quadratic_3, color="green", linestyle="solid")
    ax[1][0].plot(xR, dh_quadratic_4, color="green", linestyle="dashed")
    ax[1][0].set_title("Derivative quadratic elements")

    dh_cubic_1 = 6*(xL - x0)/(x1 - x0)**2 - 6*(xL - x0)**2 / (x1 - x0)**3
    dh_cubic_2 = -6*(x2 - xR)/(x2 - x1)**2 + 6*(x2 - xR)**2 / (x2 - x1)**3
    dh_cubic_3 = (2*(xL - x1)*(xL - x0) + (xL - x0)**2) / (x1 - x0)**2
    dh_cubic_4 = (2*(xR - x1)*(xR - x2) + (xR - x2)) / (x2 - x1)**2

    ax[1][1].plot(xL, dh_cubic_1, color="blue", linestyle="solid")
    ax[1][1].plot(xR, dh_cubic_2, color="blue", linestyle="dashed")
    ax[1][1].plot(xL, dh_cubic_3, color="green", linestyle="solid")
    ax[1][1].plot(xR, dh_cubic_4, color="green", linestyle="dashed")
    ax[1][1].set_title("Derivative cubic elements")

    fig.tight_layout()
    
    plt.show()
