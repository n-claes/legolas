import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
from utilities import plot_data, read_data

from argparse import ArgumentParser

def init_global_vars():
    global w_idx_list, current_var, plot_real
    w_idx_list = []
    current_var = 0
    plot_real = True


def select_next_variable():
    global current_var
    current_var = current_var + 1
    if current_var > 7:
        current_var = 0
    return eigenf_list[current_var]


def select_prev_variable():
    global current_var
    current_var = current_var - 1
    if current_var < 0:
        current_var = 7
    return eigenf_list[current_var]


def get_current_variable():
    return eigenf_list[current_var]


def switch_real():
    global plot_real
    if plot_real:
        plot_real = False
    else:
        plot_real = True


def find_spectrum_point_idx(x, y):
    w_clicked = x + y*1j
    idx = (np.abs(omegas - w_clicked)).argmin()
    return idx


def on_clicking(event):
    toolbar = fig3.canvas.manager.toolbar
    # Toolbar should be deactivated when doing interactive stuff
    if toolbar.mode == '':
        if event.button == 1:
            # Prevent error when clicking outside of axes window
            if event.xdata is None or event.ydata is None:
                return
            # Search nearest spectrum point to click, then plot
            idx = find_spectrum_point_idx(event.xdata, event.ydata)
            # Check for double plotting of points
            if idx not in w_idx_list:
                w_idx_list.append(idx)
                ax3.plot(np.real(omegas[idx]), np.imag(omegas[idx]), 'rx',
                         markersize=8, picker=10, label='w_point')
    fig3.canvas.draw()
    return

def on_picking(event):
    global w_idx_list
    # right-click on point removes the selected point. Radius where
    # clicking is detected depends on 'picker' argument when plotting the point
    if event.mouseevent.button == 3:
        if (hasattr(event.artist, 'get_label') and
            event.artist.get_label() == 'w_point'):
            event.artist.remove()
            # Remove point from index list as well
            idx = find_spectrum_point_idx(event.artist.get_xdata(),
                                          event.artist.get_ydata())
            w_idx_list.remove(idx)
    fig3.canvas.draw()
    return

def on_typing(event):
    global w_idx_list, plot_real

    # Pressing 'x' closes figures and finishes program
    if event.key == 'x':
        for i in plt.get_fignums():
            plt.close(i)
        return

    # Do nothing if nothing is selected
    if not w_idx_list:
        return

    if event.key == "i":
        switch_real()
        ax2.clear()
        var = get_current_variable()
        plot_data.plot_eigenfunctions(fig2, ax2, omegas, grid, eigenfunctions,
                                      var, w_idx_list, plot_real)


    # Pressing 'd' clears the current figure and selection
    if event.key == 'd':
        w_idx_list = []
        for artist in ax3.get_children():
            if (hasattr(artist, 'get_label') and
                artist.get_label() == 'w_point'):
                artist.remove()
        ax2.clear()
        fig2.set_visible(False)

    if event.key == 'down':
        var = select_prev_variable()

    if event.key == 'up':
        var = select_next_variable()

    # Pressing 'enter' plots eigenfunctions corresponding to selected points
    # Pressing 'left' or 'right' cycles through variable eigenfunctions
    if (event.key == 'enter' or event.key == 'up' or event.key == 'down' or
        event.key == "i"):
        fig2.set_visible(True)
        ax2.clear()
        var = get_current_variable()
        plot_data.plot_eigenfunctions(fig2, ax2, omegas, grid, eigenfunctions,
                                      var, w_idx_list, plot_real)

    fig2.canvas.draw()
    fig3.canvas.draw()
    return





if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument("-i", "--config_file", dest="config_file")
    args = parser.parse_args()

    config_file = args.config_file

    if config_file is None:
        raise TypeError('Configuration file could not be parsed.')

    namelist = read_data.read_config_file(config_file)

    # Retrieve values from namelist. Defaults should be provided INSIDE legolas so getting these
    # should never fail. If it does, something is wrong.
    gridpts = namelist['gridlist'].get('gridpts')
    matrix_gridpts = namelist['gridlist'].get('matrix_gridpts')
    eigenf_gridpts = namelist['gridlist'].get('ef_gridpts')
    current_equil  = namelist['equilibriumlist'].get('equilibrium_type')

    # setup the different filenames
    output_folder = namelist['filelist'].get('output_folder')
    file_ext = namelist['filelist'].get('file_extension')
    filename_grid = output_folder + namelist['filelist'].get('savename_efgrid') + file_ext
    filename_matA = output_folder + namelist['filelist'].get('savename_matrix') + '_A' + file_ext
    filename_matB = output_folder + namelist['filelist'].get('savename_matrix') + '_B' + file_ext
    filename_omegas = output_folder + namelist['filelist'].get('savename_eigenvalues') + file_ext

    eigenf_list = np.asarray(["rho", "v1", "v2", "v3", "T", "a1", "a2", "a3"])
    eigenfunctions = {}

    if namelist['savelist']['show_eigenfunctions']:
        print("")
        print("-"*50)
        print(">>> INTERACTIVE PLOTTING OF EIGENFUNCTIONS <<<")
        print("- Left-click  : select spectrum points (red cross)")
        print("- Right-click : deselect spectrum points")
        print("- Enter       : plot eigenfunctions corresponding to selected points")
        print("- Up-arrow    : cycle upwards through eigenfunction variables")
        print("- Down-arrow  : cycle downwards through eigenfunction variables")
        print("- d-key       : clears selection and eigenfunction figure")
        print("- x-key       : close figures and finish program ")
        print("-"*50)
        print("")


    # =============== READING DATA ================

    # Read eigenvalues
    omegas   = read_data.read_stream_data(filename_omegas, content_type=np.complex,
                                          rows=1, cols=matrix_gridpts)[0]

    # Read matrices
    if namelist['savelist']['show_matrices']:
        # Read matrix elements
        matrix_B = read_data.read_stream_data(filename_matB,
                                              content_type=np.float64,
                                              rows=matrix_gridpts, cols=matrix_gridpts)
        matrix_A = read_data.read_stream_data(filename_matA,
                                              content_type=np.complex,
                                              rows=matrix_gridpts, cols=matrix_gridpts)

    # Read eigenfunctions
    if namelist['savelist']['show_eigenfunctions']:
        # Read grid data
        grid = read_data.read_stream_data(filename_grid,
                                          content_type=np.float64,
                                          rows=1, cols=eigenf_gridpts)[0]

        # Every row has the eigenvector corresponding to the
        # eigenvalue at the same index in the 'omegas' array.
        for i in range(1, 8+1):
            varname = eigenf_list[i-1]
            ef_name = "output/eigenfunctions/{}_{}_eigenfunction.dat".format(i,
                                                                        varname)
            eigenfunctions[varname] = read_data.read_stream_data(ef_name,
                                                                 content_type=np.complex,
                                                                 rows=matrix_gridpts,
                                                                 cols=eigenf_gridpts)


    # ================= CREATE FIGURES FOR PLOTTING ============
    if namelist['savelist']['show_matrices']:
        fig1, ax1 = plt.subplots(1, 2, figsize=(12, 8))
        plot_data.plot_matrix(fig1, ax1, 0, matrix_A, title="Matrix A")
        plot_data.plot_matrix(fig1, ax1, 1, matrix_B, title="Matrix B", log=False)

    if namelist['savelist']['show_eigenfunctions']:
        fig2, ax2 = plt.subplots(1, 1, figsize=(12, 8))
        fig3, ax3 = plt.subplots(1, 1, figsize=(12, 8))
        # Initialise global variables in this case
        init_global_vars()
        # Connect interactive functions with figures
        fig2.canvas.mpl_connect('key_press_event', on_typing)
        fig3.canvas.mpl_connect('key_press_event', on_typing)
        fig3.canvas.mpl_connect('button_press_event', on_clicking)
        fig3.canvas.mpl_connect('pick_event', on_picking)
        fig2.set_visible(False)
    else:
        # If no eigenfunctions are required, no need to connect interactively
        fig3, ax3 = plt.subplots(1, 1, figsize=(12, 8))

    # Plot eigenvalues
    plot_data.plot_spectrum(fig3, ax3, omegas, title=current_equil)




    plt.show()

    print("\nProgram finished")
