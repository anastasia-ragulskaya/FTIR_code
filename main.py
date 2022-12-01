#!/usr/bin/env python
# img_viewer.py

import os.path

print('loading PySimpleGui')
import PySimpleGUI as sg

print('loading matplotlib')
import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt

print('loading numpy')
import numpy as np

print('loading scipy.optimize')
from scipy.optimize import curve_fit

print('loading parameters and defining a window')

matplotlib.use('TkAgg')

sg.ChangeLookAndFeel('DarkBlue')
NAME_SIZE = 14

# Define functions which are used
def name(name):
    """Create dot layout for the Window"""
    dots = NAME_SIZE - len(name) - 2
    return sg.Text(name + ' ' + '•' * dots, size=(NAME_SIZE, 1), justification='r', pad=(0, 0), font='Courier 10')


def popup_text(filename, text):
    """ Create a popup window with the information, loaded from the chosen file"""
    layout = [
        [sg.Multiline(text, size=(80, 25)), ],
    ]
    win = sg.Window(filename, layout, modal=True, finalize=True)

    while True:
        event, values = win.read()
        if event == sg.WINDOW_CLOSED:
            break
    win.close()


def draw_figure(canvas, figure):
    """ Draw the figure in the chosen canvas"""
    figure_canvas_agg = FigureCanvasTkAgg(figure, canvas)
    figure_canvas_agg.draw()
    figure_canvas_agg.get_tk_widget().pack(side='top', fill='both', expand=1.2)
    return figure_canvas_agg


def import_data(text):
    """ Import data from the .txt file. The outcome of the FTIR measurement is the .txt file with the
     wave number and absorption, written in the first and second rows, correspondingly.
     The function returns datax (wavenumber) and datay (absorption). """
    datax = []
    datay = []
    for line in text:
        datax.append(float(line.split('\t')[0]))
        datay.append(float(line.split('\t')[1]))
    return datax, datay


def make_plot(datax, datay, mode='scatter'):
    """Create a scatter or lineplot."""
    # plt.figure(figsize=(6, 5))
    # fig = Figure(figsize=(5, 4), dpi=100)
    max_data_point = max(np.array(datay)) + 0.1

    if mode == 'scatter':
        plt.scatter(np.array(datax), np.array(datay))
    else:
        plt.plot(np.array(datax), np.array(datay))
    plt.xlabel(r'$Wavenumber (cm^{-1})$')
    plt.ylabel(r'$Absorption (a.u.)$')
    plt.ylim([0, max_data_point])


def delete_figure_agg(fig_canvas_agg):
    """Delete a figure which is presently displayed in the canvas window"""
    fig_canvas_agg.get_tk_widget().forget()
    plt.close('all')


def extract_amid1_band(datax, datay):
    """Extract Amid I band region from the spectrum. It corresponds to the 1600-1700 cm-1 wave number range."""
    datax_amid1 = []
    datay_amid1 = []
    for i in range(len(datax)):
        if (datax[i] < 1700) & (datax[i] > 1600):
            datax_amid1.append(datax[i])
            datay_amid1.append(datay[i])
    return datax_amid1, datay_amid1


def subtract_background(datax_amid1, datay_amid1):
    """Perform the linear background subtraction. The background is the line, which crosses the initial and final point
    of the data. The outcome is the difference between the initial data and the background. """
    grad = (datay_amid1[0] - datay_amid1[-1]) / (datax_amid1[0] - datax_amid1[-1])
    intercept = (datay_amid1[-1] * datax_amid1[0] - datay_amid1[0] * datax_amid1[-1]) / (
            datax_amid1[0] - datax_amid1[-1])
    datay_amid_b = [(datay_amid1[i] - grad * datax_amid1[i] - intercept) for i in range(len(datax_amid1))]
    return datay_amid_b


def fit_gauss5(x, y, bounds=((0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (1, 1, 1, 1, 1, 100, 100, 100, 100, 100))):
    """Perform fit of the data with 5 gauss, which correspond to different secondary structures of the protein.
    The band positions of the secondary structures are fixed, according to the literature (dissolved in D2O):
    alpha-helix      1652 cm-1
    beta-sheet       1630 cm-1
                     1679 cm-1
    Loop             1671 cm-1
    Disordered       1645 cm-1"""

    def Gaus5(x, A1, A2, A3, A4, A5, s1, s2, s3, s4, s5):
        y = (A1 * np.exp(-(x - 1652) ** 2 / (2 * s1 ** 2)) +
             A2 * np.exp(-(x - 1630) ** 2 / (2 * s2 ** 2)) +
             A3 * np.exp(-(x - 1679) ** 2 / (2 * s3 ** 2)) +
             A4 * np.exp(-(x - 1671) ** 2 / (2 * s4 ** 2)) +
             A5 * np.exp(-(x - 1645) ** 2 / (2 * s5 ** 2))
             )
        return y

    parameters, covariance = curve_fit(Gaus5, x, y, bounds=bounds)

    error_bars = np.sqrt(np.diag(covariance))

    A1, A2, A3, A4, A5, s1, s2, s3, s4, s5 = parameters
    y_fit = (A1 * np.exp(-(x - 1652) ** 2 / (2 * s1 ** 2)) + A2 * np.exp(-(x - 1630) ** 2 / (2 * s2 ** 2)) +
             A3 * np.exp(-(x - 1679) ** 2 / (2 * s3 ** 2)) + A4 * np.exp(-(x - 1671) ** 2 / (2 * s4 ** 2)) +
             A5 * np.exp(-(x - 1645) ** 2 / (2 * s5 ** 2)))
    return y_fit, parameters, error_bars


def display_fit_results(parameters, error_bars):
    A = []
    s = []
    for i in range(len(parameters)):
        if i < 5:
            A.append('{:.4f} +/- {:.4f}'.format(parameters[i], error_bars[i]))
        else:
            s.append('{:.4f} +/- {:.4f}'.format(parameters[i], error_bars[i]))
    fit_results_table = [['\u03B1 - Helix', '1652', A[0], s[0]], ['\u03B2 - Sheet', '1630', A[1], s[1]],
                         ['\u03B2 - Sheet', '1679', A[2], s[2]], ['Loop', '1671', A[3], s[3]],
                         ['Disordered', '1645', A[4], s[4]]]
    return fit_results_table

def save_to_file_as_txt(filename, datax, datay):
    n_names = ["{}\t{}\n".format(i, j) for i,j in zip(datax, datay)]
    with open(filename, 'w') as file:
        file.writelines(n_names)

def create_file_to_save(xdata, ydata):
    data = list(np.vstack((np.array(xdata), np.array(ydata))).T)
    return data






# ----------Window Layout-----------------

file_list_column = [
    [
        sg.Text("Data Folder"),
        sg.In(size=(25, 1), enable_events=True, key="-FOLDER-"),
        sg.FolderBrowse(),
    ],
    [
        sg.Listbox(
            values=[], enable_events=True, size=(40, 4), key="-FILE LIST-"
        )
    ],
    [
        sg.Button('I. Extract Amide I band'),
        sg.Button('II. Background subtraction'),
        sg.Button('III. Fit')
    ],
    [
        name(''),
        sg.Text('Fit function of Amide I modes:', font='Courier 12')

    ],
    [
        sg.Text("Absorption (x) = A0exp(-(x-1652)\u00b2/(2s0\u00b2))+"
                "A1exp(-(x-1630)\u00b2/(2s1\u00b2))+A2exp(-(x-1679)\u00b2/(2s2\u00b2))+"
                "A3exp(-(x-1671)\u00b2/(2s3\u00b2))+A4exp(-(x-1645)\u00b2/(2s4\u00b2))", size=(50, 5),
                font='Courier 12')

    ],
    [
        name(''),
        sg.Text('Boundaries for fit parameters:', font='Courier 12')

    ],
    [
        sg.Text('                                   '),
        sg.Text('Min.', size=(10, 1)),
        sg.Text('Max.', size=(10, 1)),
        sg.Text('                        '),
        sg.Text('Min.', size=(10, 1)),
        sg.Text('Max.', size=(10, 1))

    ],
    [
        name('A0'),
        sg.Slider((0, 4), default_value=0, resolution=0.01, orientation='h', s=(10, 15), key='-A0MIN-'),
        sg.Slider((0, 4), default_value=1, resolution=0.01, orientation='h', s=(10, 15), key='-A0MAX-'),
        name('A1'),
        sg.Slider((0, 4), default_value=0, resolution=0.01, orientation='h', s=(10, 15), key='-A1MIN-'),
        sg.Slider((0, 4), default_value=1, resolution=0.01, orientation='h', s=(10, 15), key='-A1MAX-')

    ],
    [
        name('A2'),
        sg.Slider((0, 4), default_value=0, resolution=0.01, orientation='h', s=(10, 15), key='-A2MIN-'),
        sg.Slider((0, 4), default_value=1, resolution=0.01, orientation='h', s=(10, 15), key='-A2MAX-'),
        name('A3'),
        sg.Slider((0, 4), default_value=0, resolution=0.01, orientation='h', s=(10, 15), key='-A3MIN-'),
        sg.Slider((0, 4), default_value=1, resolution=0.01, orientation='h', s=(10, 15), key='-A3MAX-')

    ],
    [
        name('A4'),
        sg.Slider((0, 4), default_value=0, resolution=0.01, orientation='h', s=(10, 15), key='-A4MIN-'),
        sg.Slider((0, 4), default_value=1, resolution=0.01, orientation='h', s=(10, 15), key='-A4MAX-'),
        name('s0'),
        sg.Slider((0, 100), default_value=0, resolution=1, orientation='h', s=(10, 15), key='-S0MIN-'),
        sg.Slider((0, 100), default_value=100, resolution=1, orientation='h', s=(10, 15), key='-S0MAX-')

    ],
    [
        name('s1'),
        sg.Slider((0, 100), default_value=0, resolution=1, orientation='h', s=(10, 15), key='-S1MIN-'),
        sg.Slider((0, 100), default_value=100, resolution=1, orientation='h', s=(10, 15), key='-S1MAX-'),
        name('s2'),
        sg.Slider((0, 100), default_value=0, resolution=1, orientation='h', s=(10, 15), key='-S2MIN-'),
        sg.Slider((0, 100), default_value=100, resolution=1, orientation='h', s=(10, 15), key='-S2MAX-')

    ],
    [
        name('s3'),
        sg.Slider((0, 100), default_value=0, resolution=1, orientation='h', s=(10, 15), key='-S3MIN-'),
        sg.Slider((0, 100), default_value=100, resolution=1, orientation='h', s=(10, 15), key='-S3MAX-'),
        name('s4'),
        sg.Slider((0, 100), default_value=0, resolution=1, orientation='h', s=(10, 15), key='-S4MIN-'),
        sg.Slider((0, 100), default_value=100, resolution=1, orientation='h', s=(10, 15), key='-S4MAX-')

    ],
    [
        name(''),
        #sg.FileSaveAs(enable_events=True, key='-SAVE-')
        sg.Button('Save data after baseline correction', key='-SAVE-')

    ],

]

# For now will only show the name of the file that was chosen
image_viewer_column = [
    [sg.Text("Choose a file from list on left:")],
    [sg.Text(size=(60, 1), key="-TOUT-")],
    [sg.Canvas(key="-DATA-")],
    [
        sg.Text('Fit results:', font='Courier 12')

    ],
    [
        sg.Table([['\u03B1 - Helix', '1652', 'Nan', 'Nan'], ['\u03B2 - Sheet', '1630', 'Nan', 'Nan'],
                  ['\u03B2 - Sheet', '1679', 'Nan', 'Nan'], ['Loop', '1671', 'Nan', 'Nan'],
                  ['Disordered', '1645', 'Nan', 'Nan']],
                 headings=['Secondary structure', 'Band position in D2O, cm \u208b¹', 'Height, (a.u)',
                           'FWHM, cm\u208b¹'], max_col_width=50, def_col_width=15,
                 auto_size_columns=False,
                 # cols_justification=('left','center','right','c', 'l', 'bad'),       # Added on GitHub only as of June 2022
                 display_row_numbers=True,
                 justification='center',
                 num_rows=5,
                 alternating_row_color='dodgerblue',
                 key='-TABLE-',
                 selected_row_colors='red on yellow',
                 enable_events=True,
                 expand_x=False,
                 expand_y=True,
                 vertical_scroll_only=False)
    ],
]

# ----- Full layout -----
layout = [
    [
        sg.Column(file_list_column),
        sg.VSeperator(),
        sg.Column(image_viewer_column),
    ]
]

window = sg.Window("Absorption spectrum", layout)

# Define properties of figures

plt.rcParams['font.size'] = 14
plt.rcParams['figure.dpi'] = 100
plt.rc('legend', fontsize=14)
fig_canvas_agg = None

# Run the Event Loop
while True:
    event, values = window.read()
    if event == "Exit" or event == sg.WIN_CLOSED:
        break
    # Folder name was filled in, make a list of files in the folder
    if event == "-FOLDER-":
        folder = values["-FOLDER-"]
        try:
            # Get list of files in folder
            file_list = os.listdir(folder)
        except:
            file_list = []

        fnames = [
            f
            for f in file_list
            if os.path.isfile(os.path.join(folder, f))
               and f.lower().endswith((".txt"))
        ]
        window["-FILE LIST-"].update(fnames)
    elif event == "-FILE LIST-":  # A file was chosen from the listbox
        try:
            filename = os.path.join(
                values["-FOLDER-"], values["-FILE LIST-"][0]
            )
            window["-TOUT-"].update(filename)

            with open(filename, "rt", encoding='utf-8') as f:
                text = f.read().splitlines()
            # popup_text(filename, text)
            datax, datay = import_data(text)
            if fig_canvas_agg:
                delete_figure_agg(fig_canvas_agg)

            make_plot(datax[1:-310], datay[1:-310])  # cut the data before 1000 cm-1 (detector is not reliable there)
            fig = plt.gcf()

            fig_canvas_agg = draw_figure(window['-DATA-'].TKCanvas, fig)


        except:
            pass

    elif event == 'I. Extract Amide I band':
        try:
            datax_amid1, datay_amid1 = extract_amid1_band(datax, datay)
            if fig_canvas_agg:
                delete_figure_agg(fig_canvas_agg)
            make_plot(datax_amid1, datay_amid1)  # cut the data before 1000 cm-1 (detector is not reliable there)
            fig = plt.gcf()

            fig_canvas_agg = draw_figure(window['-DATA-'].TKCanvas, fig)


        except:
            pass

    elif event == 'II. Background subtraction':
        try:
            datay_amid1_b = subtract_background(datax_amid1, datay_amid1)
            if fig_canvas_agg:
                delete_figure_agg(fig_canvas_agg)
            make_plot(datax_amid1, datay_amid1_b)
            fig = plt.gcf()

            fig_canvas_agg = draw_figure(window['-DATA-'].TKCanvas, fig)

        except:
            pass


    elif event == 'III. Fit':
        try:

            y_fit, parameters, error_bars = fit_gauss5(np.array(datax_amid1), np.array(datay_amid1_b),
                                                       bounds=((values['-A0MIN-'], values['-A1MIN-'], values['-A2MIN-'],
                                                                 values['-A3MIN-'], values['-A4MIN-'], values['-S0MIN-'],
                                                                 values['-S1MIN-'], values['-S2MIN-'], values['-S3MIN-'],
                                                                 values['-S4MIN-']),
                                                               (values['-A0MAX-'], values['-A1MAX-'], values['-A2MAX-'],
                                                                 values['-A3MAX-'], values['-A4MAX-'], values['-S0MAX-'],
                                                                 values['-S1MAX-'], values['-S2MAX-'], values['-S3MAX-'],
                                                                 values['-S4MAX-'])))
            if fig_canvas_agg:
                delete_figure_agg(fig_canvas_agg)

            make_plot(datax_amid1, y_fit, 'lin')
            make_plot(datax_amid1, datay_amid1_b)
            fig = plt.gcf()
            fig_canvas_agg = draw_figure(window['-DATA-'].TKCanvas, fig)
            window['-TABLE-'].update(values=display_fit_results(parameters, error_bars))


        except:
            sg.popup_error_with_traceback('Fit did not converge. Change boundary conditions.')


    elif event == '-SAVE-':
        try:
            file_name = sg.popup_get_file('Save file', save_as=True)
            save_to_file_as_txt(file_name, datax_amid1, datay_amid1_b) #datay_amid1_b
        except:
            pass

window.close()
