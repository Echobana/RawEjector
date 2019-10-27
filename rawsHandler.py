#! python3
# -*- coding utf-8 -*-

'''Processes data recieved from supersponic ejector,
 presented in the form of an file'''

'''sensor_info document format:
    1st column: sensorName
        200 - p_02
        1 - p_1
        2 - p_2
        3 - p_3
        4 - p_4
        5 - p_5
        6 - p_6
        7 - p_7
        8 - p_8
        9 - p_9
        10 - p_10
        11 - p_11
        400 - p_04
        or - p_orifce (flow washer)
        100 - p_01
    2nd column: sensor status
        on - 1
        off - 0
    3rd column: lower limit of measurement (atm)
    4th column: upper limit of measurement (atm)
    5th column: location (mm)
    '''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import scipy.optimize as opt
import os


def find_max(xval, yval, n=4):
    array = np.polyfit(xval, yval, n)
    max_x = opt.fminbound(lambda x: -np.polyval(array, x), min(xval), max(xval))
    return max_x, np.polyval(array, max_x)


def atma_to_mpa(p):
    return 0.1 * p


def atmo_to_mpa(p):
    return 0.1 * p + 0.1


def area(d):
    return np.pi * d * d / 4 * pow(10, -6)


def mass_flow(mu, F, p, A, R, T):
    return mu * F * p * A / np.sqrt(R * T)


def eff(n, p4, p2, p1, k):
    return n * (pow(p4 / p2, (k - 1) / k) - 1) / (1 - pow(p4 / p1, (k - 1) / k))


def opener(file_path):
    f = open(file_path, 'r')
    df = f.readlines()

    for i in range(len(df)):
        df[i] = df[i].replace(',', '.')
        df[i] = df[i].split("\t")
        for j in range(len(df[i])):
            try:
                df[i][j] = float(df[i][j])
            except ValueError:
                df[i][j] = df[i][j]
    df = np.array(df, dtype=object)
    return df.T


def verification(sensor_info, data):
    ver_data = dict()
    for i in range(len(sensor_info)):
        if sensor_info[i][1] == 1:
            flags = []
            for j in range(len(data[0])):
                if data[i][j] > sensor_info[i][2]:
                    if data[i][j] < sensor_info[i][3]:
                        flags.append(0)
            if 1 not in flags:
                ver_data.setdefault(sensor_info[i][0], data[i])
    ver_data.setdefault('marker', data[len(data) - 1])
    return ver_data


def av(verData_dict):
    array_list = []
    for v in verData_dict.values():
        array_list.append(v)
    data_matrix = np.array(array_list)
    '''was a function:
        return list of steady points locations (indexes):
        ((1, 2, 3), (100, 101, 102)) '''
    flags = []
    for rownum in range(len(data_matrix[0])):
        if data_matrix[len(data_matrix) - 1][rownum] == 1:
            flags.append(rownum)
    markers = []
    i = 0
    while flags:
        if i + 1 < len(flags):
            if flags[i + 1] - flags[i] > 1:
                markers.append(flags[:i + 1])
                del flags[:i + 1]
                i = 0
            else:
                i += 1
        else:
            markers.append(flags[:])
            flags = []
    # returns markers
    # sensor_dict = sensor_info(sensor_path)
    data_dict = dict()
    for k in verData_dict.keys():
        data_dict.setdefault(k, [])

    for marker in markers:
        for k in data_dict.keys():
            sum = 0
            for i in marker:
                sum += verData_dict[k][i]
            data_dict[k].append(sum / len(marker))
    return data_dict


def transient(verData_dict):
    first_std, = np.where(verData_dict['marker'] == 1)
    first_std_id = first_std[0]
    trans = dict()
    for colname in verData_dict.keys():
        trans.setdefault(colname, [])
        for rownum in range(first_std_id, len(verData_dict[colname])):
            if verData_dict['marker'][rownum] == 0:
                trans[colname].append(verData_dict[colname][rownum])
    return trans


def solver(data_dict):
    ej_coeff_list = []
    comp_ratio_list = []
    eff_list = []
    p_or = list(map(atma_to_mpa, data_dict['or']))
    p_01 = list(map(atma_to_mpa, data_dict[100]))
    p_02 = list(map(atma_to_mpa, data_dict[200]))
    p_04 = list(map(atma_to_mpa, data_dict[400]))
    for i in range(len(data_dict['or'])):
        m2_cur = mass_flow(mu_orifice, area(d_orifice), p_or[i], A_k, R_air, T)
        m1_cur = mass_flow(mu_nozzle, area(d_cr), p_01[i], A_k, R_air, T)
        ej_coeff_cur = m2_cur / m1_cur
        comp_ratio_cur = p_04[i] / p_02[i]
        eff_cur = eff(ej_coeff_cur, p_04[i], p_02[i], p_01[i], k)
        ej_coeff_list.append(ej_coeff_cur)
        comp_ratio_list.append(comp_ratio_cur)
        eff_list.append(eff_cur)
    return p_or, p_01, p_04, p_02, ej_coeff_list, comp_ratio_list, eff_list


def mult_plot(graph):
    for k, v in graph.items():
        fig, ax = plt.subplots()
        # if k in maximize:
        #     max_arg, max_value = find_max(v[0], v[1])
        #     ax.scatter(max_arg, max_value, color='purple', s=80)
        plt.plot(v[0], v[1], color='black')
        plt.scatter(v[2], v[3], color='black', marker='x', s=13, label='Нестационарные точки')
        plt.scatter(v[0], v[1], color='red', marker='o', label='Осредненные стационарные точки')
        plt.xlabel(k[0])
        plt.ylabel(k[1], rotation=90)
        plt.title(k[2])
        plt.grid()
        for i in range(len(v[0])):
            plt.text(v[0][i] - 0.01 * (max(v[0]) - min(v[0])), v[1][i] + 0.06 * (max(v[1]) - min(v[1])), str(i + 1))


def length_plot(sensor_info, data_dict):
    sensor_info = sensor_info.T
    distances = dict()
    for i in range(len(sensor_info)):
        if sensor_info[i][1] == 1:
            if sensor_info[i][4] != 'None':
                distances.setdefault(sensor_info[i][0], sensor_info[i][4])
    fig_l, ax_l = plt.subplots()
    for i in range(len(data_dict[400])):
        x_val, y_val = [], []
        for k, v in data_dict.items():
            if k not in ('or', 100, 'marker') and k in (200, 400):
                x_val.append(distances[k])
                y_val.append(atma_to_mpa(data_dict[k][i]))
                # ax_l.text(distances[k], atma_to_mpa(data_dict[k][i]), 's')
            elif k not in ('or', 100, 'marker'):
                x_val.append(distances[k])
                y_val.append(atma_to_mpa(data_dict[k][i]))
                # ax_l.text(distances[k], atmo_to_mpa(data_dict[k][i]), 's')
        plt.plot(x_val, y_val)
        plt.scatter(x_val, y_val, marker=i+1)
    plt.grid()
    plt.xlabel('L, мм')
    plt.ylabel("p, МПа")


def pdf_saver(path):
    pdf = matplotlib.backends.backend_pdf.PdfPages(path)
    for fig in range(1, plt.figure().number):
        pdf.savefig(fig)
    pdf.close()


d_orifice = 3.03
d_cr = 4.28
mu_orifice = 0.98
mu_nozzle = 0.98
R_air = 287
A_k = 0.685
k = 1.4
T = 298

if __name__ == '__main__':
    data_file = os.path.abspath(r"F:\ejector_raw_files\2")
    sensor_file = os.path.abspath(r"F:\ejector_raw_files\sensor_status")

    data = opener(data_file)
    sensor_info = opener(sensor_file)

    ver_data = verification(sensor_info.T, data)

    data_dict = av(ver_data)
    p_or, p_01, p_04, p_02, ej_coeff_list, comp_ratio_list, eff_list = solver(data_dict)

    trans = transient(ver_data)
    p_or_tns, p_01_tns, p_04_tns, p_02_tns, ej_coeff_list_tns, comp_ratio_list_tns, eff_list_tns = solver(trans)

    graph = {
        (r'$p_{04}$, МПа', r'Коэффициент эжекции, k', 'Зависимость коэффициента эжекции от давления'): (p_04, ej_coeff_list, p_04_tns, ej_coeff_list_tns),
        (r'$p_{04}$, МПа', r'$p_{02}$, МПа', 'Title_2'): (p_04, p_02, p_04_tns, p_02_tns),
        (r'$p_{04}$, МПа', r'$\varepsilon$', 'Title_3'): (p_04, comp_ratio_list, p_04_tns, comp_ratio_list_tns),
        (r'$p_{04}$, МПа', r'$\eta$', 'Title_4'): (p_04, eff_list, p_04_tns, eff_list_tns)
    }  # Оси и легенды графиков 1-4

    maximize = ((r'$p_{04}$, МПа', r'$\varepsilon$', 'Title_3'), (r'$p_{04}$, МПа', r'$\eta$', 'Title_4'))

    mult_plot(graph)
    length_plot(sensor_info, data_dict)
    pdf_saver(r'./result_plot.pdf')
