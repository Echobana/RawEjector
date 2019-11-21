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
from conversion import RateGraph as RG
import sensor_class


def find_max(xval, yval, n=4):
    array = np.polyfit(xval, yval, n)
    max_x = opt.fminbound(lambda x: -np.polyval(array, x), min(xval), max(xval))
    return max_x, np.polyval(array, max_x)


def atma_to_mpa(p):
    return (0.0980665 * p) * 1000


def atmo_to_mpa(p):
    return (0.0980665 * p + 0.0980665) * 1000


def area(d):
    return np.pi * d * d / 4 * pow(10, -6)


def converter(magnitude, unit_of_measurment):
    rates = [['atma', 'atmg', (1, 1)],
             ['atma', 'MPa', (0.1, 0)],
             ['atmg', 'MPa', (0.1, 0.1)]]
    graph = RG(rates)
    coefficients = graph.get_conversion(unit_of_measurment, 'kPa')
    return coefficients[0] * (magnitude + coefficients[1])


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


def experiment_parameters(m1_list, m1_list_tns, m2_list, m2_list_tns, p_01, p_01_tns, ej_coeff_list, params):
    fig, ax = plt.subplots()
    ax.axis('off')
    m1 = concat2(m1_list, m1_list_tns)
    m2 = concat2(m2_list, m2_list_tns)
    p_01l = concat2(p_01, p_01_tns)
    m1_min, m1_max = min(m1), max(m1)
    m2_min, m2_max = min(m2), max(m2)
    p_01_min, p_01_max = min(p_01l), max(p_01l)
    text = 'Сопло: \n' \
           '$d_{{кр}}$={0:.2f} мм \n' \
           '$d_{{a}}$={1:.2f} мм\n' \
           '${{\\alpha}}^o$={2:.2f}${{^o}}$ \n' \
           '${{\\mu_{{c}}}}$={3:.2f} \n' \
           '$D_{{кс}}$={4:.2f} мм\n' \
           '$L_{{кс}}$={5:.2f} мм\n' \
           'Диффузор: \n' \
           '${{\\alpha}}_{{вых}}^o$={6:.2f} \n' \
           '$D_{{вых}}$={7:.2f} мм\n\n' \
           'Расходомерная шайба: \n' \
           '$d_{{ш}}$={8:.2f} мм\n' \
           '${{\\mu}}_{{c}}$={9:.2f} \n\n' \
           'Параметры эксперимента: \n\n' \
           'Активный газ - воздух, $T_{{01}}$={10:.2f} \n' \
           'Пассивный газ - воздух, $T_{{02}}$={11:.2f} \n' \
           '$p_{{01}}$={12:.2f}...{13:.2f} кПа \n' \
           '$m_{{1}}$={14:.2f}...{15:.2f} г/c \n' \
           '$m_{{2}}$={16:.2f}...{17:.2f} г/c \n' \
           '$k_{{cр}}$={18:.2f} \n\n' \
        .format(params[0],  # d_cr
                params[1],  # d_a
                params[2],  # alpha
                params[3],  # mu_nozzle
                params[4],  # D_cc
                params[5],  # L_cc
                params[6],  # alpha_outlet
                params[7],  # D_outlet
                params[8],  # d_orifice
                params[9],  # mu_orifice
                params[10],  # T_01
                params[11],  # T_02
                p_01_min,  # p_left
                p_01_max,  # p_right
                1000 * m1_min,  # m1_left
                1000 * m1_max,  # m1_right
                1000 * m2_min,  # m2_left
                1000 * m2_max,  # m2_right
                ej_coeff_report(ej_coeff_list)
                )
    plt.text(0., 1., text, ha='left', va='top', transform=ax.transAxes, fontsize='11')


def solver(data_dict, params):
    k = 1.4
    A_k = 0.685
    R_air = 287

    ej_coeff_list = []
    comp_ratio_list = []
    eff_list = []
    m1_list = []
    m2_list = []
    p_or = list(map(atma_to_mpa, data_dict['or']))
    p_01 = list(map(atma_to_mpa, data_dict[100]))
    p_02 = list(map(atma_to_mpa, data_dict[200]))
    p_04 = list(map(atma_to_mpa, data_dict[400]))
    for i in range(len(data_dict['or'])):
        m2_cur = mass_flow(params[9], area(params[8]), 1000 * p_or[i], A_k, R_air, params[10])
        m1_cur = mass_flow(params[3], area(params[0]), 1000 * p_01[i], A_k, R_air, params[10])
        ej_coeff_cur = m2_cur / m1_cur
        comp_ratio_cur = p_04[i] / p_02[i]
        eff_cur = eff(ej_coeff_cur, p_04[i], p_02[i], p_01[i], k)
        ej_coeff_list.append(ej_coeff_cur)
        comp_ratio_list.append(comp_ratio_cur)
        eff_list.append(eff_cur)
        m1_list.append(m1_cur)
        m2_list.append(m2_cur)
    return p_or, p_01, p_04, p_02, ej_coeff_list, comp_ratio_list, eff_list, m1_list, m2_list


def ej_coeff_report(ej_coeff_list):
    sum = 0
    for i in range(len(ej_coeff_list)):
        sum += ej_coeff_list[i]
    return sum / len(ej_coeff_list)


def comp_ratio_report(comp_ratio_list):
    return max(comp_ratio_list)


def plot_ej(data_list, title_list):
    default_plot(data_list, title_list)


def plot_p04comp(data_list, title_list):
    max_variable, max_value = data_max(data_list)
    text = r'$\varepsilon_{{max}}$={0:.2f}'.format(
        max_value) + '\n' + r'$p_{{04}}^{{\varepsilon_{{max}}}}$={0:.2f} кПа'.format(
        max_variable)
    default_plot(data_list, title_list, text=text)


def plot_p04eff(data_list, title_list):
    max_variable, max_value = data_max(data_list)
    text = r'$\eta_{{max}}$={0:.2f}'.format(
        max_value) + '\n' + r'$p_{{04}}^{{\eta_{{max}}}}$={0:.2f} кПа'.format(
        max_variable)
    default_plot(data_list, title_list, text=text)


def data_max(data_list):
    x_data = data_list[0] + data_list[2]
    y_data = data_list[1] + data_list[3]
    max_value = max(y_data)
    max_value_index = y_data.index(max_value)
    max_variable = x_data[max_value_index]
    return max_variable, max_value


def default_plot(data_list, title_list, text=''):
    fig, ax = plt.subplots()
    plt.scatter(data_list[2], data_list[3], marker='x', s=13, color='black')
    plt.scatter(data_list[0], data_list[1], marker='o', color='red')
    plt.xlabel(title_list[0])
    plt.ylabel(title_list[1])
    plt.grid()
    for i in range(len(data_list[0])):
        plt.text(data_list[0][i] - 0.01 * (max(data_list[0]) - min(data_list[0])),
                 data_list[1][i] + 0.06 * (max(data_list[1]) - min(data_list[1])), str(i + 1))
    props = dict(boxstyle='round', facecolor='wheat')
    plt.text(0.04, 0.9, text, ha='left', va='center', transform=ax.transAxes, bbox=props, fontsize='11')


def plot_p04p02(data_list, title_list):
    default_plot(data_list, title_list)


def concat2(list_1, list_2):
    return list_1 + list_2


def set_experiment_parameters(parameters_file):
    parameters = opener(parameters_file)
    params = []
    for line in parameters:
        params.append(line[1])
    # d_cr = params[0]
    # d_a = params[1]
    # alpha = params[2]
    # mu_nozzle = params[3]
    # d_cc = params[4]
    # l_cc = params[5]
    # alpha_2 = params[6]
    # d_out = params[7]
    # d_orifice = params[8]
    # mu_orifice = params[9]
    # T_01 = params[10]
    # T_02 = params[11]
    # k = 1.4
    # A_k = 0.685
    # R_air = 287
    return params


def mult_plot(graph, sensor_info, data_dict):
    plot_p04p02(graph[(r'$p_{04}$, кПа', r'$p_{02}$, кПа', 'Title_2')], (r'$p_{04}$, кПа', r'$p_{02}$, кПа', 'Title_2'))
    plot_p04comp(graph[(r'$p_{04}$, кПа', r'$\varepsilon$', 'Title_3')],
                 (r'$p_{04}$, кПа', r'$\varepsilon$', 'Title_3'))
    plot_p04eff(graph[(r'$p_{04}$, кПа', r'$\eta$', 'Title_4')], (r'$p_{04}$, кПа', r'$\eta$', 'Title_4'))
    length_plot(sensor_info, data_dict)


def throttle_plot(graph, sensor_info, data_dict):
    plot_ej(graph[(r'$p_{04}$, кПа', r'Коэффициент эжекции, k', 'Зависимость коэффициента эжекции от давления')],
            (r'$p_{04}$, кПа', r'Коэффициент эжекции, k', 'Зависимость коэффициента эжекции от давления'))
    length_plot(sensor_info, data_dict)


def length_plot(sensor_info, data_dict):
    sensor_info = sensor_info.T
    distances = dict()
    for i in range(len(sensor_info)):
        if sensor_info[i][1] == 1:
            if sensor_info[i][4] != 'None':
                distances.setdefault(sensor_info[i][0], sensor_info[i][4])
    plt.figure(figsize=(20, 10))
    for i in range(len(data_dict[400])):
        x_val, y_val = [], []
        for k, v in data_dict.items():
            if k not in ('or', 100, 'marker') and k in (200, 400):
                x_val.append(distances[k])
                y_val.append(atma_to_mpa(data_dict[k][i]))
            elif k not in ('or', 100, 'marker'):
                x_val.append(distances[k])
                y_val.append(atma_to_mpa(data_dict[k][i]))
        if i % 2 == 0:
            plt.text(1.01 * x_val[-1], 0.99 * y_val[-1], str(i + 1))
        else:
            plt.text(1.02 * x_val[-1], 0.99 * y_val[-1], str(i + 1))
        plt.scatter(x_val, y_val, marker=i + 1)
        plt.plot(x_val, y_val)
    plt.grid()
    plt.xlabel('L, мм')
    plt.ylabel("p, кПа")


def pdf_saver(path):
    pdf = matplotlib.backends.backend_pdf.PdfPages(path)
    for fig in range(1, plt.figure().number):
        pdf.savefig(fig)
    pdf.close()


if __name__ == '__main__':
    parameters_file = r'./parameters'
    data_file = os.path.abspath(r"F:\ejector_raw_files\2")
    sensor_file = os.path.abspath(r"F:\ejector_raw_files\sensor_status")

    data = opener(data_file)
    sensor_info = opener(sensor_file)
    params = set_experiment_parameters(
        parameters_file)

    ver_data = verification(sensor_info.T, data)

    data_dict = av(ver_data)
    p_or, p_01, p_04, p_02, ej_coeff_list, comp_ratio_list, eff_list, m1_list, m2_list = solver(data_dict, params)

    trans = transient(ver_data)

    p_or_tns, p_01_tns, p_04_tns, p_02_tns, \
    ej_coeff_list_tns, comp_ratio_list_tns, eff_list_tns, m1_list_tns, m2_list_tns = solver(
        trans, params)

    graph = {
        (r'$p_{04}$, кПа', r'Коэффициент эжекции, k', 'Зависимость коэффициента эжекции от давления'): (
            ej_coeff_list, p_02, ej_coeff_list_tns, p_02_tns),
        (r'$p_{04}$, кПа', r'$p_{02}$, кПа', 'Title_2'): (p_04, p_02, p_04_tns, p_02_tns),
        (r'$p_{04}$, кПа', r'$\varepsilon$', 'Title_3'): (p_04, comp_ratio_list, p_04_tns, comp_ratio_list_tns),
        (r'$p_{04}$, кПа', r'$\eta$', 'Title_4'): (p_04, eff_list, p_04_tns, eff_list_tns)
    }  # Оси и легенды графиков 1-4

    maximize = ((r'$p_{04}$, МПа', r'$\varepsilon$', 'Title_3'), (r'$p_{04}$, МПа', r'$\eta$', 'Title_4'))

    checkbox = 'throttle'
    if checkbox == 'throttle':
        experiment_parameters(m1_list, m1_list_tns, m2_list, m2_list_tns, p_01, p_01_tns, ej_coeff_list, params)
        throttle_plot(graph, sensor_info, data_dict)
    elif checkbox == 'mass_flow':
        experiment_parameters(m1_list, m1_list_tns, m2_list, m2_list_tns, p_01, p_01_tns, ej_coeff_list, params)
        mult_plot(graph, sensor_info, data_dict)
    else:
        print('did not found')

    # plt.show()
    pdf_saver(r'./result_plot.pdf')
