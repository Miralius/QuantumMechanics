import numpy as np
import matplotlib.pyplot as plt


# Обезразмеренный потенциал U(r)=-u0 exp(-r/a), где r = ax
def given_potential(x):
    return -np.exp(-x)


# Осциллятор для теста U(x) = (1/2)x²
def oscillator(x):
    return 0.5 * x ** 2


# Получаем вектор потенциалов на заданной области
def get_potential_function(a, b, n, function):
    h = (b - a) / (n - 1)
    return h, np.array([a + i * h for i in range(n)]), np.array([function(a + i * h) for i in range(n)])


# Находим u0 исходя из заданного коэффициента B. a — половина интервала решения задачи, b — коэффициент B
# h_reduced — редуцированная постоянная Планка, m1 — масса первой частицы (пусть ядро водорода)
# m2 — масса второй частицы (пусть электрона), m — их приведённая масса
# B = 2ma²u0/ħ² → u0 = Bħ²/(2ma²)
# def u0(a, b):
#    h_reduced = 1.054572 * 10 ** (-34)
#    m1 = 1.672622 * 10 ** (-27)
#    m2 = 9.109383 * 10 ** (-31)
#    m = m1 * m2 / (m1 + m2)
#    a *= 2
#    return b * h_reduced ** 2 / (2 * m * a ** 2)


# Поиск собственных векторов и значений для заданного потенциала и интервала энергии
def find_eigenvector_and_eigenvalue(e_min, e_max, h_energy, accuracy, b_coefficient, u, n, h):
    a = e_min
    vectors = np.array([])
    values = np.array([])
    indexes = np.array([])
    psi_temp = np.array([])
    discrepancy_temp = 1
    while a <= e_max:
        discrepancy, psi, m = shooting_method(b_coefficient, a, u, n, h)
        if discrepancy < accuracy:
            if discrepancy_temp > discrepancy:
                discrepancy_temp = discrepancy
                psi_temp = psi
        else:
            if len(psi_temp) != 0:
                values = np.append(values, a)
                indexes = np.append(indexes, m)
                if len(vectors) == 0:
                    vectors = psi_temp
                else:
                    vectors = np.vstack((vectors, psi_temp))
                psi_temp = np.array([])
        a += h_energy
    return vectors, values, indexes


# Ищем невязку и собственный вектор методом пристрелки
def shooting_method(b_coefficient, energy, u, n, h):
    # Ищем узел для сшивания решений
    m = n // 2
    i = 0
    while i < n - 2:
        crossing_i_plus_1 = energy - u[i + 2]
        crossing_i = energy - u[i + 1]
        if crossing_i_plus_1 * crossing_i < 0:
            m = i + 2
            break
        i += 1
    if i == n - 2:
        return 1, np.array([]), 0
    # Интегрируем уравнение Шрёдингера слева-направо
    hi = np.zeros(n)
    hi[1] = 9.999999E-10
    i = 0
    while i < m:
        second_derivative_hi_i = -b_coefficient * (energy - u[i]) * hi[i]
        hi[i + 2] = h ** 2 * second_derivative_hi_i + 2 * hi[i + 1] - hi[i]
        i += 1
    hi_left_i = hi[m]
    hi_left_i_minus_1 = hi[m - 1]
    # Интегрируем уравнение Шрёдингера справа-налево
    hi[n - 2] = 9.999999E-10
    i = n - 3
    while i >= m - 1:
        hi[i] = (2 * hi[i + 1] - hi[i + 2]) / (1 + h ** 2 * b_coefficient * (energy - u[i]))
        i -= 1
    hi_right_i = hi[m]
    hi_right_i_minus_1 = hi[m - 1]
    # Сшиваем решения
    c = hi_left_i / hi_right_i
    hi_right_i_minus_1 *= c
    i = m - 1
    while i < n:
        hi[i] *= c
        i += 1
    # Ищем максимальное по модулю решение и считаем невязку
    psi_max = np.amax(abs(hi))
    discrepancy = abs(hi_left_i_minus_1 - hi_right_i_minus_1) / psi_max
    # Нормируем решения
    a = 0
    for i in range(n):
        a += h * hi[i] ** 2
    a = a ** (-0.5)
    hi *= a
    # hi_left_i_minus_1 *= a
    # hi_right_i_minus_1 *= a
    return discrepancy, hi, m


# Отрисовка одномерных графиков
def plot(x, functions, names, y_label, common_label):
    xy = plt.subplot()
    xy.grid(which='major', color='k')
    xy.minorticks_on()
    xy.grid(which="minor", color='gray', linestyle=':')
    if len(names) == 1:
        xy.plot(x, functions, label=names[0])
    else:
        i = 0
        while i < len(functions):
            xy.plot(x, functions[i], label=names[i])
            i += 1
    xy.set_xlabel("x")
    xy.set_ylabel(y_label)
    plt.title(common_label)
    xy.legend()
    plt.show()
