from functions import *

# Несколько вводных слов. В docx файле содержится весь алгоритм, будет полезно почитать перед тем, как лезть в код,
# чтоб понять, что к чему. При B = 10 код должен выдать на нулевом уровне ε = -0.219 (у Салеева -0.218), на первом —
# ε ≈ -0.006. Далее Салеев просил проверить работоспособность кода при B = 15. На нулевом уровне я получил ε = -0.286
# (а у Салеева -0.2845, затем первый уровень пытался найти, но не успел, он мне зачел задачу. Так чтоо, вот. В коде
# могут быть ошибки, так что наезжать на меня не надо по поводу того, что шо-то ниправииильнааа считается!
# © Miralius, 12/12/2020

if __name__ == '__main__':
    # Зададим область решения задачи и кол-во точек интегрирования. Возьмем для начала гармонический осциллятор
    a = 0.001
    b = 10
    n = 1000
    h, x, u = get_potential_function(a, b, n, oscillator)
    # Задаем границы собственного значение энергии, точность
    e_min = 0
    e_max = 2
    h_energy = 5 * 10 ** (-3)
    accuracy = 5 * 10 ** (-5)
    b_coefficient = 2
    psi, energy, m = find_eigenvector_and_eigenvalue(e_min, e_max, h_energy, accuracy, b_coefficient, u, n, h)
    potentials = np.array([u])
    labels_energy, labels_number = np.array("F(x) = (1/2)x²"), np.array("F(x) = (1/2)x²")
    for i in range(len(energy)):
        print("Осциллятор, n = " + str(i) + ", собственное значение энергии = " + str("%.3f" % energy[i]))
        print("Решение сшито в узле № " + str(int(m[i])) + ", x = " + str("%.3f" % x[int(m[i])]))
        potentials = np.vstack((potentials, np.array([energy[i] for j in range(n)])))
        labels_energy = np.append(labels_energy, "ε = " + str("%.3f" % energy[i]))
        labels_number = np.append(labels_number, "n = " + str(i))
    plot(x, potentials, labels_energy, "U(x)", "График потенциала")
    plot(x, psi, np.array([labels_number[i + 1] for i in range(len(energy))]), "Ψ(x)", "Графики волновых функций")
    plot(x, np.vstack((u, psi)), labels_number, "U(x) | Ψ(x)", "Графики потенциала и волновых функций")

    # Теперь берем потенциал согласно варианту
    b = 10
    h_energy = 10 ** (-5)
    accuracy = 5 * 10 ** (-5)
    b_coefficient = 15
    h, x, u = get_potential_function(a, b, n, given_potential)
    # Задаем границы собственного значение энергии, точность
    e_min = -1
    e_max = 0
    psi, energy, m = find_eigenvector_and_eigenvalue(e_min, e_max, h_energy, accuracy, b_coefficient, u, n, h)
    potentials = np.array([u])
    labels_energy, labels_number = np.array("F(x) = -exp(-x)"), np.array("F(x) = -exp(-x)")
    for i in range(len(energy)):
        print("F(x) = -exp(-x), n = " + str(i) + ", собственное значение энергии = " + str("%.5f" % energy[i]))
        print("Решение сшито в узле № " + str(int(m[i])) + ", x = " + str("%.3f" % x[int(m[i])]))
        potentials = np.vstack((potentials, np.array([energy[i] for j in range(n)])))
        labels_energy = np.append(labels_energy, "ε = " + str("%.5f" % energy[i]))
        labels_number = np.append(labels_number, "n = " + str(i))
    plot(x, potentials, labels_energy, "U(x)", "График потенциала")
    plot(x, psi, np.array([labels_number[i + 1] for i in range(len(energy))]), "Ψ(x)", "Графики волновых функций")
    plot(x, np.vstack((u, psi)), labels_number, "U(x) | Ψ(x)", "Графики потенциала и волновых функций")
