import numpy as np

def generate_m_s_c():
    m = 10
    n = 10
    m1 = np.array([[1] * m] * n)
    for i in range(m):
        for j in range(n):
            if np.random.uniform(0, 1) > 0.3:
                m1[i][j] = 0

    print(m1)
    # row names/samples
    s1 = []
    for i in range(1, len(m1)+1):
        si = 's' + str(i)
        s1.append(si)
    print(s1)
    s1 = np.array(s1, dtype='S100')

    # column names/features
    c1 = []
    for i in range(1, len(m1)+1):
        ci = 'c' + str(i)
        c1.append(ci)
    print(c1)
    c1 = np.array(c1, dtype='S100')

    return m1, s1, c1

generate_m_s_c()