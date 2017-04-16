
def zero_matrix(h, w):
    return [[0 for c in range(w)] for r in range(h)]


def identity_matrix(h, w):
    z = zero_matrix(h, w)
    for i in range(h):
        z[i][i] = 1
    return z


def get(m, i ,j):
    return m[i-1][j-1]


def set_(m, i, j, x):
    m[i-1][j-1] = x
    return m


# assuming matrix is consistent
def dim(m):
    return (len(m), len(m[0]))


def transpose(m):
    h, w = dim(m)
    new_m = zero_matrix(w, h)

    for i in range(h):
        for j in range(w):
            set_(new_m, j+1, i+1, get(m, i+1, j+1))

    return new_m


def dot_prod(u, v):
    len_u = len(u)
    len_v = len(v)

    if len_u != len_v:
        raise ValueError('dimensions do not match')

    return sum([u[i] * v[i] for i in range(len_u)])


def magnitude(v):
    return sum([pow(c, 2) for c in v])**0.5


def vector_scale(v, s):
    return [c * s for c in v]


def unit_vector(v):
    mag = magnitude(v)
    return [c / mag for c in v]


def scalar_mult(m, s):
    return [[s * c for c in r] for r in m]


def sub_matrix(m, i, j):
    new_m = []
    new_new_m = []

    for x in range(0, len(m)):
        if x != i-1:
            new_m.append(m[x])

    trans_new_m = transpose(new_m)

    for y in range(0, len(trans_new_m)):
        if y != j-1:
            new_new_m.append(trans_new_m[y])

    return transpose(new_new_m)


def matrix_addition(a, b):

    len_a = len(a)
    len_b = len(b)

    if len_a != len_b:
        raise ValueError('dimensions do not match')
    for r in range(len_a):
        if len(a[r]) != len(b[r]):
            raise ValueError('dimensions do not match')

    res = []
    for r in range(len_a):
        new_row = []
        row1 = a[r]
        row2 = b[r]
        for i in range(len(row1)):
            new_row.append(row1[i] + row2[i])
        res.append(new_row)

    return res


def matrix_mult(a, b):
    a_h, a_w = dim(a)
    b_h, b_w = dim(b)

    if a_w != b_h:
        raise ValueError('dimensions do not match')

    res = zero_matrix(a_h, b_w)
    trans_b = transpose(b)
    for i in range(a_h):
        row = a[i]
        for j in range(len(trans_b)):
            row2 = trans_b[j]
            res[i][j] = dot_prod(row, row2)

    return res


def matrix_pow(m, p):
    h, w = dim(m)
    if h != w:
        raise ValueError('non-square matrix')
    if p == 1:
        return m
    if p % 2 == 0:
        return matrix_pow(matrix_mult(m,m), p/2)
    else:
        return matrix_mult(matrix_pow(m, p-1), m)


def append_col(m, c):
    for i in range(len(m)):
        m[i].append(c[i])
    return m


def append_row(m, r):
    m.append(r)
    return m


def augment(m1, m2):
    for r in transpose(m2):
        m1 = append_col(m1, r)
    return m1


# O(n!)
def laplace_det(m):
    h, w = dim(m)
    if h != w:
        raise ValueError('non square matrix')
    if dim(m) == (2, 2):
        return get(m, 1, 1) * get(m, 2, 2) - get(m, 1, 2) * get(m, 2, 1)
    else:
        return sum([get(m, 1, i + 1) * laplace_det(sub_matrix(m, 1, i + 1)) * (-1)**(i) for i in range(len(m))])


# O(n^3)
def det(m):
    h, w = dim(m)
    if h != w:
        raise ValueError('non square matrix')
    m = row_echelon(m)
    det = 1
    for i in range(h):
        det *= m[i][i]

    return det


def rank_null(m):
    h, w = dim(m)
    m = row_echelon(m)
    rank = 0
    for r in range(h):
        for c in range(w):
            if m[r][c] != 0:
                rank += 1
                break

    return rank, w-rank


def row_echelon(m):
    h, w = dim(m)

    for i in range(h):
        # print prettify(m)
        for j in range(i + 1, h):
            if m[i][i] != 0:
                r = float(m[j][i]) / m[i][i]
                for c in range(w):
                    m[j][c] -= m[i][c] * r
            else:
                #find pivot
                for k in range(i, h):
                    if m[k][i] != 0:
                        tmp = m[k]
                        m[k] = m[i]
                        m[i] = tmp
    return m


#BUGS:
    # doesn't work for h > w
def rre(m):
    h, w = dim(m)

    m = row_echelon(m)

    for i in range(h-1, -1, -1):
        for j in range(i+1, h):
            ratio = float(m[i][j]) / m[i+1][i+1]
            for k in range(w):
                m[i][k] -= ratio * m[j][k]


        for j in range(w):
            m[i][j] *= 1.0/m[i][i]

    #scale down
    for r in range(min(h,w)):
        for c in range(w):
            if m[r][r] != 0:
                m[r] = vector_scale(m[r], 1.0/m[r][r])

    return m


# returns x in Ax = b
def solve(A, b):
    return matrix_mult(inverse(A), b)


# do checks?
def inverse(m):
    h, w = dim(m)
    return transpose(transpose(rre(augment(m, identity_matrix(h, w))))[h:])


# cannot prettify vector
# wide matrix has issues
def prettify(m):
    res = ''
    largest = 0

    # find largest item
    for r in m:
        c = 0
        for i in r:
            str_i = str(i)
            c = len(str_i) * 2
        if c > largest:
            largest = c

    # print spaces
    for r in range(len(m)):
        s = ''
        for i in m[r]:
            str_i = str(i)
            s += str_i
            for x in range(largest - len(str_i)):
                s += ' '
        res = res + s + '\n'

    return res


if __name__ == '__main__':

    m = [[1,3],
         [5,6],
         [4,3]]

    m2 = [[1, 3, 1, 1],
          [2, 3, 5, 2],
          [5, 2, 4, 5]]

    m3 = [[1, 5, 6, 7, 1],
          [4, 3, 2, 8, 2],
          [4, 2, 5, 3, 6]]

    m1 = [[1, 3, 4, 12],
          [2, 3, 5, 54],
          [3, 10, 4, 3],
          [2, 54, 67, 9]]

    # m4 = [[random.randint(0, 10) for i in range(5)] for j in range(5)]

    # print prettify(m)
    # print 'inverse:'
    # print prettify(inverse(m))

    print prettify(rre(m))
