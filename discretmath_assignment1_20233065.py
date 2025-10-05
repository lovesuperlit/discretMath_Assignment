def input_matrix() :
    n = int(input("n을 입력 : "))
    print(f"{n} x {n} 행렬의 각 행을 입력 :")
    matrix = []
    for i in range(n) :
        row = list(map(float, input(f"{i+1}번째 행 입력 : ").split()))
        if len(row) != n :
            raise ValueError("Error")
        matrix.append(row)
    return matrix

def determinant(matrix) :
    n = len(matrix)
    if n == 1 :
        return matrix[0][0]
    if n == 2 :
        return matrix[0][0]*matrix[1][1] - matrix[0][1]*matrix[1][0]

    det = 0
    for c in range(n) :
        minor = [row[:c] + row[c+1:] for row in matrix[1:]]
        det += ((-1) ** c) * matrix[0][c] * determinant(minor)
    return det


def adjugate(matrix) :
    n = len(matrix)
    adj = []
    for i in range(n) :
        adj_row = []
        for j in range(n) :
            minor = [row[:j] + row[j+1:] for k, row in enumerate(matrix)
            if k != i]
            cofactor = ((-1) ** (i + j)) * determinant(minor)
            adj_row.append(cofactor)
        adj.append(adj_row)
    adj = [list(row) for row in zip(*adj)]
    return adj


def inverse_determinant(matrix) :
    det = determinant(matrix)
    if det == 0 :
        raise ValueError("Error")
    adj = adjugate(matrix)
    n = len(matrix)
    inv = [[adj[i][j] / det for j in range(n)] for i in range(n)]
    return inv

def inverse_gaussian_elimination(matrix) :
    n = len(matrix)
    identity = [[float(i == j) for j in range(n)] for i in range(n)]
    aug = [matrix[i] + identity[i] for i in range(n)]

    for i in range(n) :
        if aug[i][i] == 0 :
            for k in range(i + 1, n) :
                if aug[k][i] != 0 :
                    aug[i], aug[k] = aug[k], aug[i]
                    break
            else :
                raise ValueError("Error")

        pivot = aug[i][i]
        for j in range(i, 2 * n) :
            aug[i][j] /= pivot

        for k in range(i + 1, n) :
            factor = aug[k][i]
            for j in range(i, 2 * n) :
                aug[k][j] -= factor * aug[i][j]

    for i in range(n - 1, -1, -1) :
        for k in range(i - 1, -1, -1) :
            factor = aug[k][i]
            for j in range(i, 2 * n) :
                aug[k][j] -= factor * aug[i][j]

    inverse = [row[n:] for row in aug]
    return inverse

def inverse_gauss_jordan(matrix) :
    n = len(matrix)
    identity = [[float(i == j) for j in range(n)] for i in range(n)]
    aug = [matrix[i] + identity[i] for i in range(n)]

    for i in range(n) :
        if aug[i][i] == 0 :
            for k in range(i+1, n) :
                if aug[k][i] != 0 :
                    aug[i], aug[k] = aug[k], aug[i]
                    break
            else:
                raise ValueError("Error")

        pivot = aug[i][i]
        aug[i] = [x / pivot for x in aug[i]]

        for j in range(n) :
            if j != i :
                factor = aug[j][i]
                aug[j] = [aug[j][k] - factor * aug[i][k] for k in range(2 * n)]

    inverse = [row[n:] for row in aug]
    return inverse

def print_matrix(matrix, name="Matrix") :
    print(f"\n{name}:")
    for row in matrix :
        print("  ".join(f"{val:8.4f}" for val in row))


def matrix_equal(A, B, tol=1e-6) :
    n = len(A)
    for i in range(n) :
        for j in range(n) :
            if abs(A[i][j] - B[i][j]) > tol :
                return False
    return True

if __name__ == "__main__" :
    try :
        A = input_matrix()

        inv_det = inverse_determinant(A)
        inv_gauss = inverse_gaussian_elimination([row[:] for row in A])
        inv_gj = inverse_gauss_jordan([row[:] for row in A])

        print_matrix(inv_det, "행렬식 이용 역행렬")
        print_matrix(inv_gauss, "가우스 소거법 이용 역행렬")
        print_matrix(inv_gj, "가우스-조던 소거법 이용 역행렬")

        same12 = matrix_equal(inv_det, inv_gauss)
        same13 = matrix_equal(inv_det, inv_gj)
        same23 = matrix_equal(inv_gauss, inv_gj)

        print("\n<결과>")
        print(f"행렬식과 가우스 소거법: {same12}")
        print(f"행렬식과 가우스-조던 소거법: {same13}")
        print(f"가우스 소거법과 가우스-조던 소거법: {same23}")

    except ValueError as Error :
        print("\n오류:", Error)
