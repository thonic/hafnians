from sympy import symbols, expand, Matrix, prod
import numpy as np

N = 5
R = 3

# Create symbolic variables x1, x2, ..., xN
x = symbols("x1:{}".format(N + 1))
xv = Matrix(x)

# Create symbolic variables y1, y2, ..., yR
y = symbols("y1:{}".format(R + 1))
yv = Matrix(y)

# random N by N matrix with entries -4..4
t = np.random.randint(-4, 5, size=(N, R))

# apply the random transformation
xv = t * yv

# the big product of partial derivatives
polynomial = prod(xv)

# Expand the polynomial
expanded_polynomial = expand(polynomial)

# Get the coefficients and exponents
coefficients = expanded_polynomial.as_coefficients_dict()
terms = list(coefficients.keys())

for term in terms:
    print("Coefficient:", coefficients[term])
    print("Exponents:", term.as_coefficients_dict())
    print()
