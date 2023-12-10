from methods_AG1817 import *

# problem 2.1
latitude = 60
longitude = 30
height = 0

# a)
print("\n", "a)")
x, y, z = geodetic_to_geocentric(GRS80, latitude, longitude, height)
print("(x, y, z) = ", (x, y, z), "\n")

# b)
print("b)")
r, φ1, λ = rectangular_to_spherical(x, y, z)
print(r, φ1, λ, "\n")

# c)
print("c)")
a, rf = GRS80  # semi-major axis, reciprocal flattening
e2 = 1 - (1 - 1 / rf) ** 2
beta = degrees(atan(tan(radians(φ1)) / sqrt(1 - e2)))
print(beta, "\n")

# d)
print("d)")
φ = radians(latitude)
d1 = degrees(φ) - φ1
d2 = degrees(φ) - beta
print(d1, "\n", d2, "\n")

# uppg 2.2
print("2.2)")
x, y, z = 3042053.7355, 988423.1757, 5503075.2100
ellipsoid = GRS80

(latitude, longitude, height) = geocentric_to_geodetic(ellipsoid, x, y, z)
print(
    "latitude: ",
    latitude,
    "\n",
    "longitude: ",
    longitude,
    "\n",
    "height: ",
    height,
    "\n",
)

# uppg 2.3
print("2.3)")
ellipsoid = GRS80
r, φ, λ = iterative_cartesian_to_geodetic(ellipsoid, x, y, z)
print("r: ", r, "\n", "φ: ", φ, "\n", "λ: ", λ, "\n")

# 2.4
print("2.4)")
# ITRF97
print("ITRF97:")
P = np.array([[3370641.970], [711866.128], [5349796.160]])
print(P, "\n")

# GRS80
print("GRS80:")
grs = np.array([[57.395483300], [11.925395149], [43.372]])
print(grs, "\n")
φ = radians(grs[0, 0])
λ = radians(grs[1, 0])
h = radians(grs[2, 0])

# a)
print("a)")
Q = np.array(
    [
        [-sin(φ) * cos(λ), -sin(φ) * sin(λ), cos(φ)],
        [-sin(λ), cos(λ), 0],
        [cos(φ) * cos(λ), cos(φ) * sin(λ), sin(φ)],
    ]
)
print("Q is: ", "\n", Q, "\n")

# b)
print("b)")
A = np.array([[3098889.388], [1011032.696], [5463980.133]])
D = Q * (A - P)
print("point A in P-NEU: ", "\n", D, "\n")

# c)
print("c)")
vel_1 = np.array([[-0.0136], [0.0147], [0.0084]])
vel_2 = Q * vel_1
print(
    "velocity vector: ",
    "\n",
    vel_2,
    "\n",
)

# d)
print("d)")
v = sqrt((vel_2[0, 0] ** 2) + (vel_2[1, 0] ** 2))
print("velocity: ", v)
direction = degrees(atan(vel_2[1, 0] / vel_2[0, 0]))
print("direction: ", direction, "\n")
