from math import asin
from AG1817_Methods import *

# 3.1
print("3.1)", "\n")
latitude_p = 60
longitude_p = 18
ellipsoid = GRS80
a, rf = ellipsoid
e2 = 1 - (1 - 1 / rf) ** 2

print("a)")
# north --> 10m
# ds * cos(alpha) = M * d_latitude
alpha = 0
ds = 10
d_longitude = 0

M = (a * (1 - e2)) / (((1 - (e2 * (sin(radians(latitude_p))) ** 2))) ** (3 / 2))
d_latitude = 3600 * degrees(
    (ds * cos(radians(alpha))) / M
)  # 3600 to convert to minutes and seconds
print("dφ: ", d_latitude, '"')

print(
    "position: ",
    "\n",
    "latitude: ",
    latitude_p,
    "'",
    "00 '",
    d_latitude,
    '"',
    "\n",
    "longitude: ",
    longitude_p,
    "'",
    "00 '",
    d_longitude,
    '"',
    "\n",
)

print("b)")
# east --> 10m
# ds * sin(alpha) = N * cos(latitude) * d_longitude
alpha = 90
ds = 10
d_latitude = 0

N = a / sqrt(1 - (e2 * (sin(radians(latitude_p)) ** 2)))
d_longitude = 3600 * degrees(ds * sin(radians(alpha))) / (N * cos(radians(latitude_p)))
print("dλ: ", d_longitude, '"')
latitude_a = latitude_p + d_latitude
longitude_a = longitude_p + d_longitude

print(
    "position: ",
    "\n",
    "latitude: ",
    latitude_p,
    "'",
    "00 '",
    d_latitude,
    '"',
    "\n",
    "longitude: ",
    longitude_p,
    "'",
    "00 '",
    d_longitude,
    '"',
    "\n",
)

print("c)")
# 60˚ --> 10m
# compute dφ and dλ
alpha = 60
ds = 10

d_latitude = 3600 * degrees((ds * cos(radians(alpha))) / M)
d_longitude = 3600 * degrees(ds * sin(radians(alpha))) / (N * cos(radians(latitude_p)))
print("dφ: ", d_latitude, '"')
print("dλ: ", d_longitude, '"')
latitude_a = latitude_p + d_latitude
longitude_a = longitude_p + d_longitude

print(
    "position: ",
    "\n",
    "latitude: ",
    latitude_p,
    "'",
    "00 '",
    d_latitude,
    '"',
    "\n",
    "longitude: ",
    longitude_p,
    "'",
    "00 '",
    d_longitude,
    '"',
    "\n",
)

# 3.2
print("3.2)", "\n")
# calculate the geocentric latitude (φ2) and longitude (λ2) of the station WASA
φ1 = 59 + 20 / 60  # degrees and minutes (stockholm)
λ1 = 18 + 5 / 60  # degrees and minutes (stockholm)
R = 6371000.000  # earth radius
s12 = 14723000.000  # distance to WASA from Sthlm
a12 = (
    181 + 50 / 60 + 00.00000 / 3600
)  # azimuth to WASA from Sthlm in degrees + minutes + seconds

Ψ = s12 / R
φ2 = degrees(
    asin(cos(Ψ) * sin(radians(φ1)) + sin(Ψ) * cos(radians(φ1)) * cos(radians(a12)))
)
dλ = degrees(
    acos(
        (cos(radians(φ1)) * sin(radians(φ2)) - sin(Ψ) * cos(radians(a12)))
        / (sin(radians(φ1)) * cos(radians(φ2)))
    )
)

λ2 = λ1 - dλ
print("latitude (φ2): ", φ2, "\n", "longitude (λ2): ", λ2, "\n")

# 3.3
print("3.3)", "\n")
# Gauss method of mean arguments

# given data:
φ1 = 59 + (55 / 60) + 00.00000 / 3600  # N Uppsala
λ1 = 17 + (38 / 60) + 00.00000 / 3600  # E Uppsala
φ2 = 59 + (20 / 60) + 00.00000 / 3600  # E Stockholm
λ2 = 18 + (5 / 60) + 00.00000 / 3600  # E Stockholm

s12, a12, a21 = gauss_mean_arguments(GRS80, φ1, φ2, λ1, λ2)
print("avstånd: ", s12, "m", "\n", "a12", a12, "grader", "\n", "a21: ", a21, "grader")
