from methods_AG1817 import *

print("\n", "4.1)", "\n")

simple_projection("eratosthenes")

print("\n", "4.2)", "\n")

simple_projection("lambert")

print("\n", "4.3)", "\n")
latitude_0 = 59 + 21 / 60 + 0.00000 / 3600
longitude_0 = 18 + 4 / 60 + 9 / 3600
latitude = 50 + 5 / 60 + 10 / 3600
longitude = 14 + 24 / 60 + 50 / 3600

print("a)")
oblique_orthographic_azimuthal_projection(
    latitude_0, longitude_0, latitude, longitude, 6371000.000
)
print("\n", "b)")
oblique_lambert_azimuthal_projection(
    latitude_0, longitude_0, latitude, longitude, 6371000.000
)
print("\n", "c)")
oblique_equidistant_azimuthal_projection(
    latitude_0, longitude_0, latitude, longitude, 6371000.000
)
print("\n", "d)")
oblique_stereographic_azimuthal_projection(
    latitude_0, longitude_0, latitude, longitude, 6371000.000
)
print("\n", "e)")
oblique_gnomic_azimuthal_projection(
    latitude_0, longitude_0, latitude, longitude, 6371000.000
)

print("\n", "4.4)", "\n")
latitude_1 = 45
latitude_2 = 65
R = 6371000.000

print("a)")
abers_projection(
    latitude_1, latitude_2, latitude_0, longitude_0, latitude, longitude, R
)
print("\n", "b)")
lambert_conformal_conical_projection(
    latitude_1, latitude_2, latitude_0, longitude_0, latitude, longitude, R
)
