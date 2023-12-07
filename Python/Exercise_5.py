from methods_AG1817 import *


print("Exercise 5.1: \n")
latitude = 59 + 21 / 60 + 0 / 3600
longitude = 18 + 4 / 60 + 9 / 3600
latitude_0 = 0
longitude_0 = 15
R = 6371000


normal_mercator_projection_for_sphere(latitude_0, longitude_0, latitude, longitude, R)

print("\nExercise 5.2: \n")

normal_mercator_projection_for_ellipsoid(
    latitude_0, longitude_0, latitude, longitude, GRS80
)

print("\nExercise 5.3: \n")
transverse_mercator_projection_for_sphere(
    latitude_0, longitude_0, latitude, longitude, R
)

print("\nExercise 5.4: \n")
transverse_mercator_projection_for_ellipsoid(
    latitude_0, longitude_0, latitude, longitude, GRS80, 0.9996
)

