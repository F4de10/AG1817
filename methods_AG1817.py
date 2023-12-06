from math import (
    atanh,
    cos,
    atan,
    acos,
    cosh,
    degrees,
    radians,
    asin,
    sin,
    sinh,
    sqrt,
    tan,
    pi,
)
import numpy as np
import sympy as sp
import icecream as ic

# Ellipsoid parameters: semi major axis in metres, reciprocal flattening.
GRS80 = 6378137, 298.257222101
WGS84 = 6378137, 298.257223563


def geodetic_to_geocentric(ellipsoid, latitude, longitude, height):
    """
    Return geocentric (Cartesian) Coordinates x, y, z corresponding to
    the geodetic coordinates given by latitude and longitude (in
    degrees) and height above ellipsoid. The ellipsoid must be
    specified by a pair (semi-major axis, reciprocal flattening).

    Parameters:
    - ellipsoid: tuple (semi-major axis, reciprocal flattening)
    - latitude: float, latitude in degrees
    - longitude: float, longitude in degrees
    - height: float, height above ellipsoid

    Returns:
    - x: float, geocentric x-coordinate
    - y: float, geocentric y-coordinate
    - z: float, geocentric z-coordinate
    """
    φ = radians(latitude)
    λ = radians(longitude)
    sin_φ = sin(φ)
    a, rf = ellipsoid  # semi-major axis, reciprocal flattening
    e2 = 1 - (1 - 1 / rf) ** 2  # eccentricity squared
    n = a / sqrt(1 - e2 * sin_φ**2)  # prime vertical radius
    r = (n + height) * cos(φ)  # perpendicular distance from z axis
    x = r * cos(λ)
    y = r * sin(λ)
    z = (n * (1 - e2) + height) * sin_φ
    return x, y, z


def geocentric_to_geodetic(ellipsoid, x, y, z):
    """
    Convert geocentric coordinates to geodetic coordinates.

    Args:
        ellipsoid (tuple): A pair (semi-major axis, reciprocal flattening) that specifies the ellipsoid.
        x (float): The x-coordinate in geocentric coordinates.
        y (float): The y-coordinate in geocentric coordinates.
        z (float): The z-coordinate in geocentric coordinates.

    Returns:
        tuple: A tuple containing the geodetic coordinates (latitude, longitude, height) in degrees.

    """
    a, rf = ellipsoid
    e2 = 1 - (1 - 1 / rf) ** 2
    longitude = degrees(atan(y / x))
    p = sqrt((x**2) + (y**2))
    theta = atan(z / (p * sqrt(1 - e2)))
    latitude = degrees(
        atan(
            (z + (((a * e2) / (sqrt(1 - e2))) * (sin(theta) ** 3)))
            / (p - (a * e2 * cos(theta) ** 3))
        )
    )
    N = a / (sqrt(1 - e2 * sin(radians(latitude)) ** 2))
    height = (p / cos(radians(latitude))) - N

    return (latitude, longitude, height)


def spherical_to_rectangular(r, latitude, longitude):
    """
    Converts spherical coordinates to rectangular coordinates.

    Args:
        r (float): The radial distance from the origin.
        latitude (float): The latitude angle in degrees.
        longitude (float): The longitude angle in degrees.

    Returns:
        tuple: The rectangular coordinates (x, y, z).
    """
    φ = radians(latitude)
    λ = radians(longitude)
    x = r * cos(φ) * cos(λ)
    y = r * cos(φ) * sin(λ)
    z = r * sin(φ)
    return x, y, z


def rectangular_to_spherical(x, y, z):
    """
    Converts rectangular coordinates to spherical coordinates.

    Parameters:
    x (float): The x-coordinate.
    y (float): The y-coordinate.
    z (float): The z-coordinate.

    Returns:
    tuple: A tuple containing the spherical coordinates (r, φ, λ).
        r (float): The radial distance.
        φ (float): The polar angle in degrees.
        λ (float): The azimuthal angle in degrees.
    """
    r = sqrt((x**2) + (y**2) + (z**2))
    φ = degrees(atan(z / sqrt(x**2 + y**2)))
    λ = degrees(atan(y / x))
    return r, φ, λ


def iterative_Cartesian_to_geodetic(ellipsoid, x, y, z):
    """
    Converts Cartesian coordinates (x, y, z) to geodetic coordinates (latitude, longitude, height)
    using an iterative algorithm.

    Args:
        ellipsoid (tuple): A tuple containing the semi-major axis (a) and the reciprocal of the flattening (rf) of the ellipsoid.
        x (float): The x-coordinate in Cartesian coordinates.
        y (float): The y-coordinate in Cartesian coordinates.
        z (float): The z-coordinate in Cartesian coordinates.

    Returns:
        tuple: A tuple containing the latitude, longitude, and height in geodetic coordinates.
    """
    n = 1
    height_old = 0
    latitude_old = 0
    dh = 1
    dφ = 1

    a, rf = ellipsoid
    e2 = 1 - (1 - 1 / rf) ** 2
    p = sqrt((x**2) + (y**2))
    longitude_old = degrees(atan(y / x))

    i = 0
    while abs(dh) > 0.001 or abs(dφ) > 1e-09:
        latitude_new = degrees(atan(z / (p * (1 - (e2 * (n / (n + height_old)))))))
        n = a / sqrt(1 - e2 * sin(radians(latitude_new)) ** 2)
        height_new = (p / cos(radians(latitude_new))) - n
        dh = height_new - height_old
        dφ = latitude_new - latitude_old
        latitude_old = latitude_new
        height_old = height_new
        i = i + 1
        print("iteration: ", i)
        print("dh: ", dh, "\n", "dφ: ", dφ)
        print(
            "latitude: ",
            latitude_new,
            "\n",
            "longitude: ",
            longitude_old,
            "\n",
            "height: ",
            height_new,
        )
    return latitude_old, longitude_old, height_old


def M_N_W(ellipsoid, latitude):
    """
    Calculate the M, N, and W values for a given ellipsoid and latitude.

    Parameters:
    ellipsoid (tuple): A tuple containing the semi-major axis and reciprocal flattening of the ellipsoid.
    latitude (float): The latitude in degrees.

    Returns:
    tuple: A tuple containing the M, N, and W values.
    """
    a, rf = ellipsoid  # semi-major axis, reciprocal flattening
    e2 = 1 - (1 - 1 / rf) ** 2  # eccentricity squared
    φm = radians(latitude)
    W = sqrt((1 - e2 * sin(φm) ** 2))
    N = a / W
    M = (a * (1 - e2)) / W**3
    return M, N, W


def gauss_mean_arguments(ellipsoid, φ1, φ2, λ1, λ2):
    """
    Calculate the distance and azimuth between two points on an ellipsoid using the Gauss mean method.

    Args:
        ellipsoid (tuple): A tuple containing the semi-major axis and reciprocal flattening of the ellipsoid.
        φ1 (float): Latitude of the first point in degrees.
        φ2 (float): Latitude of the second point in degrees.
        λ1 (float): Longitude of the first point in degrees.
        λ2 (float): Longitude of the second point in degrees.

    Returns:
        tuple: A tuple containing the distance between the points in meters, the forward azimuth from the first point to the second point in degrees, and the backward azimuth from the second point to the first point in degrees.
    """
    a, rf = ellipsoid  # semi-major axis, reciprocal flattening
    e2 = 1 - (1 - 1 / rf) ** 2  # eccentricity squared
    e_p = sqrt(e2 / (1 - e2))
    φ_m = (radians(φ1) + radians(φ2)) / 2
    λ_m = (radians(λ1) + radians(λ2)) / 2

    W = sqrt((1 - (e2 * (sin(φ_m) ** 2))))
    N = a / W
    M = (a * (1 - e2)) / (W**3)

    n_m2 = (e_p * cos(φ_m)) ** 2
    t_m2 = tan(φ_m) ** 2
    t_m4 = t_m2**2
    v_m2 = 1 + n_m2
    v_m4 = v_m2**2

    c1 = 1 / M
    c2 = 1 / N
    c3 = 1 / 24
    c4 = (1 + n_m2 - 9 * n_m2 * t_m2) / (24 * (v_m4))
    c5 = (1 - 2 * n_m2) / 24
    c6 = (n_m2 * (t_m2 - 1 - n_m2 - 4 * n_m2 * t_m2)) / (8 * (v_m4))
    c7 = v_m2 / 12
    c8 = (3 + 5 * n_m2) / (24 * v_m2)
    c9 = 1 / 2880
    c10 = ((4 + 15 * t_m2) * (cos(φ_m) ** 2)) / 1440
    c11 = ((12 * t_m2 + t_m4) * (cos(φ_m) ** 4)) / 2880
    c12 = ((14 + 40 * t_m2 + 15 * t_m4) * (cos(φ_m) ** 4)) / 2880
    c13 = 1 / 192
    c14 = (sin(φ_m) ** 2) / 48
    c15 = ((7 - 6 * t_m2) * (cos(φ_m) ** 4)) / 1440

    dφ = radians(φ2) - radians(φ1)
    dφ_2 = dφ**2
    dλ = radians(λ2) - radians(λ1)
    dλ_2 = dλ**2

    delta_1 = ((dφ * cos(0.5 * dλ)) / c1) * (
        1
        + (c5 * (cos(φ_m) ** 2) * dλ_2)
        - (c6 * dφ_2)
        - (((c10 * dφ_2 * dλ_2) + (c12 * (dλ_2**2))))
    )
    delta_2 = (dλ * cos(φ_m) / c2) * (
        1
        - c3 * sin(φ_m) ** 2 * dλ_2
        + c4 * dφ_2
        + (c9 * dφ_2**2 - c10 * dφ_2 * dλ_2 - c11 * dλ_2**2)
    )
    d_alpha = degrees(
        sin(φ_m)
        * dλ
        * (
            1
            + c7 * cos(φ_m) ** 2 * dλ_2
            + c8 * dφ_2
            + (c13 * dφ_2**2 - c14 * dφ_2 * dλ_2 + c15 * dλ_2**2)
        )
    )

    s12 = sqrt((delta_1**2) + (delta_2**2))
    a_m = degrees(atan(delta_2 / delta_1)) + 180
    a12 = a_m - 0.5 * d_alpha
    a21 = a_m + 0.5 * d_alpha + 180

    return s12, a12, a21


def simple_projection(projection):
    """
    Perform calculations for different types of map projections.

    Args:
        projection (str): The type of projection to use. Can be "eratosthenes" or "lambert".

    Returns:
        None
    """

    R = sp.symbols("R", positive=True, real=True)
    latitude, longitude = sp.symbols("latitude longitude", real=True)

    eratosthenes = R * latitude, R * longitude
    lambert = R * sp.sin(latitude), R * longitude

    if projection == "eratosthenes":
        x, y = eratosthenes
    elif projection == "lambert":
        x, y = lambert

    # first fundamental coefficients (FFC)
    print("\na)")
    e = sp.simplify(sp.diff(x, latitude) ** 2 + sp.diff(y, latitude) ** 2)  # type: ignore
    print("\ne: ", e)

    f = sp.simplify(
        sp.diff(x, latitude) * sp.diff(x, longitude)  # type: ignore
        + sp.diff(y, latitude) * sp.diff(y, longitude)  # type: ignore
    )
    print("f: ", f)

    g = sp.simplify(sp.diff(x, longitude) ** 2 + sp.diff(y, longitude) ** 2)  # type: ignore
    print("g: ", g)

    H = sp.simplify(sp.sqrt((e * g) - (f**2)))
    H_no_abs = H.subs(
        [
            (sp.Abs(sp.cos(latitude)), sp.cos(latitude)),
            (sp.Abs(sp.sin(latitude)), sp.sin(latitude)),
        ]
    )
    print("H: ", H_no_abs)

    # deformation parameters ƞ, k, 𝛆
    print("\nb)")
    n = sp.sqrt(e) / R
    n_no_abs = n.subs(
        [
            (sp.Abs(sp.cos(latitude)), sp.cos(latitude)),
            (sp.Abs(sp.sin(latitude)), sp.sin(latitude)),
        ]
    )
    print("\nh = ", n_no_abs)

    k = sp.simplify(sp.sqrt(g) / (R * sp.cos(latitude)))
    print("k = ", k)

    epsilon = sp.simplify(H / (R**2 * sp.cos(latitude)))
    epsilon_no_abs = epsilon.subs(
        [
            (sp.Abs(sp.cos(latitude)), sp.cos(latitude)),
            (sp.Abs(sp.sin(latitude)), sp.sin(latitude)),
        ]
    )
    print("𝜺 =", epsilon_no_abs)

    # characteristics
    print("\nc)")
    if n_no_abs == k:
        print("\nh = k, conformal projection")
    else:
        print("\nh != k, non-conformal projection")

    if n_no_abs * k == 1:
        print("h * k = 1, equal-area projection")
    else:
        print("h * k != 1, non-equal-area projection")

    if n_no_abs == 1:
        print("h = 1, equidistant at meridian")
    elif k == 1:
        print("k = 1, equidistant for parallel circles")
    elif sp.simplify(epsilon_no_abs, force=True) == 1:
        print("μ = 1, equidistant for some curves (e.g., great circles)")

    # are the projections of the meridians and parallel circles vertical to each other?
    print("\nd)")
    theta_prim = degrees(acos(f / sp.sqrt(e * g)))
    print("\nthe angle between the projections of the meridians and parallel circles:")
    print("θ' = ", theta_prim)
    if theta_prim == 90:
        print("--> The meridians and parallel circles are vertical")

    return


def oblique_orthographic_azimuthal_projection(
    latitude_0, longitude_0, latitude, longitude, R
):
    """
    Performs an oblique orthographic azimuthal projection.

    Args:
        latitude_0 (float): The latitude of the center of projection.
        longitude_0 (float): The longitude of the center of projection.
        latitude (float): The latitude of the point to project.
        longitude (float): The longitude of the point to project.
        R (float): The radius of the sphere.

    Returns:
        None
    """
    x_prim = R * (
        cos(radians(latitude_0)) * sin(radians(latitude))
        - sin(radians(latitude_0))
        * cos(radians(latitude))
        * cos(radians(longitude) - radians(longitude_0))
    )
    y_prim = R * cos(radians(latitude)) * sin(radians(longitude) - radians(longitude_0))
    print("x': ", round(x_prim, 9), "\n", "y': ", round(y_prim, 9))

    r = sqrt(x_prim**2 + y_prim**2)
    print("r: ", round(r, 9))

    h = sin(radians(latitude_0)) * sin(radians(latitude)) + cos(
        radians(latitude_0)
    ) * cos(radians(latitude)) * cos(radians(longitude) - radians(longitude_0))
    k = 1
    𝜉 = sin(radians(latitude_0)) * sin(radians(latitude)) + cos(
        radians(latitude_0) * cos(radians(latitude))
    ) * cos(radians(longitude) - radians(longitude_0))
    # print("\n", "h: ", h, "\n", "k: ", k, "\n", "𝜉: ", 𝜉)
    return


def oblique_lambert_azimuthal_projection(
    latitude_0, longitude_0, latitude, longitude, R
):
    """
    Calculates the oblique Lambert azimuthal projection coordinates (x', y') and the radial distance (r)
    for a given latitude and longitude, using the oblique Lambert azimuthal projection formula.

    Parameters:
    latitude_0 (float): The latitude of the projection center in degrees.
    longitude_0 (float): The longitude of the projection center in degrees.
    latitude (float): The latitude of the point to project in degrees.
    longitude (float): The longitude of the point to project in degrees.
    R (float): The radius of the sphere or ellipsoid in the same units as the coordinates.

    Returns:
    None: The function prints the calculated coordinates (x', y') and the radial distance (r).
    """
    k = sqrt(
        2
        / (
            1
            + sin(radians(latitude_0)) * sin(radians(latitude))
            + cos(radians(latitude_0))
            * cos(radians(latitude))
            * cos(radians(longitude) - radians(longitude_0))
        )
    )
    y_prim = (
        R * k * cos(radians(latitude)) * sin(radians(longitude) - radians(longitude_0))
    )
    x_prim = (
        R
        * k
        * (
            cos(radians(latitude_0)) * sin(radians(latitude))
            - sin(radians(latitude_0))
            * cos(radians(latitude))
            * cos(radians(longitude) - radians(longitude_0))
        )
    )
    r = sqrt(x_prim**2 + y_prim**2)
    print(
        "x': ",
        round(x_prim, 9),
        "\n",
        "y': ",
        round(y_prim, 9),
        "\n",
        "r: ",
        round(r, 9),
    )


def oblique_equidistant_azimuthal_projection(
    latitude_0, longitude_0, latitude, longitude, R
):
    """
    Calculates the oblique equidistant azimuthal projection coordinates (x', y') and the radial distance (r)
    from the given latitude and longitude to the reference latitude and longitude.

    Parameters:
    - latitude_0 (float): The reference latitude in degrees.
    - longitude_0 (float): The reference longitude in degrees.
    - latitude (float): The latitude in degrees.
    - longitude (float): The longitude in degrees.
    - R (float): The radius of the sphere.

    Returns:
    None
    """
    c = acos(
        sin(radians(latitude_0)) * sin(radians(latitude))
        + cos(radians(latitude_0))
        * cos(radians(latitude))
        * cos(radians(longitude) - radians(longitude_0))
    )
    sin_lambda_prim = (
        cos(radians(latitude)) * sin(radians(longitude - longitude_0))
    ) / sin(c)
    cos_lambda_prim = (
        -cos(radians(latitude_0)) * sin(radians(latitude))
        + sin(radians(latitude_0))
        * cos(radians(latitude))
        * cos(radians(longitude - longitude_0))
    ) / sin(c)
    x_prim = -R * c * cos_lambda_prim
    y_prim = R * c * sin_lambda_prim
    r = sqrt(x_prim**2 + y_prim**2)
    print(
        "x': ",
        round(x_prim, 9),
        "\n",
        "y': ",
        round(y_prim, 9),
        "\n",
        "r: ",
        round(r, 9),
    )


def oblique_stereographic_azimuthal_projection(
    latitude_0, longitude_0, latitude, longitude, R
):
    """
    Performs an oblique stereographic azimuthal projection.

    Args:
        latitude_0 (float): Latitude of the projection center in degrees.
        longitude_0 (float): Longitude of the projection center in degrees.
        latitude (float): Latitude of the point to project in degrees.
        longitude (float): Longitude of the point to project in degrees.
        R (float): Radius of the sphere in the projection.

    Returns:
        None: The function prints the projected coordinates and the distance from the projection center.

    """
    k = 2 / (
        1
        + sin(radians(latitude_0)) * sin(radians(latitude))
        + cos(radians(latitude_0))
        * cos(radians(latitude))
        * cos(radians(longitude) - radians(longitude_0))
    )
    x_prim = (
        R
        * k
        * (
            cos(radians(latitude_0)) * sin(radians(latitude))
            - sin(radians(latitude_0))
            * cos(radians(latitude))
            * cos(radians(longitude) - radians(longitude_0))
        )
    )
    y_prim = (
        R * k * cos(radians(latitude)) * sin(radians(longitude) - radians(longitude_0))
    )
    print("x': ", round(x_prim, 4), "\n", "y': ", round(y_prim, 4))

    r = sqrt(x_prim**2 + y_prim**2)
    print("r: ", round(r, 9))


def oblique_gnomic_azimuthal_projection(
    latitude_0, longitude_0, latitude, longitude, R
):
    """
    Calculates the oblique gnomic azimuthal projection coordinates (x', y') and the radial distance (r)
    from the given latitude and longitude to the reference latitude and longitude.

    Parameters:
    - latitude_0 (float): The reference latitude in degrees.
    - longitude_0 (float): The reference longitude in degrees.
    - latitude (float): The latitude in degrees.
    - longitude (float): The longitude in degrees.
    - R (float): The radius of the sphere.

    Returns:
    None
    """
    k = 1 / (
        sin(radians(latitude_0)) * sin(radians(latitude))
        + cos(radians(latitude_0))
        * cos(radians(latitude))
        * cos(radians(longitude) - radians(longitude_0))
    )
    x_prim = (
        R
        * k
        * (
            cos(radians(latitude_0)) * sin(radians(latitude))
            - sin(radians(latitude_0))
            * cos(radians(latitude))
            * cos(radians(longitude) - radians(longitude_0))
        )
    )
    y_prim = (
        R * k * cos(radians(latitude)) * sin(radians(longitude) - radians(longitude_0))
    )
    print("x': ", round(x_prim, 4), "\n", "y': ", round(y_prim, 4))

    r = sqrt(x_prim**2 + y_prim**2)
    print("r: ", round(r, 9))


def Abers_projection(
    latitude_1, latitude_2, latitude_0, longitude_0, latitude, longitude, R
):
    """
    Perform Abers projection to convert geographic coordinates to Cartesian coordinates.

    Parameters:
    latitude_1 (float): Latitude of the first standard parallel in radians.
    latitude_2 (float): Latitude of the second standard parallel in radians.
    latitude_0 (float): Latitude of the projection origin in radians.
    longitude_0 (float): Longitude of the projection origin in radians.
    latitude (float): Latitude of the point to be projected in radians.
    longitude (float): Longitude of the point to be projected in radians.
    R (float): Radius of the Earth.

    Returns:
    None
    """
    latitude_1, latitude_2, latitude, longitude, latitude_0, longitude_0 = (
        radians(latitude_1),
        radians(latitude_2),
        radians(latitude),
        radians(longitude),
        radians(latitude_0),
        radians(longitude_0),
    )
    n = (sin(latitude_1) + sin(latitude_2)) / 2
    C = cos(latitude_1) ** 2 + (2 * n * sin(latitude_1))
    rho = (R / n) * sqrt(C - (2 * n * sin(latitude)))
    rho_0 = (R / n) * sqrt(C - (2 * n * sin(latitude_0)))
    x = rho_0 - (rho * cos(n * (longitude - longitude_0)))
    y = rho * sin(n * (longitude - longitude_0))
    r = sqrt(x**2 + y**2)
    print("x: ", round(x, 4), "\n", "y: ", round(y, 4), "\n", "r: ", round(r, 4))


def lambert_conformal_conical_projection(
    latitude_1, latitude_2, latitude_0, longitude_0, latitude, longitude, R
):
    """
    Performs the Lambert Conformal Conical Projection on a given latitude and longitude.

    Args:
        latitude_1 (float): First standard parallel latitude in radians.
        latitude_2 (float): Second standard parallel latitude in radians.
        latitude_0 (float): Latitude of the projection origin in radians.
        longitude_0 (float): Longitude of the projection origin in radians.
        latitude (float): Latitude to project in radians.
        longitude (float): Longitude to project in radians.
        R (float): Radius of the Earth.

    Returns:
        None
    """
    latitude_1, latitude_2, latitude, longitude, latitude_0, longitude_0 = (
        radians(latitude_1),
        radians(latitude_2),
        radians(latitude),
        radians(longitude),
        radians(latitude_0),
        radians(longitude_0),
    )
    n = (np.log(cos(latitude_1) * (1 / cos(latitude_2)))) / (
        np.log(
            tan((pi / 4) + (latitude_2 / 2)) * (1 / tan((pi / 4) + (latitude_1 / 2)))
        )
    )
    F = (cos(latitude_1) / n) * (tan((pi / 4) + (latitude_1 / 2)) ** n)
    rho = (R * F) / (tan((pi / 4) + (latitude / 2)) ** n)
    rho_0 = (R * F) / (tan((pi / 4) + (latitude_0 / 2)) ** n)
    x = rho_0 - (rho * cos(n * (longitude - longitude_0)))
    y = rho * sin(n * (longitude - longitude_0))
    r = sqrt(x**2 + y**2)
    print("x: ", round(x, 4), "\n", "y: ", round(y, 4), "\n", "r: ", round(r, 4))
    return x, y, r


def normal_mercator_projection_for_sphere(
    latitude_0, longitude_0, latitude, longitude, R
):
    """
    Calculates the normal Mercator projection coordinates (x', y') and the radial distance (r)
    from the given latitude and longitude to the reference latitude and longitude.

    Parameters:
    - latitude_0 (float): The reference latitude in degrees.
    - longitude_0 (float): The reference longitude in degrees.
    - latitude (float): The latitude in degrees.
    - longitude (float): The longitude in degrees.
    - R (float): The radius of the sphere.

    Returns:
    None
    """

    x = R * cos(radians(latitude_0)) * np.log(tan((pi / 4) + (radians(latitude) / 2)))
    y = R * cos(radians(latitude_0)) * radians(longitude - longitude_0)

    print("x: ", round(x, 4), "\n", "y: ", round(y, 4))
    return x, y


def normal_mercator_projection_for_ellipsoid(
    latitude_0, longitude_0, latitude, longitude, ellipsoid
):
    """
    Calculates the normal Mercator projection for an ellipsoid.

    Args:
        latitude_0 (float): The reference latitude in degrees.
        longitude_0 (float): The reference longitude in degrees.
        latitude (float): The latitude in degrees.
        longitude (float): The longitude in degrees.
        ellipsoid (tuple): A tuple containing the semi-major axis (a) and the reciprocal of the flattening (rf) of the ellipsoid.

    Returns:
        tuple: A tuple containing the projected x-coordinate, projected y-coordinate, scale factor (h), convergence factor (k), and eccentricity squared (epsilon).

    """

    latitude_0, longitude_0, latitude, longitude = (
        radians(latitude_0),
        radians(longitude_0),
        radians(latitude),
        radians(longitude),
    )

    # Input parameters
    a, rf = ellipsoid
    e2 = 1 - (1 - 1 / rf) ** 2
    e = sqrt(e2)
    C = (a * cos(latitude_0)) / (sqrt(1 - (e2 * (sin(latitude_0) ** 2))))
    x = C * np.log(
        tan(pi / 4 + latitude / 2)
        * ((1 - e * sin(latitude)) / (1 + e * sin(latitude))) ** (e / 2)
    )
    print("x: ", round(x, 4))
    y = C * (longitude - longitude_0)
    print("y: ", round(y, 4))

    h = k = (
        cos(latitude_0)
        / cos(latitude)
        * sqrt((1 - e2 * sin(latitude) ** 2) / (1 - e2 * sin(latitude_0) ** 2))
    )
    epsilon = k**2

    return x, y, h, k, epsilon


def transverse_mercator_projection_for_sphere(
    latitude_0, longitude_0, latitude, longitude, R
):
    """
    Transverse Mercator Projection for Sphere.

    Calculates the x, y, r, h, and k coordinates using the Transverse Mercator Projection formula for a sphere.

    Parameters:
    - latitude_0 (float): Reference latitude in radians.
    - longitude_0 (float): Reference longitude in radians.
    - latitude (float): Latitude in radians.
    - longitude (float): Longitude in radians.
    - R (float): Radius of the sphere.

    Returns:
    - x (float): x-coordinate.
    - y (float): y-coordinate.
    - r (float): Distance from the origin.
    - h (float): Scale factor in the x-direction.
    - k (float): Scale factor in the y-direction.
    """

    latitude_0, longitude_0, latitude, longitude = (
        radians(latitude_0),
        radians(longitude_0),
        radians(latitude),
        radians(longitude),
    )

    B = cos(latitude) * sin(longitude - longitude_0)
    x = R * (atan(tan(latitude) / cos(longitude - longitude_0)) - latitude_0)
    y = 0.5 * R * np.log((1 + B) / (1 - B))
    h = k = 1 / sqrt(1 - (B**2))
    print("x: ", round(x, 4), "\n", "y: ", round(y, 4))

    return x, y, h, k


def transverse_mercator_projection_for_ellipsoid(
    latitude_0, longitude_0, latitude, longitude, ellipsoid, scale_factor
):
    """
    Transforms geographic coordinates (latitude, longitude) to projected coordinates (x, y)
    using the Transverse Mercator projection for a given ellipsoid.

    Args:
        latitude_0 (float): Latitude of the projection origin in degrees.
        longitude_0 (float): Longitude of the projection origin in degrees.
        latitude (float): Latitude of the point to be projected in degrees.
        longitude (float): Longitude of the point to be projected in degrees.
        ellipsoid (tuple): Tuple containing the semi-major axis (a) and reciprocal flattening (rf) of the ellipsoid.
        scale_factor (float): Scale factor of the projection.

    Returns:
        tuple: A tuple containing the projected coordinates (x, y).
    """

    a, rf = ellipsoid
    e2 = 1 - (1 - 1 / rf) ** 2
    b = a * (1 - 1 / rf)
    latitude_0, latitude = (radians(latitude_0), radians(latitude))

    n = (a - b) / (a + b)
    n2 = n**2
    n3 = n2 * n
    n4 = n3 * n

    a_hat = (a / (1 + n)) * (1 + (n2 / 4) + (n4 / 64))

    beta_1 = n / 2 - (2 * n2 / 3) + (5 * n3 / 16) + (41 * n4 / 180)
    beta_2 = (13 * n2 / 48) - (3 * n3 / 5) + (557 * n4 / 1440)
    beta_3 = (61 * n3 / 240) - (103 * n4 / 140)
    beta_4 = 49561 * n4 / 161280

    A = e2
    B = ((5 * e2**2) - (e2**3)) / 6
    C = ((104 * e2**3) - (45 * e2**4)) / 120
    D = (1237 * e2**4) / 1260

    latitude_temp = latitude - (
        sin(latitude)
        * cos(latitude)
        * (A + B * sin(latitude) ** 2 + C * sin(latitude) ** 4 + D * sin(latitude) ** 6)
    )

    delta_longitude = radians(longitude - longitude_0)

    epsilon_prim = atan(tan(latitude_temp) / cos(delta_longitude))

    eta_prim = atanh(cos(latitude_temp) * sin(delta_longitude))

    x = (
        a_hat
        * scale_factor
        * (
            epsilon_prim
            + beta_1 * sin(2 * epsilon_prim) * cosh(2 * eta_prim)
            + beta_2 * sin(4 * epsilon_prim) * cosh(4 * eta_prim)
            + beta_3 * sin(6 * epsilon_prim) * cosh(6 * eta_prim)
            + beta_4 * sin(8 * epsilon_prim) * cosh(8 * eta_prim)
        )
    )
    y = (
        a_hat
        * scale_factor
        * (
            eta_prim
            + beta_1 * cos(2 * epsilon_prim) * sinh(2 * eta_prim)
            + beta_2 * cos(4 * epsilon_prim) * sinh(4 * eta_prim)
            + beta_3 * cos(6 * epsilon_prim) * sinh(6 * eta_prim)
            + beta_4 * cos(8 * epsilon_prim) * sinh(8 * eta_prim)
        )
        + 500000  # 500 km false easting
    )

    print("x: ", round(x, 4), "\n", "y: ", round(y, 4))
    return x, y


def convert_to_hours_minutes_seconds(time):
    """
    Converts a time in seconds to hours, minutes, and seconds.

    Args:
        time (float): The time in seconds.

    Returns:
        tuple: A tuple containing the hours, minutes, and seconds.
    """
    hours = int(time / 3600)
    minutes = int((time - (hours * 3600)) / 60)
    seconds = time - (hours * 3600) - (minutes * 60)
    return hours, minutes, seconds


def convert_to_degrees_minutes_seconds(degrees):
    """
    Converts degrees to degrees, minutes, and seconds.

    Args:
        degrees (float): The input degrees.

    Returns:
        tuple: A tuple containing the degrees, minutes, and seconds.
    """
    degrees = abs(degrees)
    deg_int = int(degrees)
    minutes = int((degrees - deg_int) * 60)
    seconds = (degrees - deg_int - (minutes / 60)) * 3600
    return deg_int, minutes, seconds


if __name__ == "__main__":
    print(convert_to_degrees_minutes_seconds(60.59444444444444))