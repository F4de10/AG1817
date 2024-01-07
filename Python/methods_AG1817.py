import numpy as np
import sympy as sp

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
    latitude = np.radians(latitude)
    longitude = np.radians(longitude)
    a, rf = ellipsoid  # semi-major axis, reciprocal flattening
    e2 = 1 - (1 - 1 / rf) ** 2  # eccentricity squared
    n = a / np.sqrt(1 - e2 * np.sin(latitude) ** 2)  # prime vertical radius
    r = (n + height) * np.cos(latitude)  # perpendicular distance from z axis
    x = r * np.cos(longitude)
    y = r * np.sin(longitude)
    z = (n * (1 - e2) + height) * np.sin(latitude)
    return x, y, z


def geocentric_to_geodetic(ellipsoid, x, y, z, height_0=None):
    """
    Convert geocentric coordinates to geodetic coordinates.

    Args:
        ellipsoid (tuple): A pair (semi-major axis, reciprocal flattening) that specifies the ellipsoid.
        x (float): The x-coordinate in geocentric coordinates.
        y (float): The y-coordinate in geocentric coordinates.
        z (float): The z-coordinate in geocentric coordinates.

    Returns:
        tuple: A tuple containing the geodetic coordinates (latitude, longitude, height) in np.degrees.

    """
    a, rf = ellipsoid
    e2 = 1 - (1 - 1 / rf) ** 2
    longitude = np.degrees(np.arctan(y / x))
    p = np.sqrt((x**2) + (y**2))
    theta = np.arctan(z / (p * np.sqrt(1 - e2)))

    if height_0 == 0:
        latitude = np.degrees(
            np.arctan((1 / (1 - e2)) * (z / (np.sqrt((x**2) + (y**2)))))
        )
    else:
        latitude = np.degrees(
            np.arctan(
                (z + (((a * e2) / (np.sqrt(1 - e2))) * (np.sin(theta) ** 3)))
                / (p - (a * e2 * np.cos(theta) ** 3))
            )
        )

    N = a / (np.sqrt(1 - e2 * np.sin(np.radians(latitude)) ** 2))
    height = (p / np.cos(np.radians(latitude))) - N

    return latitude, longitude, height


def spherical_to_rectangular(r, latitude, longitude):
    """
    Converts spherical coordinates to rectangular coordinates.

    Args:
        r (float): The radial distance from the origin.
        latitude (float): The latitude angle in np.degrees.
        longitude (float): The longitude angle in np.degrees.

    Returns:
        tuple: The rectangular coordinates (x, y, z).
    """
    latitude = np.radians(latitude)
    longitude = np.radians(longitude)
    x = r * np.cos(latitude) * np.cos(longitude)
    y = r * np.cos(latitude) * np.sin(longitude)
    z = r * np.sin(latitude)
    return x, y, z


def rectangular_to_spherical(x, y, z):
    """
    Converts rectangular coordinates to spherical coordinates.

    Parameters:
    x (float): The x-coordinate.
    y (float): The y-coordinate.
    z (float): The z-coordinate.

    Returns:
    tuple: A tuple containing the spherical coordinates (r, latitude, longitude).
        r (float): The radial distance.
        latitude (float): The polar angle in np.degrees.
        longitude (float): The azimuthal angle in np.degrees.
    """
    r = np.sqrt((x**2) + (y**2) + (z**2))
    latitude = np.degrees(np.arctan(z / np.sqrt(x**2 + y**2)))
    longitude = np.degrees(np.arctan(y / x))
    return r, latitude, longitude


def iterative_cartesian_to_geodetic(ellipsoid, x, y, z):
    """
    Converts Cartesian coordinates (x, y, z) to geodetic coordinates (latitude, longitude, height)
    using an iterative algorithm.

    Args:
        ellipsoid (tuple): A tuple containing the semi-major axis (a)
            and the reciprocal of the flattening (rf) of the ellipsoid.
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
    d_latitude = 1

    a, rf = ellipsoid
    e2 = 1 - (1 - 1 / rf) ** 2
    p = np.sqrt((x**2) + (y**2))
    longitude_old = np.degrees(np.arctan(y / x))

    i = 0
    while abs(dh) > 0.001 or abs(d_latitude) > 1e-09:
        latitude_new = np.degrees(
            np.arctan(z / (p * (1 - (e2 * (n / (n + height_old))))))
        )
        n = a / np.sqrt(1 - e2 * np.sin(np.radians(latitude_new)) ** 2)
        height_new = (p / np.cos(np.radians(latitude_new))) - n
        dh = height_new - height_old
        d_latitude = latitude_new - latitude_old
        latitude_old = latitude_new
        height_old = height_new
        i = i + 1
        print("iteration: ", i)
        print("dh: ", dh, "\n", "d_latitude: ", d_latitude)
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


def m_n_w(ellipsoid, latitude):
    """
    Calculate the M, N, and W values for a given ellipsoid and latitude.

    Parameters:
    ellipsoid (tuple): A tuple containing the semi-major axis and reciprocal flattening of the ellipsoid.
    latitude (float): The latitude in np.degrees.

    Returns:
    tuple: A tuple containing the M, N, and W values.
    """
    a, rf = ellipsoid  # semi-major axis, reciprocal flattening
    e2 = 1 - (1 - 1 / rf) ** 2  # eccentricity squared
    latitude_m = np.radians(latitude)
    W = np.sqrt((1 - e2 * np.sin(latitude_m) ** 2))
    N = a / W
    M = (a * (1 - e2)) / W**3
    return M, N, W


def gauss_mean_arguments(ellipsoid, latitude_1, latitude_2, longitude_1, longitude_2):
    """
    Calculate the distance and azimuth between two points on an ellipsoid using the Gauss mean method.

    Args:
        ellipsoid (tuple): A tuple containing the semi-major axis and reciprocal flattening of the ellipsoid.
        latitude_1 (float): Latitude of the first point in np.degrees.
        latitude_2 (float): Latitude of the second point in np.degrees.
        longitude_1 (float): Longitude of the first point in np.degrees.
        longitude_2 (float): Longitude of the second point in np.degrees.

    Returns:
        tuple: A tuple containing the distance between the points in meters,
            the forward azimuth from the first point to the second point in degrees,
            and the backward azimuth from the second point to the first point in np.degrees.
    """
    a, rf = ellipsoid  # semi-major axis, reciprocal flattening
    e2 = 1 - (1 - 1 / rf) ** 2  # eccentricity squared
    e_p = np.sqrt(e2 / (1 - e2))
    latitude_m = (np.radians(latitude_1) + np.radians(latitude_2)) / 2
    # longitude_m = (np.radians(longitude_1) + np.radians(longitude_2)) / 2

    W = np.sqrt((1 - (e2 * (np.sin(latitude_m) ** 2))))
    N = a / W
    M = (a * (1 - e2)) / (W**3)

    n_m2 = (e_p * np.cos(latitude_m)) ** 2
    t_m2 = np.tan(latitude_m) ** 2
    t_m4 = t_m2**2
    v_m2 = 1 + n_m2
    v_m4 = v_m2**2

    c1 = 1 / M
    c2 = 1 / N
    c3 = 1 / 24
    c4 = (1 + n_m2 - 9 * n_m2 * t_m2) / (24 * v_m4)
    c5 = (1 - 2 * n_m2) / 24
    c6 = (n_m2 * (t_m2 - 1 - n_m2 - 4 * n_m2 * t_m2)) / (8 * v_m4)
    c7 = v_m2 / 12
    c8 = (3 + 5 * n_m2) / (24 * v_m2)
    c9 = 1 / 2880
    c10 = ((4 + 15 * t_m2) * (np.cos(latitude_m) ** 2)) / 1440
    c11 = ((12 * t_m2 + t_m4) * (np.cos(latitude_m) ** 4)) / 2880
    c12 = ((14 + 40 * t_m2 + 15 * t_m4) * (np.cos(latitude_m) ** 4)) / 2880
    c13 = 1 / 192
    c14 = (np.sin(latitude_m) ** 2) / 48
    c15 = ((7 - 6 * t_m2) * (np.cos(latitude_m) ** 4)) / 1440

    d_latitude = np.radians(latitude_2) - np.radians(latitude_1)
    d_latitude_2 = d_latitude**2
    d_longitude = np.radians(longitude_2) - np.radians(longitude_1)
    d_longitude_2 = d_longitude**2

    delta_1 = ((d_latitude * np.cos(0.5 * d_longitude)) / c1) * (
        1
        + (c5 * (np.cos(latitude_m) ** 2) * d_longitude_2)
        - (c6 * d_latitude_2)
        - ((c10 * d_latitude_2 * d_longitude_2) + (c12 * (d_longitude_2**2)))
    )
    delta_2 = (d_longitude * np.cos(latitude_m) / c2) * (
        1
        - c3 * np.sin(latitude_m) ** 2 * d_longitude_2
        + c4 * d_latitude_2
        + (
            c9 * d_latitude_2**2
            - c10 * d_latitude_2 * d_longitude_2
            - c11 * d_longitude_2**2
        )
    )
    d_alpha = np.degrees(
        np.sin(latitude_m)
        * d_longitude
        * (
            1
            + c7 * np.cos(latitude_m) ** 2 * d_longitude_2
            + c8 * d_latitude_2
            + (
                c13 * d_latitude_2**2
                - c14 * d_latitude_2 * d_longitude_2
                + c15 * d_longitude_2**2
            )
        )
    )

    s12 = np.sqrt((delta_1**2) + (delta_2**2))
    a_m = np.degrees(np.arctan(delta_2 / delta_1)) + 180
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

    x, y = None, None
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
    theta_prim = sp.acos(f / sp.sqrt(e * g))
    print("\nthe angle between the projections of the meridians and parallel circles:")
    print("θ' = ", theta_prim)
    if theta_prim == 90:
        print("--> The meridians and parallel circles are vertical")

    return


# noinspection DuplicatedCode
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
        np.cos(np.radians(latitude_0)) * np.sin(np.radians(latitude))
        - np.sin(np.radians(latitude_0))
        * np.cos(np.radians(latitude))
        * np.cos(np.radians(longitude) - np.radians(longitude_0))
    )
    y_prim = (
        R
        * np.cos(np.radians(latitude))
        * np.sin(np.radians(longitude) - np.radians(longitude_0))
    )
    print("x': ", round(x_prim, 9), "\n", "y': ", round(y_prim, 9))

    r = np.sqrt(x_prim**2 + y_prim**2)
    print("r: ", round(r, 9))

    h = np.sin(np.radians(latitude_0)) * np.sin(np.radians(latitude)) + np.cos(
        np.radians(latitude_0)
    ) * np.cos(np.radians(latitude)) * np.cos(
        np.radians(longitude) - np.radians(longitude_0)
    )
    k = 1
    xi = np.sin(np.radians(latitude_0)) * np.sin(np.radians(latitude)) + np.cos(
        np.radians(latitude_0) * np.cos(np.radians(latitude))
    ) * np.cos(np.radians(longitude) - np.radians(longitude_0))
    return x_prim, y_prim, r, h, k, xi


def oblique_lambert_azimuthal_projection(
    latitude_0, longitude_0, latitude, longitude, R
):
    """
    Calculates the oblique Lambert azimuthal projection coordinates (x', y') and the radial distance (r)
    for a given latitude and longitude, using the oblique Lambert azimuthal projection formula.

    Parameters:
    latitude_0 (float): The latitude of the projection center in np.degrees.
    longitude_0 (float): The longitude of the projection center in np.degrees.
    latitude (float): The latitude of the point to project in np.degrees.
    longitude (float): The longitude of the point to project in np.degrees.
    R (float): The radius of the sphere or ellipsoid in the same units as the coordinates.

    Returns:
    None: The function prints the calculated coordinates (x', y') and the radial distance (r).
    """
    k = np.sqrt(
        2
        / (
            1
            + np.sin(np.radians(latitude_0)) * np.sin(np.radians(latitude))
            + np.cos(np.radians(latitude_0))
            * np.cos(np.radians(latitude))
            * np.cos(np.radians(longitude) - np.radians(longitude_0))
        )
    )
    y_prim = (
        R
        * k
        * np.cos(np.radians(latitude))
        * np.sin(np.radians(longitude) - np.radians(longitude_0))
    )
    x_prim = (
        R
        * k
        * (
            np.cos(np.radians(latitude_0)) * np.sin(np.radians(latitude))
            - np.sin(np.radians(latitude_0))
            * np.cos(np.radians(latitude))
            * np.cos(np.radians(longitude) - np.radians(longitude_0))
        )
    )
    r = np.sqrt(x_prim**2 + y_prim**2)
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
    - latitude_0 (float): The reference latitude in np.degrees.
    - longitude_0 (float): The reference longitude in np.degrees.
    - latitude (float): The latitude in np.degrees.
    - longitude (float): The longitude in np.degrees.
    - R (float): The radius of the sphere.

    Returns:
    None
    """
    c = np.arccos(
        np.sin(np.radians(latitude_0)) * np.sin(np.radians(latitude))
        + np.cos(np.radians(latitude_0))
        * np.cos(np.radians(latitude))
        * np.cos(np.radians(longitude) - np.radians(longitude_0))
    )
    sin_lambda_prim = (
        np.cos(np.radians(latitude)) * np.sin(np.radians(longitude - longitude_0))
    ) / np.sin(c)
    cos_lambda_prim = (
        -np.cos(np.radians(latitude_0)) * np.sin(np.radians(latitude))
        + np.sin(np.radians(latitude_0))
        * np.cos(np.radians(latitude))
        * np.cos(np.radians(longitude - longitude_0))
    ) / np.sin(c)
    x_prim = -R * c * cos_lambda_prim
    y_prim = R * c * sin_lambda_prim
    r = np.sqrt(x_prim**2 + y_prim**2)
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
        latitude_0 (float): Latitude of the projection center in np.degrees.
        longitude_0 (float): Longitude of the projection center in np.degrees.
        latitude (float): Latitude of the point to project in np.degrees.
        longitude (float): Longitude of the point to project in np.degrees.
        R (float): Radius of the sphere in the projection.

    Returns:
        None: The function prints the projected coordinates and the distance from the projection center.

    """
    k = 2 / (
        1
        + np.sin(np.radians(latitude_0)) * np.sin(np.radians(latitude))
        + np.cos(np.radians(latitude_0))
        * np.cos(np.radians(latitude))
        * np.cos(np.radians(longitude) - np.radians(longitude_0))
    )
    x_prim = (
        R
        * k
        * (
            np.cos(np.radians(latitude_0)) * np.sin(np.radians(latitude))
            - np.sin(np.radians(latitude_0))
            * np.cos(np.radians(latitude))
            * np.cos(np.radians(longitude) - np.radians(longitude_0))
        )
    )
    y_prim = (
        R
        * k
        * np.cos(np.radians(latitude))
        * np.sin(np.radians(longitude) - np.radians(longitude_0))
    )
    print("x': ", round(x_prim, 4), "\n", "y': ", round(y_prim, 4))

    r = np.sqrt(x_prim**2 + y_prim**2)
    print("r: ", round(r, 9))


def oblique_gnomic_azimuthal_projection(
    latitude_0, longitude_0, latitude, longitude, R
):
    """
    Calculates the oblique gnomic azimuthal projection coordinates (x', y') and the radial distance (r)
    from the given latitude and longitude to the reference latitude and longitude.

    Parameters:
    - latitude_0 (float): The reference latitude in np.degrees.
    - longitude_0 (float): The reference longitude in np.degrees.
    - latitude (float): The latitude in np.degrees.
    - longitude (float): The longitude in np.degrees.
    - R (float): The radius of the sphere.

    Returns:
    None
    """
    k = 1 / (
        np.sin(np.radians(latitude_0)) * np.sin(np.radians(latitude))
        + np.cos(np.radians(latitude_0))
        * np.cos(np.radians(latitude))
        * np.cos(np.radians(longitude) - np.radians(longitude_0))
    )
    x_prim = (
        R
        * k
        * (
            np.cos(np.radians(latitude_0)) * np.sin(np.radians(latitude))
            - np.sin(np.radians(latitude_0))
            * np.cos(np.radians(latitude))
            * np.cos(np.radians(longitude) - np.radians(longitude_0))
        )
    )
    y_prim = (
        R
        * k
        * np.cos(np.radians(latitude))
        * np.sin(np.radians(longitude) - np.radians(longitude_0))
    )
    print("x': ", round(x_prim, 4), "\n", "y': ", round(y_prim, 4))

    r = np.sqrt(x_prim**2 + y_prim**2)
    print("r: ", round(r, 9))


def abers_projection(
    latitude_1, latitude_2, latitude_0, longitude_0, latitude, longitude, R
):
    """
    Perform Abers projection to convert geographic coordinates to Cartesian coordinates.

    Parameters:
    latitude_1 (float): Latitude of the first standard parallel in np.radians.
    latitude_2 (float): Latitude of the second standard parallel in np.radians.
    latitude_0 (float): Latitude of the projection origin in np.radians.
    longitude_0 (float): Longitude of the projection origin in np.radians.
    latitude (float): Latitude of the point to be projected in np.radians.
    longitude (float): Longitude of the point to be projected in np.radians.
    R (float): Radius of the Earth.

    Returns:
    None
    """
    latitude_1, latitude_2, latitude, longitude, latitude_0, longitude_0 = (
        np.radians(latitude_1),
        np.radians(latitude_2),
        np.radians(latitude),
        np.radians(longitude),
        np.radians(latitude_0),
        np.radians(longitude_0),
    )
    n = (np.sin(latitude_1) + np.sin(latitude_2)) / 2
    C = np.cos(latitude_1) ** 2 + (2 * n * np.sin(latitude_1))
    rho = (R / n) * np.sqrt(C - (2 * n * np.sin(latitude)))
    rho_0 = (R / n) * np.sqrt(C - (2 * n * np.sin(latitude_0)))
    x = rho_0 - (rho * np.cos(n * (longitude - longitude_0)))
    y = rho * np.sin(n * (longitude - longitude_0))
    r = np.sqrt(x**2 + y**2)
    print("x: ", round(x, 4), "\n", "y: ", round(y, 4), "\n", "r: ", round(r, 4))


def lambert_conformal_conical_projection(
    latitude_1, latitude_2, latitude_0, longitude_0, latitude, longitude, R
):
    """
    Performs the Lambert Conformal Conical Projection on a given latitude and longitude.

    Args:
        latitude_1 (float): First standard parallel latitude in np.radians.
        latitude_2 (float): Second standard parallel latitude in np.radians.
        latitude_0 (float): Latitude of the projection origin in np.radians.
        longitude_0 (float): Longitude of the projection origin in np.radians.
        latitude (float): Latitude to project in np.radians.
        longitude (float): Longitude to project in np.radians.
        R (float): Radius of the Earth.

    Returns:
        None
    """
    latitude_1, latitude_2, latitude, longitude, latitude_0, longitude_0 = (
        np.radians(latitude_1),
        np.radians(latitude_2),
        np.radians(latitude),
        np.radians(longitude),
        np.radians(latitude_0),
        np.radians(longitude_0),
    )
    n = (np.log(np.cos(latitude_1) * (1 / np.cos(latitude_2)))) / (
        np.log(
            np.tan((np.pi / 4) + (latitude_2 / 2))
            * (1 / np.tan((np.pi / 4) + (latitude_1 / 2)))
        )
    )
    F = (np.cos(latitude_1) / n) * (np.tan((np.pi / 4) + (latitude_1 / 2)) ** n)
    rho = (R * F) / (np.tan((np.pi / 4) + (latitude / 2)) ** n)
    rho_0 = (R * F) / (np.tan((np.pi / 4) + (latitude_0 / 2)) ** n)
    x = rho_0 - (rho * np.cos(n * (longitude - longitude_0)))
    y = rho * np.sin(n * (longitude - longitude_0))
    r = np.sqrt(x**2 + y**2)
    print("x: ", round(x, 4), "\n", "y: ", round(y, 4), "\n", "r: ", round(r, 4))
    return x, y, r


def normal_mercator_projection_for_sphere(
    latitude_0, longitude_0, latitude, longitude, R
):
    """
    Calculates the normal Mercator projection coordinates (x', y') and the radial distance (r)
    from the given latitude and longitude to the reference latitude and longitude.

    Parameters:
    - latitude_0 (float): The reference latitude in np.degrees.
    - longitude_0 (float): The reference longitude in np.degrees.
    - latitude (float): The latitude in np.degrees.
    - longitude (float): The longitude in np.degrees.
    - R (float): The radius of the sphere.

    Returns:
    None
    """

    x = (
        R
        * np.cos(np.radians(latitude_0))
        * np.log(np.tan((np.pi / 4) + (np.radians(latitude) / 2)))
    )
    y = R * np.cos(np.radians(latitude_0)) * np.radians(longitude - longitude_0)

    print("x: ", round(x, 4), "\n", "y: ", round(y, 4))
    return x, y


def normal_mercator_projection_for_ellipsoid(
    latitude_0, longitude_0, latitude, longitude, ellipsoid
):
    """
    Calculates the normal Mercator projection for an ellipsoid.

    Args:
        latitude_0 (float): The reference latitude in np.degrees.
        longitude_0 (float): The reference longitude in np.degrees.
        latitude (float): The latitude in np.degrees.
        longitude (float): The longitude in np.degrees.
        ellipsoid (tuple): A tuple containing the semi-major axis (a)
            and the reciprocal of the flattening (rf) of the ellipsoid.

    Returns:
        tuple: A tuple containing the projected x-coordinate, projected y-coordinate, scale factor (h),
            convergence factor (k), and eccentricity squared (epsilon).

    """

    latitude_0, longitude_0, latitude, longitude = (
        np.radians(latitude_0),
        np.radians(longitude_0),
        np.radians(latitude),
        np.radians(longitude),
    )

    # Input parameters
    a, rf = ellipsoid
    e2 = 1 - (1 - 1 / rf) ** 2
    e = np.sqrt(e2)
    C = (a * np.cos(latitude_0)) / (np.sqrt(1 - (e2 * (np.sin(latitude_0) ** 2))))
    x = C * np.log(
        np.tan(np.pi / 4 + latitude / 2)
        * ((1 - e * np.sin(latitude)) / (1 + e * np.sin(latitude))) ** (e / 2)
    )
    print("x: ", round(x, 4))
    y = C * (longitude - longitude_0)
    print("y: ", round(y, 4))

    h = k = (
        np.cos(latitude_0)
        / np.cos(latitude)
        * np.sqrt((1 - e2 * np.sin(latitude) ** 2) / (1 - e2 * np.sin(latitude_0) ** 2))
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
    - latitude_0 (float): Reference latitude in np.radians.
    - longitude_0 (float): Reference longitude in np.radians.
    - latitude (float): Latitude in np.radians.
    - longitude (float): Longitude in np.radians.
    - R (float): Radius of the sphere.

    Returns:
    - x (float): x-coordinate.
    - y (float): y-coordinate.
    - r (float): Distance from the origin.
    - h (float): Scale factor in the x-direction.
    - k (float): Scale factor in the y-direction.
    """

    latitude_0, longitude_0, latitude, longitude = (
        np.radians(latitude_0),
        np.radians(longitude_0),
        np.radians(latitude),
        np.radians(longitude),
    )

    B = np.cos(latitude) * np.sin(longitude - longitude_0)
    x = R * (np.arctan(np.tan(latitude) / np.cos(longitude - longitude_0)) - latitude_0)
    y = 0.5 * R * np.log((1 + B) / (1 - B))
    h = k = 1 / np.sqrt(1 - (B**2))

    return x, y, h, k


def transverse_mercator_projection_for_ellipsoid(
    latitude_0, longitude_0, latitude, longitude, ellipsoid, scale_factor, false_easting
):
    """
    Transforms geographic coordinates (latitude, longitude) to projected coordinates (x, y)
    using the Transverse Mercator projection for a given ellipsoid.

    Args:
        latitude_0 (float): Latitude of the projection origin in np.degrees.
        longitude_0 (float): Longitude of the projection origin in np.degrees.
        latitude (float): Latitude of the point to be projected in np.degrees.
        longitude (float): Longitude of the point to be projected in np.degrees.
        ellipsoid (tuple): Tuple containing the semi-major axis (a) and reciprocal flattening (rf) of the ellipsoid.
        scale_factor (float): Scale factor of the projection.

    Returns:
        tuple: A tuple containing the projected coordinates (x, y).
    """

    a, rf = ellipsoid
    e2 = 1 - (1 - 1 / rf) ** 2
    b = a * (1 - 1 / rf)
    latitude_0, latitude = np.radians(latitude_0), np.radians(latitude)

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
        np.sin(latitude)
        * np.cos(latitude)
        * (
            A
            + B * np.sin(latitude) ** 2
            + C * np.sin(latitude) ** 4
            + D * np.sin(latitude) ** 6
        )
    )

    delta_longitude = np.radians(longitude - longitude_0)

    epsilon_prim = np.arctan(np.tan(latitude_temp) / np.cos(delta_longitude))

    eta_prim = np.arctanh(np.cos(latitude_temp) * np.sin(delta_longitude))

    x = (
        a_hat
        * scale_factor
        * (
            epsilon_prim
            + beta_1 * np.sin(2 * epsilon_prim) * np.cosh(2 * eta_prim)
            + beta_2 * np.sin(4 * epsilon_prim) * np.cosh(4 * eta_prim)
            + beta_3 * np.sin(6 * epsilon_prim) * np.cosh(6 * eta_prim)
            + beta_4 * np.sin(8 * epsilon_prim) * np.cosh(8 * eta_prim)
        )
    )
    y = (
        a_hat
        * scale_factor
        * (
            eta_prim
            + beta_1 * np.cos(2 * epsilon_prim) * np.sinh(2 * eta_prim)
            + beta_2 * np.cos(4 * epsilon_prim) * np.sinh(4 * eta_prim)
            + beta_3 * np.cos(6 * epsilon_prim) * np.sinh(6 * eta_prim)
            + beta_4 * np.cos(8 * epsilon_prim) * np.sinh(8 * eta_prim)
        )
        + false_easting  # 500 km false easting
    )
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
        degrees (float): The input np.degrees.

    Returns:
        tuple: A tuple containing the degrees, minutes, and seconds.
    """
    degrees = abs(degrees)
    deg_int = int(degrees)
    minutes = int((degrees - deg_int) * 60)
    seconds = (degrees - deg_int - (minutes / 60)) * 3600
    return deg_int, minutes, seconds


def R1_transformation(alpha):
    """
    Calculates the R1 transformation matrix.

    Args:
        alpha (float): The rotation angle in np.radians.

    Returns:
        numpy.ndarray: The R1 transformation matrix.
    """
    return np.array(
        [
            [1, 0, 0],
            [0, np.cos(alpha), np.sin(alpha)],
            [0, -np.sin(alpha), np.cos(alpha)],
        ]
    )


def R2_transformation(alpha):
    """
    Calculates the R2 transformation matrix.

    Args:
        alpha (float): The rotation angle in np.radians.

    Returns:
        numpy.ndarray: The R2 transformation matrix.
    """
    return np.array(
        [
            [np.cos(alpha), 0, -np.sin(alpha)],
            [0, 1, 0],
            [np.sin(alpha), 0, np.cos(alpha)],
        ]
    )


def R3_transformation(alpha):
    """
    Calculates the R3 transformation matrix.

    Args:
        alpha (float): The rotation angle in np.radians.

    Returns:
        numpy.ndarray: The R3 transformation matrix.
    """
    return np.array(
        [
            [np.cos(alpha), np.sin(alpha), 0],
            [-np.sin(alpha), np.cos(alpha), 0],
            [0, 0, 1],
        ]
    )


def R_transformation(alpha1, alpha2, alpha3):
    """
    Calculates the R transformation matrix.

    Args:
        alpha1 (float): The rotation angle in radians for the first rotation.
        alpha2 (float): The rotation angle in radians for the second rotation.
        alpha3 (float): The rotation angle in radians for the third rotation.

    Returns:
        numpy.ndarray: The R transformation matrix.
    """
    return (
        R3_transformation(alpha3)
        @ R2_transformation(alpha2)
        @ R1_transformation(alpha1)
    )


# def Rp_transformation():


if __name__ == "__main__":
    print(convert_to_degrees_minutes_seconds(60.59444444444444))
