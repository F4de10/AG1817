
from math import radians, sin, cos, tan, atanh, atan, sinh, cosh


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