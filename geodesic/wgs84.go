package geodesic

//goland:noinspection ALL
const (
	// WGS84_a is the equatorial radius of the WGS84 ellipsoid (6378137 m)
	WGS84_a = 6378137.

	// WGS84_f is the flattening of the WGS84 ellipsoid (1/298.257223563)
	WGS84_f = 1 / 298.257223563
)

var (
	// WGS84 is an instance of Geodesic which represents the Earth according to the
	// World Geodetic System. See also:
	// https://en.wikipedia.org/wiki/World_Geodetic_System
	WGS84, _ = NewGeodesic(WGS84_a, WGS84_f)
)
