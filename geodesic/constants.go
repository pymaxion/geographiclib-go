package geodesic

import (
	"math"
)

// Constants defining the WGS84 ellipsoid

//goland:noinspection ALL
const (
	// WGS84_a is the equitorial radius of the WGS84 ellipsoid (6378137 m)
	WGS84_a = 6378137.

	// WGS84_f is the flattening of the WGS84 ellipsoid (1/298.257223563)
	WGS84_f = 1 / 298.257223563
)

const (
	geodesicOrder = 6
	nA1           = geodesicOrder
	nC1           = geodesicOrder
	nC1p          = geodesicOrder
	nA2           = geodesicOrder
	nC2           = geodesicOrder
	nA3           = geodesicOrder
	nA3x          = nA3
	nC3           = geodesicOrder
	nC3x          = (nC3 * (nC3 - 1)) / 2
	nC4           = geodesicOrder
	nC4x          = (nC4 * (nC4 + 1)) / 2
	maxit1        = 20
	maxit2        = maxit1 + digits + 10
)

//goland:noinspection ALL
var (
	// WGS84 is an instance of Geodesic which represents the Earth according to the World Geodetic
	// System. See also: https://en.wikipedia.org/wiki/World_Geodetic_System
	WGS84, _ = NewGeodesic(WGS84_a, WGS84_f)
)

var (
	// using 2^-1022 instead of math.SmallestNonzeroFloat64 to maintain consistency with other
	// geographiclib implementations
	tiny    = math.Sqrt(math.Pow(2, -1022))
	epsilon = math.Nextafter(1., 2.) - 1. // https://stackoverflow.com/a/22185792/4755732
	tol0    = epsilon
	tol1    = 200 * tol0
	tol2    = math.Sqrt(tol0)
	tolb    = tol0 * tol2
	xthresh = 1000 * tol2
)
