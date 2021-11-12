package geodesic

import (
	"errors"
	"math"

	"geographiclib-go/geodesic/capabilities"
)

/*
 Geodesic is an interface representing geodesic calculations.

 The shortest path between two points on a ellipsoid at (lat1, lon1) and (lat2, lon2) is called the
 geodesic. Its length is s12 and the geodesic from point 1 to point 2 has azimuths azi1 and azi2 at
 the two end points. (The azimuth is the heading measured clockwise from north. azi2 is the
 "forward" azimuth, i.e., the heading that takes you beyond point 2 not back to point 1.)

 Given lat1, lon1, azi1, and s12, we can determine lat2, lon2, and azi2. This is the direct geodesic
 problem and its solution is given by the function Direct. (If s12 is sufficiently large that the
 geodesic wraps more than halfway around the earth, there will be another geodesic between the
 points with a smaller s12.)

 Given lat1, lon1, lat2, and lon2, we can determine azi1, azi2, and s12. This is the inverse
 geodesic problem, whose solution is given by Inverse. Usually, the solution to the inverse problem
 is unique. In cases where there are multiple solutions (all with the same s12, of course), all the
 solutions can be easily generated once a particular solution is provided.

 The standard way of specifying the direct problem is to specify the distance s12 to the second
 point. However it is sometimes useful instead to specify the arc length a12 (in degrees) on the
 auxiliary sphere. This is a mathematical construct used in solving the geodesic problems. The
 solution of the direct problem in this form is provided by ArcDirect. An arc length in excess of
 180° indicates that the geodesic is not a shortest path. In addition, the arc length between an
 equatorial crossing and the next extremum of latitude for a geodesic is 90°.

 This interface can also calculate several other quantities related to geodesics. These are:

 Reduced length:
  If we fix the first point and increase azi1 by dazi1 (radians), the second point is displaced m12
  dazi1 in the direction azi2 + 90°. The quantity m12 is called the "reduced length" and is
  symmetric under interchange of the two points. On a curved surface the reduced length obeys a
  symmetry relation, m12 + m21 = 0. On a flat surface, we have m12 = s12. The ratio s12/m12 gives
  the azimuthal scale for an azimuthal equidistant projection.

 Geodesic scale:
  Consider a reference geodesic and a second geodesic parallel to this one at point 1 and separated
  by a small distance dt. The separation of the two geodesics at point 2 is M12 dt where M12 is
  called the "geodesic scale". M21 is defined similarly (with the geodesics being parallel at point
  2). On a flat surface, we have M12 = M21 = 1. The quantity 1/M12 gives the scale of the
  Cassini-Soldner projection.

 Area:
  The area between the geodesic from point 1 to point 2 and the equation is represented by S12; it
  is the area, measured counter-clockwise, of the geodesic quadrilateral with corners (lat1,lon1),
  (0,lon1), (0,lon2), and (lat2,lon2). It can be used to compute the area of any geodesic polygon.

 The quantities m12, M12, M21, which all specify the behavior of nearby geodesics, obey addition
 rules. If points 1, 2, and 3 all lie on a single geodesic, then the following rules hold:

  s13 = s12 + s23
  a13 = a12 + a23
  S13 = S12 + S23
  m13 = m12 M23 + m23 M21
  M13 = M12 M23 - (1 - M12 M21) m23 / m12
  M31 = M32 M21 - (1 - M23 M32) m12 / m23

 The results of the geodesic calculations are bundled up into a Data object which includes the input
 parameters and all the computed results, i.e., lat1, lon1, azi1, lat2, lon2, azi2, s12, a12, m12
 (as M12Reduced), M12, M21, S12 (as S12Area).

 The functions DirectWithCapabilities, InverseWithCapabilities, and ArcDirectWithCapabilities
 include an optional final argument of type capabilities.BitMask to allow you to specify which
 results should be computed and returned. The default functions Direct, Inverse, and ArcDirect are
 equivalent to calling DirectWithCapabilities, InverseWithCapabilities, and
 ArcDirectWithCapabilities with the capabilities.Standard argument, which will compute "standard"
 geodesic results (latitudes, longitudes, azimuths, and distance). A custom capabilities.BitMask
 argument can be bitor'ed out of a combination of other capabilities.BitMask values. For example, if
 you wish to only compute the distance between two points, you would call, e.g.,

  geodesicData := WGS84.Inverse(lat1, lon1, lat2, lon2, capabilities.Distance)

 The calculations are accurate to better than 15 nm (15 nanometers) for the WGS84 ellipsoid. See
 Sec. 9 of arXiv:1102.1215v1 (https://arxiv.org/abs/1102.1215v1) for details. The algorithms used by
 this interface are based on series expansions using the flattening f as a small parameter. These
 are only accurate for |f| < 0.02; however, reasonably accurate results will be obtained for |f| <
 0.2. Here is a table of the approximate maximum error (expressed as a distance) for an ellipsoid
 with the same equatorial radius as the WGS84 ellipsoid and different values of the flattening.
  |f|      error
  0.01     25 nm
  0.02     30 nm
  0.05     10 um
  0.1     1.5 mm
  0.2     300 mm

 The algorithms are described in
  C. F. F. Karney, Algorithms for geodesics, J. Geodesy 87, 43-55 (2013)
  Link: https://doi.org/10.1007/s00190-012-0578-z
  Addenda: https://geographiclib.sourceforge.io/geod-addenda.html
*/
type Geodesic interface {
	Direct(lat1, lon1, azi1, s12 float64) Data
	DirectWithCapabilities(lat1, lon1, azi1, s12 float64, mask capabilities.BitMask) Data
	ArcDirect(lat1, lon1, azi1, a12 float64) Data
	ArcDirectWithCapabilities(lat1, lon1, azi1, a12 float64, mask capabilities.BitMask) Data
	Inverse(lat1, lon1, lat2, lon2 float64) Data
	InverseWithCapabilities(lat1, lon1, lat2, lon2 float64, mask capabilities.BitMask) Data
}

// NewGeodesic creates a new instance of Geodesic with the given equatorial radius a (in meters) and
// flattening of the ellipsoid f. Setting f = 0 gives a sphere; negative f gives a prolate
// ellipsoid.
//
// This function will return an error if either a or (1-f) is not positive.
func NewGeodesic(a, f float64) (Geodesic, error) {
	return newGeodesicImpl(a, f)
}

// geodesicImpl is an unexported implementation of the Geodesic interface
type geodesicImpl struct {
	a     float64
	f     float64
	f1    float64
	e2    float64
	ep2   float64
	b     float64
	c2    float64
	n     float64
	etol2 float64
	a3x   []float64
	c3x   []float64
	c4x   []float64
}

func newGeodesicImpl(a, f float64) (*geodesicImpl, error) {
	if !(isfinite(a) && a > 0) {
		return nil, errors.New("equatorial radius is not positive")
	}
	f1 := 1 - f
	b := a * f1
	if !(isfinite(b) && b > 0) {
		return nil, errors.New("polar semi-axis is not positive")
	}
	e2 := f * (2 - f)
	ep2 := e2 / sq(f1) // e2 / 1 - e2
	n := f / (2 - f)
	// authalic radius squared
	c2 := calculateAuthalicRadiusSquared(a, b, e2)
	// The sig12 threshold for "really short". Using the auxiliary sphere solution with dnm computed
	// at (bet1 + bet2) / 2, the relative error in the azimuth consistency check is sig12^2 * abs(f)
	// * min(1, 1-f/2) / 2. (Error measured for 1/100 < b/a < 100 and abs(f) >= 1/1000. For a given f
	// and sig12, the max error occurs for lines near the pole. If the old rule for computing dnm =
	// (dn1 + dn2)/2 is used, then the error increases by a factor of 2.) Setting this equal to
	// epsilon gives sig12 = etol2. Here 0.1 is a safety factor (error decreased by 100) and
	// max(0.001, abs(f)) stops etol2 getting too large in the nearly spherical case.
	etol2 := 0.1 * tol2 / math.Sqrt(math.Max(0.001, math.Abs(f))*math.Min(1, 1-f/2)/2)
	a3x := initA3x(n)
	c3x := initC3x(n)
	c4x := initC4x(n)
	return &geodesicImpl{
		a:     a,
		f:     f,
		f1:    f1,
		e2:    e2,
		ep2:   ep2,
		b:     b,
		c2:    c2,
		n:     n,
		etol2: etol2,
		a3x:   a3x,
		c3x:   c3x,
		c4x:   c4x,
	}, nil
}

func calculateAuthalicRadiusSquared(a, b, e2 float64) float64 {
	var multiplier float64
	if e2 == 0 {
		multiplier = 1
	} else {
		var dividend float64
		if e2 > 0 {
			dividend = atanh(math.Sqrt(e2))
		} else {
			dividend = math.Atan(math.Sqrt(-e2))
		}
		divisor := math.Sqrt(math.Abs(e2))
		multiplier = dividend / divisor
	}
	return (sq(a) + sq(b)*multiplier) / 2
}

// a3f evaluates A3
func (g *geodesicImpl) a3f(eps float64) float64 {
	return polyval(nA3-1, g.a3x, 0, eps)
}

// c3f evaluates C3 coefficients and sets elements c[1] thru c[nC3 - 1]
func (g *geodesicImpl) c3f(eps float64, c []float64) {
	mult := 1.
	o := 0
	for l := 1; l < nC3; l++ { // l is index of C3[l]
		m := nC3 - l - 1 // order of polynomial in eps
		mult *= eps
		c[l] = mult * polyval(m, g.c3x, o, eps)
		o += m + 1
	}
}

// c4f evaluates C4 coefficients and sets elements c[0] thru c[nC4 - 1]
func (g *geodesicImpl) c4f(eps float64, c []float64) {
	mult := 1.
	o := 0
	for l := 0; l < nC4; l++ { // l is index of C4[l]
		m := nC4 - l - 1 // order of polynomial in eps
		c[l] = mult * polyval(m, g.c4x, o, eps)
		o += m + 1
		mult *= eps
	}
}

/*
 Direct solves the direct geodesic problem where the length of the geodesic is specified in terms of
 distance.

  lat1: latitude of point 1 (degrees). Should be in the range [-90°, 90°].
  lon1: longitude of point 1 (degrees).
  azi1: azimuth at point 1 (degrees).
  s12: distance between point 1 and point 2 (meters); can be negative.

 The returned values Data.Lon2 and Data.Azi2 are in the range [-180°, 180°].

 If either point is at a pole, the azimuth is defined by keeping the longitude fixed, writing lat =
 ±(90° - ε), and taking the limit ε → 0+. An arc length greater that 180° signifies a geodesic which
 is not a shortest path. (For a prolate ellipsoid, an additional condition is necessary for a
 shortest path: the longitudinal extent must not exceed of 180°.)

 This function is equivalent to calling DirectWithCapabilities with capabilities.Standard.
*/
func (g *geodesicImpl) Direct(lat1, lon1, azi1, s12 float64) Data {
	return g.DirectWithCapabilities(lat1, lon1, azi1, s12, capabilities.Standard)
}

/*
 DirectWithCapabilities solves the direct geodesic problem where the length of the geodesic is
 specified in terms of distance. It also allows you to specify which results should be computed and
 returned via the capabilities.BitMask argument.

 See Direct for more details.
*/
func (g *geodesicImpl) DirectWithCapabilities(lat1, lon1, azi1, s12 float64, mask capabilities.BitMask) Data {
	panic("not yet implemented")
}

/*
 ArcDirect solves the direct geodesic problem where the length of the geodesic is specified in terms
 of arc length.

  lat1: latitude of point 1 (degrees). Should be in the range [-90°, 90°].
  lon1: longitude of point 1 (degrees).
  azi1: azimuth at point 1 (degrees).
  a12: arc length between point 1 and point 2 (meters); can be negative.

 The returned values Data.Lon2 and Data.Azi2 are in the range [-180°, 180°].

 If either point is at a pole, the azimuth is defined by keeping the longitude fixed, writing lat =
 ±(90° - ε), and taking the limit ε → 0+. An arc length greater that 180° signifies a geodesic which
 is not a shortest path. (For a prolate ellipsoid, an additional condition is necessary for a
 shortest path: the longitudinal extent must not exceed of 180°.)

 This function is equivalent to calling ArcDirectWithCapabilities with capabilities.Standard.
*/
func (g *geodesicImpl) ArcDirect(lat1, lon1, azi1, a12 float64) Data {
	return g.ArcDirectWithCapabilities(lat1, lon1, azi1, a12, capabilities.Standard)
}

/*
 ArcDirectWithCapabilities solves the direct geodesic problem where the length of the geodesic is
 specified in terms of arc length. It also allows you to specify which results should be computed
 and returned via the capabilities.BitMask argument.

 See ArcDirect for more details.
*/
func (g *geodesicImpl) ArcDirectWithCapabilities(lat1, lon1, azi1, a12 float64, mask capabilities.BitMask) Data {
	panic("not yet implemented")
}

/*
 Inverse solves the inverse geodesic problem.

  lat1 latitude of point 1 (degrees). Should be in the range [-90°, 90°].
  lon1 longitude of point 1 (degrees).
  lat2 latitude of point 2 (degrees). Should be in the range [-90°, 90°].
  lon2 longitude of point 2 (degrees).

 The returned values Data.Azi1 and Data.Azi2 are in the range [-180°, 180°].

 If either point is at a pole, the azimuth is defined by keeping the longitude fixed, writing lat =
 ±(90° - ε), taking the limit ε → 0+.

 The solution to the inverse problem is found using Newton's method. If this fails to converge (this
 is very unlikely in geodetic applications but does occur for very eccentric ellipsoids), then the
 bisection method is used to refine the solution.

 This function is equivalent to calling InverseWithCapabilities with capabilities.Standard.
*/
func (g *geodesicImpl) Inverse(lat1, lon1, lat2, lon2 float64) Data {
	return g.InverseWithCapabilities(lat1, lon1, lat2, lon2, capabilities.Standard)
}

/*
 InverseWithCapabilities solves the inverse geodesic problem. It also allows you to specify which
 results should be computed and returned via the capabilities.BitMask argument.

 See Inverse for more details.
*/
func (g *geodesicImpl) InverseWithCapabilities(lat1, lon1, lat2, lon2 float64, mask capabilities.BitMask) Data {
	solver := inverseSolver{g}
	return solver.inverse(lat1, lon1, lat2, lon2, mask)
}
