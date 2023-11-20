package geodesic

import (
	"math"

	"github.com/pymaxion/geographiclib-go/geodesic/capabilities"
)

/*
 Gnomonic is a type representing calculations for gnomonic projections centered
 at an arbitrary position C on the ellipsoid. This projection is derived in
 Section 8 of C. F. F. Karney, Algorithms for geodesics, J. Geodesy 87, 43–55
 (2013); DOI: 10.1007/s00190-012-0578-z; addenda:
 https://geographiclib.sourceforge.io/geod-addenda.html.

 The gnomonic projection of a point P on the ellipsoid is defined as follows:
 compute the geodesic line from C to P; compute the reduced length m12, geodesic
 scale M12, and ρ = m12/M12; finally, this gives the coordinates x and y of P in
 gnomonic projection with x = ρ sin azi1; y = ρ cos azi1, where azi1 is the
 azimuth of the geodesic at C. The method Forward(double, double, double, double)
 performs the forward projection and Reverse(double, double, double, double) is
 the inverse of the projection. The methods also return the azimuth azi of the
 geodesic at P and reciprocal scale rk in the azimuthal direction. The scale in
 the radial direction is 1/(rk^2).

 For a sphere, ρ reduces to a tan(s12/a), where s12 is the length of the geodesic
 from C to P, and the gnomonic projection has the property that all geodesics
 appear as straight lines. For an ellipsoid, this property holds only for
 geodesics intersecting the centers. However geodesic segments close to the center
 are approximately straight.

 Consider a geodesic segment of length l. Let T be the point on the geodesic
 (extended if necessary) closest to C, the center of the projection, and t, be
 the distance CT. To lowest order, the maximum deviation (as a true distance) of
 the corresponding gnomonic line segment (i.e., with the same end points) from
 the geodesic is:

  K(T) − K(C)) l^2 t / 32.

 where K is the Gaussian curvature.

 This result applies for any surface. For an ellipsoid of revolution, consider
 all geodesics whose end points are within a distance r of C. For a given r, the
 deviation is maximum when the latitude of C is 45°, when endpoints are a
 distance r away, and when their azimuths from the center are ± 45° or ± 135°. To
 lowest order in r and the flattening f, the deviation is f (r/2a)^3 r.

 CAUTION: The definition of this projection for a sphere is standard. However,
 there is no standard for how it should be extended to an ellipsoid. The choices
 are:
  * Declare that the projection is undefined for an ellipsoid.
  * Project to a tangent plane from the center of the ellipsoid. This causes great
  ellipses to appear as straight lines in the projection; i.e., it generalizes the
  spherical great circle to a great ellipse. This was proposed by independently by
  Bowring and Williams in 1997.
  * Project to the conformal sphere with the constant of integration chosen so
  that the values of the latitude match for the center point and perform a central
  projection onto the plane tangent to the conformal sphere at the center point.
  This causes normal sections through the center point to appear as straight lines
  in the projection; i.e., it generalizes the spherical great circle to a normal
  section. This was proposed by I. G. Letoval'tsev, Generalization of the gnomonic
  projection for a spheroid and the principal geodetic problems involved in the
  alignment of surface routes, Geodesy and Aerophotography (5), 271–274 (1963).
  * The projection given here. This causes geodesics close to the center point to
  appear as straight lines in the projection; i.e., it generalizes the spherical
  great circle to a geodesic.

 The algorithms are described in
  C. F. F. Karney, Algorithms for geodesics, J. Geodesy 87, 43-55 (2013)
  Link: https://doi.org/10.1007/s00190-012-0578-z
  Addenda: https://geographiclib.sourceforge.io/geod-addenda.html
*/
type Gnomonic struct {
	Earth *Geodesic
}

/*
 Forward performs the forward projection from geographic to gnomonic.

  lat0: latitude of center point of projection (degrees). Should be in the range [-90°, 90°].
  lon0: longitude of center point of projection (degrees). Should be in the range [−540°, 540°).
  lat - latitude of point (degrees). Should be in the range [-90°, 90°].
  lon - longitude of point (degrees). Should be in the range [−540°, 540°).

 The scale of the projection is 1/rk2 in the "radial" direction, azi clockwise
 from true north, and is 1/rk in the direction perpendicular to this. If the
 point lies "over the horizon", i.e., if rk ≤ 0, then NaNs are returned for x and
 y (the correct values are returned for azi and rk). A call to Forward followed
 by a call to Reverse will return the original (lat, lon) (to within roundoff)
 provided the point is not over the horizon.
*/
func (g *Gnomonic) Forward(lat0, lon0, lat, lon float64) GnomonicData {
	fwd := newGnomonicData()
	fwd.Lat0, fwd.Lon0, fwd.Lat, fwd.Lon = lat0, lon0, lat, lon

	caps := capabilities.Azimuth | capabilities.GeodesicScale | capabilities.ReducedLength
	inv := g.Earth.InverseWithCapabilities(lat0, lon0, lat, lon, caps)
	fwd.Azi, fwd.Rk = inv.Azi2, inv.M12

	if inv.M12 > 0 {
		rho := inv.M12Reduced / inv.M12
		x, y := sincosd(inv.Azi1)
		fwd.X, fwd.Y = rho*x, rho*y
	}
	return fwd
}

/*
 Reverse performs the reverse projection from gnomonic to geographic.

  lat0: latitude of center point of projection (degrees). Should be in the range [-90°, 90°].
  lon0: longitude of center point of projection (degrees). Should be in the range [−540°, 540°).
  x: easting of point (meters).
  y: northing of point (meters).

 The returned values GnomonicData.Lat and GnomonicData.Lon are in the range [-180°, 180°].

 The scale of the projection is 1/rk2 in the "radial" direction, azi clockwise
 from true north, and is 1/rk in the direction perpendicular to this. Even though
 all inputs should return a valid lat and lon, it's possible that the procedure
 fails to converge for very large x or y; in this case NaNs are returned for all
 the output arguments. A call to Reverse followed by a call to Forward will
 return the original (x, y) (to roundoff).
*/
func (g *Gnomonic) Reverse(lat0, lon0, x, y float64) GnomonicData {
	rev := newGnomonicData()
	rev.Lat0, rev.Lon0, rev.X, rev.Y = lat0, lon0, x, y

	azi0 := atan2d(x, y)
	rho := math.Hypot(x, y)
	a := g.Earth.EquatorialRadius()
	s := a * math.Atan(rho/a)

	little := rho <= a
	if !little {
		rho = 1 / rho
	}

	caps := capabilities.Latitude | capabilities.Longitude | capabilities.Azimuth | capabilities.DistanceIn |
		capabilities.ReducedLength | capabilities.GeodesicScale
	line := g.Earth.LineWithCapabilities(lat0, lon0, azi0, caps)

	var pos Data
	numIt := 10
	trip := 0
	tripEpsilon := 0.01 * math.Sqrt(epsilon)

	for count := numIt; count > 0; count-- {
		pos = line.PositionWithCapabilities(s, caps)
		if trip > 0 {
			break
		}

		var ds float64
		if little {
			ds = ((pos.M12Reduced / pos.M12) - rho) * pos.M12 * pos.M12
		} else {
			ds = (rho - (pos.M12 / pos.M12Reduced)) * pos.M12Reduced * pos.M12Reduced
		}
		s -= ds

		if math.Abs(ds) < tripEpsilon*a {
			trip++
		}
	}

	if trip == 0 {
		return rev
	}

	rev.Lat, rev.Lon, rev.Azi, rev.Rk = pos.Lat2, pos.Lon2, pos.Azi1, pos.M12
	return rev
}

// GnomonicData describes the results of gnomonic projection. This is used to
// return the results for a gnomonic projection of a point (lat, lon) given a
// center point of projection (lat0, lon0). The returned GnomonicData structs
// always include the parameters provided to Gnomonic.Forward and
// Gnomonic.Reverse and it always includes the fields x, y, azi. and rk.
type GnomonicData struct {
	Lat0 float64
	Lon0 float64
	Lat  float64
	Lon  float64
	X    float64
	Y    float64
	Azi  float64
	Rk   float64
}

func newGnomonicData() GnomonicData {
	return GnomonicData{
		Lat0: math.NaN(),
		Lon0: math.NaN(),
		Lat:  math.NaN(),
		Lon:  math.NaN(),
		X:    math.NaN(),
		Y:    math.NaN(),
		Azi:  math.NaN(),
		Rk:   math.NaN(),
	}
}
