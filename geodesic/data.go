package geodesic

import "math"

// Data describes the characteristics of a geodesic between point 1 (Lat1, Lon1)
// and point 2 (Lat2, Lon2). Fields that have not been set will be filled by
// math.NaN. Data will always include the field A12.
type Data struct {

	// Lat1 is the latitude of point 1 in degrees
	Lat1 float64

	// Lon1 is the longitude of point 1 in degrees
	Lon1 float64

	// Azi1 is the azimuth at point 1 in degrees
	Azi1 float64

	// Lat2 is the latitude of point 2 in degrees
	Lat2 float64

	// Lon2 is the longitude of point 2 in degrees
	Lon2 float64

	// Azi2 is the azimuth at point 2 in degrees
	Azi2 float64

	// S12 is the distance between point 1 and point 2 in meters
	S12 float64

	// A12 is the arc length on the auxiliary sphere between point 1 and point 2 in
	// degrees
	A12 float64

	// M12Reduced is the reduced length of the geodesic in meters.
	//
	// Reduced length:
	//  If we fix the first point and increase azi1 by dazi1 (radians), the second
	//  point is displaced m12 dazi1 in the direction azi2 + 90 degrees. The quantity
	//  m12 is called the "reduced length" and is symmetric under interchange of the
	//  two points. On a curved surface the reduced length obeys a symmetry relation,
	//  m12 + m21 = 0. On a flat surface, we have m12 = s12. The ratio s12/m12 gives
	//  the azimuthal scale for an azimuthal equidistant projection.
	M12Reduced float64

	// M12 is the geodesic scale of point 2 relative to point 1 (dimensionless)
	//
	// Geodesic scale:
	//  Consider a reference geodesic and a second geodesic parallel to this one at
	//  point 1 and separated by a small distance dt. The separation of the two
	//  geodesics at point 2 is M12 dt where M12 is called the "geodesic scale". M21
	//  is defined similarly (with the geodesics being parallel at point 2). On a flat
	//  surface, we have M12 = M21 = 1. The quantity 1/M12 gives the scale of the
	//  Cassini-Soldner projection.
	M12 float64

	// M21 is the geodesic scale of point 1 relative to point 2 (dimensionless)
	//
	// Geodesic scale:
	//  Consider a reference geodesic and a second geodesic parallel to this one at
	//  point 1 and separated by a small distance dt. The separation of the two
	//  geodesics at point 2 is M12 dt where M12 is called the "geodesic scale". M21
	//  is defined similarly (with the geodesics being parallel at point 2). On a flat
	//  surface, we have M12 = M21 = 1. The quantity 1/M12 gives the scale of the
	//  Cassini-Soldner projection.
	M21 float64

	// S12Area is the area under the geodesic (meters²)
	//
	// Area:
	//  The area between the geodesic from point 1 to point 2 and the equation is
	//  represented by S12Area; it is the area, measured counter-clockwise, of the
	//  geodesic quadrilateral with corners (lat1,lon1), (0,lon1), (0,lon2), and
	//  (lat2,lon2). It can be used to compute the area of any geodesic polygon.
	S12Area float64
}

func newData() Data {
	return Data{
		Lat1:       math.NaN(),
		Lon1:       math.NaN(),
		Azi1:       math.NaN(),
		Lat2:       math.NaN(),
		Lon2:       math.NaN(),
		Azi2:       math.NaN(),
		S12:        math.NaN(),
		A12:        math.NaN(),
		M12Reduced: math.NaN(),
		M12:        math.NaN(),
		M21:        math.NaN(),
		S12Area:    math.NaN(),
	}
}
