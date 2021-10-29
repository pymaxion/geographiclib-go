package geographiclib

import "math"

// GeodesicData describes the characteristics of a geodesic between point 1 (Lat1, Lon1) and point 2
// (Lat2, Lon2). Fields that have not been set will be filled by math.NaN. GeodesicData will always
// include the field A12.
type GeodesicData interface {

	// Lat1 is the latitude of point 1 in degrees
	Lat1() float64

	// Lon1 is the longitude of point 1 in degrees
	Lon1() float64

	// Azi1 is the azimuth at point 1 in degrees
	Azi1() float64

	// Lat2 is the latitude of point 2 in degrees
	Lat2() float64

	// Lon2 is the longitude of point 2 in degrees
	Lon2() float64

	// Azi2 is the azimuth at point 2 in degrees
	Azi2() float64

	// S12 is the distance between point 1 and point 2 in meters
	S12() float64

	// A12 is the arc length on the auxiliary sphere between point 1 and point 2 in degrees
	A12() float64

	// M12Reduced is the reduced length of the geodesic in meters
	M12Reduced() float64

	// M12 is the geodesic scale of point 2 relative to point 1 (dimensionless)
	M12() float64

	// M21 is the geodesic scale of point 1 relative to point 2 (dimensionless)
	M21() float64

	// S12Area is the area under the geodesic (meters²)
	S12Area() float64
}

type geodesicDataImpl struct {
	lat1       float64
	lon1       float64
	azi1       float64
	lat2       float64
	lon2       float64
	azi2       float64
	s12        float64
	a12        float64
	m12Reduced float64
	m12        float64
	m21        float64
	s12Area    float64
}

func newGeodesicDataImpl() *geodesicDataImpl {
	return &geodesicDataImpl{
		lat1:       math.NaN(),
		lon1:       math.NaN(),
		azi1:       math.NaN(),
		lat2:       math.NaN(),
		lon2:       math.NaN(),
		azi2:       math.NaN(),
		s12:        math.NaN(),
		a12:        math.NaN(),
		m12Reduced: math.NaN(),
		m12:        math.NaN(),
		m21:        math.NaN(),
		s12Area:    math.NaN(),
	}
}

func (g *geodesicDataImpl) Lat1() float64 {
	return g.lat1
}

func (g *geodesicDataImpl) Lon1() float64 {
	return g.lon1
}

func (g *geodesicDataImpl) Azi1() float64 {
	return g.azi1
}

func (g *geodesicDataImpl) Lat2() float64 {
	return g.lat2
}

func (g *geodesicDataImpl) Lon2() float64 {
	return g.lon2
}

func (g *geodesicDataImpl) Azi2() float64 {
	return g.azi2
}

func (g *geodesicDataImpl) S12() float64 {
	return g.s12
}

func (g *geodesicDataImpl) A12() float64 {
	return g.a12
}

func (g *geodesicDataImpl) M12Reduced() float64 {
	return g.m12Reduced
}

func (g *geodesicDataImpl) M12() float64 {
	return g.m12
}

func (g *geodesicDataImpl) M21() float64 {
	return g.m21
}

func (g *geodesicDataImpl) S12Area() float64 {
	return g.s12Area
}
