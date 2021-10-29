package geographiclib

// This file defines bit masks for what geodesic calculations to do.
//
// These masks do double duty. They specify (via the outmask parameter) which results to return in
// the GeodesicData struct returned by the general functions Geodesic.Direct(double, double, double,
// double, int) and Geodesic.Inverse(double, double, double, double, int) routines. They also
// signify (via the caps parameter) to the NewGeodesicLine(Geodesic, double, double, double, int)
// function and to Geodesic.Line(double, double, double, int) what capabilities should be included
// in the GeodesicLine struct.
const (
	capNone = 0
	capC1   = 1 << 0
	capC1p  = 1 << 1
	capC2   = 1 << 2
	capC3   = 1 << 3
	capC4   = 1 << 4
	capAll  = 0x1F
	capMask = capAll
	outAll  = 0x7F80
	outMask = 0xFF80 // Include LongUnroll

	// None specifies: no capabilities, no output.
	None = 0

	// Latitude specifies: calculate latitude lat2. (It's not necessary to include this as a
	// capability to GeodesicLine because this is included by default.)
	Latitude = 1<<7 | capNone

	// Longitude specifies: calculate longitude lon2.
	Longitude = 1<<8 | capC3

	// Azimuth specifies: calculate azimuths azi1 and azi2. (It's not necessary to include this as a
	// capability to GeodesicLine because this is included by default.)
	Azimuth = 1<<9 | capNone

	// Distance specifies: calculate distance s12.
	Distance = 1<<10 | capC1

	// All of the above, the "standard" output and capabilities.
	Standard = Latitude | Longitude | Azimuth | Distance

	// DistanceIn specifies: allow distance s12 to be used as input in the direct geodesic problem.
	DistanceIn = 1<<11 | capC1 | capC1p

	// ReducedLength specifies: calculate reduced length m12.
	ReducedLength = 1<<12 | capC1 | capC2

	// GeodesicScale specifies: calculate geodesic scales M12 and M21.
	GeodesicScale = 1<<13 | capC1 | capC2

	// Area specifies: calculate area S12.
	Area = 1<<14 | capC4

	// All capabilities, calculate everything. (LongUnroll is not included in this mask.)
	All = outAll | capAll

	// LongUnroll specifies: unroll lon2.
	LongUnroll = 1 << 15
)
