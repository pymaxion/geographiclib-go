package capabilities

// Mask represents a set of geodesic calculations to perform as an integer
// bitmask.
//
// When used as an argument to Geodesic.DirectWithCapabilities and
// Geodesic.InverseWithCapabilities, Mask specifies which results to return in
// the GeodesicData struct.
//
// When used as an argument to NewGeodesicLineWithCapabilities and
// Geodesic.LineWithCapabilities, Mask specifies which capabilities should be
// included in the GeodesicLine struct.
type Mask int

const (
	none    Mask = 0
	C1      Mask = 1 << 0
	C1p     Mask = 1 << 1
	C2      Mask = 1 << 2
	C3      Mask = 1 << 3
	C4      Mask = 1 << 4
	all     Mask = 0x1F
	mask         = all
	outAll  Mask = 0x7F80
	OutMask Mask = 0xFF80 // Include LongUnroll

	// None specifies: no capabilities, no output.
	None Mask = 0

	// Latitude specifies: calculate latitude lat2. (It's not necessary to include
	// this as a capability to GeodesicLine because this is included by default.)
	Latitude = 1<<7 | none

	// Longitude specifies: calculate longitude lon2.
	Longitude = 1<<8 | C3

	// Azimuth specifies: calculate azimuths azi1 and azi2. (It's not necessary to
	// include this as a capability to GeodesicLine because this is included by
	// default.)
	Azimuth = 1<<9 | none

	// Distance specifies: calculate distance s12.
	Distance = 1<<10 | C1

	// Standard specifies the default output and capabilities (latitudes, longitudes,
	// azimuths, and distance)
	Standard = Latitude | Longitude | Azimuth | Distance

	// DistanceIn specifies: allow distance s12 to be used as input in the direct
	// geodesic problem.
	DistanceIn = 1<<11 | C1 | C1p

	// ReducedLength specifies: calculate reduced length m12.
	ReducedLength = 1<<12 | C1 | C2

	// GeodesicScale specifies: calculate geodesic scales M12 and M21.
	GeodesicScale = 1<<13 | C1 | C2

	// Area specifies: calculate area S12.
	Area = 1<<14 | C4

	// All capabilities, calculate everything. (LongUnroll is not included in this
	// caps.)
	All = outAll | all

	// LongUnroll specifies: unroll lon2.
	LongUnroll Mask = 1 << 15
)
