package capabilities

// BitMask represents a set of geodesic calculations to perform as an integer bitmask.
//
// When used as an argument to Geodesic.DirectWithCapabilities and Geodesic.InverseWithCapabilities,
// BitMask specifies which results to return in the GeodesicData struct.
//
// When used as an argument to NewGeodesicLineWithCapabilities and Geodesic.LineWithCapabilities,
// BitMask specifies which capabilities should be included in the GeodesicLine struct.
type BitMask int

const (
	none    BitMask = 0
	C1      BitMask = 1 << 0
	C1p     BitMask = 1 << 1
	C2      BitMask = 1 << 2
	C3      BitMask = 1 << 3
	C4      BitMask = 1 << 4
	all     BitMask = 0x1F
	mask            = all
	outAll  BitMask = 0x7F80
	OutMask BitMask = 0xFF80 // Include LongUnroll

	// None specifies: no capabilities, no output.
	None BitMask = 0

	// Latitude specifies: calculate latitude lat2. (It's not necessary to include this as a
	// capability to GeodesicLine because this is included by default.)
	Latitude = 1<<7 | none

	// Longitude specifies: calculate longitude lon2.
	Longitude = 1<<8 | C3

	// Azimuth specifies: calculate azimuths azi1 and azi2. (It's not necessary to include this as a
	// capability to GeodesicLine because this is included by default.)
	Azimuth = 1<<9 | none

	// Distance specifies: calculate distance s12.
	Distance = 1<<10 | C1

	// Standard specifies the default output and capabilities (latitudes, longitudes, azimuths, and
	// distance)
	Standard = Latitude | Longitude | Azimuth | Distance

	// DistanceIn specifies: allow distance s12 to be used as input in the direct geodesic problem.
	DistanceIn = 1<<11 | C1 | C1p

	// ReducedLength specifies: calculate reduced length m12.
	ReducedLength = 1<<12 | C1 | C2

	// GeodesicScale specifies: calculate geodesic scales M12 and M21.
	GeodesicScale = 1<<13 | C1 | C2

	// Area specifies: calculate area S12.
	Area = 1<<14 | C4

	// All capabilities, calculate everything. (LongUnroll is not included in this mask.)
	All = outAll | all

	// LongUnroll specifies: unroll lon2.
	LongUnroll BitMask = 1 << 15
)
