package capabilities

import (
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestCapabilities(t *testing.T) {
	testCases := []struct {
		name          string
		mask          BitMask
		expectedValue int
	}{
		{"None", None, 0},
		{"Latitude", Latitude, 128},
		{"Longitude", Longitude, 264},
		{"Azimuth", Azimuth, 512},
		{"Distance", Distance, 1025},
		{"Standard", Standard, 1929},
		{"DistanceIn", DistanceIn, 2051},
		{"ReducedLength", ReducedLength, 4101},
		{"GeodesicScale", GeodesicScale, 8197},
		{"Area", Area, 16400},
		{"All", All, 32671},
		{"LongUnroll", LongUnroll, 32768},
	}

	for _, testCase := range testCases {
		t.Run(testCase.name, func(t *testing.T) {
			assert.Equal(t, testCase.expectedValue, int(testCase.mask))
		})
	}
}
