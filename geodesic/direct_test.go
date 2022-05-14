package geodesic

import (
	"fmt"
	"testing"

	"github.com/pymaxion/geographiclib-go/geodesic/capabilities"

	"github.com/stretchr/testify/assert"
)

func TestDirect(t *testing.T) {
	for i, tt := range commonTestCases {
		t.Run(fmt.Sprintf("common test case #%d", i), func(t *testing.T) {
			r := WGS84.DirectWithCapabilities(tt.Lat1, tt.Lon1, tt.Azi1, tt.S12, capabilities.All|capabilities.LongUnroll)
			assert.InDelta(t, tt.Lat1, r.Lat1, 1e-13)
			assert.InDelta(t, tt.Lon1, r.Lon1, 1e-13)
			assert.InDelta(t, tt.Lat2, r.Lat2, 1e-13)
			assert.InDelta(t, tt.Lon2, r.Lon2, 1e-13)
			assert.InDelta(t, tt.Azi1, r.Azi1, 1e-13)
			assert.InDelta(t, tt.Azi2, r.Azi2, 1e-13)
			assert.InDelta(t, tt.S12, r.S12, 1e-8)
			assert.InDelta(t, tt.A12, r.A12, 1e-13)
			assert.InDelta(t, tt.M12Reduced, r.M12Reduced, 1e-8)
			assert.InDelta(t, tt.M12, r.M12, 1e-15)
			assert.InDelta(t, tt.M21, r.M21, 1e-15)
			assert.InDelta(t, tt.S12Area, r.S12Area, 0.1)
		})
	}

	for _, testCase := range []geodSolve{
		geodSolve1,
		geodSolve5,
		geodSolve15,
		geodSolve17,
		geodSolve28,
		geodSolve61,
		geodSolve73,
		geodSolve84,
	} {
		t.Run(testCase.String(), testCase.logic)
	}

	t.Run("azimuths = +/-0 and +/-180 for the direct problem", func(t *testing.T) {
		testCases := []struct {
			azi1 float64
			lon2 float64
			azi2 float64
		}{
			{0, +180, +180},
			{minusZero, -180, -180},
			{+180, +180, 0},
			{-180, -180, minusZero},
		}

		for _, tt := range testCases {
			t.Run(fmt.Sprintf("azi1: %.3f", tt.azi1), func(t *testing.T) {
				r := WGS84.DirectWithCapabilities(0, 0, tt.azi1, 15e6, capabilities.Standard|capabilities.LongUnroll)
				assert.True(t, equiv(tt.lon2, r.Lon2))
				assert.True(t, equiv(tt.azi2, r.Azi2))
			})
		}
	})
}

func TestArcDirect(t *testing.T) {
	for i, tt := range commonTestCases {
		t.Run(fmt.Sprintf("common test case #%d", i), func(t *testing.T) {
			r := WGS84.ArcDirectWithCapabilities(tt.Lat1, tt.Lon1, tt.Azi1, tt.A12, capabilities.All|capabilities.LongUnroll)
			assert.InDelta(t, tt.Lat1, r.Lat1, 1e-13)
			assert.InDelta(t, tt.Lon1, r.Lon1, 1e-13)
			assert.InDelta(t, tt.Lat2, r.Lat2, 1e-13)
			assert.InDelta(t, tt.Lon2, r.Lon2, 1e-13)
			assert.InDelta(t, tt.Azi1, r.Azi1, 1e-13)
			assert.InDelta(t, tt.Azi2, r.Azi2, 1e-13)
			assert.InDelta(t, tt.S12, r.S12, 1e-8)
			assert.InDelta(t, tt.A12, r.A12, 1e-13)
			assert.InDelta(t, tt.M12Reduced, r.M12Reduced, 1e-8)
			assert.InDelta(t, tt.M12, r.M12, 1e-15)
			assert.InDelta(t, tt.M21, r.M21, 1e-15)
			assert.InDelta(t, tt.S12Area, r.S12Area, 0.1)
		})
	}
}
