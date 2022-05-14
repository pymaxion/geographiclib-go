package geodesic

import (
	"fmt"
	"testing"

	"github.com/pymaxion/geographiclib-go/geodesic/capabilities"
	"github.com/stretchr/testify/assert"
)

func TestInverse(t *testing.T) {
	for i, tt := range commonTestCases {
		t.Run(fmt.Sprintf("common test case #%d", i), func(t *testing.T) {
			r := WGS84.InverseWithCapabilities(tt.Lat1, tt.Lon1, tt.Lat2, tt.Lon2, capabilities.All|capabilities.LongUnroll)
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
		geodSolve0,
		geodSolve2,
		geodSolve4,
		geodSolve6,
		geodSolve9,
		geodSolve10,
		geodSolve11,
		geodSolve12,
		geodSolve14,
		geodSolve26,
		geodSolve29,
		geodSolve33,
		geodSolve55,
		geodSolve59,
		geodSolve74,
		geodSolve76,
		geodSolve78,
		geodSolve80,
		geodSolve92,
		geodSolve94,
		geodSolve96,
	} {
		t.Run(testCase.String(), testCase.logic)
	}

	t.Run("azimuth with coincident point on equator", func(t *testing.T) {
		testCases := []struct {
			lat1 float64
			lat2 float64
			azi  float64
		}{
			{0, minusZero, 180},
			{minusZero, 0, 0},
		}

		for _, tt := range testCases {
			t.Run(fmt.Sprintf("lat1: %.3f, lat2: %.3f, azi: %.3f", tt.lat1, tt.lat2, tt.azi), func(t *testing.T) {
				r := WGS84.Inverse(tt.lat1, 0, tt.lat2, 0)
				assert.True(t, equiv(tt.azi, r.Azi1))
				assert.True(t, equiv(tt.azi, r.Azi2))
			})
		}
	})

	t.Run("Does the nearly antipodal equatorial solution go north or south?", func(t *testing.T) {
		testCases := []struct {
			lat1 float64
			lat2 float64
			azi1 float64
			azi2 float64
		}{
			{0, 0, 56, 124},
			{minusZero, minusZero, 124, 56},
		}

		for _, tt := range testCases {
			t.Run(fmt.Sprintf("lat1: %.3f, lat2: %.3f", tt.lat1, tt.lat2), func(t *testing.T) {
				r := WGS84.Inverse(tt.lat1, 0, tt.lat2, 179.5)
				assert.InDelta(t, tt.azi1, r.Azi1, 1)
				assert.InDelta(t, tt.azi2, r.Azi2, 1)
			})
		}
	})

	t.Run("Does the exact antipodal equatorial path go N/S + E/W?", func(t *testing.T) {
		testCases := []struct {
			lat1 float64
			lat2 float64
			lon2 float64
			azi1 float64
			azi2 float64
		}{
			{0, 0, +180, 0, +180},
			{minusZero, minusZero, +180, +180, 0},
			{0, 0, -180, minusZero, -180},
			{minusZero, minusZero, -180, -180, minusZero},
		}

		for _, tt := range testCases {
			t.Run(fmt.Sprintf("lat1: %.3f, lat2: %.3f, lon2: %.3f", tt.lat1, tt.lat2, tt.lon2), func(t *testing.T) {
				r := WGS84.Inverse(tt.lat1, 0, tt.lat2, tt.lon2)
				assert.True(t, equiv(tt.azi1, r.Azi1))
				assert.True(t, equiv(tt.azi2, r.Azi2))
			})
		}
	})

	t.Run("Antipodal points on the equator with prolate ellipsoid", func(t *testing.T) {
		testCases := []struct {
			lon2 float64
			azi  float64
		}{
			{+180, +90},
			{-180, -90},
		}

		g, _ := NewGeodesic(6.4e6, -1/300.0)
		for _, tt := range testCases {
			t.Run(fmt.Sprintf("lon2: %.3f", tt.lon2), func(t *testing.T) {
				r := g.Inverse(0, 0, 0, tt.lon2)
				assert.True(t, equiv(tt.azi, r.Azi1))
				assert.True(t, equiv(tt.azi, r.Azi2))
			})
		}
	})
}
