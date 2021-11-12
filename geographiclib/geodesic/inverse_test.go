package geodesic

import (
	"fmt"
	"math"
	"testing"

	caps "geographiclib-go/geographiclib/geodesic/capabilities"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

func TestInverse(t *testing.T) {
	for i, tt := range commonTestCases {
		t.Run(fmt.Sprintf("common test case #%d", i), func(t *testing.T) {
			r := WGS84.InverseWithCapabilities(tt.Lat1(), tt.Lon1(), tt.Lat2(), tt.Lon2(), caps.All|caps.LongUnroll)
			assert.InDelta(t, tt.Lat1(), r.Lat1(), 1e-13)
			assert.InDelta(t, tt.Lon1(), r.Lon1(), 1e-13)
			assert.InDelta(t, tt.Lat2(), r.Lat2(), 1e-13)
			assert.InDelta(t, tt.Lon2(), r.Lon2(), 1e-13)
			assert.InDelta(t, tt.Azi1(), r.Azi1(), 1e-13)
			assert.InDelta(t, tt.Azi2(), r.Azi2(), 1e-13)
			assert.InDelta(t, tt.S12(), r.S12(), 1e-8)
			assert.InDelta(t, tt.A12(), r.A12(), 1e-13)
			assert.InDelta(t, tt.M12Reduced(), r.M12Reduced(), 1e-8)
			assert.InDelta(t, tt.M12(), r.M12(), 1e-15)
			assert.InDelta(t, tt.M21(), r.M21(), 1e-15)
			assert.InDelta(t, tt.S12Area(), r.S12Area(), 0.1)
		})
	}
}

func TestGeodSolveCases(t *testing.T) {
	t.Run("GeodSolve0", func(t *testing.T) {
		r := WGS84.Inverse(40.6, -73.8, 49.01666667, 2.55)
		assert.InDelta(t, 53.47022, r.Azi1(), 0.5e-5)
		assert.InDelta(t, 111.59367, r.Azi2(), 0.5e-5)
		assert.InDelta(t, 5853226, r.S12(), 0.5)
	})

	// Check fix for antipodal prolate bug found 2010-09-04
	t.Run("GeodSolve2", func(t *testing.T) {
		geod, err := NewGeodesic(6.4e6, -1/150.0)
		require.Nil(t, err)

		r := geod.Inverse(0.07476, 0, -0.07476, 180)
		assert.InDelta(t, 90.00078, r.Azi1(), 0.5e-5)
		assert.InDelta(t, 90.00078, r.Azi2(), 0.5e-5)
		assert.InDelta(t, 20106193, r.S12(), 0.5)

		r = geod.Inverse(0.1, 0, -0.1, 180)
		assert.InDelta(t, 90.00105, r.Azi1(), 0.5e-5)
		assert.InDelta(t, 90.00105, r.Azi2(), 0.5e-5)
		assert.InDelta(t, 20106193, r.S12(), 0.5)
	})

	// Check fix for short line bug found 2010-05-21
	t.Run("GeodSolve4", func(t *testing.T) {
		r := WGS84.Inverse(36.493349428792, 0, 36.49334942879201, .0000008)
		assert.InDelta(t, 0.072, r.S12(), 0.5e-3)
	})

	// Check fix for volatile sbet12a bug found 2011-06-25 (gcc 4.4.4 x86 -O3). Found again on
	// 2012-03-27 with tdm-mingw32 (g++ 4.6.1).
	t.Run("GeodSolve6", func(t *testing.T) {
		r := WGS84.Inverse(88.202499451857, 0, -88.202499451857, 179.981022032992859592)
		assert.InDelta(t, 20003898.214, r.S12(), 0.5e-3)

		r = WGS84.Inverse(89.262080389218, 0, -89.262080389218, 179.992207982775375662)
		assert.InDelta(t, 20003925.854, r.S12(), 0.5e-3)

		r = WGS84.Inverse(89.333123580033, 0, -89.333123580032997687, 179.99295812360148422)
		assert.InDelta(t, 20003926.881, r.S12(), 0.5e-3)
	})

	// Check fix for volatile x bug found 2011-06-25 (gcc 4.4.4 x86 -O3)
	t.Run("GeodSolve9", func(t *testing.T) {
		r := WGS84.Inverse(56.320923501171, 0, -56.320923501171, 179.664747671772880215)
		assert.InDelta(t, 19993558.287, r.S12(), 0.5e-3)
	})

	// Check fix for adjust tol1_ bug found 2011-06-25 (Visual Studio 10 rel + debug)
	t.Run("GeodSolve10", func(t *testing.T) {
		r := WGS84.Inverse(52.784459512564, 0, -52.784459512563990912, 179.634407464943777557)
		assert.InDelta(t, 19991596.095, r.S12(), 0.5e-3)
	})

	// Check fix for bet2 = -bet1 bug found 2011-06-25 (Visual Studio 10 rel + debug)
	t.Run("GeodSolve11", func(t *testing.T) {
		r := WGS84.Inverse(48.522876735459, 0, -48.52287673545898293, 179.599720456223079643)
		assert.InDelta(t, 19989144.774, r.S12(), 0.5e-3)
	})

	// Check fix for inverse geodesics on extreme prolate/oblate ellipsoids Reported 2012-08-29
	// Stefan Guenther <stefan.gunther@embl.de>; fixed 2012-10-07
	t.Run("GeodSolve12", func(t *testing.T) {
		geod, err := NewGeodesic(89.8, -1.83)
		require.Nil(t, err)

		r := geod.Inverse(0, 0, -10, 160)
		assert.InDelta(t, 120.27, r.Azi1(), 1e-2)
		assert.InDelta(t, 105.15, r.Azi2(), 1e-2)
		assert.InDelta(t, 266.7, r.S12(), 1e-1)
	})

	// Check fix for inverse ignoring lon12 = nan
	t.Run("GeodSolve14", func(t *testing.T) {
		r := WGS84.Inverse(0, 0, 1, math.NaN())
		assert.True(t, math.IsNaN(r.Azi1()))
		assert.True(t, math.IsNaN(r.Azi2()))
		assert.True(t, math.IsNaN(r.S12()))
	})

	t.Run("GeodSolve26", func(t *testing.T) {
		geod, err := NewGeodesic(6.4e6, 0)
		require.Nil(t, err)

		r := geod.InverseWithCapabilities(1, 2, 3, 4, caps.Area)
		assert.InDelta(t, 49911046115.0, r.S12Area(), 0.5)
	})
}
