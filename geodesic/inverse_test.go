package geodesic

import (
	"fmt"
	"testing"

	caps "geographiclib-go/geodesic/capabilities"

	"github.com/stretchr/testify/assert"
)

func TestInverse_CommonCases(t *testing.T) {
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

func TestInverse_GeodSolveCases(t *testing.T) {
	testCases := []geodSolve{
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
	}
	for _, testCase := range testCases {
		t.Run(testCase.String(), testCase.logic)
	}
}
