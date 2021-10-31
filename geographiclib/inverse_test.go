package geographiclib

import (
	"fmt"
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestInverse(t *testing.T) {
	for i, tt := range testCases {
		t.Run(fmt.Sprintf("test case #%d", i), func(t *testing.T) {
			r := WGS84.Inverse(tt.lat1, tt.lon1, tt.lat2, tt.lon2, All|LongUnroll)

			assert.InDelta(t, tt.lat1, r.Lat1(), 1e-13)
			assert.InDelta(t, tt.lon1, r.Lon1(), 1e-13)
			assert.InDelta(t, tt.lat2, r.Lat2(), 1e-13)
			assert.InDelta(t, tt.lon2, r.Lon2(), 1e-13)
			assert.InDelta(t, tt.azi1, r.Azi1(), 1e-13)
			assert.InDelta(t, tt.azi2, r.Azi2(), 1e-13)
			assert.InDelta(t, tt.s12, r.S12(), 1e-8)
			assert.InDelta(t, tt.a12, r.A12(), 1e-13)
			assert.InDelta(t, tt.m12, r.M12Reduced(), 1e-8)
			assert.InDelta(t, tt.M12, r.M12(), 1e-15)
			assert.InDelta(t, tt.M21, r.M21(), 1e-15)
			assert.InDelta(t, tt.S12, r.S12Area(), 0.1)
		})
	}
}
