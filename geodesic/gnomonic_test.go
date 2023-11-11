package geodesic

import (
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestGnomonicForward(t *testing.T) {
	t.Run("forward calc", func(t *testing.T) {
		lat0, lon0 := 48+50/60.0, 2+20/60.0 // Paris
		lat, lon := 50.9, 1.8               // Calais
		g := Gnomonic{Earth: WGS84}
		r := g.Forward(lat0, lon0, lat, lon)
		assert.InDelta(t, -37543.7, r.X, 0.05)
		assert.InDelta(t, 230103, r.Y, 0.25)
	})
}

func TestGnomonicReverse(t *testing.T) {
	t.Run("reverse calc", func(t *testing.T) {
		lat0, lon0 := 48+50/60.0, 2+20/60.0 // Paris
		x, y := -38e3, 230e3                // Calais
		g := Gnomonic{Earth: WGS84}
		r := g.Reverse(lat0, lon0, x, y)
		assert.InDelta(t, 50.899, r.Lat, 0.0005)
		assert.InDelta(t, 1.79353, r.Lon, 0.000005)
	})
}
