package readme

import (
	"fmt"
	"geographiclib-go/geodesic"
	"geographiclib-go/geodesic/capabilities"
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestWellingtonToSalamanca(t *testing.T) {
	r := geodesic.WGS84.Inverse(-41.32, 174.81, 40.96, -5.50)
	fmt.Printf("The distance is %.3f m.\n", r.S12)

	assert.InDelta(t, 19959679.267, r.S12, 1e-3)
}

func TestPointSouthwestOfPerth(t *testing.T) {
	r := geodesic.WGS84.Direct(-32.06, 115.74, 225, 20000e3)
	fmt.Printf("The position is (%.8f, %.8f).\n", r.Lat2, r.Lon2)

	assert.InDelta(t, 32.11195529, r.Lat2, 1e-8)
	assert.InDelta(t, -63.95925278, r.Lon2, 1e-8)
}

func TestJFKToLHR(t *testing.T) {
	r := geodesic.WGS84.InverseWithCapabilities(40.6, -73.8, 51.6, -0.5, capabilities.Area)
	fmt.Printf("The area is %.1f m^2.\n", r.S12Area)

	assert.InDelta(t, 40041368848742.5, r.S12Area, 1e-1)
}

func TestWaypointsBeijingToSanFrancisco(t *testing.T) {
	l := geodesic.WGS84.InverseLine(40.1, 116.6, 37.6, -122.4)
	ds := 1000e3
	n := int(math.Ceil(l.Distance() / ds))
	rs := make([]string, 0, n)
	for i := 0; i < n+1; i++ {
		if i == 0 {
			fmt.Println("distance latitude longitude azimuth")
		}
		s := math.Min(ds*float64(i), l.Distance())
		r := l.PositionWithCapabilities(s, capabilities.Standard|capabilities.LongUnroll)
		fmt.Printf("%.0f %.5f %.5f %.5f\n", r.S12, r.Lat2, r.Lon2, r.Azi2)
		rs = append(rs, fmt.Sprintf("%.0f %.5f %.5f %.5f", r.S12, r.Lat2, r.Lon2, r.Azi2))
	}

	assert.Len(t, rs, 11)
	assert.Equal(t, "0 40.10000 116.60000 42.91642", rs[0])
	assert.Equal(t, "1000000 46.37321 125.44903 48.99365", rs[1])
	assert.Equal(t, "2000000 51.78786 136.40751 57.29433", rs[2])
	assert.Equal(t, "3000000 55.92437 149.93825 68.24573", rs[3])
	assert.Equal(t, "4000000 58.27452 165.90776 81.68242", rs[4])
	assert.Equal(t, "5000000 58.43499 183.03167 96.29014", rs[5])
	assert.Equal(t, "6000000 56.37430 199.26948 109.99924", rs[6])
	assert.Equal(t, "7000000 52.45769 213.17327 121.33210", rs[7])
	assert.Equal(t, "8000000 47.19436 224.47209 129.98619", rs[8])
	assert.Equal(t, "9000000 41.02145 233.58294 136.34359", rs[9])
	assert.Equal(t, "9513998 37.60000 237.60000 138.89027", rs[10])
}

func TestWaypointsBeijingToSanFranciscoExpressedInArcLen(t *testing.T) {
	l := geodesic.WGS84.InverseLine(40.1, 116.6, 37.6, -122.4)
	n := int(math.Ceil(l.Arc()))
	da := l.Arc() / float64(n)
	rs := make([]string, 0, n)
	for i := 0; i < n+1; i++ {
		if i == 0 {
			fmt.Println("latitude longitude")
		}
		a := da * float64(i)
		r := l.ArcPositionWithCapabilities(a, capabilities.Latitude|capabilities.Longitude|capabilities.LongUnroll)
		fmt.Printf("%.5f %.5f\n", r.Lat2, r.Lon2)
		rs = append(rs, fmt.Sprintf("%.5f %.5f", r.Lat2, r.Lon2))
	}

	assert.Len(t, rs, 87)
	assert.Equal(t, "40.10000 116.60000", rs[0])
	assert.Equal(t, "40.82573 117.49243", rs[1])
	assert.Equal(t, "41.54435 118.40447", rs[2])
	assert.Equal(t, "42.25551 119.33686", rs[3])
	assert.Equal(t, "42.95886 120.29036", rs[4])
	assert.Equal(t, "43.65403 121.26575", rs[5])
	assert.Equal(t, "44.34062 122.26380", rs[6])
	// ...
	assert.Equal(t, "39.82385 235.05331", rs[83])
	assert.Equal(t, "39.08884 235.91990", rs[84])
	assert.Equal(t, "38.34746 236.76857", rs[85])
	assert.Equal(t, "37.60000 237.60000", rs[86])
}

func TestAreaOfAntarctica(t *testing.T) {
	p := geodesic.NewPolygonArea(geodesic.WGS84, false)
	antarctica := [][]float64{
		{-63.1, -58}, {-72.9, -74}, {-71.9, -102}, {-74.9, -102}, {-74.3, -131},
		{-77.5, -163}, {-77.4, 163}, {-71.7, 172}, {-65.9, 140}, {-65.7, 113},
		{-66.6, 88}, {-66.9, 59}, {-69.8, 25}, {-70.0, -4}, {-71.0, -14},
		{-77.3, -33}, {-77.9, -46}, {-74.7, -61},
	}

	for _, pnt := range antarctica {
		p.AddPoint(pnt[0], pnt[1])
	}
	r := p.Compute(false, true)
	fmt.Printf("Perimeter/area of Antarctica are %.3f m / %.1f m^2.\n", r.Perimeter, r.Area)

	assert.InDelta(t, 16831067.893, r.Perimeter, 1e-3)
	assert.InDelta(t, 13662703680020.1, r.Area, 1e-1)
}
