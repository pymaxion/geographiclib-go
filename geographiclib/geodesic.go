package geographiclib

import (
	"errors"
	"math"
)

type Geodesic interface {
	Inverse(lat1, lon1, lat2, lon2 float64, outmask int) GeodesicData
}

func NewGeodesic(a, f float64) (Geodesic, error) {
	return newGeodesicImpl(a, f)
}

type geodesicImpl struct {
	a     float64
	f     float64
	f1    float64
	e2    float64
	ep2   float64
	b     float64
	c2    float64
	n     float64
	etol2 float64
	a3x   []float64
	c3x   []float64
	c4x   []float64
}

func newGeodesicImpl(a, f float64) (*geodesicImpl, error) {
	if !(isfinite(a) && a > 0) {
		return nil, errors.New("equatorial radius is not positive")
	}
	f1 := 1 - f
	b := a * f1
	if !(isfinite(b) && b > 0) {
		return nil, errors.New("polar semi-axis is not positive")
	}
	e2 := f * (2 - f)
	ep2 := e2 / sq(f1) // e2 / 1 - e2
	n := f / (2 - f)
	// authalic radius squared
	c2 := calculateAuthalicRadiusSquared(a, b, e2)
	// The sig12 threshold for "really short". Using the auxiliary sphere solution with dnm computed
	// at (bet1 + bet2) / 2, the relative error in the azimuth consistency check is sig12^2 * abs(f)
	// * min(1, 1-f/2) / 2. (Error measured for 1/100 < b/a < 100 and abs(f) >= 1/1000. For a given f
	// and sig12, the max error occurs for lines near the pole. If the old rule for computing dnm =
	// (dn1 + dn2)/2 is used, then the error increases by a factor of 2.) Setting this equal to
	// epsilon gives sig12 = etol2. Here 0.1 is a safety factor (error decreased by 100) and
	// max(0.001, abs(f)) stops etol2 getting too large in the nearly spherical case.
	etol2 := 0.1 * tol2 / math.Sqrt(math.Max(0.001, math.Abs(f))*math.Min(1, 1-f/2)/2)
	a3x := initA3x(n)
	c3x := initC3x(n)
	c4x := initC4x(n)
	return &geodesicImpl{
		a:     a,
		f:     f,
		f1:    f1,
		e2:    e2,
		ep2:   ep2,
		b:     b,
		c2:    c2,
		n:     n,
		etol2: etol2,
		a3x:   a3x,
		c3x:   c3x,
		c4x:   c4x,
	}, nil
}

func calculateAuthalicRadiusSquared(a, b, e2 float64) float64 {
	var multiplier float64
	if e2 == 0 {
		multiplier = 1
	} else {
		var dividend float64
		if e2 > 0 {
			dividend = atanh(math.Sqrt(e2))
		} else {
			dividend = math.Atan(math.Sqrt(-e2))
		}
		divisor := math.Sqrt(math.Abs(e2))
		multiplier = dividend / divisor
	}
	return (sq(a) + sq(b)*multiplier) / 2
}

// a3f evaluates A3
func (g *geodesicImpl) a3f(eps float64) float64 {
	return polyval(nA3-1, g.a3x, 0, eps)
}

// c3f evaluates C3 coefficients and sets elements c[1] thru c[nC3 - 1]
func (g *geodesicImpl) c3f(eps float64, c []float64) {
	mult := 1.
	o := 0
	for l := 1; l < nC3; l++ { // l is index of C3[l]
		m := nC3 - l - 1 // order of polynomial in eps
		mult *= eps
		c[l] = mult * polyval(m, g.c3x, o, eps)
		o += m + 1
	}
}

// c4f evaluates C4 coefficients and sets elements c[0] thru c[nC4 - 1]
func (g *geodesicImpl) c4f(eps float64, c []float64) {
	mult := 1.
	o := 0
	for l := 0; l < nC4; l++ { // l is index of C4[l]
		m := nC4 - l - 1 // order of polynomial in eps
		c[l] = mult * polyval(m, g.c4x, o, eps)
		o += m + 1
		mult *= eps
	}
}

func (g *geodesicImpl) Inverse(lat1, lon1, lat2, lon2 float64, outmask int) GeodesicData {
	solver := inverseSolver{g}
	return solver.inverse(lat1, lon1, lat2, lon2, outmask)
}
