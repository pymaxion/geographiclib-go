package geographiclib

import (
	"errors"
	"math"
)

type Geodesic interface {
	F1() float64

	Inverse(lat1, lon1, lat2, lon2 float64, outmask int) GeodesicData
}

func NewGeodesic(a, f float64) (Geodesic, error) {
	return newGeodesicImpl(a, f)
}

const (
	geodesicOrder = 6
	nA1           = geodesicOrder
	nC1           = geodesicOrder
	nC1p          = geodesicOrder
	nA2           = geodesicOrder
	nC2           = geodesicOrder
	nA3           = geodesicOrder
	nA3x          = nA3
	nC3           = geodesicOrder
	nC3x          = (nC3 * (nC3 - 1)) / 2
	nC4           = geodesicOrder
	nC4x          = (nC4 * (nC4 + 1)) / 2
	maxit1        = 20
	maxit2        = maxit1 + digits + 10
)

var (
	// using 2^-1022 instead of math.SmallestNonzeroFloat64 to maintain consistency with other
	// geographiclib implementations
	tiny    = math.Sqrt(math.Pow(2, -1022))
	epsilon = math.Nextafter(1., 2.) - 1. // https://stackoverflow.com/a/22185792/4755732
	tol0    = epsilon
	tol1    = 200 * tol0
	tol2    = math.Sqrt(tol0)
	tolb    = tol0 * tol2
	xthresh = 1000 * tol2
)

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

func initA3x(n float64) []float64 {
	a3x := make([]float64, 0, nA3x)
	coeff := []float64{
		// A3, coeff of eps^5, polynomial in n of order 0
		-3, 128,
		// A3, coeff of eps^4, polynomial in n of order 1
		-2, -3, 64,
		// A3, coeff of eps^3, polynomial in n of order 2
		-1, -3, -1, 16,
		// A3, coeff of eps^2, polynomial in n of order 2
		3, -1, -2, 8,
		// A3, coeff of eps^1, polynomial in n of order 1
		1, -1, 2,
		// A3, coeff of eps^0, polynomial in n of order 0
		1, 1,
	}
	o, k := 0, 0
	for j := nA3 - 1; j >= 0; j-- { // coeff of eps^j
		m := min(nA3-j-1, j) // order of polynomial in n
		a3x[k] = polyval(m, coeff, o, n) / coeff[o+m+1]
		k++
		o += m + 2
	}
	return a3x
}

func initC3x(n float64) []float64 {
	c3x := make([]float64, 0, nC3x)
	coeff := []float64{
		// C3[1], coeff of eps^5, polynomial in n of order 0
		3, 128,
		// C3[1], coeff of eps^4, polynomial in n of order 1
		2, 5, 128,
		// C3[1], coeff of eps^3, polynomial in n of order 2
		-1, 3, 3, 64,
		// C3[1], coeff of eps^2, polynomial in n of order 2
		-1, 0, 1, 8,
		// C3[1], coeff of eps^1, polynomial in n of order 1
		-1, 1, 4,
		// C3[2], coeff of eps^5, polynomial in n of order 0
		5, 256,
		// C3[2], coeff of eps^4, polynomial in n of order 1
		1, 3, 128,
		// C3[2], coeff of eps^3, polynomial in n of order 2
		-3, -2, 3, 64,
		// C3[2], coeff of eps^2, polynomial in n of order 2
		1, -3, 2, 32,
		// C3[3], coeff of eps^5, polynomial in n of order 0
		7, 512,
		// C3[3], coeff of eps^4, polynomial in n of order 1
		-10, 9, 384,
		// C3[3], coeff of eps^3, polynomial in n of order 2
		5, -9, 5, 192,
		// C3[4], coeff of eps^5, polynomial in n of order 0
		7, 512,
		// C3[4], coeff of eps^4, polynomial in n of order 1
		-14, 7, 512,
		// C3[5], coeff of eps^5, polynomial in n of order 0
		21, 2560,
	}
	o, k := 0, 0
	for l := 1; l < nC3; l++ { // l is index of C3[l]
		for j := nC3 - 1; j >= l; j-- { // coeff of eps^j
			m := min(nC3-j-1, j) // order of polynomial in n
			c3x[k] = polyval(m, coeff, o, n) / coeff[o+m+1]
			k++
			o += m + 2
		}
	}
	return c3x
}

func initC4x(n float64) []float64 {
	c4x := make([]float64, 0, nC4x)
	coeff := []float64{
		// C4[0], coeff of eps^5, polynomial in n of order 0
		97, 15015,
		// C4[0], coeff of eps^4, polynomial in n of order 1
		1088, 156, 45045,
		// C4[0], coeff of eps^3, polynomial in n of order 2
		-224, -4784, 1573, 45045,
		// C4[0], coeff of eps^2, polynomial in n of order 3
		-10656, 14144, -4576, -858, 45045,
		// C4[0], coeff of eps^1, polynomial in n of order 4
		64, 624, -4576, 6864, -3003, 15015,
		// C4[0], coeff of eps^0, polynomial in n of order 5
		100, 208, 572, 3432, -12012, 30030, 45045,
		// C4[1], coeff of eps^5, polynomial in n of order 0
		1, 9009,
		// C4[1], coeff of eps^4, polynomial in n of order 1
		-2944, 468, 135135,
		// C4[1], coeff of eps^3, polynomial in n of order 2
		5792, 1040, -1287, 135135,
		// C4[1], coeff of eps^2, polynomial in n of order 3
		5952, -11648, 9152, -2574, 135135,
		// C4[1], coeff of eps^1, polynomial in n of order 4
		-64, -624, 4576, -6864, 3003, 135135,
		// C4[2], coeff of eps^5, polynomial in n of order 0
		8, 10725,
		// C4[2], coeff of eps^4, polynomial in n of order 1
		1856, -936, 225225,
		// C4[2], coeff of eps^3, polynomial in n of order 2
		-8448, 4992, -1144, 225225,
		// C4[2], coeff of eps^2, polynomial in n of order 3
		-1440, 4160, -4576, 1716, 225225,
		// C4[3], coeff of eps^5, polynomial in n of order 0
		-136, 63063,
		// C4[3], coeff of eps^4, polynomial in n of order 1
		1024, -208, 105105,
		// C4[3], coeff of eps^3, polynomial in n of order 2
		3584, -3328, 1144, 315315,
		// C4[4], coeff of eps^5, polynomial in n of order 0
		-128, 135135,
		// C4[4], coeff of eps^4, polynomial in n of order 1
		-2560, 832, 405405,
		// C4[5], coeff of eps^5, polynomial in n of order 0
		128, 99099,
	}
	o, k := 0, 0
	for l := 0; l < nC4; l++ { // l is index of C4[l]
		for j := nC4 - 1; j >= l; j-- { // coeff of eps^j
			m := nC4 - j - 1 // order of polynomial in n
			c4x[k] = polyval(m, coeff, o, n) / coeff[o+m+1]
			k++
			o += m + 2
		}
	}
	return c4x
}

func (g *geodesicImpl) F1() float64 {
	return g.f1
}

func (g *geodesicImpl) Inverse(lat1, lon1, lat2, lon2 float64, outmask int) GeodesicData {
	solver := inverseSolver{g}
	return solver.inverse(lat1, lon1, lat2, lon2, outmask)
}
