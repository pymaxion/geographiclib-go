package geodesic

import (
	"math"
)

const (
	// digits represents the number of binary digits in the fraction of a double
	// precision number. This is equivalent to C++'s numeric_limits<double>::digits.
	digits        = 53
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
	epsilon = math.Nextafter(1., 2.) - 1. // https://stackoverflow.com/a/22185792/4755732
	// tiny is an underflow guard. We require tiny * epsilon > 0 and tiny + epsilon
	// == epsilon. Note that we are using 2^-1022 here instead of
	// math.SmallestNonzeroFloat64 to maintain consistency with other geographiclib
	// implementations.
	tiny    = math.Sqrt(math.Pow(2, -1022))
	tol0    = epsilon
	tol1    = 200 * tol0
	tol2    = math.Sqrt(tol0)
	tolb    = tol0 * tol2 // Check on bisection interval
	xthresh = 1000 * tol2
)

// sq squares a number (and avoids the overhead of math.Pow)
func sq(x float64) float64 {
	return x * x
}

// atanh calculates the inverse hyperbolic tangent of x. This is defined in terms
// of log1p(x) in order to maintain accuracy near x = 0. In addition, the odd
// parity of the function is enforced.
func atanh(x float64) float64 {
	y := math.Abs(x) // enforce odd parity
	y = math.Log1p(2*y/(1-y)) / 2
	if x > 0 {
		return y
	} else if x < 0 {
		return -y
	} else {
		return x
	}
}

// norm normalizes a sine/cosine pair (i.e., make sin²(x) + cos²(x) = 1)
func norm(sinx, cosx float64) (float64, float64) {
	// math.Hypot is inaccurate for 3.[89]. Problem reported by agdhruv
	// https://github.com/geopy/geopy/issues/466; see
	// https://bugs.python.org/issue43088
	// Visual Studio 2015 32-bit has a similar problem.
	r := math.Sqrt(sq(sinx) + sq(cosx))
	return sinx / r, cosx / r
}

// sum returns the error-free sum (s, t) of two numbers (u, v) where s =
// round(u+v) and t = u+v-s. See D. E. Knuth, TAOCP, Vol 2, 4.2.2, Theorem B.
func sum(u, v float64) (float64, float64) {
	s := u + v
	up := s - v
	vpp := s - up
	up -= u
	vpp -= v
	t := s
	if s != 0 {
		t = 0 - (up + vpp)
	}
	return s, t
}

// polyval evaluates a polynomial using Horner's method, where:
//  n = order of the polynomial
//  p = coefficient slice (of size n + s + 1 or more)
//  s = starting index for the slice
//  x = the variable
func polyval(n int, p []float64, s int, x float64) float64 {
	var y float64
	if n < 0 {
		y = 0
	} else {
		y = p[s]
	}
	for n > 0 {
		n--
		s++
		y = y*x + p[s]
	}
	return y
}

// angRound coarsens a value close to zero. The makes the smallest gap in x =
// 1/16 - nextafter(1/16, 0) = 1/2^57 for reals = 0.7 pm on the earth if x is an
// angle in degrees. (This is about 1000 times more resolution than we get with
// angles around 90 degrees.) We use this to avoid having to deal with near
// singular cases when x is non-zero but tiny (e.g., 1.0e-200). Note that tiny
// negative numbers get converted to -0.
func angRound(x float64) float64 {
	if x == 0 {
		return 0
	}
	const z = 1 / 16.
	y := math.Abs(x)
	// The compiler mustn't "simplify" z - (z - y) to y
	if y < z {
		y = z - (z - y)
	}
	return math.Copysign(y, x)
}

// angNormalize normalizes an angle in degrees to the range [-180, 180)
func angNormalize(x float64) float64 {
	y := math.Remainder(x, 360)
	if math.Abs(y) == 180 {
		return math.Copysign(180, x)
	}
	return y
}

// latFix replaces latitudes outside the range [-90, 90] by math.NaN
func latFix(x float64) float64 {
	if math.Abs(x) > 90 {
		return math.NaN()
	} else {
		return x
	}
}

// angDiff calculates the exact difference of two angles in degrees, reduced to
// [-180, 180]. More specifically, this function computes z = y - x exactly,
// reduced to [-180, 180], and then sets z = d + e where d is the nearest
// representable number to z and e is the truncation error. If z = ±0° or ±180°,
// then the sign of d is given by the sign of y - x. The maximum absolute value
// of e is 2^-26.
func angDiff(x, y float64) (float64, float64) {
	t1, t2 := sum(math.Remainder(-x, 360), math.Remainder(y, 360))
	d, e := sum(math.Remainder(t1, 360), t2)
	if d == 0 || math.Abs(d) == 180 {
		// If e == 0, take sign from y - x
		// else (e != 0, implies d = +/-180), d and e must have opposite signs
		d = math.Copysign(d, ternary(e == 0, y-x, -e))
	}
	return d, e
}

// deg2rad converts degrees to radians
func deg2rad(d float64) float64 {
	return d * math.Pi / 180.0
}

// rad2deg converts radians to degrees
func rad2deg(r float64) float64 {
	return r * 180.0 / math.Pi
}

// sincosd calculates the sine and cosine of x in degrees. The results obey
// exactly the elementary properties of the trigonometric functions, e.g., sin 9
// = cos 81 = - sin 123456789.
func sincosd(x float64) (float64, float64) {
	// In order to minimize round-off errors, this function exactly reduces the
	// argument to the range [-45, 45] before converting it to radians.
	r := math.Mod(x, 360)
	var q int
	if math.IsNaN(r) {
		q = 0
	} else {
		q = int(math.Round(r / 90))
	}
	r -= 90 * float64(q)
	r = deg2rad(r)
	s, c := math.Sincos(r)

	var sinx, cosx float64
	switch q & 3 {
	case 0:
		sinx, cosx = s, c
		break
	case 1:
		sinx, cosx = c, -s
		break
	case 2:
		sinx, cosx = -s, -c
		break
	default: // case 3
		sinx, cosx = -c, s
	}

	if sinx == 0 {
		sinx = math.Copysign(sinx, x)
	}
	return sinx, 0 + cosx
}

// sincosde evaluates the sine and cosine function in degrees with reduced
// arguments, plus a correction specified in degrees. This is a variant of
// sincosd allowing a correction to the angle to be supplied. x must be in
// [-180°, 180°] and t is assumed to be a *small* correction. angRound is applied
// to the reduced angle to prevent problems with x + t being extremely close, but
// not exactly equal, to one of the four cardinal directions.
func sincosde(x, t float64) (float64, float64) {
	// In order to minimize round-off errors, this function exactly reduces the
	// argument to the range [-45, 45] before converting it to radians.
	var q int
	if math.IsNaN(x) {
		q = 0
	} else {
		q = int(math.Round(x / 90))
	}
	r := x - 90*float64(q)

	// now abs(r) <= 45
	r = deg2rad(angRound(r + t))
	s, c := math.Sincos(r)

	var sinx, cosx float64
	switch q & 3 {
	case 0:
		sinx, cosx = s, c
		break
	case 1:
		sinx, cosx = c, -s
		break
	case 2:
		sinx, cosx = -s, -c
		break
	default: // case 3
		sinx, cosx = -c, s
	}

	if sinx == 0 {
		sinx = math.Copysign(sinx, x)
	}
	return sinx, 0 + cosx

}

// atan2d calculates atan2(y, x) with the result in degrees where y = the sine of
// the angle and x = the cosine of the angle. The result is in the range [-180
// 180]. N.B., atan2d(±0, -1) = +180.
func atan2d(y, x float64) float64 {
	// In order to minimize round-off errors, this function rearranges the arguments
	// so that result of atan2 is in the range [-pi/4, pi/4] before converting it to
	// degrees and mapping the result to the correct quadrant.
	q := 0
	if math.Abs(y) > math.Abs(x) {
		q = 2
		x, y = y, x
	}
	if x < 0 {
		x = -x
		q++
	}
	ang := rad2deg(math.Atan2(y, x))
	switch q {
	// Note that atan2d(-0.0, 1.0) will return -0.  However, we expect that
	// atan2d will not be called with y = -0.  If need be, include
	//
	//   case 0: ang = 0 + ang; break;
	//
	// and handle mpfr as in AngRound.
	case 1:
		ang = math.Copysign(180, y) - ang
		break
	case 2:
		ang = 90 - ang
		break
	case 3:
		ang = -90 + ang
		break
	default:
		break
	}
	return ang
}

// isfinite tests if x is finite
func isfinite(x float64) bool {
	return !math.IsNaN(x) && !math.IsInf(x, 0)
}

// min returns the minimum of two ints
func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

// ternary is a simple, naive implementation of a ternary operator
func ternary(predicate bool, valIfTrue, valIfFalse float64) float64 {
	if predicate {
		return valIfTrue
	}
	return valIfFalse
}

// sinCosSeries evaluates a trig series using Clenshaw summation
func sinCosSeries(sinp bool, sinx, cosx float64, c []float64) float64 {
	// Evaluate:
	//  y = sinp ? sum(c[i] * sin( 2*i    * x), i, 1, n) :
	//             sum(c[i] * cos((2*i+1) * x), i, 0, n-1)
	// using Clenshaw summation. N.B. c[0] is unused for sin series Approx operation
	// count = (n + 5) mult and (2 * n + 2) add
	k := len(c) // Point to one beyond last element
	n := k - int(ternary(sinp, 1, 0))
	ar := 2 * (cosx - sinx) * (cosx + sinx) // 2 * cos(2 * x)
	y1 := 0.                                // accumulators for sum
	y0 := 0.
	if n&1 != 0 {
		k -= 1
		y0 = c[k]
	}
	// Now n is even
	n /= 2
	for n != 0 { // while n--:
		n -= 1
		// Unroll loop x 2, so accumulators return to their original role
		k -= 1
		y1 = ar*y0 - y1 + c[k]
		k -= 1
		y0 = ar*y1 - y0 + c[k]
	}

	if sinp {
		return 2 * sinx * cosx * y0 // sin(2 * x) * y0
	} else {
		return cosx * (y0 - y1) // cos(x) * (y0 - y1)
	}
}

func initA3x(n float64) []float64 {
	a3x := make([]float64, nA3x)
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
	c3x := make([]float64, nC3x)
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
	c4x := make([]float64, nC4x)
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

// a1m1f calculates the scale factor A1-1 = mean value of (d/dsigma)I1 - 1
func a1m1f(eps float64) float64 {
	coeff := []float64{
		// (1-eps)*A1-1, polynomial in eps2 of order 3
		1, 4, 64, 0, 256,
	}
	m := nA1 / 2
	t := polyval(m, coeff, 0, sq(eps)) / coeff[m+1]
	return (t + eps) / (1 - eps)
}

// c1f computes the coefficients C1[l] in the Fourier expansion of B1
func c1f(eps float64, c []float64) {
	coeff := []float64{
		// C1[1]/eps^1, polynomial in eps2 of order 2
		-1, 6, -16, 32,
		// C1[2]/eps^2, polynomial in eps2 of order 2
		-9, 64, -128, 2048,
		// C1[3]/eps^3, polynomial in eps2 of order 1
		9, -16, 768,
		// C1[4]/eps^4, polynomial in eps2 of order 1
		3, -5, 512,
		// C1[5]/eps^5, polynomial in eps2 of order 0
		-7, 1280,
		// C1[6]/eps^6, polynomial in eps2 of order 0
		-7, 2048,
	}
	eps2 := sq(eps)
	d := eps
	o := 0
	for l := 1; l <= nC1; l++ { // l is index of C1p[l]
		m := (nC1 - l) / 2 // order of polynomial in eps^2
		c[l] = d * polyval(m, coeff, o, eps2) / coeff[o+m+1]
		o += m + 2
		d *= eps
	}
}

// c1pf computes the coefficients C1p[l] in the Fourier expansion of B1p
func c1pf(eps float64, c []float64) {
	coeff := []float64{
		// C1p[1]/eps^1, polynomial in eps2 of order 2
		205, -432, 768, 1536,
		// C1p[2]/eps^2, polynomial in eps2 of order 2
		4005, -4736, 3840, 12288,
		// C1p[3]/eps^3, polynomial in eps2 of order 1
		-225, 116, 384,
		// C1p[4]/eps^4, polynomial in eps2 of order 1
		-7173, 2695, 7680,
		// C1p[5]/eps^5, polynomial in eps2 of order 0
		3467, 7680,
		// C1p[6]/eps^6, polynomial in eps2 of order 0
		38081, 61440,
	}
	eps2 := sq(eps)
	d := eps
	o := 0
	for l := 1; l <= nC1p; l++ { // l is index of C1p[l]
		m := (nC1p - l) / 2 // order of polynomial in eps^2
		c[l] = d * polyval(m, coeff, o, eps2) / coeff[o+m+1]
		o += m + 2
		d *= eps
	}
}

// a2m1f calculates the scale factor A2-1 = mean value of (d/dsigma)I2 - 1
func a2m1f(eps float64) float64 {
	coeff := []float64{
		// (eps+1)*A2-1, polynomial in eps2 of order 3
		-11, -28, -192, 0, 256,
	}
	m := nA2 / 2
	t := polyval(m, coeff, 0, sq(eps)) / coeff[m+1]
	return (t - eps) / (1 + eps)
}

// c2f calculates the coefficients C2[l] in the Fourier expansion of B2
func c2f(eps float64, c []float64) {
	coeff := []float64{
		// C2[1]/eps^1, polynomial in eps2 of order 2
		1, 2, 16, 32,
		// C2[2]/eps^2, polynomial in eps2 of order 2
		35, 64, 384, 2048,
		// C2[3]/eps^3, polynomial in eps2 of order 1
		15, 80, 768,
		// C2[4]/eps^4, polynomial in eps2 of order 1
		7, 35, 512,
		// C2[5]/eps^5, polynomial in eps2 of order 0
		63, 1280,
		// C2[6]/eps^6, polynomial in eps2 of order 0
		77, 2048,
	}

	eps2 := sq(eps)
	d := eps
	o := 0
	for l := 1; l <= nC2; l++ { // l is index of C2[l]
		m := (nC2 - l) / 2 // order of polynomial in eps^2
		c[l] = d * polyval(m, coeff, o, eps2) / coeff[o+m+1]
		o += m + 2
		d *= eps
	}
}
