package geographiclib

import "math"

const (
	// digits represents the number of binary digits in the fraction of a double precision number.
	// This is equivalent to C++'s numeric_limits<double>::digits.
	digits = 53
)

// sq squares a number (and avoids the overhead of math.Pow)
func sq(x float64) float64 {
	return x * x
}

// atanh calculates the inverse hyperbolic tangent of x. This is defined in terms of log1p(x) in
// order to maintain accuracy near x = 0. In addition, the odd parity of the function is enforced.
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
	r := math.Hypot(sinx, cosx)
	return sinx / r, cosx / r
}

// sum returns the error-free sum (s, t) of two numbers (u, v) where s = round(u+v) and t = u+v-s.
// See D. E. Knuth, TAOCP, Vol 2, 4.2.2, Theorem B.
func sum(u, v float64) (float64, float64) {
	s := u + v
	up := s - v
	vpp := s - up
	up -= u
	vpp -= v
	t := -(up + vpp)
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

// angRound coarsens a value close to zero. The makes the smallest gap in x = 1/16 - nextafter(1/16,
// 0) = 1/2^57 for reals = 0.7 pm on the earth if x is an angle in degrees. (This is about 1000
// times more resolution than we get with angles around 90 degrees.) We use this to avoid having to
// deal with near singular cases when x is non-zero but tiny (e.g., 1.0e-200). Note that tiny
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
	if x < 0 {
		return -y
	} else {
		return y
	}
}

// remainder calculates the remainder of x/y in the range [-y/2, y/2]. The range of x is
// unrestricted, but y must be positive.
func remainder(x, y float64) float64 {
	x = math.Mod(x, y)
	if x < -y/2 {
		return x + y
	} else if x < y/2 {
		return x
	} else {
		return x - y
	}
}

// angNormalize normalizes an angle in degrees to the range [-180, 180)
func angNormalize(x float64) float64 {
	x = remainder(x, 360)
	if x == -180 {
		return 180
	} else {
		return x
	}
}

// latFix replaces latitudes outside the range [-90, 90] by math.NaN
func latFix(x float64) float64 {
	if math.Abs(x) > 90 {
		return math.NaN()
	} else {
		return x
	}
}

// angDiff calculates the exact difference of two angles in degrees, reduced to (-180, 180]. More
// specifically, this function z = y - x exactly, reduced to (-180, 180], and then sets z = d + e
// where d is the nearest representable number to z and e is the truncation error. If d = -180, then
// e > 0; if d = 180, then e <= 0.
func angDiff(x, y float64) (float64, float64) {
	d, t := sum(angNormalize(-x), angNormalize(y))
	d = angNormalize(d)
	if d == 180 && t > 0 {
		return sum(-180, t)
	} else {
		return sum(d, t)
	}
}

// deg2rad converts degrees to radians
func deg2rad(d float64) float64 {
	return d * math.Pi / 180.0
}

// rad2deg converts radians to degrees
func rad2deg(r float64) float64 {
	return r * 180.0 / math.Pi
}

// sincosd calculates the sine and cosine of x in degrees. The results obey exactly the elementary
// properties of the trigonometric functions, e.g., sin 9 = cos 81 = - sin 123456789.
func sincosd(x float64) (float64, float64) {
	// In order to minimize round-off errors, this function exactly reduces the argument to the
	//  range [-45, 45] before converting it to radians.
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

	// remove the minus sign on -0.0 except for sin(-0.0).
	if x != 0 {
		sinx += 0
		cosx += 0
	}
	return sinx, cosx
}

// atan2d calculates atan2(y, x) with the result in degrees where y = the sine of the angle and x =
// the cosine of the angle. The result is in the range (-180 180]. N.B., atan2d(±0, -1) = +180;
// atan2d(-ε, -1) = -180, for ε positive and tiny; atan2d(±0, 1) = ±0.
func atan2d(y, x float64) float64 {
	// In order to minimize round-off errors, this function rearranges the arguments so that result
	// of atan2 is in the range [-pi/4, pi/4] before converting it to degrees and mapping the result
	// to the correct quadrant.
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
		if y >= 0 {
			ang = 180 - ang
		} else {
			ang = -180 - ang
		}
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
	return math.Abs(x) <= math.MaxFloat64
}

// min returns the minimum of two ints
func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

// sign returns 1 if x is non-negative, else -1
func sign(x float64) float64 {
	if x >= 0 {
		return 1
	} else {
		return -1
	}
}
