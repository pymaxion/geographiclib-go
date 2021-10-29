package geographiclib

// accumulator is an accumulator for sums.
//
// This allow many double precision numbers to be added together with twice the normal precision.
// Thus the effective precision of the sum is 106 bits or about 32 decimal places.
//
// The implementation follows J. R. Shewchuk, "Adaptive Precision Floating-Point Arithmetic and Fast
// Robust Geometric Predicates" (https://doi.org/10.1007/PL00009321), Discrete & Computational
// Geometry 18(3) 305-363 (1997).
//
// In the documentation of the struct's methods, sum stands for the value currently held
// in the accumulator.
type accumulator struct {
	s float64
	t float64
}

// newAccumulator constructs an accumulator from a float64
func newAccumulator(y float64) *accumulator {
	return &accumulator{
		s: y,
		t: 0,
	}
}

// set sets the value to a float64
func (a *accumulator) set(y float64) {
	a.s = y
	a.t = 0
}

// sum returns the value held in the accumulator
func (a *accumulator) sum() float64 {
	return a.s
}

// sumWith returns the result of adding a number to sum (but doesn't change sum)
func (a *accumulator) sumWith(y float64) float64 {
	b := *a
	b.add(y)
	return b.sum()
}

// add adds a number to the accumulator
func (a *accumulator) add(y float64) {
	// Here's Shewchuk's solution...
	y, u := sum(y, a.t) // accumulate stating at least significant end
	a.s, a.t = sum(y, a.s)
	// Start is s, t decreasing and non-adjacent. Sum is now (s + t + u) exactly with s, t, u
	// non-adjacent and in decreasing order (except for possible zeros). The following code tries to
	// normalize the result. Ideally, we want s = round(s+t+u) and u = round(s+t+u - s). The
	// following does an approximate job (and maintains the decreasing non-adjacent property). Here
	// are two "failures" using 3-bit floats:
	//
	// Case 1: s is not equal to round(s+t+u) -- off by 1 ulp
	// [12, -1] - 8 -> [4, 0, -1] -> [4, -1] = 3 should be [3, 0] = 3
	//
	// Case 2: s+t is not as close to s+t+u as it shold be
	// [64, 5] + 4 -> [64, 8, 1] -> [64,  8] = 72 (off by 1)
	//                    should be [80, -7] = 73 (exact)
	//
	// "Fixing" these problems is probably not worth the expense. The representation inevitably leads
	// to small errors in the accumulated values. The additional errors illustrated here amount to 1
	// ulp of the less significant word during each addition to the Accumulator and an additional
	// possible error of 1 ulp in the reported sum.
	//
	// Incidentally, the "ideal" representation described above is not canonical, because s =
	// round(s + t) may not be true. For example, with 3-bit floats:
	//
	// [128, 16] + 1 -> [160, -16] -- 160 = round(145).
	// But [160, 0] - 16 -> [128, 16] -- 128 = round(144).
	//
	if a.s == 0 { // This implies t == 0,
		a.s = u // so result is u
	} else {
		a.t += u // otherwise just accumulate u to t.
	}
}

// negate negates an accumulator
func (a *accumulator) negate() {
	a.s *= -1
	a.t *= -1
}

// remainder calculates the remainder of the sum on division by y
func (a *accumulator) remainder(y float64) {
	a.s = remainder(a.s, y)
	a.add(0) // renormalize
}
