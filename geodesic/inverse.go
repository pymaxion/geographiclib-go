package geodesic

import (
	"math"

	"github.com/pymaxion/geographiclib-go/geodesic/capabilities"
)

type inverseSolver struct {
	*Geodesic
}

func (s *inverseSolver) inverse(lat1, lon1, lat2, lon2 float64, caps capabilities.Mask) Data {
	caps &= capabilities.OutMask
	ir := s.genInverse(lat1, lon1, lat2, lon2, caps)
	if (caps & capabilities.Azimuth) != 0 {
		ir.Azi1, ir.Azi2 = atan2d(ir.salp1, ir.calp1), atan2d(ir.salp2, ir.calp2)
	}
	return ir.Data
}

type inverseResult struct {
	Data
	salp1 float64
	calp1 float64
	salp2 float64
	calp2 float64
}

func (s *inverseSolver) genInverse(lat1, lon1, lat2, lon2 float64, caps capabilities.Mask) inverseResult {
	r := inverseResult{
		Data:  newData(),
		salp1: math.NaN(),
		calp1: math.NaN(),
		salp2: math.NaN(),
		calp2: math.NaN(),
	}

	// Compute longitude difference (AngDiff does this carefully). Result is in
	// [-180, 180] but -180 is only for west-going geodesics. 180 is for east-going
	// and meridional geodesics.
	lat1, lat2 = latFix(lat1), latFix(lat2)
	r.Lat1, r.Lat2 = lat1, lat2

	// If really close to the equator, treat as on equator.
	lat1, lat2 = angRound(lat1), angRound(lat2)
	lon12, lon12s := angDiff(lon1, lon2)
	if (caps & capabilities.LongUnroll) != 0 {
		r.Lon1, r.Lon2 = lon1, (lon1+lon12)+lon12s
	} else {
		r.Lon1, r.Lon2 = angNormalize(lon1), angNormalize(lon2)
	}

	// Make longitude difference positive.
	lonSign := math.Copysign(1, lon12)
	lon12 *= lonSign
	lon12s *= lonSign
	lam12 := deg2rad(lon12)
	// Calculate sincos of lon12 + error (this applies angRound internally)
	slam12, clam12 := sincosde(lon12, lon12s)
	lon12s = (180 - lon12) - lon12s // the supplementary longitude difference

	// Swap points so that point with higher (abs) latitude is point 1
	// If one latitude is a nan, then it becomes lat1.
	swapp := ternary(math.Abs(lat1) < math.Abs(lat2) || lat2 != lat2, -1, 1)
	if swapp < 0 {
		lonSign *= -1
		lat1, lat2 = lat2, lat1
	}
	// Make lat1 <= 0
	latSign := math.Copysign(1, -lat1)
	lat1 *= latSign
	lat2 *= latSign
	// Now we have
	//
	//     0 <= lon12 <= 180
	//     -90 <= lat1 <= 0
	//     lat1 <= lat2 <= -lat1
	//
	// lonSign, swapp, latSign register the transformation to bring the coordinates
	// to this canonical form. In all cases, 1 means no change was made. We make
	// these transformations so that there are few cases to check, e.g., on verifying
	// quadrants in atan2. In addition, this enforces some symmetries in the results
	// returned.

	sbet1, cbet1 := sincosd(lat1)
	sbet1 *= s.f1
	// Ensure cbet1 = +epsilon at poles; doing the fix on beta means that sig12
	// will be <= 2*tiny for two points at the same pole.
	sbet1, cbet1 = norm(sbet1, cbet1)
	cbet1 = math.Max(tiny, cbet1)
	sbet2, cbet2 := sincosd(lat2)
	sbet2 *= s.f1
	// Ensure cbet2 = +epsilon at poles
	sbet2, cbet2 = norm(sbet2, cbet2)
	cbet2 = math.Max(tiny, cbet2)

	// If cbet1 < -sbet1, then cbet2 - cbet1 is a sensitive measure of the |bet1| -
	// |bet2|. Alternatively (cbet1 >= -sbet1), abs(sbet2) + sbet1 is a better
	// measure. This logic is used in assigning calp2 in Lambda12. Sometimes these
	// quantities vanish and in that case we force bet2 = +/- bet1 exactly. An
	// example where this is necessary is the inverse problem 48.522876735459 0
	// -48.52287673545898293 179.599720456223079643 which failed with Visual Studio
	// 10 (Release and Debug)
	if cbet1 < -sbet1 {
		if cbet2 == cbet1 {
			sbet2 = math.Copysign(sbet1, sbet2)
		}
	} else {
		if math.Abs(sbet2) == -sbet1 {
			cbet2 = cbet1
		}
	}
	dn1 := math.Sqrt(1 + s.ep2*sq(sbet1))
	dn2 := math.Sqrt(1 + s.ep2*sq(sbet2))

	s12x, m12x, sig12 := math.NaN(), math.NaN(), math.NaN()
	c1a := make([]float64, nC1+1)
	c2a := make([]float64, nC2+1)
	c3a := make([]float64, nC3)

	meridian := lat1 == -90 || slam12 == 0
	if meridian {
		// Endpoints are on a single full meridian, so the geodesic might lie on a meridian.
		r.calp1, r.salp1 = clam12, slam12 // Head to the target longitude
		r.calp2, r.salp2 = 1, 0           // At the target we're heading north
		// tan(bet) = tan(sig) * cos(alp)
		ssig1, csig1 := sbet1, r.calp1*cbet1
		ssig2, csig2 := sbet2, r.calp2*cbet2
		// sig12 = sig2 - sig1 (N.B., math.Max(+0.0, -0.0) -> +0.0)
		sig12 = math.Atan2(math.Max(0.0, csig1*ssig2-ssig1*csig2), csig1*csig2+ssig1*ssig2)

		lCaps := caps | capabilities.Distance | capabilities.ReducedLength
		lr := s.lengths(s.n, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2, cbet1, cbet2, lCaps, c1a, c2a)
		s12x, m12x, r.M12, r.M21 = lr.s12b, lr.m12b, lr.M12, lr.M21
		// Add the check for sig12 since zero length geodesics might yield m12 < 0. Test
		// case was:
		//
		//    echo 20.001 0 20.001 0 | GeodSolve -i
		//
		// In fact, we will have sig12 > pi/2 for meridional geodesic which is not a
		// shortest path.
		if sig12 < 1 || m12x >= 0 {
			// Need at least 2, to handle 90 0 90 180
			if sig12 < 3*tiny ||
				// Prevent negative s12 or m12 for short lines
				(sig12 < tol0 && (s12x < 0 || m12x < 0)) {
				sig12, m12x, s12x = 0, 0, 0
			}
			m12x *= s.b
			s12x *= s.b
			r.A12 = rad2deg(sig12)
		} else {
			// m12 < 0, i.e., prolate and too close to anti-podal
			meridian = false
		}
	}

	omg12, somg12, comg12 := math.NaN(), 2., math.NaN()
	if !meridian && sbet1 == 0 && // and sbet2 == 0
		// Mimic the way Lambda12 works with calp1 = 0
		(s.f <= 0 || lon12s >= s.f*180) {
		// Geodesic runs along equator
		r.calp1, r.calp2 = 0, 0
		r.salp1, r.salp2 = 1, 1
		s12x = s.a * lam12
		sig12 = lam12 / s.f1
		omg12 = sig12
		m12x = s.b * math.Sin(sig12)
		if (caps & capabilities.GeodesicScale) != 0 {
			r.M12 = math.Cos(sig12)
			r.M21 = r.M12
		}
		r.A12 = lon12 / s.f1
	} else if !meridian {
		// Now point1 and point2 belong within a hemisphere bounded by a
		// meridian and geodesic is neither meridional or equatorial.

		// Figure a starting point for Newton's method
		sr := s.inverseStart(sbet1, cbet1, dn1, sbet2, cbet2, dn2, lam12, slam12, clam12, c1a, c2a)
		sig12, r.salp1, r.calp1, r.salp2, r.calp2 = sr.sig12, sr.salp1, sr.calp1, sr.salp2, sr.calp2

		if sig12 >= 0 {
			// Short lines (InverseStart sets salp2, calp2, dnm)
			s12x = sig12 * s.b * sr.dnm
			m12x = sq(sr.dnm) * s.b * math.Sin(sig12/sr.dnm)
			if (caps & capabilities.GeodesicScale) != 0 {
				r.M12 = math.Cos(sig12 / sr.dnm)
				r.M21 = r.M12
			}
			r.A12 = rad2deg(sig12)
			omg12 = lam12 / (s.f1 * sr.dnm)
		} else {
			// Newton's method. This is a straightforward solution of f(alp1) =
			// lambda12(alp1) - lam12 = 0 with one wrinkle. f(alp) has exactly one root in
			// the interval (0, pi) and its derivative is positive at the root. Thus f(alp)
			// is positive for alp > alp1 and negative for alp < alp1. During the course of
			// the iteration, a range (alp1a, alp1b) is maintained which brackets the root
			// and with each evaluation of f(alp) the range is shrunk, if possible. Newton's
			// method is restarted whenever the derivative of f is negative (because the new
			// value of alp1 is then further from the solution) or if the new estimate of
			// alp1 lies outside (0,pi); in this case, the new starting guess is taken to be
			// (alp1a + alp1b) / 2.
			numit := 0
			// Bracketing range
			salp1a, calp1a := tiny, 1.
			salp1b, calp1b := tiny, -1.
			ssig1, csig1, ssig2, csig2, eps, domg12 := math.NaN(), math.NaN(), math.NaN(), math.NaN(), math.NaN(), math.NaN()

			for tripn, tripb := false, false; numit < maxit2; numit++ {
				// the WGS84 logic set: mean = 1.47, sr = 1.25, max = 16
				// WGS84 and random input: mean = 2.85, sr = 0.60
				l12r := s.lambda12(sbet1, cbet1, dn1, sbet2, cbet2, dn2, r.salp1, r.calp1, slam12, clam12, numit < maxit1, c1a, c2a, c3a)
				v, dv := l12r.lam12, l12r.dlam12
				sig12, ssig1, csig1, ssig2, csig2, eps, domg12 = l12r.sig12, l12r.ssig1, l12r.csig1, l12r.ssig2, l12r.csig2, l12r.eps, l12r.domg12
				r.salp2, r.calp2 = l12r.salp2, l12r.calp2

				// Reversed logic to allow escape with NaNs
				if tripb || !(math.Abs(v) >= ternary(tripn, 8, 1)*tol0) {
					break
				}
				// Update bracketing values
				if v > 0 && (numit > maxit1 || r.calp1/r.salp1 > calp1b/salp1b) {
					salp1b, calp1b = r.salp1, r.calp1
				} else if v < 0 && (numit > maxit1 || r.calp1/r.salp1 < calp1a/salp1a) {
					salp1a, calp1a = r.salp1, r.calp1
				}
				if numit < maxit1 && dv > 0 {
					dalp1 := -v / dv
					sdalp1, cdalp1 := math.Sincos(dalp1)
					nsalp1 := r.salp1*cdalp1 + r.calp1*sdalp1
					if nsalp1 > 0 && math.Abs(dalp1) < math.Pi {
						r.calp1 = r.calp1*cdalp1 - r.salp1*sdalp1
						r.salp1 = nsalp1
						r.salp1, r.calp1 = norm(r.salp1, r.calp1)
						// In some regimes we don't get quadratic convergence because
						// slope -> 0.  So use convergence conditions based on epsilon
						// instead of sqrt(epsilon).
						tripn = math.Abs(v) <= 16*tol0
						continue
					}
				}

				// Either dV was not positive or updated value was outside legal range. Use the
				// midpoint of the bracket as the next estimate. This mechanism is not needed for
				// the WGS84 ellipsoid, but it does catch problems with more eccentric
				// ellipsoids. Its efficacy is such for the WGS84 logic set with the starting
				// guess set to alp1 = 90deg:
				//  the WGS84 logic set: mean = 5.21, sr = 3.93, max = 24
				//  WGS84 and random input: mean = 4.74, sr = 0.99
				r.salp1, r.calp1 = (salp1a+salp1b)/2, (calp1a+calp1b)/2
				r.salp1, r.calp1 = norm(r.salp1, r.calp1)
				tripn = false
				tripb = math.Abs(salp1a-r.salp1)+(calp1a-r.calp1) < tolb ||
					math.Abs(r.salp1-salp1b)+(r.calp1-calp1b) < tolb
			}

			// Ensure that the reduced length and geodesic scale are computed in
			// a "canonical" way, with the I2 integral.
			var lMask capabilities.Mask
			if (caps & (capabilities.ReducedLength | capabilities.GeodesicScale)) != 0 {
				lMask = capabilities.OutMask | capabilities.Distance
			} else {
				lMask = capabilities.OutMask | capabilities.None
			}
			lr := s.lengths(eps, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2, cbet1, cbet2, lMask, c1a, c2a)
			s12x, m12x, r.M12, r.M21 = lr.s12b, lr.m12b, lr.M12, lr.M21
			m12x *= s.b
			s12x *= s.b
			r.A12 = rad2deg(sig12)
			if (caps & capabilities.Area) != 0 {
				// omg12 = lam12 - domg12
				sdomg12, cdomg12 := math.Sincos(domg12)
				somg12 = slam12*cdomg12 - clam12*sdomg12
				comg12 = clam12*cdomg12 + slam12*sdomg12
			}
		}
	} // end elif not meridian

	if (caps & capabilities.Distance) != 0 {
		r.S12 = 0 + s12x // Convert -0 to 0
	}
	if (caps & capabilities.ReducedLength) != 0 {
		r.M12Reduced = 0 + m12x // Convert -0 to 0
	}
	if (caps & capabilities.Area) != 0 {
		// From lambda12: sin(alp1) * cos(bet1) = sin(alp0)
		salp0 := r.salp1 * cbet1
		calp0 := math.Hypot(r.calp1, r.salp1*sbet1) // calp0 > 0
		var alp12 float64
		if calp0 != 0 && salp0 != 0 {
			// From lambda12: tan(bet) = tan(sig) * cos(alp)
			ssig1, csig1 := sbet1, r.calp1*cbet1
			ssig2, csig2 := sbet2, r.calp2*cbet2
			k2 := sq(calp0) * s.ep2
			eps := k2 / (2*(1+math.Sqrt(1+k2)) + k2)
			// Multiplier = a^2 * e^2 * cos(alpha0) * sin(alpha0).
			A4 := sq(s.a) * calp0 * salp0 * s.e2
			ssig1, csig1 = norm(ssig1, csig1)
			ssig2, csig2 = norm(ssig2, csig2)
			c4a := make([]float64, nC4)
			s.c4f(eps, c4a)

			B41, B42 := sinCosSeries(false, ssig1, csig1, c4a), sinCosSeries(false, ssig2, csig2, c4a)
			r.S12Area = A4 * (B42 - B41)
		} else {
			// Avoid problems with indeterminate sig1, sig2 on equator
			r.S12Area = 0
		}

		if !meridian && somg12 == 2 {
			somg12, comg12 = math.Sincos(omg12)
		}

		if !meridian && comg12 > -0.7071 && // Long difference not too big
			sbet2-sbet1 < 1.75 { // Lat difference not too big
			// Use tan(Gamma/2) = tan(omg12/2)
			// * (tan(bet1/2)+tan(bet2/2))/(1+tan(bet1/2)*tan(bet2/2))
			// with tan(x/2) = sin(x)/(1+cos(x))
			domg12 := 1 + comg12
			dbet1, dbet2 := 1+cbet1, 1+cbet2
			alp12 = 2 * math.Atan2(somg12*(sbet1*dbet2+sbet2*dbet1), domg12*(sbet1*sbet2+dbet1*dbet2))
		} else {
			// alp12 = alp2 - alp1, used in atan2 so no need to normalize
			salp12, calp12 := r.salp2*r.calp1-r.calp2*r.salp1, r.calp2*r.calp1+r.salp2*r.salp1
			// The right thing appears to happen if alp1 = +/-180 and alp2 = 0, viz salp12 =
			// -0 and alp12 = -180. However this depends on the sign being attached to 0
			// correctly. The following ensures the correct behavior.
			if salp12 == 0 && calp12 < 0 {
				salp12, calp12 = tiny*r.calp1, -1
			}
			alp12 = math.Atan2(salp12, calp12)
		}

		r.S12Area += s.c2 * alp12
		r.S12Area *= swapp * lonSign * latSign
		// Convert -0 to 0
		r.S12Area += 0
	}

	// Convert calp, salp to azimuth accounting for lonsign, swapp, latsign.
	if swapp < 0 {
		r.salp2, r.salp1 = r.salp1, r.salp2
		r.calp2, r.calp1 = r.calp1, r.calp2
		if (caps & capabilities.GeodesicScale) != 0 {
			r.M21, r.M12 = r.M12, r.M21
		}
	}

	r.salp1 *= swapp * lonSign
	r.calp1 *= swapp * latSign
	r.salp2 *= swapp * lonSign
	r.calp2 *= swapp * latSign
	return r
}

type lengthsResult struct {
	s12b float64
	m12b float64
	m0   float64
	M12  float64
	M21  float64
}

func (s *inverseSolver) lengths(eps, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2, cbet1, cbet2 float64, caps capabilities.Mask, c1a, c2a []float64) lengthsResult {
	r := lengthsResult{
		s12b: math.NaN(),
		m12b: math.NaN(),
		m0:   math.NaN(),
		M12:  math.NaN(),
		M21:  math.NaN(),
	}

	// Return m12b = (reduced length)/b; also calculate s12b = distance/b, and m0 =
	// coefficient of secular term in expression for reduced length.
	caps &= capabilities.OutMask
	var m0x, j12, a1, a2 float64
	if (caps & (capabilities.Distance | capabilities.ReducedLength | capabilities.GeodesicScale)) != 0 {
		a1 = a1m1f(eps)
		c1f(eps, c1a)
		if (caps & (capabilities.ReducedLength | capabilities.GeodesicScale)) != 0 {
			a2 = a2m1f(eps)
			c2f(eps, c2a)
			m0x = a1 - a2
			a2 = 1 + a2
		}
		a1 = 1 + a1
	}

	if (caps & capabilities.Distance) != 0 {
		b1 := sinCosSeries(true, ssig2, csig2, c1a) - sinCosSeries(true, ssig1, csig1, c1a)
		// Missing a factor of b
		r.s12b = a1 * (sig12 + b1)
		if (caps & (capabilities.ReducedLength | capabilities.GeodesicScale)) != 0 {
			b2 := sinCosSeries(true, ssig2, csig2, c2a) - sinCosSeries(true, ssig1, csig1, c2a)
			j12 = m0x*sig12 + (a1*b1 - a2*b2)
		}
	} else if (caps & (capabilities.ReducedLength | capabilities.GeodesicScale)) != 0 {
		// Assume here that nC1 >= nC2
		for l := 1; l <= nC2; l++ {
			c2a[l] = a1*c1a[l] - a2*c2a[l]
		}
		j12 = m0x*sig12 + (sinCosSeries(true, ssig2, csig2, c2a) - sinCosSeries(true, ssig1, csig1, c2a))
	}

	if (caps & capabilities.ReducedLength) != 0 {
		r.m0 = m0x
		// Missing a factor of b.
		// Add parens around (csig1 * ssig2) and (ssig1 * csig2) to ensure
		// accurate cancellation in the case of coincident points.
		r.m12b = dn2*(csig1*ssig2) - dn1*(ssig1*csig2) - csig1*csig2*j12
	}

	if (caps & capabilities.GeodesicScale) != 0 {
		csig12 := csig1*csig2 + ssig1*ssig2
		t := s.ep2 * (cbet1 - cbet2) * (cbet1 + cbet2) / (dn1 + dn2)
		r.M12 = csig12 + (t*ssig2-csig2*j12)*ssig1/dn1
		r.M21 = csig12 - (t*ssig1-csig1*j12)*ssig2/dn2
	}
	return r
}

type startResult struct {
	sig12 float64
	salp1 float64
	calp1 float64
	salp2 float64
	calp2 float64
	dnm   float64
}

// inverseStart returns a starting point for Newton's method in salp1 and calp1
// (function value is -1). If Newton's method doesn't need to be used, return
// also salp2 and calp2 and function value is sig12.
func (s *inverseSolver) inverseStart(sbet1, cbet1, dn1, sbet2, cbet2, dn2, lam12, slam12, clam12 float64, c1a, c2a []float64) startResult {
	r := startResult{
		sig12: -1,
		salp1: math.NaN(),
		calp1: math.NaN(),
		salp2: math.NaN(),
		calp2: math.NaN(),
		dnm:   math.NaN(),
	}

	// bet12 = bet2 - bet1 in [0, pi); bet12a = bet2 + bet1 in (-pi, 0]
	sbet12 := sbet2*cbet1 - cbet2*sbet1
	cbet12 := cbet2*cbet1 + sbet2*sbet1
	sbet12a := sbet2*cbet1 + cbet2*sbet1

	var somg12, comg12 float64
	shortline := cbet12 >= 0 && sbet12 < 0.5 && cbet2*lam12 < 0.5
	if shortline {
		sbetm2 := sq(sbet1 + sbet2)
		// sin((bet1+bet2)/2)^2
		// =  (sbet1 + sbet2)^2 / ((sbet1 + sbet2)^2 + (cbet1 + cbet2)^2)
		sbetm2 /= sbetm2 + sq(cbet1+cbet2)
		r.dnm = math.Sqrt(1 + s.ep2*sbetm2)
		omg12 := lam12 / (s.f1 * r.dnm)
		somg12, comg12 = math.Sincos(omg12)
	} else {
		somg12, comg12 = slam12, clam12
	}

	r.salp1 = cbet2 * somg12
	if comg12 >= 0 {
		r.calp1 = sbet12 + cbet2*sbet1*sq(somg12)/(1+comg12)
	} else {
		r.calp1 = sbet12a - cbet2*sbet1*sq(somg12)/(1-comg12)
	}
	ssig12 := math.Hypot(r.salp1, r.calp1)
	csig12 := sbet1*sbet2 + cbet1*cbet2*comg12

	if shortline && ssig12 < s.etol2 {
		// really short lines
		r.salp2 = cbet1 * somg12
		var t float64
		if comg12 >= 0 {
			t = sq(somg12) / (1 + comg12)
		} else {
			t = 1 - comg12
		}
		r.calp2 = sbet12 - cbet1*sbet2*t
		r.salp2, r.calp2 = norm(r.salp2, r.calp2)
		// Set return value
		r.sig12 = math.Atan2(ssig12, csig12)
	} else if math.Abs(s.n) > 0.1 || // Skip astroid calc if too eccentric
		csig12 >= 0 ||
		ssig12 >= 6*math.Abs(s.n)*math.Pi*sq(cbet1) {
		// Nothing to do, zeroth order spherical approximation is OK
	} else {
		// Scale lam12 and bet2 to x, y coordinate system where antipodal point
		// is at origin and singular point is at y = 0, x = -1.
		var x, y, lamscale, betscale float64
		lam12x := math.Atan2(-slam12, -clam12) // lam12 - pi
		if s.f >= 0 {                          // In fact f == 0 does not get here
			// x = dlong, y = dlat
			k2 := sq(sbet1) * s.ep2
			eps := k2 / (2*(1+math.Sqrt(1+k2)) + k2)
			lamscale = s.f * cbet1 * s.a3f(eps) * math.Pi
			betscale = lamscale * cbet1
			x = lam12x / lamscale
			y = sbet12a / betscale
		} else { // _f < 0
			// x = dlat, y = dlong
			cbet12a := cbet2*cbet1 - sbet2*sbet1
			bet12a := math.Atan2(sbet12a, cbet12a)
			// In the case of lon12 = 180, this repeats a calculation made in
			// Inverse.
			ls := s.lengths(s.n, math.Pi+bet12a, sbet1, -cbet1, dn1, sbet2, cbet2, dn2, cbet1, cbet2, capabilities.ReducedLength, c1a, c2a)
			t := cbet1 * cbet2 * ls.m0 * math.Pi
			x = -1 + ls.m12b/t
			if x < -0.01 {
				betscale = sbet12a / x
			} else {
				betscale = -s.f * sq(cbet1) * math.Pi
			}
			lamscale = betscale / cbet1
			y = lam12x / lamscale
		}

		if y > -tol1 && x > -1-xthresh {
			// strip near cut
			if s.f >= 0 {
				r.salp1 = math.Min(1.0, -x)
				r.calp1 = -math.Sqrt(1 - sq(r.salp1))
			} else {
				r.calp1 = math.Max(ternary(x > -tol1, 0., -1.), x)
				r.salp1 = math.Sqrt(1 - sq(r.calp1))
			}
		} else {
			// Estimate alp1, by solving the astroid problem.
			//
			// Could estimate alpha1 = theta + pi/2, directly, i.e.,
			//   calp1 = y/k; salp1 = -x/(1+k);  for _f >= 0
			//   calp1 = x/(1+k); salp1 = -y/k;  for _f < 0 (need to check)
			//
			// However, it's better to estimate omg12 from astroid and use spherical formula
			// to compute alp1. This reduces the mean number of Newton iterations for astroid
			// cases from 2.24 (min 0, max 6) to 2.12 (min 0 max 5). The changes in the
			// number of iterations are as follows:
			//
			// change percent
			//    1       5
			//    0      78
			//   -1      16
			//   -2       0.6
			//   -3       0.04
			//   -4       0.002
			//
			// The histogram of iterations is (m = number of iterations estimating alp1
			// directly, n = number of iterations estimating via omg12, total number of
			// trials = 148605):
			//
			//  iter    m      n
			//    0   148    186
			//    1 13046  13845
			//    2 93315 102225
			//    3 36189  32341
			//    4  5396      7
			//    5   455      1
			//    6    56      0
			//
			// Because omg12 is near pi, estimate work with omg12a = pi - omg12
			k := s.astroid(x, y)
			var t float64
			if s.f >= 0 {
				t = -x * k / (1 + k)
			} else {
				t = -y * (1 + k) / k
			}
			omg12a := lamscale * t
			somg12, comg12 = math.Sincos(omg12a)
			comg12 *= -1
			// Update spherical estimate of alp1 using omg12 instead of lam12
			r.salp1 = cbet2 * somg12
			r.calp1 = sbet12a - cbet2*sbet1*sq(somg12)/(1-comg12)
		}
	}

	// Sanity check on starting guess. Backwards check allows NaN through.
	if !(r.salp1 <= 0) {
		r.salp1, r.calp1 = norm(r.salp1, r.calp1)
	} else {
		r.salp1, r.calp1 = 1, 0
	}
	return r
}

// astroid solves the astroid equation k^4+2*k^3-(x^2+y^2-1)*k^2-2*y^2*k-y^2 = 0
// for positive root k. This solution is adapted from Geocentric::Reverse.
func (s *inverseSolver) astroid(x, y float64) float64 {
	var k float64
	p, q := sq(x), sq(y)
	r := (p + q - 1) / 6
	if !(q == 0 && r <= 0) {
		// Avoid possible division by zero when r = 0 by multiplying equations for s and
		// t by r^3 and r, resp.
		S := p * q / 4 // S = r^3 * s
		r2 := sq(r)
		r3 := r * r2
		// The discriminant of the quadratic equation for T3. This is zero on the evolute
		// curve p^(1/3)+q^(1/3) = 1
		disc := S * (S + 2*r3)
		u := r
		if disc >= 0 {
			T3 := S + r3
			// Pick the sign on the sqrt to maximize abs(T3). This minimizes loss of
			// precision due to cancellation. The result is unchanged because of the way the
			// T is used in definition of u.
			j := math.Sqrt(disc)
			T3 += ternary(T3 < 0, -j, j) // T3 = (r * t)^3
			// N.B. cbrt always returns the root.  cbrt(-8) = -2.
			T := math.Cbrt(T3) // T = r * t
			// T can be zero; but then r2 / T -> 0.
			if T != 0 {
				j = r2 / T
			} else {
				j = 0
			}
			u += T + j
		} else {
			// T is complex, but the way u is defined the result is .
			ang := math.Atan2(math.Sqrt(-disc), -(S + r3))
			// There are three possible cube roots. We choose the root which avoids
			// cancellation. Note that disc < 0 implies that r < 0.
			u += 2 * r * math.Cos(ang/3)
		}
		v := math.Sqrt(sq(u) + q) // guaranteed positive
		// Avoid loss of accuracy when u < 0.
		var uv float64
		if u < 0 {
			uv = q / (v - u)
		} else {
			uv = u + v // u+v, guaranteed positive
		}
		w := (uv - q) / (2 * v) // positive?
		// Rearrange expression for k to avoid loss of accuracy due to subtraction.
		// Division by 0 not possible because uv > 0, w >= 0.
		k = uv / (math.Sqrt(uv+sq(w)) + w) // guaranteed positive
	} else { // q == 0 && r <= 0
		// y = 0 with |x| <= 1.  Handle this case directly.
		// for y small, positive root is k = abs(y)/sqrt(1-x^2)
		k = 0
	}
	return k
}

type lambda12Result struct {
	lam12  float64
	salp2  float64
	calp2  float64
	sig12  float64
	ssig1  float64
	csig1  float64
	ssig2  float64
	csig2  float64
	eps    float64
	domg12 float64
	dlam12 float64
}

func (s *inverseSolver) lambda12(sbet1, cbet1, dn1, sbet2, cbet2, dn2, salp1, calp1, slam120, clam120 float64, diffp bool, c1a, c2a, c3a []float64) lambda12Result {
	r := lambda12Result{}
	if sbet1 == 0 && calp1 == 0 {
		// Break degeneracy of equatorial line. This case has already been handled.
		calp1 = -tiny
	}

	// sin(alp1) * cos(bet1) = sin(alp0)
	salp0 := salp1 * cbet1
	calp0 := math.Hypot(calp1, salp1*sbet1) // calp0 > 0
	// tan(bet1) = tan(sig1) * cos(alp1)
	// tan(omg1) = sin(alp0) * tan(sig1) = tan(omg1)=tan(alp1)*sin(bet1)
	r.ssig1 = sbet1
	somg1 := salp0 * sbet1
	comg1 := calp1 * cbet1
	r.csig1 = comg1
	r.ssig1, r.csig1 = norm(r.ssig1, r.csig1)
	// norm(somg1, comg1) -- don't need to normalize!

	// Enforce symmetries in the case abs(bet2) = -bet1.  Need to be careful
	// about this case, since this can yield singularities in the Newton
	// iteration.
	// sin(alp2) * cos(bet2) = sin(alp0)
	if cbet2 != cbet1 {
		r.salp2 = salp0 / cbet2
	} else {
		r.salp2 = salp1
	}
	// calp2 = sqrt(1 - sq(salp2))
	//       = sqrt(sq(calp0) - sq(sbet2)) / cbet2
	// and subst for calp0 and rearrange to give (choose positive sqrt to give alp2
	// in [0, pi/2]).
	if cbet2 != cbet1 || math.Abs(sbet2) != -sbet1 {
		var t float64
		if cbet1 < -sbet1 {
			t = (cbet2 - cbet1) * (cbet1 + cbet2)
		} else {
			t = (sbet1 - sbet2) * (sbet1 + sbet2)
		}
		r.calp2 = math.Sqrt(sq(calp1*cbet1)+t) / cbet2
	} else {
		r.calp2 = math.Abs(calp1)
	}
	// tan(bet2) = tan(sig2) * cos(alp2)
	// tan(omg2) = sin(alp0) * tan(sig2).
	r.ssig2 = sbet2
	somg2 := salp0 * sbet2
	comg2 := r.calp2 * cbet2
	r.csig2 = comg2
	r.ssig2, r.csig2 = norm(r.ssig2, r.csig2)
	// norm(somg2, comg2); -- don't need to normalize!

	// sig12 = sig2 - sig1, limit to [0, pi]
	y := r.csig1*r.ssig2 - r.ssig1*r.csig2
	x := r.csig1*r.csig2 + r.ssig1*r.ssig2
	r.sig12 = math.Atan2(math.Max(0.0, y), x)
	// omg12 = omg2 - omg1, limit to [0, pi]
	somg12 := math.Max(0.0, comg1*somg2-somg1*comg2)
	comg12 := comg1*comg2 + somg1*somg2
	// eta = omg12 - lam120
	y = somg12*clam120 - comg12*slam120
	x = comg12*clam120 + somg12*slam120
	eta := math.Atan2(y, x)
	k2 := sq(calp0) * s.ep2
	r.eps = k2 / (2*(1+math.Sqrt(1+k2)) + k2)
	s.c3f(r.eps, c3a)
	B312 := sinCosSeries(true, r.ssig2, r.csig2, c3a) - sinCosSeries(true, r.ssig1, r.csig1, c3a)
	r.domg12 = -s.f * s.a3f(r.eps) * salp0 * (r.sig12 + B312)
	r.lam12 = eta + r.domg12

	if diffp {
		if r.calp2 == 0 {
			r.dlam12 = -2 * s.f1 * dn1 / sbet1
		} else {
			lr := s.lengths(r.eps, r.sig12, r.ssig1, r.csig1, dn1, r.ssig2, r.csig2, dn2, cbet1, cbet2, capabilities.ReducedLength, c1a, c2a)
			r.dlam12 = lr.m12b
			r.dlam12 *= s.f1 / (r.calp2 * cbet2)
		}
	} else {
		r.dlam12 = math.NaN()
	}
	return r
}
